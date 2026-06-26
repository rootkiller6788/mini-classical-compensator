/**
 * @file cascade_performance.c
 * @brief Performance evaluation and metrics for cascade control systems
 *
 * L2 -- Core Concepts: Performance metrics (ISE, IAE, ITAE), disturbance
 *       rejection quantification, cascade efficiency ratio.
 * L5 -- Computational Methods: Step response analysis, performance indices
 *       computation, settling time and overshoot extraction.
 * L6 -- Canonical Systems: Performance assessment for standard cascade
 *       configurations (DC motor, reactor, flow-pressure, level-tank).
 * L7 -- Applications: Comparative analysis of cascade vs. single-loop,
 *       demonstrating the disturbance rejection advantage.
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * L5: Step Response Performance Metrics
 * ========================================================================== */

int cascade_compute_performance(const CascadeTF *cl_tf,
                                 double t_final, int n_pts,
                                 CascadePerformance *perf) {
    /* Compute comprehensive performance metrics from step response.
     *
     * ISE = integral(e^2 dt) ? penalizes large errors heavily
     * IAE = integral(|e| dt) ? linear penalty, easier to minimize
     * ITAE = integral(t*|e| dt) ? penalizes late errors more
     *
     * Settling time: time for response to stay within 2% of final value
     * Overshoot: (y_max - y_ss) / y_ss
     * Rise time: 10% to 90% of final value
     *
     * Reference: Astrom & Hagglund, "PID Controllers" (1995), Ch. 2. */
    if (!cl_tf || !perf || t_final <= 0 || n_pts < 10) return -1;

    double *t = (double *)malloc((size_t)n_pts * sizeof(double));
    double *y = (double *)malloc((size_t)n_pts * sizeof(double));
    if (!t || !y) { free(t); free(y); return -1; }

    /* Generate time vector and step response */
    for (int i = 0; i < n_pts; i++) {
        t[i] = (double)i * t_final / (double)(n_pts - 1);
        y[i] = cascade_step_response(cl_tf, t[i], 15);
    }

    double y_ss = y[n_pts - 1];  /* Assume steady-state reached */
    if (fabs(y_ss) < 1e-10) y_ss = 1.0;

    /* Compute performance indices */
    double ise = 0.0, iae = 0.0, itae = 0.0;
    double y_max = y[0];
    int peak_idx = 0;
    double y_10 = 0.1 * y_ss, y_90 = 0.9 * y_ss;
    double t_rise_start = -1.0, t_rise_end = -1.0;
    double t_settle = t_final;

    for (int i = 0; i < n_pts; i++) {
        double e = 1.0 - y[i];
        double dt = (i > 0) ? t[i] - t[i - 1] : t_final / (double)n_pts;

        ise += e * e * dt;
        iae += fabs(e) * dt;
        itae += t[i] * fabs(e) * dt;

        /* Track peak */
        if (y[i] > y_max) {
            y_max = y[i];
            peak_idx = i;
        }

        /* Rise time detection */
        if (t_rise_start < 0 && y[i] >= y_10) t_rise_start = t[i];
        if (t_rise_end < 0 && y[i] >= y_90) t_rise_end = t[i];

        /* Settling time: last time response leaves 2% band */
        double band = 0.02 * y_ss;
        if (fabs(y[i] - y_ss) > band) t_settle = t[i];
    }

    /* Overshoot */
    double overshoot = 0.0;
    if (y_max > y_ss && y_ss > 1e-10) {
        overshoot = (y_max - y_ss) / fabs(y_ss);
    }

    /* Decay ratio (ratio of second peak to first peak) */
    double decay_ratio = 0.0;
    if (peak_idx > 0 && peak_idx < n_pts - 1) {
        /* Simple estimate from overshoot and damping */
        double zeta_est = -log(overshoot + 1e-10) /
                          sqrt(M_PI * M_PI + log(overshoot + 1e-10) * log(overshoot + 1e-10));
        decay_ratio = overshoot * overshoot;  /* For second-order systems */
        (void)zeta_est;
    }

    /* Fill performance struct */
    perf->outer_rise_time = (t_rise_end > 0 && t_rise_start > 0)
                            ? t_rise_end - t_rise_start : t_final;
    perf->outer_settle_time = t_settle;
    perf->outer_overshoot = overshoot;
    perf->outer_peak_time = (peak_idx > 0) ? t[peak_idx] : 0.0;
    perf->outer_ess = fabs(1.0 - y_ss);
    perf->outer_decay_ratio = decay_ratio;

    /* These would be filled by comparing with inner loop separately */
    perf->inner_rise_time = 0.0;
    perf->inner_settle_time = 0.0;
    perf->inner_overshoot = 0.0;
    perf->inner_ess = 0.0;

    perf->dist_rejection_ratio = 0.0;
    perf->cascade_efficiency = 0.0;
    perf->robustness_margin = 1.0;

    free(t); free(y);
    return 0;
}

/* ==========================================================================
 * L2: Disturbance Rejection Quantification
 * ========================================================================== */

int cascade_compare_single_vs_cascade(const CascadeTF *inner_model,
                                       const CascadeTF *outer_model,
                                       const CascadePID *ci,
                                       const CascadePID *co,
                                       double t_final, int n_pts,
                                       double *improvement_factor) {
    /* Quantify the advantage of cascade over single-loop control.
     *
     * Setup: A disturbance d enters the inner loop.
     *   Cascade: inner loop rejects d before it affects outer PV significantly
     *   Single-loop: only outer controller responds, much slower
     *
     * Improvement factor = ISE_single / ISE_cascade
     * Typical values: 2-10x improvement for well-designed cascade
     *
     * Theorem (Disturbance Rejection of Cascade):
     * For disturbance d entering inner PV:
     *   y_cascade/d = Go * S_i * S_o
     *   y_single/d = Go * S_single
     * where S_i = 1/(1+Ci*Gi) is inner sensitivity.
     * Since |S_i(jw)| << 1 for w << wb_i, cascade provides order-of-magnitude
     * improvement for disturbances within inner loop bandwidth. */
    if (!inner_model || !outer_model || !ci || !co || !improvement_factor)
        return -1;

    /* Cascade: inner CL rejects disturbance */
    CascadeTF Gi_cl;
    cascade_inner_closed_loop(inner_model, ci, &Gi_cl);

    CascadeTF cascade_cl;
    cascade_overall_closed_loop(&Gi_cl, co, outer_model, &cascade_cl);

    /* Single-loop: only outer controller, no inner loop */
    CascadeTF co_tf;
    cascade_pid_to_tf(co, &co_tf);
    CascadePoly num_sl = cascade_poly_mul(&co_tf.num, &outer_model->num);
    CascadePoly den_sl = cascade_poly_mul(&co_tf.den, &outer_model->den);

    /* Close outer loop directly (simplified single-loop) */
    CascadeTF single_cl_tf;
    CascadePoly den_cl_sl = cascade_poly_add(&den_sl, &num_sl);
    single_cl_tf.num = num_sl;
    single_cl_tf.den = den_cl_sl;
    single_cl_tf.gain = 1.0;
    single_cl_tf.delay = 0.0;
    single_cl_tf.has_delay = 0;

    /* Compute ISE for both */
    double ise_cascade = 0.0, ise_single = 0.0;
    for (int i = 0; i < n_pts; i++) {
        double t = (double)i * t_final / (double)(n_pts - 1);
        /* Disturbance step response (simplified: y ~ 1 - step_response) */
        double y_c = cascade_step_response(&cascade_cl, t, 10);
        double y_s = cascade_step_response(&single_cl_tf, t, 10);
        double dt = t_final / (double)n_pts;
        ise_cascade += (1.0 - y_c) * (1.0 - y_c) * dt;
        ise_single += (1.0 - y_s) * (1.0 - y_s) * dt;
    }

    *improvement_factor = (ise_cascade > 1e-15)
                          ? ise_single / ise_cascade : 1.0;

    cascade_tf_free(&Gi_cl);
    cascade_tf_free(&cascade_cl);
    cascade_tf_free(&co_tf);
    cascade_poly_free(&den_sl);
    cascade_tf_free(&single_cl_tf);
    return 0;
}

/* ==========================================================================
 * L5: Performance Index Computation
 * ========================================================================== */

double cascade_compute_ise(const double *error, const double *time, int n) {
    /* ISE = Integral of Squared Error
     * Numerically integrated using trapezoidal rule. */
    if (!error || !time || n < 2) return -1.0;
    double ise = 0.0;
    for (int i = 1; i < n; i++) {
        double dt = time[i] - time[i - 1];
        double e_avg = (error[i] * error[i] + error[i - 1] * error[i - 1]) / 2.0;
        ise += e_avg * dt;
    }
    return ise;
}

double cascade_compute_iae(const double *error, const double *time, int n) {
    /* IAE = Integral of Absolute Error */
    if (!error || !time || n < 2) return -1.0;
    double iae = 0.0;
    for (int i = 1; i < n; i++) {
        double dt = time[i] - time[i - 1];
        double e_avg = (fabs(error[i]) + fabs(error[i - 1])) / 2.0;
        iae += e_avg * dt;
    }
    return iae;
}

double cascade_compute_itae(const double *error, const double *time, int n) {
    /* ITAE = Integral of Time-weighted Absolute Error
     * Penalizes errors that persist at later times. */
    if (!error || !time || n < 2) return -1.0;
    double itae = 0.0;
    for (int i = 1; i < n; i++) {
        double dt = time[i] - time[i - 1];
        double e_avg = (time[i] * fabs(error[i]) +
                        time[i - 1] * fabs(error[i - 1])) / 2.0;
        itae += e_avg * dt;
    }
    return itae;
}

double cascade_compute_tv(const double *control, int n) {
    /* Total Variation of control signal.
     * TV = sum |u(k) - u(k-1)|
     * Measures control effort / actuator wear.
     * Higher TV means more aggressive control action. */
    if (!control || n < 2) return 0.0;
    double tv = 0.0;
    for (int i = 1; i < n; i++) {
        tv += fabs(control[i] - control[i - 1]);
    }
    return tv;
}

/* ==========================================================================
 * L5: Settling Time and Overshoot Extraction
 * ========================================================================== */

int cascade_extract_rise_time(const double *t, const double *y, int n,
                               double y_final, double *rise_time) {
    /* 10%-90% rise time extraction from step response data. */
    if (!t || !y || !rise_time || n < 3 || fabs(y_final) < 1e-10)
        return -1;

    double y10 = 0.10 * y_final;
    double y90 = 0.90 * y_final;
    double t10 = -1.0, t90 = -1.0;

    for (int i = 1; i < n; i++) {
        if (t10 < 0 && y[i] >= y10 && y[i - 1] < y10) {
            /* Linear interpolation */
            double frac = (y10 - y[i - 1]) / (y[i] - y[i - 1] + 1e-15);
            t10 = t[i - 1] + frac * (t[i] - t[i - 1]);
        }
        if (t90 < 0 && y[i] >= y90 && y[i - 1] < y90) {
            double frac = (y90 - y[i - 1]) / (y[i] - y[i - 1] + 1e-15);
            t90 = t[i - 1] + frac * (t[i] - t[i - 1]);
        }
        if (t10 >= 0 && t90 >= 0) break;
    }

    if (t10 >= 0 && t90 >= 0) {
        *rise_time = t90 - t10;
        return 0;
    }
    return -1;
}

int cascade_extract_settling_time(const double *t, const double *y, int n,
                                   double y_final, double band_pct,
                                   double *settle_time) {
    /* 2% (or band_pct%) settling time: last time response leaves the band. */
    if (!t || !y || !settle_time || n < 3) return -1;

    double band = (band_pct / 100.0) * fabs(y_final);
    if (band < 1e-10) band = 0.02 * fabs(y_final);
    if (band < 1e-10) band = 0.01;

    double t_settle = t[n - 1];
    for (int i = n - 1; i >= 0; i--) {
        if (fabs(y[i] - y_final) > band) {
            t_settle = t[i];
            break;
        }
    }
    *settle_time = t_settle;
    return 0;
}

double cascade_extract_overshoot(const double *y, int n, double y_final) {
    /* Maximum overshoot as fraction of final value. */
    if (!y || n < 2 || fabs(y_final) < 1e-10) return 0.0;
    double y_max = y[0];
    for (int i = 0; i < n; i++) {
        if (y[i] > y_max) y_max = y[i];
    }
    if (y_max > y_final) return (y_max - y_final) / fabs(y_final);
    return 0.0;
}

/* ==========================================================================
 * L6: Application-Specific Performance Assessment
 * ========================================================================== */

int cascade_assess_dc_motor(const CascadeDCMotor *motor,
                             CascadePerformance *perf) {
    /* Assess cascade performance for DC motor position/velocity control.
     *
     * Typical cascade structure:
     *   Inner: PI velocity control (current/torque loop may be innermost)
     *   Outer: P position control
     *
     * Performance expectations:
     *   - Velocity loop bandwidth: 50-500 Hz for servo motors
     *   - Position loop bandwidth: 5-50 Hz
     *   - Bandwidth ratio: typically 5-10
     *   - Overshoot: <5% for precision positioning */
    if (!motor || !perf) return -1;

    /* Compute velocity loop bandwidth (inner) */
    double tau_e = motor->L / (motor->R > 0.01 ? motor->R : 1.0);
    double tau_m = motor->J * motor->R /
                   (motor->Kb * motor->Kt + motor->B * motor->R + 1e-10);

    /* Inner bandwidth estimate: 1/(tau_e + tau_m) */
    double wb_inner = 1.0 / (tau_e + tau_m + 1e-6);
    perf->inner_settle_time = 4.0 / (wb_inner + 1e-6);
    perf->inner_rise_time = 2.2 / (wb_inner + 1e-6);
    perf->inner_overshoot = 0.02;  /* Well-tuned PI gives ~2% overshoot */
    perf->inner_ess = 0.0;  /* PI control: zero steady-state error */

    /* Outer position loop: typically 5-10x slower */
    double wb_outer = wb_inner / 5.0;
    perf->outer_settle_time = 4.0 / (wb_outer + 1e-6);
    perf->outer_rise_time = 2.2 / (wb_outer + 1e-6);
    perf->outer_overshoot = 0.05;
    perf->outer_ess = 0.0;

    /* Disturbance rejection improvement */
    perf->dist_rejection_ratio = 5.0;
    perf->cascade_efficiency = perf->dist_rejection_ratio / 5.0;
    perf->robustness_margin = 1.5;  /* Ms ~ 1.5 */

    return 0;
}

int cascade_assess_reactor(const CascadeReactor *rx,
                            CascadePerformance *perf) {
    /* Assess cascade performance for jacketed reactor temperature control.
     *
     * Cascade: jacket flow -> jacket temperature (inner)
     *          jacket temperature -> reactor temperature (outer)
     *
     * Typical bandwidth ratio: 3-5 (limited by heat transfer dynamics)
     * Key challenge: reactor may be open-loop unstable at some operating
     * points, requiring stabilizing outer controller. */
    if (!rx || !perf) return -1;

    /* Jacket time constant (inner) */
    double F_j_safe = (rx->F_j > 1e-10) ? rx->F_j : 1e-10;
    double tau_j = rx->V_j / F_j_safe;

    /* Reactor time constant (outer) */
    double F_in_safe = (rx->F_in > 1e-10) ? rx->F_in : 1e-10;
    double tau_r = rx->V_r / F_in_safe;

    perf->inner_settle_time = 4.0 * tau_j;
    perf->inner_rise_time = 2.2 * tau_j;
    perf->outer_settle_time = 4.0 * tau_r;
    perf->outer_rise_time = 2.2 * tau_r;
    perf->inner_overshoot = 0.05;
    perf->outer_overshoot = 0.10;
    perf->dist_rejection_ratio = 3.0;
    perf->cascade_efficiency = 3.0 / 4.0;
    perf->robustness_margin = 1.6;
    return 0;
}

int cascade_assess_level_tank(const CascadeLevelTank *tank,
                               CascadePerformance *perf) {
    /* Assess cascade performance for level-on-flow surge tank control.
     *
     * Cascade: outflow valve -> outflow rate (inner)
     *          outflow rate -> tank level (outer, integrating)
     *
     * Since the outer process is integrating (1/(A*s)), the cascade
     * inner loop is critical for disturbance rejection. Without cascade,
     * level control is very sluggish due to the integration.
     *
     * Typical bandwidth ratio: 5-10 (higher because outer is integrating) */
    if (!tank || !perf) return -1;

    /* Flow loop (inner): fast, ~seconds */
    double tau_flow = 0.5;
    perf->inner_settle_time = 4.0 * tau_flow;
    perf->inner_rise_time = 2.2 * tau_flow;
    perf->inner_overshoot = 0.02;
    perf->inner_ess = 0.0;

    /* Level loop (outer): integrating, ~minutes */
    double A = (tank->tank_area > 0.01) ? tank->tank_area : 1.0;
    /* Level time constant for P-control: tau_level ~ A * level_range / flow */
    double tau_level = A * (tank->max_level - tank->min_level) /
                       (tank->pump_max_flow + 1e-10);
    if (tau_level < 1.0) tau_level = 10.0;

    perf->outer_settle_time = 4.0 * tau_level;
    perf->outer_rise_time = 2.2 * tau_level;
    perf->outer_overshoot = 0.0;  /* Integrating: no overshoot with P-only */
    perf->outer_ess = 0.0;  /* Integrating: zero error for step with P control */

    perf->dist_rejection_ratio = 8.0;  /* High benefit for integrating outer */
    perf->cascade_efficiency = 8.0 / 10.0;
    perf->robustness_margin = 1.8;
    return 0;
}
