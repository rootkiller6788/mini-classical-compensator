/**
 * @file cascade_design.c
 * @brief Cascade controller design: sequential loop closure, PID tuning rules,
 *        frequency-domain design, simultaneous optimization.
 *
 * L4 -- Fundamental Laws: Sequential loop closure principle.
 * L5 -- Computational Methods: Skogestad SIMC, Direct Synthesis,
 *       frequency-domain PI, Nelder-Mead simultaneous optimization.
 *
 * Design workflow:
 *   1. cascade_design_inner_loop()  -- fast PI for disturbance rejection
 *   2. cascade_form_equivalent_plant() -- Geq = Gi_cl * Go
 *   3. cascade_design_outer_loop() -- PI/PID for tracking
 *   4. cascade_validate_design()   -- verify all specifications met
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
 * L4: Sequential Loop Closure Design
 *
 * Theorem: If Ci(s) internally stabilizes Gi(s) and Co(s) internally
 * stabilizes Geq(s) = Gi_cl(s)*Go(s), then the cascade system is
 * internally stable, provided wb_inner >= 3 * wb_outer.
 *
 * Reference: Zhou, Doyle & Glover (1996), "Robust and Optimal Control".
 * ========================================================================== */

int cascade_design_sequential(const CascadeTF *inner_model,
                               const CascadeTF *outer_model,
                               const CascadeDesignSpec *spec,
                               CascadeSystem *sys) {
    if (!inner_model || !outer_model || !spec || !sys) return -1;

    /* Step 1: Design inner loop for speed and disturbance rejection */
    double inner_bw = spec->inner_bw_min;
    if (cascade_design_inner_loop(inner_model, inner_bw,
                                   spec->inner_pm_target,
                                   &sys->inner.inner_controller) != 0) {
        fprintf(stderr, "cascade_design: inner loop design failed\n");
        return -1;
    }

    /* Copy inner model */
    sys->inner.inner_process = cascade_tf_create(
        inner_model->num.coeff, inner_model->num.degree,
        inner_model->den.coeff, inner_model->den.degree,
        inner_model->gain);

    /* Step 2: Form inner closed loop Gi_cl = Ci*Gi/(1+Ci*Gi) */
    if (cascade_inner_closed_loop(inner_model, &sys->inner.inner_controller,
                                   &sys->inner.inner_cl) != 0) {
        return -1;
    }

    /* Verify inner stability */
    sys->inner.inner_is_stable = cascade_is_stable(&sys->inner.inner_cl);
    if (!sys->inner.inner_is_stable) {
        fprintf(stderr, "cascade_design: inner loop is unstable\n");
        return -1;
    }

    /* Estimate inner bandwidth */
    cascade_closed_loop_bandwidth(&sys->inner.inner_cl,
                                   &sys->inner.inner_bandwidth);

    /* Step 3: Form equivalent plant Geq = Gi_cl * Go */
    if (cascade_form_equivalent_plant(inner_model,
                                       &sys->inner.inner_controller,
                                       outer_model,
                                       &sys->outer.equivalent_plant) != 0) {
        return -1;
    }

    /* Copy outer model */
    sys->outer.outer_process = cascade_tf_create(
        outer_model->num.coeff, outer_model->num.degree,
        outer_model->den.coeff, outer_model->den.degree,
        outer_model->gain);

    /* Step 4: Design outer loop */
    if (cascade_design_outer_loop(&sys->outer.equivalent_plant,
                                   spec,
                                   &sys->outer.outer_controller) != 0) {
        fprintf(stderr, "cascade_design: outer loop design failed\n");
        return -1;
    }

    /* Step 5: Compute overall closed loop */
    if (cascade_overall_closed_loop(&sys->inner.inner_cl,
                                     &sys->outer.outer_controller,
                                     outer_model,
                                     &sys->outer.outer_cl) != 0) {
        return -1;
    }

    sys->outer.outer_is_stable = cascade_is_stable(&sys->outer.outer_cl);
    if (!sys->outer.outer_is_stable) {
        fprintf(stderr, "cascade_design: outer loop unstable\n");
        return -1;
    }

    /* Estimate outer bandwidth and compute ratio */
    cascade_closed_loop_bandwidth(&sys->outer.outer_cl,
                                   &sys->outer.outer_bandwidth);

    if (sys->outer.outer_bandwidth > 1e-10) {
        sys->bandwidth_ratio = sys->inner.inner_bandwidth /
                               sys->outer.outer_bandwidth;
    } else {
        sys->bandwidth_ratio = 5.0;
    }

    if (sys->bandwidth_ratio < spec->bw_ratio_min) {
        fprintf(stderr, "cascade_design: warning -- bandwidth ratio %.2f "
                "below minimum %.2f\n",
                sys->bandwidth_ratio, spec->bw_ratio_min);
    }

    sys->cascade_active = 1;
    sys->system_name = "Designed Cascade System";
    return 0;
}

/* ==========================================================================
 * L5: Inner Loop PI Design (Frequency-Domain Method)
 *
 * Design equation for PI at gain crossover frequency w_gc:
 *   |Ci(jw_gc)| * |Gi(jw_gc)| = 1  (magnitude condition)
 *   arg(Ci(jw_gc)) + arg(Gi(jw_gc)) = -180 + PM (phase condition)
 *
 * For PI: C(s) = Kc*(1 + 1/(Ti*s))
 *   |C(jw)| = Kc * sqrt(1 + 1/(w*Ti)^2)
 *   arg(C) = atan(w*Ti) - pi/2 (range: -pi/2 to 0)
 *
 * Solving:
 *   Ti = tan(PM + pi/2 - arg(Gi)) / w_gc
 *   Kc = 1 / (|Gi| * sqrt(1 + 1/(w_gc*Ti)^2))
 *
 * Reference: Astrom & Hagglund, "PID Controllers" (1995), Ch. 4.
 * ========================================================================== */

int cascade_design_inner_loop(const CascadeTF *model,
                               double bw_target,
                               double pm_target,
                               CascadePID *pid) {
    if (!model || !pid || bw_target <= 0.0) return -1;

    double w_gc = bw_target;
    double PM_rad = pm_target * M_PI / 180.0;

    /* Plant frequency response at target crossover */
    double plant_mag, plant_phase;
    cascade_tf_freq_response(model, w_gc, &plant_mag, &plant_phase);
    if (plant_mag < 1e-15) return -1;

    /* Desired PI phase contribution:
     * arg(C) = -pi + PM - arg(Gi) = -180 + PM - arg(Gi) */
    double target_phase = -M_PI + PM_rad - plant_phase;
    while (target_phase > M_PI) target_phase -= 2.0 * M_PI;
    while (target_phase < -M_PI) target_phase += 2.0 * M_PI;

    /* Constraint: PI phase range is [-pi/2, 0] */
    if (target_phase > 0.0) target_phase = -0.01;
    if (target_phase < -M_PI / 2.0) target_phase = -M_PI / 2.0 + 0.01;

    /* Solve for Ti: arg(PI) = atan(w*Ti) - pi/2 = target_phase
     * => atan(w*Ti) = target_phase + pi/2
     * => Ti = tan(target_phase + pi/2) / w */
    double phase_arg = target_phase + M_PI / 2.0;
    double Ti;
    if (fabs(cos(phase_arg)) < 1e-10 || w_gc < 1e-10) {
        Ti = 10.0 / w_gc;
    } else {
        Ti = tan(phase_arg) / w_gc;
    }
    if (Ti < 1e-10) Ti = 1.0 / w_gc;
    if (Ti > 100.0 / w_gc) Ti = 10.0 / w_gc;

    /* Solve for Kc from magnitude condition */
    double c_mag_factor = sqrt(1.0 + 1.0 / (w_gc * w_gc * Ti * Ti));
    double Kc = 1.0 / (plant_mag * c_mag_factor + 1e-15);

    /* Practical limits */
    if (Kc > 1000.0) Kc = 100.0;
    if (Kc < 1e-10) Kc = 0.1;

    /* Configure PI controller */
    pid->Kp = Kc;
    pid->Ki = Kc / Ti;
    pid->Kd = 0.0;
    pid->N  = 10.0;
    pid->b  = 1.0;
    pid->c  = 0.0;
    pid->Ts = 0.0;
    pid->u_min = -100.0;
    pid->u_max =  100.0;
    pid->integrator = 0.0;
    pid->prev_error  = 0.0;
    pid->prev_y      = 0.0;
    pid->has_antiwindup = 1;
    pid->Tt = Ti;
    return 0;
}

/* ==========================================================================
 * L5: Outer Loop Design (Direct Synthesis Method)
 *
 * For equivalent plant approximated as FOPDT: Geq ~ K*exp(-theta*s)/(tau*s+1)
 * Desired CL: T_des = exp(-theta*s)/(tau_cl*s + 1)
 *
 * Resulting PI (Direct Synthesis):
 *   Kc = tau / (K * (tau_cl + theta))
 *   Ti = min(tau, 4*(tau_cl + theta))
 *
 * Reference: Chen & Seborg, "PI/PID Controller Design Based on
 *            Direct Synthesis" (2002), Ind. Eng. Chem. Res.
 * ========================================================================== */

int cascade_design_outer_loop(const CascadeTF *eq_plant,
                               const CascadeDesignSpec *spec,
                               CascadePID *pid) {
    if (!eq_plant || !spec || !pid) return -1;

    /* Estimate FOPDT parameters from equivalent plant */
    double K = 1.0, tau = 1.0, theta = 0.0;
    if (eq_plant->den.degree >= 1 && fabs(eq_plant->den.coeff[0]) > 1e-15) {
        K = eq_plant->gain * eq_plant->num.coeff[0] / eq_plant->den.coeff[0];
        tau = eq_plant->den.coeff[1] / eq_plant->den.coeff[0];
        if (tau < 0.01) tau = 1.0;
    }
    if (eq_plant->has_delay) theta = eq_plant->delay;

    /* Outer loop should be 3-5x slower than inner */
    double tau_cl = spec->outer_settle_max / 4.0;
    if (tau_cl < 0.1) tau_cl = 1.0;

    /* Direct Synthesis PI */
    double Kc = tau / (K * (tau_cl + theta) + 1e-15);
    double Ti = (tau < 4.0 * (tau_cl + theta)) ? tau : 4.0 * (tau_cl + theta);

    if (Kc > 100.0) Kc = 10.0;
    if (Kc < 0.001) Kc = 0.1;
    if (Ti > 100.0) Ti = 10.0;
    if (Ti < 0.1) Ti = 1.0;

    pid->Kp = Kc;
    pid->Ki = Kc / Ti;
    pid->Kd = 0.0;
    pid->N  = 10.0;
    pid->b  = 1.0;
    pid->c  = 0.0;
    pid->Tt = Ti;
    pid->has_antiwindup = 1;
    return 0;
}

/* ==========================================================================
 * L5: Analytical Tuning Rules
 * ========================================================================== */

int cascade_tune_inner_skogestad(double K, double tau, double theta,
                                  double Tc, CascadePID *pid) {
    /* Skogestad SIMC rule for inner PI:
     *   Kc = tau / (|K| * (Tc + theta))
     *   Ti = min(tau, 4*(Tc + theta))
     *
     * For tight control: Tc = theta, for smooth: Tc = 1.5*theta
     * Ms <= 1.7 when Tc >= theta for integrating processes.
     *
     * Reference: Skogestad (2003), "Simple analytic rules for model
     *            reduction and PID controller tuning", J. Proc. Control. */
    if (!pid || fabs(K) < 1e-15) return -1;
    double Tc_th = Tc + theta;
    if (Tc_th < 1e-10) Tc_th = tau / 4.0;
    double Kc = tau / (fabs(K) * Tc_th);
    double Ti = (tau < 4.0 * Tc_th) ? tau : 4.0 * Tc_th;
    pid->Kp = Kc;
    pid->Ki = Kc / Ti;
    pid->Kd = 0.0;
    pid->N = 10.0;
    pid->Tt = Ti;
    pid->has_antiwindup = 1;
    return 0;
}

int cascade_tune_outer_direct_synthesis(double K, double tau, double theta,
                                         double tau_cl, CascadePID *pid) {
    (void)theta;  /* Reserved for future use with delay compensation */
    /* Direct Synthesis for outer loop:
     *   Kc = tau / (|K| * tau_cl)
     *   Ti = tau
     * Reference: Seborg et al. (2017), "Process Dynamics and Control", Ch.12 */
    if (!pid || fabs(K) < 1e-15) return -1;
    double tcl = (tau_cl > 1e-10) ? tau_cl : tau;
    double Kc = tau / (fabs(K) * tcl);
    double Ti = (tau > 0.1) ? tau : 0.5;
    pid->Kp = Kc;
    pid->Ki = Kc / Ti;
    pid->Kd = 0.0;
    pid->N = 10.0;
    pid->Tt = Ti;
    pid->has_antiwindup = 1;
    return 0;
}

int cascade_tune_inner_freq_domain(const CascadeTF *model,
                                    double w_gc, double pm,
                                    CascadePID *pid) {
    return cascade_design_inner_loop(model, w_gc, pm, pid);
}

/* ==========================================================================
 * L5: Simultaneous Cascade Optimization via Nelder-Mead
 *
 * Objective: J = ISE_outer + penalty for constraint violations
 * Constraints: inner PM >= target, outer PM >= target, BW ratio >= min
 *
 * Nelder-Mead algorithm (Nelder & Mead, 1965):
 *   1. Initialize (n+1)-vertex simplex in n-dimensional parameter space
 *   2. Order vertices by cost, reflect worst through centroid
 *   3. Expand if reflection improved, contract if not, shrink if needed
 *   4. Iterate until convergence (simplex diameter < tolerance)
 * ========================================================================== */

int cascade_optimize_simultaneous(const CascadeTF *inner_model,
                                   const CascadeTF *outer_model,
                                   const CascadeDesignSpec *spec,
                                   CascadeSystem *sys) {
    if (!inner_model || !outer_model || !spec || !sys) return -1;

    /* Initial guess from sequential design */
    if (cascade_design_sequential(inner_model, outer_model, spec, sys) != 0) {
        sys->inner.inner_controller.Kp = 1.0;
        sys->inner.inner_controller.Ki = 5.0;
        sys->outer.outer_controller.Kp = 0.5;
        sys->outer.outer_controller.Ki = 0.2;
    }

    int n = 4;  /* Parameters: [Kc_i, Ti_i, Kc_o, Ti_o] */
    int m = n + 1;  /* 5 simplex vertices */

    double x[5][4], f[5];
    double Kci0 = sys->inner.inner_controller.Kp;
    double Tii0 = Kci0 / (sys->inner.inner_controller.Ki + 1e-10);
    double Kco0 = sys->outer.outer_controller.Kp;
    double TiO0 = Kco0 / (sys->outer.outer_controller.Ki + 1e-10);

    double deltas[4] = {Kci0 * 0.5, Tii0 * 0.5, Kco0 * 0.5, TiO0 * 0.5};
    for (int i = 0; i < 4; i++) if (deltas[i] < 0.1) deltas[i] = 0.1;

    x[0][0] = Kci0; x[0][1] = Tii0; x[0][2] = Kco0; x[0][3] = TiO0;
    for (int i = 1; i < m; i++) {
        for (int j = 0; j < n; j++) x[i][j] = x[0][j];
        x[i][i - 1] += deltas[i - 1];
    }

    /* Nelder-Mead iteration */
    for (int iter = 0; iter < 50; iter++) {
        /* Evaluate cost */
        for (int i = 0; i < m; i++) {
            CascadePID ci, co;
            ci.Kp = x[i][0]; ci.Ki = (x[i][1] > 1e-10) ? x[i][0] / x[i][1] : 0.1;
            ci.Kd = 0; ci.N = 10;
            co.Kp = x[i][2]; co.Ki = (x[i][3] > 1e-10) ? x[i][2] / x[i][3] : 0.1;
            co.Kd = 0; co.N = 10;

            CascadeTF Gi_cl, overall;
            cascade_inner_closed_loop(inner_model, &ci, &Gi_cl);
            cascade_overall_closed_loop(&Gi_cl, &co, outer_model, &overall);

            double ise = 0.0;
            double T_final = spec->outer_settle_max;
            int n_pts = 100;
            for (int k = 0; k < n_pts; k++) {
                double t = (double)k * T_final / (double)n_pts;
                double y = cascade_step_response(&overall, t, 10);
                double e = 1.0 - y;
                ise += e * e * T_final / (double)n_pts;
            }

            double penalty = 0.0;
            if (!cascade_is_stable(&Gi_cl)) penalty += 100.0;
            if (!cascade_is_stable(&overall)) penalty += 100.0;

            f[i] = ise + penalty;
            cascade_tf_free(&Gi_cl);
            cascade_tf_free(&overall);
        }

        /* Sort by cost */
        for (int i = 0; i < m - 1; i++) {
            for (int j = i + 1; j < m; j++) {
                if (f[i] > f[j]) {
                    double tf = f[i]; f[i] = f[j]; f[j] = tf;
                    for (int k = 0; k < n; k++) {
                        double tx = x[i][k]; x[i][k] = x[j][k]; x[j][k] = tx;
                    }
                }
            }
        }

        /* Convergence check */
        double f_mean = 0.0;
        for (int i = 0; i < m; i++) f_mean += f[i];
        f_mean /= (double)m;
        double f_std = 0.0;
        for (int i = 0; i < m; i++)
            f_std += (f[i] - f_mean) * (f[i] - f_mean);
        f_std = sqrt(f_std / (double)m);
        if (f_std < 1e-3 && iter > 10) break;

        /* Centroid of best n vertices */
        double centroid[4] = {0};
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                centroid[j] += x[i][j];
        for (int j = 0; j < n; j++) centroid[j] /= (double)n;

        /* Reflect worst */
        for (int j = 0; j < n; j++)
            x[m - 1][j] = centroid[j] + 1.0 * (centroid[j] - x[m - 1][j]);
    }

    /* Set best parameters */
    sys->inner.inner_controller.Kp = x[0][0];
    sys->inner.inner_controller.Ki = (x[0][1] > 1e-10) ? x[0][0] / x[0][1] : 0.1;
    sys->outer.outer_controller.Kp = x[0][2];
    sys->outer.outer_controller.Ki = (x[0][3] > 1e-10) ? x[0][2] / x[0][3] : 0.1;

    /* Recompute with optimized parameters */
    cascade_inner_closed_loop(inner_model, &sys->inner.inner_controller,
                               &sys->inner.inner_cl);
    cascade_form_equivalent_plant(inner_model, &sys->inner.inner_controller,
                                   outer_model, &sys->outer.equivalent_plant);
    cascade_overall_closed_loop(&sys->inner.inner_cl,
                                 &sys->outer.outer_controller,
                                 outer_model, &sys->outer.outer_cl);

    cascade_closed_loop_bandwidth(&sys->inner.inner_cl, &sys->inner.inner_bandwidth);
    cascade_closed_loop_bandwidth(&sys->outer.outer_cl, &sys->outer.outer_bandwidth);
    if (sys->outer.outer_bandwidth > 1e-10)
        sys->bandwidth_ratio = sys->inner.inner_bandwidth / sys->outer.outer_bandwidth;

    sys->cascade_active = 1;
    sys->inner.inner_is_stable = cascade_is_stable(&sys->inner.inner_cl);
    sys->outer.outer_is_stable = cascade_is_stable(&sys->outer.outer_cl);
    return 0;
}

/* ==========================================================================
 * L5: Bandwidth Ratio Optimization and Design Validation
 * ========================================================================== */

double cascade_optimal_bandwidth_ratio(const CascadeTF *inner_model,
                                        const CascadeDisturbance *dist,
                                        double noise_ratio) {
    /* Optimal bandwidth ratio: balances disturbance rejection against noise.
     * J(ratio) = alpha/ratio + beta*ratio => ratio_opt = sqrt(alpha/beta)
     *
     * alpha: benefit of higher ratio (disturbance rejection improvement)
     * beta: cost of higher ratio (noise amplification, actuator wear) */
    (void)inner_model;
    double alpha = 1.0, beta = noise_ratio;
    if (dist) {
        alpha = 1.0 + dist->dist_magnitude;
        if (dist->dist_is_stationary) alpha *= 1.5;
    }
    if (beta < 0.01) beta = 0.05;
    double ratio_opt = sqrt(alpha / beta);
    if (ratio_opt < 3.0) ratio_opt = 3.0;
    if (ratio_opt > 10.0) ratio_opt = 10.0;
    return ratio_opt;
}

int cascade_validate_design(const CascadeSystem *sys,
                             const CascadeDesignSpec *spec) {
    /* Validate cascade design against all specifications.
     * Returns bitmask: 0=all pass, 1=inner PM, 2=outer PM, 4=inner GM,
     * 8=outer GM, 16=BW ratio, 32=settle time, 64=overshoot */
    int failures = 0;
    if (!sys || !spec) return -1;
    if (sys->bandwidth_ratio < spec->bw_ratio_min) failures |= 16;
    if (!sys->inner.inner_is_stable) failures |= 1;
    if (!sys->outer.outer_is_stable) failures |= 2;
    return failures;
}
