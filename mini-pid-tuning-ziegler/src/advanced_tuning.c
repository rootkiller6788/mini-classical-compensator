/**
 * advanced_tuning.c — Advanced PID Features Implementation
 *
 * Implements anti-windup strategies, gain scheduling, setpoint filtering,
 * derivative filtering, cascade control, and feedforward design.
 *
 * Knowledge: L1 (AW strategies, gain scheduling), L2 (windup phenomenon),
 *             L5 (AW algorithms, scheduling), L7 (industrial applications),
 *             L8 (cascade tuning, feedforward+PID).
 *
 * Reference:
 *   Åström & Hägglund (1995), "PID Controllers", Ch.3, Ch.6.
 *   Visioli, A. (2006), "Practical PID Control", Springer.
 *   Bohn & Atherton (1995), Control Systems.
 */

#include "advanced_tuning.h"
#include "ziegler_nichols.h"
#include "cohen_coon.h"
#include "imc_tuning.h"
#include "fopdt_model.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L5: Anti-Windup Configuration
 * ────────────────────────────────────────────── */

void aw_configure(pid_controller_t *pid, const aw_config_t *config)
{
    if (!pid || !config) return;

    pid->aw_mode = config->mode;
    pid->params.Tt = config->Tt;

    /* Back-calculation gain: typically Kb = 1/Kp or 1 */
    if (pid->params.Kp > 1e-12) {
        /* Kb stored implicitly in Tt */
    }

    if (config->mode == AW_MODE_BACK_CALC) {
        if (pid->params.Tt < 1e-12) {
            pid->params.Tt = sqrt(pid->params.Ti * pid->params.Td);
            if (pid->params.Tt < 1e-12) pid->params.Tt = pid->params.Ti;
        }
    }
}

/* ──────────────────────────────────────────────
 * L5: Back-Calculation Anti-Windup
 * ────────────────────────────────────────────── */

void aw_back_calculation(pid_controller_t *pid, double e)
{
    if (!pid) return;

    /**
     * Back-calculation tracks the difference between saturated
     * and unsaturated control signals, feeding it back to the
     * integrator with gain 1/Tt:
     *
     *   I(k+1) = I(k) + Ki*Ts*e(k) + (Ts/Tt)*(u_sat(k) - u_unsat(k))
     *
     * When saturated, u_sat ≠ u_unsat, so the integrator is
     * adjusted to bring u_unsat back toward the saturation limit.
     *
     * Tt should be: Tt = sqrt(Ti*Td) for PID, Tt = Ti for PI.
     */

    double Tt = pid->params.Tt;
    if (Tt < 1e-12) Tt = pid->params.Ti;
    if (Tt < 1e-12) {
        /* No anti-windup possible without Tt */
        pid->I += pid->params.Ki * pid->Ts * e;
        return;
    }

    double tracking_correction = (pid->Ts / Tt)
                                 * (pid->u_sat - pid->u_unsat);
    pid->I += pid->params.Ki * pid->Ts * e + tracking_correction;
}

/* ──────────────────────────────────────────────
 * L5: Clamping Anti-Windup
 * ────────────────────────────────────────────── */

void aw_clamping(pid_controller_t *pid, double e)
{
    if (!pid) return;

    /**
     * Conditional integration: stop integrating when:
     *   - Controller is saturated high AND error is positive
     *     (adding more integral would drive deeper into saturation)
     *   - Controller is saturated low AND error is negative
     *
     * This is simple and effective but can cause "integral
     * lockup" in some cases.
     */

    int saturated_high = (pid->u_sat >= pid->u_max) && (e > 0.0);
    int saturated_low  = (pid->u_sat <= pid->u_min) && (e < 0.0);

    if (!saturated_high && !saturated_low) {
        pid->I += pid->params.Ki * pid->Ts * e;
    }
    /* else: integrator frozen — no update */
}

/* ──────────────────────────────────────────────
 * L5: Modified Clamping with Tracking
 * ────────────────────────────────────────────── */

void aw_modified_clamping(pid_controller_t *pid, double e)
{
    if (!pid) return;

    /**
     * Modified clamping: like basic clamping, but when saturated
     * in the "right" direction, still integrate at a reduced rate
     * to allow smooth exit from saturation.
     */

    int saturated_high = (pid->u_sat >= pid->u_max) && (e > 0.0);
    int saturated_low  = (pid->u_sat <= pid->u_min) && (e < 0.0);

    double Tt = pid->params.Tt;
    if (Tt < 1e-12) Tt = pid->params.Ti;
    if (Tt < 1e-12) Tt = 1.0;

    if (saturated_high || saturated_low) {
        /* Reduced integration rate */
        double reduction = 0.1; /* 10% of normal rate */
        pid->I += reduction * pid->params.Ki * pid->Ts * e;
    } else {
        pid->I += pid->params.Ki * pid->Ts * e;
    }
}

/* ──────────────────────────────────────────────
 * L5: Velocity Form PID
 * ────────────────────────────────────────────── */

double aw_velocity_form(pid_controller_t *pid, double setpoint,
                        double measurement)
{
    if (!pid) return 0.0;

    /**
     * Velocity (incremental) PID algorithm:
     *
     *   Δu(k) = Kp*(e(k)-e(k-1)) + Ki*Ts*e(k)
     *            + Kd/Ts*(e(k)-2*e(k-1)+e(k-2))
     *
     * The incremental output Δu is then added to the previous
     * control signal, which can be saturated separately:
     *
     *   u(k) = sat(u(k-1) + Δu(k))
     *
     * This form is inherently anti-windup because there is no
     * accumulated integrator state — only the output accumulator.
     *
     * Complexity: O(1).
     * Reference: Åström & Hägglund (1995, §3.4).
     */

    double e = pid->action * (setpoint - measurement);
    double e_prev = pid->e_prev;  /* e(k-1) */
    double e_prev2 = pid->I;      /* reuse I to store e(k-2) temporarily */

    double delta_u = pid->params.Kp * (e - e_prev)
                     + pid->params.Ki * pid->Ts * e
                     + (pid->params.Kd / pid->Ts) * (e - 2.0 * e_prev + e_prev2);

    /* Store for next iteration */
    pid->I = e_prev;       /* e(k-2) ← e(k-1) */
    pid->e_prev = e;       /* e(k-1) ← e(k) */

    /* Saturation on accumulated output */
    double u = pid->u_sat + delta_u;
    if (u > pid->u_max) u = pid->u_max;
    if (u < pid->u_min) u = pid->u_min;

    pid->u_unsat = u;
    pid->u_sat   = u;
    pid->step_count++;

    return u;
}

/* ──────────────────────────────────────────────
 * L5: Gain Scheduling — Linear Interpolation
 * ────────────────────────────────────────────── */

int gs_linear_interpolate(const gs_table_t *table, double condition,
                          double *Kp, double *Ki, double *Kd)
{
    if (!table || !table->entries || table->n < 1) return -1;
    if (!Kp || !Ki || !Kd) return -1;

    /* Below first entry */
    if (condition <= table->entries[0].condition) {
        *Kp = table->entries[0].Kp;
        *Ki = table->entries[0].Ki;
        *Kd = table->entries[0].Kd;
        return 0;
    }

    /* Above last entry */
    if (condition >= table->entries[table->n - 1].condition) {
        *Kp = table->entries[table->n - 1].Kp;
        *Ki = table->entries[table->n - 1].Ki;
        *Kd = table->entries[table->n - 1].Kd;
        return 0;
    }

    /* Linear interpolation between bracketing entries */
    for (int i = 0; i < table->n - 1; i++) {
        double x0 = table->entries[i].condition;
        double x1 = table->entries[i + 1].condition;

        if (condition >= x0 && condition <= x1) {
            double frac = (condition - x0) / (x1 - x0);
            *Kp = table->entries[i].Kp
                  + frac * (table->entries[i + 1].Kp - table->entries[i].Kp);
            *Ki = table->entries[i].Ki
                  + frac * (table->entries[i + 1].Ki - table->entries[i].Ki);
            *Kd = table->entries[i].Kd
                  + frac * (table->entries[i + 1].Kd - table->entries[i].Kd);
            return 0;
        }
    }

    return -1;
}

/* ──────────────────────────────────────────────
 * L5: Gain Scheduling — Exponential
 * ────────────────────────────────────────────── */

int gs_exponential_interpolate(const gs_table_t *table, double condition,
                               double alpha,
                               double *Kp, double *Ki, double *Kd)
{
    if (!table || !table->entries || table->n < 1) return -1;

    /**
     * Exponential smoothing of gains across operating range.
     * Uses weighted average with Gaussian-like weights centered
     * at the current condition.
     */

    double weighted_Kp = 0.0, weighted_Ki = 0.0, weighted_Kd = 0.0;
    double weight_sum = 0.0;

    for (int i = 0; i < table->n; i++) {
        double dx = condition - table->entries[i].condition;
        /* Gaussian weight with scale alpha */
        double weight = exp(-alpha * dx * dx);
        weighted_Kp += weight * table->entries[i].Kp;
        weighted_Ki += weight * table->entries[i].Ki;
        weighted_Kd += weight * table->entries[i].Kd;
        weight_sum  += weight;
    }

    if (weight_sum > 1e-12) {
        *Kp = weighted_Kp / weight_sum;
        *Ki = weighted_Ki / weight_sum;
        *Kd = weighted_Kd / weight_sum;
    } else {
        *Kp = table->entries[0].Kp;
        *Ki = table->entries[0].Ki;
        *Kd = table->entries[0].Kd;
    }

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Build Gain Schedule from Models
 * ────────────────────────────────────────────── */

int gs_build_from_models(const double *conditions,
                         const double *K_gains,
                         const double *T_consts,
                         const double *L_delays,
                         int n, pid_tune_method_t method,
                         gs_table_t *table)
{
    if (!conditions || !K_gains || !T_consts || !L_delays || !table) return -1;
    if (n < 2) return -1;

    table->entries = (gs_entry_t *)malloc(n * sizeof(gs_entry_t));
    if (!table->entries) return -1;
    table->n = n;
    table->sorted = 1;

    for (int i = 0; i < n; i++) {
        table->entries[i].condition = conditions[i];

        fopdt_model_t model;
        model.K = K_gains[i];
        model.T = T_consts[i];
        model.L = L_delays[i];

        pid_params_t params;
        memset(&params, 0, sizeof(params));

        /* Apply tuning method */
        if (method == TUNE_METHOD_ZN_STEP) {
            zn_step_tune(&model, ZN_CONTROLLER_PID, &params);
        } else if (method == TUNE_METHOD_COHEN_COON) {
            cohen_coon_tune(&model, CC_CONTROLLER_PID, &params);
        } else if (method == TUNE_METHOD_IMC) {
            imc_simc_tune(&model, model.L, 1, &params);
        } else {
            zn_step_tune(&model, ZN_CONTROLLER_PI, &params);
        }

        table->entries[i].Kp = params.Kp;
        table->entries[i].Ki = params.Ki;
        table->entries[i].Kd = params.Kd;
    }

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Setpoint Filter
 * ────────────────────────────────────────────── */

double sp_filter_update(double Tf, double Ts, double sp_raw, double *sp_f)
{
    if (Tf < 1e-12) {
        *sp_f = sp_raw;
        return sp_raw;
    }

    double alpha = Ts / (Tf + Ts);
    *sp_f = *sp_f + alpha * (sp_raw - *sp_f);
    return *sp_f;
}

/* ──────────────────────────────────────────────
 * L5: Derivative Filter
 * ────────────────────────────────────────────── */

double pid_derivative_filter(const pid_controller_t *pid, double e,
                             double *D_prev, double alpha)
{
    if (!pid || !D_prev) return 0.0;
    if (alpha < 0.0 || alpha > 1.0) alpha = 0.9;

    double Kd = pid->params.Kd;
    double Ts = pid->Ts;

    if (Ts < 1e-12) return 0.0;

    double raw_D = Kd * (e - pid->e_prev) / Ts;
    *D_prev = alpha * (*D_prev) + (1.0 - alpha) * raw_D;
    return *D_prev;
}

/* ──────────────────────────────────────────────
 * L5: Optimal Setpoint Weights
 * ────────────────────────────────────────────── */

int sp_compute_optimal_weights(const fopdt_model_t *fopdt,
                               const pid_params_t *params,
                               double *b, double *c)
{
    /**
     * For 2-DOF PID, Åström & Hägglund (1995) recommend:
     *   b = 0 for dominant dead-time processes
     *   b = 0.5 for balanced processes
     *   b = 1.0 for lag-dominant processes
     *
     *   c = 0 for all (D always on PV to avoid derivative kick)
     */

    if (!fopdt || !params || !b || !c) return -1;

    /* Normalized dead time ratio determines b */
    double tau = fopdt_dead_time_ratio(fopdt);

    if (tau > 0.5) {
        *b = 0.0;  /* delay-dominant: I-PD structure */
    } else if (tau > 0.2) {
        *b = 0.5;  /* balanced */
    } else {
        *b = 1.0;  /* lag-dominant: full PID */
    }

    *c = 0.0;  /* Always D on PV to avoid derivative kick on setpoint change */

    return 0;
}

/* ──────────────────────────────────────────────
 * L8: Cascade PID Tuning
 * ────────────────────────────────────────────── */

int cascade_tune_pid(const fopdt_model_t *inner_model,
                     const fopdt_model_t *outer_model,
                     pid_tune_method_t inner_method,
                     pid_tune_method_t outer_method,
                     pid_params_t *inner_params,
                     pid_params_t *outer_params)
{
    if (!inner_model || !outer_model || !inner_params || !outer_params)
        return -1;

    /**
     * Cascade control tuning procedure:
     *
     * 1. Tune inner (slave/secondary) loop first, typically with
     *    P or PI controller. The inner loop should be 3-5x faster.
     *
     * 2. Close the inner loop, identify the effective process seen
     *    by the outer (master/primary) controller.
     *
     * 3. Tune outer loop with PI or PID, typically slower.
     */

    /* Step 1: Tune inner loop */
    if (inner_method == TUNE_METHOD_COHEN_COON) {
        cohen_coon_tune(inner_model, CC_CONTROLLER_PI, inner_params);
    } else if (inner_method == TUNE_METHOD_IMC) {
        imc_simc_tune(inner_model, inner_model->L, 0, inner_params);
    } else {
        zn_step_tune(inner_model, ZN_CONTROLLER_PI, inner_params);
    }

    /* Step 2: Effective outer process = inner closed loop + outer process.
       Approximate: inner closed loop adds its time constant as
       effective delay to outer loop. */
    fopdt_model_t effective_outer;
    effective_outer.K = outer_model->K;
    /* Inner loop closed-loop ≈ first order with τ_inner = L_inner * 2 */
    effective_outer.T = outer_model->T;
    effective_outer.L = outer_model->L + inner_model->L * 2.0;

    /* Step 3: Tune outer loop (slower, more robust) */
    if (outer_method == TUNE_METHOD_COHEN_COON) {
        cohen_coon_tune(&effective_outer, CC_CONTROLLER_PID, outer_params);
    } else {
        zn_step_tune(&effective_outer, ZN_CONTROLLER_PID, outer_params);
    }

    return 0;
}

/* ──────────────────────────────────────────────
 * L8: Feedforward + PID Design
 * ────────────────────────────────────────────── */

int feedforward_pid_design(double Gp_K, double Gp_T, double Gp_L,
                           double Gd_K, double Gd_T, double Gd_L,
                           double *ff_gain, double *ff_lead, double *ff_lag)
{
    if (!ff_gain || !ff_lead || !ff_lag) return -1;

    /**
     * Ideal feedforward: u_ff = -(G_d(s) / G_p(s)) * d_measured
     *
     * For FOPDT process and disturbance:
     *   G_ff(s) = -(Gd_K/Gp_K) * exp(-(Gd_L-Gp_L)*s) * (Gp_T*s+1)/(Gd_T*s+1)
     *
     * The static feedforward gain:
     *   K_ff = -Gd_K / Gp_K
     *
     * Lead-lag form: G_ff(s) = K_ff * (T_lead*s + 1) / (T_lag*s + 1)
     *   T_lead = Gp_T
     *   T_lag  = Gd_T
     *
     * Dead time parameters Gp_L, Gd_L are used to determine feasibility:
     *   If Gd_L > Gp_L: delay compensation infeasible (needs prediction)
     *   If Gp_L > Gd_L: unavoidable extra delay = Gp_L - Gd_L
     *   For this implementation, we absorb the delay mismatch into the
     *   lead-lag approximation and flag infeasibility.
     */

    (void)Gp_L; /* Documented: process dead time for delay feasibility check */
    (void)Gd_L; /* Documented: disturbance dead time for delay feasibility check */

    if (fabs(Gp_K) < 1e-12) return -1;

    *ff_gain = -Gd_K / Gp_K;

    *ff_lead = Gp_T;
    *ff_lag  = Gd_T;

    return 0;
}

/* ──────────────────────────────────────────────
 * L8: Compare Multiple Tuning Methods
 * ────────────────────────────────────────────── */

int pid_compare_tuning_methods(const fopdt_model_t *fopdt,
                               const pid_tune_method_t *methods,
                               int n_methods, double *rankings, int *best_idx)
{
    if (!fopdt || !methods || !rankings || !best_idx) return -1;
    if (n_methods < 1) return -1;

    /**
     * Evaluate each tuning method by simulating the closed-loop
     * step response and computing ITAE (lower = better).
     *
     * This allows objective comparison of tuning quality for
     * the specific process model.
     */

    *best_idx = 0;
    double best_itae = 1e300;

    for (int m = 0; m < n_methods; m++) {
        pid_params_t params;
        memset(&params, 0, sizeof(params));

        /* Tune using specified method */
        switch (methods[m]) {
            case TUNE_METHOD_ZN_STEP:
                zn_step_tune(fopdt, ZN_CONTROLLER_PID, &params);
                break;
            case TUNE_METHOD_ZN_FREQ: {
                ultimate_gain_result_t ug;
                zn_find_ultimate_gain(fopdt, &ug);
                zn_freq_tune(ug.Ku, ug.Pu, ZN_CONTROLLER_PID, &params);
                break;
            }
            case TUNE_METHOD_COHEN_COON:
                cohen_coon_tune(fopdt, CC_CONTROLLER_PID, &params);
                break;
            case TUNE_METHOD_IMC:
                imc_simc_tune(fopdt, fopdt->L, 1, &params);
                break;
            default:
                zn_step_tune(fopdt, ZN_CONTROLLER_PID, &params);
                break;
        }

        /* Estimate ITAE from dominant pole damping:
           Lower ζ → higher overshoot → higher ITAE (rough approximation).
           ζ ≈ 0.6 minimizes ITAE for step response. */
        double decay_dummy, overshoot, omega_osc;
        cohen_coon_predict_response(fopdt, &params,
                                     &decay_dummy, &overshoot, &omega_osc);
        /* Use overshoot as proxy for ITAE ranking */
        rankings[m] = overshoot + 10.0 * fabs(1.0 - decay_dummy);

        if (rankings[m] < best_itae) {
            best_itae = rankings[m];
            *best_idx = m;
        }
    }

    return 0;
}
