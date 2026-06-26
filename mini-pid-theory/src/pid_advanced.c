
/**
 * pid_advanced.c - Advanced PID Features
 *
 * L8: Advanced Topics
 *   - Cascade control
 *   - Feedforward control
 *   - Gain scheduling
 *   - Setpoint filtering/ramping
 *   - Explicit anti-windup with tracking time
 *   - Setpoint weighting optimization
 *   - Bumpless parameter change
 *   - Smith Predictor for deadtime-dominant processes
 */

#include "mini-pid-theory.h"
#include "pid_advanced.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --------------------------------------------------------------------------
 * Cascade Control
 *
 * Outer loop (primary) controls the main PV.
 * Inner loop (secondary) controls an intermediate PV for fast disturbance rejection.
 *
 * The inner loop must be at least 3-5x faster than the outer loop.
 * Tuning procedure:
 *   1. Tune inner loop first (with outer in manual)
 *   2. Tune outer loop with inner in auto
 * -------------------------------------------------------------------------- */

int pid_cascade_init(pid_cascade_t *cascade,
                     const pid_params_t *outer, const pid_params_t *inner)
{
    if (!cascade || !outer || !inner) return -1;

    memcpy(&cascade->outer_params, outer, sizeof(pid_params_t));
    memcpy(&cascade->inner_params, inner, sizeof(pid_params_t));

    pid_init(outer, &cascade->outer_state);
    pid_init(inner, &cascade->inner_state);

    cascade->inner_setpoint_max = inner->umax;
    cascade->inner_setpoint_min = inner->umin;

    return 0;
}

double pid_cascade_compute(pid_cascade_t *cascade,
                           double ysp_outer, double y_outer,
                           double y_inner, double dt)
{
    if (!cascade) return 0.0;

    /* Outer loop: compute inner setpoint */
    double ysp_inner = pid_compute(&cascade->outer_params, &cascade->outer_state,
                                   ysp_outer, y_outer, dt);

    /* Clamp inner setpoint to valid range */
    if (ysp_inner > cascade->inner_setpoint_max)
        ysp_inner = cascade->inner_setpoint_max;
    if (ysp_inner < cascade->inner_setpoint_min)
        ysp_inner = cascade->inner_setpoint_min;

    /* Inner loop: compute control signal */
    double u = pid_compute(&cascade->inner_params, &cascade->inner_state,
                           ysp_inner, y_inner, dt);

    return u;
}

/* --------------------------------------------------------------------------
 * Feedforward Control
 *
 * Adds a term based on measured disturbance d:
 *   u_total = u_pid + u_ff
 *
 * Static feedforward: u_ff = K_ff * d
 * Dynamic feedforward: u_ff(s) = K_ff * (T_lead*s + 1)/(T_lag*s + 1) * d(s)
 *
 * The ideal feedforward gain for disturbance rejection:
 *   K_ff = -G_d(0) / G(0)
 * where G_d(s) is the disturbance transfer function.
 * -------------------------------------------------------------------------- */

void pid_feedforward_init(pid_feedforward_t *ff, double K_ff,
                          double T_lead, double T_lag)
{
    if (!ff) return;
    ff->K_ff = K_ff;
    ff->T_lead = T_lead;
    ff->T_lag = T_lag;
    ff->ff_saturation = 1e10;
}

double pid_feedforward_compute(const pid_feedforward_t *ff,
                               double disturbance, double dt,
                               double *ff_state)
{
    if (!ff || !ff_state) return 0.0;

    /* Static feedforward */
    double u_ff = ff->K_ff * disturbance;

    /* Dynamic lead-lag filter if time constants are specified */
    if (ff->T_lead > 1e-10 || ff->T_lag > 1e-10) {
        /* Discrete first-order lead-lag:
         * u_ff[k] = alpha*u_ff[k-1] + beta*d[k] + gamma*d[k-1]
         * where alpha = T_lag/(T_lag+Ts), beta = K_ff*(T_lead+Ts)/(T_lag+Ts),
         *       gamma = -K_ff*T_lead/(T_lag+Ts)
         */
        double Ts = (dt > 0.0) ? dt : 0.01;
        double alpha = ff->T_lag / (ff->T_lag + Ts);
        double beta = ff->K_ff * (ff->T_lead + Ts) / (ff->T_lag + Ts);
        double gamma = -ff->K_ff * ff->T_lead / (ff->T_lag + Ts);

        double prev_state = ff_state[0];
        double prev_dist = ff_state[1];

        u_ff = alpha * prev_state + beta * disturbance + gamma * prev_dist;

        ff_state[0] = u_ff;
        ff_state[1] = disturbance;
    }

    /* Saturation */
    if (u_ff > ff->ff_saturation) u_ff = ff->ff_saturation;
    if (u_ff < -ff->ff_saturation) u_ff = -ff->ff_saturation;

    return u_ff;
}

/* --------------------------------------------------------------------------
 * Gain Scheduling
 *
 * PID gains are interpolated from a lookup table based on a scheduling
 * variable z (e.g., production rate, valve position, measured load).
 * -------------------------------------------------------------------------- */

int pid_gain_schedule_init(pid_gain_schedule_t *gs,
                           const double *z_bp, const double *Kc_bp,
                           const double *Ti_bp, const double *Td_bp,
                           int N)
{
    if (!gs || !z_bp || !Kc_bp || !Ti_bp || !Td_bp || N < 2) return -1;

    gs->N = N;
    gs->z_breakpoints = (double*)malloc(N * sizeof(double));
    gs->Kc_table = (double*)malloc(N * sizeof(double));
    gs->Ti_table = (double*)malloc(N * sizeof(double));
    gs->Td_table = (double*)malloc(N * sizeof(double));

    if (!gs->z_breakpoints || !gs->Kc_table || !gs->Ti_table || !gs->Td_table) {
        pid_gain_schedule_free(gs);
        return -1;
    }

    memcpy(gs->z_breakpoints, z_bp, N * sizeof(double));
    memcpy(gs->Kc_table, Kc_bp, N * sizeof(double));
    memcpy(gs->Ti_table, Ti_bp, N * sizeof(double));
    memcpy(gs->Td_table, Td_bp, N * sizeof(double));
    gs->own_memory = 1;

    return 0;
}

void pid_gain_schedule_free(pid_gain_schedule_t *gs)
{
    if (!gs) return;
    if (gs->own_memory) {
        free(gs->z_breakpoints);
        free(gs->Kc_table);
        free(gs->Ti_table);
        free(gs->Td_table);
    }
    memset(gs, 0, sizeof(*gs));
}

int pid_gain_schedule_lookup(const pid_gain_schedule_t *gs, double z,
                             double *Kc, double *Ti, double *Td)
{
    if (!gs || !Kc || !Ti || !Td || gs->N < 2) return -1;

    const double *zbp = gs->z_breakpoints;

    /* Below first breakpoint: use first value */
    if (z <= zbp[0]) {
        *Kc = gs->Kc_table[0];
        *Ti = gs->Ti_table[0];
        *Td = gs->Td_table[0];
        return 0;
    }

    /* Above last breakpoint: use last value */
    if (z >= zbp[gs->N - 1]) {
        int last = gs->N - 1;
        *Kc = gs->Kc_table[last];
        *Ti = gs->Ti_table[last];
        *Td = gs->Td_table[last];
        return 0;
    }

    /* Binary search for interval */
    int lo = 0, hi = gs->N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (z < zbp[mid]) hi = mid;
        else lo = mid;
    }

    /* Linear interpolation */
    double frac = (z - zbp[lo]) / (zbp[hi] - zbp[lo]);
    *Kc = gs->Kc_table[lo] + frac * (gs->Kc_table[hi] - gs->Kc_table[lo]);
    *Ti = gs->Ti_table[lo] + frac * (gs->Ti_table[hi] - gs->Ti_table[lo]);
    *Td = gs->Td_table[lo] + frac * (gs->Td_table[hi] - gs->Td_table[lo]);

    return 0;
}

/* --------------------------------------------------------------------------
 * Setpoint Filter / Rate Limiter
 *
 * Smooths setpoint changes to reduce overshoot and derivative kick.
 *
 * First-order filter: ysp_f(s) = ysp(s) / (T_f*s + 1)
 * Rate limiter: |d(ysp)/dt| <= ramp_rate
 * -------------------------------------------------------------------------- */

void pid_sp_filter_init(pid_setpoint_filter_t *filt, double T_f, double ramp_rate)
{
    if (!filt) return;
    filt->T_f = T_f;
    filt->ramp_rate = ramp_rate;
    filt->ysp_current = 0.0;
}

double pid_sp_filter_compute(pid_setpoint_filter_t *filt,
                             double ysp_target, double dt)
{
    if (!filt) return ysp_target;

    double Ts = (dt > 0.0) ? dt : 0.01;
    double ysp_out = ysp_target;

    /* First-order filter */
    if (filt->T_f > 1e-10) {
        double alpha = Ts / (filt->T_f + Ts);
        ysp_out = filt->ysp_current + alpha * (ysp_target - filt->ysp_current);
    }

    /* Rate limiter */
    if (filt->ramp_rate > 1e-10) {
        double max_step = filt->ramp_rate * Ts;
        double diff = ysp_out - filt->ysp_current;
        if (diff > max_step) {
            ysp_out = filt->ysp_current + max_step;
        } else if (diff < -max_step) {
            ysp_out = filt->ysp_current - max_step;
        }
    }

    filt->ysp_current = ysp_out;
    return ysp_out;
}

/* --------------------------------------------------------------------------
 * Anti-windup Tracking Time Constant
 *
 * Recommended Tt values for back-calculation anti-windup:
 *   Tt = sqrt(Ti * Td)   (Astrom & Hagglund default)
 *   Tt = Ti/2            (more rapid reset)
 *   Tt = Ti              (slower, smoother reset)
 *
 * Returns the recommended Tt.
 * -------------------------------------------------------------------------- */

double pid_antiwindup_tracking_time(const pid_params_t *params)
{
    if (!params) return 1.0;
    double Tt = sqrt(params->Ti * params->Td + 1e-10);
    if (Tt < 1e-10) Tt = params->Ti * 0.5;
    if (Tt < 1e-10) Tt = 1.0;
    return Tt;
}

/* --------------------------------------------------------------------------
 * Setpoint Weighting Optimization
 *
 * Finds optimal b and c weights to minimize IAE subject to Ms <= Ms_max.
 *
 * b: proportional setpoint weight [0, 1]
 *    b = 0: P-term on measurement only (no proportional kick)
 *    b = 1: full P-term on setpoint error (classic)
 *
 * c: derivative setpoint weight [0, 1]
 *    c = 0: D-term on measurement only (recommended)
 *    c = 1: full D-term on setpoint error (causes kick)
 *
 * Generally, reducing b reduces overshoot at the cost of slower response.
 * Optimal b typically 0.3-0.7 for most processes.
 *
 * This function uses grid search over b and c.
 * -------------------------------------------------------------------------- */

int pid_optimize_setpoint_weights(double K, double tau, double theta,
                                  const pid_params_t *base_params,
                                  double Ms_max,
                                  double *b_opt, double *c_opt,
                                  double *iae_opt)
{
    if (!base_params || !b_opt || !c_opt || !iae_opt) return -1;

    double best_iae = 1e308;
    double best_b = 1.0, best_c = 0.0;
    int found = 0;

    /* Grid search: b in 0.1 steps, c in 0.25 steps */
    for (int ib = 0; ib <= 10; ib++) {
        double b = ib * 0.1;
        for (int ic = 0; ic <= 4; ic++) {
            double c = ic * 0.25;

            pid_params_t test_params = *base_params;
            test_params.b = b;
            test_params.c = c;
            test_params.form = PID_FORM_2DOF;

            /* Check robustness constraint */
            double Ms, w_peak;
            pid_max_sensitivity(K, tau, theta, &test_params, &Ms, &w_peak);

            if (Ms <= Ms_max) {
                /* Evaluate performance */
                pid_performance_t perf;
                pid_evaluate_fopdt(K, tau, theta, &test_params,
                                   1.0, 5.0*(tau+theta), 500, &perf);

                if (perf.iae < best_iae) {
                    best_iae = perf.iae;
                    best_b = b;
                    best_c = c;
                    found = 1;
                }
            }
        }
    }

    if (found) {
        *b_opt = best_b;
        *c_opt = best_c;
        *iae_opt = best_iae;
        return 0;
    }
    return -1;
}

/* --------------------------------------------------------------------------
 * Bumpless Parameter Change
 *
 * When PID gains change online (gain scheduling, adaptive tuning),
 * adjust integrator to avoid output discontinuity.
 *
 * Old output: u_old = Kc_old * (b_old*ysp - y) + I_old + D_old
 * New output: u_new = Kc_new * (b_new*ysp - y) + I_new + D_new
 *
 * We want u_new = u_old, so:
 *   I_new = u_old - Kc_new * (b_new*ysp - y) - D_new
 *
 * Note: D_new ~ 0 at steady state (derivative of constant ~ 0)
 * -------------------------------------------------------------------------- */

void pid_bumpless_param_change(const pid_params_t *old_params,
                               const pid_params_t *new_params,
                               pid_state_t *state,
                               double ysp, double y)
{
    if (!old_params || !new_params || !state) return;

    /* Compute what the new P and D contributions would be */
    double P_new = new_params->Kc * (new_params->b * ysp - y);
    double D_new = state->prev_deriv; /* Approximate, derivative state preserved */

    /* Back-calculate integrator to maintain output continuity */
    double u_prev = state->prev_output;
    state->integral = u_prev - P_new - D_new;

    /* The first pid_compute call will fine-tune this */
}

/* --------------------------------------------------------------------------
 * Smith Predictor for Deadtime-Dominant Processes
 *
 * When theta/tau > 1, standard PID performs poorly because the deadtime
 * limits the achievable bandwidth. The Smith Predictor uses an internal
 * model to predict what the delay-free output would be.
 *
 * Architecture:
 *   y'(t) = model output without delay (G0(s))
 *   Feedback signal to PID = y'(t) + (y(t) - y'(t-theta))
 *
 * If model matches plant perfectly, PID controls G0(s) directly,
 * effectively removing the deadtime from the feedback path.
 *
 * Reference: Smith, O.J.M., "A Controller to Overcome Dead Time", ISA, 1959
 * -------------------------------------------------------------------------- */

int pid_smith_init(pid_smith_predictor_t *smith,
                   const pid_params_t *params,
                   double K_model, double tau_model, double theta_model,
                   double Ts)
{
    if (!smith || !params || Ts <= 0.0) return -1;

    memcpy(&smith->pid, params, sizeof(pid_params_t));
    pid_init(params, &smith->pid_state);

    smith->K_model = K_model;
    smith->tau_model = tau_model;
    smith->theta_model = theta_model;

    /* Allocate delay buffer for Smith predictor */
    smith->delay_length = (int)ceil(theta_model / Ts) + 1;
    smith->delay_buffer = (double*)calloc(smith->delay_length, sizeof(double));
    smith->delay_index = 0;

    if (!smith->delay_buffer) return -1;
    return 0;
}

void pid_smith_free(pid_smith_predictor_t *smith)
{
    if (!smith) return;
    free(smith->delay_buffer);
    smith->delay_buffer = NULL;
    smith->delay_length = 0;
}

double pid_smith_compute(pid_smith_predictor_t *smith,
                         double ysp, double y, double dt)
{
    if (!smith) return 0.0;

    double Ts = (dt > 0.0) ? dt : smith->pid.Ts;

    /* Internal model without delay:
     * Simulate G0(s) = K_model/(tau_model*s+1) */
    /* We need the current model input, which is the PID output (stored in state) */
    double u_model = smith->pid_state.prev_output;

    /* Delay-free model output */
    static double y0_model = 0.0;
    double alpha = exp(-Ts / smith->tau_model);
    y0_model = alpha * y0_model + smith->K_model * (1.0 - alpha) * u_model;

    /* Delayed model output (from delay buffer) */
    double y0_delayed = smith->delay_buffer[smith->delay_index];

    /* Store current delay-free output in buffer */
    smith->delay_buffer[smith->delay_index] = y0_model;
    smith->delay_index = (smith->delay_index + 1) % smith->delay_length;

    /* Smith predictor feedback signal:
     * y_fb = y0_model + (y - y0_delayed)
     * This is the prediction of what y would be without deadtime */
    double y_fb = y0_model + (y - y0_delayed);

    /* Compute PID using the predicted feedback */
    double u = pid_compute(&smith->pid, &smith->pid_state, ysp, y_fb, Ts);

    return u;
}
