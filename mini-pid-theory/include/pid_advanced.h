
#ifndef PID_ADVANCED_H
#define PID_ADVANCED_H

#include "mini-pid-theory.h"

#ifdef __cplusplus
extern "C" {
#endif

/* L8: Advanced PID features - anti-windup, cascade, feedforward, etc. */

/* Cascade control structure.
 * Outer loop (primary): controls the main process variable
 * Inner loop (secondary): controls an intermediate variable for faster disturbance rejection
 *
 * Architecture:
 *   ysp_outer -> C_outer(s) -> ysp_inner -> C_inner(s) -> G_inner(s) -> y_inner -> G_outer(s) -> y_outer
 *                    ^                                                             |
 *                    |_____________ feedback outer ________________________________|
 */
typedef struct {
    pid_params_t outer_params;
    pid_params_t inner_params;
    pid_state_t  outer_state;
    pid_state_t  inner_state;
    double       inner_setpoint_max;
    double       inner_setpoint_min;
} pid_cascade_t;

int pid_cascade_init(pid_cascade_t *cascade,
                     const pid_params_t *outer, const pid_params_t *inner);
double pid_cascade_compute(pid_cascade_t *cascade,
                           double ysp_outer, double y_outer,
                           double y_inner, double dt);

/* Feedforward control structure.
 * Adds a feedforward term based on measured disturbance d:
 *   u_total = u_pid + K_ff * d
 * where K_ff = -G_d(0)/G(0) for perfect steady-state disturbance rejection
 */
typedef struct {
    double K_ff;         /* Static feedforward gain */
    double T_lead;       /* Lead time constant for dynamic FF */
    double T_lag;        /* Lag time constant for dynamic FF */
    double ff_saturation; /* Feedforward signal saturation limit */
} pid_feedforward_t;

void pid_feedforward_init(pid_feedforward_t *ff, double K_ff,
                          double T_lead, double T_lag);
double pid_feedforward_compute(const pid_feedforward_t *ff,
                               double disturbance, double dt,
                               double *ff_state);

/* Gain scheduling structure.
 * PID gains Kc, Ti, Td are scheduled as functions of a scheduling variable z
 * (e.g., production rate, operating point, measured load).
 *
 * Linear interpolation between breakpoints in a table.
 */
typedef struct {
    double *z_breakpoints;  /* Scheduling variable breakpoints, length N */
    double *Kc_table;       /* Gain at each breakpoint */
    double *Ti_table;       /* Integral time at each breakpoint */
    double *Td_table;       /* Derivative time at each breakpoint */
    int     N;              /* Number of breakpoints */
    int     own_memory;     /* 1 if struct owns the arrays */
} pid_gain_schedule_t;

int pid_gain_schedule_init(pid_gain_schedule_t *gs,
                           const double *z_bp, const double *Kc_bp,
                           const double *Ti_bp, const double *Td_bp,
                           int N);
void pid_gain_schedule_free(pid_gain_schedule_t *gs);
int pid_gain_schedule_lookup(const pid_gain_schedule_t *gs, double z,
                             double *Kc, double *Ti, double *Td);

/* Setpoint filter / ramp structure.
 * Smooths step changes in setpoint to reduce overshoot:
 *   ysp_filtered(s) = ysp(s) / (T_f*s + 1)
 * Or rate-limited ramp: |d(ysp)/dt| <= ramp_rate
 */
typedef struct {
    double T_f;          /* First-order filter time constant [s] */
    double ramp_rate;    /* Maximum rate of change [units/s] (0 = no limit) */
    double ysp_current;  /* Current filtered setpoint value (state) */
} pid_setpoint_filter_t;

void pid_sp_filter_init(pid_setpoint_filter_t *filt, double T_f, double ramp_rate);
double pid_sp_filter_compute(pid_setpoint_filter_t *filt,
                             double ysp_target, double dt);

/* Explicit anti-windup with configurable tracking time constant.
 *
 * The back-calculation method feeds the saturation error back to the integrator:
 *   I[k] = I[k-1] + (Kc*Ts/Ti)*e[k] + (Ts/Tt)*(u_sat[k-1] - u_ideal[k-1])
 *
 * Tt (tracking time constant) determines how fast the integrator resets.
 * Common choices: Tt = sqrt(Ti*Td)  or  Tt = Ti/2  or  Tt = Ti
 *
 * Reference: Astrom & Hagglund (1995), Section 3.5
 */
double pid_antiwindup_tracking_time(const pid_params_t *params);

/* Setpoint weighting optimization.
 * For given plant model and PID tuning, find optimal b and c weights
 * that minimize IAE subject to Ms <= Ms_max constraint.
 *
 * Uses golden-section search over b in [0,1] and c in [0,1].
 *
 * Returns optimal b, c values.
 */
int pid_optimize_setpoint_weights(double K, double tau, double theta,
                                  const pid_params_t *base_params,
                                  double Ms_max,
                                  double *b_opt, double *c_opt,
                                  double *iae_opt);

/* Bumpless parameter change.
 * When PID gains are changed online (e.g., via gain scheduling),
 * the integrator and derivative states must be adjusted to avoid bumps.
 *
 * Adjusts I[k] so that with new parameters, the output remains continuous:
 *   I_new = u_prev - Kc_new * (b_new*ysp - y) - D_new
 */
void pid_bumpless_param_change(const pid_params_t *old_params,
                               const pid_params_t *new_params,
                               pid_state_t *state,
                               double ysp, double y);

/* PID with Smith Predictor for dominant deadtime processes.
 *
 * G(s) = G0(s) * exp(-theta*s)
 * The Smith Predictor uses an internal model G0(s) to predict the
 * delay-free output y'(t), allowing the PID to control G0(s) directly.
 *
 * This is especially beneficial when theta/tau > 1.
 *
 * Reference: Smith (1957), Astrom & Hagglund Section 6.4
 */
typedef struct {
    pid_params_t pid;
    pid_state_t  pid_state;
    double       K_model;      /* Internal model gain */
    double       tau_model;    /* Internal model time constant */
    double       theta_model;  /* Internal model deadtime */
    /* Delay line for storing past model inputs */
    double      *delay_buffer;
    int          delay_length;
    int          delay_index;
} pid_smith_predictor_t;

int pid_smith_init(pid_smith_predictor_t *smith,
                   const pid_params_t *params,
                   double K_model, double tau_model, double theta_model,
                   double Ts);
void pid_smith_free(pid_smith_predictor_t *smith);
double pid_smith_compute(pid_smith_predictor_t *smith,
                         double ysp, double y, double dt);

#ifdef __cplusplus
}
#endif
#endif /* PID_ADVANCED_H */
