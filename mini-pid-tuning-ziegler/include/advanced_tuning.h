#ifndef ADVANCED_TUNING_H
#define ADVANCED_TUNING_H

/**
 * advanced_tuning.h — Advanced PID Features and Extensions
 *
 * Knowledge coverage:
 *   L1: Anti-windup strategies, gain scheduling, setpoint filtering
 *   L2: Integrator windup phenomenon and bumpless transfer concept
 *   L4: Circle criterion and describing function for anti-windup stability
 *   L5: Advanced PID algorithms: velocity form, cascade, override
 *   L7: Industrial applications: HVAC, chemical reactors, servo drives
 *   L8: Adaptive gain scheduling, Fuzzy-PID, fractional-order PID
 *
 * Reference:
 *   Åström, K.J. & Hägglund, T. (1995) Ch.6 "Automatic Tuning and Adaptation".
 *   Visioli, A. (2006) "Practical PID Control", Springer.
 *   Bohn, C. & Atherton, D.P. (1995)
 *   "An Analysis Package Comparing PID Anti-Windup Strategies", Control Systems.
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: Anti-Windup Structures
 * ────────────────────────────────────────────── */

/**
 * Anti-windup tracking configuration.
 */
typedef struct {
    antiwindup_mode_t mode;
    double Tt;          /* Tracking time constant (back-calculation) */
    double Kb;          /* Back-calculation gain (default = 1/Kp) */
    int    use_conditional; /* 1 = conditional integration enabled */
    double clamp_high;  /* Upper clamp threshold (% of saturation) */
    double clamp_low;   /* Lower clamp threshold */
} aw_config_t;

/**
 * Bumpless transfer configuration.
 *
 * When switching between manual/auto or between different controllers,
 * bumpless transfer prevents abrupt changes in control signal.
 */
typedef struct {
    int    source_mode;  /* 0 = tracking, 1 = balancing, 2 = reset */
    double tracking_time; /* Time constant for tracking mode */
    int    balanced;      /* Internal: 1 when balanced */
} bumpless_config_t;

/* ──────────────────────────────────────────────
 * L5: Anti-Windup Implementations
 * ────────────────────────────────────────────── */

/**
 * Configure anti-windup for a PID controller.
 *
 * Sets up the anti-windup mechanism based on the selected strategy.
 * The back-calculation method modifies the integrator by:
 *   ΔI = (Kp/Ti)*e + (1/Tt)*(u_sat - u_unsat)
 *
 * @param pid    Controller.
 * @param config Anti-windup configuration.
 */
void aw_configure(pid_controller_t *pid, const aw_config_t *config);

/**
 * Back-calculation anti-windup: the standard industrial solution.
 *
 * When the control signal saturates, the integrator is "back-calculated"
 * to prevent windup. The integral update becomes:
 *
 *   I(k+1) = I(k) + Ki*Ts*e(k) + (Ts/Tt)*(u_sat(k) - u_unsat(k))
 *
 * where Tt is the tracking time constant (typically Tt = sqrt(Ti*Td)
 * for PID, or Tt = Ti for PI).
 *
 * Complexity: O(1) per call.
 * Reference: Fertik & Ross (1967), Åström & Hägglund (1995, §3.5).
 *
 * @param pid Controller state.
 * @param e   Current error.
 */
void aw_back_calculation(pid_controller_t *pid, double e);

/**
 * Clamping (conditional integration) anti-windup.
 *
 * The integrator is only updated when:
 *   - The control signal is not saturated, OR
 *   - The error has the opposite sign to the saturation direction
 *     (i.e., the integral action would help bring the output out of
 *      saturation rather than deeper into it).
 *
 * @param pid Controller state.
 * @param e   Current error.
 */
void aw_clamping(pid_controller_t *pid, double e);

/**
 * Modified clamping with tracking time constant.
 */
void aw_modified_clamping(pid_controller_t *pid, double e);

/**
 * Velocity-form PID: inherently anti-windup (no integral accumulator).
 *
 * Instead of computing absolute control signal, compute its increment:
 *   Δu(k) = Kp*(e(k)-e(k-1)) + Ki*Ts*e(k)
 *            + (Kd/Ts)*(e(k)-2*e(k-1)+e(k-2))
 *
 * The absolute control is then saturating accumulator:
 *   u(k) = sat(u(k-1) + Δu(k))
 *
 * This form has no explicit integrator to wind up.
 *
 * @param pid    Controller state (e_prev stores e(k-1); I stores e(k-2)).
 * @param setpoint Current setpoint.
 * @param measurement Current measurement.
 * @return        Control increment Δu.
 */
double aw_velocity_form(pid_controller_t *pid, double setpoint,
                        double measurement);

/* ──────────────────────────────────────────────
 * L1: Gain Scheduling
 * ────────────────────────────────────────────── */

/**
 * Gain schedule entry: (operating condition, PID gains).
 */
typedef struct {
    double condition;  /* Operating condition variable value */
    double Kp;         /* Proportional gain at this condition */
    double Ki;         /* Integral gain */
    double Kd;         /* Derivative gain */
} gs_entry_t;

/**
 * Gain schedule table.
 */
typedef struct {
    gs_entry_t *entries; /* Array of entries [length n] */
    int         n;       /* Number of entries */
    int         sorted;  /* 1 = sorted by condition ascending */
} gs_table_t;

/**
 * Linear interpolation gain scheduling.
 *
 * Given a sorted gain schedule table, interpolates PID gains
 * at the current operating condition.
 *
 * @param table        Gain schedule table (must be sorted).
 * @param condition    Current operating condition.
 * @param Kp, Ki, Kd   [out] Interpolated gains.
 * @return             0 on success, -1 if out of table range.
 */
int gs_linear_interpolate(const gs_table_t *table, double condition,
                          double *Kp, double *Ki, double *Kd);

/**
 * Exponential interpolation gain scheduling (smoother transitions).
 *
 * Uses exponential weighting for smoother gain transitions near
 * breakpoints compared to linear interpolation.
 *
 * @param table        Gain schedule table.
 * @param condition    Current operating condition.
 * @param alpha        Smoothing factor (0 = nearest, 1 = fully weighted avg).
 * @param Kp, Ki, Kd   [out] Interpolated gains.
 * @return             0 on success.
 */
int gs_exponential_interpolate(const gs_table_t *table, double condition,
                               double alpha,
                               double *Kp, double *Ki, double *Kd);

/**
 * Build a gain schedule table from process models at different
 * operating points using a specified tuning method.
 *
 * @param conditions   Array of operating points [length n].
 * @param K_gains      Process gains at each condition.
 * @param T_consts     Process time constants at each condition.
 * @param L_delays     Process dead times at each condition.
 * @param n            Number of operating points.
 * @param method       Tuning method to apply at each point.
 * @param table        [out] Constructed gain schedule.
 * @return             0 on success.
 */
int gs_build_from_models(const double *conditions,
                         const double *K_gains,
                         const double *T_consts,
                         const double *L_delays,
                         int n, pid_tune_method_t method,
                         gs_table_t *table);

/* ──────────────────────────────────────────────
 * L2: Setpoint and Measurement Filtering
 * ────────────────────────────────────────────── */

/**
 * First-order setpoint filter: reduces proportional kick on step changes.
 *
 *   F_sp(s) = 1 / (T_f * s + 1)
 *
 * The filter time constant T_f is typically:
 *   T_f = 0.1 * Td  (for aggressive)
 *   T_f = 0.5 * Td  (for conservative)
 *
 * @param Tf     Filter time constant.
 * @param Ts     Sampling period.
 * @param sp_raw Raw setpoint value.
 * @param sp_f   [in/out] Filtered setpoint state (initial = sp_raw).
 * @return       Filtered setpoint.
 */
double sp_filter_update(double Tf, double Ts, double sp_raw, double *sp_f);

/**
 * Derivative filter: first-order low-pass on the derivative term.
 *
 * The derivative term with filter becomes:
 *   D(s) = Kd * Td * s / (1 + s * Td / N) * E(s)
 *
 * Discrete implementation (backward Euler):
 *   D(k) = α * D(k-1) + (1-α) * Kd * (e(k) - e(k-1)) / Ts
 *   where α = Td / (Td + N*Ts)
 *
 * @param pid   PID controller (params must have N and Td set).
 * @param e     Current error.
 * @param D_prev [in/out] Previous filtered derivative value.
 * @param alpha Precomputed filter coefficient.
 * @return      Filtered derivative term.
 */
double pid_derivative_filter(const pid_controller_t *pid, double e,
                             double *D_prev, double alpha);

/**
 * Compute optimal setpoint weights using the Åström-Hägglund method.
 *
 * For 2-DOF PID, the setpoint weights b and c reduce overshoot
 * without affecting disturbance rejection.
 *
 * Recommended values (Åström & Hägglund, 1995):
 *   b = 0 for pure I-PD (no proportional kick)
 *   b = 1 for full PID on error
 *   c = 0 for D on PV only (derivative of measurement, not error)
 *
 * This function computes b and c that minimize ITAE for a given FOPDT.
 *
 * @param fopdt Process model.
 * @param params PID parameters (tuned for disturbance rejection).
 * @param b     [out] Optimal proportional setpoint weight.
 * @param c     [out] Optimal derivative setpoint weight.
 * @return      0 on success.
 */
int sp_compute_optimal_weights(const fopdt_model_t *fopdt,
                               const pid_params_t *params,
                               double *b, double *c);

/* ──────────────────────────────────────────────
 * L8: Cascade Control Tuning
 * ────────────────────────────────────────────── */

/**
 * Cascade PID tuning: tune inner (secondary) and outer (primary) loops
 * sequentially according to standard cascade design procedure.
 *
 * Procedure:
 *   1. Tune inner loop (fast) using any method.
 *   2. Close inner loop; identify effective process for outer loop.
 *   3. Tune outer loop (slower) using the effective process.
 *
 * The inner loop should be 3-5x faster than the outer loop.
 *
 * @param inner_model     Inner loop process model.
 * @param outer_model     Outer loop process model.
 * @param inner_method    Tuning method for inner loop.
 * @param outer_method    Tuning method for outer loop.
 * @param inner_params    [out] Inner loop PID parameters.
 * @param outer_params    [out] Outer loop PID parameters.
 * @return                0 on success.
 */
int cascade_tune_pid(const fopdt_model_t *inner_model,
                     const fopdt_model_t *outer_model,
                     pid_tune_method_t inner_method,
                     pid_tune_method_t outer_method,
                     pid_params_t *inner_params,
                     pid_params_t *outer_params);

/**
 * Feedforward + PID combined design.
 *
 * The feedforward term is:
 *   u_ff = -(G_d(s) / G_p(s)) * d_measured
 *
 * This function computes the static feedforward gain and lead-lag
 * parameters for dynamic compensation.
 *
 * @param Gp_K    Process gain.
 * @param Gp_T    Process time constant.
 * @param Gp_L    Process dead time.
 * @param Gd_K    Disturbance gain.
 * @param Gd_T    Disturbance time constant.
 * @param Gd_L    Disturbance dead time.
 * @param ff_gain [out] Feedforward gain K_ff.
 * @param ff_lead [out] Feedforward lead time constant.
 * @param ff_lag  [out] Feedforward lag time constant.
 * @return        0 on success.
 */
int feedforward_pid_design(double Gp_K, double Gp_T, double Gp_L,
                           double Gd_K, double Gd_T, double Gd_L,
                           double *ff_gain, double *ff_lead, double *ff_lag);

/**
 * Compare multiple PID tuning methods on the same process model
 * and rank them by performance criteria.
 *
 * @param fopdt    Process model.
 * @param methods  Array of tuning methods to compare.
 * @param n_methods Number of methods.
 * @param rankings [out] Array of ITAE values (lower = better).
 * @param best_idx [out] Index of best method.
 * @return         0 on success.
 */
int pid_compare_tuning_methods(const fopdt_model_t *fopdt,
                               const pid_tune_method_t *methods,
                               int n_methods, double *rankings, int *best_idx);

#ifdef __cplusplus
}
#endif

#endif /* ADVANCED_TUNING_H */
