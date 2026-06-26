#ifndef PID_TUNING_H
#define PID_TUNING_H

/**
 * pid_tuning.h — Core PID Controller Structures and Tuning Framework
 *
 * Knowledge coverage:
 *   L1 Definitions: PID controller forms (parallel, series, ideal),
 *                   tuning parameters Kp/Ki/Kd, Ti/Td, N-filter
 *   L2 Concepts:    feedback error, proportional/derivative kick,
 *                   reset windup, setpoint tracking vs regulation
 *   L3 Math:        PID transfer function, Laplace/Fourier forms
 *
 * Reference: Åström & Hägglund, "PID Controllers: Theory, Design, and
 *            Tuning" (ISA, 1995), Chapter 2-3.
 */

#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: PID Controller Forms & Parameter Structures
 * ────────────────────────────────────────────── */

/** PID controller form enumeration.
 *
 *  IDEAL (standard/ISA):  C(s) = Kp*(1 + 1/(Ti*s) + Td*s)
 *  PARALLEL (independent): C(s) = Kp + Ki/s + Kd*s
 *  SERIES (interacting):   C(s) = Kp'*(1 + 1/(Ti'*s))*(1 + Td'*s)
 *  DERIVATIVE_ON_PV:       D-term acts on measurement only (not error)
 *  PROPORTIONAL_ON_PV:     P-term acts on measurement only (I-PD structure)
 */
typedef enum {
    PID_FORM_IDEAL = 0,
    PID_FORM_PARALLEL,
    PID_FORM_SERIES,
    PID_FORM_DERIVATIVE_ON_PV,
    PID_FORM_PROPORTIONAL_ON_PV,
    PID_FORM_2DOF
} pid_form_t;

/** Tuning method enumeration.
 *
 *  Each method represents a distinct tuning philosophy documented
 *  in the control theory literature.
 */
typedef enum {
    TUNE_METHOD_ZN_STEP = 0,       /* Ziegler-Nichols step response */
    TUNE_METHOD_ZN_FREQ,           /* Ziegler-Nichols frequency response */
    TUNE_METHOD_COHEN_COON,        /* Cohen-Coon (1953) */
    TUNE_METHOD_TYREUS_LUYBEN,     /* Tyreus-Luyben (1992) */
    TUNE_METHOD_IMC,               /* Internal Model Control based */
    TUNE_METHOD_AMIGO,             /* Åström & Hägglund AMIGO */
    TUNE_METHOD_LAMBDA,            /* Lambda tuning (Dahlin/Smith) */
    TUNE_METHOD_CHIEN_HRONES,      /* Chien-Hrones-Reswick */
    TUNE_METHOD_GM_PM,             /* Gain/Phase margin specification */
    TUNE_METHOD_RELAY_AUTOTUNE,    /* Åström-Hägglund relay auto-tuner */
    TUNE_METHOD_MANUAL             /* User-specified parameters */
} pid_tune_method_t;

/** PID Controller action direction */
typedef enum {
    PID_DIRECT_ACTING = 1,   /* Output increases with error (normal) */
    PID_REVERSE_ACTING = -1  /* Output decreases with error (e.g., cooling) */
} pid_action_t;

/** Back-calculation anti-windup mode */
typedef enum {
    AW_MODE_NONE = 0,         /* No anti-windup protection */
    AW_MODE_BACK_CALC,        /* Back-calculation (feedback from saturated output) */
    AW_MODE_CLAMPING,         /* Conditional integration (clamping) */
    AW_MODE_VELOCITY,         /* Velocity form (inherently anti-windup) */
    AW_MODE_MODIFIED_CLAMPING /* Tracking time constant based */
} antiwindup_mode_t;

/**
 * L1: Core PID parameter structure.
 *
 * Stores gains in multiple representations for convenience.
 * Internal consistency: Kp always defined. For IDEAL form:
 *   Ki = Kp / Ti  (if Ti > 0)
 *   Kd = Kp * Td  (if Td > 0)
 */
typedef struct {
    /* Primary gains — always valid */
    double Kp;           /* Proportional gain [output/error] */
    double Ki;           /* Integral gain [output/(error*time)] */
    double Kd;           /* Derivative gain [output*time/error] */

    /* Time-constant representation (IDEAL form) */
    double Ti;           /* Integral time constant [time], Ti > 0 means active */
    double Td;           /* Derivative time constant [time], Td > 0 means active */

    /* Derivative filter: N determines filter bandwidth */
    double N;            /* Derivative filter factor (typical 2-20) */

    /* Setpoint weighting (2-DOF structure) */
    double b;            /* Proportional setpoint weight (0 ≤ b ≤ 1) */
    double c;            /* Derivative setpoint weight (0 ≤ c ≤ 1) */

    /* Anti-windup tracking time constant */
    double Tt;           /* Tracking time constant for back-calculation */
} pid_params_t;

/**
 * L1: PID Controller state (run-time data).
 *
 * Stores integrator state, last error for derivative, saturation flags.
 */
typedef struct {
    /* Controller parameters */
    pid_params_t  params;

    /* Configuration */
    pid_form_t    form;
    pid_action_t  action;
    antiwindup_mode_t aw_mode;

    /* State variables */
    double I;          /* Integral accumulator */
    double e_prev;     /* Previous error (for derivative) */
    double y_prev;     /* Previous measurement (for D-on-PV) */
    double D_filt;     /* Filtered derivative term state */
    double u_sat;      /* Saturated control signal */
    double u_unsat;    /* Unsaturated (ideal) control signal */

    /* Limits */
    double u_min;      /* Output lower saturation limit */
    double u_max;      /* Output upper saturation limit */

    /* Sample time */
    double Ts;         /* Sampling period [time] */
    int    step_count; /* Number of controller updates */
} pid_controller_t;

/**
 * L1: Process model for tuning: First-Order Plus Dead Time (FOPDT).
 *
 *   G(s) = K * exp(-L*s) / (T*s + 1)
 *
 * This is the most fundamental process model in PID tuning.
 *   K  = static gain
 *   T  = time constant (lag)
 *   L  = apparent dead time (delay)
 */
typedef struct {
    double K;   /* Static gain [(output)/(input)] */
    double T;   /* Time constant [time] */
    double L;   /* Apparent dead time [time] */
} fopdt_model_t;

/**
 * L1: Second-Order Plus Dead Time (SOPDT) model.
 *
 *   G(s) = K * exp(-L*s) / ( (T1*s+1)*(T2*s+1) )
 *   or
 *   G(s) = K * exp(-L*s) * ωₙ² / (s² + 2ζωₙs + ωₙ²)
 */
typedef struct {
    double K;    /* Static gain */
    double T1;   /* First time constant */
    double T2;   /* Second time constant */
    double L;    /* Apparent dead time */
    /* Alternative oscilatory representation */
    double zeta; /* Damping ratio */
    double wn;   /* Natural frequency [rad/time] */
    int    use_osc; /* 1 = use zeta/wn form, 0 = T1/T2 form */
} sopdt_model_t;

/**
 * L1: Integrating Plus Dead Time (IPDT) model.
 *
 *   G(s) = Kv * exp(-L*s) / s
 *
 * Common in level control and slow process applications.
 */
typedef struct {
    double Kv;  /* Velocity gain [1/time] */
    double L;   /* Dead time [time] */
} ipdt_model_t;

/**
 * L1: Step-response experiment data for Ziegler-Nichols identification.
 */
typedef struct {
    double *time;     /* Time vector [length N] */
    double *input;    /* Input step signal */
    double *output;   /* Process output response */
    int     N;        /* Number of data points */
    double  Ts;       /* Sample period */
    double  step_mag; /* Magnitude of the input step */
} step_response_data_t;

/** Ultimate gain/period experiment result (ZN frequency method) */
typedef struct {
    double Ku;  /* Ultimate gain (gain at stability limit) */
    double Pu;  /* Ultimate period [time] */
    int    converged;  /* 1 if experiment identified valid Ku/Pu */
} ultimate_gain_result_t;

/**
 * L2: PID Performance Metrics.
 *
 * Used to evaluate closed-loop response quality.
 */
typedef struct {
    double IAE;    /* Integral of Absolute Error */
    double ISE;    /* Integral of Squared Error */
    double ITAE;   /* Integral of Time-weighted Absolute Error */
    double ITSE;   /* Integral of Time-weighted Squared Error */
    double overshoot_pct; /* Percent overshoot */
    double settling_time; /* Settling time (to within 2%) */
    double rise_time;     /* Rise time (10% to 90%) */
    double steady_state_error; /* Final steady-state error */
    int    oscillation_count; /* Number of sign changes */
} pid_perf_t;

/**
 * L3: Transfer function representation of the PID controller.
 *
 * Parallel canonical form:
 *   C(s) = Kp + Ki/s + Kd*s
 *
 * With derivative filter (realizable):
 *   C(s) = Kp + Ki/s + (Kd*s) / (1 + s*Td/N)
 *
 * Frequency response:
 *   C(jω) = Kp + Ki/(jω) + jω*Kd / (1 + jω*Td/N)
 */
typedef struct {
    double Kp, Ki, Kd, Td, N;
} pid_tf_params_t;

/* ──────────────────────────────────────────────
 * L5: PID Controller Runtime Functions
 * ────────────────────────────────────────────── */

/**
 * Initialize a PID controller structure with sensible defaults.
 *
 * Complexity: O(1). Theorem reference: Åström & Hägglund (1995) Ch.3.
 *
 * @param pid     Controller to initialize.
 * @param Kp      Proportional gain.
 * @param Ki      Integral gain (parallel form).
 * @param Kd      Derivative gain (parallel form).
 * @param Ts      Sampling period.
 */
void pid_init(pid_controller_t *pid, double Kp, double Ki, double Kd, double Ts);

/**
 * Convert ideal-form parameters to parallel-form gains.
 *
 * Given Kp, Ti, Td (IDEAL form), compute Kp, Ki, Kd (PARALLEL form).
 *   Ki = Kp / Ti
 *   Kd = Kp * Td
 *
 * @param Kp  Proportional gain (in/out for consistency)
 * @param Ti  Integral time constant
 * @param Td  Derivative time constant
 * @param Ki  [out] Integral gain
 * @param Kd_out [out] Derivative gain
 */
void pid_ideal_to_parallel(double Kp, double Ti, double Td,
                           double *Ki, double *Kd_out);

/**
 * Convert parallel-form gains back to ideal-form time constants.
 *
 * @param Kp   Proportional gain
 * @param Ki   Integral gain
 * @param Kd   Derivative gain
 * @param Ti   [out] Integral time constant
 * @param Td   [out] Derivative time constant
 */
void pid_parallel_to_ideal(double Kp, double Ki, double Kd,
                           double *Ti, double *Td);

/**
 * Convert series (interacting) form to parallel gains.
 *
 * Series:  C(s) = Kp' * (1 + 1/(Ti'*s)) * (1 + Td'*s)
 *          = Kp' * (1 + Td'/Ti' + 1/(Ti'*s) + Td'*s)
 *
 * @param Kp_s  Series proportional gain
 * @param Ti_s  Series integral time
 * @param Td_s  Series derivative time
 * @param Kp_p  [out] Parallel Kp
 * @param Ki_p  [out] Parallel Ki
 * @param Kd_p  [out] Parallel Kd
 */
void pid_series_to_parallel(double Kp_s, double Ti_s, double Td_s,
                            double *Kp_p, double *Ki_p, double *Kd_p);

/**
 * Compute PID controller output for one time step.
 *
 * Implements the velocity-form PID algorithm with anti-windup
 * and derivative filtering. This is the core computational method.
 *
 * Complexity: O(1) per call.
 * Reference: Åström & Hägglund (1995), Algorithm 3.1.
 *
 * @param pid        Controller state (updated in-place).
 * @param setpoint   Desired value.
 * @param measurement Current measured value.
 * @return           Control signal (saturated to [u_min, u_max]).
 */
double pid_update(pid_controller_t *pid, double setpoint, double measurement);

/**
 * Reset the PID controller's integrator and state.
 *
 * Complexity: O(1). Used during bumpless transfer or mode changes.
 */
void pid_reset(pid_controller_t *pid);

/**
 * Compute PID frequency response at given frequency ω.
 *
 * Evaluates the PID transfer function C(jω) = P(jω) + I(jω) + D(jω).
 * The derivative filter is included: D(s) = Kd*s / (1 + s*Td/N).
 *
 * @param params  PID parameters.
 * @param omega   Angular frequency [rad/time].
 * @param mag     [out] Magnitude |C(jω)|.
 * @param phase   [out] Phase ∠C(jω) [radians].
 */
void pid_freq_response(const pid_params_t *params, double omega,
                       double *mag, double *phase);

/**
 * Compute the PID controller transfer function coefficients as rational
 * polynomial ratio: C(s) = (a2*s² + a1*s + a0) / (s*(1 + s*Td/N))
 *
 * @param p     PID parameters.
 * @param a2    [out] s² coefficient of numerator.
 * @param a1    [out] s¹ coefficient of numerator.
 * @param a0    [out] s⁰ coefficient of numerator.
 */
void pid_transfer_function_coeffs(const pid_params_t *p,
                                  double *a2, double *a1, double *a0);

/**
 * Evaluate PID closed-loop performance from step-response data.
 *
 * Computes IAE, ISE, ITAE, ITSE, overshoot, settling time.
 * Complexity: O(N) where N = number of data points.
 *
 * @param t      Time vector [length N].
 * @param y      Output response [length N].
 * @param sp     Setpoint (constant step).
 * @param N      Number of points.
 * @param perf   [out] Computed performance metrics.
 */
void pid_eval_performance(const double *t, const double *y,
                          double sp, int N, pid_perf_t *perf);

/**
 * Set controller output limits and anti-windup configuration.
 *
 * @param pid   Controller.
 * @param mode  Anti-windup strategy.
 * @param u_min Lower saturation bound.
 * @param u_max Upper saturation bound.
 */
void pid_set_limits(pid_controller_t *pid, antiwindup_mode_t mode,
                    double u_min, double u_max);

/**
 * Configure setpoint weighting for 2-DOF structure.
 *
 * The controller equation becomes:
 *   u = Kp*(b*sp - y) + Ki*∫(sp - y)dt + Kd*(c*dsp/dt - dy/dt)
 *
 * b = 1, c = 1 → standard PID (error feedback)
 * b = 0, c = 0 → I-PD (integral only on error)
 *
 * @param pid Controller.
 * @param b   Proportional setpoint weight (default 1.0).
 * @param c   Derivative setpoint weight (default 0.0 for D-on-PV).
 */
void pid_set_setpoint_weights(pid_controller_t *pid, double b, double c);

/* ──────────────────────────────────────────────
 * L5: PID Tuning Rule Application Functions
 * ────────────────────────────────────────────── */

/**
 * Apply a pre-computed tuning rule result to a PID controller.
 *
 * Fills the pid_params_t structure and configures form appropriately.
 *
 * @param pid       Controller to configure.
 * @param p         Filled parameter structure.
 * @param form      Controller form (IDEAL or PARALLEL).
 * @param Ts        Sampling period.
 */
void pid_apply_tuning(pid_controller_t *pid, const pid_params_t *p,
                      pid_form_t form, double Ts);

/**
 * Convert between PID controller forms analytically.
 *
 * Handles conversion between IDEAL, PARALLEL, and SERIES forms
 * with correct algebraic relationships.
 *
 * @param src        Source parameters.
 * @param src_form   Source controller form.
 * @param dst        [out] Destination parameters.
 * @param dst_form   Destination controller form.
 */
void pid_convert_form(const pid_params_t *src, pid_form_t src_form,
                      pid_params_t *dst, pid_form_t dst_form);

/**
 * Derive derivative filter constant N from desired filter bandwidth.
 *
 * The derivative term becomes: D(s) = Kd*s / (1 + s*Td/N)
 * The cutoff frequency of this filter is ωc = N / Td.
 * N = 2 to 20 is typical; N = ∞ means no filter.
 *
 * @param Td      Derivative time constant.
 * @param cutoff  Desired filter cutoff frequency [rad/time].
 * @return        Filter constant N.
 */
double pid_derive_filter_N(double Td, double cutoff);

/**
 * Compute the minimum recommended sampling period for a PID controller
 * based on the process dynamics (Nyquist-Shannon implied by ω_bw).
 *
 * Rule of thumb: Ts ≤ (0.1 to 0.3) * (L + T), but exact bound depends
 * on closed-loop bandwidth.
 *
 * @param omega_bw  Closed-loop bandwidth estimate [rad/time].
 * @return          Recommended Ts.
 */
double pid_recommend_sampling(double omega_bw);

#ifdef __cplusplus
}
#endif

#endif /* PID_TUNING_H */
