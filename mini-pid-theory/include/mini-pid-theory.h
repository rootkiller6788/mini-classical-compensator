
#ifndef MINI_PID_THEORY_H
#define MINI_PID_THEORY_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* L1: PID form enumeration */
typedef enum {
    PID_FORM_STANDARD = 0,
    PID_FORM_PARALLEL = 1,
    PID_FORM_SERIES   = 2,
    PID_FORM_2DOF     = 3
} pid_form_t;

typedef enum {
    PID_ACTION_DIRECT  = 0,
    PID_ACTION_REVERSE = 1
} pid_action_t;

typedef enum {
    PID_DERIV_ON_ERROR       = 0,
    PID_DERIV_ON_MEASUREMENT = 1,
    PID_DERIV_ON_FILTERED    = 2
} pid_deriv_mode_t;

typedef enum {
    PID_ERROR_ABSOLUTE      = 0,
    PID_ERROR_SQUARED       = 1,
    PID_ERROR_ITAE_WEIGHTED = 2,
    PID_ERROR_SATURATED     = 3
} pid_error_type_t;

typedef enum {
    PID_WINDUP_NONE          = 0,
    PID_WINDUP_CLAMPING      = 1,
    PID_WINDUP_BACK_CALC     = 2,
    PID_WINDUP_INCREMENTAL   = 3
} pid_antiwindup_t;

/* L1: PID parameter structure - Standard (ISA) form
 * u(t) = Kc * [ e(t) + (1/Ti)*integral(e)dt + Td*de/dt ]
 * Parallel: Kp=Kc, Ki=Kc/Ti, Kd=Kc*Td
 */
typedef struct {
    double Kc, Ti, Td, Ts, N;
    double b, c;
    double umin, umax, deadband, Tt;
    pid_form_t        form;
    pid_action_t      action;
    pid_deriv_mode_t  deriv_mode;
    pid_antiwindup_t  antiwindup;
} pid_params_t;

/* L1: PID state structure */
typedef struct {
    double integral, prev_error, prev_error2;
    double prev_measurement, prev_output, prev_deriv;
    double filtered_error, saturated_output;
    int    saturation_count;
    double integrated_absolute_error;
    double integrated_squared_error;
    double integrated_time_iae;
    double output_variance, total_variation;
    uint32_t sample_count;
    uint8_t  is_initialized, is_saturated;
    uint8_t  derivative_spike, manual_mode;
    double _reserved[4];
} pid_state_t;

/* L1: 2DOF PID parameters */
typedef struct {
    double b, c, F_sp, F_ff, K_ff;
} pid_2dof_params_t;

/* L1: Performance metrics */
typedef struct {
    double iae, ise, itae, itse, tv_u;
    double overshoot, settling_time, rise_time;
    double steady_state_error;
    double gain_margin, phase_margin;
} pid_performance_t;

/* L2: Core API */
int pid_init(const pid_params_t *params, pid_state_t *state);
double pid_compute(const pid_params_t *params, pid_state_t *state,
                   double ysp, double y, double dt);
void pid_reset(pid_state_t *state);
void pid_bumpless_transfer(const pid_params_t *params, pid_state_t *state,
                           double manual_output);

/* L3: Frequency domain */
double pid_freq_magnitude(const pid_params_t *params, double omega);
double pid_freq_phase(const pid_params_t *params, double omega);
void pid_loop_transfer(double K, double tau, double theta,
                       const pid_params_t *params, double omega,
                       double *mag, double *phase);
void pid_sensitivity(double K, double tau, double theta,
                     const pid_params_t *params, double omega,
                     double *mag, double *phase);
void pid_complementary_sensitivity(double K, double tau, double theta,
                                   const pid_params_t *params, double omega,
                                   double *mag, double *phase);
int pid_max_sensitivity(double K, double tau, double theta,
                        const pid_params_t *params,
                        double *Ms, double *w_peak);

/* L4: Stability */
int pid_stability_routh(double K, double tau, const pid_params_t *params);
int pid_stability_margins(double K, double tau, double theta,
                          const pid_params_t *params,
                          double *gm, double *pm, double *w180, double *w0db);

/* L5: Computation */
int pid_closed_loop_poles(double K, double tau, double theta,
                          const pid_params_t *params,
                          double *poles, int *npoles);
int pid_evaluate_fopdt(double K, double tau, double theta,
                       const pid_params_t *params,
                       double ysp, double tfinal, int nsteps,
                       pid_performance_t *perf);
int pid_convert_form(const pid_params_t *src_params, pid_form_t src_form,
                     pid_params_t *dst_params, pid_form_t dst_form);
int pid_relay_feedback(double K, double tau, double theta, double d,
                       double *Ku, double *Tu);

#ifdef __cplusplus
}
#endif
#endif /* MINI_PID_THEORY_H */
