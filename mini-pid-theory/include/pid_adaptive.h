
#ifndef PID_ADAPTIVE_H
#define PID_ADAPTIVE_H

#include "mini-pid-theory.h"
#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* L8: Adaptive PID - online adaptation of PID parameters. */

/* Auto-tuning state machine states.
 *
 * AUTO_TUNE_IDLE:    Normal PID operation, not tuning
 * AUTO_TUNE_RELAY:   Relay feedback experiment in progress
 * AUTO_TUNE_ANALYZE: Analyzing oscillation data to extract Ku, Tu
 * AUTO_TUNE_COMPUTE: Computing PID parameters from Ku, Tu
 * AUTO_TUNE_DONE:    Tuning complete, switching to new parameters
 * AUTO_TUNE_FAILED:  Tuning failed - process too noisy or no oscillation
 */
typedef enum {
    AUTO_TUNE_IDLE    = 0,
    AUTO_TUNE_RELAY   = 1,
    AUTO_TUNE_ANALYZE = 2,
    AUTO_TUNE_COMPUTE = 3,
    AUTO_TUNE_DONE    = 4,
    AUTO_TUNE_FAILED  = 5
} autotune_state_t;

/* Auto-tuning data structure.
 * Manages the Astrom-Hagglund relay feedback auto-tuning experiment.
 */
typedef struct {
    autotune_state_t state;
    double relay_amplitude;
    double hysteresis;
    int    max_cycles;
    double Ts;
    /* Measured data */
    double y_prev;
    double y_peak_high;
    double y_peak_low;
    double t_last_crossing;
    double oscillation_amplitude;
    double oscillation_period;
    int    half_cycle_count;
    int    samples_this_cycle;
    double prev_relay_output;
    /* Results */
    double Ku;
    double Tu;
} autotune_t;

void autotune_init(autotune_t *at, double relay_amplitude, double hysteresis,
                   int max_cycles, double Ts);
int autotune_update(autotune_t *at, double y, double ysp);
double autotune_get_control(const autotune_t *at);
int autotune_get_results(const autotune_t *at, double *Ku, double *Tu);

/* Adaptive PID using MIT rule (Model Reference Adaptive Control).
 *
 * Reference model: Gm(s) = omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
 * Adaptation law (MIT rule): dKc/dt = -gamma * e * ym * sign(K)
 * where e = y - ym (error between process and reference model)
 *
 * Reference: Astrom & Wittenmark "Adaptive Control" (1995), Section 5.3
 */
typedef struct {
    double omega_n;    /* Reference model natural frequency [rad/s] */
    double zeta;       /* Reference model damping ratio */
    double gamma;      /* Adaptation gain */
    double Kc_min;     /* Minimum allowed gain */
    double Kc_max;     /* Maximum allowed gain */
    /* Reference model state */
    double ym;         /* Model output */
    double dym;        /* Model derivative */
    double y_prev;     /* Process output previous step */
    double Kc_current; /* Current adapted gain */
} pid_mit_adaptive_t;

void pid_mit_init(pid_mit_adaptive_t *adapt, double omega_n, double zeta,
                  double gamma, double Kc0, double Kc_min, double Kc_max);
double pid_mit_update(pid_mit_adaptive_t *adapt, double ysp, double y, double dt);

/* Extremum-seeking auto-tuner.
 * Perturbs PID parameters with sinusoidal dither and uses gradient
 * estimation to optimize a cost function (e.g., IAE over a moving window).
 *
 * Based on: Krstic & Wang (2000), "Stability of extremum seeking feedback"
 */
typedef struct {
    double a;            /* Dither amplitude */
    double omega;        /* Dither frequency [rad/s] */
    double gamma;        /* Gradient descent gain */
    double Kc_base;      /* Current base gain */
    double Ti_base;      /* Current base integral time */
    double cost_prev;    /* Previous cost function value */
    double phase;        /* Current dither phase */
    double highpass_prev;/* Highpass filter state */
    double cost_integral;/* Moving cost integral */
    int    window_size;  /* Cost evaluation window [samples] */
    int    sample_idx;   /* Current sample in window */
    double *error_window;/* Ring buffer of recent errors */
} pid_extremum_seeking_t;

int pid_es_init(pid_extremum_seeking_t *es, double a, double omega,
                double gamma, double Kc0, double Ti0, int window_size);
void pid_es_free(pid_extremum_seeking_t *es);
int pid_es_update(pid_extremum_seeking_t *es, double error, double dt,
                  double *Kc_out, double *Ti_out);

/* Iterative Feedback Tuning (IFT).
 *
 * Model-free tuning method that uses closed-loop experiments to estimate
 * the gradient of a cost criterion with respect to PID parameters.
 *
 * IFT performs three experiments per iteration:
 *   1. Normal operation: collect r1(t), y1(t)
 *   2. Special experiment: r2(t) = r(t) - y1(t), collect y2(t)
 *   3. Gradient estimation from correlation of y1, y2 signals
 *
 * Reference: Hjalmarsson et al. (1998), "Iterative feedback tuning"
 */
typedef struct {
    double Kc, Ti, Td;       /* Current parameters */
    double dJ_dKc, dJ_dTi, dJ_dTd; /* Gradient estimates */
    double step_size;        /* Gradient descent step size */
    int    iteration;
    int    experiment_phase; /* 0=normal, 1=gradient Kc, 2=gradient Ti, 3=gradient Td */
    /* Buffers for experiment data */
    double *r1_buf, *y1_buf, *r2_buf, *y2_buf;
    int    buf_size;
    int    buf_idx;
} pid_ift_tuner_t;

int pid_ift_init(pid_ift_tuner_t *ift, double Kc0, double Ti0, double Td0,
                 double step_size, int buf_size);
void pid_ift_free(pid_ift_tuner_t *ift);
int pid_ift_iteration(pid_ift_tuner_t *ift);

/* Self-Tuning Regulator (STR) for PID.
 *
 * Uses recursive least squares (RLS) to estimate a FOPDT model online,
 * then re-tunes the PID using SIMC rules.
 *
 * This combines system identification with PID tuning in a closed loop.
 */
typedef struct {
    /* RLS parameter estimates: y[k] + a1*y[k-1] = b0*u[k-d] + b1*u[k-d-1] */
    double a1, b0, b1;
    double P[3][3];   /* Covariance matrix */
    double lambda;     /* Forgetting factor (0.95-0.995) */
    double theta_est[3]; /* [a1, b0, b1] */
    double phi[3];     /* Regression vector */
    int    delay;      /* Estimated delay in samples */
    /* PID re-tuning parameters */
    double K_model, tau_model, theta_model;
    double tau_c;      /* Desired closed-loop time constant */
    int    param_update_interval; /* Re-tune every N samples */
    int    sample_count;
} pid_str_tuner_t;

void pid_str_init(pid_str_tuner_t *str, double lambda, int delay,
                  double tau_c, int update_interval);
int pid_str_update(pid_str_tuner_t *str, double y, double u_prev, double Ts,
                   pid_params_t *params);

#ifdef __cplusplus
}
#endif
#endif /* PID_ADAPTIVE_H */
