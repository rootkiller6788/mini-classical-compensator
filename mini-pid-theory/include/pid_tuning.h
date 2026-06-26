
#ifndef PID_TUNING_H
#define PID_TUNING_H

#include "mini-pid-theory.h"

#ifdef __cplusplus
extern "C" {
#endif

/* L5: Tuning method enumeration.
 *
 * ZIEGLER_NICHOLS_STEP: Open-loop step response method (1942)
 *   - Fit FOPDT model to process reaction curve
 *   - Tune based on a = K*theta/tau (normalized deadtime)
 *
 * ZIEGLER_NICHOLS_ULTIMATE: Closed-loop ultimate sensitivity method
 *   - Increase Kc until sustained oscillation at Ku, Tu
 *   - Kc = 0.6*Ku, Ti = 0.5*Tu, Td = 0.125*Tu
 *
 * COHEN_COON: Open-loop FOPDT-based tuning (1953)
 *   - Designed for 1/4 decay ratio
 *   - Uses K, tau, theta directly from reaction curve
 *
 * AMIGO: Approximate M-constrained Integral Gain Optimization
 *   - Astrom & Hagglund (2004) robust tuning rules
 *   - Maximizes integral gain subject to Ms <= Ms_target
 *
 * IMC: Internal Model Control based tuning
 *   - Lambda tuning for FOPDT: Kc = tau/(K*(lambda+theta))
 *   - Ti = tau, Td = 0 (PI) or Td = tau*theta/(2*tau+theta) (PID)
 *
 * CHIEN_HRONES_RESNICK: Setpoint vs disturbance response tuning
 *   - Two sets of rules: one for fast SP response, one for fast LD rejection
 *
 * TYREUS_LUYBEN: Conservative tuning for integrator plants
 *   - Kc = 0.45*Ku, Ti = 2.2*Tu, Td = Tu/6.3
 *
 * SKOGESTAD_IMC: SIMC tuning rules (2003)
 *   - Simple, analytically derived, handles wide range of processes
 */
typedef enum {
    TUNE_ZN_STEP       = 0,
    TUNE_ZN_ULTIMATE   = 1,
    TUNE_COHEN_COON    = 2,
    TUNE_AMIGO         = 3,
    TUNE_IMC_LAMBDA    = 4,
    TUNE_CHIEN_HRONES  = 5,
    TUNE_TYREUS_LUYBEN = 6,
    TUNE_SIMC           = 7,
    TUNE_CUSTOM         = 8
} pid_tuning_method_t;

/* Tuning specification structure */
typedef struct {
    pid_tuning_method_t method;
    double K;           /* Plant static gain */
    double tau;         /* Plant dominant time constant [s] */
    double theta;       /* Plant apparent deadtime [s] */
    double Ku;          /* Ultimate gain (for ZN-ultimate, TL) */
    double Tu;          /* Ultimate period [s] (for ZN-ultimate, TL) */
    double lambda;      /* IMC closed-loop time constant */
    double Ms_target;   /* Target maximum sensitivity (AMIGO, default 1.4) */
    double decay_ratio; /* Target decay ratio (C-C, default 0.25) */
    int    controller_type; /* 1=PI only, 2=PID full */
} pid_tuning_spec_t;

/* L5: Main tuning function - compute PID parameters from plant model */
int pid_tune(const pid_tuning_spec_t *spec, pid_params_t *params);

/* L5: Ziegler-Nichols step response method.
 * Computes Kc, Ti, Td from FOPDT model: G(s) = K*exp(-theta*s)/(tau*s+1)
 *
 * PI:  Kc = 0.9/(a) * (1 + 0.33*theta/tau),  Ti = 3.33*theta/(1 + 0.33*theta/tau)
 * PID: Kc = 1.2/a,                            Ti = 2.0*theta,  Td = 0.5*theta
 * where a = K*theta/tau
 *
 * Reference: Ziegler & Nichols (1942), Ogata Section 10.4
 */
int pid_tune_zn_step(double K, double tau, double theta,
                     int controller_type, pid_params_t *params);

/* L5: Ziegler-Nichols ultimate sensitivity method.
 * PI:  Kc = 0.45*Ku,  Ti = 0.83*Tu
 * PID: Kc = 0.6*Ku,   Ti = 0.5*Tu,  Td = 0.125*Tu
 */
int pid_tune_zn_ultimate(double Ku, double Tu,
                         int controller_type, pid_params_t *params);

/* L5: Cohen-Coon open-loop tuning.
 * PI:  Kc = (1/K)*(tau/theta)*(0.9 + theta/(12*tau)),
 *      Ti = theta*(30 + 3*theta/tau)/(9 + 20*theta/tau)
 * PID: Kc = (1/K)*(tau/theta)*(4/3 + theta/(4*tau)),
 *      Ti = theta*(32 + 6*theta/tau)/(13 + 8*theta/tau),
 *      Td = theta*4/(11 + 2*theta/tau)
 *
 * Reference: Cohen & Coon (1953)
 */
int pid_tune_cohen_coon(double K, double tau, double theta,
                        int controller_type, pid_params_t *params);

/* L5: AMIGO tuning rules (Astrom & Hagglund, 2004).
 * Robust PID tuning that maximizes integral gain subject to Ms constraint.
 *
 * PI:  Kc = (0.15/K)*(1 + (0.35*tau/(tau+0.35*theta))),
 *      Ti = 0.35*theta + 13*theta*tau^2/(tau^2+12*theta*tau+7*theta^2)
 * PID: Kc = (1/K)*(0.2 + 0.45*tau/theta),
 *      Ti = 0.4*theta + 0.8*tau,
 *      Td = 0.5*theta*tau/(0.3*theta + tau)
 *
 * Ms_target = 1.4 (default, conservative)
 */
int pid_tune_amigo(double K, double tau, double theta, double Ms_target,
                   int controller_type, pid_params_t *params);

/* L5: IMC-based lambda tuning.
 * PI:  Kc = tau/(K*(lambda+theta)),  Ti = tau
 * PID: Kc = (2*tau+theta)/(2*K*(lambda+theta)),
 *      Ti = tau + theta/2,  Td = tau*theta/(2*tau+theta)
 *
 * lambda >= 0.1*tau typically; larger lambda = more robust
 * Reference: Rivera, Morari & Skogestad (1986)
 */
int pid_tune_imc_lambda(double K, double tau, double theta, double lambda,
                        int controller_type, pid_params_t *params);

/* L5: SIMC tuning (Skogestad IMC, 2003).
 * PI:  Kc = tau/(K*(tau_c + theta)),  Ti = min(tau, 4*(tau_c+theta))
 * PID: same as PI with Td = theta/2 for theta-dominant processes
 * tau_c = desired closed-loop time constant
 *
 * Recommended: tau_c = theta for tight control, tau_c = 1.5*theta for robust
 */
int pid_tune_simc(double K, double tau, double theta, double tau_c,
                  int controller_type, pid_params_t *params);

/* L5: Tyreus-Luyben conservative tuning for integrating processes.
 * PI:  Kc = 0.45*Ku,  Ti = 2.2*Tu
 * PID: Kc = 0.45*Ku,  Ti = 2.2*Tu,  Td = Tu/6.3
 *
 * More conservative than ZN; designed for chemical process stability
 * Reference: Tyreus & Luyben (1992)
 */
int pid_tune_tyreus_luyben(double Ku, double Tu,
                           int controller_type, pid_params_t *params);

/* L5: Chien-Hrones-Reswick tuning for setpoint tracking or disturbance rejection.
 * mode = 0: fastest response without overshoot (0% OS)
 * mode = 1: fastest response with 20% overshoot
 * tuning_for = 0: setpoint tracking, 1: disturbance rejection
 *
 * Reference: Chien, Hrones & Reswick (1952)
 */
int pid_tune_chien_hrones(double K, double tau, double theta,
                          int mode, int tuning_for,
                          int controller_type, pid_params_t *params);

/* L5: Compute normalized deadtime (controllability index).
 * theta_norm = theta / (theta + tau) in [0, 1]
 * 0 = delay-free (easy), 1 = pure deadtime (hard)
 */
double pid_normalized_deadtime(double tau, double theta);

/* L5: Determine if PID (vs PI) is recommended based on normalized deadtime.
 * Returns 2 for PID, 1 for PI
 * Rule: theta/(theta+tau) > 0.2 => PID recommended
 */
int pid_recommend_controller_type(double tau, double theta);

/* L5: Tuning robustness analysis.
 * Given tuned parameters and plant model, estimate:
 *   dKc_max  = maximum gain increase before instability (% of Kc)
 *   dtheta_max = maximum deadtime increase before instability (% of theta)
 * These give the "fragility" of the tuning
 */
void pid_tuning_robustness(const pid_params_t *params,
                           double K, double tau, double theta,
                           double *dKc_max, double *dtheta_max);

#ifdef __cplusplus
}
#endif
#endif /* PID_TUNING_H */
