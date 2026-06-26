
/**
 * pid_tuning.c - PID Tuning Methods
 *
 * Implements classical and modern PID tuning rules.
 * Each tuning method is an independent knowledge point.
 *
 * L5: Tuning Algorithms
 *   - Ziegler-Nichols step response (1942)
 *   - Ziegler-Nichols ultimate sensitivity (1942)
 *   - Cohen-Coon (1953)
 *   - AMIGO (Astrom & Hagglund, 2004)
 *   - IMC Lambda tuning (Rivera, Morari, Skogestad, 1986)
 *   - SIMC (Skogestad, 2003)
 *   - Tyreus-Luyben (1992)
 *   - Chien-Hrones-Reswick (1952)
 */

#include "mini-pid-theory.h"
#include "pid_tuning.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --------------------------------------------------------------------------
 * Normalized Deadtime (controllability index)
 *
 * theta_n = theta / (tau + theta)   in [0, 1]
 * 0   = pure lag (easy to control)
 * 0.5 = balanced
 * 1   = pure deadtime (very hard to control)
 * -------------------------------------------------------------------------- */

double pid_normalized_deadtime(double tau, double theta)
{
    if (tau < 0.0) tau = 0.0;
    if (theta < 0.0) theta = 0.0;
    double denom = tau + theta;
    if (denom < 1e-10) return 0.0;
    return theta / denom;
}

/* --------------------------------------------------------------------------
 * Controller type recommendation.
 * Returns 2 for PID, 1 for PI.
 * Rule of thumb: if normalized deadtime > 0.2, use PID (not just PI).
 * -------------------------------------------------------------------------- */

int pid_recommend_controller_type(double tau, double theta)
{
    double theta_n = pid_normalized_deadtime(tau, theta);
    return (theta_n > 0.2) ? 2 : 1;
}

/* --------------------------------------------------------------------------
 * Ziegler-Nichols Step Response Method (1942)
 *
 * Based on the process reaction curve from an open-loop step test.
 * Fit a FOPDT model: G(s) = K*exp(-theta*s)/(tau*s+1)
 * Then compute PID parameters based on a = K*theta/tau.
 *
 * P:   Kc = 1/a
 * PI:  Kc = 0.9/a,  Ti = 3.33*theta
 * PID: Kc = 1.2/a,  Ti = 2.0*theta,  Td = 0.5*theta
 *
 * This method typically gives a decay ratio of 1/4.
 * Reference: Ziegler & Nichols, Trans. ASME, 1942
 * -------------------------------------------------------------------------- */

int pid_tune_zn_step(double K, double tau, double theta,
                     int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0 || tau < 0.0 || theta < 0.0) return -1;

    if (theta < 1e-10) theta = tau * 0.1; /* Avoid division by zero, assume 10% deadtime */

    double a = K * theta / tau;
    if (a < 1e-10) return -1;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = theta / 10.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    switch (controller_type) {
        case 0: /* P only */
            params->Kc = 1.0 / a;
            params->Ti = 1e308;
            params->Td = 0.0;
            break;
        case 1: /* PI */
            params->Kc = 0.9 / a;
            params->Ti = 3.33 * theta;
            params->Td = 0.0;
            break;
        case 2: /* PID */
        default:
            params->Kc = 1.2 / a;
            params->Ti = 2.0 * theta;
            params->Td = 0.5 * theta;
            break;
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * Ziegler-Nichols Ultimate Sensitivity Method (1942)
 *
 * The closed-loop method:
 * 1. Set Ti = infinity, Td = 0
 * 2. Increase Kc until sustained oscillation at Ku (ultimate gain)
 * 3. Measure oscillation period Tu (ultimate period)
 *
 * P:   Kc = 0.5*Ku
 * PI:  Kc = 0.45*Ku, Ti = 0.83*Tu
 * PID: Kc = 0.6*Ku,  Ti = 0.5*Tu,  Td = 0.125*Tu
 *
 * This is the most widely taught ZN method in control textbooks.
 * Reference: Ziegler & Nichols, Trans. ASME, 1942
 * -------------------------------------------------------------------------- */

int pid_tune_zn_ultimate(double Ku, double Tu,
                         int controller_type, pid_params_t *params)
{
    if (!params || Ku <= 0.0 || Tu <= 0.0) return -1;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = Tu / 100.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    switch (controller_type) {
        case 0: /* P only */
            params->Kc = 0.5 * Ku;
            params->Ti = 1e308;
            params->Td = 0.0;
            break;
        case 1: /* PI */
            params->Kc = 0.45 * Ku;
            params->Ti = 0.83 * Tu;
            params->Td = 0.0;
            break;
        case 2: /* PID */
        default:
            params->Kc = 0.6 * Ku;
            params->Ti = 0.5 * Tu;
            params->Td = 0.125 * Tu;
            break;
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * Cohen-Coon Open-Loop Tuning (1953)
 *
 * Designed for 1/4 decay ratio, similar to ZN step but uses more
 * information from the FOPDT model.
 *
 * PI:
 *   Kc = (1/K)*(tau/theta)*(0.9 + theta/(12*tau))
 *   Ti = theta*(30 + 3*theta/tau)/(9 + 20*theta/tau)
 *
 * PID:
 *   Kc = (1/K)*(tau/theta)*(4/3 + theta/(4*tau))
 *   Ti = theta*(32 + 6*theta/tau)/(13 + 8*theta/tau)
 *   Td = theta*4/(11 + 2*theta/tau)
 *
 * Reference: Cohen & Coon, Trans. ASME, 1953
 * -------------------------------------------------------------------------- */

int pid_tune_cohen_coon(double K, double tau, double theta,
                        int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0 || tau <= 0.0 || theta <= 0.0) return -1;

    double ratio = theta / tau;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = theta / 10.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    if (controller_type <= 1) { /* PI */
        params->Kc = (1.0/K) * (tau/theta) * (0.9 + theta/(12.0*tau));
        params->Ti = theta * (30.0 + 3.0*ratio) / (9.0 + 20.0*ratio);
        params->Td = 0.0;
    } else { /* PID */
        params->Kc = (1.0/K) * (tau/theta) * (4.0/3.0 + theta/(4.0*tau));
        params->Ti = theta * (32.0 + 6.0*ratio) / (13.0 + 8.0*ratio);
        params->Td = theta * 4.0 / (11.0 + 2.0*ratio);
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * AMIGO Tuning Rules (Astrom & Hagglund, 2004)
 *
 * Approximate M-constrained Integral Gain Optimization.
 * Maximizes integral gain Ki = Kc/Ti subject to Ms <= Ms_max.
 *
 * PI:
 *   Kc = (0.15/K) * (1 + 0.35*tau/(tau+0.35*theta))
 *   Ti = 0.35*theta + 13*theta*tau^2/(tau^2+12*theta*tau+7*theta^2)
 *
 * PID:
 *   Kc = (1/K) * (0.2 + 0.45*tau/theta)
 *   Ti = 0.4*theta + 0.8*tau
 *   Td = 0.5*theta*tau/(0.3*theta + tau)
 *
 * These formulas are for Ms = 1.4 (default, conservative).
 * For Ms = 2.0, multiply Kc by ~1.3.
 *
 * Reference: Astrom & Hagglund, Control Engineering Practice, 2004
 * -------------------------------------------------------------------------- */

int pid_tune_amigo(double K, double tau, double theta, double Ms_target,
                   int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0) return -1;
    if (tau < 0.0) tau = 0.0;
    if (theta < 1e-10) theta = tau * 0.1;

    /* AMIGO formulas are calibrated for Ms = 1.4 */
    double Ms_factor = (Ms_target > 0.0) ? 1.4 / Ms_target : 1.0;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = theta / 10.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    if (controller_type <= 1) { /* PI */
        params->Kc = (0.15/K) * (1.0 + 0.35*tau/(tau + 0.35*theta)) * Ms_factor;
        params->Ti = 0.35*theta + 13.0*theta*tau*tau/(tau*tau + 12.0*theta*tau + 7.0*theta*theta);
        params->Td = 0.0;
    } else { /* PID */
        params->Kc = (1.0/K) * (0.2 + 0.45*tau/theta) * Ms_factor;
        params->Ti = 0.4*theta + 0.8*tau;
        params->Td = 0.5*theta*tau/(0.3*theta + tau);
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * IMC Lambda Tuning (Rivera, Morari & Skogestad, 1986)
 *
 * Internal Model Control based tuning. The IMC filter time constant lambda
 * determines the closed-loop speed.
 *
 * PI (for FOPDT):
 *   Kc = tau / (K * (lambda + theta))
 *   Ti = tau
 *
 * PID (for FOPDT):
 *   Kc = (2*tau + theta) / (2*K * (lambda + theta))
 *   Ti = tau + theta/2
 *   Td = tau*theta / (2*tau + theta)
 *
 * lambda = desired closed-loop time constant.
 *   lambda = theta:  tight control
 *   lambda = 3*theta: robust control
 *
 * Reference: Rivera, Morari & Skogestad, IEC Process Des. Dev., 1986
 * -------------------------------------------------------------------------- */

int pid_tune_imc_lambda(double K, double tau, double theta, double lambda,
                        int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0) return -1;
    if (lambda < 1e-10) lambda = theta;
    if (tau < 0.0) tau = 0.0;
    if (theta < 0.0) theta = 0.0;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = (tau + theta) / 100.0;
    if (params->Ts < 0.001) params->Ts = 0.001;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    if (controller_type <= 1) { /* PI */
        params->Kc = tau / (K * (lambda + theta));
        params->Ti = tau;
        params->Td = 0.0;
    } else { /* PID */
        params->Kc = (2.0*tau + theta) / (2.0 * K * (lambda + theta));
        params->Ti = tau + theta/2.0;
        params->Td = tau * theta / (2.0*tau + theta + 1e-10);
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * SIMC Tuning (Skogestad IMC, 2003)
 *
 * Simple IMC tuning rules that are analytically derived but easy to use.
 * Widely adopted in process industry.
 *
 * PI:
 *   Kc = tau / (K * (tau_c + theta))
 *   Ti = min(tau, 4*(tau_c + theta))
 *
 * PID (for tau > 8*theta, cascade form recommended instead):
 *   Kc same as PI
 *   Ti = min(tau, 4*(tau_c + theta))
 *   Td = theta/2
 *
 * tau_c = desired closed-loop time constant.
 *   tau_c = theta:           tight control
 *   tau_c = 1.5*theta:       robust control
 *
 * Key insight: Ti should be exactly tau for "perfect" IMC, but bounded
 * by 4*(tau_c+theta) for robustness when tau_c is small.
 *
 * Reference: Skogestad, J. Process Control, 2003
 * -------------------------------------------------------------------------- */

int pid_tune_simc(double K, double tau, double theta, double tau_c,
                  int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0) return -1;
    if (tau_c < 1e-10) tau_c = theta;
    if (tau < 0.0) tau = 0.0;
    if (theta < 0.0) theta = 0.0;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = (tau + theta) / 100.0;
    if (params->Ts < 0.001) params->Ts = 0.001;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    params->Kc = tau / (K * (tau_c + theta));
    params->Ti = (tau < 4.0*(tau_c + theta)) ? tau : 4.0*(tau_c + theta);

    if (controller_type >= 2 && theta > 1e-10) {
        params->Td = theta / 2.0;
    } else {
        params->Td = 0.0;
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * Tyreus-Luyben Conservative Tuning (1992)
 *
 * More conservative than ZN ultimate. Designed for chemical processes
 * where stability is paramount.
 *
 * PI:  Kc = 0.45*Ku,  Ti = 2.2*Tu
 * PID: Kc = 0.45*Ku,  Ti = 2.2*Tu,  Td = Tu/6.3
 *
 * Note: Same Kc as ZN PI but much larger Ti.
 * The larger Ti reduces integral action, making the loop more stable.
 *
 * Reference: Tyreus & Luyben, IEC Research, 1992
 * -------------------------------------------------------------------------- */

int pid_tune_tyreus_luyben(double Ku, double Tu,
                           int controller_type, pid_params_t *params)
{
    if (!params || Ku <= 0.0 || Tu <= 0.0) return -1;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = Tu / 100.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    params->Kc = 0.45 * Ku;
    params->Ti = 2.2 * Tu;

    if (controller_type >= 2) {
        params->Td = Tu / 6.3;
    } else {
        params->Td = 0.0;
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * Chien-Hrones-Reswick Tuning (1952)
 *
 * Two sets of rules: one optimized for setpoint tracking, one for
 * load disturbance rejection. Two modes per set: 0% overshoot vs 20% overshoot.
 *
 * Setpoint tracking, 0% overshoot:
 *   P:  Kc = 0.3/a
 *   PI: Kc = 0.35/a, Ti = 1.2*tau
 *   PID: Kc = 0.6/a, Ti = 1.0*tau, Td = 0.5*theta
 *
 * Setpoint tracking, 20% overshoot:
 *   P:  Kc = 0.7/a
 *   PI: Kc = 0.6/a, Ti = 1.0*tau
 *   PID: Kc = 0.95/a, Ti = 1.4*tau, Td = 0.47*theta
 *
 * Disturbance rejection, 0% overshoot:
 *   P:  Kc = 0.3/a
 *   PI: Kc = 0.6/a, Ti = 4.0*theta
 *   PID: Kc = 0.95/a, Ti = 2.4*theta, Td = 0.42*theta
 *
 * Disturbance rejection, 20% overshoot:
 *   P:  Kc = 0.7/a
 *   PI: Kc = 0.7/a, Ti = 2.3*theta
 *   PID: Kc = 1.2/a, Ti = 2.0*theta, Td = 0.42*theta
 *
 * where a = K*theta/tau.
 *
 * Reference: Chien, Hrones & Reswick, Trans. ASME, 1952
 * -------------------------------------------------------------------------- */

int pid_tune_chien_hrones(double K, double tau, double theta,
                          int mode, int tuning_for,
                          int controller_type, pid_params_t *params)
{
    if (!params || K <= 0.0 || tau < 0.0 || theta < 0.0) return -1;
    if (theta < 1e-10) theta = tau * 0.1;

    double a = K * theta / tau;
    if (a < 1e-10) return -1;

    memset(params, 0, sizeof(*params));
    params->form = PID_FORM_STANDARD;
    params->action = PID_ACTION_DIRECT;
    params->deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params->antiwindup = PID_WINDUP_BACK_CALC;
    params->N = 10.0;
    params->Ts = theta / 10.0;
    params->b = 1.0;
    params->c = 0.0;
    params->umin = -1e10;
    params->umax = 1e10;

    /* mode: 0=0% OS, 1=20% OS */
    /* tuning_for: 0=setpoint, 1=disturbance */

    if (tuning_for == 0) {
        /* Setpoint tracking */
        if (mode == 0) {
            /* 0% overshoot */
            if (controller_type == 0) {
                params->Kc = 0.3/a; params->Ti = 1e308; params->Td = 0.0;
            } else if (controller_type == 1) {
                params->Kc = 0.35/a; params->Ti = 1.2*tau; params->Td = 0.0;
            } else {
                params->Kc = 0.6/a; params->Ti = tau; params->Td = 0.5*theta;
            }
        } else {
            /* 20% overshoot */
            if (controller_type == 0) {
                params->Kc = 0.7/a; params->Ti = 1e308; params->Td = 0.0;
            } else if (controller_type == 1) {
                params->Kc = 0.6/a; params->Ti = tau; params->Td = 0.0;
            } else {
                params->Kc = 0.95/a; params->Ti = 1.4*tau; params->Td = 0.47*theta;
            }
        }
    } else {
        /* Disturbance rejection */
        if (mode == 0) {
            /* 0% overshoot */
            if (controller_type == 0) {
                params->Kc = 0.3/a; params->Ti = 1e308; params->Td = 0.0;
            } else if (controller_type == 1) {
                params->Kc = 0.6/a; params->Ti = 4.0*theta; params->Td = 0.0;
            } else {
                params->Kc = 0.95/a; params->Ti = 2.4*theta; params->Td = 0.42*theta;
            }
        } else {
            /* 20% overshoot */
            if (controller_type == 0) {
                params->Kc = 0.7/a; params->Ti = 1e308; params->Td = 0.0;
            } else if (controller_type == 1) {
                params->Kc = 0.7/a; params->Ti = 2.3*theta; params->Td = 0.0;
            } else {
                params->Kc = 1.2/a; params->Ti = 2.0*theta; params->Td = 0.42*theta;
            }
        }
    }

    params->Tt = sqrt(params->Ti * params->Td + 1e-10);
    return 0;
}

/* --------------------------------------------------------------------------
 * Universal tuning dispatcher
 * -------------------------------------------------------------------------- */

int pid_tune(const pid_tuning_spec_t *spec, pid_params_t *params)
{
    if (!spec || !params) return -1;

    switch (spec->method) {
        case TUNE_ZN_STEP:
            return pid_tune_zn_step(spec->K, spec->tau, spec->theta,
                                    spec->controller_type, params);
        case TUNE_ZN_ULTIMATE:
            return pid_tune_zn_ultimate(spec->Ku, spec->Tu,
                                        spec->controller_type, params);
        case TUNE_COHEN_COON:
            return pid_tune_cohen_coon(spec->K, spec->tau, spec->theta,
                                       spec->controller_type, params);
        case TUNE_AMIGO:
            return pid_tune_amigo(spec->K, spec->tau, spec->theta,
                                  spec->Ms_target, spec->controller_type, params);
        case TUNE_IMC_LAMBDA:
            return pid_tune_imc_lambda(spec->K, spec->tau, spec->theta,
                                       spec->lambda, spec->controller_type, params);
        case TUNE_SIMC:
            return pid_tune_simc(spec->K, spec->tau, spec->theta,
                                 spec->lambda, spec->controller_type, params);
        case TUNE_TYREUS_LUYBEN:
            return pid_tune_tyreus_luyben(spec->Ku, spec->Tu,
                                          spec->controller_type, params);
        case TUNE_CHIEN_HRONES:
            return pid_tune_chien_hrones(spec->K, spec->tau, spec->theta,
                                         0, 0, spec->controller_type, params);
        default:
            return -1;
    }
}

/* --------------------------------------------------------------------------
 * Tuning Robustness Analysis
 *
 * Estimates how much the plant model can change before the closed loop
 * becomes unstable. This measures the "fragility" of a tuning.
 *
 * dKc_max: Maximum percentage increase in K*Kc product before instability
 * dtheta_max: Maximum percentage increase in deadtime before instability
 *
 * Uses Routh-Hurwitz to find the stability boundary.
 * -------------------------------------------------------------------------- */

void pid_tuning_robustness(const pid_params_t *params,
                           double K, double tau, double theta,
                           double *dKc_max, double *dtheta_max)
{
    if (!params || !dKc_max || !dtheta_max) return;

    /* Find maximum gain: binary search on K*Kc factor */
    double K_factor_lo = 1.0, K_factor_hi = 100.0;
    pid_params_t test_params = *params;

    if (pid_stability_routh(K, tau, &test_params) != 1) {
        *dKc_max = 0.0;
        *dtheta_max = 0.0;
        return;
    }

    for (int i = 0; i < 40; i++) {
        double mid = (K_factor_lo + K_factor_hi) / 2.0;
        test_params.Kc = params->Kc * mid;
        if (pid_stability_routh(K, tau, &test_params) == 1) {
            K_factor_lo = mid;
        } else {
            K_factor_hi = mid;
        }
    }
    *dKc_max = (K_factor_lo - 1.0) * 100.0;

    /* Deadtime robustness: Routh doesn't directly include theta,
     * but we can estimate via phase margin equivalently.
     * For FOPDT, deadtime affects the phase lag linearly.
     * dtheta_max ~ PM(deg) * pi/180 * tau_total / (tau_total * w_gc)
     * Simplified: dtheta_max ~ PM_current / w_gc / theta * 100%
     */
    double gm, pm, w180, w0db;
    pid_stability_margins(K, tau, theta, params, &gm, &pm, &w180, &w0db);
    if (pm > 0.0 && w0db > 1e-10) {
        *dtheta_max = (pm * M_PI / 180.0) / (w0db * theta + 1e-10) * 100.0;
        if (*dtheta_max > 500.0) *dtheta_max = 500.0;
    } else {
        *dtheta_max = 0.0;
    }
}
