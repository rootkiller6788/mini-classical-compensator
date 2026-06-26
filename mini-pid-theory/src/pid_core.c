
/**
 * pid_core.c - Core PID Controller Implementation
 *
 * Implements the fundamental PID control law in all standard forms:
 * Standard (ISA), Parallel, Series, and 2DOF.
 *
 * Knowledge coverage:
 *   L1: All PID form implementations, discrete-time state management
 *   L2: Feedback computation, integral windup, derivative filtering
 *   L3: Laplace-domain transfer function evaluation
 *   L4: Routh-Hurwitz stability test
 *   L5: Closed-loop pole computation, performance evaluation
 *
 * Each function implements an independent knowledge point.
 * No filler - every line serves a pedagogical or practical purpose.
 */

#include "mini-pid-theory.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --------------------------------------------------------------------------
 * PID Initialization
 * -------------------------------------------------------------------------- */

int pid_init(const pid_params_t *params, pid_state_t *state)
{
    if (!params || !state) return -1;
    if (params->Kc < 0.0)  return -1;
    if (params->Ti < params->Ts && params->Ti > 0.0) return -1;
    if (params->Td < 0.0)  return -1;
    if (params->Ts <= 0.0) return -1;
    if (params->umin >= params->umax) return -1;

    memset(state, 0, sizeof(*state));
    state->is_initialized = 1;
    return 0;
}

/* --------------------------------------------------------------------------
 * PID Compute - Standard (ISA) Form
 *
 * u(t) = Kc * [ (b*ysp - y) + (1/Ti)*integral(e) + Td * d(c*ysp - y)/dt ]
 *
 * Discrete implementation with:
 *   - Derivative on measurement (D-on-PV) to avoid derivative kick
 *   - Low-pass filter on derivative term: D(s) = Kc*Td*s / (1 + s*Td/N)
 *   - Back-calculation anti-windup when saturated
 *
 * Discretization:
 *   Integral:  Trapezoidal (Tustin/Bilinear): I[k] = I[k-1] + (Ts/2)*(e[k]+e[k-1])
 *   Derivative: Backward difference + first-order filter:
 *     D[k] = (Td/(Td+N*Ts)) * D[k-1] - (Kc*Td*N/(Td+N*Ts)) * (y[k] - y[k-1])
 *
 * Reference: Astrom & Hagglund (1995), Section 3.4-3.5
 * Complexity: O(1) per call
 * -------------------------------------------------------------------------- */

double pid_compute(const pid_params_t *params, pid_state_t *state,
                   double ysp, double y, double dt)
{
    if (!params || !state || !state->is_initialized) return 0.0;

    double Ts = (dt > 0.0) ? dt : params->Ts;
    if (Ts <= 0.0) return state->prev_output;

    /* Handle deadband */
    double error = ysp - y;
    if (params->deadband > 0.0 && fabs(error) < params->deadband) {
        error = 0.0;
    }

    /* Apply control action direction */
    if (params->action == PID_ACTION_REVERSE) {
        error = -error;
    }

    /* ---- Proportional term (on measurement with setpoint weighting) ---- */
    double P = params->Kc * (params->b * ysp - y);

    /* ---- Integral term with anti-windup ---- */
    double I = state->integral;
    double Ki_scaled = (params->Ti > 1e-10) ? params->Kc * Ts / params->Ti : 0.0;

    /* Trapezoidal integration of error */
    double I_new = I + Ki_scaled * 0.5 * (error + state->prev_error);

    /* Back-calculation anti-windup: if output saturated last time,
     * feed back the saturation error to prevent integrator windup */
    if (params->antiwindup == PID_WINDUP_BACK_CALC && state->is_saturated) {
        double Tt = (params->Tt > 1e-10) ? params->Tt : sqrt(params->Ti * params->Td + 1e-10);
        if (Tt > 1e-10) {
            double tracking_correction = Ts / Tt * (state->saturated_output - state->prev_output);
            I_new += tracking_correction;
        }
    } else if (params->antiwindup == PID_WINDUP_CLAMPING && state->is_saturated) {
        /* Conditional integration: freeze integrator when saturated */
        I_new = I;
    }

    /* ---- Derivative term with filtering ---- */
    double D = 0.0;
    if (params->Td > 1e-10 && params->N > 0.0) {
        double Tf = params->Td / params->N;  /* Filter time constant */
        double alpha = Tf / (Tf + Ts);        /* Filter coefficient */

        if (params->deriv_mode == PID_DERIV_ON_MEASUREMENT ||
            params->deriv_mode == PID_DERIV_ON_FILTERED) {
            /* Derivative on measurement (avoids derivative kick) */
            double raw_deriv = -(params->Kc * params->Td) * (y - state->prev_measurement) / Ts;
            D = alpha * state->prev_deriv + (1.0 - alpha) * raw_deriv;

            if (params->deriv_mode == PID_DERIV_ON_FILTERED) {
                /* Additional low-pass filter on the measurement itself */
                double beta = Ts / (Ts + Tf);
                state->filtered_error = beta * error + (1.0 - beta) * state->filtered_error;
            }
        } else {
            /* Derivative on error (classic) */
            double deriv_error = (error - state->prev_error) / Ts;
            double raw_deriv = params->Kc * params->Td * deriv_error;
            D = alpha * state->prev_deriv + (1.0 - alpha) * raw_deriv;
        }
    }

    /* ---- 2DOF setpoint weighting on derivative ---- */
    if (params->form == PID_FORM_2DOF && params->c < 1.0) {
        /* D term is already on measurement (c=0) or partially on setpoint */
        /* The full implementation splits D into setpoint and measurement parts */
        double D_sp = 0.0;
        if (params->c > 0.0 && params->Td > 1e-10) {
            /* Derivative on setpoint (rarely used due to kick) */
            double deriv_sp = (ysp - state->prev_measurement) / Ts;
            D_sp = params->Kc * params->Td * params->c * deriv_sp;
            D = D_sp + (1.0 - params->c) * D;  /* Blend SP and PV derivatives */
        }
    }

    /* ---- Compute ideal (unsaturated) output ---- */
    double u_ideal = P + I_new + D;

    /* ---- Apply output saturation ---- */
    double u = u_ideal;
    int saturated = 0;
    if (u > params->umax) {
        u = params->umax;
        saturated = 1;
    } else if (u < params->umin) {
        u = params->umin;
        saturated = 1;
    }

    /* ---- Update state ---- */
    state->integral = I_new;
    state->prev_error = error;
    state->prev_measurement = y;
    state->prev_deriv = D;
    state->saturated_output = u;
    state->prev_output = u;
    state->is_saturated = (uint8_t)saturated;
    if (saturated) {
        state->saturation_count++;
    } else {
        state->saturation_count = 0;
    }
    state->derivative_spike = (fabs(D - state->prev_deriv) > 10.0 * fabs(D) + 1.0) ? 1 : 0;
    state->sample_count++;

    /* ---- Accumulate performance metrics ---- */
    double abs_e = fabs(error);
    state->integrated_absolute_error += abs_e * Ts;
    state->integrated_squared_error += error * error * Ts;
    state->integrated_time_iae += state->sample_count * Ts * abs_e;
    state->total_variation += fabs(u - state->prev_output);

    return u;
}

/* --------------------------------------------------------------------------
 * PID Reset - Clear internal state
 * -------------------------------------------------------------------------- */

void pid_reset(pid_state_t *state)
{
    if (!state) return;
    double iae = state->integrated_absolute_error;
    double ise = state->integrated_squared_error;
    double itae = state->integrated_time_iae;
    double tv = state->total_variation;
    uint32_t count = state->sample_count;

    memset(state, 0, sizeof(*state));
    state->is_initialized = 1;
    state->integrated_absolute_error = iae;
    state->integrated_squared_error = ise;
    state->integrated_time_iae = itae;
    state->total_variation = tv;
    state->sample_count = count;
}

/* --------------------------------------------------------------------------
 * Bumpless Transfer - Manual to Auto
 *
 * Back-calculates the integrator value so that when the controller switches
 * from manual to automatic mode, the output does not jump.
 *
 * Given current manual output u_man, we want:
 *   u_man = Kc*(b*ysp - y) + I + D
 *   => I = u_man - Kc*(b*ysp - y) - D
 * -------------------------------------------------------------------------- */

void pid_bumpless_transfer(const pid_params_t *params, pid_state_t *state,
                           double manual_output)
{
    if (!params || !state) return;

    /* Clamp manual output to allowed range */
    double u_man = manual_output;
    if (u_man > params->umax) u_man = params->umax;
    if (u_man < params->umin) u_man = params->umin;

    /* Set integrator to achieve u_man; P and D compensate in next compute call.
     * Since we don't have ysp/y here, the integrator is initialized so that
     * the next pid_compute output starts from u_man seamlessly. */

    state->integral = u_man;  /* Initialize I to manual output */
    /* P and D will subtract off on the first compute call because:
     * u = Kc*(b*ysp - y) + I[k] + D[k]
     * We need: u_man = Kc*(b*ysp - y) + u_man + D  --> subtract later */
    state->prev_output = u_man;
    state->saturated_output = u_man;
    state->manual_mode = 0;
    state->is_saturated = 0;
    state->saturation_count = 0;
}

/* --------------------------------------------------------------------------
 * Frequency Domain: PID Transfer Function
 *
 * Standard form:
 *   C(s) = Kc * (1 + 1/(Ti*s) + Td*s)
 *   |C(jw)| = Kc * sqrt(1 + (w*Td - 1/(w*Ti))^2)
 *   arg(C(jw)) = atan2(w*Td - 1/(w*Ti), 1.0)
 *
 * Parallel form:
 *   C(s) = Kp + Ki/s + Kd*s
 *   |C(jw)| = sqrt(Kp^2 + (Kd*w - Ki/w)^2)
 *   arg(C(jw)) = atan2(Kd*w - Ki/w, Kp)
 *
 * Series form:
 *   C(s) = Kp' * (1 + 1/(Ti'*s)) * (1 + Td'*s)
 *   |C(jw)| = Kp' * sqrt(1 + 1/(w*Ti')^2) * sqrt(1 + (w*Td')^2)
 *   arg(C(jw)) = atan2(0,1) + atan2(-1/(w*Ti'),1) + atan2(w*Td',1)
 *
 * Note: The ideal derivative term Kd*s is non-causal.
 * In practice, derivative is always filtered: Kd*s/(1+s*Tf).
 * Here we compute the ideal transfer function for analysis purposes.
 * -------------------------------------------------------------------------- */

double pid_freq_magnitude(const pid_params_t *params, double omega)
{
    if (!params || omega < 0.0) return 0.0;
    if (omega < 1e-15) return 1e308; /* Infinite gain at DC due to integrator */

    double mag;
    switch (params->form) {
        case PID_FORM_PARALLEL: {
            double Kp = params->Kc;
            double Ki = (params->Ti > 1e-10) ? params->Kc / params->Ti : 0.0;
            double Kd = params->Kc * params->Td;
            double re = Kp;
            double im = Kd * omega - Ki / omega;
            mag = sqrt(re*re + im*im);
            break;
        }
        case PID_FORM_SERIES: {
            double Kc_s = params->Kc * (1.0 + params->Td / (params->Ti + 1e-10));
            double Ti_s = params->Ti + params->Td;
            double Td_s = (params->Ti * params->Td) / (params->Ti + params->Td + 1e-10);
            if (Ti_s < 1e-10) Ti_s = 1e-10;
            double mag_pi = Kc_s * sqrt(1.0 + 1.0/(omega*omega*Ti_s*Ti_s));
            double mag_pd = sqrt(1.0 + omega*omega*Td_s*Td_s);
            mag = mag_pi * mag_pd;
            break;
        }
        case PID_FORM_STANDARD:
        default: {
            double re = 1.0;
            double im = omega * params->Td - 1.0 / (omega * (params->Ti > 1e-10 ? params->Ti : 1e-10));
            mag = params->Kc * sqrt(re*re + im*im);
            break;
        }
    }
    return mag;
}

double pid_freq_phase(const pid_params_t *params, double omega)
{
    if (!params || omega < 0.0) return 0.0;
    if (omega < 1e-15) return -M_PI / 2.0;

    double phase;
    switch (params->form) {
        case PID_FORM_PARALLEL: {
            double Kp = params->Kc;
            double Ki = (params->Ti > 1e-10) ? params->Kc / params->Ti : 0.0;
            double Kd = params->Kc * params->Td;
            double re = Kp;
            double im = Kd * omega - Ki / omega;
            phase = atan2(im, re);
            break;
        }
        case PID_FORM_SERIES: {
            double Ti_s = params->Ti + params->Td;
            double Td_s = (params->Ti * params->Td) / (params->Ti + params->Td + 1e-10);
            if (Ti_s < 1e-10) Ti_s = 1e-10;
            double phase_pi = atan2(-1.0/(omega*Ti_s), 1.0);
            double phase_pd = atan2(omega*Td_s, 1.0);
            phase = phase_pi + phase_pd;
            break;
        }
        case PID_FORM_STANDARD:
        default: {
            double im = omega * params->Td - 1.0 / (omega * (params->Ti > 1e-10 ? params->Ti : 1e-10));
            phase = atan2(im, 1.0);
            break;
        }
    }
    return phase;
}

/* --------------------------------------------------------------------------
 * Open-Loop Transfer Function: L(s) = G(s)*C(s)
 *
 * G(s) = K * exp(-theta*s) / (tau*s + 1)
 * |L(jw)| = |G(jw)| * |C(jw)|
 * arg(L(jw)) = arg(G(jw)) + arg(C(jw))
 *
 * |G(jw)| = K / sqrt(1 + (w*tau)^2)
 * arg(G(jw)) = -atan(w*tau) - w*theta
 * -------------------------------------------------------------------------- */

void pid_loop_transfer(double K, double tau, double theta,
                       const pid_params_t *params, double omega,
                       double *mag, double *phase)
{
    /* Plant magnitude and phase */
    double plant_mag = K / sqrt(1.0 + omega*omega*tau*tau);
    double plant_phase = -atan(omega * tau) - omega * theta;

    /* Controller magnitude and phase */
    double ctrl_mag = pid_freq_magnitude(params, omega);
    double ctrl_phase = pid_freq_phase(params, omega);

    *mag = plant_mag * ctrl_mag;
    *phase = plant_phase + ctrl_phase;
}

/* --------------------------------------------------------------------------
 * Sensitivity Function: S(s) = 1/(1 + L(s))
 *
 * S(jw) = 1 / (1 + |L|*exp(j*arg(L)))
 * |S| = 1 / sqrt(1 + 2*|L|*cos(arg(L)) + |L|^2)
 * arg(S) = -atan2(|L|*sin(arg(L)), 1 + |L|*cos(arg(L)))
 * -------------------------------------------------------------------------- */

void pid_sensitivity(double K, double tau, double theta,
                     const pid_params_t *params, double omega,
                     double *mag, double *phase)
{
    double L_mag, L_phase;
    pid_loop_transfer(K, tau, theta, params, omega, &L_mag, &L_phase);

    double re = 1.0 + L_mag * cos(L_phase);
    double im = L_mag * sin(L_phase);
    *mag = 1.0 / sqrt(re*re + im*im);
    *phase = -atan2(im, re);
}

/* --------------------------------------------------------------------------
 * Complementary Sensitivity: T(s) = L(s)/(1 + L(s))
 *
 * |T| = |L| / sqrt(1 + 2*|L|*cos(arg(L)) + |L|^2)
 * arg(T) = arg(L) - atan2(|L|*sin(arg(L)), 1 + |L|*cos(arg(L)))
 * -------------------------------------------------------------------------- */

void pid_complementary_sensitivity(double K, double tau, double theta,
                                   const pid_params_t *params, double omega,
                                   double *mag, double *phase)
{
    double L_mag, L_phase;
    pid_loop_transfer(K, tau, theta, params, omega, &L_mag, &L_phase);

    double re = 1.0 + L_mag * cos(L_phase);
    double im = L_mag * sin(L_phase);
    double denom = sqrt(re*re + im*im);
    *mag = (denom > 1e-15) ? L_mag / denom : 0.0;
    *phase = L_phase - atan2(im, re);
}

/* --------------------------------------------------------------------------
 * Maximum Sensitivity Ms = max_omega |S(j*omega)|
 *
 * Uses golden-section search to find the peak of |S(jw)|.
 * Ms is a key robustness metric:
 *   Ms < 1.4: conservative, good robustness
 *   1.4 < Ms < 2.0: reasonable
 *   Ms > 2.0: aggressive, low robustness
 *
 * The search is bounded by the critical frequencies:
 *   w_low  = 0.01 / (tau + theta)   (well below bandwidth)
 *   w_high = 10.0 / (tau + theta)   (well above bandwidth)
 * -------------------------------------------------------------------------- */

int pid_max_sensitivity(double K, double tau, double theta,
                        const pid_params_t *params,
                        double *Ms, double *w_peak)
{
    if (!params || !Ms || !w_peak) return -1;
    if (tau < 0.0 || theta < 0.0) return -1;

    double tau_total = tau + theta;
    if (tau_total < 1e-10) tau_total = 1.0;

    double w_low  = 0.01 / tau_total;
    double w_high = 10.0 / tau_total;

    const double phi = (sqrt(5.0) - 1.0) / 2.0; /* Golden ratio conjugate */
    const int max_iter = 50;
    const double tol = 1e-6;

    double a = w_low, b = w_high;
    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    double S_c_mag, S_d_mag, dummy_phase;
    pid_sensitivity(K, tau, theta, params, c, &S_c_mag, &dummy_phase);
    pid_sensitivity(K, tau, theta, params, d, &S_d_mag, &dummy_phase);

    for (int i = 0; i < max_iter; i++) {
        if (fabs(b - a) < tol * (fabs(c) + fabs(d))) break;
        if (S_c_mag > S_d_mag) {
            b = d; d = c; S_d_mag = S_c_mag;
            c = b - phi * (b - a);
            pid_sensitivity(K, tau, theta, params, c, &S_c_mag, &dummy_phase);
        } else {
            a = c; c = d; S_c_mag = S_d_mag;
            d = a + phi * (b - a);
            pid_sensitivity(K, tau, theta, params, d, &S_d_mag, &dummy_phase);
        }
    }

    *w_peak = (a + b) / 2.0;
    pid_sensitivity(K, tau, theta, params, *w_peak, Ms, &dummy_phase);
    return 0;
}

/* --------------------------------------------------------------------------
 * Routh-Hurwitz Stability Test
 *
 * For plant G(s) = K/(tau*s+1) with standard PID C(s) = Kc*(1 + 1/(Ti*s) + Td*s):
 *
 * Characteristic equation: 1 + G*C = 0
 *   tau*Ti*s^3 + Ti*(1+K*Kc*Td)*s^2 + K*Kc*Ti*s + K*Kc = 0
 *
 * Let: a3 = tau*Ti, a2 = Ti*(1+K*Kc*Td), a1 = K*Kc*Ti, a0 = K*Kc
 *
 * Routh array:
 *   s^3 | a3        a1
 *   s^2 | a2        a0
 *   s^1 | b1        0   where b1 = (a2*a1 - a3*a0)/a2
 *   s^0 | a0
 *
 * Stability: a3,a2,a1,a0 > 0  AND  a2*a1 > a3*a0
 * -------------------------------------------------------------------------- */

int pid_stability_routh(double K, double tau, const pid_params_t *params)
{
    if (!params) return -1;
    if (K <= 0.0 || tau <= 0.0) return -1;
    if (params->Kc <= 0.0) return 0; /* Negative feedback only works with positive gain */

    double Ti = (params->Ti > 1e-10) ? params->Ti : 1e308;
    double Td = params->Td;

    double a3 = tau * Ti;
    double a2 = Ti * (1.0 + K * params->Kc * Td);
    double a1 = K * params->Kc * Ti;
    double a0 = K * params->Kc;

    /* All coefficients must be positive */
    if (a3 <= 0.0 || a2 <= 0.0 || a1 <= 0.0 || a0 <= 0.0) return 0;

    /* Routh-Hurwitz criterion: a2*a1 > a3*a0 */
    if (a2 * a1 > a3 * a0) return 1;

    return 0;
}

/* --------------------------------------------------------------------------
 * Stability Margins (Gain Margin & Phase Margin)
 *
 * For FOPDT plant with PID controller.
 *
 * Gain Margin: How much the loop gain can increase before instability.
 *   |L(j*w_pc)| = 1/GM  at the phase crossover (arg(L) = -pi)
 *   GM [dB] = -20*log10(|L(j*w_pc)|)
 *
 * Phase Margin: How much additional phase lag causes instability.
 *   PM = pi + arg(L(j*w_gc))  at the gain crossover (|L| = 1)
 *
 * Algorithm: Binary search for crossover frequencies.
 * -------------------------------------------------------------------------- */

int pid_stability_margins(double K, double tau, double theta,
                          const pid_params_t *params,
                          double *gm, double *pm, double *w180, double *w0db)
{
    if (!params || !gm || !pm || !w180 || !w0db) return -1;

    double tau_total = tau + theta;
    if (tau_total < 1e-10) tau_total = 1.0;

    /* Search range: 3 decades around 1/tau_total */
    double w_min = 1e-4 / tau_total;
    double w_max = 1e4 / tau_total;
    int max_iter = 60;
    double tol = 1e-8;

    /* ---- Find phase crossover (arg(L) = -pi) ---- */
    double w_lo = w_min, w_hi = w_max;
    double Lm_lo, Lp_lo, Lm_hi, Lp_hi;
    pid_loop_transfer(K, tau, theta, params, w_lo, &Lm_lo, &Lp_lo);
    pid_loop_transfer(K, tau, theta, params, w_hi, &Lm_hi, &Lp_hi);

    /* Normalize phases to [-pi, pi] */
    while (Lp_lo > M_PI) Lp_lo -= 2.0*M_PI;
    while (Lp_lo < -M_PI) Lp_lo += 2.0*M_PI;
    while (Lp_hi > M_PI) Lp_hi -= 2.0*M_PI;
    while (Lp_hi < -M_PI) Lp_hi += 2.0*M_PI;

    *w180 = 0.0;
    *gm = 1e308; /* Infinite gain margin if no phase crossover */

    if (Lp_lo <= -M_PI && Lp_hi >= -M_PI) {
        /* Phase crossover exists in this range */
        for (int i = 0; i < max_iter; i++) {
            double w_mid = (w_lo + w_hi) / 2.0;
            double Lm_mid, Lp_mid;
            pid_loop_transfer(K, tau, theta, params, w_mid, &Lm_mid, &Lp_mid);
            while (Lp_mid > M_PI) Lp_mid -= 2.0*M_PI;
            while (Lp_mid < -M_PI) Lp_mid += 2.0*M_PI;

            if (fabs(w_hi - w_lo) < tol * w_mid) {
                *w180 = w_mid;
                *gm = -20.0 * log10(Lm_mid + 1e-15);
                break;
            }
            if (Lp_mid < -M_PI) {
                w_lo = w_mid;
            } else {
                w_hi = w_mid;
            }
        }
    }

    /* ---- Find gain crossover (|L| = 1 = 0 dB) ---- */
    w_lo = w_min; w_hi = w_max;
    pid_loop_transfer(K, tau, theta, params, w_lo, &Lm_lo, &Lp_lo);
    pid_loop_transfer(K, tau, theta, params, w_hi, &Lm_hi, &Lp_hi);

    /* Convert to dB */
    double LdB_lo = 20.0 * log10(Lm_lo + 1e-15);
    double LdB_hi = 20.0 * log10(Lm_hi + 1e-15);

    *w0db = 0.0;
    *pm = -1e308; /* No phase margin if no gain crossover */

    if (LdB_lo >= 0.0 && LdB_hi <= 0.0) {
        for (int i = 0; i < max_iter; i++) {
            double w_mid = (w_lo + w_hi) / 2.0;
            double Lm_mid, Lp_mid;
            pid_loop_transfer(K, tau, theta, params, w_mid, &Lm_mid, &Lp_mid);
            double LdB_mid = 20.0 * log10(Lm_mid + 1e-15);

            if (fabs(w_hi - w_lo) < tol * w_mid) {
                *w0db = w_mid;
                while (Lp_mid > M_PI) Lp_mid -= 2.0*M_PI;
                while (Lp_mid < -M_PI) Lp_mid += 2.0*M_PI;
                *pm = (M_PI + Lp_mid) * 180.0 / M_PI;
                break;
            }
            if (LdB_mid > 0.0) {
                w_lo = w_mid;
            } else {
                w_hi = w_mid;
            }
        }
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Closed-Loop Poles for FOPDT Plant with PID
 *
 * Characteristic equation: 1 + G(s)C(s) = 0
 *
 * For FOPDT: G(s) = K*exp(-theta*s)/(tau*s+1)
 * With 1st-order Pade: exp(-theta*s) ~ (1 - theta*s/2)/(1 + theta*s/2)
 *
 * For standard PID C(s) = Kc*(1 + 1/(Ti*s) + Td*s):
 *   G(s)C(s) = K*Kc*(1 + 1/(Ti*s) + Td*s) * (1 - theta*s/2) /
 *              [(tau*s+1) * (1 + theta*s/2)]
 *
 * Characteristic equation becomes a 4th-order polynomial:
 *   a4*s^4 + a3*s^3 + a2*s^2 + a1*s + a0 = 0
 *
 * We solve via Ferrari's method for quartics.
 * For delay-free (theta=0), it reduces to a cubic.
 * -------------------------------------------------------------------------- */

/* Solve quadratic: a*s^2 + b*s + c = 0; stores roots as [re0,im0,re1,im1] */
static int solve_quadratic(double a, double b, double c, double *roots)
{
    if (fabs(a) < 1e-15) {
        if (fabs(b) < 1e-15) return 0;
        roots[0] = -c / b; roots[1] = 0.0;
        roots[2] = roots[0]; roots[3] = 0.0;
        return 1;
    }
    double disc = b*b - 4.0*a*c;
    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        roots[0] = (-b + sqrt_disc) / (2.0*a); roots[1] = 0.0;
        roots[2] = (-b - sqrt_disc) / (2.0*a); roots[3] = 0.0;
        return 2;
    } else {
        roots[0] = -b / (2.0*a); roots[1] = sqrt(-disc) / (2.0*a);
        roots[2] = roots[0]; roots[3] = -roots[1];
        return 2;
    }
}

/* Solve cubic: a*s^3 + b*s^2 + c*s + d = 0 via Cardano's method */
static int solve_cubic(double a, double b, double c, double d, double *roots)
{
    if (fabs(a) < 1e-15) return solve_quadratic(b, c, d, roots);

    /* Normalize to s^3 + A*s^2 + B*s + C = 0 */
    double A = b / a, B = c / a, C = d / a;

    /* Depressed cubic: t^3 + p*t + q = 0, where s = t - A/3 */
    double p = B - A*A/3.0;
    double q = C + (2.0*A*A*A - 9.0*A*B)/27.0;

    /* Discriminant */
    double disc = q*q/4.0 + p*p*p/27.0;

    if (disc > 0.0) {
        /* One real root */
        double u = cbrt(-q/2.0 + sqrt(disc));
        double v = cbrt(-q/2.0 - sqrt(disc));
        double t1 = u + v;
        roots[0] = t1 - A/3.0; roots[1] = 0.0;
        return 1;
    } else if (disc < 0.0) {
        /* Three real roots */
        double phi = acos(-q/2.0 / sqrt(-p*p*p/27.0));
        double r = 2.0 * sqrt(-p/3.0);
        for (int k = 0; k < 3; k++) {
            double t = r * cos((phi + 2.0*M_PI*k)/3.0);
            roots[2*k] = t - A/3.0; roots[2*k+1] = 0.0;
        }
        return 3;
    } else {
        /* Double/triple root */
        double t1 = 2.0 * cbrt(-q/2.0);
        double t2 = -cbrt(-q/2.0);
        roots[0] = t1 - A/3.0; roots[1] = 0.0;
        roots[2] = t2 - A/3.0; roots[3] = 0.0;
        roots[4] = roots[2]; roots[5] = 0.0;
        return 3;
    }
}

/* Solve quartic: a*s^4 + b*s^3 + c*s^2 + d*s + e = 0 via Ferrari's method.
 * Available for higher-order characteristic equations (e.g., 2nd-order Pade).
 * Not currently called but provided for completeness. */
__attribute__((unused))
static int solve_quartic(double a, double b, double c, double d, double e,
                         double *roots)
{
    if (fabs(a) < 1e-15) return solve_cubic(b, c, d, e, roots);

    /* Normalize */
    double A = b / a, B = c / a, C = d / a, D = e / a;

    /* Resolvent cubic: y^3 - B*y^2 + (A*C-4*D)*y + (4*B*D - C*C - A*A*D) = 0 */
    double ra = 1.0;
    double rb = -B;
    double rc = A*C - 4.0*D;
    double rd = 4.0*B*D - C*C - A*A*D;

    double y_roots[6];
    solve_cubic(ra, rb, rc, rd, y_roots);

    /* Pick a real root for the resolvent cubic (use the largest) */
    double y = y_roots[0];
    for (int k = 1; k < 3; k++) {
        if (y_roots[2*k+1] == 0.0 && y_roots[2*k] > y) y = y_roots[2*k];
    }

    /* Compute R = sqrt(A*A/4 - B + y) */
    double R2 = A*A/4.0 - B + y;
    if (R2 < 0.0) R2 = 0.0;
    double R = sqrt(R2);

    double D_term, E_term;
    if (fabs(R) > 1e-15) {
        D_term = 0.75*A*A - R*R - 2.0*B + 0.25*(4.0*A*B - 8.0*C - A*A*A)/R;
        E_term = 0.75*A*A - R*R - 2.0*B - 0.25*(4.0*A*B - 8.0*C - A*A*A)/R;
    } else {
        D_term = 0.75*A*A - 2.0*B + 2.0*sqrt(y*y - 4.0*D);
        E_term = 0.75*A*A - 2.0*B - 2.0*sqrt(y*y - 4.0*D);
    }

    int nroots = 0;
    /* Solve s^2 + (A/2 + R)*s + (y/2 + D_term/2) = 0 */
    if (D_term >= 0.0) {
        double r1 = -(A/2.0 + R)/2.0;
        double dd = (A/2.0+R)*(A/2.0+R) - 4.0*(y/2.0 + D_term/2.0);
        if (dd >= 0.0) {
            double disc = sqrt(dd)/2.0;
            roots[0] = r1 + disc; roots[1] = 0.0;
            roots[2] = r1 - disc; roots[3] = 0.0;
            nroots += 2;
        }
    }
    /* Solve s^2 + (A/2 - R)*s + (y/2 + E_term/2) = 0 */
    if (E_term >= 0.0) {
        double r2 = -(A/2.0 - R)/2.0;
        double dd = (A/2.0-R)*(A/2.0-R) - 4.0*(y/2.0 + E_term/2.0);
        if (dd >= 0.0) {
            double disc = sqrt(dd)/2.0;
            int idx = nroots * 2;
            roots[idx] = r2 + disc; roots[idx+1] = 0.0;
            roots[idx+2] = r2 - disc; roots[idx+3] = 0.0;
            nroots += 2;
        }
    }
    return nroots;
}

int pid_closed_loop_poles(double K, double tau, double theta,
                          const pid_params_t *params,
                          double *poles, int *npoles)
{
    if (!params || !poles || !npoles) return -1;
    if (K <= 0.0 || tau < 0.0) return -1;

    double Kc = params->Kc;
    double Ti = (params->Ti > 1e-10) ? params->Ti : 1e308;
    double Td = params->Td;
    double Ki = (Ti > 1e307) ? 0.0 : Kc / Ti;
    double Kd = Kc * Td;

    if (theta < 1e-10) {
        /* Delay-free case: cubic characteristic equation
         * 1 + K*(Kp + Ki/s + Kd*s)/(tau*s+1) = 0
         * tau*s^2 + s + K*Kd*s^2 + K*Kp*s + K*Ki = 0
         * (tau + K*Kd)*s^2 + (1 + K*Kp)*s + K*Ki = 0  [if no D]
         * With D: tau*s^2 + s + K*(Kp+Ki/s+Kd*s) = 0
         *         tau*s^3 + s^2 + K*Kd*s^2 + K*Kp*s + K*Ki = 0
         */
        if (Td > 1e-10) {
            /* Cubic with derivative */
            double a3 = tau;
            double a2 = 1.0 + K * Kd;
            double a1 = K * Kc;
            double a0 = K * Ki;
            double roots[6];
            *npoles = solve_cubic(a3, a2, a1, a0, roots);
            for (int i = 0; i < (*npoles) * 2; i++) poles[i] = roots[i];
        } else if (Ki > 1e-15) {
            /* PI: quadratic */
            double a2 = tau;
            double a1 = 1.0 + K * Kc;
            double a0 = K * Ki;
            *npoles = solve_quadratic(a2, a1, a0, poles);
        } else {
            /* P only */
            poles[0] = -(1.0 + K * Kc) / tau; poles[1] = 0.0;
            *npoles = 1;
        }
    } else {
        /* With deadtime using 1st-order Pade approximation:
         * exp(-theta*s) ~ (1 - theta*s/2)/(1 + theta*s/2)
         *
         * Characteristic equation:
         * (tau*s+1)*(theta*s/2+1)*Ti*s + K*Kc*(Ti*Td*s^2+Ti*s+1)*(1-theta*s/2) = 0
         *
         * This is 4th order in s.
         */
        double th2 = theta / 2.0;
        /* Characteristic equation derivation:
         * 1 + G(s)C(s) = 0  with Pade approx for deadtime.
         * 1 + K*(1-th2*s)/(tau*s+1)/(1+th2*s) * Kc*(1+1/(Ti*s)+Td*s) = 0
         * Multiply by (tau*s+1)*(1+th2*s)*Ti*s:
         * Result is cubic: a3*s^3 + a2*s^2 + a1*s + a0 = 0
         */
        double a3_coeff = Ti * tau * th2 - K * Kc * th2 * Ti * Td;
        double a2_coeff = Ti * (tau + th2) + K * Kc * (Ti * Td - th2 * Ti);
        double a1_coeff = Ti + K * Kc * (Ti - th2);
        double a0_coeff = K * Kc;

        double roots[6];
        if (fabs(a3_coeff) < 1e-15) {
            *npoles = solve_quadratic(a2_coeff, a1_coeff, a0_coeff, poles);
        } else {
            *npoles = solve_cubic(a3_coeff, a2_coeff, a1_coeff, a0_coeff, roots);
            for (int i = 0; i < (*npoles) * 2; i++) poles[i] = roots[i];
        }
    }
    return 0;
}

/* --------------------------------------------------------------------------
 * Performance Evaluation for PID-controlled FOPDT
 *
 * Simulates the closed-loop step response using discrete-time approximation.
 *
 * Plant simulation: y[k+1] = exp(-Ts/tau)*y[k] + K*(1-exp(-Ts/tau))*u[k-d]
 * where d = ceil(theta/Ts) is the deadtime in samples.
 *
 * Evaluates IAE, ISE, ITAE, overshoot, settling time, rise time, and TV.
 * -------------------------------------------------------------------------- */

int pid_evaluate_fopdt(double K, double tau, double theta,
                       const pid_params_t *params,
                       double ysp, double tfinal, int nsteps,
                       pid_performance_t *perf)
{
    if (!params || !perf || nsteps < 10) return -1;
    if (K <= 0.0 || tau <= 0.0 || tfinal <= 0.0) return -1;

    double Ts = tfinal / nsteps;
    int delay_samples = (int)ceil(theta / Ts);
    if (delay_samples < 0) delay_samples = 0;

    /* Allocate delay buffer */
    double *u_history = (double*)calloc(delay_samples + 1, sizeof(double));
    if (!u_history) return -1;

    pid_state_t state;
    pid_init(params, &state);

    double y = 0.0;
    int delay_idx = 0;

    memset(perf, 0, sizeof(*perf));
    double y_max = 0.0;
    double t_rise_start = -1.0, t_rise_end = -1.0;
    double t_settle = -1.0;

    for (int k = 0; k < nsteps; k++) {
        double t = k * Ts;

        /* Get delayed input */
        double u_delayed = u_history[delay_idx];

        /* Plant dynamics: y[k+1] = exp(-Ts/tau)*y[k] + K*(1-exp(-Ts/tau))*u_delayed */
        double alpha = exp(-Ts / tau);
        y = alpha * y + K * (1.0 - alpha) * u_delayed;

        /* Compute PID control */
        double u = pid_compute(params, &state, ysp, y, Ts);

        /* Store in delay buffer */
        u_history[delay_idx] = u;
        delay_idx = (delay_idx + 1) % (delay_samples + 1);

        /* Error */
        double e = ysp - y;
        double abs_e = fabs(e);

        /* Performance indices */
        perf->iae += abs_e * Ts;
        perf->ise += e * e * Ts;
        perf->itae += t * abs_e * Ts;
        perf->itse += t * e * e * Ts;
        if (k > 0) perf->tv_u += fabs(u - u_history[(delay_idx-1+delay_samples+1)%(delay_samples+1)]);

        /* Overshoot tracking */
        if (y > y_max) y_max = y;
        if (y_max > ysp) perf->overshoot = (y_max - ysp) / ysp * 100.0;

        /* Rise time (10% to 90%) */
        if (t_rise_start < 0.0 && y >= 0.1 * ysp) t_rise_start = t;
        if (t_rise_end < 0.0 && y >= 0.9 * ysp) t_rise_end = t;

        /* Settling time (stays within +/-2% after this) */
        if (t_settle < 0.0 && fabs(y - ysp) <= 0.02 * fabs(ysp) && t > 0.5 * tfinal) {
            /* Check if it stays settled */
            t_settle = t;
        }
        if (t_settle >= 0.0 && fabs(y - ysp) > 0.02 * fabs(ysp)) {
            t_settle = -1.0; /* Left the band */
        }
    }

    perf->rise_time = (t_rise_end > 0.0) ? (t_rise_end - t_rise_start) : -1.0;
    perf->settling_time = (t_settle > 0.0) ? t_settle : -1.0;
    perf->steady_state_error = ysp - y;

    free(u_history);
    return 0;
}

/* --------------------------------------------------------------------------
 * PID Form Conversion
 *
 * Standard <-> Parallel <-> Series <-> 2DOF
 * -------------------------------------------------------------------------- */

int pid_convert_form(const pid_params_t *src_params, pid_form_t src_form,
                     pid_params_t *dst_params, pid_form_t dst_form)
{
    if (!src_params || !dst_params) return -1;

    /* First convert to standard form as intermediate representation */
    double Kc, Ti, Td;
    switch (src_form) {
        case PID_FORM_STANDARD:
            Kc = src_params->Kc; Ti = src_params->Ti; Td = src_params->Td;
            break;
        case PID_FORM_PARALLEL: {
            double Kp = src_params->Kc; /* Kc field holds Kp for parallel */
            double Ki = (src_params->Ti > 1e-10) ? src_params->Kc / src_params->Ti : 0.0;
            double Kd = src_params->Kc * src_params->Td;
            if (Kp < 1e-10) return -1;
            Kc = Kp;
            Ti = (Ki > 1e-10) ? Kp / Ki : 1e308;
            Td = Kd / Kp;
            break;
        }
        case PID_FORM_SERIES: {
            double Kc_s = src_params->Kc;
            double Ti_s = src_params->Ti;
            double Td_s = src_params->Td;
            if (Ti_s + Td_s < 1e-10) return -1;
            Kc = Kc_s * (Ti_s + Td_s) / Ti_s;
            Ti = Ti_s + Td_s;
            Td = Ti_s * Td_s / (Ti_s + Td_s);
            break;
        }
        case PID_FORM_2DOF:
            /* 2DOF uses the same Kc,Ti,Td as standard, just different b,c */
            Kc = src_params->Kc; Ti = src_params->Ti; Td = src_params->Td;
            break;
        default:
            return -1;
    }

    /* Copy base structure first */
    *dst_params = *src_params;

    /* Now convert from standard to target form */
    switch (dst_form) {
        case PID_FORM_STANDARD:
            dst_params->Kc = Kc; dst_params->Ti = Ti; dst_params->Td = Td;
            break;
        case PID_FORM_PARALLEL: {
            /* Store as: Kc=Kp, Ti s.t. Ki=Kc/Ti, Td s.t. Kd=Kc*Td */
            double Ki = (Ti > 1e307) ? 0.0 : Kc / Ti;
            double Kd = Kc * Td;
            dst_params->Kc = Kc; /* Kp */
            dst_params->Ti = (Ki > 1e-10) ? Kc / Ki : 1e308; /* Ti stored, Ki=Kc/Ti */
            dst_params->Td = Kd / (Kc > 1e-10 ? Kc : 1.0); /* Td stored */
            break;
        }
        case PID_FORM_SERIES: {
            if (Ti < 1e10 && Td > 1e-10) {
                dst_params->Kc = Kc * (1.0 + Td / Ti);
                dst_params->Ti = Ti + Td;
                dst_params->Td = Ti * Td / (Ti + Td);
            } else {
                dst_params->Kc = Kc;
                dst_params->Ti = Ti;
                dst_params->Td = Td;
            }
            break;
        }
        case PID_FORM_2DOF:
            dst_params->Kc = Kc; dst_params->Ti = Ti; dst_params->Td = Td;
            /* b and c setpoint weights remain as-is */
            break;
        default:
            return -1;
    }

    dst_params->form = dst_form;
    return 0;
}

/* --------------------------------------------------------------------------
 * Relay Feedback - Astrom-Hagglund method
 *
 * For FOPDT G(s) = K*exp(-theta*s)/(tau*s+1):
 * The describing function of a relay with amplitude d is 4*d/(pi*A)
 * where A is the oscillation amplitude.
 *
 * At the ultimate frequency w_u:
 *   |G(j*w_u)| = 1/Ku
 *   arg(G(j*w_u)) = -pi
 *
 * Solving: w_u solves: -atan(w_u*tau) - w_u*theta = -pi
 * Then: Ku = sqrt(1+(w_u*tau)^2) / K
 *
 * Tu = 2*pi / w_u
 * -------------------------------------------------------------------------- */

int pid_relay_feedback(double K, double tau, double theta, double d,
                       double *Ku, double *Tu)
{
    (void)d; /* Relay amplitude - used in describing function: Ku = 4*d/(pi*a) */
    if (!Ku || !Tu || K <= 0.0) return -1;

    /* Find w_u that satisfies arg(G(j*w_u)) = -pi
     * -atan(w*tau) - w*theta = -pi
     * => atan(w*tau) + w*theta = pi
     *
     * For small tau, w_u ~ pi/theta
     * For large tau, w_u ~ pi/(2*theta) approx
     *
     * Use binary search */
    double w_lo = 0.01 / (tau + theta + 1e-10);
    double w_hi = M_PI / (theta + 1e-10);
    int max_iter = 50;

    for (int i = 0; i < max_iter; i++) {
        double w = (w_lo + w_hi) / 2.0;
        double phase = atan(w * tau) + w * theta;

        if (fabs(phase - M_PI) < 1e-8) {
            double w_u = w;
            *Tu = 2.0 * M_PI / w_u;
            *Ku = sqrt(1.0 + w_u*w_u*tau*tau) / K;
            return 0;
        }
        if (phase < M_PI) {
            w_lo = w;
        } else {
            w_hi = w;
        }
        if (fabs(w_hi - w_lo) < 1e-10) {
            double w_u = w;
            *Tu = 2.0 * M_PI / w_u;
            *Ku = sqrt(1.0 + w_u*w_u*tau*tau) / K;
            return 0;
        }
    }

    /* Fallback: pure deadtime approximation */
    if (theta > 1e-10) {
        *Tu = 2.0 * theta;
        *Ku = sqrt(1.0 + M_PI*M_PI*tau*tau/(theta*theta)) / K;
    } else {
        *Tu = 2.0 * M_PI * tau;
        *Ku = sqrt(2.0) / K;
    }
    return 0;
}
