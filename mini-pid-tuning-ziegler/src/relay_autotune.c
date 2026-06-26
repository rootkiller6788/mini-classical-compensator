/**
 * relay_autotune.c — Åström-Hägglund Relay Auto-Tuning Implementation
 *
 * Implements the relay feedback auto-tuning method for automatic
 * identification of ultimate gain Ku and ultimate period Pu.
 *
 * Knowledge: L1 (relay auto-tuner definition), L2 (describing function),
 *             L4 (oscillation condition via describing function),
 *             L5 (simulation and experiment algorithms),
 *             L7 (industrial auto-tuners: Honeywell UDC, Foxboro EXACT).
 *
 * Reference:
 *   Åström, K.J. & Hägglund, T. (1984), Automatica, 20(5), 645-651.
 *   Åström, K.J. & Hägglund, T. (1995), "PID Controllers", ISA, Ch.8.
 */

#include "relay_autotune.h"
#include "ziegler_nichols.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L2: Describing Function of Relay Nonlinearity
 * ────────────────────────────────────────────── */

void relay_describing_function(relay_type_t type, double amplitude,
                               double hysteresis, double input_amplitude,
                               double *mag, double *phase)
{
    /**
     * Ideal relay (on-off with no hysteresis):
     *   Output: u = +d if e > 0, u = -d if e < 0
     *   DF: N(a) = 4*d / (π * a)    [real, zero phase]
     *
     * Relay with hysteresis (±ε):
     *   Output: switches to +d when e > +ε
     *           switches to -d when e < -ε
     *   DF: N(a) = (4*d/(π*a)) * sqrt(1 - (ε/a)²)
     *              - j * (4*d*ε) / (π * a²)
     *
     * The negative imaginary part represents phase lag.
     *
     * Reference: Gelb & Vander Velde (1968), "Multiple-Input
     *            Describing Functions".
     */

    double d = amplitude;
    double a = input_amplitude;
    double eps = hysteresis;

    if (a < 1e-12) {
        *mag   = 0.0;
        *phase = 0.0;
        return;
    }

    if (type == RELAY_IDEAL || eps < 1e-12) {
        *mag   = (4.0 * d) / (M_PI * a);
        *phase = 0.0;
    } else {
        double ratio = eps / a;
        if (ratio > 1.0) {
            /* No switching: input amplitude below hysteresis */
            *mag   = 0.0;
            *phase = 0.0;
        } else {
            double re = (4.0 * d) / (M_PI * a) * sqrt(1.0 - ratio * ratio);
            double im = -(4.0 * d * eps) / (M_PI * a * a);
            *mag   = sqrt(re * re + im * im);
            *phase = atan2(im, re);
        }
    }
}

/* ──────────────────────────────────────────────
 * L5: Extract Ku/Pu from Oscillation
 * ────────────────────────────────────────────── */

int relay_extract_ultimate(const relay_config_t *config,
                           double osc_amplitude, double osc_period,
                           relay_result_t *result)
{
    if (!config || !result) return -1;
    if (osc_amplitude < 1e-12 || osc_period < 1e-12) return -1;

    double d = config->amplitude;
    double a = osc_amplitude;

    /* Describing function magnitude */
    double N_mag = 0.0, N_phase = 0.0;
    relay_describing_function(config->type, d, config->hysteresis,
                              a, &N_mag, &N_phase);

    /* Oscillation condition: N(a) * G(jω_osc) = -1
       → |G(jω_osc)| = 1 / |N(a)|
       → Ku = |N(a)|  (since Ku * |G| = 1 at oscillation) */

    result->Ku = N_mag;
    result->Pu = osc_period;
    result->oscillation_ampl = a;
    result->phase_cross_freq = 2.0 * M_PI / osc_period;
    result->cycles = 1;
    result->converged = 1;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Design Relay Experiment
 * ────────────────────────────────────────────── */

void relay_design_experiment(double approx_K, double approx_Tsum,
                             double u_range, relay_config_t *config)
{
    if (!config) return;

    (void)approx_K; /* Approximate gain estimated for oscillation amplitude prediction */

    /**
     * Design guidelines:
     *   d ≈ 0.05 * (u_max - u_min)  → small enough to not disrupt
     *   ε ≈ 0.01 * (setpoint range)  → small enough for accuracy
     *   Ts ← fast enough to capture oscillations
     *   max_duration ← long enough for 3-5 cycles to develop
     */

    config->amplitude  = 0.05 * u_range;
    config->hysteresis = 0.01 * u_range;

    /* Estimate ultimate period generously */
    double approx_Pu = 2.0 * approx_Tsum;
    if (approx_Pu < 0.1) approx_Pu = 0.1;

    config->Ts = approx_Pu / 50.0;  /* 50 samples per estimated period */
    if (config->Ts > 0.1) config->Ts = 0.1;

    config->max_duration = 10.0 * approx_Pu;
    config->max_cycles   = 5;
    config->settle_tol   = 0.02; /* 2% period convergence */
    config->type = RELAY_HYSTERESIS;
}

/* ──────────────────────────────────────────────
 * L5: Simulate Relay Feedback on FOPDT Model
 * ────────────────────────────────────────────── */

int relay_simulate_fopdt(const fopdt_model_t *fopdt,
                         const relay_config_t *config,
                         relay_result_t *result)
{
    if (!fopdt || !config || !result) return -1;

    /**
     * Emulate relay feedback loop digitally:
     *
     *   e(t) = -y(t)  (setpoint = 0, relay on error)
     *   u(t) = relay(e(t))
     *   y(t) = FOPDT response to u(t)
     *
     * The relay creates sustained oscillations at the process
     * phase crossover frequency (ω_180), from which Ku and Pu
     * are obtained.
     */

    double K = fopdt->K;
    double T = fopdt->T;
    double L = fopdt->L;
    double d = config->amplitude;
    double eps = config->hysteresis;

    if (K < 1e-12 || T < 1e-12) return -1;

    /* Simulation setup */
    int N_steps = (int)(config->max_duration / config->Ts);
    if (N_steps < 100) N_steps = 500;
    if (N_steps > 50000) N_steps = 50000;

    double Ts = config->Ts;

    /* State for FOPDT simulation (Euler with delay buffer) */
    int delay_steps = (int)(L / Ts);
    if (delay_steps < 0) delay_steps = 0;
    if (delay_steps > N_steps) delay_steps = N_steps;

    double *u_history = (double *)calloc(N_steps + 1, sizeof(double));
    double *y_history = (double *)calloc(N_steps + 1, sizeof(double));
    if (!u_history || !y_history) {
        free(u_history);
        free(y_history);
        return -1;
    }

    double y = 0.0;
    int relay_state = 1; /* +d */
    int crossing_count = 0;
    double last_crossing_time = 0.0;
    double period_sum = 0.0;
    int period_count = 0;
    double max_y = 0.0, min_y = 0.0;

    for (int k = 0; k < N_steps; k++) {
        double t = k * Ts;

        /* Relay action */
        double u;
        if (config->type == RELAY_HYSTERESIS) {
            if (y > eps) {
                relay_state = -1; /* Negative output */
            } else if (y < -eps) {
                relay_state = 1;  /* Positive output */
            }
            u = relay_state * d;
        } else {
            /* Ideal relay */
            u = (y < 0) ? d : -d;
        }

        u_history[k] = u;

        /* FOPDT: T * dy/dt + y = K * u(t-L) */
        int delay_idx = k - delay_steps;
        double u_delayed = (delay_idx >= 0) ? u_history[delay_idx] : 0.0;

        double dy = (Ts / T) * (K * u_delayed - y);
        y += dy;

        y_history[k] = y;

        /* Track oscillation */
        if (y > max_y) max_y = y;
        if (y < min_y) min_y = y;

        /* Zero-crossing detection for period measurement */
        if (k > 0 && y * y_history[k-1] < 0) {
            if (crossing_count > 1) {
                double period = t - last_crossing_time;
                /* Accept if period seems stable */
                if (crossing_count > 5) {
                    period_sum += period;
                    period_count++;
                }
            }
            last_crossing_time = t;
            crossing_count++;
        }
    }

    if (period_count < 2 || max_y - min_y < 1e-12) {
        free(u_history);
        free(y_history);
        result->converged = 0;
        return -1;
    }

    double avg_period = period_sum / period_count;
    double osc_ampl = (max_y - min_y) / 2.0;

    free(u_history);
    free(y_history);

    return relay_extract_ultimate(config, osc_ampl, avg_period, result);
}

/* ──────────────────────────────────────────────
 * L5: Full Auto-Tune Cycle
 * ────────────────────────────────────────────── */

int relay_autotune_complete(const fopdt_model_t *fopdt,
                            const relay_config_t *config,
                            pid_params_t *params)
{
    if (!fopdt || !config || !params) return -1;

    relay_result_t relay_res;
    memset(&relay_res, 0, sizeof(relay_res));

    if (relay_simulate_fopdt(fopdt, config, &relay_res) != 0) {
        return -1;
    }

    if (!relay_res.converged) return -1;

    return zn_freq_tune(relay_res.Ku, relay_res.Pu,
                        ZN_CONTROLLER_PID, params);
}

/* ──────────────────────────────────────────────
 * L5: Adaptive Re-Tune Trigger
 * ────────────────────────────────────────────── */

int relay_needs_retune(double current_Ku, double current_Pu,
                       double observed_Ku, double observed_Pu,
                       double K_threshold, double P_threshold)
{
    if (current_Ku < 1e-12 || current_Pu < 1e-12) return 1;

    double K_change = fabs(observed_Ku - current_Ku) / current_Ku;
    double P_change = fabs(observed_Pu - current_Pu) / current_Pu;

    return (K_change > K_threshold || P_change > P_threshold) ? 1 : 0;
}

/* ──────────────────────────────────────────────
 * L5: Variable Hysteresis Relay Method
 * ────────────────────────────────────────────── */

int relay_variable_hysteresis(const fopdt_model_t *fopdt,
                              double base_amplitude,
                              double eps_min, double eps_max, double eps_step,
                              double max_ampl, relay_result_t *result)
{
    if (!fopdt || !result) return -1;

    /**
     * Progressively increase hysteresis to reduce oscillation amplitude.
     * This protects the process from large excursions.
     */

    relay_config_t config;
    memset(&config, 0, sizeof(config));
    config.amplitude  = base_amplitude;
    config.Ts         = 0.01;
    config.max_duration = 100.0;
    config.max_cycles   = 5;
    config.settle_tol   = 0.02;
    config.type = RELAY_HYSTERESIS;

    int attempts = 0;
    double eps = eps_min;
    while (eps <= eps_max && attempts < 20) {
        config.hysteresis = eps;

        relay_result_t temp_res;
        memset(&temp_res, 0, sizeof(temp_res));

        if (relay_simulate_fopdt(fopdt, &config, &temp_res) == 0
            && temp_res.converged) {
            if (temp_res.oscillation_ampl <= max_ampl) {
                /* Found acceptable oscillation */
                *result = temp_res;
                result->converged = 1;
                return 0;
            }
        }
        eps *= eps_step;
        attempts++;
    }

    result->converged = 0;
    return -1;
}
