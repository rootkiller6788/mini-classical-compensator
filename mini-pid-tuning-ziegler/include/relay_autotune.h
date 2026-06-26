#ifndef RELAY_AUTOTUNE_H
#define RELAY_AUTOTUNE_H

/**
 * relay_autotune.h — Åström-Hägglund Relay Auto-Tuning (1984)
 *
 * Knowledge coverage:
 *   L1: Relay auto-tuning definition — automated identification of Ku, Pu
 *   L2: Describing function analysis for relay nonlinearity
 *   L4: Oscillation condition: N(a)*G(jω_osc) = -1
 *   L5: Relay experiment implementation and parameter extraction
 *   L6: Application to industrial PID auto-tuners (Honeywell, Foxboro)
 *   L8: Advanced: adaptive PID with periodic relay updates
 *
 * Reference:
 *   Åström, K.J. & Hägglund, T. (1984)
 *   "Automatic Tuning of Simple Regulators with Specifications on
 *    Phase and Amplitude Margins", Automatica, 20(5), 645-651.
 *
 *   Åström, K.J. & Hägglund, T. (1995)
 *   "PID Controllers: Theory, Design, and Tuning", ISA.
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: Relay Auto-Tuner Configuration
 * ────────────────────────────────────────────── */

/**
 * Relay types for auto-tuning experiments.
 *
 * IDEAL:      u = +d if e > 0, u = -d if e < 0
 * HYSTERESIS: u = +d if e > ε, u = -d if e < -ε (dead zone ±ε)
 * SATURATING: u = ±min(d, |e|) — limits the relay amplitude
 */
typedef enum {
    RELAY_IDEAL = 0,
    RELAY_HYSTERESIS,
    RELAY_SATURATING
} relay_type_t;

/**
 * L1: Relay experiment configuration.
 */
typedef struct {
    double amplitude;      /* Relay output amplitude d (> 0) */
    double hysteresis;     /* Hysteresis width ε (≥ 0; 0 = ideal relay) */
    double Ts;             /* Sampling interval */
    double max_duration;   /* Maximum experiment duration [time] */
    int    max_cycles;     /* Maximum number of oscillation cycles */
    double settle_tol;     /* Convergence tolerance for period detection */
    relay_type_t type;     /* Relay type */
} relay_config_t;

/**
 * L1: Relay experiment result — contains identified parameters.
 */
typedef struct {
    double Ku;             /* Ultimate gain */
    double Pu;             /* Ultimate period [time] */
    double oscillation_ampl; /* Output oscillation amplitude */
    int    cycles;         /* Number of cycles observed */
    double phase_cross_freq; /* Phase crossover frequency ω_180 (rad/time) */
    int    converged;      /* 1 = experiment converged successfully */
} relay_result_t;

/* ──────────────────────────────────────────────
 * L2: Describing Function Analysis
 * ────────────────────────────────────────────── */

/**
 * Compute the describing function (DF) of the relay nonlinearity.
 *
 * For an ideal relay with amplitude d, input amplitude a:
 *
 *   N(a) = 4*d / (π * a)
 *
 * This is the fundamental harmonic approximation. The oscillation
 * condition is:
 *
 *   N(a) * G(jω_osc) = -1
 *
 * For a relay with hysteresis width ε:
 *
 *   N(a) = (4*d / (π*a)) * sqrt(1 - (ε/a)²) - j*(4*d*ε) / (π*a²)
 *
 * Complexity: O(1).
 *
 * @param type            Relay type.
 * @param amplitude       Relay amplitude d.
 * @param hysteresis      Hysteresis width ε.
 * @param input_amplitude Estimated input amplitude a.
 * @param mag             [out] |N(a)|.
 * @param phase           [out] ∠N(a) [radians].
 */
void relay_describing_function(relay_type_t type, double amplitude,
                               double hysteresis, double input_amplitude,
                               double *mag, double *phase);

/**
 * Extract Ku and Pu from relay experiment oscillation data.
 *
 * Given measured oscillation amplitude a and period P_osc:
 *
 * Ideal relay:
 *   Ku = 4*d / (π * a)
 *   Pu = P_osc
 *
 * Hysteresis relay:
 *   Ku = (4*d / (π*a)) * sqrt(1 - (ε/a)²)
 *   Pu = compensates for additional phase lag from hysteresis.
 *
 * @param config         Relay configuration.
 * @param osc_amplitude  Measured oscillation amplitude (a).
 * @param osc_period     Measured oscillation period (P_osc).
 * @param result         [out] Computed Ku and Pu.
 * @return               0 on success.
 */
int relay_extract_ultimate(const relay_config_t *config,
                           double osc_amplitude, double osc_period,
                           relay_result_t *result);

/* ──────────────────────────────────────────────
 * L5: Relay Auto-Tuning Algorithms
 * ────────────────────────────────────────────── */

/**
 * Design relay experiment parameters based on process characteristics.
 *
 * This function determines appropriate relay amplitude d and hysteresis ε
 * given rough process knowledge (approximate gain, time scale).
 *
 * Rules of thumb (Åström & Hägglund, 1995):
 *   d ≈ 0.05 * (u_max - u_min)   (5% of actuator range)
 *   ε ≈ 0.01 * (setpoint range)   (1% of measurement range)
 *
 * @param approx_K    Approximate process gain estimate.
 * @param approx_Tsum Approximate sum of time constants (L+T).
 * @param u_range     Actuator output range = u_max - u_min.
 * @param config      [out] Recommended relay configuration.
 */
void relay_design_experiment(double approx_K, double approx_Tsum,
                             double u_range, relay_config_t *config);

/**
 * Simulate relay feedback experiment for a given FOPDT process.
 *
 * This is the core auto-tuning algorithm: it emulates the relay
 * feedback loop digitally to find Ku and Pu without physically
 * running the experiment (useful for preliminary analysis).
 *
 * Algorithm:
 *   1. Apply relay feedback to process model.
 *   2. Wait for steady oscillation to develop.
 *   3. Measure oscillation period (Pu) and amplitude (a).
 *   4. Compute Ku = 4*d / (π * a).
 *
 * Complexity: O(N_steps) where N = max_cycles * (Pu / Ts).
 *
 * @param fopdt   Process model.
 * @param config  Relay configuration.
 * @param result  [out] Identified Ku, Pu.
 * @return        0 on success.
 */
int relay_simulate_fopdt(const fopdt_model_t *fopdt,
                         const relay_config_t *config,
                         relay_result_t *result);

/**
 * Complete auto-tune cycle: relay experiment → Ku/Pu → PID tuning.
 *
 * This implements the full Åström-Hägglund automatic tuning procedure:
 *   1. Run relay experiment → obtain Ku, Pu.
 *   2. From Ku, Pu → compute PID parameters using
 *      modified Z-N with setpoint weighting.
 *   3. Assign controller gains.
 *
 * @param fopdt   Process model.
 * @param config  Relay experiment configuration.
 * @param params  [out] Auto-tuned PID parameters.
 * @return        0 on success.
 */
int relay_autotune_complete(const fopdt_model_t *fopdt,
                            const relay_config_t *config,
                            pid_params_t *params);

/**
 * Adaptive relay re-tuning trigger: determines if process has changed
 * sufficiently to warrant a new auto-tune cycle.
 *
 * Monitors the difference between predicted and actual closed-loop
 * oscillation characteristics.
 *
 * @param current_Ku    Currently stored ultimate gain.
 * @param current_Pu    Currently stored ultimate period.
 * @param observed_Ku   Recently observed ultimate gain.
 * @param observed_Pu   Recently observed ultimate period.
 * @param K_threshold   Fractional change in Ku to trigger (e.g., 0.20 = 20%).
 * @param P_threshold   Fractional change in Pu to trigger.
 * @return              1 if re-tuning is recommended, 0 otherwise.
 */
int relay_needs_retune(double current_Ku, double current_Pu,
                       double observed_Ku, double observed_Pu,
                       double K_threshold, double P_threshold);

/**
 * Relay with variable hysteresis: progressively increase hysteresis
 * to reduce oscillation amplitude (safer for sensitive processes).
 *
 * Starts with ε = ε_min and doubles ε until oscillation amplitude
 * is within acceptable bounds.
 *
 * @param fopdt    Process model.
 * @param base_amplitude Base relay amplitude.
 * @param eps_min  Minimum hysteresis.
 * @param eps_max  Maximum hysteresis.
 * @param eps_step Hysteresis increment factor (> 1.0).
 * @param max_ampl Maximum allowed oscillation amplitude.
 * @param result   [out] Final Ku, Pu.
 * @return         0 on success.
 */
int relay_variable_hysteresis(const fopdt_model_t *fopdt,
                              double base_amplitude,
                              double eps_min, double eps_max, double eps_step,
                              double max_ampl, relay_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* RELAY_AUTOTUNE_H */
