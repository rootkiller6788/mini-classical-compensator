/**
 * @file lag_frequency.h
 * @brief Frequency-domain analysis: Bode, Nyquist, stability margins
 *
 * Tools for analyzing the frequency response of lag-compensated systems.
 * Central to the classical control design methodology.
 *
 * L4: Nyquist stability criterion, gain/phase margin computation
 * L5: Bode plot construction, frequency response evaluation
 *
 * Course alignment:
 *   MIT 6.302 — frequency response methods for compensator analysis
 *   Stanford ENGR105 — Bode/Nyquist analysis for controller design
 *   Berkeley ME232 — frequency-domain loop analysis
 *   Caltech CDS 110 — stability margins in frequency domain
 *   ETH 151-0591 — Frequenzgangdarstellung
 *   Cambridge 3F2 — frequency response methods
 *   Tsinghua 自动控制原理 — 频域分析
 *
 * Textbook: Ogata Ch. 8; Franklin/Powell/Emami-Naeini Ch. 6
 */

#ifndef LAG_FREQUENCY_H
#define LAG_FREQUENCY_H

#include "lag_types.h"
#include "lag_compensator.h"

/* ==========================================================================
 * L5: Bode plot computation
 * ========================================================================== */

/**
 * Compute the complete Bode plot for a transfer function.
 *
 * Generates magnitude (dB) and phase (degrees) over a logarithmic
 * frequency range. The frequency range is automatically determined
 * to span at least 2 decades below and above the corner frequencies.
 *
 * Algorithm:
 *   1. Determine frequency range from TF poles/zeros
 *   2. Generate logarithmically spaced frequency points
 *   3. Evaluate magnitude and phase at each point
 *   4. Detect gain and phase crossover frequencies
 *   5. Compute stability margins
 *
 * Complexity: O(n_points) where n_points is number of frequency samples
 *
 * @param tf          transfer function to analyze
 * @param n_points    number of frequency points (suggest 500-2000)
 * @param[out] bode   allocated BodeData struct (caller frees points)
 * @return 0 on success, negative on failure
 */
int lag_compute_bode(const LagTransferFunction *tf, int n_points,
                     LagBodeData *bode);

/**
 * Compute Bode plot for the open-loop transfer function L(s) = G_c(s)*G(s).
 *
 * This is the key analysis for stability margin verification.
 *
 * Complexity: O(n_points)
 */
int lag_compute_open_loop_bode(const LagCompensator *lag,
                                const LagTransferFunction *plant,
                                int n_points, LagBodeData *bode);

/**
 * Compute Bode plot for the closed-loop transfer function T(s) = L/(1+L).
 *
 * Complexity: O(n_points)
 */
int lag_compute_closed_loop_bode(const LagCompensator *lag,
                                  const LagTransferFunction *plant,
                                  int n_points, LagBodeData *bode);

/**
 * Compute Bode plot for the sensitivity function S(s) = 1/(1+L).
 *
 * Sensitivity analysis is important for disturbance rejection assessment.
 *
 * Complexity: O(n_points)
 */
int lag_compute_sensitivity_bode(const LagCompensator *lag,
                                  const LagTransferFunction *plant,
                                  int n_points, LagBodeData *bode);

/* ==========================================================================
 * L5: Nyquist plot computation
 * ========================================================================== */

/**
 * Compute Nyquist plot data for a transfer function.
 *
 * Evaluates G(jw) for w in [0, inf) and its conjugate for w in (-inf, 0],
 * producing a closed contour in the complex plane.
 *
 * L4 Theorem (Nyquist Stability Criterion):
 *   Z = N + P
 *   where Z = # of unstable closed-loop poles
 *         N = # of clockwise encirclements of -1
 *         P = # of unstable open-loop poles
 *
 * Algorithm:
 *   1. Evaluate G(jw) for positive frequencies
 *   2. Mirror for negative frequencies (conjugate)
 *   3. Count encirclements of the -1 point using winding number
 *   4. Determine closed-loop stability from N + P
 *
 * Complexity: O(n_points)
 *
 * @param tf           open-loop transfer function
 * @param n_points     number of frequency points
 * @param P            number of unstable open-loop poles
 * @param[out] nyquist allocated NyquistData struct
 * @return 0 on success
 */
int lag_compute_nyquist(const LagTransferFunction *tf, int n_points,
                        int P, LagNyquistData *nyquist);

/* ==========================================================================
 * L4: Stability margin computation
 * ========================================================================== */

/**
 * Find the gain crossover frequency w_gc where |G(jw)| = 1 (0 dB).
 *
 * Uses bisection search on the Bode magnitude data.
 * If no crossover exists (gain always < 1), returns -1.
 * If multiple crossovers exist, returns the highest one.
 *
 * Complexity: O(n_points) for one pass through frequency data
 *
 * @param bode  Bode plot data for the system
 * @return gain crossover frequency in rad/s, or -1 if none
 */
double lag_find_gain_crossover(const LagBodeData *bode);

/**
 * Find the phase crossover frequency w_pc where angle(G(jw)) = -180 deg.
 *
 * Uses linear interpolation between frequency points.
 * If no crossover exists, returns -1.
 *
 * Complexity: O(n_points)
 */
double lag_find_phase_crossover(const LagBodeData *bode);

/**
 * Compute the phase margin at the gain crossover frequency.
 *
 * Definition: PM = 180 + angle(L(j*w_gc))  [in degrees]
 *
 * Theorem (Stability): For minimum-phase systems,
 * closed-loop stability requires PM > 0.
 *
 * Complexity: O(n_points) (calls lag_find_gain_crossover)
 *
 * @param bode  Bode plot data for open-loop transfer function
 * @return phase margin in degrees
 */
double lag_compute_phase_margin(const LagBodeData *bode);

/**
 * Compute the gain margin at the phase crossover frequency.
 *
 * Definition: GM = -20*log10(|L(j*w_pc)|)  [in dB]
 *
 * Theorem (Stability): For minimum-phase systems,
 * closed-loop stability requires GM > 0 dB.
 *
 * Complexity: O(n_points)
 *
 * @param bode  Bode plot data for open-loop transfer function
 * @return gain margin in dB
 */
double lag_compute_gain_margin(const LagBodeData *bode);

/**
 * Compute both phase and gain margins in one call.
 *
 * Complexity: O(n_points)
 */
void lag_compute_stability_margins(const LagBodeData *bode,
                                    double *phase_margin_deg,
                                    double *gain_margin_db);

/* ==========================================================================
 * L2: Bandwidth computation
 * ========================================================================== */

/**
 * Compute the closed-loop bandwidth.
 *
 * Definition: frequency w_bw where |T(jw)| = -3 dB (0.707)
 * T(s) = L(s) / (1 + L(s)) is the complementary sensitivity.
 *
 * The bandwidth indicates the frequency range over which
 * the system can track reference signals effectively.
 *
 * A lag compensator typically reduces bandwidth because
 * it attenuates high frequencies by factor beta.
 *
 * Complexity: O(n_points)
 *
 * @param bode  Bode plot data for closed-loop transfer function
 * @return bandwidth in rad/s, or -1 if not found
 */
double lag_compute_bandwidth(const LagBodeData *bode);

/* ==========================================================================
 * L4: Steady-state error from frequency response
 * ========================================================================== */

/**
 * Compute steady-state error from the low-frequency behavior.
 *
 * Using the Final Value Theorem:
 *   Kp = lim_{w->0} |L(jw)|                    (position error constant)
 *   Kv = lim_{w->0} w * |L(jw)|                (velocity error constant)
 *   Ka = lim_{w->0} w^2 * |L(jw)|              (acceleration error constant)
 *
 * Then:
 *   e_ss_step = 1 / (1 + Kp)
 *   e_ss_ramp = 1 / Kv
 *   e_ss_parabolic = 1 / Ka
 *
 * Complexity: O(1) (uses lowest-frequency data point)
 *
 * @param bode     Bode plot data for open-loop system
 * @param ess_type type of input for error computation
 * @return steady-state error (dimensionless or in appropriate units)
 */
double lag_compute_steady_state_error(const LagBodeData *bode,
                                       LagESSType ess_type);

/**
 * Estimate the error constant directly from the low-frequency asymptote.
 *
 * Extrapolates the -20*N dB/dec slope back to w = 1 rad/s
 * to estimate Kv, Ka, etc.
 *
 * Complexity: O(1)
 */
double lag_estimate_error_constant(const LagBodeData *bode,
                                    LagESSType ess_type);

/* ==========================================================================
 * L4: Frequency-domain performance metrics
 * ========================================================================== */

/**
 * Compute the resonant peak magnitude M_p from Bode data.
 *
 * M_p = max_{w} |T(jw)|
 *
 * Related to damping ratio: M_p = 1/(2*zeta*sqrt(1-zeta^2)) for 2nd order
 *
 * Complexity: O(n_points)
 */
double lag_compute_resonant_peak(const LagBodeData *closed_loop_bode);

/**
 * Compute the resonant frequency w_r where |T(jw)| is maximum.
 *
 * Complexity: O(n_points)
 */
double lag_compute_resonant_frequency(const LagBodeData *closed_loop_bode);

#endif /* LAG_FREQUENCY_H */