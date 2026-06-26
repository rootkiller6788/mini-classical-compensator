/**
 * @file lead_frequency.h
 * @brief Frequency Domain Analysis for Lead Compensator Design
 *
 * Bode plot computation, phase/gain margin calculation, and frequency
 * response analysis tools essential for lead compensator design.
 *
 * L3 - Mathematical Structures: Complex frequency response
 * L4 - Fundamental Laws: Nyquist stability criterion on compensated loop
 * L5 - Computational Methods: Bode/Nyquist data generation
 *
 * Reference: Ogata Ch.8, Dorf&Bishop Ch.10, Franklin et al. Ch.6
 * MIT 6.302, Berkeley ME132, Cambridge 3F2
 */

#ifndef LEAD_FREQUENCY_H
#define LEAD_FREQUENCY_H

#include "lead_compensator.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L3 - Transfer Function Evaluation
 * ========================================================================= */

/** Evaluate transfer function at complex frequency s.
 *  G(s) = gain * num(s) / den(s).
 *  Uses Horner's method for polynomial evaluation.
 *  Complexity: O(order) */
lead_complex_t lead_tf_evaluate(const lead_tf_t *tf, lead_complex_t s);

/** Evaluate transfer function magnitude at s = j*omega.
 *  Complexity: O(order) */
double lead_tf_magnitude(const lead_tf_t *tf, double omega);

/** Evaluate transfer function phase at s = j*omega (radians).
 *  Complexity: O(order) */
double lead_tf_phase(const lead_tf_t *tf, double omega);

/** Evaluate plant G(s) at complex s.
 *  Complexity: O(order) */
lead_complex_t lead_system_evaluate(const lead_system_t *sys, lead_complex_t s);

/** Evaluate plant magnitude |G(jw)|.
 *  Complexity: O(order) */
double lead_system_magnitude(const lead_system_t *sys, double omega);

/** Evaluate plant phase angle(G(jw)) in radians.
 *  Complexity: O(order) */
double lead_system_phase(const lead_system_t *sys, double omega);

/* =========================================================================
 * L5 - Bode Plot Generation
 * ========================================================================= */

/** Generate Bode plot data for a transfer function.
 *  Sweeps frequency from f_min to f_max with logarithmic spacing.
 *  Complexity: O(N * order) where N is number of frequency points */
void lead_bode_compute(const lead_tf_t *tf, double f_min, double f_max,
                        int num_points, lead_bode_data_t *bode);

/** Generate Bode data for compensated open-loop C(s)*G(s).
 *  Complexity: O(N * order) */
void lead_bode_compensated(const lead_compensator_t *comp,
                            const lead_system_t *sys,
                            double f_min, double f_max, int num_points,
                            lead_bode_data_t *bode);

/** Generate Bode data for the lead compensator alone.
 *  Complexity: O(N) */
void lead_bode_compensator_only(const lead_compensator_t *comp,
                                 double f_min, double f_max, int num_points,
                                 lead_bode_data_t *bode);

/* =========================================================================
 * L4-L5 - Stability Margins (Bode/Nyquist)
 * ========================================================================= */

/** Compute gain crossover frequency: |G(jw_gc)| = 1 (0 dB).
 *  Uses bisection on magnitude in dB.
 *  Complexity: O(log(range/tol)) */
double lead_find_gain_crossover(const lead_system_t *sys, double w_min,
                                 double w_max);

/** Compute phase crossover frequency: angle(G(jw_pc)) = -180 deg.
 *  Uses bisection and angle interpolation.
 *  Complexity: O(log(range/tol)) */
double lead_find_phase_crossover(const lead_system_t *sys, double w_min,
                                  double w_max);

/** Compute phase margin from Bode plot.
 *  PM = 180 + angle(G(jw_gc)) in degrees.
 *  Complexity: O(log(range/tol)) */
double lead_compute_phase_margin(const lead_system_t *sys);

/** Compute gain margin from Bode plot.
 *  GM = -20*log10(|G(jw_pc)|) in dB.
 *  Complexity: O(log(range/tol)) */
double lead_compute_gain_margin(const lead_system_t *sys);

/** Compute phase margin of compensated system C(s)*G(s).
 *  Complexity: O(log(range/tol)) */
double lead_compensated_phase_margin(const lead_compensator_t *comp,
                                      const lead_system_t *sys);

/** Compute gain margin of compensated system C(s)*G(s).
 *  Complexity: O(log(range/tol)) */
double lead_compensated_gain_margin(const lead_compensator_t *comp,
                                     const lead_system_t *sys);

/** Classify phase margin status.
 *  Complexity: O(1) */
lead_pm_status_t lead_classify_pm(double phase_margin_deg);

/* =========================================================================
 * L4 - Nyquist Stability Criterion (Compensated Loop)
 * ========================================================================= */

/** Compute Nyquist encirclements count for compensated open loop.
 *  Evaluates C(jw)*G(jw) along the Nyquist contour, counts
 *  encirclements of the -1 point.
 *
 *  L4 Theorem: N = Z - P, where Z = N + P closed-loop RHP poles.
 *
 *  Complexity: O(N * order) */
int lead_nyquist_encirclements(const lead_compensator_t *comp,
                                const lead_system_t *sys, int num_points);

/** Check closed-loop stability via Nyquist criterion.
 *  Returns true if compensated system is stable.
 *  Complexity: O(N * order) */
bool lead_is_stable_nyquist(const lead_compensator_t *comp,
                             const lead_system_t *sys);

/** Count open-loop RHP poles of the compensated system.
 *  Uses Routh-Hurwitz criterion.
 *  Complexity: O(order^2) */
int lead_count_rhp_poles(const lead_system_t *sys);

/** Routh-Hurwitz stability test on a polynomial.
 *  L4 Theorem: Number of sign changes in first column = RHP roots.
 *  a[0] + a[1]*s + a[2]*s^2 + ... + a[n]*s^n = 0.
 *  Complexity: O(n^2) */
int lead_routh_hurwitz(const double *coeffs, int order);

/** Compute the closed-loop characteristic polynomial:
 *  1 + C(s)*G(s) = 0  =>  den_C*den_G + num_C*num_G = 0.
 *  Complexity: O(order^2) */
void lead_closed_loop_polynomial(const lead_compensator_t *comp,
                                  const lead_system_t *sys,
                                  double *cl_poly, int *cl_order);

#ifdef __cplusplus
}
#endif

#endif /* LEAD_FREQUENCY_H */
