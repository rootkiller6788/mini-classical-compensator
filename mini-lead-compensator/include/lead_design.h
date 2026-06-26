/**
 * @file lead_design.h
 * @brief Lead Compensator Design Methods (L5 Computational Methods)
 *
 * Covers frequency-domain (Bode), root-locus, and analytic design.
 * Reference: Ogata Ch.7, Dorf&Bishop Ch.10, Franklin et al. Ch.6
 * MIT 6.302, Stanford ENGR105, Berkeley ME132, Caltech CDS 110
 */

#ifndef LEAD_DESIGN_H
#define LEAD_DESIGN_H

#include "lead_compensator.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L5 - Frequency Domain Design (Bode Method, Ogata Sec 7-4)
 * ========================================================================= */

/** Design lead compensator via frequency-domain (Bode) method.
 *  Steps: 1) K from e_ss, 2) PM of KG(s), 3) phi_add = PM_des-PM_exist+eps,
 *  4) alpha=(1-sin(phi_add))/(1+sin(phi_add)),
 *  5) find w where |KG(jw)|=-10*log10(alpha) dB,
 *  6) T=1/(w*sqrt(alpha)), 7) K_c=K/alpha, 8) verify PM.
 *  Complexity: O(N * order), Course: MIT 6.302 */
bool lead_design_frequency(const lead_system_t *system,
                            const lead_specs_t *specs,
                            lead_design_result_t *result);

/** Compute required gain K from steady-state error specification.
 *  Type-0: K from position error constant.
 *  Type-1: K from velocity error constant.
 *  Complexity: O(1) */
double lead_design_gain_from_error(const lead_system_t *system, double e_ss);

/** Find frequency where |KG(jw)| in dB crosses -10*log10(alpha).
 *  Uses bisection. This becomes the new gain crossover frequency.
 *  Complexity: O(log(range/tol)) */
double lead_find_new_crossover(const lead_system_t *system, double alpha,
                                double w_min, double w_max);

/* =========================================================================
 * L5 - Root Locus Design
 * ========================================================================= */

/** Design lead compensator via root-locus method.
 *  Uses Evans angle condition at desired pole to compute angle deficiency,
 *  then places lead zero/pole to supply the needed angle.
 *  Complexity: O(order), Course: Stanford ENGR105 */
bool lead_design_root_locus(const lead_system_t *system,
                             lead_complex_t s_desired,
                             lead_design_result_t *result);

/** Compute angle deficiency at point s: what the compensator must supply
 *  to make s part of the compensated root locus.
 *  Complexity: O(order) */
double lead_angle_deficiency(const lead_system_t *system, lead_complex_t s);

/** Place lead zero and pole to provide angle phi_c at point s.
 *  Uses bisector method for geometric placement.
 *  Complexity: O(1) */
bool lead_place_zero_pole_angle(lead_complex_t s, double phi_c,
                                 double *z_c_out, double *p_c_out);

/* =========================================================================
 * L5 - Analytic Design
 * ========================================================================= */

/** Analytic design from time-domain specs (PO%, ts -> zeta, wn -> PM, w_gc).
 *  Complexity: O(log range) */
bool lead_design_analytic(const lead_specs_t *specs,
                           const lead_system_t *system,
                           lead_design_result_t *result);

/** Compute dominant closed-loop pole from (zeta, wn).
 *  s_d = -zeta*wn +/- j*wn*sqrt(1-zeta^2).
 *  Complexity: O(1) */
lead_complex_t lead_dominant_pole(double zeta, double wn, bool upper_half);

/* =========================================================================
 * L2 - Performance Prediction
 * ========================================================================= */

/** Predict step response metrics from phase margin and crossover.
 *  Uses 2nd-order approximation. Complexity: O(1) */
void lead_predict_performance(double phase_margin, double crossover_freq,
                               lead_performance_t *perf);

/** Verify compensator design by computing closed-loop metrics.
 *  Complexity: O(order) */
void lead_verify_design(const lead_compensator_t *compensator,
                         const lead_system_t *system,
                         lead_performance_t *perf);

/** Check if design meets all specs. Complexity: O(1) */
bool lead_check_specs(const lead_performance_t *perf,
                       const lead_specs_t *specs);

/** PM-zeta relationship (accuracy improved over linear rule).
 *  PM = atan(2*zeta/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)).
 *  Complexity: O(1) */
double lead_zeta_to_pm(double zeta);

/** Required bandwidth from (zeta, settling time).
 *  w_BW = 4/(zeta*ts)*sqrt(sqrt(1+4*zeta^4)-2*zeta^2).
 *  Complexity: O(1) */
double lead_bandwidth_from_zeta_ts(double zeta, double ts);

/** Approximate crossover from bandwidth: w_gc ~ w_BW / 1.6.
 *  Complexity: O(1) */
double lead_crossover_from_bandwidth(double bandwidth);

#ifdef __cplusplus
}
#endif

#endif /* LEAD_DESIGN_H */
