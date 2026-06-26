/**
 * @file lag_design.h
 * @brief Lag compensator design algorithms and procedures
 *
 * Systematic design methods for synthesizing a lag compensator to meet
 * steady-state error specifications while preserving stability margins.
 *
 * Design procedure (Ogata, Ch. 7):
 *   1. Determine the required error constant from ESS specification
 *   2. Compute beta = required_Kv / current_Kv
 *   3. Determine gain crossover frequency of uncompensated system
 *   4. Place corner frequency 1/T at w_gc/10 (one decade below)
 *   5. Compute T from corner frequency placement
 *   6. Verify phase margin of compensated system
 *
 * L5: All design algorithms are computational implementations of
 *     engineering design procedures documented in control textbooks.
 *
 * Course alignment:
 *   MIT 6.302 / Stanford ENGR105 / Berkeley ME232 / Caltech CDS 110
 *   ETH 151-0591 / Cambridge 3F2 / Tsinghua 自动控制原理
 */

#ifndef LAG_DESIGN_H
#define LAG_DESIGN_H

#include "lag_types.h"
#include "lag_compensator.h"

/* ==========================================================================
 * L5: Steady-state error based design
 * ========================================================================== */

/**
 * Design a lag compensator to satisfy a steady-state error requirement.
 *
 * Theorem (Final Value Theorem): For unity feedback with open-loop TF
 * G_ol(s), the steady-state error to a step input is:
 *   e_ss = 1 / (1 + Kp)  where Kp = lim_{s->0} G_ol(s)
 * For ramp input:  e_ss = 1/Kv where Kv = lim_{s->0} s*G_ol(s)
 * For parabolic:   e_ss = 1/Ka where Ka = lim_{s->0} s^2*G_ol(s)
 *
 * The lag compensator multiplies the error constant by Kc.
 *
 * Algorithm: compute required Kc, set beta=Kc, place corner freq,
 * verify design.
 *
 * Complexity: O(1)
 *
 * @param plant       plant transfer function G(s)
 * @param spec        design specification (ESS target, type)
 * @param[out] result designed lag compensator
 * @return 0 on success, negative on failure
 */
int lag_design_for_steady_state_error(const LagTransferFunction *plant,
                                       const LagDesignSpec *spec,
                                       LagCompensator *result);

/**
 * Design a lag compensator from error constant improvement factor.
 *
 * Given current error constant K_cur and required K_req,
 * compute beta = K_req / K_cur. Construct:
 *   G_c(s) = beta * (T*s + 1) / (beta*T*s + 1)
 * with T chosen to place zero corner one decade below w_gc.
 *
 * Complexity: O(1)
 *
 * @param K_current     current error constant of uncompensated system
 * @param K_required    required error constant
 * @param w_gc_current  current gain crossover frequency (rad/s)
 * @param[out] result   designed lag compensator
 * @return 0 on success, negative on invalid parameters
 */
int lag_design_from_error_constants(double K_current, double K_required,
                                     double w_gc_current,
                                     LagCompensator *result);

/* ==========================================================================
 * L5: Phase-margin-preserving design
 * ========================================================================== */

/**
 * Design a lag compensator that preserves a specified phase margin.
 *
 * The phase lag introduced by the compensator at the gain crossover
 * frequency must not reduce the phase margin below target.
 * This is achieved by placing both corner frequencies well below w_gc.
 *
 * Rule of thumb: 1/T <= w_gc / 10 (zero corner at least 1 decade below w_gc)
 *
 * The phase contribution at w_gc is approximately:
 *   phi_c(w_gc) = arctan(w_gc*T) - arctan(w_gc*beta*T)
 *
 * When 1/T = w_gc/10: phi_c ~= 5.71 - arctan(0.1*beta) degrees
 * This is small enough to preserve phase margin for moderate beta.
 *
 * Algorithm:
 *   1. Measure uncompensated phase margin
 *   2. Determine allowable phase reduction
 *   3. Choose T such that phase contribution at w_gc < allowable reduction
 *   4. Set beta from ESS requirement
 *   5. Verify and iterate if needed
 *
 * Complexity: O(log(beta)) due to frequency search
 */
int lag_design_for_phase_margin(const LagTransferFunction *plant,
                                 const LagDesignSpec *spec,
                                 LagCompensator *result);

/* ==========================================================================
 * L5: Bandwidth-constrained design
 * ========================================================================== */

/**
 * Design a lag compensator with bandwidth constraints.
 *
 * The lag compensator reduces high-frequency gain by factor beta,
 * which may lower the gain crossover frequency. This method selects
 * beta and T to keep the bandwidth within specified limits.
 *
 * Complexity: O(log(w_max/w_min))
 */
int lag_design_for_bandwidth(const LagTransferFunction *plant,
                              const LagDesignSpec *spec,
                              LagCompensator *result);

/* ==========================================================================
 * L5: Root-locus-based design
 * ========================================================================== */

/**
 * Design a lag compensator using root locus considerations.
 *
 * The lag compensator adds a dipole (pole-zero pair close together)
 * near the origin. This increases the error constant without
 * significantly changing the root locus away from the origin.
 *
 * The ratio beta = z_c/p_c determines error constant improvement.
 * The dipole must be placed close to the origin to minimize
 * its effect on the transient response.
 *
 * Design steps:
 *   1. Find dominant closed-loop poles from transient spec
 *   2. Compute required error constant improvement (beta)
 *   3. Place zero at |z_c| = |dominant_pole| / 100 (near origin)
 *   4. Place pole at |p_c| = |z_c| / beta (closer to origin)
 *   5. Adjust Kc for root locus preservation
 *
 * Complexity: O(1)
 */
int lag_design_root_locus(const LagTransferFunction *plant,
                           double dominant_pole, double beta_required,
                           LagCompensator *result);

/* ==========================================================================
 * L5: Bode method — full iterative frequency-domain design
 * ========================================================================== */

/**
 * Complete frequency-domain lag compensator design (Bode method).
 *
 * Standard textbook procedure (Ogata Sec. 7-5):
 *
 * Step 1: Determine open-loop gain K to satisfy ESS requirement.
 * Step 2: Draw Bode plot of uncompensated system with gain K.
 * Step 3: Find w_gc and phase margin of uncompensated system.
 * Step 4: Determine required phase addition (allow 5-12 deg safety).
 * Step 5: Find frequency where phase = -180 + PM_target + safety.
 *         This becomes the new target w_gc.
 * Step 6: Place zero corner 1/T at w_gc/10.
 * Step 7: Determine beta from attenuation needed at w_gc.
 * Step 8: Compute pole corner = 1/(beta*T).
 * Step 9: Verify design, iterate if necessary.
 *
 * Complexity: O(n_freq) where n_freq ~= 1000 frequency points
 *
 * @param plant     plant transfer function G(s)
 * @param spec      complete design specification
 * @param[out] result  designed lag compensator
 * @return 0 on success, negative on failure
 */
int lag_design_bode_method(const LagTransferFunction *plant,
                            const LagDesignSpec *spec,
                            LagCompensator *result);

/* ==========================================================================
 * L5: Design verification and helpers
 * ========================================================================== */

/**
 * Verify that a designed lag compensator meets all specifications.
 *
 * Checks: phase margin, gain margin, steady-state error, bandwidth.
 *
 * @return bitmask: bit0=PM, bit1=GM, bit2=ESS, bit3=bandwidth
 *         (15 = all passed)
 *
 * Complexity: O(n_freq) for frequency response sweep
 */
int lag_verify_design(const LagCompensator *lag,
                       const LagTransferFunction *plant,
                       const LagDesignSpec *spec);

/**
 * Compute recommended safety margin (extra phase) based on beta.
 *
 * Heuristic: safety_margin_deg = 5 + 10*log10(beta)
 * Clamped to [5, 15] degrees for typical beta in [1, 10].
 *
 * Complexity: O(1)
 */
double lag_recommended_safety_margin(double beta);

/**
 * Initialize a design specification with reasonable defaults.
 * PM=45deg, GM=10dB, ESS=2%, step input, one decade safety margin.
 * Complexity: O(1)
 */
LagDesignSpec lag_design_spec_default(void);

/** Set ESS specification. Complexity: O(1) */
void lag_design_spec_set_ess(LagDesignSpec *spec,
                              LagESSType type, double ess_target);

/** Set PM/GM specification. Complexity: O(1) */
void lag_design_spec_set_pm(LagDesignSpec *spec,
                             double pm_target, double gm_target);

#endif /* LAG_DESIGN_H */