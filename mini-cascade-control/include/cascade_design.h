/**
 * @file cascade_design.h
 * @brief Cascade controller design and synthesis API
 *
 * L4 --- Fundamental Laws: Sequential loop closure theorem.
 * L5 --- Computational Methods: Sequential tuning, optimization.
 *
 * Design Procedure (Sequential Loop Closure):
 *   1. Design inner (secondary) controller Ci(s) for fast disturbance rejection
 *   2. Close inner loop: Gi_cl(s) = Ci(s)*Gi(s) / (1 + Ci(s)*Gi(s))
 *   3. Form equivalent plant: Geq(s) = Gi_cl(s) * Go(s)
 *   4. Design outer (primary) controller Co(s) for Geq(s)
 *   5. Verify bandwidth separation: wbi >= 3 * wbo
 *
 * Course alignment:
 *   MIT 6.302 -- loop shaping for cascade
 *   Stanford ENGR207B -- sequential design
 *   Berkeley ME232 -- cascade tuning rules
 *   ETH 151-0591 -- Kaskadenregler-Entwurf
 *   Tsinghua -- Automatic Control cascade design
 */

#ifndef CASCADE_DESIGN_H
#define CASCADE_DESIGN_H

#include "cascade_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L4: Sequential Loop Closure -- Fundamental Theorem
 * ========================================================================== */

/**
 * Theorem (Sequential Loop Closure): If Ci(s) internally stabilizes Gi(s)
 * and Co(s) internally stabilizes Geq(s) = Gi_cl(s) * Go(s), then the
 * complete cascade system is internally stable provided that:
 *   (a) Geq(s) has no unstable hidden modes (minimal realization), and
 *   (b) wb_inner >= 3 * wb_outer (bandwidth separation condition).
 *
 * Reference: Zhou, Doyle & Glover, "Robust and Optimal Control" (1996),
 *            Ch. 17 on nested controller architectures.
 */

/** Design a cascade controller using sequential loop closure.
 *  @param inner_model  Transfer function of inner process Gi(s)
 *  @param outer_model  Transfer function of outer process Go(s)
 *  @param spec         Design specifications (PM, GM, BW targets)
 *  @param sys          [out] Fully designed cascade system
 *  @return 0 on success, -1 if design fails specifications */
int cascade_design_sequential(const CascadeTF *inner_model,
                               const CascadeTF *outer_model,
                               const CascadeDesignSpec *spec,
                               CascadeSystem *sys);

/** Design inner loop controller to meet bandwidth and margin targets.
 *  @param model     Inner process Gi(s)
 *  @param bw_target Desired inner closed-loop bandwidth (rad/s)
 *  @param pm_target Desired phase margin (deg)
 *  @param pid       [out] Designed PI controller
 *  @return 0 on success */
int cascade_design_inner_loop(const CascadeTF *model,
                               double bw_target,
                               double pm_target,
                               CascadePID *pid);

/** Design outer loop controller for the equivalent plant.
 *  @param eq_plant Equivalent plant Geq(s) = Gi_cl(s) * Go(s)
 *  @param spec     Design specifications
 *  @param pid      [out] Designed PID controller
 *  @return 0 on success */
int cascade_design_outer_loop(const CascadeTF *eq_plant,
                               const CascadeDesignSpec *spec,
                               CascadePID *pid);

/* ==========================================================================
 * L5: PID Tuning Rules for Cascade Control
 * ========================================================================== */

/** Inner-loop PI tuning using Skogestad SIMC rule.
 *  For Gi(s) = K*exp(-theta*s)/(tau*s+1):
 *    Kc = tau / (K * (Tc + theta))
 *    Ti = min(tau, 4*(Tc + theta))
 *  where Tc is desired closed-loop time constant.
 *
 *  Reference: Skogestad, "Simple analytic rules for model reduction
 *  and PID controller tuning" (2003), J. Process Control.
 *
 *  @param K     Process gain
 *  @param tau   Time constant (s)
 *  @param theta Time delay (s)
 *  @param Tc    Desired closed-loop time constant (s)
 *  @param pid   [out] PI controller (Kc, Ti set, Kd=0)
 *  @return 0 */
int cascade_tune_inner_skogestad(double K, double tau, double theta,
                                  double Tc, CascadePID *pid);

/** Outer-loop PID tuning using Direct Synthesis.
 *  Reference: Chen & Seborg, "PI/PID controller design based on
 *  direct synthesis" (2002), Ind. Eng. Chem. Res.
 *
 *  @param K      Equivalent plant gain
 *  @param tau    Dominant time constant (s)
 *  @param theta  Time delay (s)
 *  @param tau_cl Desired closed-loop time constant (s)
 *  @param pid    [out] PID controller
 *  @return 0 */
int cascade_tune_outer_direct_synthesis(double K, double tau, double theta,
                                         double tau_cl, CascadePID *pid);

/** Frequency-domain inner loop PI design.
 *  Solves: |Ci(jw_gc)*Gi(jw_gc)| = 1, phase = -180 + PM
 *  @param model  Inner process TF
 *  @param w_gc   Target gain crossover frequency (rad/s)
 *  @param pm     Target phase margin (deg)
 *  @param pid    [out] PI controller
 *  @return 0 on success, -1 if impossible */
int cascade_tune_inner_freq_domain(const CascadeTF *model,
                                    double w_gc, double pm,
                                    CascadePID *pid);

/* ==========================================================================
 * L5: Cascade Design Optimization
 * ========================================================================== */

/** Simultaneous cascade optimization minimizing ISE.
 *  Uses Nelder-Mead simplex over (Kc_i, Ti_i, Kc_o, Ti_o).
 *  @param inner_model Inner process model
 *  @param outer_model Outer process model
 *  @param spec        Design specifications
 *  @param sys         [out] Optimized system
 *  @return 0 on success */
int cascade_optimize_simultaneous(const CascadeTF *inner_model,
                                   const CascadeTF *outer_model,
                                   const CascadeDesignSpec *spec,
                                   CascadeSystem *sys);

/** Compute optimal bandwidth ratio for cascade control.
 *  @param inner_model Inner process model
 *  @param dist        Disturbance specification
 *  @param noise_ratio Inner/outer noise variance ratio
 *  @return Optimal bandwidth ratio (typ. 3-10) */
double cascade_optimal_bandwidth_ratio(const CascadeTF *inner_model,
                                        const CascadeDisturbance *dist,
                                        double noise_ratio);

/** Validate cascade design against specifications.
 *  @param sys   Designed system
 *  @param spec  Target specifications
 *  @return 0 if all met, bitmask of failures otherwise */
int cascade_validate_design(const CascadeSystem *sys,
                             const CascadeDesignSpec *spec);

#ifdef __cplusplus
}
#endif

#endif /* CASCADE_DESIGN_H */
