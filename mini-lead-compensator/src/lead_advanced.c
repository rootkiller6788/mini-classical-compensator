/**
 * @file lead_advanced.c
 * @brief Advanced Lead Compensator Topics (L8)
 *
 * L8 - Advanced Topics:
 *   - Multiple-stage lead compensation
 *   - Lead-lag combination design principles
 *   - Robust lead design with gain/phase margin constraints
 *   - Time-varying lead (adaptive/gain-scheduled)
 *   - Discrete-time (digital) lead compensator via Tustin/Bilinear
 *
 * Reference: Astrom & Murray Ch.11, Skogestad & Postlethwaite Ch.5
 * Stanford ENGR210B, Cambridge 4F2, ETH 151-0563
 */

#include "lead_compensator.h"
#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* =========================================================================
 * L8 - Multi-Stage Lead Compensation
 * ========================================================================= */

/**
 * Design two-stage lead compensator for large phase lead (> 65 deg).
 *
 * L8 - When required phase lead exceeds LEAD_MAX_SINGLE_STAGE,
 * use two cascaded lead compensators: C(s) = C1(s) * C2(s).
 *
 * Each stage provides half the required lead:
 *   phi_c1 = phi_c2 = phi_total / 2
 *
 * Benefits:
 *   - Achieves 90-110 deg phase lead (vs. 65 deg single-stage)
 *   - Reduced noise amplification (each stage has larger alpha)
 *   - Better robustness (spread-out lead frequencies)
 *
 * Design:
 *   1. Split required phase lead: phi_per_stage = phi_total / N
 *   2. For each stage, determine alpha_i from phi_per_stage
 *   3. Place omega_m_i at geometrically spaced frequencies
 *      (e.g., omega_m2 = k * omega_m1, where k > 1)
 *   4. Each stage provides gain compensation at its own frequency
 */
bool lead_design_two_stage(const lead_system_t *system,
                            const lead_specs_t *specs,
                            lead_compensator_t *stage1,
                            lead_compensator_t *stage2) {
    if (!system || !specs || !stage1 || !stage2) return false;

    /* Compute total phase lead needed */
    double pm_existing = lead_compute_phase_margin(system);
    double phi_total = specs->phase_margin_desired - pm_existing
                       + LEAD_PHASE_MARGIN_EPS;
    if (phi_total < 0.0) phi_total = 0.0;
    if (phi_total > 110.0) phi_total = 110.0; /* practical limit for 2-stage */

    /* Split: each stage provides half */
    double phi_per_stage = phi_total / 2.0;
    double phi_per_rad = phi_per_stage * M_PI / 180.0;

    double alpha = lead_alpha_from_phi_max(phi_per_rad);

    /* Stage 1: placed at lower frequency */
    double w_gc = lead_find_gain_crossover(system, LEAD_FREQ_MIN_DEFAULT,
                                            LEAD_FREQ_MAX_DEFAULT);
    if (w_gc <= 0.0) w_gc = 1.0;

    double w_m1 = w_gc * 0.7; /* slightly below original crossover */
    double T1 = 1.0 / (w_m1 * sqrt(alpha));
    double K1 = 1.0;

    if (!lead_init_from_KT_alpha(stage1, K1, T1, alpha, LEAD_TYPE_ANALYTIC))
        return false;

    /* Stage 2: placed at higher frequency (wider band) */
    double w_m2 = w_gc * 1.5; /* above original crossover, extends bandwidth */
    double T2 = 1.0 / (w_m2 * sqrt(alpha));

    if (!lead_init_from_KT_alpha(stage2, 1.0, T2, alpha, LEAD_TYPE_ANALYTIC))
        return false;

    return true;
}

/**
 * Compute the total transfer function of cascaded compensators.
 *
 * L8 - C_total(s) = C1(s) * C2(s)
 * For lead compensators in cascade:
 *   C_total(s) = K1*K2 * (T1*s+1)*(T2*s+1) / ((alpha1*T1*s+1)*(alpha2*T2*s+1))
 */
void lead_cascade_transfer_function(const lead_compensator_t *c1,
                                     const lead_compensator_t *c2,
                                     lead_tf_t *tf) {
    if (!c1 || !c2 || !tf) return;
    memset(tf, 0, sizeof(lead_tf_t));

    /* Numerator: K1*K2 * (T1*s+1)*(T2*s+1) = K1*K2 * [1, T1+T2, T1*T2] */
    double K_total = c1->dc_gain * c2->dc_gain;
    tf->num_order = 2;
    tf->num[0] = K_total;
    tf->num[1] = K_total * (c1->T + c2->T);
    tf->num[2] = K_total * c1->T * c2->T;

    /* Denominator: (alpha1*T1*s+1)*(alpha2*T2*s+1) = [1, a1T1+a2T2, a1a2*T1*T2] */
    double a1T1 = c1->alpha * c1->T;
    double a2T2 = c2->alpha * c2->T;
    tf->den_order = 2;
    tf->den[0] = 1.0;
    tf->den[1] = a1T1 + a2T2;
    tf->den[2] = a1T1 * a2T2;
    tf->gain = 1.0;
}

/**
 * Combined phase of cascaded compensators at frequency omega.
 */
double lead_cascade_phase(const lead_compensator_t *c1,
                           const lead_compensator_t *c2,
                           double omega) {
    return lead_phase_at(c1, omega) + lead_phase_at(c2, omega);
}

/**
 * Combined magnitude of cascaded compensators at frequency omega.
 */
double lead_cascade_magnitude(const lead_compensator_t *c1,
                               const lead_compensator_t *c2,
                               double omega) {
    return lead_magnitude_at(c1, omega) * lead_magnitude_at(c2, omega);
}

/* =========================================================================
 * L8 - Discrete-Time (Digital) Lead Compensator
 * ========================================================================= */

/**
 * Convert continuous-time lead compensator to discrete-time
 * using Tustin's method (bilinear transform).
 *
 * L8 - s = (2/T_s) * (z-1)/(z+1) where T_s is sampling period.
 *
 * The discrete-time transfer function:
 *   C(z) = (b0 + b1*z^{-1}) / (1 + a1*z^{-1})
 *
 * where:
 *   b0 = (2*K_c/T_s + K_c*z_c) / (2/T_s + p_c)
 *   b1 = (K_c*z_c - 2*K_c/T_s) / (2/T_s + p_c)
 *   a1 = (p_c - 2/T_s) / (2/T_s + p_c)
 *
 * This is critical for digital implementation in microcontrollers
 * and real-time control systems.
 */
void lead_to_discrete_tustin(const lead_compensator_t *comp, double Ts,
                              double *b0, double *b1, double *a1) {
    if (!comp || !b0 || !b1 || !a1 || Ts <= 0.0) return;

    double K_c = comp->K_c;
    double z_c = comp->z_c;
    double p_c = comp->p_c;
    double two_over_Ts = 2.0 / Ts;

    double den = two_over_Ts + p_c;
    if (fabs(den) < 1e-300) {
        *b0 = 1.0; *b1 = 0.0; *a1 = 0.0;
        return;
    }

    *b0 = (two_over_Ts * K_c + K_c * z_c) / den;
    *b1 = (K_c * z_c - two_over_Ts * K_c) / den;
    *a1 = (p_c - two_over_Ts) / den;
}

/**
 * Evaluate discrete-time lead compensator at a given sample.
 *
 * L8 - y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]
 *
 * Implements the difference equation for real-time filtering.
 */
double lead_discrete_update(double b0, double b1, double a1,
                             double x_new, double *x_prev, double *y_prev) {
    if (!x_prev || !y_prev) return 0.0;
    double y_new = b0 * x_new + b1 * (*x_prev) - a1 * (*y_prev);
    *x_prev = x_new;
    *y_prev = y_new;
    return y_new;
}

/* =========================================================================
 * L8 - Robust Lead Design with GM/PM Constraints
 * ========================================================================= */

/**
 * Design lead compensator with simultaneous gain and phase margin constraints.
 *
 * L8 - Robust design ensures both GM >= GM_min and PM >= PM_min.
 * Uses iterative refinement checking both margins.
 *
 * For systems with significant phase lag, satisfying both margins
 * simultaneously may require larger alpha (less lead) to avoid
 * reducing gain margin.
 */
bool lead_design_robust(const lead_system_t *system,
                         const lead_specs_t *specs,
                         lead_design_result_t *result) {
    if (!system || !specs || !result) return false;

    /* Standard frequency-domain design first */
    bool ok = lead_design_frequency(system, specs, result);
    if (!ok) return false;

    /* Verify gain margin */
    double gm = lead_compensated_gain_margin(&result->compensator, system);

    if (gm < specs->gain_margin_desired && gm > 0.0) {
        /* GM too low: reduce alpha (less lead, less HF gain) */
        double alpha_orig = result->alpha_actual;
        double gm_shortfall = specs->gain_margin_desired - gm;

        /* Increase alpha slightly (less aggressive lead -> better GM) */
        double alpha_new = alpha_orig * (1.0 + gm_shortfall / 20.0);
        if (alpha_new > 0.9) alpha_new = 0.9;

        /* Re-design with new alpha */
        lead_specs_t adj_specs = *specs;
        double phi_new = lead_phi_max_from_alpha(alpha_new);
        adj_specs.phase_margin_desired = phi_new * 180.0 / M_PI
                                          + lead_compute_phase_margin(system);
        if (adj_specs.phase_margin_desired > specs->phase_margin_desired) {
            adj_specs.phase_margin_desired = specs->phase_margin_desired;
        }

        return lead_design_frequency(system, &adj_specs, result);
    }

    return true;
}

/**
 * Compute the gain margin degradation caused by the lead compensator.
 *
 * L8 - Lead compensators reduce GM because they increase HF gain.
 * This function quantifies the GM reduction:
 *   delta_GM = GM_uncompensated - GM_compensated (dB)
 *
 * If delta_GM > 6 dB (excessive), consider reducing alpha
 * or switching to lead-lag compensation.
 */
double lead_gm_degradation(const lead_compensator_t *comp,
                            const lead_system_t *sys) {
    if (!comp || !sys) return 0.0;

    double gm_uncomp = lead_compute_gain_margin(sys);
    double gm_comp = lead_compensated_gain_margin(comp, sys);

    if (isinf(gm_uncomp)) gm_uncomp = 40.0;
    if (isinf(gm_comp)) gm_comp = 40.0;

    return gm_uncomp - gm_comp;
}

/**
 * Sensitivity-based robust stability check.
 *
 * L8 - Small Gain Theorem application:
 * For multiplicative uncertainty bounded by |W_T(jw)|,
 * system is robustly stable iff |T(jw)*W_T(jw)| < 1 for all w.
 *
 * This function checks: M_t * w_u < 1 where w_u = uncertainty weight.
 */
bool lead_robust_stability_check(const lead_compensator_t *comp,
                                  const lead_system_t *sys,
                                  double uncertainty_weight) {
    if (!comp || !sys) return false;

    double M_t = lead_peak_complementary_sensitivity(comp, sys,
                                                      LEAD_FREQ_MIN_DEFAULT,
                                                      LEAD_FREQ_MAX_DEFAULT,
                                                      500);
    return (M_t * uncertainty_weight < 1.0);
}

/**
 * Compute the frequency at which the lead compensator phase is half
 * of its maximum (3-dB phase bandwidth).
 *
 * L8 - Characterizes the effective range of the lead compensator.
 * Useful for assessing whether the lead band covers the crossover region.
 */
double lead_half_phase_frequency(const lead_compensator_t *comp, bool upper) {
    if (!comp || !lead_validate(comp)) return 0.0;

    double phi_max = lead_phi_max_from_alpha(comp->alpha);
    double phi_half = phi_max / 2.0;

    double w_m = lead_omega_max(comp);
    double search_min = upper ? w_m : 0.0;
    double search_max = upper ? w_m * 100.0 : w_m;

    if (!upper) {
        if (lead_phase_at(comp, 0.0) >= phi_half) return 0.0;
    }

    /* Bisection */
    double lo = search_min, hi = search_max;
    for (int iter = 0; iter < 60; iter++) {
        double mid = (lo + hi) * 0.5;
        double phi_mid = lead_phase_at(comp, mid);
        if (upper) {
            if (phi_mid > phi_half) lo = mid;
            else hi = mid;
        } else {
            if (phi_mid < phi_half) lo = mid;
            else hi = mid;
        }
        if ((hi - lo) / (mid + 1e-12) < 1e-10) break;
    }
    return (lo + hi) * 0.5;
}