/**
 * @file lead_design.c
 * @brief Lead Compensator Design Methods Implementation
 *
 * L5 - Computational Methods:
 *   - Frequency-domain (Bode) design (Ogata Sec 7-4)
 *   - Root-locus design (Evans angle condition)
 *   - Analytic design from time-domain specs
 *
 * Reference: Ogata Ch.7, Dorf&Bishop Ch.10, Franklin et al. Ch.6
 * MIT 6.302, Stanford ENGR105, Berkeley ME132, Caltech CDS 110
 */

#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* =========================================================================
 * L5 - Frequency-Domain (Bode) Design
 * ========================================================================= */

/**
 * Design lead compensator via Bode method (Ogata 7-4).
 *
 * Algorithm:
 * 1. Determine K from steady-state error requirement
 * 2. Compute PM of KG(s)
 * 3. phi_add = PM_desired - PM_existing + epsilon (5-12 deg)
 * 4. Cap phi_add at LEAD_MAX_SINGLE_STAGE (65 deg)
 * 5. alpha = (1 - sin(phi_add)) / (1 + sin(phi_add))
 * 6. Find w_new where |KG(jw)|_dB = -10*log10(alpha)
 * 7. T = 1/(w_new * sqrt(alpha)), K_c = K/alpha
 * 8. Verify PM of C(s)G(s), iterate if needed
 */
bool lead_design_frequency(const lead_system_t *system,
                            const lead_specs_t *specs,
                            lead_design_result_t *result) {
    if (!system || !specs || !result) return false;

    memset(result, 0, sizeof(lead_design_result_t));
    result->specs = *specs;
    result->converged = false;

    /* Step 1: Determine gain K from steady-state error */
    double K = lead_design_gain_from_error(system, specs->steady_state_error);
    lead_system_t sys_scaled = *system;
    sys_scaled.tf.gain *= K;

    /* Step 2: Compute existing PM */
    double pm_existing = lead_compute_phase_margin(&sys_scaled);
    double w_gc_old = lead_find_gain_crossover(&sys_scaled,
                                                LEAD_FREQ_MIN_DEFAULT,
                                                LEAD_FREQ_MAX_DEFAULT);
    result->crossover_freq_old = w_gc_old;

    /* Step 3: Additional phase lead needed */
    double epsilon = LEAD_PHASE_MARGIN_EPS;
    double phi_m_add_deg = specs->phase_margin_desired - pm_existing + epsilon;
    if (phi_m_add_deg < 0.0) phi_m_add_deg = 0.0;
    if (phi_m_add_deg > LEAD_MAX_SINGLE_STAGE)
        phi_m_add_deg = LEAD_MAX_SINGLE_STAGE;
    double phi_m_add = phi_m_add_deg * M_PI / 180.0;

    /* Step 4: Compute alpha */
    double alpha = lead_alpha_from_phi_max(phi_m_add);
    if (alpha < LEAD_ALPHA_MIN) alpha = LEAD_ALPHA_MIN;

    /* Step 5: Find new crossover frequency */
    double w_new = lead_find_new_crossover(&sys_scaled, alpha,
                                            LEAD_FREQ_MIN_DEFAULT,
                                            LEAD_FREQ_MAX_DEFAULT);
    if (w_new <= 0.0) w_new = w_gc_old;

    /* Step 6-7: Compute T and K_c */
    double T = 1.0 / (w_new * sqrt(alpha));
    double K_c = K / alpha;

    /* Step 8: Create compensator */
    lead_compensator_t comp;
    if (!lead_init(&comp, K_c, 1.0/T, 1.0/(alpha*T), LEAD_TYPE_ANALYTIC))
        return false;

    /* Step 9: Iterative refinement */
    double pm_achieved = 0.0;
    bool converged = false;
    int iter = 0;

    for (iter = 0; iter < LEAD_MAX_ITERATIONS; iter++) {
        pm_achieved = lead_compensated_phase_margin(&comp, system);
        if (pm_achieved >= specs->phase_margin_desired - 0.1) {
            converged = true;
            break;
        }
        double shortfall = specs->phase_margin_desired - pm_achieved;
        if (shortfall < 0.5) break;

        phi_m_add_deg += shortfall * 0.5;
        if (phi_m_add_deg > LEAD_MAX_SINGLE_STAGE)
            phi_m_add_deg = LEAD_MAX_SINGLE_STAGE;
        phi_m_add = phi_m_add_deg * M_PI / 180.0;
        alpha = lead_alpha_from_phi_max(phi_m_add);
        if (alpha < LEAD_ALPHA_MIN) alpha = LEAD_ALPHA_MIN;

        w_new = lead_find_new_crossover(&sys_scaled, alpha,
                                         LEAD_FREQ_MIN_DEFAULT,
                                         LEAD_FREQ_MAX_DEFAULT);
        if (w_new <= 0.0) w_new = w_gc_old * 1.5;
        T = 1.0 / (w_new * sqrt(alpha));
        K_c = K / alpha;

        if (!lead_init(&comp, K_c, 1.0/T, 1.0/(alpha*T), LEAD_TYPE_ANALYTIC))
            break;
    }

    result->compensator = comp;
    result->achieved_phase_margin = pm_achieved;
    result->achieved_gain_margin = lead_compensated_gain_margin(&comp, system);
    result->crossover_freq_new = lead_open_loop_crossover(&comp, system);
    result->phase_lead_added = phi_m_add_deg;
    result->alpha_actual = alpha;
    result->T_actual = T;
    result->iterations = iter;
    result->converged = converged;

    return converged || (fabs(pm_achieved - specs->phase_margin_desired) < 2.0);
}
/**
 * Determine gain K from steady-state error specification.
 * L5 - For unity feedback: e_ss = 1/(1+K_p) for step input.
 * K = (1-e_ss)/(e_ss*|plant_dc_gain|)
 */
double lead_design_gain_from_error(const lead_system_t *system, double e_ss) {
    if (!system || e_ss <= 0.0) return 1.0;
    double num_dc = system->tf.num[0];
    double den_dc = system->tf.den[0];
    if (fabs(den_dc) < 1e-300) den_dc = 1e-6;
    double plant_dc_gain = system->tf.gain * num_dc / den_dc;
    double K = (1.0 - e_ss) / (e_ss * fabs(plant_dc_gain));
    if (K < 0.01) K = 0.01;
    if (K > 1e6) K = 1e6;
    return K;
}

/**
 * Find frequency where |KG(jw)|_dB = -10*log10(alpha).
 * L5 - Log-spaced scan + bisection. Returns new gain crossover frequency.
 */
double lead_find_new_crossover(const lead_system_t *system, double alpha,
                                double w_min, double w_max) {
    if (!system || alpha <= 0.0 || alpha >= 1.0) return 0.0;
    double target_db = -10.0 * log10(alpha);
    double best_w = 0.0, best_diff = INFINITY;
    int n_scan = 200;

    for (int i = 0; i < n_scan; i++) {
        double frac = (double)i / (double)(n_scan - 1);
        double w = w_min * pow(w_max / w_min, frac);
        double num_re = 0.0, num_im = 0.0;
        { double sr = 1.0, si = 0.0;
          for (int j = 0; j <= system->tf.num_order; j++) {
            num_re += system->tf.num[j] * sr;
            num_im += system->tf.num[j] * si;
            double nr = -si * w, ni = sr * w;
            sr = nr; si = ni;
          }
        }
        double nms = num_re*num_re + num_im*num_im;
        double den_re = 0.0, den_im = 0.0;
        { double sr = 1.0, si = 0.0;
          for (int j = 0; j <= system->tf.den_order; j++) {
            den_re += system->tf.den[j] * sr;
            den_im += system->tf.den[j] * si;
            double nr = -si * w, ni = sr * w;
            sr = nr; si = ni;
          }
        }
        double dms = den_re*den_re + den_im*den_im;
        double mag = 1.0;
        if (dms > 1e-300 && nms > 1e-300)
            mag = system->tf.gain * sqrt(nms / dms);
        double mag_db = 20.0 * log10(mag);
        double diff = fabs(mag_db - target_db);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    if (best_w > 0.0 && best_diff > 0.01) {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_min) lo = w_min;
        if (hi > w_max) hi = w_max;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double nr = 0.0, ni = 0.0;
            { double sr = 1.0, si = 0.0;
              for (int j = 0; j <= system->tf.num_order; j++) {
                nr += system->tf.num[j] * sr;
                ni += system->tf.num[j] * si;
                double tr = -si*mid, ti = sr*mid;
                sr = tr; si = ti;
              }
            }
            double dr = 0.0, di = 0.0;
            { double sr = 1.0, si = 0.0;
              for (int j = 0; j <= system->tf.den_order; j++) {
                dr += system->tf.den[j] * sr;
                di += system->tf.den[j] * si;
                double tr = -si*mid, ti = sr*mid;
                sr = tr; si = ti;
              }
            }
            double nms2 = nr*nr+ni*ni, dms2 = dr*dr+di*di;
            double mag_mid = 1.0;
            if (dms2 > 1e-300 && nms2 > 1e-300)
                mag_mid = system->tf.gain * sqrt(nms2/dms2);
            if (20.0*log10(mag_mid) > target_db) lo = mid;
            else hi = mid;
            if ((hi-lo)/mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }
    return best_w;
}

/* =========================================================================
 * L5 - Root Locus Design
 * ========================================================================= */

/**
 * Design lead compensator via root-locus method.
 *
 * Given a desired closed-loop pole location s_d, compute the angle
 * deficiency and place lead zero/pole to supply the needed angle.
 *
 * L5 - Evans angle condition: sum_angle(s-z_i) - sum_angle(s-p_i) = -180
 * The compensator must contribute: phi_c = -180 - angle(G(s_d))
 *
 * For a lead compensator C(s) = K_c*(s+z_c)/(s+p_c):
 *   angle from zero at -z_c: atan2(im, re + z_c)
 *   angle from pole at -p_c: -atan2(im, re + p_c)
 *   Total: atan2(im, re+z_c) - atan2(im, re+p_c) = phi_c
 */
bool lead_design_root_locus(const lead_system_t *system,
                             lead_complex_t s_desired,
                             lead_design_result_t *result) {
    if (!system || !result) return false;
    memset(result, 0, sizeof(lead_design_result_t));

    /* Compute angle deficiency */
    double phi_c = lead_angle_deficiency(system, s_desired);
    if (phi_c <= 0.0) {
        /* No lead needed or deficiency in wrong direction */
        return false;
    }

    /* Place zero and pole */
    double z_c, p_c;
    if (!lead_place_zero_pole_angle(s_desired, phi_c, &z_c, &p_c))
        return false;

    /* Compute gain K_c from magnitude condition */
    /* |K_c*(s_d+z_c)/(s_d+p_c) * G(s_d)| = 1 */
    lead_complex_t s_plus_zc = {s_desired.re + z_c, s_desired.im};
    lead_complex_t s_plus_pc = {s_desired.re + p_c, s_desired.im};
    double num_mag = sqrt(s_plus_zc.re*s_plus_zc.re + s_plus_zc.im*s_plus_zc.im);
    double den_mag = sqrt(s_plus_pc.re*s_plus_pc.re + s_plus_pc.im*s_plus_pc.im);

    /* |G(s_d)| via Horner */
    double G_re = 0.0, G_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.num_order; i++) {
        G_re += system->tf.num[i] * sr;
        G_im += system->tf.num[i] * si;
        double nr = sr*s_desired.re - si*s_desired.im;
        double ni = sr*s_desired.im + si*s_desired.re;
        sr = nr; si = ni;
      }
    }
    double numG_mag = sqrt(G_re*G_re + G_im*G_im) * system->tf.gain;

    double D_re = 0.0, D_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.den_order; i++) {
        D_re += system->tf.den[i] * sr;
        D_im += system->tf.den[i] * si;
        double nr = sr*s_desired.re - si*s_desired.im;
        double ni = sr*s_desired.im + si*s_desired.re;
        sr = nr; si = ni;
      }
    }
    double denG_mag = sqrt(D_re*D_re + D_im*D_im);

    double G_mag = (denG_mag > 1e-300) ? numG_mag / denG_mag : INFINITY;
    double K_c = den_mag / (num_mag * G_mag);

    lead_compensator_t comp;
    if (!lead_init(&comp, K_c, z_c, p_c, LEAD_TYPE_ANALYTIC))
        return false;

    result->compensator = comp;
    result->alpha_actual = comp.alpha;
    result->T_actual = comp.T;
    result->converged = true;
    return true;
}

/**
 * Compute angle deficiency at a point s for root-locus design.
 *
 * L5 - angle deficiency = -M_PI - angle_G_s
 * where angle_G_s = atan2(Im{G(s)}, Re{G(s)})
 *
 * If deficiency > 0, need phase lead.
 * If deficiency < 0, need phase lag.
 */
double lead_angle_deficiency(const lead_system_t *system, lead_complex_t s) {
    if (!system) return 0.0;

    /* Evaluate G(s) */
    double G_re = 0.0, G_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.num_order; i++) {
        G_re += system->tf.num[i] * sr;
        G_im += system->tf.num[i] * si;
        double nr = sr*s.re - si*s.im;
        double ni = sr*s.im + si*s.re;
        sr = nr; si = ni;
      }
    }
    G_re *= system->tf.gain;
    G_im *= system->tf.gain;

    double D_re = 0.0, D_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.den_order; i++) {
        D_re += system->tf.den[i] * sr;
        D_im += system->tf.den[i] * si;
        double nr = sr*s.re - si*s.im;
        double ni = sr*s.im + si*s.re;
        sr = nr; si = ni;
      }
    }

    double den_mag_sq = D_re*D_re + D_im*D_im;
    double angle_G = 0.0;
    if (den_mag_sq > 1e-300) {
        double val_re = (G_re*D_re + G_im*D_im) / den_mag_sq;
        double val_im = (G_im*D_re - G_re*D_im) / den_mag_sq;
        angle_G = atan2(val_im, val_re);
    }

    /* For standard negative-feedback root locus:
     * Required: angle(G(s)) = -pi (mod 2*pi)
     * Deficiency: phi_needed = -pi - angle_G (bring into (-pi, pi]) */
    double phi_needed = -M_PI - angle_G;
    while (phi_needed <= -M_PI) phi_needed += 2.0*M_PI;
    while (phi_needed > M_PI) phi_needed -= 2.0*M_PI;

    /* For lead design, we need positive phase */
    if (phi_needed < 0.0) phi_needed += 2.0*M_PI;

    return phi_needed;
}

/**
 * Place lead zero and pole to provide angle phi_c at point s.
 *
 * L5 - Bisector method (Ogata 7-6):
 * Draw line from s to the real axis at angle theta = (pi - phi_c) / 2
 * from the horizontal. The intersection gives zero and pole locations.
 *
 * From geometry at s = x + j*y:
 *   Let theta_z = atan2(y, x + z_c)   [angle contributed by zero]
 *   Let theta_p = atan2(y, x + p_c)   [angle contributed by pole]
 *   Requirement: theta_z - theta_p = phi_c
 *
 * Using bisector: draw angle bisector from s to real axis.
 * Bisector angle = angle(s) + phi_c/2 (measured from positive real)
 * The intersection z_c is at distance d along the bisector scaled properly.
 */
bool lead_place_zero_pole_angle(lead_complex_t s, double phi_c,
                                 double *z_c_out, double *p_c_out) {
    if (!z_c_out || !p_c_out) return false;
    if (phi_c <= 0.0 || phi_c >= M_PI) return false;
    if (s.im <= 0.0) return false; /* must be in upper half-plane */

    /* Bisector method: angle between zero and pole vectors from s */
    double gamma = (M_PI - phi_c) / 2.0;
    double theta_s = atan2(s.im, -s.re); /* angle from s to origin */

    /* Place zero: extend from s toward real axis at angle (theta_s + gamma) */
    double theta_z = theta_s + gamma;
    /* Intersection with real axis: z_c = s.im / tan(theta_z) - s.re */
    double tan_theta_z = tan(theta_z);
    if (fabs(tan_theta_z) < 1e-12) return false;
    double z_c = s.im / tan_theta_z - s.re;
    if (z_c <= 0.0) z_c = fabs(z_c) * 0.1;

    /* Place pole: angle from s to pole is (theta_s - gamma) */
    double theta_p = theta_s - gamma;
    double tan_theta_p = tan(theta_p);
    if (fabs(tan_theta_p) < 1e-12) return false;
    double p_c = s.im / tan_theta_p - s.re;
    if (p_c <= 0.0) p_c = fabs(p_c);

    /* Ensure lead condition: |z_c| < |p_c| */
    if (z_c >= p_c) {
        /* Swap: use alternative placement */
        double tmp = z_c;
        z_c = p_c * 0.3;
        p_c = tmp * 3.0;
    }

    *z_c_out = z_c;
    *p_c_out = p_c;
    return true;
}

/* =========================================================================
 * L5 - Analytic Design from Time-Domain Specs
 * ========================================================================= */

/**
 * Analytic design: PO% and ts -> zeta, wn -> PM, w_gc -> design.
 *
 * L5 - Uses 2nd-order approximations to convert time-domain specs
 * to frequency-domain parameters, then applies Bode design method.
 */
bool lead_design_analytic(const lead_specs_t *specs,
                           const lead_system_t *system,
                           lead_design_result_t *result) {
    if (!specs || !system || !result) return false;

    /* Convert time-domain specs to frequency-domain */
    double zeta = lead_zeta_from_overshoot(specs->max_overshoot_pct);
    double wn = lead_wn_from_settling(zeta, specs->settling_time);
    double pm_desired = lead_damping_to_pm(zeta);
    /* w_gc derived by Bode method internally from PM and plant */
    (void)lead_crossover_from_bandwidth(
        lead_bandwidth_from_zeta_wn(zeta, wn));

    /* Build frequency-domain specs */
    lead_specs_t freq_specs;
    memset(&freq_specs, 0, sizeof(freq_specs));
    freq_specs.phase_margin_desired = pm_desired;
    freq_specs.gain_margin_desired = 10.0; /* typical */
    freq_specs.steady_state_error = specs->steady_state_error;
    freq_specs.use_frequency_domain = true;
    freq_specs.method = LEAD_METHOD_ANALYTIC;

    return lead_design_frequency(system, &freq_specs, result);
}

/**
 * Dominant closed-loop pole from (zeta, wn).
 *
 * s_d = -zeta*wn +/- j*wn*sqrt(1 - zeta^2)
 *
 * L3 - Standard 2nd-order pole location formula.
 */
lead_complex_t lead_dominant_pole(double zeta, double wn, bool upper_half) {
    lead_complex_t pole;
    pole.re = -zeta * wn;
    pole.im = wn * sqrt(1.0 - zeta * zeta);
    if (!upper_half) pole.im = -pole.im;
    return pole;
}

/* =========================================================================
 * L2 - Performance Prediction and Verification
 * ========================================================================= */

/**
 * Predict step response metrics from PM and crossover frequency.
 * Uses 2nd-order dominant-pole approximation.
 */
void lead_predict_performance(double phase_margin, double crossover_freq,
                               lead_performance_t *perf) {
    if (!perf) return;
    memset(perf, 0, sizeof(lead_performance_t));

    double zeta = lead_pm_to_damping(phase_margin);
    double wn = crossover_freq; /* approximate: w_gc �� wn for 2nd-order */

    perf->damping_ratio = zeta;
    perf->natural_freq = wn;
    perf->percent_overshoot = lead_overshoot_from_zeta(zeta);
    perf->settling_time = lead_settling_time_2pct(zeta, wn);
    perf->peak_time = lead_peak_time(zeta, wn);
    perf->rise_time = lead_rise_time(zeta, wn);
    perf->bandwidth = lead_bandwidth_from_zeta_wn(zeta, wn);
    perf->steady_state_error = 0.0; /* depends on system type */
}

/**
 * Verify compensator design by computing compensated closed-loop metrics.
 * Uses the 2nd-order approximation from achieved PM and crossover.
 */
void lead_verify_design(const lead_compensator_t *compensator,
                         const lead_system_t *system,
                         lead_performance_t *perf) {
    if (!compensator || !system || !perf) return;

    double pm = lead_compensated_phase_margin(compensator, system);
    double w_gc = lead_open_loop_crossover(compensator, system);
    lead_predict_performance(pm, w_gc, perf);
}

/**
 * Check if achieved performance meets all specifications.
 */
bool lead_check_specs(const lead_performance_t *perf,
                       const lead_specs_t *specs) {
    if (!perf || !specs) return false;

    if (perf->percent_overshoot > specs->max_overshoot_pct + 1.0)
        return false;
    if (perf->settling_time > specs->settling_time * 1.05)
        return false;
    if (perf->steady_state_error > specs->steady_state_error * 1.01)
        return false;
    return true;
}

/**
 * Accurate PM-zeta relationship for 2nd-order systems.
 * PM = atan(2*zeta/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)).
 */
double lead_zeta_to_pm(double zeta) {
    return lead_damping_to_pm(zeta);
}

/**
 * Bandwidth from (zeta, settling time).
 * w_BW = 4/(zeta*ts)*sqrt(sqrt(1+4*zeta^4)-2*zeta^2).
 */
double lead_bandwidth_from_zeta_ts(double zeta, double ts) {
    if (zeta <= 0.0 || ts <= 0.0) return 0.0;
    double wn = 4.0 / (zeta * ts);
    return lead_bandwidth_from_zeta_wn(zeta, wn);
}

/**
 * Approximate crossover from bandwidth: w_gc ~ w_BW / 1.6.
 */
double lead_crossover_from_bandwidth(double bandwidth) {
    return bandwidth / 1.6;
}

/**
 * Compute the required compensator DC gain K to achieve a desired
 * velocity error constant K_v for a type-1 system.
 *
 * L5 - For type-1: K_v = lim_{s->0} s*C(s)*G(s) = K * K_v_plant
 * => K = K_v_desired / K_v_plant
 */
double lead_gain_from_Kv(double K_v_desired, const lead_system_t *system) {
    if (!system || K_v_desired <= 0.0) return 1.0;

    /* K_v_plant = lim_{s->0} s*G(s) = num[0]*gain / den[1] */
    double num_0 = system->tf.num[0];
    double den_1 = system->tf.den[1]; /* coefficient of s in denominator */
    if (fabs(den_1) < 1e-300) return 1.0;

    double Kv_plant = system->tf.gain * num_0 / den_1;
    double K = K_v_desired / fabs(Kv_plant);
    if (K < 0.01) K = 0.01;
    return K;
}
