/**
 * @file lead_compensator.c
 * @brief Lead Compensator Core Implementation
 *
 * L1-L4: compensator creation, validation, complex evaluation,
 * phase/magnitude, maximum phase lead theorem, PM/zeta mapping.
 */

#include "lead_compensator.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* =========================================================================
 * L3 - Complex arithmetic utilities (internal use only)
 * ========================================================================= */

static lead_complex_t lc_div(lead_complex_t a, lead_complex_t b) {
    double d = b.re * b.re + b.im * b.im;
    if (d < 1e-300) {
        lead_complex_t r = {INFINITY, INFINITY};
        return r;
    }
    lead_complex_t r = {(a.re * b.re + a.im * b.im) / d,
                        (a.im * b.re - a.re * b.im) / d};
    return r;
}

/* =========================================================================
 * L1 - Compensator Initialization
 * ========================================================================= */

/**
 * Initialize from pole-zero: C(s) = K_c * (s + z_c) / (s + p_c)
 * Conditions: K_c > 0, z_c > 0, p_c > 0, z_c < p_c.
 * L1 - Core definition.
 */
bool lead_init(lead_compensator_t *comp, double K_c, double z_c,
               double p_c, lead_type_t type) {
    if (!comp) return false;
    if (K_c <= 0.0 || z_c <= 0.0 || p_c <= 0.0) return false;
    if (z_c >= p_c) return false;

    double alpha = z_c / p_c;
    if (alpha < LEAD_ALPHA_MIN || alpha >= 1.0) return false;

    comp->K_c     = K_c;
    comp->z_c     = z_c;
    comp->p_c     = p_c;
    comp->T       = 1.0 / z_c;
    comp->alpha   = alpha;
    comp->dc_gain = K_c * alpha;
    comp->type    = type;

    return true;
}

/**
 * Initialize from K-T-alpha: C(s) = K * (T*s + 1) / (alpha*T*s + 1)
 * L1 - Alternative parameterization.
 */
bool lead_init_from_KT_alpha(lead_compensator_t *comp, double K, double T,
                              double alpha, lead_type_t type) {
    if (!comp) return false;
    if (K <= 0.0 || T <= 0.0) return false;
    if (alpha < LEAD_ALPHA_MIN || alpha >= 1.0) return false;

    comp->K_c     = K / alpha;
    comp->z_c     = 1.0 / T;
    comp->p_c     = 1.0 / (alpha * T);
    comp->T       = T;
    comp->alpha   = alpha;
    comp->dc_gain = K;
    comp->type    = type;

    return true;
}

void lead_copy(lead_compensator_t *dest, const lead_compensator_t *src) {
    if (!dest || !src) return;
    memcpy(dest, src, sizeof(lead_compensator_t));
}

/**
 * Validate lead compensator.
 * L1 - Checks all consistency conditions.
 */
bool lead_validate(const lead_compensator_t *comp) {
    if (!comp) return false;
    if (comp->K_c <= 0.0) return false;
    if (comp->z_c <= 0.0 || comp->p_c <= 0.0) return false;
    if (comp->z_c >= comp->p_c) return false;
    if (comp->alpha < LEAD_ALPHA_MIN || comp->alpha >= 1.0) return false;
    if (comp->T <= 0.0) return false;
    double alpha_check = comp->z_c / comp->p_c;
    if (fabs(alpha_check - comp->alpha) > 1e-10) return false;
    return true;
}

/**
 * Convert compensator to transfer function.
 * C(s) = K_c*(s+z_c)/(s+p_c). num = [K_c*z_c, K_c], den = [p_c, 1].
 */
void lead_to_transfer_function(const lead_compensator_t *comp, lead_tf_t *tf) {
    if (!comp || !tf) return;
    memset(tf, 0, sizeof(lead_tf_t));
    tf->num_order = 1;
    tf->num[0] = comp->K_c * comp->z_c;
    tf->num[1] = comp->K_c;
    tf->den_order = 1;
    tf->den[0] = comp->p_c;
    tf->den[1] = 1.0;
    tf->gain = 1.0;
}

/* =========================================================================
 * L3 - Complex Evaluation
 * ========================================================================= */

lead_complex_t lead_evaluate(const lead_compensator_t *comp, lead_complex_t s) {
    lead_complex_t zero_val = {0.0, 0.0};
    if (!comp) return zero_val;

    lead_complex_t num = {s.re + comp->z_c, s.im};
    lead_complex_t den = {s.re + comp->p_c, s.im};

    double den_mag_sq = den.re * den.re + den.im * den.im;
    if (den_mag_sq < 1e-300) {
        lead_complex_t inf_val = {INFINITY, INFINITY};
        return inf_val;
    }

    lead_complex_t scaled_num = {comp->K_c * num.re, comp->K_c * num.im};
    return lc_div(scaled_num, den);
}

lead_complex_t lead_evaluate_jw(const lead_compensator_t *comp, double omega) {
    lead_complex_t s = {0.0, omega};
    return lead_evaluate(comp, s);
}

/* =========================================================================
 * L2 - Phase and Magnitude
 * ========================================================================= */

/**
 * Phase: phi(w) = atan(wT) - atan(alpha*wT)
 * L2 - Core concept.
 */
double lead_phase_at(const lead_compensator_t *comp, double omega) {
    if (!comp) return 0.0;
    double wT = omega * comp->T;
    double awT = comp->alpha * omega * comp->T;
    return atan(wT) - atan(awT);
}

/**
 * Magnitude: |C(jw)| = K * sqrt(1+(wT)^2) / sqrt(1+(alpha*wT)^2)
 * L2 - Core concept.
 */
double lead_magnitude_at(const lead_compensator_t *comp, double omega) {
    if (!comp) return 0.0;
    double wT = omega * comp->T;
    double awT = comp->alpha * omega * comp->T;
    double num = sqrt(1.0 + wT * wT);
    double den = sqrt(1.0 + awT * awT);
    if (den < 1e-300) return INFINITY;
    return comp->dc_gain * num / den;
}

/* =========================================================================
 * L4 - Maximum Phase Lead Theorem
 * ========================================================================= */

/**
 * Theorem: phi_max = arcsin((1-alpha)/(1+alpha))
 *
 * Proof: Let x = wT. phi(x) = atan(x) - atan(alpha*x).
 * d(phi)/dx = 0 => x = 1/sqrt(alpha) => omega_m = 1/(T*sqrt(alpha)).
 * tan(phi_max) = (1-alpha)/(2*sqrt(alpha)).
 * sin(phi_max) = (1-alpha)/(1+alpha).
 */
double lead_phi_max_from_alpha(double alpha) {
    if (alpha <= 0.0 || alpha >= 1.0) return 0.0;
    return asin((1.0 - alpha) / (1.0 + alpha));
}

/**
 * Inverse: alpha = (1 - sin(phi_m)) / (1 + sin(phi_m))
 */
double lead_alpha_from_phi_max(double phi_m) {
    if (phi_m <= 0.0) return 1.0;
    if (phi_m >= M_PI / 2.0) return LEAD_ALPHA_MIN;
    double sin_phi = sin(phi_m);
    double alpha = (1.0 - sin_phi) / (1.0 + sin_phi);
    if (alpha < LEAD_ALPHA_MIN) alpha = LEAD_ALPHA_MIN;
    if (alpha > 1.0) alpha = 1.0;
    return alpha;
}

/**
 * omega_m = 1/(T * sqrt(alpha))
 */
double lead_omega_max(const lead_compensator_t *comp) {
    if (!comp || !lead_validate(comp)) return 0.0;
    return 1.0 / (comp->T * sqrt(comp->alpha));
}

/**
 * |C(j*omega_m)| = K / sqrt(alpha)
 */
double lead_magnitude_at_omega_max(const lead_compensator_t *comp) {
    if (!comp || !lead_validate(comp)) return 0.0;
    return comp->dc_gain / sqrt(comp->alpha);
}

/* =========================================================================
 * L2 - Compensated System
 * ========================================================================= */

double lead_compensated_phase(const lead_compensator_t *comp,
                               const lead_system_t *system, double omega) {
    if (!comp || !system) return 0.0;
    double phi_c = lead_phase_at(comp, omega);

    double num_re = 0.0, num_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.num_order; i++) {
        num_re += system->tf.num[i] * sr;
        num_im += system->tf.num[i] * si;
        double nr = -si * omega, ni = sr * omega;
        sr = nr; si = ni;
      }
    }
    double den_re = 0.0, den_im = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.den_order; i++) {
        den_re += system->tf.den[i] * sr;
        den_im += system->tf.den[i] * si;
        double nr = -si * omega, ni = sr * omega;
        sr = nr; si = ni;
      }
    }
    double dms = den_re*den_re + den_im*den_im;
    double phi_g = 0.0;
    if (dms > 1e-300) {
        double Gr = system->tf.gain*(num_re*den_re + num_im*den_im)/dms;
        double Gi = system->tf.gain*(num_im*den_re - num_re*den_im)/dms;
        phi_g = atan2(Gi, Gr);
    }
    return phi_c + phi_g;
}

double lead_compensated_magnitude(const lead_compensator_t *comp,
                                   const lead_system_t *system, double omega) {
    if (!comp || !system) return 0.0;

    double mag_c = lead_magnitude_at(comp, omega);

    double nr = 0.0, ni = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.num_order; i++) {
        nr += system->tf.num[i] * sr;
        ni += system->tf.num[i] * si;
        double tr = -si*omega, ti = sr*omega;
        sr = tr; si = ti;
      }
    }
    double nms = nr*nr + ni*ni;

    double dr = 0.0, di = 0.0;
    { double sr = 1.0, si = 0.0;
      for (int i = 0; i <= system->tf.den_order; i++) {
        dr += system->tf.den[i] * sr;
        di += system->tf.den[i] * si;
        double tr = -si*omega, ti = sr*omega;
        sr = tr; si = ti;
      }
    }
    double dms = dr*dr + di*di;

    double mag_g = 0.0;
    if (dms > 1e-300) mag_g = system->tf.gain * sqrt(nms/dms);
    else if (nms < 1e-300) mag_g = 0.0;
    else mag_g = INFINITY;

    return mag_c * mag_g;
}

/* =========================================================================
 * L1 - DC and HF Gain
 * ========================================================================= */

double lead_dc_gain(const lead_compensator_t *comp) {
    if (!comp) return 0.0;
    return comp->dc_gain;
}

double lead_hf_gain(const lead_compensator_t *comp) {
    if (!comp) return 0.0;
    return comp->K_c;
}

/* =========================================================================
 * L3 - PM/Damping/PO/ts Relationships
 * ========================================================================= */

double lead_pm_to_damping(double phase_margin_deg) {
    if (phase_margin_deg <= 0.0) return 0.0;
    if (phase_margin_deg >= 90.0) return 1.0;
    double pm_rad = phase_margin_deg * M_PI / 180.0;
    double cos_pm = cos(pm_rad);
    if (cos_pm <= 0.0) return 1.0;
    double zeta = sin(pm_rad) / (2.0 * sqrt(cos_pm));
    if (zeta > 1.0) zeta = 1.0;
    if (zeta < 0.0) zeta = 0.0;
    return zeta;
}

double lead_damping_to_pm(double zeta) {
    if (zeta <= 0.0) return 0.0;
    if (zeta >= 1.0) return 90.0;
    double z2 = zeta * zeta;
    double inner = sqrt(1.0 + 4.0 * z2 * z2);
    double denom = sqrt(inner - 2.0 * z2);
    if (denom < 1e-300) return 90.0;
    return atan(2.0 * zeta / denom) * 180.0 / M_PI;
}

double lead_overshoot_from_zeta(double zeta) {
    if (zeta <= 0.0) return 100.0;
    if (zeta >= 1.0) return 0.0;
    return 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
}

double lead_zeta_from_overshoot(double po_pct) {
    if (po_pct <= 0.0) return 1.0;
    if (po_pct >= 100.0) return 0.0;
    double ln_po = log(po_pct / 100.0);
    return -ln_po / sqrt(M_PI * M_PI + ln_po * ln_po);
}

double lead_settling_time_2pct(double zeta, double wn) {
    if (zeta <= 0.0 || wn <= 0.0) return INFINITY;
    return 4.0 / (zeta * wn);
}

double lead_wn_from_settling(double zeta, double ts) {
    if (zeta <= 0.0 || ts <= 0.0) return INFINITY;
    return 4.0 / (zeta * ts);
}

double lead_peak_time(double zeta, double wn) {
    if (zeta >= 1.0 || wn <= 0.0) return INFINITY;
    return M_PI / (wn * sqrt(1.0 - zeta * zeta));
}

double lead_rise_time(double zeta, double wn) {
    if (zeta >= 1.0 || wn <= 0.0) return INFINITY;
    double beta = acos(zeta);
    return (M_PI - beta) / (wn * sqrt(1.0 - zeta * zeta));
}

double lead_damped_frequency(double zeta, double wn) {
    if (zeta >= 1.0) return 0.0;
    return wn * sqrt(1.0 - zeta * zeta);
}

double lead_bandwidth_from_zeta_wn(double zeta, double wn) {
    if (wn <= 0.0) return 0.0;
    if (zeta >= 1.0) return wn * (sqrt(2.0) - 1.0);
    double z2 = zeta * zeta;
    double inner = sqrt(2.0 - 4.0 * z2 + 4.0 * z2 * z2);
    return wn * sqrt(1.0 - 2.0 * z2 + inner);
}

/* =========================================================================
 * L2 - Frequency Response and Utilities
 * ========================================================================= */

void lead_freq_response(const lead_compensator_t *comp, double omega,
                         double *mag_out, double *phase_out) {
    if (!comp) {
        if (mag_out) *mag_out = 0.0;
        if (phase_out) *phase_out = 0.0;
        return;
    }
    if (mag_out) *mag_out = lead_magnitude_at(comp, omega);
    if (phase_out) *phase_out = lead_phase_at(comp, omega);
}

double lead_efficiency(const lead_compensator_t *comp, double omega) {
    if (!comp || !lead_validate(comp)) return 0.0;
    double phi_max = lead_phi_max_from_alpha(comp->alpha);
    if (phi_max < 1e-10) return 0.0;
    return lead_phase_at(comp, omega) / phi_max;
}

double lead_phase_bandwidth(const lead_compensator_t *comp, double threshold_ratio) {
    if (!comp || !lead_validate(comp)) return 0.0;
    if (threshold_ratio <= 0.0 || threshold_ratio >= 1.0) return 0.0;
    double phi_max = lead_phi_max_from_alpha(comp->alpha);
    double phi_target = phi_max * threshold_ratio;
    double w_m = lead_omega_max(comp);
    double w_lo = 0.0, w_hi = w_m;
    for (int iter = 0; iter < 60; iter++) {
        double w_mid = (w_lo + w_hi) * 0.5;
        if (lead_phase_at(comp, w_mid) < phi_target) w_lo = w_mid;
        else w_hi = w_mid;
        if ((w_hi - w_lo) < 1e-10) break;
    }
    double w_lower = w_lo;
    w_lo = w_m; w_hi = w_m * 100.0;
    for (int iter = 0; iter < 60; iter++) {
        double w_mid = (w_lo + w_hi) * 0.5;
        if (lead_phase_at(comp, w_mid) > phi_target) w_lo = w_mid;
        else w_hi = w_mid;
        if ((w_hi - w_lo) < 1e-10) break;
    }
    return w_hi - w_lower;
}

void lead_print(const lead_compensator_t *comp) {
    if (!comp) {
        printf("lead_compensator_t: NULL\n");
        return;
    }
    printf("Lead Compensator:\n");
    printf("  C(s) = %.6g * (s + %.6g) / (s + %.6g)\n",
           comp->K_c, comp->z_c, comp->p_c);
    printf("  C(s) = %.6g * (%.6g*s + 1) / (%.6g*s + 1)\n",
           comp->dc_gain, comp->T, comp->alpha * comp->T);
    printf("  alpha = %.6g, T = %.6g, DC gain = %.6g, HF gain = %.6g\n",
           comp->alpha, comp->T, comp->dc_gain, comp->K_c);
    printf("  omega_m = %.6g rad/s, phi_max = %.6g deg\n",
           lead_omega_max(comp),
           lead_phi_max_from_alpha(comp->alpha) * 180.0 / M_PI);
    printf("  |C(j*omega_m)| = %.6g (%.6g dB)\n",
           lead_magnitude_at_omega_max(comp),
           20.0 * log10(lead_magnitude_at_omega_max(comp)));
}