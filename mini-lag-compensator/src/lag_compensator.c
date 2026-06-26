/**
 * @file lag_compensator.c
 * @brief Core lag compensator implementation: constructors, evaluation,
 *        accessors, validation, transfer function conversion.
 *
 * Knowledge coverage:
 *   L1: All constructor variants, accessor functions
 *   L2: Steady-state error improvement calculation
 *   L3: s-domain complex evaluation, magnitude/phase computation
 *   L4: Stability and minimum-phase checks
 *   L5: Validation, string formatting, comparison
 *
 * Textbook: Ogata, "Modern Control Engineering" (2010)
 *           Franklin, Powell, Emami-Naeini, "Feedback Control of Dynamic Systems"
 */

#include "lag_compensator.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ==========================================================================
 * L1: Constructors
 * ========================================================================== */

LagCompensator lag_create(double Kc, double T, double beta) {
    LagCompensator lag;
    lag.Kc = Kc;
    lag.T = T;
    lag.beta = beta;
    /* zero at s = -1/T */
    lag.zero = -1.0 / T;
    /* pole at s = -1/(beta*T); closer to origin than zero (|pole| < |zero|) */
    lag.pole = -1.0 / (beta * T);
    /* DC gain: G_c(0) = Kc * (0+1)/(0+1) = Kc */
    lag.dc_gain = Kc;
    /* HF gain: G_c(inf) = Kc / beta */
    lag.hf_gain = Kc / beta;
    /* Maximum phase lag frequency: w_max = 1/(T*sqrt(beta)) */
    lag.max_lag_freq = 1.0 / (T * sqrt(beta));
    /* Maximum phase lag magnitude:
     * |phi_max| = arctan(sqrt(beta)) - arctan(1/sqrt(beta))
     *           = arcsin((beta-1)/(beta+1))  [positive value = lag magnitude]
     */
    double sqrt_beta = sqrt(beta);
    lag.max_phase_lag_rad = atan(sqrt_beta) - atan(1.0 / sqrt_beta);
    /* Ensure positive (magnitude of lag) */
    if (lag.max_phase_lag_rad < 0) {
        lag.max_phase_lag_rad = -lag.max_phase_lag_rad;
    }
    return lag;
}

LagCompensator lag_create_from_pole_zero(double Kc, double zero, double pole) {
    /* Input: G_c(s) = Kc * (s + zero) / (s + pole)
     *        with zero = -z_c (z_c > 0), pole = -p_c (p_c > 0)
     *        lag condition: p_c < z_c => |pole| < |zero|
     *
     * Convert to (Kc, T, beta) form:
     *   T = 1/z_c = -1/zero
     *   beta = z_c/p_c = zero/pole
     *
     * Verify: G_c(s) = Kc * (s + z_c) / (s + p_c)
     *                = Kc * (z_c*(s/z_c + 1)) / (p_c*(s/p_c + 1))
     *                = Kc*(z_c/p_c) * (s/z_c + 1) / (s/p_c + 1)
     *                = Kc*beta * (T*s + 1) / (beta*T*s + 1)
     * Wait - there's a subtlety. Let me re-derive:
     *
     * G_c(s) = Kc * (s + z_c) / (s + p_c)
     *        = Kc * z_c * (s/z_c + 1) / [p_c * (s/p_c + 1)]
     *        = Kc * (z_c/p_c) * (s/z_c + 1) / (s/p_c + 1)
     *        = Kc * beta * (T*s + 1) / (beta*T*s + 1)
     * where T = 1/z_c, beta = z_c/p_c.
     *
     * But the standard form is G_c(s) = Kc' * (T*s + 1) / (beta*T*s + 1)
     * with DC gain = Kc'. Comparing:
     * At s=0: standard gives Kc', pole-zero gives Kc*z_c/p_c = Kc*beta.
     * So Kc' = Kc * beta.
     *
     * Therefore the compensator parameters are:
     *   T = 1/z_c
     *   beta = z_c/p_c
     *   Kc' = Kc * beta
     *
     * But wait, in our standard form, Kc IS the DC gain.
     * So if we want DC gain = Kc_target, we set:
     *   internal_Kc = Kc_target
     * and the pole-zero form DC gain = internal_Kc * (z_c/p_c) = Kc_target * beta
     *
     * There are two conventions. We adopt:
     *   lag_create(Kc, T, beta) => G_c(s) = Kc * (Ts+1)/(beta*Ts+1)
     *   DC gain = Kc, HF gain = Kc/beta
     *
     * For pole-zero input G_c(s) = Kc_in * (s+z)/(s+p):
     *   DC gain = Kc_in * z/p
     *   So effective Kc = Kc_in * z/p = Kc_in * beta
     *
     * We compute T = 1/z, beta = z/p, Kc_eff = Kc_in * beta,
     * then call lag_create(Kc_eff, T, beta).
     */
    double z = -zero;  /* zero = -z, so z > 0 */
    double p = -pole;  /* pole = -p, so p > 0 */
    double z_over_p = z / p;  /* beta */
    double T_from_z = 1.0 / z;
    double Kc_eff = Kc * z_over_p;
    return lag_create(Kc_eff, T_from_z, z_over_p);
}

LagCompensator lag_create_from_corners(double Kc, double omega_z, double omega_p) {
    /* Given G_c(s) = Kc * (s + w_z) / (s + w_p) with w_p < w_z.
     * This is equivalent to pole at -w_p, zero at -w_z.
     * Converts to: T = 1/w_z, beta = w_z/w_p.
     * DC gain = Kc * w_z/w_p = Kc * beta.
     */
    double beta = omega_z / omega_p;
    double T = 1.0 / omega_z;
    double Kc_eff = Kc * beta;
    return lag_create(Kc_eff, T, beta);
}

LagCompensator lag_create_identity(void) {
    /* Identity compensator: G_c(s) = 1.
     * Achieved with Kc=1, any T>0, beta=1 (no separation).
     * beta=1 means pole = zero, giving pole-zero cancellation.
     */
    LagCompensator lag;
    lag.Kc = 1.0;
    lag.T = 1.0;
    lag.beta = 1.0;
    lag.zero = -1.0;
    lag.pole = -1.0;
    lag.dc_gain = 1.0;
    lag.hf_gain = 1.0;
    lag.max_lag_freq = 1.0;
    lag.max_phase_lag_rad = 0.0;  /* no lag when beta=1 */
    return lag;
}

/* ==========================================================================
 * L1: Accessors — all O(1)
 * ========================================================================== */

double lag_get_zero(const LagCompensator *lag) {
    return lag->zero;
}

double lag_get_pole(const LagCompensator *lag) {
    return lag->pole;
}

double lag_get_dc_gain(const LagCompensator *lag) {
    return lag->dc_gain;
}

double lag_get_hf_gain(const LagCompensator *lag) {
    return lag->hf_gain;
}

double lag_get_beta(const LagCompensator *lag) {
    return lag->beta;
}

double lag_get_time_constant(const LagCompensator *lag) {
    return lag->T;
}

double lag_get_zero_pole_ratio(const LagCompensator *lag) {
    return lag->beta;
}

/* ==========================================================================
 * L2: Steady-state error improvement
 * ========================================================================== */

double lag_ess_improvement(const LagCompensator *lag, LagESSType ess_type) {
    /* Theorem: For unity feedback, the lag compensator increases
     * the low-frequency open-loop gain by factor Kc.
     *
     * For all input types (step/ramp/parabolic), the error constants
     * (Kp, Kv, Ka) are multiplied by Kc, and steady-state error is
     * reduced by factor Kc.
     *
     * e_ss_new = e_ss_old / Kc
     * improvement factor = e_ss_old / e_ss_new = Kc
     *
     * The factor is the same regardless of input type because the
     * compensator affects only the DC gain, not the system type.
     */
    (void)ess_type;  /* same improvement factor for all types */
    return lag->Kc;
}

/* ==========================================================================
 * L3: s-domain evaluation
 * ========================================================================== */

LagComplex lag_eval_s(const LagCompensator *lag, LagComplex s) {
    /* Evaluate G_c(s) = Kc * (T*s + 1) / (beta*T*s + 1)
     *
     * Let s = sigma + j*w.
     *
     * Numerator:   Kc * (T*(sigma+jw) + 1) = Kc*((T*sigma+1) + j*(T*w))
     * Denominator: beta*T*(sigma+jw) + 1 = (beta*T*sigma+1) + j*(beta*T*w)
     *
     * Result = num/den = (num_re + j*num_im) / (den_re + j*den_im)
     *        = (num_re*den_re + num_im*den_im + j*(num_im*den_re - num_re*den_im))
     *          / (den_re^2 + den_im^2)
     */
    double T = lag->T;
    double beta = lag->beta;
    double Kc = lag->Kc;

    /* Numerator: Kc * (T*s + 1) */
    double num_re = Kc * (T * s.re + 1.0);
    double num_im = Kc * (T * s.im);

    /* Denominator: beta*T*s + 1 */
    double den_re = beta * T * s.re + 1.0;
    double den_im = beta * T * s.im;

    /* Complex division */
    double den_mag2 = den_re * den_re + den_im * den_im;

    LagComplex result;
    if (den_mag2 < 1e-300) {
        /* Pole encountered: s = -1/(beta*T) + j*0.
         * Return infinity with correct sign. */
        result.re = (num_re * den_re >= 0) ? INFINITY : -INFINITY;
        result.im = (num_im * den_re >= 0) ? INFINITY : -INFINITY;
    } else {
        result.re = (num_re * den_re + num_im * den_im) / den_mag2;
        result.im = (num_im * den_re - num_re * den_im) / den_mag2;
    }
    return result;
}

double lag_eval_magnitude(const LagCompensator *lag, double omega) {
    /* |G_c(jw)| = Kc * sqrt(1 + (w*T)^2) / sqrt(1 + (w*beta*T)^2)
     *
     * Derivation:
     *   G_c(jw) = Kc * (1 + j*w*T) / (1 + j*w*beta*T)
     *   |G_c| = Kc * |1+j*w*T| / |1+j*w*beta*T|
     *         = Kc * sqrt(1+w^2*T^2) / sqrt(1+w^2*beta^2*T^2)
     */
    double wT = omega * lag->T;
    double w_beta_T = omega * lag->beta * lag->T;
    return lag->Kc * sqrt(1.0 + wT * wT) / sqrt(1.0 + w_beta_T * w_beta_T);
}

double lag_eval_phase(const LagCompensator *lag, double omega) {
    /* phi(w) = arctan(w*T) - arctan(w*beta*T)
     *
     * Since beta > 1, the second term dominates and phi(w) < 0 for w > 0.
     * This is the defining characteristic of a lag compensator:
     * it produces phase lag (negative phase).
     */
    double wT = omega * lag->T;
    double w_beta_T = omega * lag->beta * lag->T;
    return atan(wT) - atan(w_beta_T);
}

LagFreqPoint lag_eval_frequency(const LagCompensator *lag, double omega) {
    /* Complete frequency response evaluation at a single frequency.
     * Returns magnitude (linear and dB), phase (rad and deg),
     * and real/imaginary parts for Nyquist analysis.
     */
    LagFreqPoint fp;
    fp.omega = omega;

    double wT = omega * lag->T;
    double w_beta_T = omega * lag->beta * lag->T;

    /* G_c(jw) = Kc * (1 + j*wT) / (1 + j*w_beta_T) */
    double den_mag2 = 1.0 + w_beta_T * w_beta_T;

    /* Real part: Re = Kc * (1 + w^2*T*beta*T) / (1 + w^2*beta^2*T^2)
     * Imag part: Im = Kc * (w*T - w*beta*T) / (1 + w^2*beta^2*T^2)
     */
    fp.real_part = lag->Kc * (1.0 + wT * w_beta_T) / den_mag2;
    fp.imag_part = lag->Kc * (wT - w_beta_T) / den_mag2;

    /* Magnitude */
    fp.magnitude = sqrt(fp.real_part * fp.real_part +
                        fp.imag_part * fp.imag_part);
    if (fp.magnitude < 1e-300) {
        fp.magnitude_db = -INFINITY;
    } else {
        fp.magnitude_db = 20.0 * log10(fp.magnitude);
    }

    /* Phase */
    fp.phase_rad = atan2(fp.imag_part, fp.real_part);
    fp.phase_deg = fp.phase_rad * 180.0 / M_PI;

    return fp;
}

/* ==========================================================================
 * L2: Corner frequencies and characteristic values
 * ========================================================================== */

void lag_get_corner_frequencies(const LagCompensator *lag,
                                double *omega_low, double *omega_high) {
    /* Corner frequencies are where the Bode asymptotes change:
     *   w_low  = 1/(beta*T)  -- pole corner (lower frequency)
     *   w_high = 1/T         -- zero corner (higher frequency)
     *
     * Below w_low: flat asymptote at 20*log10(Kc)
     * Between w_low and w_high: -20 dB/dec slope (lag region)
     * Above w_high: flat asymptote at 20*log10(Kc/beta)
     */
    *omega_low = 1.0 / (lag->beta * lag->T);
    *omega_high = 1.0 / lag->T;
}

double lag_get_max_lag_frequency(const LagCompensator *lag) {
    /* Frequency at which phase lag is maximum.
     *
     * Theorem: Differentiating phi(w) = arctan(wT) - arctan(w*beta*T)
     * and setting d(phi)/dw = 0 gives:
     *   w_max = 1 / (T * sqrt(beta))
     *
     * This is the geometric mean of the two corner frequencies:
     *   w_max = sqrt(w_low * w_high)
     *        = sqrt(1/(beta*T) * 1/T) = 1/(T*sqrt(beta))
     */
    return lag->max_lag_freq;
}

double lag_get_max_phase_lag(const LagCompensator *lag) {
    /* Magnitude of the maximum phase lag.
     *
     * At w_max = 1/(T*sqrt(beta)):
     *   phi_max = arctan(1/sqrt(beta)) - arctan(sqrt(beta))
     *
     * Using the identity: arctan(x) - arctan(y) = arctan((x-y)/(1+xy))
     *   phi_max = arctan((1/sqrt(beta) - sqrt(beta))/(1 + 1))
     *           = arctan((1-beta)/(2*sqrt(beta)))
     *
     * Alternatively: sin(phi_max) = (1-beta)/(1+beta)
     *                phi_max = arcsin((1-beta)/(1+beta))
     *
     * Since beta > 1, phi_max is negative (phase lag).
     * We return the magnitude (positive value).
     */
    return lag->max_phase_lag_rad;
}

/* ==========================================================================
 * L2: Asymptotic behavior
 * ========================================================================== */

double lag_low_freq_asymptote_db(const LagCompensator *lag) {
    /* As w -> 0: |G_c(jw)| -> Kc, so asymptote = 20*log10(Kc). */
    if (lag->Kc <= 0) return -INFINITY;
    return 20.0 * log10(lag->Kc);
}

double lag_high_freq_asymptote_db(const LagCompensator *lag) {
    /* As w -> inf: |G_c(jw)| -> Kc/beta, so asymptote = 20*log10(Kc/beta). */
    double hf_gain = lag->Kc / lag->beta;
    if (hf_gain <= 0) return -INFINITY;
    return 20.0 * log10(hf_gain);
}

/* ==========================================================================
 * L3: Transfer function conversion
 * ========================================================================== */

LagTransferFunction lag_to_transfer_function(const LagCompensator *lag) {
    /* G_c(s) = Kc * (T*s + 1) / (beta*T*s + 1)
     *        = (Kc*T*s + Kc) / (beta*T*s + 1)
     *
     * Numerator:  n1*s + n0 = Kc*T*s + Kc
     *   n0 = Kc, n1 = Kc*T
     * Denominator: d1*s + d0 = beta*T*s + 1
     *   d0 = 1, d1 = beta*T
     */
    LagTransferFunction tf;

    /* Numerator: order 1 */
    tf.numerator.order = 1;
    tf.numerator.coeff = (double*)malloc(2 * sizeof(double));
    if (tf.numerator.coeff) {
        tf.numerator.coeff[0] = lag->Kc;        /* constant term */
        tf.numerator.coeff[1] = lag->Kc * lag->T; /* s term */
    }

    /* Denominator: order 1 */
    tf.denominator.order = 1;
    tf.denominator.coeff = (double*)malloc(2 * sizeof(double));
    if (tf.denominator.coeff) {
        tf.denominator.coeff[0] = 1.0;               /* constant term */
        tf.denominator.coeff[1] = lag->beta * lag->T; /* s term */
    }

    /* DC gain */
    if (tf.denominator.coeff && tf.denominator.coeff[0] != 0.0) {
        tf.dc_gain = tf.numerator.coeff[0] / tf.denominator.coeff[0];
    } else {
        tf.dc_gain = INFINITY;
    }

    return tf;
}

LagCompensator lag_from_first_order_tf(double b1, double b0, double a1, double a0) {
    /* Convert G(s) = (b1*s + b0) / (a1*s + a0) to lag compensator form.
     *
     * Requirements:
     *   1. a0 > 0, a1 > 0  =>  stable pole at -a0/a1 < 0
     *   2. b0 > 0, b1 > 0  =>  minimum-phase zero at -b0/b1 < 0
     *   3. a0/a1 < b0/b1   =>  |pole| < |zero| (lag condition)
     *
     * The lag compensator form: G_c(s) = Kc * (T*s + 1) / (beta*T*s + 1)
     *
     * Rewrite G(s):
     *   G(s) = (b0/a0) * (b1/b0*s + 1) / (a1/a0*s + 1)
     *
     * Identify: Kc = b0/a0, T = b1/b0, beta*T = a1/a0
     * So: beta = (a1/a0) / (b1/b0) = (a1*b0)/(a0*b1)
     *
     * But wait, we need beta > 1 for a lag compensator.
     * beta = (a1/a0)/(b1/b0) = (a1*b0)/(a0*b1)
     *
     * The lag condition |pole| < |zero| means:
     *   a0/a1 < b0/b1  =>  a0*b1 < a1*b0  =>  (a1*b0)/(a0*b1) > 1
     *
     * So beta > 1 is equivalent to the lag condition. Good.
     */

    /* Check basic validity */
    if (a0 <= 0 || b0 <= 0) {
        /* Cannot form valid lag compensator from these params.
         * Return identity as fallback. */
        return lag_create_identity();
    }

    double Kc = b0 / a0;
    double T = b1 / b0;
    double beta_T = a1 / a0;
    double beta = beta_T / T;

    /* Ensure beta > 1 (lag condition) */
    if (beta <= 1.0) {
        beta = 1.0 + 1e-6;  /* clamp to near-identity */
    }

    return lag_create(Kc, T, beta);
}

/* ==========================================================================
 * L4: Stability and minimum-phase checks
 * ========================================================================== */

int lag_is_stable(const LagCompensator *lag) {
    /* The lag compensator is stable if its pole is in the left half-plane.
     * pole = -1/(beta*T), so pole < 0 iff beta*T > 0.
     * Since beta > 1 > 0 and T > 0, the compensator is always stable
     * for valid parameters.
     */
    return (lag->pole < 0.0) ? 1 : 0;
}

int lag_is_minimum_phase(const LagCompensator *lag) {
    /* The compensator is minimum-phase if its zero is in the LHP.
     * zero = -1/T, so zero < 0 iff T > 0.
     * Always true for valid parameters.
     */
    return (lag->zero < 0.0) ? 1 : 0;
}

/* ==========================================================================
 * L5: Validation
 * ========================================================================== */

int lag_validate(const LagCompensator *lag) {
    /* Comprehensive parameter validation.
     *
     * Check order is important: structural validity first,
     * then numerical validity, then domain constraints.
     */

    /* Check for NaN */
    if (isnan(lag->Kc) || isnan(lag->T) || isnan(lag->beta)) {
        return -5;
    }

    /* Check for infinity */
    if (isinf(lag->Kc) || isinf(lag->T) || isinf(lag->beta)) {
        return -5;
    }

    /* Kc must be positive */
    if (lag->Kc <= 0.0) {
        return -1;
    }

    /* T must be positive (ensures zero in LHP) */
    if (lag->T <= 0.0) {
        return -2;
    }

    /* beta must be > 1 (lag condition; beta=1 is identity, not lag) */
    if (lag->beta <= 1.0) {
        return -3;
    }

    /* Verify pole-zero ordering: |pole| < |zero| (pole closer to origin)
     * This is the defining characteristic of a lag compensator. */
    if (fabs(lag->pole) >= fabs(lag->zero)) {
        return -4;
    }

    /* Verify internal consistency: pole and zero should be negative real */
    if (lag->pole >= 0.0 || lag->zero >= 0.0) {
        return -4;
    }

    /* Verify derived values are consistent */
    double expected_pole = -1.0 / (lag->beta * lag->T);
    double expected_zero = -1.0 / lag->T;
    if (fabs(lag->pole - expected_pole) > 1e-9 * fabs(expected_pole)) {
        return -4;
    }
    if (fabs(lag->zero - expected_zero) > 1e-9 * fabs(expected_zero)) {
        return -4;
    }

    return 0;  /* valid */
}

const char* lag_validate_error_string(int error_code) {
    switch (error_code) {
        case 0:   return "Valid";
        case -1:  return "Kc must be positive";
        case -2:  return "T must be positive";
        case -3:  return "beta must be greater than 1";
        case -4:  return "Pole-zero ordering violates lag condition (|pole| < |zero| required)";
        case -5:  return "Non-finite values detected (NaN or Inf)";
        default:  return "Unknown validation error";
    }
}

/* ==========================================================================
 * L5: Utility functions
 * ========================================================================== */

void lag_to_string(const LagCompensator *lag, char *buf, int bufsz) {
    /* Format: G_c(s) = Kc*(T*s+1)/(beta*T*s+1) */
    if (buf && bufsz > 0) {
        snprintf(buf, bufsz,
                 "G_c(s) = %.4g*(%.4g*s+1)/(%.4g*%.4g*s+1) = %.4g*(%.4g*s+1)/(%.4g*s+1)",
                 lag->Kc, lag->T, lag->beta, lag->T,
                 lag->Kc, lag->T, lag->beta * lag->T);
    }
}

int lag_equals(const LagCompensator *a, const LagCompensator *b, double tol) {
    /* Compare Kc, T, beta independently.
     * Two compensators are equal if all three defining parameters match.
     */
    if (fabs(a->Kc - b->Kc) > tol) return 0;
    if (fabs(a->T - b->T) > tol) return 0;
    if (fabs(a->beta - b->beta) > tol) return 0;
    return 1;
}