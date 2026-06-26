/**
 * @file lead_frequency.c
 * @brief Frequency Domain Analysis for Lead Compensator Design
 *
 * L3 - Complex transfer function evaluation (Horner's method)
 * L4 - Stability margins, Nyquist criterion, Routh-Hurwitz
 * L5 - Bode plot generation, crossover detection
 *
 * Reference: Ogata Ch.8, Dorf&Bishop Ch.10, Franklin et al. Ch.6
 * MIT 6.302, Berkeley ME132, Cambridge 3F2, ETH 151-0591
 */

#include "lead_frequency.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* =========================================================================
 * Forward declarations of internal helpers
 * ========================================================================= */

static void eval_poly_complex(const double *coeffs, int order,
                               lead_complex_t s, lead_complex_t *result);
static double eval_poly_mag_sq(const double *coeffs, int order, double omega);

/* =========================================================================
 * L3 - Transfer Function Evaluation
 * ========================================================================= */

/**
 * Evaluate polynomial at complex s using Horner's method.
 * p(s) = coeffs[0] + coeffs[1]*s + coeffs[2]*s^2 + ... + coeffs[n]*s^n
 *
 * Horner: p(s) = coeffs[n]; for i=n-1..0: p = p*s + coeffs[i]
 * L3 - Polynomial evaluation, O(n).
 */
static void eval_poly_complex(const double *coeffs, int order,
                               lead_complex_t s, lead_complex_t *result) {
    result->re = 0.0;
    result->im = 0.0;

    if (order < 0) return;

    /* Horner: start from highest power */
    for (int i = order; i >= 0; i--) {
        /* result = result * s + coeffs[i] */
        double new_re = result->re * s.re - result->im * s.im + coeffs[i];
        double new_im = result->re * s.im + result->im * s.re;
        result->re = new_re;
        result->im = new_im;
    }
}

/**
 * Evaluate transfer function at complex s.
 * G(s) = gain * num(s) / den(s)
 * L3 - Rational function evaluation, O(max(num_order, den_order)).
 */
lead_complex_t lead_tf_evaluate(const lead_tf_t *tf, lead_complex_t s) {
    lead_complex_t zero_val = {0.0, 0.0};
    if (!tf) return zero_val;

    lead_complex_t num_val, den_val;
    eval_poly_complex(tf->num, tf->num_order, s, &num_val);
    eval_poly_complex(tf->den, tf->den_order, s, &den_val);

    double den_mag_sq = den_val.re * den_val.re + den_val.im * den_val.im;
    if (den_mag_sq < 1e-300) {
        lead_complex_t inf_val = {INFINITY, INFINITY};
        return inf_val;
    }

    /* G(s) = gain * num(s) / den(s) */
    lead_complex_t result;
    double scale = tf->gain / den_mag_sq;
    result.re = scale * (num_val.re * den_val.re + num_val.im * den_val.im);
    result.im = scale * (num_val.im * den_val.re - num_val.re * den_val.im);
    return result;
}

/**
 * Evaluate magnitude of polynomial at s = j*omega.
 * |p(jw)|^2 = Re(p(jw))^2 + Im(p(jw))^2
 * L3 - Frequency-domain magnitude, O(n).
 */
static double eval_poly_mag_sq(const double *coeffs, int order, double omega) {
    double re = 0.0, im = 0.0;

    for (int i = 0; i <= order; i++) {
        double s_re = 1.0, s_im = 0.0;
        /* Compute (j*w)^i */
        for (int k = 0; k < i; k++) {
            double nr = -s_im * omega;
            double ni = s_re * omega;
            s_re = nr;
            s_im = ni;
        }
        re += coeffs[i] * s_re;
        im += coeffs[i] * s_im;
    }

    return re * re + im * im;
}

/** Magnitude of transfer function |G(jw)| */
double lead_tf_magnitude(const lead_tf_t *tf, double omega) {
    if (!tf) return 0.0;
    double num_msq = eval_poly_mag_sq(tf->num, tf->num_order, omega);
    double den_msq = eval_poly_mag_sq(tf->den, tf->den_order, omega);
    if (den_msq < 1e-300) return INFINITY;
    return tf->gain * sqrt(num_msq / den_msq);
}

/** Phase of transfer function angle(G(jw)) in radians */
double lead_tf_phase(const lead_tf_t *tf, double omega) {
    if (!tf) return 0.0;
    lead_complex_t s = {0.0, omega};
    lead_complex_t val = lead_tf_evaluate(tf, s);
    return atan2(val.im, val.re);
}

/** Evaluate plant G(s) at complex s */
lead_complex_t lead_system_evaluate(const lead_system_t *sys, lead_complex_t s) {
    if (!sys) {
        lead_complex_t z = {0.0, 0.0};
        return z;
    }
    return lead_tf_evaluate(&sys->tf, s);
}

double lead_system_magnitude(const lead_system_t *sys, double omega) {
    if (!sys) return 0.0;
    return lead_tf_magnitude(&sys->tf, omega);
}

double lead_system_phase(const lead_system_t *sys, double omega) {
    if (!sys) return 0.0;
    return lead_tf_phase(&sys->tf, omega);
}
/* =========================================================================
 * L5 - Bode Plot Generation
 * ========================================================================= */

/**
 * Generate Bode plot data with logarithmic frequency spacing.
 * L5 - Bode plot is the foundation of frequency-domain design.
 */
void lead_bode_compute(const lead_tf_t *tf, double f_min, double f_max,
                        int num_points, lead_bode_data_t *bode) {
    if (!tf || !bode) return;
    if (num_points > LEAD_MAX_FREQ_POINTS) num_points = LEAD_MAX_FREQ_POINTS;
    if (num_points < 2) return;

    memset(bode, 0, sizeof(lead_bode_data_t));
    bode->freq_min = f_min;
    bode->freq_max = f_max;
    bode->log_spacing = true;
    bode->num_points = num_points;

    double log_fmin = log10(f_min);
    double log_fmax = log10(f_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_fmin + frac * (log_fmax - log_fmin));

        lead_complex_t s = {0.0, w};
        lead_complex_t val = lead_tf_evaluate(tf, s);

        bode->points[i].frequency = w;
        bode->points[i].real_part = val.re;
        bode->points[i].imag_part = val.im;

        double mag = sqrt(val.re * val.re + val.im * val.im);
        bode->points[i].magnitude = mag;
        bode->points[i].magnitude_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);

        double phase = atan2(val.im, val.re);
        bode->points[i].phase_rad = phase;
        bode->points[i].phase_deg = phase * 180.0 / M_PI;
    }
}

/** Bode plot of compensated open-loop C(s)*G(s) */
void lead_bode_compensated(const lead_compensator_t *comp,
                            const lead_system_t *sys,
                            double f_min, double f_max, int num_points,
                            lead_bode_data_t *bode) {
    if (!comp || !sys || !bode) return;
    if (num_points > LEAD_MAX_FREQ_POINTS) num_points = LEAD_MAX_FREQ_POINTS;
    if (num_points < 2) return;

    memset(bode, 0, sizeof(lead_bode_data_t));
    bode->freq_min = f_min;
    bode->freq_max = f_max;
    bode->log_spacing = true;
    bode->num_points = num_points;

    double log_fmin = log10(f_min);
    double log_fmax = log10(f_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_fmin + frac * (log_fmax - log_fmin));

        double mag = lead_compensated_magnitude(comp, sys, w);
        double phase = lead_compensated_phase(comp, sys, w);

        bode->points[i].frequency = w;
        bode->points[i].magnitude = mag;
        bode->points[i].magnitude_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        bode->points[i].phase_rad = phase;
        bode->points[i].phase_deg = phase * 180.0 / M_PI;
        bode->points[i].real_part = mag * cos(phase);
        bode->points[i].imag_part = mag * sin(phase);
    }
}

/** Bode plot of lead compensator alone */
void lead_bode_compensator_only(const lead_compensator_t *comp,
                                 double f_min, double f_max, int num_points,
                                 lead_bode_data_t *bode) {
    if (!comp || !bode) return;
    if (num_points > LEAD_MAX_FREQ_POINTS) num_points = LEAD_MAX_FREQ_POINTS;
    if (num_points < 2) return;

    memset(bode, 0, sizeof(lead_bode_data_t));
    bode->freq_min = f_min;
    bode->freq_max = f_max;
    bode->log_spacing = true;
    bode->num_points = num_points;

    double log_fmin = log10(f_min);
    double log_fmax = log10(f_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_fmin + frac * (log_fmax - log_fmin));

        double mag = lead_magnitude_at(comp, w);
        double phase = lead_phase_at(comp, w);

        bode->points[i].frequency = w;
        bode->points[i].magnitude = mag;
        bode->points[i].magnitude_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        bode->points[i].phase_rad = phase;
        bode->points[i].phase_deg = phase * 180.0 / M_PI;
        bode->points[i].real_part = mag * cos(phase);
        bode->points[i].imag_part = mag * sin(phase);
    }
}

/* =========================================================================
 * L4-L5 - Stability Margins
 * ========================================================================= */

/**
 * Find gain crossover frequency: |KG(jw_gc)| = 1 (0 dB).
 * L4 - Gain crossover is where magnitude crosses unity.
 * Uses log-spaced scan + bisection refinement.
 */
double lead_find_gain_crossover(const lead_system_t *sys, double w_min,
                                 double w_max) {
    if (!sys) return 0.0;

    double best_w = 0.0, best_diff = INFINITY;
    int n_scan = 200;

    for (int i = 0; i < n_scan; i++) {
        double frac = (double)i / (double)(n_scan - 1);
        double w = w_min * pow(w_max / w_min, frac);
        double mag = lead_system_magnitude(sys, w);
        double mag_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        double diff = fabs(mag_db);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    /* Bisection refinement */
    if (best_w > 0.0 && best_diff > 0.1) {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_min) lo = w_min;
        if (hi > w_max) hi = w_max;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double mag_mid = lead_system_magnitude(sys, mid);
            if (mag_mid > 1.0) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    return best_w;
}

/**
 * Find phase crossover frequency: angle(KG(jw_pc)) = -180 deg.
 * L4 - Phase crossover is where Nyquist curve crosses negative real axis.
 */
double lead_find_phase_crossover(const lead_system_t *sys, double w_min,
                                  double w_max) {
    if (!sys) return 0.0;

    double best_w = 0.0, best_diff = INFINITY;
    int n_scan = 300;

    for (int i = 0; i < n_scan; i++) {
        double frac = (double)i / (double)(n_scan - 1);
        double w = w_min * pow(w_max / w_min, frac);
        double phase = lead_system_phase(sys, w);
        double phase_deg = phase * 180.0 / M_PI;
        /* Normalize to [-360, 0] range */
        while (phase_deg > 0.0) phase_deg -= 360.0;
        while (phase_deg < -360.0) phase_deg += 360.0;
        double diff = fabs(phase_deg + 180.0);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    /* Bisection refinement */
    if (best_w > 0.0 && best_diff > 0.5) {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_min) lo = w_min;
        if (hi > w_max) hi = w_max;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double phase = lead_system_phase(sys, mid);
            double phase_deg = phase * 180.0 / M_PI;
            while (phase_deg > 0.0) phase_deg -= 360.0;
            while (phase_deg < -360.0) phase_deg += 360.0;
            if (phase_deg > -180.0) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    return best_w;
}

/**
 * Compute phase margin: PM = 180 + angle(KG(jw_gc)) in degrees.
 *
 * L4 - Positive PM => stable (for minimum-phase systems).
 * PM > 30 deg typically required for adequate damping.
 * PM > 45 deg preferred for robust design.
 * PM > 60 deg for low overshoot.
 */
double lead_compute_phase_margin(const lead_system_t *sys) {
    if (!sys) return 0.0;

    double w_gc = lead_find_gain_crossover(sys, LEAD_FREQ_MIN_DEFAULT,
                                            LEAD_FREQ_MAX_DEFAULT);
    if (w_gc <= 0.0) return 0.0;

    double phase_rad = lead_system_phase(sys, w_gc);
    double phase_deg = phase_rad * 180.0 / M_PI;

    /* PM = 180 + phase(G(jw_gc)) */
    double pm = 180.0 + phase_deg;

    /* Normalize to [-180, 180] */
    while (pm > 180.0) pm -= 360.0;
    while (pm < -180.0) pm += 360.0;

    return pm;
}

/**
 * Compute gain margin: GM = -20*log10(|KG(jw_pc)|) in dB.
 *
 * L4 - Positive GM => stable (for minimum-phase systems).
 * GM > 6 dB typically required.
 */
double lead_compute_gain_margin(const lead_system_t *sys) {
    if (!sys) return INFINITY;

    double w_pc = lead_find_phase_crossover(sys, LEAD_FREQ_MIN_DEFAULT,
                                             LEAD_FREQ_MAX_DEFAULT);
    if (w_pc <= 0.0) return INFINITY;

    double mag = lead_system_magnitude(sys, w_pc);
    if (mag < 1e-300) return INFINITY;

    double gm_db = -20.0 * log10(mag);
    return gm_db;
}

/**
 * Phase margin of compensated system C(s)*G(s).
 * L4 - Verifies that lead compensator achieves the target PM.
 */
double lead_compensated_phase_margin(const lead_compensator_t *comp,
                                      const lead_system_t *sys) {
    if (!comp || !sys) return 0.0;

    /* Find new gain crossover: |C(jw)*G(jw)| = 1 */
    double w_lo = LEAD_FREQ_MIN_DEFAULT, w_hi = LEAD_FREQ_MAX_DEFAULT;
    double best_w = 0.0, best_diff = INFINITY;

    for (int i = 0; i < 300; i++) {
        double frac = (double)i / 299.0;
        double w = w_lo * pow(w_hi / w_lo, frac);
        double mag = lead_compensated_magnitude(comp, sys, w);
        double mag_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        double diff = fabs(mag_db);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    if (best_w <= 0.0) return 0.0;

    /* Refine with bisection */
    {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_lo) lo = w_lo;
        if (hi > w_hi) hi = w_hi;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double mag_mid = lead_compensated_magnitude(comp, sys, mid);
            if (mag_mid > 1.0) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    double phase_rad = lead_compensated_phase(comp, sys, best_w);
    double phase_deg = phase_rad * 180.0 / M_PI;
    double pm = 180.0 + phase_deg;
    while (pm > 180.0) pm -= 360.0;
    while (pm < -180.0) pm += 360.0;

    return pm;
}

/**
 * Gain margin of compensated system C(s)*G(s).
 */
double lead_compensated_gain_margin(const lead_compensator_t *comp,
                                     const lead_system_t *sys) {
    if (!comp || !sys) return INFINITY;

    /* Find phase crossover for compensated system */
    double w_lo = LEAD_FREQ_MIN_DEFAULT, w_hi = LEAD_FREQ_MAX_DEFAULT;
    double best_w = 0.0, best_diff = INFINITY;

    for (int i = 0; i < 300; i++) {
        double frac = i / 299.0;
        double w = w_lo * pow(w_hi / w_lo, frac);
        double phase = lead_compensated_phase(comp, sys, w);
        double phase_deg = phase * 180.0 / M_PI;
        while (phase_deg > 0.0) phase_deg -= 360.0;
        while (phase_deg < -360.0) phase_deg += 360.0;
        double diff = fabs(phase_deg + 180.0);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    if (best_w <= 0.0 || best_diff > 5.0) return INFINITY;

    /* Bisection refinement */
    {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_lo) lo = w_lo;
        if (hi > w_hi) hi = w_hi;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double phase = lead_compensated_phase(comp, sys, mid);
            double phase_deg = phase * 180.0 / M_PI;
            while (phase_deg > 0.0) phase_deg -= 360.0;
            while (phase_deg < -360.0) phase_deg += 360.0;
            if (phase_deg > -180.0) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    double mag = lead_compensated_magnitude(comp, sys, best_w);
    if (mag < 1e-300) return INFINITY;
    return -20.0 * log10(mag);
}

/** Classify phase margin status */
lead_pm_status_t lead_classify_pm(double phase_margin_deg) {
    if (phase_margin_deg > 5.0) return LEAD_PM_STABLE;
    if (phase_margin_deg > -5.0 && phase_margin_deg <= 5.0)
        return LEAD_PM_MARGINAL;
    if (isinf(phase_margin_deg) || isnan(phase_margin_deg))
        return LEAD_PM_UNDEFINED;
    return LEAD_PM_UNSTABLE;
}

/* =========================================================================
 * L4 - Nyquist Stability Criterion
 * ========================================================================= */

/**
 * Count encirclements of -1 point by Nyquist contour.
 *
 * L4 - Nyquist criterion: N = Z - P
 *   N: net CCW encirclements of -1 by Nyquist plot
 *   Z: number of closed-loop RHP poles
 *   P: number of open-loop RHP poles
 *
 * Stable closed-loop system => Z = 0 => N = -P (P CW encirclements)
 * This function evaluates L(jw) = C(jw)*G(jw) along the positive
 * imaginary axis and counts crossings of the negative real axis
 * to the left of -1.
 */
int lead_nyquist_encirclements(const lead_compensator_t *comp,
                                const lead_system_t *sys, int num_points) {
    if (!comp || !sys || num_points < 10) return 0;

    int encirclements = 0;
    double prev_re = 0.0, prev_im = 0.0;
    bool prev_valid = false;

    double log_wmin = log10(LEAD_FREQ_MIN_DEFAULT);
    double log_wmax = log10(LEAD_FREQ_MAX_DEFAULT);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_wmin + frac * (log_wmax - log_wmin));

        double mag = lead_compensated_magnitude(comp, sys, w);
        double phase = lead_compensated_phase(comp, sys, w);
        double re = mag * cos(phase);
        double im = mag * sin(phase);

        if (prev_valid) {
            /* Check crossing of negative real axis */
            if (prev_im * im < 0.0) {
                /* Interpolate crossing point */
                double t = -prev_im / (im - prev_im);
                double cross_re = prev_re + t * (re - prev_re);
                /* If crossing is to the left of -1, count it */
                if (cross_re < -1.0) {
                    /* Determine direction: upward (+1) or downward (-1) */
                    if (prev_im < 0.0 && im > 0.0) {
                        encirclements++; /* CCW */
                    } else if (prev_im > 0.0 && im < 0.0) {
                        encirclements--; /* CW */
                    }
                }
            }
        }

        prev_re = re;
        prev_im = im;
        prev_valid = true;
    }

    return encirclements;
}

/**
 * Check closed-loop stability via Nyquist criterion.
 *
 * L4 - For minimum-phase systems (P=0):
 * Stable iff Nyquist plot does NOT encircle -1 (N=0).
 */
bool lead_is_stable_nyquist(const lead_compensator_t *comp,
                             const lead_system_t *sys) {
    if (!comp || !sys) return false;

    int P = lead_count_rhp_poles(sys);
    int N = lead_nyquist_encirclements(comp, sys, 500);
    int Z = N + P;

    return (Z == 0);
}

/* =========================================================================
 * L4 - Routh-Hurwitz Stability Criterion
 * ========================================================================= */

/**
 * Count open-loop RHP poles using Routh-Hurwitz.
 *
 * L4 - For polynomial a[0] + a[1]*s + ... + a[n]*s^n = 0,
 * construct Routh array and count sign changes in first column.
 */
int lead_routh_hurwitz(const double *coeffs, int order) {
    if (!coeffs || order <= 0) return 0;
    if (order > LEAD_MAX_ORDER) order = LEAD_MAX_ORDER;

    /* Build Routh array */
    int rows = order + 1;
    /* Dynamically allocate with VLA-like approach using static max */
    double routh[LEAD_MAX_ORDER + 1][LEAD_MAX_ORDER + 1];
    memset(routh, 0, sizeof(routh));

    /* First two rows */
    for (int j = 0; j <= order / 2; j++) {
        int idx = order - 2 * j;
        routh[0][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }
    for (int j = 0; j <= (order - 1) / 2; j++) {
        int idx = order - 1 - 2 * j;
        routh[1][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }

    /* Compute remaining rows */
    for (int i = 2; i < rows; i++) {
        for (int j = 0; j < (rows + 1) / 2; j++) {
            double a = routh[i - 2][0];
            double b = routh[i - 2][j + 1];
            double c = routh[i - 1][0];
            double d = routh[i - 1][j + 1];

            if (fabs(c) < 1e-300) {
                /* Zero in first column: epsilon method */
                c = 1e-10;
            }

            routh[i][j] = (c * b - a * d) / c;
        }
    }

    /* Count sign changes in first column */
    int sign_changes = 0;
    double prev = routh[0][0];
    for (int i = 1; i < rows; i++) {
        double curr = routh[i][0];
        /* Treat very small values as zero */
        if (fabs(curr) < 1e-12) curr = 0.0;

        if (prev * curr < 0.0) sign_changes++;
        if (fabs(curr) > 1e-12) prev = curr;
    }

    return sign_changes;
}

/**
 * Count open-loop RHP poles of compensated system.
 * Uses Routh-Hurwitz on the denominator polynomial.
 */
int lead_count_rhp_poles(const lead_system_t *sys) {
    if (!sys) return 0;
    return lead_routh_hurwitz(sys->tf.den, sys->tf.den_order);
}

/**
 * Compute closed-loop characteristic polynomial:
 *   1 + C(s)*G(s) = 0  =>  den_C * den_G + num_C * num_G = 0
 *
 * L4 - Characteristic polynomial determines closed-loop stability.
 *
 * C(s) = num_C / den_C = K_c*(s+z_c)/(s+p_c)
 *   num_C = [K_c*z_c, K_c], den_C = [p_c, 1]
 *
 * G(s) = num_G / den_G = gain * num(s) / den(s)
 *
 * den_CL = den_C * den_G + num_C * num_G
 *
 * This uses discrete convolution for polynomial multiplication.
 */
void lead_closed_loop_polynomial(const lead_compensator_t *comp,
                                  const lead_system_t *sys,
                                  double *cl_poly, int *cl_order) {
    if (!comp || !sys || !cl_poly || !cl_order) return;

    /* Get compensator numerator and denominator */
    lead_tf_t comp_tf;
    lead_to_transfer_function(comp, &comp_tf);

    /* Polynomial multiplication: den_C * den_G */
    int order_dc_dg = comp_tf.den_order + sys->tf.den_order;
    double den_prod[2 * LEAD_MAX_ORDER + 2];
    memset(den_prod, 0, sizeof(den_prod));
    for (int i = 0; i <= comp_tf.den_order; i++) {
        for (int j = 0; j <= sys->tf.den_order; j++) {
            den_prod[i + j] += comp_tf.den[i] * sys->tf.den[j];
        }
    }

    /* Polynomial multiplication: num_C * num_G */
    int order_nc_ng = comp_tf.num_order + sys->tf.num_order;
    double num_prod[2 * LEAD_MAX_ORDER + 2];
    memset(num_prod, 0, sizeof(num_prod));
    for (int i = 0; i <= comp_tf.num_order; i++) {
        for (int j = 0; j <= sys->tf.num_order; j++) {
            num_prod[i + j] += comp_tf.num[i] * sys->tf.num[j] * sys->tf.gain;
        }
    }

    /* Sum: den_CL = den_prod + num_prod */
    int max_order = (order_dc_dg > order_nc_ng) ? order_dc_dg : order_nc_ng;
    for (int i = 0; i <= max_order; i++) {
        cl_poly[i] = den_prod[i] + num_prod[i];
    }
    *cl_order = max_order;
}
