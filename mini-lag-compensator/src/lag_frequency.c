/**
 * @file lag_frequency.c
 * @brief Frequency-domain analysis: Bode plots, Nyquist plots, stability margins
 *
 * L4: Nyquist stability criterion, gain/phase margin theorems
 * L5: Bode plot computation, frequency response evaluation, bandwidth
 *
 * Course: MIT 6.302, Stanford ENGR105, Berkeley ME232,
 *         Caltech CDS 110, ETH 151-0591, Cambridge 3F2, Tsinghua
 *
 * Textbook: Ogata Ch. 8; Franklin/Powell/Emami-Naeini Ch. 6
 */

#include "lag_frequency.h"
#include "lag_compensator.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==========================================================================
 * Helper: Evaluate transfer function magnitude and phase at frequency w
 * ========================================================================== */

/**
 * Evaluate a general rational transfer function G(s) at s = j*w.
 *
 * G(jw) = N(jw) / D(jw)
 *       = Sum_{k=0}^n n_k*(jw)^k / Sum_{k=0}^m d_k*(jw)^k
 *
 * This function handles arbitrary-order TFs by evaluating
 * numerator and denominator polynomials as complex numbers.
 *
 * Complexity: O(max(order_num, order_den))
 */
static void eval_tf_at_freq(const LagTransferFunction *tf, double omega,
                            double *mag, double *phase_rad,
                            double *real_part, double *imag_part) {
    /* Evaluate N(jw) and D(jw) as complex numbers */
    double N_re = 0.0, N_im = 0.0;
    double D_re = 0.0, D_im = 0.0;

    /* Numerator: N(jw) = Sum n_k * (jw)^k */
    for (int k = 0; k <= tf->numerator.order; k++) {
        double coeff = tf->numerator.coeff[k];
        /* (jw)^k = w^k * j^k */
        int mod4 = k % 4;
        double wk = pow(omega, k);
        switch (mod4) {
            case 0: N_re += coeff * wk; break;             /* 1 */
            case 1: N_im += coeff * wk; break;             /* j */
            case 2: N_re -= coeff * wk; break;             /* -1 */
            case 3: N_im -= coeff * wk; break;             /* -j */
        }
    }

    /* Denominator: D(jw) = Sum d_k * (jw)^k */
    for (int k = 0; k <= tf->denominator.order; k++) {
        double coeff = tf->denominator.coeff[k];
        int mod4 = k % 4;
        double wk = pow(omega, k);
        switch (mod4) {
            case 0: D_re += coeff * wk; break;
            case 1: D_im += coeff * wk; break;
            case 2: D_re -= coeff * wk; break;
            case 3: D_im -= coeff * wk; break;
        }
    }

    /* G(jw) = N(jw) / D(jw) */
    double den_mag2 = D_re * D_re + D_im * D_im;
    if (den_mag2 < 1e-300) {
        /* Pole on the imaginary axis */
        *real_part = INFINITY;
        *imag_part = INFINITY;
        *mag = INFINITY;
        *phase_rad = 0.0;
        return;
    }

    *real_part = (N_re * D_re + N_im * D_im) / den_mag2;
    *imag_part = (N_im * D_re - N_re * D_im) / den_mag2;
    *mag = sqrt((*real_part) * (*real_part) + (*imag_part) * (*imag_part));
    *phase_rad = atan2(*imag_part, *real_part);
}

/* ==========================================================================
 * L5: Bode plot computation
 * ========================================================================== */

int lag_compute_bode(const LagTransferFunction *tf, int n_points,
                     LagBodeData *bode) {
    /* Compute Bode plot over a logarithmically spaced frequency range.
     *
     * The frequency range is automatically determined:
     * - Start: 2 decades below the lowest corner frequency
     * - End:   2 decades above the highest corner frequency
     *
     * With logarithmic spacing, each point is:
     *   w[i] = w_min * (w_max/w_min)^(i/(n-1))
     */
    if (!tf || !bode || n_points < 10) return -1;

    /* Determine frequency range from TF poles/zeros.
     * For a general TF, find all corner frequencies.
     * Simplified: use 0.01 to 100 rad/s as default, adjusted by TF parameters. */
    double w_min = 0.01;
    double w_max = 100.0;

    /* If we have coefficients, estimate corner frequencies from the
     * highest and lowest non-zero coefficient ratios. */
    if (tf->denominator.order >= 1 && tf->denominator.coeff) {
        /* The smallest corner frequency is approximately
         * the ratio of smallest to next coefficient. */
        if (tf->denominator.coeff[0] != 0 && tf->denominator.coeff[1] != 0) {
            double wc_low = fabs(tf->denominator.coeff[0] /
                                  tf->denominator.coeff[1]);
            if (wc_low > 0 && wc_low < w_min) w_min = wc_low / 10.0;
            if (wc_low > w_max/10.0) w_max = wc_low * 100.0;
        }
    }

    /* Allocate frequency points */
    bode->points = (LagFreqPoint*)malloc((size_t)n_points * sizeof(LagFreqPoint));
    if (!bode->points) return -2;

    bode->num_points = n_points;
    bode->freq_min = w_min;
    bode->freq_max = w_max;

    double log_w_min = log10(w_min);
    double log_w_max = log10(w_max);

    /* Initialize crossover tracking */
    bode->gain_crossover = -1.0;
    bode->phase_crossover = -1.0;
    bode->phase_margin = 0.0;
    bode->gain_margin = 0.0;
    bode->dc_gain_db = -INFINITY;

    double prev_mag_db = -INFINITY;
    double prev_phase_deg = 0.0;

    for (int i = 0; i < n_points; i++) {
        double log_w = log_w_min + (log_w_max - log_w_min) * i / (n_points - 1);
        double w = pow(10.0, log_w);

        LagFreqPoint *fp = &bode->points[i];
        fp->omega = w;

        eval_tf_at_freq(tf, w, &fp->magnitude, &fp->phase_rad,
                        &fp->real_part, &fp->imag_part);

        /* Magnitude in dB */
        if (fp->magnitude < 1e-300) {
            fp->magnitude_db = -200.0;  /* effectively -inf */
        } else {
            fp->magnitude_db = 20.0 * log10(fp->magnitude);
        }

        /* Phase in degrees */
        fp->phase_deg = fp->phase_rad * 180.0 / M_PI;

        /* Store DC gain (first point approximation) */
        if (i == 0) {
            bode->dc_gain_db = fp->magnitude_db;
        }

        /* Detect gain crossover: |G(jw)| crosses 0 dB */
        if (i > 0 && prev_mag_db * fp->magnitude_db <= 0.0) {
            /* Linear interpolation for more accurate crossover */
            double t = prev_mag_db / (prev_mag_db - fp->magnitude_db);
            double log_w_gc = log10(bode->points[i-1].omega) * (1.0 - t) +
                              log10(fp->omega) * t;
            bode->gain_crossover = pow(10.0, log_w_gc);
        }

        /* Detect phase crossover: phase crosses -180 degrees */
        if (i > 0) {
            double p_prev = prev_phase_deg;
            double p_curr = fp->phase_deg;
            /* Normalize phase to [-360, 0] range for crossover detection */
            while (p_prev > 0) p_prev -= 360.0;
            while (p_prev < -360.0) p_prev += 360.0;
            while (p_curr > 0) p_curr -= 360.0;
            while (p_curr < -360.0) p_curr += 360.0;

            if ((p_prev + 180.0) * (p_curr + 180.0) <= 0.0) {
                /* Crossed -180 */
                double t = fabs(p_prev + 180.0) /
                           (fabs(p_prev + 180.0) + fabs(p_curr + 180.0));
                double log_w_pc = log10(bode->points[i-1].omega) * (1.0 - t) +
                                  log10(fp->omega) * t;
                bode->phase_crossover = pow(10.0, log_w_pc);
            }
        }

        prev_mag_db = fp->magnitude_db;
        prev_phase_deg = fp->phase_deg;
    }

    /* Compute stability margins from crossover data */
    if (bode->gain_crossover > 0) {
        /* Find phase at gain crossover by interpolation */
        double log_w_gc = log10(bode->gain_crossover);
        for (int i = 1; i < n_points; i++) {
            double log_w_prev = log10(bode->points[i-1].omega);
            double log_w_curr = log10(bode->points[i].omega);
            if (log_w_gc >= log_w_prev && log_w_gc <= log_w_curr) {
                double t = (log_w_gc - log_w_prev) / (log_w_curr - log_w_prev);
                double phase_at_gc = bode->points[i-1].phase_deg * (1.0 - t) +
                                     bode->points[i].phase_deg * t;
                bode->phase_margin = 180.0 + phase_at_gc;
                break;
            }
        }
    }

    if (bode->phase_crossover > 0) {
        /* Find magnitude at phase crossover */
        double log_w_pc = log10(bode->phase_crossover);
        for (int i = 1; i < n_points; i++) {
            double log_w_prev = log10(bode->points[i-1].omega);
            double log_w_curr = log10(bode->points[i].omega);
            if (log_w_pc >= log_w_prev && log_w_pc <= log_w_curr) {
                double t = (log_w_pc - log_w_prev) / (log_w_curr - log_w_prev);
                double mag_db_at_pc = bode->points[i-1].magnitude_db * (1.0 - t) +
                                      bode->points[i].magnitude_db * t;
                bode->gain_margin = -mag_db_at_pc;
                break;
            }
        }
    }

    return 0;
}

/* ==========================================================================
 * Helper: Compute open-loop transfer function G_c(s)*G(s)
 * ========================================================================== */

/**
 * Build the open-loop transfer function L(s) = G_c(s) * G(s).
 * For a first-order compensator and general plant:
 *   L(s) = [Kc*(T*s+1)/(beta*T*s+1)] * [N_p(s)/D_p(s)]
 *
 * This function allocates memory for the result; caller must free.
 *
 * Complexity: O(order_num * order_den)
 */
static int build_open_loop_tf(const LagCompensator *lag,
                               const LagTransferFunction *plant,
                               LagTransferFunction *L) {
    /* L(s) = G_c(s) * G(s) = (Kc*(Ts+1)) / (beta*Ts+1) * Np(s)/Dp(s)
     * Numerator:   Kc*(Ts+1) * Np(s)
     * Denominator: (beta*Ts+1) * Dp(s)
     *
     * G_c numerator: Kc*T*s + Kc = n_c1*s + n_c0
     * G_c denominator: beta*T*s + 1 = d_c1*s + d_c0
     *
     * Product numerator: (n_c1*s + n_c0) * Np(s)
     * Product denominator: (d_c1*s + d_c0) * Dp(s)
     */
    int np_order = plant->numerator.order;
    int dp_order = plant->denominator.order;

    /* Compensator numerator: order 1, coeffs [Kc, Kc*T] */
    double nc[2] = {lag->Kc, lag->Kc * lag->T};
    /* Compensator denominator: order 1, coeffs [1, beta*T] */
    double dc[2] = {1.0, lag->beta * lag->T};

    /* Convolve polynomials */
    int L_num_order = np_order + 1;
    int L_den_order = dp_order + 1;

    L->numerator.order = L_num_order;
    L->numerator.coeff = (double*)calloc((size_t)(L_num_order + 1), sizeof(double));
    L->denominator.order = L_den_order;
    L->denominator.coeff = (double*)calloc((size_t)(L_den_order + 1), sizeof(double));

    if (!L->numerator.coeff || !L->denominator.coeff) {
        free(L->numerator.coeff);
        free(L->denominator.coeff);
        return -1;
    }

    /* Numerator convolution: nc * Np */
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= np_order; j++) {
            L->numerator.coeff[i + j] += nc[i] * plant->numerator.coeff[j];
        }
    }

    /* Denominator convolution: dc * Dp */
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= dp_order; j++) {
            L->denominator.coeff[i + j] += dc[i] * plant->denominator.coeff[j];
        }
    }

    /* DC gain */
    if (L->denominator.coeff[0] != 0.0) {
        L->dc_gain = L->numerator.coeff[0] / L->denominator.coeff[0];
    } else {
        L->dc_gain = INFINITY;
    }

    return 0;
}

/* ==========================================================================
 * L5: Open-loop and closed-loop Bode
 * ========================================================================== */

int lag_compute_open_loop_bode(const LagCompensator *lag,
                                const LagTransferFunction *plant,
                                int n_points, LagBodeData *bode) {
    /* Compute Bode plot for L(s) = G_c(s) * G(s). */
    if (!lag || !plant || !bode) return -1;

    LagTransferFunction L;
    int ret = build_open_loop_tf(lag, plant, &L);
    if (ret != 0) return ret;

    ret = lag_compute_bode(&L, n_points, bode);

    /* Clean up */
    free(L.numerator.coeff);
    free(L.denominator.coeff);

    return ret;
}

int lag_compute_closed_loop_bode(const LagCompensator *lag,
                                  const LagTransferFunction *plant,
                                  int n_points, LagBodeData *bode) {
    /* T(s) = L(s) / (1 + L(s))
     *
     * Since computing the closed-loop TF exactly requires polynomial
     * addition and division, we evaluate it pointwise:
     *   T(jw) = L(jw) / (1 + L(jw))
     *
     * We first build L(s), then evaluate T(jw) at each frequency.
     */
    if (!lag || !plant || !bode || n_points < 10) return -1;

    LagTransferFunction L;
    int ret = build_open_loop_tf(lag, plant, &L);
    if (ret != 0) return ret;

    /* We need to compute Bode for T(s) = L/(1+L).
     * Strategy: compute L(jw) at each frequency, then T = L/(1+L). */

    /* First compute L(jw) data */
    ret = lag_compute_bode(&L, n_points, bode);
    if (ret != 0) {
        free(L.numerator.coeff);
        free(L.denominator.coeff);
        return ret;
    }

    /* Transform each point: T = L / (1+L) */
    for (int i = 0; i < bode->num_points; i++) {
        LagFreqPoint *fp = &bode->points[i];
        double L_re = fp->real_part;
        double L_im = fp->imag_part;

        /* T = L / (1+L) = (L_re + j*L_im) / ((1+L_re) + j*L_im) */
        double den_re = 1.0 + L_re;
        double den_im = L_im;
        double den_mag2 = den_re * den_re + den_im * den_im;

        if (den_mag2 < 1e-300) {
            fp->real_part = INFINITY;
            fp->imag_part = INFINITY;
            fp->magnitude = INFINITY;
            fp->magnitude_db = INFINITY;
            fp->phase_rad = 0.0;
            fp->phase_deg = 0.0;
        } else {
            fp->real_part = (L_re * den_re + L_im * den_im) / den_mag2;
            fp->imag_part = (L_im * den_re - L_re * den_im) / den_mag2;
            fp->magnitude = sqrt(fp->real_part * fp->real_part +
                                 fp->imag_part * fp->imag_part);
            fp->magnitude_db = 20.0 * log10(fp->magnitude);
            fp->phase_rad = atan2(fp->imag_part, fp->real_part);
            fp->phase_deg = fp->phase_rad * 180.0 / M_PI;
        }
    }

    free(L.numerator.coeff);
    free(L.denominator.coeff);
    return 0;
}

int lag_compute_sensitivity_bode(const LagCompensator *lag,
                                  const LagTransferFunction *plant,
                                  int n_points, LagBodeData *bode) {
    /* S(s) = 1 / (1 + L(s))
     *
     * S(jw) = 1 / (1 + L(jw))
     */
    if (!lag || !plant || !bode || n_points < 10) return -1;

    /* We use the same approach as closed-loop but with S = 1/(1+L) */
    LagTransferFunction L;
    int ret = build_open_loop_tf(lag, plant, &L);
    if (ret != 0) return ret;

    ret = lag_compute_bode(&L, n_points, bode);
    if (ret != 0) {
        free(L.numerator.coeff);
        free(L.denominator.coeff);
        return ret;
    }

    /* Transform: S = 1/(1+L) */
    for (int i = 0; i < bode->num_points; i++) {
        LagFreqPoint *fp = &bode->points[i];
        double L_re = fp->real_part;
        double L_im = fp->imag_part;

        /* S = 1 / (1+L) = 1 / ((1+L_re) + j*L_im) */
        double den_re = 1.0 + L_re;
        double den_im = L_im;
        double den_mag2 = den_re * den_re + den_im * den_im;

        if (den_mag2 < 1e-300) {
            fp->real_part = INFINITY;
            fp->imag_part = INFINITY;
            fp->magnitude = INFINITY;
            fp->magnitude_db = INFINITY;
            fp->phase_rad = 0.0;
            fp->phase_deg = 0.0;
        } else {
            fp->real_part = den_re / den_mag2;
            fp->imag_part = -den_im / den_mag2;
            fp->magnitude = 1.0 / sqrt(den_mag2);
            fp->magnitude_db = -10.0 * log10(den_mag2);
            fp->phase_rad = atan2(fp->imag_part, fp->real_part);
            fp->phase_deg = fp->phase_rad * 180.0 / M_PI;
        }
    }

    free(L.numerator.coeff);
    free(L.denominator.coeff);
    return 0;
}

/* ==========================================================================
 * L5: Nyquist plot computation
 * ========================================================================== */

int lag_compute_nyquist(const LagTransferFunction *tf, int n_points,
                        int P, LagNyquistData *nyquist) {
    /* Compute Nyquist contour data.
     *
     * The Nyquist D-contour traverses:
     *   1. Positive imaginary axis: s = jw for w in (0, inf)
     *   2. Infinite semicircle in RHP (mapped to origin for proper TFs)
     *   3. Negative imaginary axis: s = -jw for w in (inf, 0)
     *
     * We evaluate G(jw) at logarithmically spaced frequencies.
     * The conjugate symmetry gives G(-jw) = conj(G(jw)).
     *
     * Encirclement counting: We use the winding number around -1.
     * N = (1/2pi) * delta_arg(G(jw) + 1) as w goes from 0 to inf.
     *
     * Closed-loop stability: Z = N + P, stable if Z = 0.
     */
    if (!tf || !nyquist || n_points < 20) return -1;

    /* Allocate points for the full contour (positive + negative frequencies) */
    int total_points = 2 * n_points;
    nyquist->points = (LagFreqPoint*)malloc((size_t)total_points *
                                             sizeof(LagFreqPoint));
    if (!nyquist->points) return -2;

    nyquist->num_points = total_points;

    /* Determine frequency range */
    double w_min = 0.01, w_max = 100.0;
    double log_w_min = log10(w_min);
    double log_w_max = log10(w_max);

    /* Evaluate for positive frequencies */
    for (int i = 0; i < n_points; i++) {
        double log_w = log_w_min + (log_w_max - log_w_min) * i / (n_points - 1);
        double w = pow(10.0, log_w);

        LagFreqPoint *fp = &nyquist->points[i];
        fp->omega = w;

        eval_tf_at_freq(tf, w, &fp->magnitude, &fp->phase_rad,
                        &fp->real_part, &fp->imag_part);
        fp->magnitude_db = 20.0 * log10(fp->magnitude > 1e-300 ? fp->magnitude : 1e-300);
        fp->phase_deg = fp->phase_rad * 180.0 / M_PI;
    }

    /* Mirror for negative frequencies using conjugate symmetry:
     * G(-jw) = conj(G(jw)) */
    for (int i = 0; i < n_points; i++) {
        int pos_idx = n_points - 1 - i;  /* reverse order */
        int neg_idx = n_points + i;

        nyquist->points[neg_idx].omega = -nyquist->points[pos_idx].omega;
        nyquist->points[neg_idx].real_part = nyquist->points[pos_idx].real_part;
        nyquist->points[neg_idx].imag_part = -nyquist->points[pos_idx].imag_part;
        nyquist->points[neg_idx].magnitude = nyquist->points[pos_idx].magnitude;
        nyquist->points[neg_idx].magnitude_db = nyquist->points[pos_idx].magnitude_db;
        nyquist->points[neg_idx].phase_rad = -nyquist->points[pos_idx].phase_rad;
        nyquist->points[neg_idx].phase_deg = -nyquist->points[pos_idx].phase_deg;
    }

    /* Count encirclements of -1 point.
     * The Nyquist contour encircles -1 if the phase of G(jw) + 1
     * changes by multiples of 2*pi as w traverses the contour.
     *
     * Simplified: check if curve crosses negative real axis to the left of -1. */
    int crossings = 0;
    for (int i = 0; i < n_points - 1; i++) {
        double x1 = nyquist->points[i].real_part;
        double y1 = nyquist->points[i].imag_part;
        double x2 = nyquist->points[i+1].real_part;
        double y2 = nyquist->points[i+1].imag_part;

        /* Check if segment crosses negative real axis left of -1 */
        if (y1 * y2 < 0) {  /* opposite signs -> crosses real axis */
            /* Interpolate crossing point */
            double t = -y1 / (y2 - y1);
            double x_cross = x1 + t * (x2 - x1);
            if (x_cross < -1.0 && y2 > y1) {
                crossings++;  /* clockwise crossing */
            } else if (x_cross < -1.0 && y2 < y1) {
                crossings--;  /* counter-clockwise */
            }
        }
    }

    nyquist->encirclements = crossings;
    nyquist->is_stable = (P + nyquist->encirclements == 0) ? 1 : 0;

    return 0;
}

/* ==========================================================================
 * L4: Stability margin computation
 * ========================================================================== */

double lag_find_gain_crossover(const LagBodeData *bode) {
    if (!bode || bode->num_points < 2) return -1.0;

    /* Find frequency where magnitude crosses 0 dB.
     * Look for sign change in magnitude_db. */
    for (int i = 1; i < bode->num_points; i++) {
        double mag1 = bode->points[i-1].magnitude_db;
        double mag2 = bode->points[i].magnitude_db;
        if (mag1 * mag2 <= 0.0 && mag1 > mag2) {
            /* Interpolate */
            double t = mag1 / (mag1 - mag2);
            double log_w1 = log10(bode->points[i-1].omega);
            double log_w2 = log10(bode->points[i].omega);
            double log_w = log_w1 + t * (log_w2 - log_w1);
            return pow(10.0, log_w);
        }
    }
    return -1.0;
}

double lag_find_phase_crossover(const LagBodeData *bode) {
    if (!bode || bode->num_points < 2) return -1.0;

    for (int i = 1; i < bode->num_points; i++) {
        double p1 = bode->points[i-1].phase_deg;
        double p2 = bode->points[i].phase_deg;

        /* Normalize */
        while (p1 > 0) p1 -= 360.0;
        while (p1 < -360.0) p1 += 360.0;
        while (p2 > 0) p2 -= 360.0;
        while (p2 < -360.0) p2 += 360.0;

        if ((p1 + 180.0) * (p2 + 180.0) <= 0.0 && p2 < p1) {
            double t = fabs(p1 + 180.0) / (fabs(p1 + 180.0) + fabs(p2 + 180.0));
            double log_w1 = log10(bode->points[i-1].omega);
            double log_w2 = log10(bode->points[i].omega);
            double log_w = log_w1 + t * (log_w2 - log_w1);
            return pow(10.0, log_w);
        }
    }
    return -1.0;
}

double lag_compute_phase_margin(const LagBodeData *bode) {
    double w_gc = lag_find_gain_crossover(bode);
    if (w_gc < 0) return 0.0;

    /* Find phase at w_gc by interpolation */
    double log_w_gc = log10(w_gc);
    for (int i = 1; i < bode->num_points; i++) {
        double log_w1 = log10(bode->points[i-1].omega);
        double log_w2 = log10(bode->points[i].omega);
        if (log_w_gc >= log_w1 && log_w_gc <= log_w2) {
            double t = (log_w_gc - log_w1) / (log_w2 - log_w1);
            double phase = bode->points[i-1].phase_deg * (1.0 - t) +
                           bode->points[i].phase_deg * t;
            return 180.0 + phase;
        }
    }
    return 0.0;
}

double lag_compute_gain_margin(const LagBodeData *bode) {
    double w_pc = lag_find_phase_crossover(bode);
    if (w_pc < 0) return INFINITY;

    /* Find magnitude at w_pc by interpolation */
    double log_w_pc = log10(w_pc);
    for (int i = 1; i < bode->num_points; i++) {
        double log_w1 = log10(bode->points[i-1].omega);
        double log_w2 = log10(bode->points[i].omega);
        if (log_w_pc >= log_w1 && log_w_pc <= log_w2) {
            double t = (log_w_pc - log_w1) / (log_w2 - log_w1);
            double mag_db = bode->points[i-1].magnitude_db * (1.0 - t) +
                            bode->points[i].magnitude_db * t;
            return -mag_db;
        }
    }
    return INFINITY;
}

void lag_compute_stability_margins(const LagBodeData *bode,
                                    double *phase_margin_deg,
                                    double *gain_margin_db) {
    *phase_margin_deg = lag_compute_phase_margin(bode);
    *gain_margin_db = lag_compute_gain_margin(bode);
}

/* ==========================================================================
 * L2: Bandwidth and performance metrics
 * ========================================================================== */

double lag_compute_bandwidth(const LagBodeData *bode) {
    /* Find frequency where closed-loop magnitude = -3 dB (0.707).
     *
     * For T(s) = L/(1+L), the bandwidth is where |T| = -3 dB.
     * This is approximately the same as the 0 dB crossover of L(s)
     * for well-behaved systems.
     *
     * We search the Bode data for the -3 dB point. */
    if (!bode || bode->num_points < 2) return -1.0;

    /* If this is closed-loop data, find where magnitude = -3 dB */
    double target_db = -3.0;
    for (int i = 1; i < bode->num_points; i++) {
        double mag1 = bode->points[i-1].magnitude_db;
        double mag2 = bode->points[i].magnitude_db;
        if ((mag1 - target_db) * (mag2 - target_db) <= 0.0) {
            double t = (mag1 - target_db) / (mag1 - mag2);
            double log_w1 = log10(bode->points[i-1].omega);
            double log_w2 = log10(bode->points[i].omega);
            double log_w = log_w1 + t * (log_w2 - log_w1);
            return pow(10.0, log_w);
        }
    }
    return -1.0;
}

double lag_compute_steady_state_error(const LagBodeData *bode,
                                       LagESSType ess_type) {
    /* Compute steady-state error from low-frequency Bode data.
     *
     * Using Final Value Theorem and low-frequency asymptote.
     *
     * For the first data point (lowest frequency), we estimate
     * the error constant from the magnitude. */
    if (!bode || bode->num_points < 1) return INFINITY;

    /* Use lowest frequency point to estimate DC behavior */
    double mag_lf = bode->points[0].magnitude;

    double ess;
    switch (ess_type) {
        case LAG_ESS_STEP:
            ess = 1.0 / (1.0 + mag_lf);
            break;
        case LAG_ESS_RAMP:
            /* For type-0 systems, Kv = 0 so ess = inf.
             * For type-1, Kv = w*|G| at low frequencies.
             * Estimate from slope if possible. */
            if (bode->num_points >= 2) {
                double mag1 = bode->points[0].magnitude;
                double mag2 = bode->points[1].magnitude;
                double w1 = bode->points[0].omega;
                double w2 = bode->points[1].omega;
                double slope = (20.0*log10(mag2) - 20.0*log10(mag1)) /
                               (log10(w2) - log10(w1));
                if (fabs(slope + 20.0) < 5.0) {
                    /* -20 dB/dec => type 1 system */
                    double Kv_est = w1 * mag1;
                    ess = 1.0 / Kv_est;
                } else {
                    ess = INFINITY;
                }
            } else {
                ess = 1.0 / mag_lf;
            }
            break;
        case LAG_ESS_PARABOLIC:
            /* Similar slope analysis for type-2 systems */
            ess = INFINITY;  /* Most systems need type >= 2 */
            break;
        default:
            ess = INFINITY;
    }
    return ess;
}

double lag_estimate_error_constant(const LagBodeData *bode,
                                    LagESSType ess_type) {
    if (!bode || bode->num_points < 1) return 0.0;

    double mag_lf = bode->points[0].magnitude;
    double w_lf = bode->points[0].omega;

    switch (ess_type) {
        case LAG_ESS_STEP:  return mag_lf;
        case LAG_ESS_RAMP:  return w_lf * mag_lf;
        case LAG_ESS_PARABOLIC: return w_lf * w_lf * mag_lf;
        default: return 0.0;
    }
}

double lag_compute_resonant_peak(const LagBodeData *closed_loop_bode) {
    if (!closed_loop_bode || closed_loop_bode->num_points < 1) return 0.0;

    double M_p = 0.0;
    for (int i = 0; i < closed_loop_bode->num_points; i++) {
        if (closed_loop_bode->points[i].magnitude > M_p) {
            M_p = closed_loop_bode->points[i].magnitude;
        }
    }
    return M_p;
}

double lag_compute_resonant_frequency(const LagBodeData *closed_loop_bode) {
    if (!closed_loop_bode || closed_loop_bode->num_points < 1) return -1.0;

    double M_p = 0.0;
    double w_r = -1.0;
    for (int i = 0; i < closed_loop_bode->num_points; i++) {
        if (closed_loop_bode->points[i].magnitude > M_p) {
            M_p = closed_loop_bode->points[i].magnitude;
            w_r = closed_loop_bode->points[i].omega;
        }
    }
    return w_r;
}