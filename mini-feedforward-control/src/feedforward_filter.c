/**
 * @file feedforward_filter.c
 * @brief Prefilter, FIR/IIR filter, Butterworth, notch filter implementations.
 *
 * L3-L5: Digital and analog filter design for feedforward signal conditioning.
 * Implements prefiltering for reference smoothing, low-pass filtering
 * for properness enforcement, and notch filtering for resonance suppression.
 */

#include "feedforward_filter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L5: FIR Filter
 * ========================================================================== */

void fir_init(FIRFilter *fir, const double *coeff, int n_taps)
{
    if (!fir || !coeff || n_taps <= 0) return;

    fir->n_taps = n_taps;
    fir->coeff = (double *)malloc(n_taps * sizeof(double));
    fir->buffer = (double *)calloc(n_taps, sizeof(double));
    if (fir->coeff) memcpy(fir->coeff, coeff, n_taps * sizeof(double));
    fir->buf_idx = 0;
}

double fir_process(FIRFilter *fir, double x)
{
    if (!fir || !fir->coeff || !fir->buffer) return x;

    /* Insert new sample into circular buffer */
    fir->buffer[fir->buf_idx] = x;

    /* Compute output: y = sum(b[k] * x[n-k]) */
    double y = 0.0;
    for (int k = 0; k < fir->n_taps; k++) {
        int idx = (fir->buf_idx - k + fir->n_taps) % fir->n_taps;
        y += fir->coeff[k] * fir->buffer[idx];
    }

    /* Advance circular buffer index */
    fir->buf_idx = (fir->buf_idx + 1) % fir->n_taps;

    return y;
}

void fir_reset(FIRFilter *fir)
{
    if (fir && fir->buffer) {
        memset(fir->buffer, 0, fir->n_taps * sizeof(double));
        fir->buf_idx = 0;
    }
}

void fir_free(FIRFilter *fir)
{
    if (fir) {
        free(fir->coeff);
        free(fir->buffer);
        fir->coeff = NULL;
        fir->buffer = NULL;
        fir->n_taps = 0;
    }
}

void fir_design_moving_average(int window_size, FIRFilter *fir)
{
    if (window_size <= 0) window_size = 1;
    double *coeff = (double *)malloc(window_size * sizeof(double));
    if (!coeff) return;

    for (int i = 0; i < window_size; i++) {
        coeff[i] = 1.0 / (double)window_size;
    }

    fir_init(fir, coeff, window_size);
    free(coeff);
}

void fir_design_lowpass(double cutoff, int n_taps, FIRFilter *fir)
{
    if (cutoff <= 0.0 || cutoff >= 0.5 || n_taps <= 0) return;

    /* Design low-pass FIR filter using Hamming window method.
     * Ideal impulse response: h_d[n] = 2*fc * sinc(2*fc*(n-M/2))
     * where sinc(x) = sin(pi*x)/(pi*x), fc = cutoff frequency (normalized).
     * Hamming window: w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1)) */

    double *h = (double *)malloc(n_taps * sizeof(double));
    if (!h) return;

    int M = n_taps - 1;
    double fc = cutoff;

    for (int n = 0; n < n_taps; n++) {
        double t = (double)(n - M / 2);
        if (fabs(t) < 1e-15) {
            h[n] = 2.0 * fc; /* sinc(0) = 1 */
        } else {
            h[n] = 2.0 * fc * sin(2.0 * M_PI * fc * t) / (M_PI * t);
        }
        /* Hamming window */
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * n / M);
        h[n] *= w;
    }

    /* Normalize for unity DC gain */
    double sum = 0.0;
    for (int n = 0; n < n_taps; n++) sum += h[n];
    if (fabs(sum) > 1e-15) {
        for (int n = 0; n < n_taps; n++) h[n] /= sum;
    }

    fir_init(fir, h, n_taps);
    free(h);
}

/* ==========================================================================
 * L5: IIR Filter
 * ========================================================================== */

void iir_init(IIRFilter *iir, const double *b, const double *a, int order)
{
    if (!iir || !b || !a || order <= 0) return;

    iir->order = order;
    iir->b = (double *)malloc((order + 1) * sizeof(double));
    iir->a = (double *)malloc(order * sizeof(double));
    iir->state = (double *)calloc(order + 1, sizeof(double));

    if (iir->b && iir->a) {
        /* Normalize by a[0] if not 1.0; here we assume a[0]=1.0.
         * Store b[0..order] and a[1..order] */
        for (int i = 0; i <= order; i++) iir->b[i] = b[i];
        for (int i = 0; i < order; i++) iir->a[i] = a[i + 1];
    }
}

double iir_process(IIRFilter *iir, double x)
{
    if (!iir || !iir->b || !iir->a || !iir->state) return x;

    /* Direct Form II Transposed:
     * y = b[0]*x + state[0]
     * state[0] = b[1]*x + state[1] - a[0]*y
     * state[1] = b[2]*x + state[2] - a[1]*y
     * ...
     * More precisely:
     *   y[n] = b0*x[n] + d0[n-1]
     *   d0[n] = b1*x[n] - a1*y[n] + d1[n-1]
     *   d1[n] = b2*x[n] - a2*y[n] + d2[n-1]
     *   ...
     */

    double y = iir->b[0] * x + iir->state[0];

    /* Update states */
    for (int i = 0; i < iir->order; i++) {
        double d_new = (i + 1 <= iir->order ? iir->b[i + 1] * x : 0.0)
                       - (i < iir->order ? iir->a[i] * y : 0.0)
                       + (i + 1 < iir->order ? iir->state[i + 1] : 0.0);

        if (i + 1 < iir->order) {
            iir->state[i] = iir->state[i + 1]; /* Shift */
        }
        /* Store */
        iir->state[i] = d_new;
    }

    return y;
}

void iir_reset(IIRFilter *iir)
{
    if (iir && iir->state) {
        memset(iir->state, 0, (iir->order + 1) * sizeof(double));
    }
}

void iir_free(IIRFilter *iir)
{
    if (iir) {
        free(iir->b);
        free(iir->a);
        free(iir->state);
        iir->b = NULL;
        iir->a = NULL;
        iir->state = NULL;
        iir->order = 0;
    }
}

/* ==========================================================================
 * L5: Butterworth Low-Pass Filter Design
 * ========================================================================== */

void butterworth_lp_design(int order, double cutoff, ButterworthLP *lp)
{
    if (!lp || order <= 0 || order > 8) return;

    lp->order = order;
    lp->cutoff = cutoff;

    /* Butterworth polynomial coefficients for orders 1-8.
     * Normalized for cutoff = 1 rad/s.
     * Coefficients: den(s) = s^n + b_{n-1}*s^{n-1} + ... + b1*s + 1
     * where b_k = product-of-cos formula for Butterworth. */

    /* Tabulated Butterworth polynomial coefficients (s^n + ... + 1) */
    static const double b_coeffs[8][8] = {
        {1.000000000000000},                          /* n=1 */
        {1.414213562373095, 1.000000000000000},       /* n=2 */
        {2.000000000000000, 2.000000000000000, 1.000000000000000}, /* n=3 */
        {2.613125929752753, 3.414213562373095, 2.613125929752753, 1.000000000000000}, /* n=4 */
        {3.236067977499790, 5.236067977499790, 5.236067977499790, 3.236067977499790, 1.000000000000000}, /* n=5 */
        {3.863703305156273, 7.464101615137754, 9.141620172685641, 7.464101615137754, 3.863703305156273, 1.000000000000000}, /* n=6 */
        {4.493959207434935, 10.097834679044615, 14.591793886021997, 14.591793886021997, 10.097834679044615, 4.493959207434935, 1.000000000000000}, /* n=7 */
        {5.125830895483816, 13.137071184544074, 21.846150969204733, 25.688355931996508, 21.846150969204733, 13.137071184544074, 5.125830895483816, 1.000000000000000}  /* n=8 */
    };

    lp->den_order = order;
    lp->num[0] = 1.0;
    for (int i = 1; i < 10; i++) lp->num[i] = 0.0;

    /* Denormalize: multiply s^k term by (1/wc)^k
     * den_coeff[k] = b_coeffs[n-1][k] * (1/wc)^(n-k)
     * where k=0 corresponds to constant term, k=n corresponds to s^n */
    lp->den[0] = 1.0; /* constant term in normalized form */
    for (int k = 0; k < order; k++) {
        double exponent = (double)(order - k);
        lp->den[order - k] = b_coeffs[order - 1][k] * pow(1.0 / cutoff, exponent);
    }
    lp->den[0] = 1.0; /* s^0 term = 1 always */

    /* Reorder: den = [a0, a1, ..., an] where a0=1, a1=b1/wc, ..., an=1/wc^n */
    for (int k = 0; k < order; k++) {
        lp->den[k + 1] = b_coeffs[order - 1][k] / pow(cutoff, k + 1);
    }
    lp->den[0] = 1.0;
}

void butterworth_poles(int order, double cutoff,
                       double *poles, int *n_poles)
{
    if (!poles || !n_poles || order <= 0) return;

    /* Butterworth poles: p_k = wc * exp(j*pi*(0.5 + (2k+1)/(2n)))
     * for k = 0, 1, ..., n-1
     * Return real and imaginary parts interleaved. */
    int count = 0;
    for (int k = 0; k < order; k++) {
        double theta = M_PI * (0.5 + (2.0 * k + 1.0) / (2.0 * order));
        double real = cutoff * cos(theta);
        double imag = cutoff * sin(theta);

        poles[count * 2] = real;
        poles[count * 2 + 1] = imag;
        count++;
    }
    *n_poles = count;
}

/* ==========================================================================
 * L5: Notch Filter Design
 * ========================================================================== */

void notch_design(double wn, double depth, double width, NotchFilter *nf)
{
    if (!nf || wn <= 0.0) return;

    nf->wn = wn;

    /* Notch filter: H(s) = (s^2 + 2*zeta_z*wn*s + wn^2) / (s^2 + 2*zeta_p*wn*s + wn^2)
     * Depth at w=wn: |H(j*wn)| = zeta_z / zeta_p
     * Width: related to zeta_p (larger zeta_p = wider notch) */

    if (depth < 0.001) depth = 0.001; /* > 0 for stability */
    if (depth > 0.999) depth = 0.999;

    /* Choose zeta_z to achieve desired depth:
     * depth = zeta_z / zeta_p => zeta_z = depth * zeta_p */
    double zeta_p = width;
    if (zeta_p < 0.05) zeta_p = 0.05;
    if (zeta_p > 2.0) zeta_p = 2.0;

    nf->zeta_z = depth * zeta_p;
    nf->zeta_p = zeta_p;
    nf->depth = depth;
}

void notch_biquad_design(double wn, double zeta, double Ts,
                         double b[3], double a[3])
{
    /* Design discrete-time notch via bilinear transform:
     * Analog notch: H(s) = (s^2 + wn^2) / (s^2 + 2*zeta*wn*s + wn^2)
     * Bilinear: s = (2/Ts) * (z-1)/(z+1)
     *
     * Resulting biquad:
     *   b0 = 1 + (wn*Ts/2)^2
     *   b1 = 2*(1 - (wn*Ts/2)^2)
     *   b2 = b0
     *   a0 = 1 + zeta*wn*Ts + (wn*Ts/2)^2
     *   a1 = 2*(1 - (wn*Ts/2)^2)
     *   a2 = 1 - zeta*wn*Ts + (wn*Ts/2)^2
     */

    if (Ts <= 0.0) { b[0]=1; b[1]=0; b[2]=0; a[0]=1; a[1]=0; a[2]=0; return; }

    double wT = wn * Ts / 2.0;
    double wT2 = wT * wT;

    b[0] = 1.0 + wT2;
    b[1] = 2.0 * (1.0 - wT2);
    b[2] = b[0];

    a[0] = 1.0 + zeta * wn * Ts + wT2;
    a[1] = 2.0 * (1.0 - wT2);
    a[2] = 1.0 - zeta * wn * Ts + wT2;

    /* Normalize so a[0] = 1 */
    double a0_inv = 1.0 / a[0];
    b[0] *= a0_inv;
    b[1] *= a0_inv;
    b[2] *= a0_inv;
    a[1] *= a0_inv;
    a[2] *= a0_inv;
    a[0] = 1.0;
}

/* ==========================================================================
 * L3: Prefilter for Reference Smoothing
 * ========================================================================== */

void prefilter_velocity_limit(double v_max, double dt,
                              const double *r_in, int n,
                              double *r_out)
{
    if (!r_in || !r_out || n <= 0) return;

    r_out[0] = r_in[0];
    for (int i = 1; i < n; i++) {
        double delta = r_in[i] - r_out[i - 1];
        double max_step = v_max * dt;
        if (delta > max_step) {
            r_out[i] = r_out[i - 1] + max_step;
        } else if (delta < -max_step) {
            r_out[i] = r_out[i - 1] - max_step;
        } else {
            r_out[i] = r_in[i];
        }
    }
}

void prefilter_scurve(double v_max, double a_max, double dt,
                      double target, double *r_out, int *n_out)
{
    /* S-curve trajectory generation:
     * Phase 1: constant acceleration (jerk = 0 assumed, or constant jerk)
     * Phase 2: constant velocity
     * Phase 3: constant deceleration
     *
     * Compute time for each phase. */

    if (!r_out || !n_out || target == 0.0 || dt <= 0.0) {
        if (n_out) *n_out = 0;
        return;
    }

    /* Distance required for acceleration to v_max and back:
     * d_accel = v_max^2 / a_max (each of accel and decel)
     * Total distance for full trapezoid = 2 * v_max^2 / a_max + v_max * t_const */
    double sign = (target > 0.0) ? 1.0 : -1.0;
    double dist = fabs(target);
    double t_accel = v_max / a_max;
    double d_accel = 0.5 * a_max * t_accel * t_accel; /* = v_max^2 / (2*a_max) */

    int phase = 0;
    int max_pts = 10000;
    int idx = 0;
    double pos = 0.0, vel = 0.0;

    while (idx < max_pts && fabs(pos) < dist) {
        double remain = dist - fabs(pos);

        if (phase == 0) {
            /* Acceleration phase */
            if (vel < v_max && remain > d_accel) {
                vel += a_max * dt;
                if (vel > v_max) vel = v_max;
            } else {
                phase = 1; /* Switch to constant velocity */
            }
        }

        if (phase == 1) {
            /* Constant velocity phase */
            if (remain <= d_accel) {
                phase = 2; /* Switch to deceleration */
            }
        }

        if (phase == 2) {
            /* Deceleration phase */
            double stop_dist = vel * vel / (2.0 * a_max);
            if (remain <= stop_dist) {
                vel -= a_max * dt;
                if (vel < 0.0) vel = 0.0;
            }
        }

        pos += vel * dt;
        if (fabs(pos) > dist) pos = sign * dist;

        r_out[idx] = sign * pos;
        idx++;

        if (fabs(pos) >= dist - 1e-9 || vel <= 1e-9) break;
    }

    /* Ensure final position is reached */
    if (idx < max_pts) {
        r_out[idx - 1] = target;
    }

    *n_out = idx;
}

void prefilter_exp_smooth(double alpha, const double *r_in,
                          int n, double *r_out)
{
    if (!r_in || !r_out || n <= 0) return;
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    r_out[0] = r_in[0];
    for (int i = 1; i < n; i++) {
        r_out[i] = alpha * r_in[i] + (1.0 - alpha) * r_out[i - 1];
    }
}
