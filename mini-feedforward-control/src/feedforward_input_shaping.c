/**
 * @file feedforward_input_shaping.c
 * @brief Input shaping implementations: ZV, ZVD, EI, two-mode shapers.
 *
 * L5 --- Computational methods for vibration-free motion control.
 * Implements the Singer-Seering input shaping framework (1990).
 *
 * Key references:
 *   Smith, O.J.M. "Posicast Control of Damped Oscillatory Systems" (1957)
 *   Singer, N.C. & Seering, W.P. "Preshaping Command Inputs to Reduce
 *     System Vibration" (1990), ASME J. Dynamic Systems.
 *   Singhose, W.E. et al. "Extra-Insensitive Input Shapers" (1994)
 */

#include "feedforward_input_shaping.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L5: ZV Shaper Design (2 impulses)
 * ========================================================================== */

void shaper_design_zv(double wn, double zeta, InputShaper *shaper)
{
    if (!shaper || wn <= 0.0) return;

    double wd = wn * sqrt(1.0 - zeta * zeta);
    if (wd < 1e-9) wd = 1e-9;

    /* K = exp(-zeta*pi / sqrt(1-zeta^2)) */
    double K = exp(-zeta * M_PI / sqrt(1.0 - zeta * zeta));

    shaper->n_imp = 2;
    shaper->impulses = (Impulse *)malloc(2 * sizeof(Impulse));
    if (!shaper->impulses) {
        shaper->n_imp = 0;
        return;
    }

    /* A1 = 1/(1+K), t1 = 0 */
    shaper->impulses[0].amplitude = 1.0 / (1.0 + K);
    shaper->impulses[0].time = 0.0;

    /* A2 = K/(1+K), t2 = pi/wd (half period of damped natural frequency) */
    shaper->impulses[1].amplitude = K / (1.0 + K);
    shaper->impulses[1].time = M_PI / wd;

    shaper->total_time = shaper->impulses[1].time;
    shaper->wn = wn;
    shaper->zeta = zeta;
}

/* ==========================================================================
 * L5: ZVD Shaper Design (3 impulses)
 * ========================================================================== */

void shaper_design_zvd(double wn, double zeta, InputShaper *shaper)
{
    if (!shaper || wn <= 0.0) return;

    double wd = wn * sqrt(1.0 - zeta * zeta);
    if (wd < 1e-9) wd = 1e-9;

    double K = exp(-zeta * M_PI / sqrt(1.0 - zeta * zeta));
    double K2 = K * K;

    double denom = 1.0 + 2.0 * K + K2;

    shaper->n_imp = 3;
    shaper->impulses = (Impulse *)malloc(3 * sizeof(Impulse));
    if (!shaper->impulses) { shaper->n_imp = 0; return; }

    shaper->impulses[0].amplitude = 1.0 / denom;
    shaper->impulses[0].time = 0.0;

    shaper->impulses[1].amplitude = 2.0 * K / denom;
    shaper->impulses[1].time = M_PI / wd;

    shaper->impulses[2].amplitude = K2 / denom;
    shaper->impulses[2].time = 2.0 * M_PI / wd;

    shaper->total_time = shaper->impulses[2].time;
    shaper->wn = wn;
    shaper->zeta = zeta;
}

/* ==========================================================================
 * L5: EI Shaper Design (Extra-Insensitive, 3 impulses)
 * ========================================================================== */

void shaper_design_ei(double wn, double zeta, double v_tol,
                      InputShaper *shaper)
{
    if (!shaper || wn <= 0.0) return;

    /* EI shaper places the zero-vibration condition at two frequencies
     * symmetric about the design frequency, widening the insensitive band.
     *
     * The amplitudes and times are computed to keep vibration below V_tol
     * over the widest possible frequency range.
     *
     * Approximate formulas from Singhose (1994):
     *   For V_tol = 0.05 (5%): A1=0.25, A2=0.5, A3=0.25
     *     t1=0, t2=0.5*Td, t3=Td
     *   For general V_tol: compute via optimization (simplified here).
     */

    if (v_tol <= 0.0) v_tol = 0.05;
    if (v_tol > 1.0) v_tol = 1.0;

    double wd = wn * sqrt(1.0 - zeta * zeta);
    if (wd < 1e-9) wd = 1e-9;
    double Td = 2.0 * M_PI / wd; /* damped period */

    /* Simplified EI: use ZVD as base, adjust for V_tol */
    double K = exp(-zeta * M_PI / sqrt(1.0 - zeta * zeta));

    /* For EI, the time spacing is fractionally different from ZVD.
     * t2 ? 0.5*Td (vs 0.5*Td for ZVD which also has t2 = Td/2)
     * EI uses optimized times that create an insensitive "notch". */
    double scale = 1.0 - 0.2 * v_tol; /* heuristic scaling */

    shaper->n_imp = 3;
    shaper->impulses = (Impulse *)malloc(3 * sizeof(Impulse));
    if (!shaper->impulses) { shaper->n_imp = 0; return; }

    /* EI amplitudes */
    double A1 = 1.0 / (1.0 + 2.0 * K + K * K);
    double A2 = 2.0 * K / (1.0 + 2.0 * K + K * K);
    double A3 = K * K / (1.0 + 2.0 * K + K * K);

    /* Adjust for wider insensitivity: spread amplitudes */
    double spread = 0.1 * (1.0 - v_tol);
    A1 = 0.25 + spread * 0.25;
    A2 = 0.50 - spread * 0.50;
    A3 = 0.25 + spread * 0.25;

    shaper->impulses[0].amplitude = A1;
    shaper->impulses[0].time = 0.0;

    shaper->impulses[1].amplitude = A2;
    shaper->impulses[1].time = 0.5 * Td * scale;

    shaper->impulses[2].amplitude = A3;
    shaper->impulses[2].time = Td * scale;

    shaper->total_time = shaper->impulses[2].time;
    shaper->wn = wn;
    shaper->zeta = zeta;
}

/* ==========================================================================
 * L5: Two-Mode Shaper
 * ========================================================================== */

void shaper_design_two_mode(double wn1, double zeta1,
                            double wn2, double zeta2,
                            ShaperType type, InputShaper *shaper)
{
    if (!shaper) return;

    /* Design individual shapers, then convolve */
    InputShaper s1, s2;
    memset(&s1, 0, sizeof(s1));
    memset(&s2, 0, sizeof(s2));

    if (type == SHAPER_ZV) {
        shaper_design_zv(wn1, zeta1, &s1);
        shaper_design_zv(wn2, zeta2, &s2);
    } else {
        shaper_design_zvd(wn1, zeta1, &s1);
        shaper_design_zvd(wn2, zeta2, &s2);
    }

    /* Convolve: each impulse in s1 with each impulse in s2 */
    int n = s1.n_imp * s2.n_imp;
    shaper->n_imp = n;
    shaper->impulses = (Impulse *)calloc(n, sizeof(Impulse));
    if (!shaper->impulses) {
        shaper->n_imp = 0;
        shaper_free(&s1);
        shaper_free(&s2);
        return;
    }

    int idx = 0;
    for (int i = 0; i < s1.n_imp; i++) {
        for (int j = 0; j < s2.n_imp; j++) {
            shaper->impulses[idx].amplitude =
                s1.impulses[i].amplitude * s2.impulses[j].amplitude;
            shaper->impulses[idx].time =
                s1.impulses[i].time + s2.impulses[j].time;
            idx++;
        }
    }

    /* Sort by time (simple bubble sort, n <= 9 typically) */
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (shaper->impulses[j].time < shaper->impulses[i].time) {
                Impulse tmp = shaper->impulses[i];
                shaper->impulses[i] = shaper->impulses[j];
                shaper->impulses[j] = tmp;
            }
        }
    }

    /* Merge coincident impulses */
    int merged = 0;
    for (int i = 0; i < n - 1; i++) {
        if (fabs(shaper->impulses[i].time - shaper->impulses[i+1].time) < 1e-10) {
            shaper->impulses[i].amplitude += shaper->impulses[i+1].amplitude;
            /* Shift remaining */
            for (int j = i + 1; j < n - 1; j++) {
                shaper->impulses[j] = shaper->impulses[j + 1];
            }
            n--;
            i--;
            merged++;
        }
    }
    shaper->n_imp = n;

    shaper->total_time = shaper->impulses[n - 1].time;
    /* Average wn for display */
    shaper->wn = (wn1 + wn2) / 2.0;
    shaper->zeta = (zeta1 + zeta2) / 2.0;

    shaper_free(&s1);
    shaper_free(&s2);
}

/* ==========================================================================
 * L5: Negative ZV Shaper
 * ========================================================================== */

void shaper_design_negative_zv(double wn, double zeta, InputShaper *shaper)
{
    if (!shaper || wn <= 0.0) return;

    /* Negative shaper uses negative impulse(s) to reduce duration.
     * For ZV: 3 impulses (+, -, +) with shorter total time.
     * Time is reduced from Td/2 to Td/3 approximately. */

    double wd = wn * sqrt(1.0 - zeta * zeta);
    if (wd < 1e-9) wd = 1e-9;
    double Td = 2.0 * M_PI / wd;

    double K = exp(-zeta * M_PI / sqrt(1.0 - zeta * zeta));

    shaper->n_imp = 3;
    shaper->impulses = (Impulse *)malloc(3 * sizeof(Impulse));
    if (!shaper->impulses) { shaper->n_imp = 0; return; }

    /* Negative ZV shaper
     * A1 = 1/(1+K-K^2), t1 = 0
     * A2 = -K/(1+K-K^2), t2 = Td/3
     * A3 = K^2/(1+K-K^2), t3 = 2*Td/3 */
    double K2 = K * K;
    double denom = 1.0 + K - K2;
    if (fabs(denom) < 1e-15) denom = 1e-15;

    shaper->impulses[0].amplitude = 1.0 / denom;
    shaper->impulses[0].time = 0.0;

    shaper->impulses[1].amplitude = -K / denom;
    shaper->impulses[1].time = Td / 3.0;

    shaper->impulses[2].amplitude = K2 / denom;
    shaper->impulses[2].time = 2.0 * Td / 3.0;

    shaper->total_time = shaper->impulses[2].time;
    shaper->wn = wn;
    shaper->zeta = zeta;
}

/* ==========================================================================
 * L5: Shaper Application and Analysis
 * ========================================================================== */

void shaper_apply(const InputShaper *shaper,
                  const double *ref_in, int n_in, double dt,
                  double *ref_out, int *n_out)
{
    if (!shaper || !ref_in || !ref_out || !n_out) return;

    /* Convolution: y[k] = sum_j A_j * r(t_k - t_j)
     * This is a discrete-time convolution of the impulse sequence
     * with the input reference trajectory. */
    int extra = (int)ceil(shaper->total_time / dt) + 1;
    *n_out = n_in + extra;

    /* Initialize output to zero */
    for (int k = 0; k < *n_out; k++) {
        ref_out[k] = 0.0;
    }

    /* For each shaper impulse, add weighted delayed reference */
    for (int j = 0; j < shaper->n_imp; j++) {
        double Aj = shaper->impulses[j].amplitude;
        double delay = shaper->impulses[j].time;
        int delay_samples = (int)(delay / dt + 0.5); /* round to nearest */

        for (int k = 0; k < n_in; k++) {
            int out_idx = k + delay_samples;
            if (out_idx < *n_out) {
                ref_out[out_idx] += Aj * ref_in[k];
            }
        }
    }
}

double shaper_residual_vibration(const InputShaper *shaper,
                                 double r, double zeta)
{
    if (!shaper || shaper->n_imp == 0) return 1.0;

    /* Frequency ratio r = w / wn
     * wd_actual = r*wn * sqrt(1-zeta^2) = r*wd_design
     * V = exp(-zeta*wn*t_last) * sqrt(C^2 + S^2) */

    double w = r * shaper->wn;
    double wd = w * sqrt(1.0 - zeta * zeta);
    double C = 0.0, S = 0.0;

    for (int i = 0; i < shaper->n_imp; i++) {
        double Ai = shaper->impulses[i].amplitude;
        double ti = shaper->impulses[i].time;
        double decay = exp(zeta * shaper->wn * ti);
        C += Ai * decay * cos(wd * ti);
        S += Ai * decay * sin(wd * ti);
    }

    double t_last = shaper->impulses[shaper->n_imp - 1].time;
    double envelope = exp(-zeta * shaper->wn * t_last);

    return envelope * sqrt(C * C + S * S);
}

void shaper_sensitivity_curve(const InputShaper *shaper,
                              double zeta, double r_min, double r_max,
                              int npts, double *r_out, double *v_out)
{
    if (!shaper || !r_out || !v_out || npts <= 0) return;

    double dr = (r_max - r_min) / (double)(npts - 1);
    for (int i = 0; i < npts; i++) {
        double r = r_min + dr * i;
        r_out[i] = r;
        v_out[i] = shaper_residual_vibration(shaper, r, zeta);
    }
}

double shaper_insensitivity_range(const InputShaper *shaper,
                                  double zeta, double v_tol,
                                  double *w_low, double *w_high)
{
    if (!shaper || shaper->n_imp == 0) {
        if (w_low) *w_low = 1.0;
        if (w_high) *w_high = 1.0;
        return 0.0;
    }

    /* Search for frequency range where V < v_tol.
     * Scan from r=0.1 to r=2.0 in fine steps. */
    double low = 0.0, high = 0.0;
    int found_low = 0;

    for (int i = 0; i < 10000; i++) {
        double r = 0.1 + 1.9 * i / 9999.0;
        double V = shaper_residual_vibration(shaper, r, zeta);

        if (V <= v_tol && !found_low) {
            low = r;
            found_low = 1;
        }
        if (V <= v_tol && found_low) {
            high = r;
        }
    }

    if (w_low) *w_low = low;
    if (w_high) *w_high = high;
    return high - low;
}

void shaper_optimal_min_time(double wn, double zeta,
                             int max_amps, double v_tol,
                             double amp_limit, InputShaper *shaper)
{
    if (!shaper || max_amps < 2) return;
    (void)amp_limit; /* simplified: ignore amplitude constraint */

    /* For minimum time with V_tol constraint, start with ZV and add
     * impulses until vibration tolerance is met across desired range.
     * This is a simplified version of the optimization. */

    /* Start with ZVD as it has wider insensitivity */
    if (max_amps >= 3) {
        shaper_design_zvd(wn, zeta, shaper);
    } else {
        shaper_design_zv(wn, zeta, shaper);
    }

    /* Check if it meets tolerance at design frequency */
    double V = shaper_residual_vibration(shaper, 1.0, zeta);
    (void)V; /* ZV/ZVD should have V=0 at r=1 */
    (void)v_tol;

    /* For wider insensitivity, would iteratively add more impulses.
     * This simplified version uses ZVD as the practical optimum. */
    if (shaper->n_imp < max_amps && max_amps >= 4) {
        /* Could add more impulses here for wider insensitivity */
    }

    shaper->wn = wn;
    shaper->zeta = zeta;
}

/* ==========================================================================
 * Utility
 * ========================================================================== */

void shaper_free(InputShaper *shaper)
{
    if (shaper && shaper->impulses) {
        free(shaper->impulses);
        shaper->impulses = NULL;
        shaper->n_imp = 0;
    }
}

void shaper_to_string(const InputShaper *shaper, char *buf, int bufsz)
{
    if (!shaper || !buf || bufsz <= 0) return;
    int pos = 0;
    pos += snprintf(buf + pos, bufsz - pos,
                    "Shaper: %d impulses, total_time=%.4f s, wn=%.2f, zeta=%.3f\n",
                    shaper->n_imp, shaper->total_time, shaper->wn, shaper->zeta);
    for (int i = 0; i < shaper->n_imp && pos < bufsz; i++) {
        pos += snprintf(buf + pos, bufsz - pos,
                        "  Impulse %d: A=%.6f, t=%.6f\n",
                        i, shaper->impulses[i].amplitude,
                        shaper->impulses[i].time);
    }
}
