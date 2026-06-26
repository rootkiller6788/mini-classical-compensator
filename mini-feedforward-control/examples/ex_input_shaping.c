/**
 * @file ex_input_shaping.c
 * @brief Example: Input shaping for crane anti-sway control.
 *
 * Demonstrates L5-L6: ZV, ZVD, and EI shapers applied to a gantry crane
 * to eliminate payload swing during point-to-point motion.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "feedforward_input_shaping.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Input Shaping for Crane Anti-Sway Control ===\n\n");

    /* Overhead crane parameters:
     * Cable length L = 3.0 m
     * Payload mass m = 500 kg
     * Natural frequency: wn = sqrt(g/L) = sqrt(9.81/3.0) = 1.808 rad/s
     * Period: T = 2*pi/wn = 3.474 s
     * Damping: very low, zeta ? 0.005 */
    double g = 9.81;
    double L = 3.0;
    double wn = sqrt(g / L);
    double zeta = 0.005;
    double T = 2.0 * M_PI / wn;

    printf("Crane parameters:\n");
    printf("  Cable length L = %.2f m\n", L);
    printf("  Natural frequency wn = %.3f rad/s\n", wn);
    printf("  Period T = %.3f s\n", T);
    printf("  Damping zeta = %.4f\n\n", zeta);

    /* Design three shapers */
    InputShaper shaper_zv, shaper_zvd, shaper_ei;
    memset(&shaper_zv, 0, sizeof(shaper_zv));
    memset(&shaper_zvd, 0, sizeof(shaper_zvd));
    memset(&shaper_ei, 0, sizeof(shaper_ei));

    shaper_design_zv(wn, zeta, &shaper_zv);
    shaper_design_zvd(wn, zeta, &shaper_zvd);
    shaper_design_ei(wn, zeta, 0.05, &shaper_ei);

    /* Print shaper details */
    char buf[512];

    printf("--- ZV Shaper (2 impulses) ---\n");
    shaper_to_string(&shaper_zv, buf, sizeof(buf));
    printf("%s\n", buf);

    printf("--- ZVD Shaper (3 impulses) ---\n");
    shaper_to_string(&shaper_zvd, buf, sizeof(buf));
    printf("%s\n", buf);

    printf("--- EI Shaper (3 impulses, 5%% tolerance) ---\n");
    shaper_to_string(&shaper_ei, buf, sizeof(buf));
    printf("%s\n", buf);

    /* Compare shaper durations */
    printf("Shader duration comparison:\n");
    printf("  ZV:  %.4f s (%.1f%% of period)\n",
           shaper_zv.total_time, 100.0 * shaper_zv.total_time / T);
    printf("  ZVD: %.4f s (%.1f%% of period)\n",
           shaper_zvd.total_time, 100.0 * shaper_zvd.total_time / T);
    printf("  EI:  %.4f s (%.1f%% of period)\n",
           shaper_ei.total_time, 100.0 * shaper_ei.total_time / T);

    /* Sensitivity analysis */
    printf("\n--- Sensitivity Analysis (Vibration vs Frequency Ratio) ---\n");
    printf("  r       V_zv(r)     V_zvd(r)    V_ei(r)\n");
    printf("  ------  ----------  ----------  ----------\n");

    int n_pts = 21;
    for (int i = 0; i < n_pts; i++) {
        double r = 0.5 + 1.5 * i / (n_pts - 1.0);
        double V_zv = shaper_residual_vibration(&shaper_zv, r, zeta);
        double V_zvd = shaper_residual_vibration(&shaper_zvd, r, zeta);
        double V_ei = shaper_residual_vibration(&shaper_ei, r, zeta);
        printf("  %.3f   %.6f    %.6f    %.6f\n", r, V_zv, V_zvd, V_ei);
    }

    /* Insensitivity range */
    double w_low, w_high;
    double insens_zv = shaper_insensitivity_range(&shaper_zv, zeta, 0.05,
                                                   &w_low, &w_high);
    double insens_zvd = shaper_insensitivity_range(&shaper_zvd, zeta, 0.05,
                                                    &w_low, &w_high);
    double insens_ei = shaper_insensitivity_range(&shaper_ei, zeta, 0.05,
                                                   &w_low, &w_high);

    printf("\n5%% Insensitivity range (frequency ratio):\n");
    printf("  ZV:  %.3f  (%.3f - %.3f)\n", insens_zv,
           (insens_zv > 0) ? 1.0 - insens_zv/2 : 1.0,
           (insens_zv > 0) ? 1.0 + insens_zv/2 : 1.0);
    printf("  ZVD: %.3f\n", insens_zvd);
    printf("  EI:  %.3f\n", insens_ei);

    /* Apply shapers to a step command (trolley move 2 meters) */
    printf("\n--- Application: 2m trolley move ---\n");

    double dt = 0.01;
    int n_in = 300;  /* 3 seconds */
    double *ref_in = (double *)malloc(n_in * sizeof(double));
    for (int i = 0; i < n_in; i++) {
        ref_in[i] = (i >= 10) ? 2.0 : 0.0;  /* Step at t=0.1s */
    }

    double *ref_zv = (double *)malloc((n_in + 200) * sizeof(double));
    double *ref_zvd = (double *)malloc((n_in + 200) * sizeof(double));
    int n_zv, n_zvd;

    shaper_apply(&shaper_zv, ref_in, n_in, dt, ref_zv, &n_zv);
    shaper_apply(&shaper_zvd, ref_in, n_in, dt, ref_zvd, &n_zvd);

    printf("  Unshaped:  %d samples, rise at t=0.1s\n", n_in);
    printf("  ZV shaped: %d samples, delay=%.2fs\n", n_zv, shaper_zv.total_time);
    printf("  ZVD shaped: %d samples, delay=%.2fs\n", n_zvd, shaper_zvd.total_time);

    /* Print first few shaped outputs */
    printf("\n  First 20 shaped outputs (ZVD):\n  ");
    for (int i = 0; i < 20 && i < n_zvd; i++) {
        printf("%.3f ", ref_zvd[i]);
    }
    printf("\n");

    free(ref_in);
    free(ref_zv);
    free(ref_zvd);
    shaper_free(&shaper_zv);
    shaper_free(&shaper_zvd);
    shaper_free(&shaper_ei);

    printf("\n=== Done ===\n");
    return 0;
}
