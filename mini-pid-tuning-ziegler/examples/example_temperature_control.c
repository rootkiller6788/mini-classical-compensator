/**
 * example_temperature_control.c — Temperature Loop PID Tuning
 *
 * Demonstrates end-to-end PID tuning for a temperature control loop
 * using multiple methods and comparing results.
 *
 * Knowledge: L6 (temperature canonical system), L7 (process control).
 */

#include "pid_tuning.h"
#include "ziegler_nichols.h"
#include "cohen_coon.h"
#include "imc_tuning.h"
#include "gain_margin_tuning.h"
#include "advanced_tuning.h"
#include "application_tuning.h"
#include "fopdt_model.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("Temperature Control Loop PID Tuning Demo\n");
    printf("========================================\n\n");

    /* A typical thermal process: */
    /*   G(s) = 1.8 * exp(-2.0*s) / (25.0*s + 1) */
    /*   K = 1.8 °C/%, T = 25 min, L = 2 min */
    fopdt_model_t furnace;
    furnace.K = 1.8;
    furnace.T = 25.0;
    furnace.L = 2.0;

    printf("Process: FOPDT G(s) = %.2f*exp(-%.1f*s)/(%.1f*s+1)\n",
           furnace.K, furnace.L, furnace.T);
    printf("Normalized dead time ratio: τ = %.3f\n\n",
           fopdt_dead_time_ratio(&furnace));

    /* Method 1: Ziegler-Nichols Step Response */
    printf("--- Method 1: Ziegler-Nichols Step Response (1942) ---\n");
    pid_params_t zn_params;
    zn_step_tune(&furnace, ZN_CONTROLLER_PID, &zn_params);
    printf("  Kp = %.3f, Ti = %.3f, Td = %.3f\n",
           zn_params.Kp, zn_params.Ti, zn_params.Td);

    double overshoot_zn, decay_zn, omega_zn;
    cohen_coon_predict_response(&furnace, &zn_params,
                                 &decay_zn, &overshoot_zn, &omega_zn);
    printf("  Predicted overshoot: %.1f%%, decay ratio: %.3f\n\n",
           overshoot_zn, decay_zn);

    /* Method 2: Cohen-Coon */
    printf("--- Method 2: Cohen-Coon (1953) ---\n");
    pid_params_t cc_params;
    cohen_coon_tune(&furnace, CC_CONTROLLER_PID, &cc_params);
    printf("  Kp = %.3f, Ti = %.3f, Td = %.3f\n",
           cc_params.Kp, cc_params.Ti, cc_params.Td);

    double overshoot_cc, decay_cc, omega_cc;
    cohen_coon_predict_response(&furnace, &cc_params,
                                 &decay_cc, &overshoot_cc, &omega_cc);
    printf("  Predicted overshoot: %.1f%%, decay ratio: %.3f\n\n",
           overshoot_cc, decay_cc);

    /* Method 3: SIMC (Skogestad IMC) */
    printf("--- Method 3: SIMC/IMC (Skogestad 2003) ---\n");
    pid_params_t simc_params;
    imc_simc_tune(&furnace, furnace.L * 1.5, 1, &simc_params);
    printf("  Kp = %.3f, Ti = %.3f, Td = %.3f\n",
           simc_params.Kp, simc_params.Ti, simc_params.Td);

    /* Check margins for SIMC */
    margin_spec_t simc_margins;
    gm_pm_compute_margins(&furnace, &simc_params, &simc_margins);
    printf("  Gain margin: %.1f (%.1f dB), Phase margin: %.1f°\n\n",
           simc_margins.Am, 20.0*log10(simc_margins.Am),
           simc_margins.Pm * 180.0 / M_PI);

    /* Method 4: Gain/Phase Margin specification */
    printf("--- Method 4: Gm=3, Pm=50° design ---\n");
    pid_params_t gm_params;
    gm_pm_tune_pid(&furnace, 3.0, 50.0 * M_PI / 180.0, 0, &gm_params);
    printf("  Kp = %.3f, Ti = %.3f, Td = %.3f\n",
           gm_params.Kp, gm_params.Ti, gm_params.Td);

    /* Compare methods */
    printf("\n--- Method Comparison (estimated ITAE ranking) ---\n");
    pid_tune_method_t methods[] = {
        TUNE_METHOD_ZN_STEP,
        TUNE_METHOD_COHEN_COON,
        TUNE_METHOD_IMC,
        TUNE_METHOD_GM_PM
    };
    const char *names[] = { "Z-N Step", "Cohen-Coon", "IMC/SIMC", "Gm/Pm" };
    double rankings[4];
    int best;

    pid_compare_tuning_methods(&furnace, methods, 4, rankings, &best);
    for (int i = 0; i < 4; i++) {
        printf("  %-15s score: %.3f %s\n",
               names[i], rankings[i],
               (i == best) ? "★ BEST" : "");
    }

    /* Application-specific recommendation */
    printf("\n--- Application-Specific Tuning ---\n");
    pid_params_t app_params;
    app_tune_temperature(furnace.K, furnace.T, furnace.L, &app_params);
    printf("  Temperature-optimized: Kp=%.3f, Ti=%.3f, Td=%.3f\n",
           app_params.Kp, app_params.Ti, app_params.Td);

    printf("\nDemo complete.\n");
    return 0;
}
