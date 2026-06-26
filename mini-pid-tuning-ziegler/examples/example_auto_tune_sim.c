/**
 * example_auto_tune_sim.c — Relay Auto-Tuning Simulation
 *
 * Demonstrates the Åström-Hägglund relay auto-tuning procedure:
 *   relay experiment → Ku/Pu identification → PID tuning.
 *
 * Knowledge: L5 (relay auto-tuning), L6 (process identification),
 *             L7 (industrial auto-tuners).
 */

#include "pid_tuning.h"
#include "ziegler_nichols.h"
#include "relay_autotune.h"
#include "fopdt_model.h"
#include "gain_margin_tuning.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("Relay Auto-Tuning Simulation\n");
    printf("============================\n\n");

    /* Define three different processes */
    fopdt_model_t processes[] = {
        { 1.0, 10.0, 1.0 },  /* Lag-dominant: T >> L */
        { 1.0,  5.0, 2.0 },  /* Balanced */
        { 1.0,  2.0, 3.0 },  /* Delay-dominant: L > T/2 */
    };
    const char *desc[] = {
        "Lag-dominant (T >> L)",
        "Balanced (T ≈ 2.5*L)",
        "Delay-dominant (L > T/2)"
    };

    for (int p = 0; p < 3; p++) {
        printf("Process %d: %s\n", p + 1, desc[p]);
        printf("  G(s) = %.1f * exp(-%.1f*s) / (%.1f*s + 1)\n",
               processes[p].K, processes[p].L, processes[p].T);
        printf("  τ = %.3f\n\n", fopdt_dead_time_ratio(&processes[p]));

        /* Design relay experiment */
        relay_config_t config;
        relay_design_experiment(processes[p].K,
                                processes[p].T + processes[p].L,
                                20.0, &config);
        printf("  Relay: amplitude=%.3f, hysteresis=%.3f, Ts=%.3f\n",
               config.amplitude, config.hysteresis, config.Ts);

        /* Run relay experiment simulation */
        relay_result_t result;
        int ret = relay_simulate_fopdt(&processes[p], &config, &result);
        if (ret == 0 && result.converged) {
            printf("  Identified: Ku = %.3f, Pu = %.3f sec\n",
                   result.Ku, result.Pu);

            /* Apply Z-N frequency response PID tuning */
            pid_params_t zn_pid;
            zn_freq_tune(result.Ku, result.Pu, ZN_CONTROLLER_PID, &zn_pid);
            printf("  Z-N PID: Kp=%.3f, Ti=%.3f, Td=%.3f\n",
                   zn_pid.Kp, zn_pid.Ti, zn_pid.Td);

            /* Alternative: Custom conservative tuning */
            pid_params_t cons_pid;
            zn_freq_tune_custom(result.Ku, result.Pu, ZN_CONTROLLER_PID,
                                0.3, 0.5, 0.125, &cons_pid);
            printf("  Conservative PID: Kp=%.3f, Ti=%.3f, Td=%.3f\n",
                   cons_pid.Kp, cons_pid.Ti, cons_pid.Td);

            /* Compute achieved margins */
            margin_spec_t margins;
            gm_pm_compute_margins(&processes[p], &zn_pid, &margins);
            printf("  Z-N margins: Gm=%.1f (%.1f dB), Pm=%.1f°\n",
                   margins.Am, 20.0*log10(margins.Am > 0.01 ? margins.Am : 0.01),
                   margins.Pm * 180.0 / M_PI);
        } else {
            printf("  ⚠ Relay experiment did not converge\n");
        }
        printf("\n");
    }

    /* Demonstrate variable hysteresis for safe operation */
    printf("--- Variable Hysteresis Demonstration ---\n");
    fopdt_model_t sensitive = { 0.5, 15.0, 3.0 };
    relay_result_t vh_result;
    int ret = relay_variable_hysteresis(&sensitive, 2.0,
                                         0.05, 0.5, 2.0, 0.3, &vh_result);
    if (ret == 0) {
        printf("  Variable hysteresis: Ku=%.3f, Pu=%.3f\n",
               vh_result.Ku, vh_result.Pu);
    } else {
        printf("  Variable hysteresis: target amplitude not met\n");
    }

    /* Retune trigger test */
    printf("\n--- Retune Trigger Test ---\n");
    int needs_retune = relay_needs_retune(3.0, 5.0, 4.0, 5.2, 0.2, 0.1);
    printf("  Ku: 3.0→4.0 (33%%), Pu: 5.0→5.2 (4%%)\n");
    printf("  Retune recommended: %s\n", needs_retune ? "YES" : "NO");

    printf("\nRelay auto-tuning demo complete.\n");
    return 0;
}
