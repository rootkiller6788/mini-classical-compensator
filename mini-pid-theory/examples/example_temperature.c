
/**
 * example_temperature.c - PID Temperature Control (L6/L7 Application)
 *
 * Temperature control of a heated tank:
 *   G(s) = K * exp(-theta*s) / (tau*s + 1)
 *
 * Typical values for a small industrial heater:
 *   K = 2.0 degC/%  (2 degrees per percent heater output)
 *   tau = 120 s     (2 minute thermal time constant)
 *   theta = 15 s    (15 second measurement/actuation delay)
 *
 * This is a lag-dominant process with moderate deadtime.
 * SIMC tuning is well-suited for this application.
 */
#include <stdio.h>
#include <math.h>
#include "../include/mini-pid-theory.h"
#include "../include/pid_tuning.h"

int main(void)
{
    /* Process parameters */
    double K_heater = 2.0;    /* degC per % output */
    double tau_heat = 120.0;  /* seconds */
    double theta_heat = 15.0; /* seconds */

    printf("=== PID Temperature Control of Heated Tank ===\n");
    printf("Process: K=%.1f degC/%%, tau=%.0f s, theta=%.0f s\n", K_heater, tau_heat, theta_heat);
    printf("Normalized deadtime = %.3f\n\n", pid_normalized_deadtime(tau_heat, theta_heat));

    /* Compare multiple tuning methods */
    pid_params_t params_zn, params_cc, params_simc;
    pid_performance_t perf_zn, perf_cc, perf_simc;

    /* Ziegler-Nichols Step Response */
    pid_tune_zn_step(K_heater, tau_heat, theta_heat, 2, &params_zn);
    pid_evaluate_fopdt(K_heater, tau_heat, theta_heat, &params_zn, 60.0, 600.0, 600, &perf_zn);

    /* Cohen-Coon */
    pid_tune_cohen_coon(K_heater, tau_heat, theta_heat, 2, &params_cc);
    pid_evaluate_fopdt(K_heater, tau_heat, theta_heat, &params_cc, 60.0, 600.0, 600, &perf_cc);

    /* SIMC (tau_c = theta for tight control) */
    pid_tune_simc(K_heater, tau_heat, theta_heat, theta_heat, 2, &params_simc);
    pid_evaluate_fopdt(K_heater, tau_heat, theta_heat, &params_simc, 60.0, 600.0, 600, &perf_simc);

    printf("Tuning Comparison (Setpoint: 60 degC, process initially at 20 degC):\n");
    printf("%-25s %10s %10s %10s %10s %10s\n", "Method", "Kc", "Ti[s]", "Td[s]", "IAE", "OS[%]");
    printf("%-25s %10.3f %10.1f %10.1f %10.1f %10.1f\n",
           "Ziegler-Nichols Step", params_zn.Kc, params_zn.Ti, params_zn.Td,
           perf_zn.iae, perf_zn.overshoot);
    printf("%-25s %10.3f %10.1f %10.1f %10.1f %10.1f\n",
           "Cohen-Coon", params_cc.Kc, params_cc.Ti, params_cc.Td,
           perf_cc.iae, perf_cc.overshoot);
    printf("%-25s %10.3f %10.1f %10.1f %10.1f %10.1f\n",
           "SIMC", params_simc.Kc, params_simc.Ti, params_simc.Td,
           perf_simc.iae, perf_simc.overshoot);

    /* Stability margins for SIMC tuning */
    double gm, pm, w180, w0db;
    pid_stability_margins(K_heater, tau_heat, theta_heat, &params_simc, &gm, &pm, &w180, &w0db);
    printf("\nSIMC Stability: GM=%.1f dB, PM=%.1f deg, BW=%.4f rad/s\n", gm, pm, w0db);

    return 0;
}
