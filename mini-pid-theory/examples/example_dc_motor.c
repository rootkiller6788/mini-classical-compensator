
/**
 * example_dc_motor.c - PID Control of a DC Motor (L6 Canonical System)
 *
 * DC motor transfer function:
 *   G(s) = K / (s * (tau*s + 1))
 * where:
 *   K = speed constant
 *   tau = mechanical time constant
 *
 * The PID controls the motor speed to a desired setpoint.
 * Tuning uses the Ziegler-Nichols ultimate sensitivity method
 * after finding Ku and Tu via relay feedback.
 */
#include <stdio.h>
#include <math.h>
#include "../include/mini-pid-theory.h"
#include "../include/pid_tuning.h"

int main(void)
{
    /* DC motor parameters (typical small servo motor) */
    double K_motor = 10.0;     /* Speed constant [rad/s/V] */
    double tau_motor = 0.05;   /* Mechanical time constant [s] */

    printf("=== PID Control of DC Motor ===\n");
    printf("Motor: K=%.1f rad/s/V, tau=%.2f s\n\n", K_motor, tau_motor);

    /* Step 1: Find ultimate gain and period via relay feedback */
    double Ku, Tu;
    pid_relay_feedback(K_motor, tau_motor, 0.0, 1.0, &Ku, &Tu);
    printf("Relay Auto-Tuning Results:\n");
    printf("  Ultimate Gain Ku = %.3f\n", Ku);
    printf("  Ultimate Period Tu = %.3f s\n\n", Tu);

    /* Step 2: Apply Ziegler-Nichols ultimate method tuning */
    pid_params_t params;
    pid_tune_zn_ultimate(Ku, Tu, 2, &params);
    printf("ZN Tuning Results:\n");
    printf("  Kc = %.3f, Ti = %.3f s, Td = %.3f s\n\n", params.Kc, params.Ti, params.Td);

    /* Step 3: Simulate step response */
    pid_performance_t perf;
    pid_evaluate_fopdt(K_motor, tau_motor, 0.0, &params, 1000.0, 0.5, 500, &perf);

    printf("Step Response Performance (1000 RPM step):\n");
    printf("  IAE = %.3f\n", perf.iae);
    printf("  ISE = %.3f\n", perf.ise);
    printf("  Overshoot = %.1f%%\n", perf.overshoot);
    printf("  Rise Time = %.3f s\n", perf.rise_time);
    printf("  Settling Time = %.3f s\n", perf.settling_time);
    printf("  Steady-State Error = %.3f RPM\n", perf.steady_state_error);
    printf("  Control Effort (TV) = %.3f\n", perf.tv_u);

    /* Step 4: Stability margins */
    double gm, pm, w180, w0db;
    pid_stability_margins(K_motor, tau_motor, 0.0, &params, &gm, &pm, &w180, &w0db);
    printf("\nStability Margins:\n");
    printf("  Gain Margin = %.1f dB\n", gm);
    printf("  Phase Margin = %.1f deg\n", pm);
    printf("  Bandwidth = %.3f rad/s\n", w0db);

    return 0;
}
