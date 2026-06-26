
/**
 * example_flow_control.c - PID Flow Control (L6/L7 Application)
 *
 * Flow control loop:
 *   G(s) = K / (tau*s + 1)   (fast first-order process)
 *
 * Flow loops are typically very fast (tau < 1 second) and
 * noise-dominated. Derivative action amplifies measurement noise,
 * so PI control is standard. For noisy measurements, the derivative
 * filter coefficient N should be increased (N = 2-5).
 *
 * This example demonstrates:
 *   1. PI tuning for a fast flow loop
 *   2. Effect of derivative filter on noise amplification
 *   3. Setpoint filtering for smooth operation
 */
#include <stdio.h>
#include <math.h>
#include "../include/mini-pid-theory.h"
#include "../include/pid_tuning.h"
#include "../include/pid_advanced.h"

int main(void)
{
    /* Flow loop parameters (typical industrial flow control) */
    double K_flow = 0.8;     /* (L/s) per % valve */
    double tau_flow = 0.5;   /* seconds (fast response) */
    double theta_flow = 0.2; /* seconds (small transport delay) */

    printf("=== PID Flow Control ===\n");
    printf("Flow loop: K=%.1f (L/s)/%%, tau=%.2f s, theta=%.2f s\n\n",
           K_flow, tau_flow, theta_flow);

    /* PI tuning via IMC (lambda = 3*theta for noise-tolerant response) */
    pid_params_t params;
    double lambda = 3.0 * theta_flow;
    pid_tune_imc_lambda(K_flow, tau_flow, theta_flow, lambda, 1, &params);

    printf("IMC-PI Tuning (lambda=%.1f s):\n", lambda);
    printf("  Kc = %.3f, Ti = %.2f s\n", params.Kc, params.Ti);

    /* PI-only for flow control (no derivative due to noise) */
    params.Td = 0.0;
    params.deriv_mode = PID_DERIV_ON_MEASUREMENT;
    params.antiwindup = PID_WINDUP_BACK_CALC;

    /* Setpoint filter for smooth flow changes */
    pid_setpoint_filter_t sp_filter;
    pid_sp_filter_init(&sp_filter, 2.0, 0.5); /* 2s filter, 0.5 (L/s)/s ramp */

    /* Simulate */
    pid_state_t state;
    pid_init(&params, &state);
    double ysp = 10.0; /* Target flow: 10 L/s */
    double y = 0.0;    /* Initial flow */
    double u = 0.0;    /* Valve position (%) */

    printf("\nFlow control simulation (1 second sample period):\n");
    printf("Time[s]   SP[L/s]   PV[L/s]   Error   Valve[%%]\n");

    double Ts = 0.05; /* 50ms control period */
    for (int k = 0; k < 100; k++) {
        double t = k * Ts;
        /* Filter the setpoint */
        double ysp_filtered = pid_sp_filter_compute(&sp_filter, ysp, Ts);

        /* Flow measurement (simulated) */
        double alpha = exp(-Ts / tau_flow);
        y = alpha * y + K_flow * (1.0 - alpha) * u;

        /* PID control */
        u = pid_compute(&params, &state, ysp_filtered, y, Ts);

        if (k % 20 == 0) {
            printf("%6.2f    %7.2f   %7.2f   %6.2f   %7.2f\n",
                   t, ysp_filtered, y, ysp_filtered - y, u);
        }
    }

    printf("\nFinal flow: %.2f L/s (target: %.1f L/s, error: %.2f%%)\n",
           y, ysp, (ysp - y) / ysp * 100.0);
    printf("IAE accumulated: %.3f\n", state.integrated_absolute_error);

    return 0;
}
