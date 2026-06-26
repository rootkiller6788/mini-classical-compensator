
/**
 * example_level_control.c - PID Level Control of a Tank (L6/L7 Application)
 *
 * Tank level control:
 *   G(s) = K / s   (integrating process)
 *
 * Level is the integral of (inflow - outflow). Controlling level
 * with PID requires special care because the integrating nature
 * means P-only control yields zero steady-state error for step
 * disturbances, but PI is needed for ramp disturbances.
 *
 * For integrating processes, the Tyreus-Luyben tuning rules
 * are recommended over standard ZN.
 */
#include <stdio.h>
#include <math.h>
#include "../include/mini-pid-theory.h"
#include "../include/pid_tuning.h"

int main(void)
{
    /* Tank parameters */
    double A_tank = 2.0;      /* Cross-sectional area [m^2] */
    double K_valve = 0.05;    /* Valve flow coefficient [m^3/s per %] */
    double K_level = 1.0 / A_tank;  /* Level gain [m/m^3] */

    /* Effective integrator gain: K/s where K = K_valve*K_level */
    double K_int = K_valve * K_level;

    printf("=== PID Level Control of Tank ===\n");
    printf("Tank area: %.1f m^2, Valve gain: %.3f m^3/s/%%\n", A_tank, K_valve);
    printf("Integrator gain K = %.4f m/s/%%\n\n", K_int);

    /* For integrating processes, approximate as FOPDT with large tau */
    /* Using tau_approx = 10x the desired closed-loop time constant */
    double tau_approx = 100.0;
    double theta_approx = 1.0;  /* Small measurement delay */

    /* Apply Tyreus-Luyben (more conservative than ZN) */
    /* First need Ku, Tu. For integrator: Ku ~ pi/(2*K*theta), Tu ~ 4*theta */
    double Ku_est = 3.1416 / (2.0 * K_int * theta_approx);
    double Tu_est = 4.0 * theta_approx;

    pid_params_t params_tl, params_zn;
    pid_tune_tyreus_luyben(Ku_est, Tu_est, 2, &params_tl);
    pid_tune_zn_ultimate(Ku_est, Tu_est, 2, &params_zn);

    printf("Tyreus-Luyben Tuning (Conservative):\n");
    printf("  Kc = %.4f, Ti = %.1f s, Td = %.1f s\n", params_tl.Kc, params_tl.Ti, params_tl.Td);

    printf("Ziegler-Nichols Tuning (Aggressive):\n");
    printf("  Kc = %.4f, Ti = %.1f s, Td = %.1f s\n\n", params_zn.Kc, params_zn.Ti, params_zn.Td);

    /* Simulate with FOPDT approximation */
    pid_performance_t perf_tl, perf_zn;
    pid_evaluate_fopdt(K_int * tau_approx, tau_approx, theta_approx,
                       &params_tl, 2.0, 200.0, 500, &perf_tl);
    pid_evaluate_fopdt(K_int * tau_approx, tau_approx, theta_approx,
                       &params_zn, 2.0, 200.0, 500, &perf_zn);

    printf("Performance (2m level step):\n");
    printf("%-20s %10s %10s %10s %10s\n", "Method", "IAE", "OS[%]", "Ts[s]", "TV");
    printf("%-20s %10.3f %10.1f %10.1f %10.3f\n",
           "Tyreus-Luyben", perf_tl.iae, perf_tl.overshoot,
           perf_tl.settling_time, perf_tl.tv_u);
    printf("%-20s %10.3f %10.1f %10.1f %10.3f\n",
           "Ziegler-Nichols", perf_zn.iae, perf_zn.overshoot,
           perf_zn.settling_time, perf_zn.tv_u);

    /* Key insight: TL gives smoother control (less TV) at expense of speed */
    if (perf_tl.tv_u < perf_zn.tv_u) {
        printf("\nResult: Tyreus-Luyben provides %.1f%% smoother control action "
               "(important for valve wear)\n",
               (1.0 - perf_tl.tv_u/perf_zn.tv_u) * 100.0);
    }

    return 0;
}
