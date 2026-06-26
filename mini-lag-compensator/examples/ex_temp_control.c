/**
 * @file ex_temp_control.c
 * @brief Temperature control with lag compensator
 * L6: Canonical system -- thermal process (FOPDT) temperature regulation.
 */

#include "lag_compensator.h"
#include "lag_design.h"
#include <stdio.h>
#include <math.h>

extern int lag_design_temperature_control(double K_plant, double tau_plant,
    double theta_plant, double temp_accuracy_req, LagCompensator *comp);

int main(void) {
    printf("Temperature Control with Lag Compensator\n");
    printf("========================================\n\n");

    double K_reactor = 2.0;
    double tau_reactor = 120.0;
    double theta_reactor = 10.0;

    printf("Chemical Reactor Model:\n");
    printf("  G(s) = %.1f * exp(-%.0f*s) / (%.0f*s + 1)\n",
           K_reactor, theta_reactor, tau_reactor);
    printf("  (First-Order Plus Dead Time)\n\n");

    double ess_old = 1.0 / (1.0 + K_reactor);
    printf("Uncompensated steady-state error: %.2f%%\n\n", ess_old * 100.0);

    LagCompensator comp;
    int ret = lag_design_temperature_control(K_reactor, tau_reactor,
                                              theta_reactor, 0.005, &comp);
    if (ret != 0) { printf("Design returned code %d\n", ret); }

    printf("Designed Lag Compensator:\n");
    char buf[256];
    lag_to_string(&comp, buf, sizeof(buf));
    printf("  %s\n", buf);
    printf("  Kc=%.2f, beta=%.2f, T=%.2f s\n",
           lag_get_dc_gain(&comp), lag_get_beta(&comp),
           lag_get_time_constant(&comp));

    double w_low, w_high;
    lag_get_corner_frequencies(&comp, &w_low, &w_high);
    double w_theta = 1.0 / theta_reactor;
    printf("\n  Corner: w_low=%.4f, w_high=%.4f rad/s\n", w_low, w_high);
    printf("  Dead-time frequency: %.4f rad/s\n", w_theta);
    if (w_high < w_theta) {
        printf("  [OK] Corners below dead-time frequency\n");
    } else {
        printf("  [WARNING] Corners may interact with dead time\n");
    }

    double ess_new = 1.0 / (1.0 + K_reactor * lag_get_dc_gain(&comp));
    printf("\nCompensated ESS: %.4f%% (%.1fx improvement)\n\n",
           ess_new * 100.0, ess_old / ess_new);

    printf("Phase at key frequencies:\n");
    double freqs[] = {0.001, 0.01, 0.1};
    for (int i = 0; i < 3; i++) {
        double phi = lag_eval_phase(&comp, freqs[i]);
        printf("  w=%.3f: phi=%.2f deg\n", freqs[i], phi*180.0/M_PI);
    }

    printf("\nTemperature Control Design Complete\n");
    return 0;
}