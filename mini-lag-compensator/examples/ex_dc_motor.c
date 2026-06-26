/**
 * @file ex_dc_motor.c
 * @brief Lag compensator design for DC motor speed control
 * L6: Canonical system -- DC motor speed regulation.
 */

#include "lag_compensator.h"
#include "lag_design.h"
#include "lag_frequency.h"
#include <stdio.h>
#include <math.h>

extern int lag_design_dc_motor_speed(const LagDCMotorParams *params,
    double speed_accuracy_req, LagCompensator *comp);
extern int lag_simulate_step_first_order(const LagCompensator *lag,
    double K_plant, double tau_plant, double t_end, int n_points,
    LagStepResponse *resp);
extern void lag_step_response_free(LagStepResponse *resp);

int main(void) {
    printf("DC Motor Speed Control with Lag Compensator\n");
    printf("========================================\n\n");

    LagDCMotorParams motor;
    motor.R = 2.0; motor.L = 0.005; motor.Kb = 0.05; motor.Kt = 0.05;
    motor.J = 0.0001; motor.B = 0.00001;
    motor.rated_speed = 500.0; motor.rated_torque = 0.1;
    motor.supply_voltage = 12.0;

    double denom = motor.B * motor.R + motor.Kt * motor.Kb;
    double Km = motor.Kt / denom;
    double tau_m = motor.J * motor.R / denom;

    printf("Motor: Km=%.4f, tau_m=%.4f s\n", Km, tau_m);
    double open_loop_dc = Km * motor.supply_voltage;
    double ess_uncomp = 1.0 / (1.0 + open_loop_dc);
    printf("Uncompensated ESS: %.2f%%\n\n", ess_uncomp * 100.0);

    LagCompensator comp;
    int ret = lag_design_dc_motor_speed(&motor, 0.02, &comp);
    if (ret != 0) { printf("Design failed!\n"); return 1; }

    char buf[256];
    lag_to_string(&comp, buf, sizeof(buf));
    printf("Designed: %s\n", buf);
    printf("  Kc=%.2f beta=%.2f T=%.4f\n",
           lag_get_dc_gain(&comp), lag_get_beta(&comp),
           lag_get_time_constant(&comp));

    double ess_comp = 1.0 / (1.0 + open_loop_dc * lag_get_dc_gain(&comp));
    printf("Compensated ESS: %.4f%% (%.1fx improvement)\n\n",
           ess_comp * 100.0, ess_uncomp / ess_comp);

    LagStepResponse resp;
    ret = lag_simulate_step_first_order(&comp, Km, tau_m, 0.2, 100, &resp);
    if (ret == 0) {
        printf("Step: y_ss=%.4f Ts=%.4f OS=%.1f%% ess=%.4f\n",
               resp.final_value, resp.settling_time,
               resp.overshoot_pct, resp.steady_state_error);
        lag_step_response_free(&resp);
    }

    printf("\nDC Motor Lag Compensation Complete\n");
    return 0;
}