/**
 * @file ex_position_servo.c
 * @brief Position servo control with lag compensator
 * L6: Canonical system -- type-1 position servo tracking accuracy.
 */

#include "lag_compensator.h"
#include "lag_design.h"
#include "lag_frequency.h"
#include <stdio.h>
#include <math.h>

extern int lag_design_position_servo(double K_plant, double tau_plant,
    double position_error_req, LagCompensator *comp);
extern int lag_simulate_step_second_order(const LagCompensator *lag,
    double K_plant, double zeta, double wn, double t_end, int n_points,
    LagStepResponse *resp);
extern void lag_step_response_free(LagStepResponse *resp);

int main(void) {
    printf("Position Servo with Lag Compensator\n");
    printf("========================================\n\n");

    double K_servo = 10.0;
    double tau_servo = 0.05;
    double zeta_servo = 0.3;
    double wn_servo = 50.0;

    printf("Position Servo Model:\n");
    printf("  G(s) = %.1f / (s * (%.3f*s + 1))\n", K_servo, tau_servo);
    printf("  (Type-1 system with integrator)\n\n");

    double ess_uncomp = 1.0 / K_servo;
    printf("Uncompensated:\n");
    printf("  Kv = %.2f, ramp error = %.4f (%.1f%%)\n\n",
           K_servo, ess_uncomp, ess_uncomp * 100.0);

    LagCompensator comp;
    double accuracy_req = 0.002;
    int ret = lag_design_position_servo(K_servo, tau_servo,
                                         accuracy_req, &comp);
    if (ret != 0) { printf("Design returned code %d\n", ret); }

    printf("Designed Lag Compensator:\n");
    char buf[256];
    lag_to_string(&comp, buf, sizeof(buf));
    printf("  %s\n", buf);

    double Kv_comp = K_servo * lag_get_dc_gain(&comp);
    double ess_comp = 1.0 / Kv_comp;
    printf("\nCompensated:\n");
    printf("  Kv = %.2f (%.1fx), ramp error = %.6f\n",
           Kv_comp, lag_get_dc_gain(&comp), ess_comp);

    double w_gc = K_servo;
    double phi_at_gc = lag_eval_phase(&comp, w_gc);
    printf("\nPhase at gain crossover (w=%.1f):\n", w_gc);
    printf("  phi = %.2f deg, |G_c| = %.4f\n",
           phi_at_gc * 180.0 / M_PI, lag_eval_magnitude(&comp, w_gc));
    if (fabs(phi_at_gc) < 0.1) {
        printf("  [OK] Negligible phase impact\n");
    } else {
        printf("  [NOTE] %.1f deg phase reduction\n",
               fabs(phi_at_gc)*180.0/M_PI);
    }

    printf("\nFrequency sweep (compensator):\n");
    printf("  w (rad/s)    |G_c| (dB)    phase (deg)\n");
    double freqs[] = {0.1, 1.0, 10.0, 100.0};
    for (int i = 0; i < 4; i++) {
        LagFreqPoint fp = lag_eval_frequency(&comp, freqs[i]);
        printf("  %8.2f    %8.2f    %8.2f\n",
               fp.omega, fp.magnitude_db, fp.phase_deg);
    }

    LagStepResponse resp;
    ret = lag_simulate_step_second_order(&comp, K_servo, zeta_servo,
                                          wn_servo, 0.5, 200, &resp);
    if (ret == 0) {
        printf("\nStep Response (compensated):\n");
        printf("  y_ss=%.4f OS=%.1f%% Ts=%.4f ess=%.6f\n",
               resp.final_value, resp.overshoot_pct,
               resp.settling_time, resp.steady_state_error);
        lag_step_response_free(&resp);
    }

    printf("\nPosition Servo Design Complete\n");
    return 0;
}