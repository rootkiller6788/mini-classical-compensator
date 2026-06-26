/**
 * @file example_lead_dc_motor.c
 * @brief Lead compensation for DC motor velocity control (L6)
 *
 * Plant: G(s) = 50 / ((s+5)*(s+10)) — DC motor transfer function
 * Design lead compensator to achieve PM >= 55 deg.
 *
 * Application: Manufacturing automation, conveyor belt speed control.
 * Purdue ECE 602, Georgia Tech AE 6530, Tsinghua Auto Control
 */

#include "lead_compensator.h"
#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <stdio.h>
#include <string.h>

int main(void) {
    printf("=== Lead Compensator: DC Motor Velocity Control ===\n\n");

    /* DC motor: G(s) = 50 / ((s+5)(s+10)) = 50 / (s^2 + 15s + 50) */
    lead_system_t motor;
    memset(&motor, 0, sizeof(motor));
    motor.tf.num[0] = 50.0;
    motor.tf.num_order = 0;
    motor.tf.den[0] = 50.0;
    motor.tf.den[1] = 15.0;
    motor.tf.den[2] = 1.0;
    motor.tf.den_order = 2;
    motor.tf.gain = 1.0;
    strcpy(motor.name, "DC Motor Vel");

    /* Uncompensated analysis */
    double pm_raw = lead_compute_phase_margin(&motor);
    double gm_raw = lead_compute_gain_margin(&motor);
    printf("Uncompensated: PM = %.1f deg, GM = %.1f dB\n\n", pm_raw, gm_raw);

    /* Design specs */
    lead_specs_t specs;
    memset(&specs, 0, sizeof(specs));
    specs.phase_margin_desired = 55.0;
    specs.gain_margin_desired = 10.0;
    specs.steady_state_error = 0.02;
    specs.max_overshoot_pct = 15.0;
    specs.settling_time = 1.0;
    specs.use_frequency_domain = true;
    specs.method = LEAD_METHOD_FREQUENCY;

    /* Design */
    lead_design_result_t result;
    lead_design_frequency(&motor, &specs, &result);

    /* Display compensator */
    lead_print(&result.compensator);

    /* Compensated metrics */
    double pm_new = lead_compensated_phase_margin(&result.compensator, &motor);
    double gm_new = lead_compensated_gain_margin(&result.compensator, &motor);
    double w_gc = lead_open_loop_crossover(&result.compensator, &motor);
    double bw = lead_closed_loop_bandwidth(&result.compensator, &motor);

    printf("\nCompensated: PM = %.1f deg, GM = %.1f dB\n", pm_new, gm_new);
    printf("Gain Crossover: %.3f rad/s\n", w_gc);
    printf("Closed-loop BW: %.3f rad/s\n", bw);

    /* Sensitivity at key frequencies */
    double S_lf = lead_sensitivity(&result.compensator, &motor, 0.1);
    double S_mf = lead_sensitivity(&result.compensator, &motor, w_gc);
    double S_hf = lead_sensitivity(&result.compensator, &motor, 100.0);
    printf("\nSensitivity |S(jw)|:\n");
    printf("  Low freq (0.1):  %.4f (%.1f dB)\n", S_lf, 20*log10(S_lf));
    printf("  At w_gc:         %.4f (%.1f dB)\n", S_mf, 20*log10(S_mf));
    printf("  High freq (100): %.4f (%.1f dB)\n", S_hf, 20*log10(S_hf));

    /* Noise analysis */
    double noise_hf = lead_noise_amplification(&result.compensator, &motor, 500.0);
    printf("\nNoise amplification at 500 rad/s: %.4f\n", noise_hf);

    /* Stability check */
    bool stable = lead_is_stable_nyquist(&result.compensator, &motor);
    printf("Nyquist stability: %s\n", stable ? "STABLE" : "UNSTABLE");

    printf("\n=== DC Motor Lead Compensation Complete ===\n");
    return 0;
}