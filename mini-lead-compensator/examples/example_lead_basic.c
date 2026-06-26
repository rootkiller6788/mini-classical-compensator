/**
 * @file example_lead_basic.c
 * @brief Basic lead compensator design and analysis example
 *
 * L6 - Classic problem: Lead compensation for a DC motor position servo.
 * Plant: G(s) = 10 / (s*(s+2))
 * Design lead compensator for PM >= 45 deg, steady-state error <= 5%.
 *
 * MIT 6.302, Stanford ENGR105, Berkeley ME132
 */

#include "lead_compensator.h"
#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("=== Lead Compensator Design: DC Motor Position Servo ===\n\n");

    /* Define plant: G(s) = 10 / (s*(s+2)) = 10 / (s^2 + 2s) */
    lead_system_t plant;
    memset(&plant, 0, sizeof(plant));
    plant.tf.num[0] = 10.0;
    plant.tf.num_order = 0;
    plant.tf.den[0] = 0.0;   /* constant term */
    plant.tf.den[1] = 2.0;   /* s term */
    plant.tf.den[2] = 1.0;   /* s^2 term */
    plant.tf.den_order = 2;
    plant.tf.gain = 1.0;
    strcpy(plant.name, "DC Motor Servo");

    printf("Plant: G(s) = 10 / (s*(s+2))\n\n");

    /* Analyze uncompensated system */
    printf("--- Uncompensated System ---\n");
    double pm_uncomp = lead_compute_phase_margin(&plant);
    double gm_uncomp = lead_compute_gain_margin(&plant);
    double w_gc_uncomp = lead_find_gain_crossover(&plant,
                                                   LEAD_FREQ_MIN_DEFAULT,
                                                   LEAD_FREQ_MAX_DEFAULT);
    double w_pc_uncomp = lead_find_phase_crossover(&plant,
                                                    LEAD_FREQ_MIN_DEFAULT,
                                                    LEAD_FREQ_MAX_DEFAULT);

    printf("Phase Margin:      %.2f deg\n", pm_uncomp);
    printf("Gain Margin:       %.2f dB\n", gm_uncomp);
    printf("Gain Crossover:    %.4f rad/s\n", w_gc_uncomp);
    printf("Phase Crossover:   %.4f rad/s\n\n", w_pc_uncomp);

    /* Design specifications */
    lead_specs_t specs;
    memset(&specs, 0, sizeof(specs));
    specs.phase_margin_desired = 45.0;
    specs.gain_margin_desired = 12.0;
    specs.steady_state_error = 0.05;
    specs.max_overshoot_pct = 20.0;
    specs.settling_time = 2.0;
    specs.use_frequency_domain = true;
    specs.method = LEAD_METHOD_FREQUENCY;

    printf("Design Specifications:\n");
    printf("  Phase Margin >= %.0f deg\n", specs.phase_margin_desired);
    printf("  Gain Margin >= %.0f dB\n", specs.gain_margin_desired);
    printf("  Steady-state Error <= %.0f%%\n", specs.steady_state_error * 100.0);
    printf("  Max Overshoot <= %.0f%%\n", specs.max_overshoot_pct);
    printf("  Settling Time <= %.1f s\n\n", specs.settling_time);

    /* Design lead compensator */
    lead_design_result_t result;
    bool ok = lead_design_frequency(&plant, &specs, &result);

    if (!ok) {
        printf("Design did not fully converge.\n\n");
    }

    printf("--- Designed Lead Compensator ---\n");
    lead_print(&result.compensator);
    printf("\n");

    printf("--- Compensated System ---\n");
    double pm_comp = result.achieved_phase_margin;
    double gm_comp = result.achieved_gain_margin;
    double w_gc_comp = result.crossover_freq_new;

    printf("Phase Margin:      %.2f deg (was %.2f deg)\n", pm_comp, pm_uncomp);
    printf("Gain Margin:       %.2f dB (was %.2f dB)\n", gm_comp, gm_uncomp);
    printf("Gain Crossover:    %.4f rad/s (was %.4f rad/s)\n",
           w_gc_comp, w_gc_uncomp);
    printf("Phase Lead Added:  %.2f deg\n", result.phase_lead_added);
    printf("Alpha:             %.4f\n", result.alpha_actual);
    printf("T:                 %.4f s\n", result.T_actual);
    printf("Iterations:        %d\n", result.iterations);
    printf("Converged:         %s\n\n", result.converged ? "Yes" : "No");

    /* Predict performance */
    lead_performance_t perf;
    lead_predict_performance(pm_comp, w_gc_comp, &perf);
    printf("--- Predicted Closed-Loop Performance ---\n");
    printf("Damping Ratio:     %.3f\n", perf.damping_ratio);
    printf("Natural Frequency: %.3f rad/s\n", perf.natural_freq);
    printf("Percent Overshoot: %.1f%%\n", perf.percent_overshoot);
    printf("Settling Time:     %.3f s\n", perf.settling_time);
    printf("Rise Time:         %.3f s\n", perf.rise_time);
    printf("Peak Time:         %.3f s\n", perf.peak_time);
    printf("Bandwidth:         %.3f rad/s\n\n", perf.bandwidth);

    /* Robustness analysis */
    printf("--- Robustness Analysis ---\n");
    double M_s = lead_peak_sensitivity(&result.compensator, &plant,
                                       0.01, 100.0, 500);
    double M_t = lead_peak_complementary_sensitivity(&result.compensator, &plant,
                                                      0.01, 100.0, 500);
    double M_m = lead_modulus_margin(&result.compensator, &plant,
                                      0.01, 100.0, 500);
    double delay_m = lead_delay_margin(&result.compensator, &plant);

    printf("Peak Sensitivity (M_s):          %.3f\n", M_s);
    printf("Peak Complementary Sens (M_t):   %.3f\n", M_t);
    printf("Modulus Margin:                  %.4f\n", M_m);
    printf("Delay Margin:                    %.6f s\n", delay_m);
    printf("Excessive HF Gain:               %s\n",
           lead_excessive_hf_gain(&result.compensator, 20.0) ? "YES" : "No");

    /* Check specs */
    bool specs_met = lead_check_specs(&perf, &specs);
    printf("\n=== Specifications %s ===\n", specs_met ? "MET" : "NOT MET");

    return specs_met ? 0 : 1;
}