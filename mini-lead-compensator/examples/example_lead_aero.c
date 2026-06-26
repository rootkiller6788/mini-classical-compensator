/**
 * @file example_lead_aero.c
 * @brief Lead compensation for aerospace pitch control (L7 Application)
 *
 * Application: Aircraft pitch angle control.
 * Plant: Short-period approximation G(s) = 1.5*(s+0.5)/(s*(s^2+0.65s+2.15))
 *
 * L7 - Real aerospace application:
 * Boeing 747 pitch dynamics approximated at Mach 0.8, 40,000 ft.
 *
 * Reference: Bryson, A.E. (1994). "Control of Spacecraft and Aircraft."
 * Stanford ENGR105, Georgia Tech AE 6530
 */

#include "lead_compensator.h"
#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(void) {
    printf("=== Lead Compensator: Aircraft Pitch Control (Boeing 747) ===\n\n");

    /* Pitch angle plant at flight condition */
    lead_system_t aircraft;
    memset(&aircraft, 0, sizeof(aircraft));
    /* G(s) = 1.5*(s+0.5) / (s*(s^2 + 0.65s + 2.15))
     *     = 1.5*s + 0.75 / (s^3 + 0.65s^2 + 2.15s) */
    aircraft.tf.num[0] = 0.75;
    aircraft.tf.num[1] = 1.5;
    aircraft.tf.num_order = 1;
    aircraft.tf.den[0] = 0.0;
    aircraft.tf.den[1] = 2.15;
    aircraft.tf.den[2] = 0.65;
    aircraft.tf.den[3] = 1.0;
    aircraft.tf.den_order = 3;
    aircraft.tf.gain = 1.0;
    strcpy(aircraft.name, "Boeing 747 Pitch");

    /* Uncompensated analysis */
    printf("Aircraft pitch dynamics (short-period approximation)\n");
    printf("G(s) = 1.5(s+0.5) / [s(s^2 + 0.65s + 2.15)]\n\n");

    double pm_raw = lead_compute_phase_margin(&aircraft);
    double gm_raw = lead_compute_gain_margin(&aircraft);
    printf("Uncompensated: PM = %.1f deg, GM = %.1f dB\n\n", pm_raw, gm_raw);

    /* Aerodynamic constraints: need PM >= 60 deg for gust rejection */
    lead_specs_t specs;
    memset(&specs, 0, sizeof(specs));
    specs.phase_margin_desired = 60.0;
    specs.gain_margin_desired = 12.0;
    specs.steady_state_error = 0.01;
    specs.max_overshoot_pct = 10.0;
    specs.settling_time = 3.0;
    specs.use_frequency_domain = true;
    specs.method = LEAD_METHOD_FREQUENCY;

    printf("Requirements:\n");
    printf("  PM >= 60 deg (gust disturbance rejection)\n");
    printf("  GM >= 12 dB (structural mode safety margin)\n");

    /* Design */
    lead_design_result_t result;
    lead_design_frequency(&aircraft, &specs, &result);

    /* Display */
    lead_print(&result.compensator);

    /* Verify */
    double pm_new = lead_compensated_phase_margin(&result.compensator, &aircraft);
    double gm_new = lead_compensated_gain_margin(&result.compensator, &aircraft);
    double w_gc = lead_open_loop_crossover(&result.compensator, &aircraft);

    printf("\nCompensated Performance:\n");
    printf("  PM = %.1f deg (requirement: >= %.0f deg)\n",
           pm_new, specs.phase_margin_desired);
    printf("  GM = %.1f dB (requirement: >= %.0f dB)\n",
           gm_new, specs.gain_margin_desired);
    printf("  Crossover = %.3f rad/s\n", w_gc);

    /* Robustness for flight envelope uncertainty */
    double M_s = lead_peak_sensitivity(&result.compensator, &aircraft,
                                       0.01, 100.0, 500);
    double delay_m = lead_delay_margin(&result.compensator, &aircraft);
    printf("  Peak Sensitivity M_s = %.3f (requirement: < 2.0)\n", M_s);
    printf("  Delay Margin = %.4f s (actuator delay tolerance)\n", delay_m);

    /* Control effort analysis */
    double max_eff = lead_max_control_effort(&result.compensator, &aircraft,
                                              0.01, 100.0, 500);
    printf("  Max Control Effort = %.3f\n", max_eff);

    /* Stability */
    bool stable = lead_is_stable_nyquist(&result.compensator, &aircraft);
    int rhp_cl = lead_count_unstable_cl_poles(&result.compensator, &aircraft);
    printf("  Nyquist Stable: %s\n", stable ? "Yes" : "No");
    printf("  RHP CL Poles: %d\n\n", rhp_cl);

    printf("=== Aerospace Lead Compensation Complete ===\n");
    return 0;
}