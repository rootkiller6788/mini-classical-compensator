#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/ratio_control.h"
#include "../include/ratio_station.h"
#include "../include/ratio_cascade.h"

/**
 * @file example_cascade_ratio.c
 * @brief Cascade Ratio Control — pH Neutralization with Chemical Dosing
 *
 * Demonstrates cascade ratio control where an outer pH controller
 * sets a reagent flow ratio that an inner flow controller tracks.
 *
 * Architecture:
 *   Primary:  pH controller (slow, nonlinear)
 *     → produces reagent demand signal
 *   Ratio station: demand → flow ratio
 *   Secondary: Reagent flow controller (fast, linear)
 *     → tracks ratio-adjusted flow setpoint
 *
 * L6: pH neutralization — canonical nonlinear process
 * L7: Water treatment application (NHS, Detroit Water)
 *
 * MIT 6.302: Cascade Control §9
 * Cambridge 4F2: Process control with cascade
 */

#define SIM_STEPS 40

static void print_header(const char *t) {
    printf("\n==========================================================\n");
    printf("  %s\n", t);
    printf("==========================================================\n\n");
}

int main(void) {
    print_header("pH Neutralization Cascade Ratio Control (L6/L7)");

    /* Step 1: Initialize cascade ratio controller */
    cascade_ratio_t cr;
    int rc = cascade_ratio_init(&cr, CASCADE_MODE_SIMPLE, 1,
                                 "pHC-100/FIC-101");
    if (rc != 0) {
        printf("ERROR: cascade init failed: %d\n", rc);
        return 1;
    }
    printf("[Init] Cascade controller: %s\n\n", "pHC-100/FIC-101");

    /* Step 2: Configure ratio station */
    ratio_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.ratio_gain = 0.5;   /* L/min reagent per % controller output */
    cfg.ratio_bias = 0.5;   /* minimum reagent flow */
    cfg.ratio_min  = 0.0;
    cfg.ratio_max  = 100.0;
    cfg.action     = RATIO_ACTION_MULTIPLY;
    cascade_ratio_configure(&cr, &cfg, NULL);
    printf("[Ratio] Reagent flow = %.2f * output + %.2f L/min\n",
           cfg.ratio_gain, cfg.ratio_bias);

    /* Step 3: Tune controllers */
    /* Primary (pH) — slow: Kp=0.5, Ki=0.02 */
    cascade_ratio_tune_primary(&cr, 0.5, 0.02, 0.0, 0.0, 100.0);
    /* Secondary (flow) — fast: Kp=2.0, Ki=0.5 */
    cascade_ratio_tune_secondary(&cr, 0, 2.0, 0.5, 0.0, 0.0, 100.0);

    /* Put secondary in cascade mode */
    cascade_ratio_switch_secondary_to_cascade(&cr, 0);

    printf("[Tuning] Primary: Kp=0.5, Ki=0.02 (slow pH loop)\n");
    printf("[Tuning] Secondary: Kp=2.0, Ki=0.5 (fast flow loop)\n\n");

    /* Step 4: Simulate pH neutralization */
    /* Stream: acidic influent, pH 4.5, needs NaOH dosing to pH 7.0 */
    double pH_setpoint = 7.0;
    double pH_actual = 4.5;    /* starting acidic */
    double reagent_flow = 0.0;
    double influent_flow = 100.0; /* L/min */
    double dt = 1.0; /* 1 second */

    printf("Time  | PrimSP | pH_act | PrimOut | SecSP | ReagFlow\n");
    printf("------+--------+--------+---------+-------+----------\n");

    for (int step = 0; step < SIM_STEPS; step++) {
        /* Disturbance: influent flow changes */
        if (step == 15) influent_flow = 120.0; /* 20% increase */
        if (step == 30) influent_flow = 80.0;  /* drop back */

        cascade_ratio_set_primary_sp(&cr, pH_setpoint);
        cascade_ratio_execute(&cr, pH_actual, reagent_flow, 0.0, dt);

        double primary_out = cascade_ratio_get_primary_output(&cr);
        double secondary_sp = cascade_ratio_get_secondary_sp(&cr, 0);
        (void)cascade_ratio_get_secondary_output(&cr, 0);

        /* Simulate reagent flow tracking (fast) */
        double tau_flow = 3.0; /* flow loop time constant */
        reagent_flow += (dt / tau_flow) * (secondary_sp - reagent_flow);

        /* Simulate pH response (nonlinear, simplified) */
        /* pH = 7 - (7 - 4.5) * exp(-k * reagent_flow / influent_flow) */
        double k = 0.3;
        double neutralization = exp(-k * reagent_flow / influent_flow);
        pH_actual = 7.0 - 2.5 * neutralization;
        /* Clamp to realistic range */
        if (pH_actual > 14.0) pH_actual = 14.0;
        if (pH_actual < 1.0)  pH_actual = 1.0;

        printf("%5.0f | %6.1f | %6.2f | %7.1f | %5.1f | %8.2f",
               step * dt, pH_setpoint, pH_actual,
               primary_out, secondary_sp, reagent_flow);

        if (step == 15) printf(" ← influent +20%%");
        if (step == 30) printf(" ← influent -33%%");
        printf("\n");
    }

    /* Step 5: Results analysis */
    printf("\n[Results]\n");
    printf("  Final pH: %.2f (target: %.1f)\n", pH_actual, pH_setpoint);
    printf("  Final reagent flow: %.2f L/min\n", reagent_flow);
    printf("  Steady-state error: %.3f pH units\n",
           fabs(pH_actual - pH_setpoint));

    if (fabs(pH_actual - 7.0) < 0.5) {
        printf("  Status: pH within acceptable range ✓\n");
    } else {
        printf("  Status: pH still converging (slow loop)\n");
    }

    /* Demonstrates cascade ratio principle:
     * - Primary pH loop corrects slowly but eliminates offset
     * - Secondary flow loop tracks fast, rejects flow disturbances
     * - Ratio station translates pH demand to flow setpoint
     */
    printf("\n[Knowledge Takeaway]\n");
    printf("  Cascade ratio control separates time scales:\n");
    printf("    Primary (pH):    slow, nonlinear → ratio demand\n");
    printf("    Ratio station:   demand → flow setpoint\n");
    printf("    Secondary (flow): fast, linear → rejects disturbances\n");
    printf("  This is fundamental in all process industries.\n");

    printf("\n==========================================================\n");
    printf("  pH Cascade Ratio Example — COMPLETE\n");
    printf("==========================================================\n\n");
    return 0;
}
