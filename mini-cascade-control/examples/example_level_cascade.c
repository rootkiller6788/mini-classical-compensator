/**
 * @file example_level_cascade.c
 * @brief End-to-end example: Level-on-flow cascade for surge tank
 *
 * L6 -- Canonical Systems: Level tank with flow cascade.
 * L7 -- Applications: Boiler drum level, distillation column, buffer tank.
 *
 * The level process is integrating (1/(A*s)), making it inherently
 * difficult to control. Cascade with inner flow loop dramatically
 * improves disturbance rejection.
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Level-on-Flow Cascade Control for Surge Tank ===\n\n");

    CascadeLevelTank tank = {
        .tank_area = 3.0, .max_level = 4.0, .min_level = 1.0,
        .inflow_nominal = 0.02, .outflow_valve_Cv = 8.0,
        .pump_max_flow = 0.05, .level_setpoint = 2.5
    };

    printf("Tank parameters: Area=%.1f m2, Level range=%.1f-%.1f m\n",
           tank.tank_area, tank.min_level, tank.max_level);
    printf("Nominal inflow=%.3f m3/s\n", tank.inflow_nominal);

    CascadeTF flow_tf, level_tf;
    cascade_create_level_tank_model(&tank, &flow_tf, &level_tf);

    printf("\nFlow TF (inner): Gi(s) = %.3f/(%.1f*s+1)\n",
           flow_tf.num.coeff[0], flow_tf.den.coeff[1]/flow_tf.den.coeff[0]);
    printf("Level TF (outer): Go(s) = %.3f/s (integrating!)\n",
           level_tf.num.coeff[0]);

    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    spec.inner_pm_target = 60.0;
    spec.outer_pm_target = 50.0;
    spec.bw_ratio_min = 5.0;
    spec.outer_settle_max = 120.0;

    CascadeSystem sys;
    int ret = cascade_design_sequential(&flow_tf, &level_tf, &spec, &sys);
    if (ret != 0) { printf("ERROR: Design failed!\n"); return 1; }

    printf("\nDesigned controllers:\n");
    printf("  Inner (flow) PI: Kc=%.3f, Ti=%.1f s\n",
           sys.inner.inner_controller.Kp,
           sys.inner.inner_controller.Kp/(sys.inner.inner_controller.Ki+1e-10));
    printf("  Outer (level) PI: Kc=%.3f, Ti=%.1f s\n",
           sys.outer.outer_controller.Kp,
           sys.outer.outer_controller.Kp/(sys.outer.outer_controller.Ki+1e-10));
    printf("  Bandwidth ratio: %.2f\n", sys.bandwidth_ratio);

    CascadePerformance perf;
    cascade_assess_level_tank(&tank, &perf);
    printf("\nEstimated performance:\n");
    printf("  Inner settle: %.1f s (flow loop)\n", perf.inner_settle_time);
    printf("  Outer settle: %.0f s (level loop, integrating)\n", perf.outer_settle_time);
    printf("  Dist rejection improvement: %.1fx over single-loop\n",
           perf.dist_rejection_ratio);

    printf("\nKey insight: Without cascade, a 10%% inflow disturbance\n");
    printf("would cause level to drift for minutes. The inner flow loop\n");
    printf("rejects flow disturbances within seconds.\n");

    cascade_tf_free(&flow_tf);
    cascade_tf_free(&level_tf);
    cascade_system_free(&sys);

    printf("\n=== Level cascade example complete ===\n");
    return 0;
}
