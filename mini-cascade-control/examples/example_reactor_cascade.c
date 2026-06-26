/**
 * @file example_reactor_cascade.c
 * @brief End-to-end example: Jacketed CSTR temperature cascade control
 *
 * L6 -- Canonical Systems: Chemical reactor with jacket/reactor temperature cascade.
 * L7 -- Applications: Process control in chemical/pharmaceutical industry.
 *
 * Demonstrates cascade control for an exothermic CSTR where jacket temperature
 * is the inner (fast) loop and reactor temperature is the outer (slow) loop.
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Jacketed CSTR Temperature Cascade Control ===\n\n");

    CascadeReactor rx = {
        .V_r = 2.0, .V_j = 0.3, .rho = 1000.0, .Cp = 4200.0,
        .UA = 800.0, .dH = -60000.0, .k0 = 5e9, .Ea = 75000.0,
        .R_gas = 8.314, .T_in = 298.0, .T_j_in = 280.0,
        .F_in = 0.0005, .F_j = 0.002
    };

    printf("Reactor parameters:\n");
    printf("  V_r=%.1f m3, V_j=%.1f m3, UA=%.0f W/K\n", rx.V_r, rx.V_j, rx.UA);
    printf("  F_in=%.4f m3/s, F_j=%.4f m3/s\n", rx.F_in, rx.F_j);

    CascadeTF jacket_tf, reactor_tf;
    cascade_create_reactor_model(&rx, &jacket_tf, &reactor_tf);

    printf("\nJacket TF (inner): Gi(s) = %.3f/(%.1f*s+1)\n",
           jacket_tf.num.coeff[0], jacket_tf.den.coeff[1]/jacket_tf.den.coeff[0]);
    printf("Reactor TF (outer): Go(s) = %.3f/(%.1f*s+1)\n",
           reactor_tf.num.coeff[0], reactor_tf.den.coeff[1]/reactor_tf.den.coeff[0]);

    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    spec.inner_pm_target = 55.0;
    spec.outer_pm_target = 40.0;
    spec.bw_ratio_min = 3.0;

    CascadeSystem sys;
    int ret = cascade_design_sequential(&jacket_tf, &reactor_tf, &spec, &sys);
    if (ret != 0) { printf("ERROR: Design failed!\n"); return 1; }

    printf("\nDesigned controllers:\n");
    printf("  Inner (jacket) PI: Kc=%.3f, Ti=%.1f s\n",
           sys.inner.inner_controller.Kp,
           sys.inner.inner_controller.Kp/(sys.inner.inner_controller.Ki+1e-10));
    printf("  Outer (reactor) PI: Kc=%.3f, Ti=%.1f s\n",
           sys.outer.outer_controller.Kp,
           sys.outer.outer_controller.Kp/(sys.outer.outer_controller.Ki+1e-10));
    printf("  Bandwidth ratio: %.2f\n", sys.bandwidth_ratio);

    CascadePerformance perf;
    cascade_assess_reactor(&rx, &perf);
    printf("\nEstimated performance:\n");
    printf("  Inner settle: %.0f s\n", perf.inner_settle_time);
    printf("  Outer settle: %.0f s\n", perf.outer_settle_time);
    printf("  Dist rejection improvement: %.1fx\n", perf.dist_rejection_ratio);

    printf("\nStep response of overall closed loop:\n");
    CascadeTF outer_cl;
    cascade_overall_closed_loop(&sys.inner.inner_cl,
                                 &sys.outer.outer_controller,
                                 &reactor_tf, &outer_cl);

    printf("  t(s)      T_response\n");
    for (int i = 0; i <= 20; i++) {
        double t = (double)i * 20.0;
        double y = cascade_step_response(&outer_cl, t, 15);
        printf("  %6.0f     %.4f\n", t, y);
    }

    cascade_tf_free(&outer_cl);
    cascade_tf_free(&jacket_tf);
    cascade_tf_free(&reactor_tf);
    cascade_system_free(&sys);

    printf("\n=== Reactor cascade example complete ===\n");
    return 0;
}
