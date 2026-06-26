/**
 * @file example_dc_motor_cascade.c
 * @brief End-to-end example: DC motor position/velocity cascade control
 *
 * L6 -- Canonical Systems: DC motor with cascade position/velocity control.
 * L7 -- Applications: Servo motor control (robotics, CNC, aerospace actuators).
 *
 * This example demonstrates:
 *   1. DC motor model creation from physical parameters
 *   2. Inner velocity loop PI design
 *   3. Outer position loop P design
 *   4. Cascade system simulation with step response
 *   5. Performance evaluation
 *
 * Physical parameters based on a typical small servo motor.
 *
 * Course alignment:
 *   MIT 6.302 -- DC motor control
 *   Stanford ENGR105 -- Servo control
 *   Berkeley ME232 -- Electromechanical systems
 *   Purdue ME 575 -- Industrial motion control
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== DC Motor Position/Velocity Cascade Control ===\n\n");

    /* Define DC motor parameters (small servo motor) */
    CascadeDCMotor motor = {
        .R = 2.0,      /* Armature resistance, Ohm */
        .L = 0.002,    /* Armature inductance, H */
        .Kb = 0.05,    /* Back-EMF constant, V*s/rad */
        .Kt = 0.05,    /* Torque constant, N*m/A */
        .J = 0.0001,   /* Rotor inertia, kg*m^2 */
        .B = 0.00005,  /* Viscous friction, N*m*s/rad */
        .gear_ratio = 1.0
    };

    printf("Motor parameters:\n");
    printf("  R=%.3f Ohm, L=%.4f H, Kt=%.3f Nm/A, J=%.6f kg*m2\n",
           motor.R, motor.L, motor.Kt, motor.J);

    /* Create motor transfer functions */
    CascadeTF vel_tf, pos_tf;
    cascade_create_dc_motor_model(&motor, &vel_tf, &pos_tf);

    printf("\nVelocity TF (inner loop): Gi(s) = Kt/(J*L*s^2 + ...)\n");
    printf("  deg(num)=%d, deg(den)=%d\n", vel_tf.num.degree, vel_tf.den.degree);

    printf("Position TF (outer loop): Go(s) = 1/s\n");
    printf("  deg(num)=%d, deg(den)=%d\n", pos_tf.num.degree, pos_tf.den.degree);

    /* Design specifications */
    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    spec.inner_pm_target = 65.0;   /* Higher PM for velocity loop */
    spec.inner_bw_min = 50.0;      /* rad/s -- fast velocity control */
    spec.outer_settle_max = 0.5;   /* 0.5s settling for position */

    printf("\nDesign specifications:\n");
    printf("  Inner PM target: %.0f deg, BW min: %.0f rad/s\n",
           spec.inner_pm_target, spec.inner_bw_min);
    printf("  Outer settle max: %.1f s\n", spec.outer_settle_max);

    /* Design cascade controller */
    CascadeSystem sys;
    int ret = cascade_design_sequential(&vel_tf, &pos_tf, &spec, &sys);
    if (ret != 0) {
        printf("ERROR: Cascade design failed!\n");
        cascade_tf_free(&vel_tf);
        cascade_tf_free(&pos_tf);
        return 1;
    }

    printf("\nDesigned controllers:\n");
    printf("  Inner PI: Kp=%.3f, Ki=%.3f\n",
           sys.inner.inner_controller.Kp, sys.inner.inner_controller.Ki);
    printf("  Outer PI: Kp=%.3f, Ki=%.3f\n",
           sys.outer.outer_controller.Kp, sys.outer.outer_controller.Ki);
    printf("  Bandwidth ratio: %.2f\n", sys.bandwidth_ratio);

    /* Assess performance */
    CascadePerformance perf;
    cascade_assess_dc_motor(&motor, &perf);

    printf("\nEstimated performance:\n");
    printf("  Inner settle: %.3f s, overshoot: %.1f%%\n",
           perf.inner_settle_time, perf.inner_overshoot * 100.0);
    printf("  Outer settle: %.3f s, overshoot: %.1f%%\n",
           perf.outer_settle_time, perf.outer_overshoot * 100.0);
    printf("  Dist rejection improvement: %.1fx\n", perf.dist_rejection_ratio);
    printf("  Robustness margin (Ms eqv): %.2f\n", perf.robustness_margin);

    /* Step response simulation */
    printf("\nStep response of outer closed loop:\n");
    printf("  t(s)      y(t)\n");
    printf("  ---------  --------\n");

    CascadeTF outer_cl;
    cascade_overall_closed_loop(&sys.inner.inner_cl,
                                 &sys.outer.outer_controller,
                                 &pos_tf, &outer_cl);

    for (int i = 0; i <= 20; i++) {
        double t = (double)i * 0.05;  /* 0 to 1 second */
        double y = cascade_step_response(&outer_cl, t, 15);
        printf("  %.3f      %.4f\n", t, y);
    }

    /* Compute ISE */
    double ise = 0.0;
    int n_pts = 200;
    for (int i = 0; i < n_pts; i++) {
        double t = (double)i * 1.0 / (double)(n_pts - 1);
        double y = cascade_step_response(&outer_cl, t, 15);
        double e = 1.0 - y;
        double dt = 1.0 / (double)n_pts;
        ise += e * e * dt;
    }
    printf("\n  ISE = %.6f\n", ise);

    /* Cleanup */
    cascade_tf_free(&outer_cl);
    cascade_tf_free(&vel_tf);
    cascade_tf_free(&pos_tf);
    cascade_system_free(&sys);

    printf("\n=== Example complete ===\n");
    return 0;
}
