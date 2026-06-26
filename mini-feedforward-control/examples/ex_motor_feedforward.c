/**
 * @file ex_motor_feedforward.c
 * @brief Example: DC motor position control with velocity/acceleration feedforward.
 *
 * Demonstrates L6 canonical problem: servo motor tracking with 2-DOF control.
 * Compares feedback-only vs feedback+feedforward performance.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "feedforward_core.h"
#include "feedforward_design.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * DC motor parameters (realistic values for a small servo):
 *   J  = 0.0001  kg*m^2  (rotor inertia)
 *   b  = 0.00005 N*m*s/rad (viscous friction)
 *   Kt = 0.05    N*m/A  (torque constant)
 *   Ke = 0.05    V*s/rad (back-EMF constant, ~Kt in SI)
 *   R  = 2.0     Ohm    (winding resistance)
 *   L  = 0.001   H      (winding inductance)
 *
 * Simplified model (neglecting inductance):
 *   P(s) = K / (s * (tau*s + 1))
 *   where K = Kt / (R*b + Ke*Kt) ? 1/Ke for b small
 *         tau = R*J / (R*b + Ke*Kt) ? R*J/(Ke*Kt)
 */

int main(void)
{
    printf("=== DC Motor Position Control with Feedforward ===\n\n");

    /* Motor parameters */
    double J  = 0.0001;
    double b  = 0.00005;
    double Kt = 0.05;
    double Ke = 0.05;
    double R  = 2.0;

    /* Derived parameters */
    double denom = R * b + Ke * Kt;
    double K = Kt / denom;
    double tau = R * J / denom;

    printf("Motor model: P(s) = %.4f / (s * (%.4f*s + 1))\n", K, tau);
    printf("  DC gain (velocity): K = %.4f rad/s/V\n", K);
    printf("  Time constant: tau = %.6f s\n\n", tau);

    /* Position transfer function: P(s) = K/(s*(tau*s+1)) */
    double num_p[] = {K};
    double den_p[] = {0.0, 1.0, tau};  /* s*(tau*s+1) = s + tau*s^2 */
    TransferFn plant = tf_create(num_p, 0, den_p, 2, 1.0);

    /* Feedback controller: PD control C(s) = Kp + Kd*s
     * Designed for bandwidth of ~50 rad/s, damping ~0.7 */
    double Kp = 250.0;  /* Proportional gain */
    double Kd = 30.0;   /* Derivative gain */
    double nc[] = {Kp, Kd};
    double dc[] = {1.0, 0.01};  /* pseudo-derivative with LPF */
    TransferFn fb_ctrl = tf_create(nc, 1, dc, 1, 1.0);

    /* Feedforward: velocity + acceleration
     * F(s) = Kv_ff*s + Ka_ff*s^2 (with LPF for realizability)
     * Kv_ff ? 1/K (ideal velocity FF)
     * Ka_ff ? tau/K (acceleration FF, compensates for inertia) */
    double Kv_ff = 1.0 / K;
    double Ka_ff = tau / K;
    printf("Feedforward gains:\n");
    printf("  Velocity FF:     Kv_ff = %.4f V/(rad/s)\n", Kv_ff);
    printf("  Acceleration FF: Ka_ff = %.6f V/(rad/s^2)\n", Ka_ff);

    /* Build feedforward filter: F(s) = (Kv_ff*s + Ka_ff*s^2) / (s+wc)^2
     * with wc high enough for properness */
    double wc = 500.0;  /* Filter bandwidth */
    double nf[] = {0.0, Kv_ff, Ka_ff};
    double df[] = {wc*wc, 2.0*wc, 1.0};
    TransferFn ff_ref = tf_create(nf, 2, df, 2, 1.0);

    /* Evaluate */
    printf("\n--- Reference Tracking Simulation ---\n");

    /* Desired trajectory: ramp to pi/2 (90 degrees) in 0.5s, hold */
    double t_final = 1.0;
    double dt = 0.001;
    int n = (int)(t_final / dt) + 1;

    /* Simulation state */
    double pos_fb = 0.0, vel_fb = 0.0;
    double pos_ff = 0.0, vel_ff = 0.0;
    double u_fb = 0.0, u_ff = 0.0;
    double target = M_PI / 2.0;
    double rise_time = 0.3;

    /* Trajectory generation: smooth S-curve */
    double r, r_dot, r_ddot;

    double ise_fb = 0.0, ise_ff = 0.0;
    double max_u_fb = 0.0, max_u_ff = 0.0;

    for (int k = 0; k < n; k++) {
        double t = k * dt;

        /* Desired trajectory: 5th-order polynomial (minimum jerk) */
        if (t < rise_time) {
            double tau_n = t / rise_time;
            r = target * (10.0*pow(tau_n,3) - 15.0*pow(tau_n,4) + 6.0*pow(tau_n,5));
            r_dot = target * (30.0*pow(tau_n,2) - 60.0*pow(tau_n,3) + 30.0*pow(tau_n,4)) / rise_time;
            r_ddot = target * (60.0*tau_n - 180.0*pow(tau_n,2) + 120.0*pow(tau_n,3)) / (rise_time*rise_time);
        } else {
            r = target;
            r_dot = 0.0;
            r_ddot = 0.0;
        }

        if (k % 200 == 0) {
            printf("  t=%.3f: r=%.4f, r_dot=%.4f, r_ddot=%.4f\n", t, r, r_dot, r_ddot);
        }

        /* --- Feedback only --- */
        double err_fb = r - pos_fb;
        u_fb = Kp * err_fb - Kd * vel_fb;
        if (fabs(u_fb) > max_u_fb) max_u_fb = fabs(u_fb);

        /* Plant dynamics: d(pos)/dt = vel, d(vel)/dt = (K*u - vel)/tau */
        double acc_fb = (K * u_fb - vel_fb) / tau;
        vel_fb += acc_fb * dt;
        pos_fb += vel_fb * dt;

        double err_curr_fb = r - pos_fb;
        ise_fb += err_curr_fb * err_curr_fb * dt;

        /* --- With feedforward --- */
        double err_ff = r - pos_ff;
        u_ff = Kp * err_ff - Kd * vel_ff + Kv_ff * r_dot + Ka_ff * r_ddot;
        if (fabs(u_ff) > max_u_ff) max_u_ff = fabs(u_ff);

        double acc_ff = (K * u_ff - vel_ff) / tau;
        vel_ff += acc_ff * dt;
        pos_ff += vel_ff * dt;

        double err_curr_ff = r - pos_ff;
        ise_ff += err_curr_ff * err_curr_ff * dt;

        /* Print result at end */
        if (k == n - 1) {
            printf("\n  Final state:\n");
            printf("    FB only:  pos=%.6f rad, ISE=%.6f\n", pos_fb, ise_fb);
            printf("    FB + FF:  pos=%.6f rad, ISE=%.6f\n", pos_ff, ise_ff);
        }
    }

    printf("\n--- Performance Summary ---\n");
    printf("  Feedback only:\n");
    printf("    Position error (final): %.6f rad\n", target - pos_fb);
    printf("    ISE:                   %.6f\n", ise_fb);
    printf("    Max control:           %.4f V\n", max_u_fb);
    printf("\n  Feedback + Feedforward:\n");
    printf("    Position error (final): %.6f rad\n", target - pos_ff);
    printf("    ISE:                   %.6f\n", ise_ff);
    printf("    Max control:           %.4f V\n", max_u_ff);
    printf("    ISE improvement:       %.1f%%\n",
           100.0 * (1.0 - ise_ff / (ise_fb + 1e-15)));

    tf_free(&plant);
    tf_free(&fb_ctrl);
    tf_free(&ff_ref);

    printf("\n=== Done ===\n");
    return 0;
}
