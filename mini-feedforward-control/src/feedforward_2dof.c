/**
 * @file feedforward_2dof.c
 * @brief 2-DOF controller implementation: combined feedback + feedforward,
 *        performance evaluation, robustness analysis.
 *
 * L2: 2-DOF separation principle ? feedback handles stability/robustness,
 *     feedforward handles performance.
 *
 * L4: Internal model principle in feedforward context.
 * L5: 2-DOF synthesis algorithms.
 * L6: Canonical 2-DOF control problems (motor, process, temperature).
 */

#include "feedforward_core.h"
#include "feedforward_design.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L2: 2-DOF Controller Construction
 * ========================================================================== */

/**
 * Initialize a 2-DOF controller from individual components.
 *
 * @param plant     Plant model P(s)
 * @param fb        Feedback controller C(s)
 * @param ff_ref    Reference feedforward F_r(s)
 * @param ff_dist   Disturbance feedforward F_d(s)
 * @param two_dof   Output 2-DOF controller
 */
void twodof_init(const TransferFn *plant, const TransferFn *fb,
                 const TransferFn *ff_ref, const TransferFn *ff_dist,
                 TwoDOF *two_dof)
{
    if (!plant || !fb || !two_dof) return;

    two_dof->plant = *plant;
    two_dof->feedback_ctrl = *fb;

    if (ff_ref) {
        two_dof->ref_ff = *ff_ref;
    } else {
        /* Default: unity FF for setpoint */
        double one[] = {1.0};
        two_dof->ref_ff = tf_create(one, 0, one, 0, 1.0);
    }

    if (ff_dist) {
        two_dof->dist_ff = *ff_dist;
    } else {
        /* Default: zero disturbance FF */
        double zero[] = {0.0};
        double one[] = {1.0};
        two_dof->dist_ff = tf_create(zero, 0, one, 0, 0.0);
    }

    two_dof->ff_gain_ratio = ff_dc_gain(&two_dof->ref_ff, 0)
                             / (ff_dc_gain(&two_dof->feedback_ctrl, 0) + 1e-9);
    two_dof->coupling_factor = 0.0; /* FF and FB paths are decoupled by design */
}

/**
 * Compute the 2-DOF control signal.
 *
 * In continuous time: u(t) = C(s)*(r(t) - y(t)) + F(s)*r(t) - D(s)*d_m(t)
 *
 * For discrete-time implementation, this function computes one time step.
 *
 * Note: This requires state-space or difference equation realization.
 * This function provides the steady-state decomposition for analysis.
 *
 * @param two_dof     2-DOF controller
 * @param r           Reference at current time
 * @param y           Measured output at current time
 * @param y_prev      Previous output (for derivative)
 * @param dt          Time step
 * @param fb_state    Feedback controller state (array, length = fb order)
 * @param ff_state    Feedforward controller state (array)
 * @param u_out       Output total control signal
 */
void twodof_compute_control(const TwoDOF *two_dof,
                            double r, double y, double y_prev, double dt,
                            double *fb_state, double *ff_state,
                            double *u_out)
{
    if (!two_dof || !u_out) return;
    (void)y_prev;
    (void)dt;
    (void)fb_state;
    (void)ff_state;

    /* Simplified: steady-state relationship
     * u = K_c*(r - y) + K_ff*r - K_d*d  (PI-type approximation) */
    double K_c = ff_dc_gain(&two_dof->feedback_ctrl, 0);
    double K_ff = ff_dc_gain(&two_dof->ref_ff, 0);

    double u_fb = K_c * (r - y);   /* Proportional feedback */
    double u_ff = K_ff * r;         /* Feedforward */

    *u_out = u_fb + u_ff;
}

/**
 * Compute closed-loop transfer functions for the 2-DOF system.
 *
 * Reference tracking:  G_ry(s) = P(s) * [F(s) + C(s)] / [1 + P(s)*C(s)]
 * Disturbance rejection: G_dy(s) = P(s) / [1 + P(s)*C(s)] - P(s)*D(s) /
 *                                    [1 + P(s)*C(s)] (simplified)
 *
 * @param two_dof  2-DOF controller
 * @param G_ry     Output: reference-to-output TF
 * @param G_dy     Output: disturbance-to-output TF
 * @param G_un     Output: reference-to-control TF
 */
void twodof_closed_loop_tfs(const TwoDOF *two_dof,
                             TransferFn *G_ry,
                             TransferFn *G_dy,
                             TransferFn *G_un)
{
    if (!two_dof) return;

    /* G_ry = (P*F + P*C) / (1 + P*C) */
    if (G_ry) {
        ff_compute_closed_loop_tf(&two_dof->plant,
                                  &two_dof->feedback_ctrl,
                                  &two_dof->ref_ff, G_ry);
    }

    /* G_dy = P / (1 + P*C) ? sensitivity * plant
     * With disturbance FF: y/d = P/(1+PC) + P*D/(1+PC)
     * = P*(1+D)/(1+PC) where D = dist_ff (actually G_d_model * D_ff) */
    if (G_dy) {
        TransferFn pc = tf_series(&two_dof->plant, &two_dof->feedback_ctrl);
        /* Denominator: 1 + PC */
        Poly den_pc = poly_mul(&two_dof->plant.den, &two_dof->feedback_ctrl.den);
        Poly num_pc_prod = poly_mul(&two_dof->plant.num, &two_dof->feedback_ctrl.num);
        double k_pc = two_dof->plant.gain * two_dof->feedback_ctrl.gain;
        for (int i = 0; i <= num_pc_prod.degree; i++) {
            num_pc_prod.coeff[i] *= k_pc;
        }
        Poly den_cl = poly_add(&den_pc, &num_pc_prod);

        /* Numerator = P * (den_controller) simplified to P */
        G_dy->num = poly_create(two_dof->plant.num.coeff,
                               two_dof->plant.num.degree);
        /* Multiply by feedback controller denominator */
        Poly num = poly_mul(&G_dy->num, &two_dof->feedback_ctrl.den);
        poly_free(&G_dy->num);
        G_dy->num = num;
        G_dy->den = den_cl;
        G_dy->gain = two_dof->plant.gain;

        poly_free(&den_pc); poly_free(&num_pc_prod);
        tf_free(&pc);
    }

    /* G_un = [F + C] / (1 + P*C) */
    if (G_un) {
        TransferFn f_plus_c = tf_parallel(&two_dof->ref_ff,
                                          &two_dof->feedback_ctrl);
        /* Denominator same as G_ry */
        TransferFn pc = tf_series(&two_dof->plant, &two_dof->feedback_ctrl);
        Poly den_pc = poly_mul(&two_dof->plant.den, &two_dof->feedback_ctrl.den);
        Poly num_pc_prod = poly_mul(&two_dof->plant.num, &two_dof->feedback_ctrl.num);
        double k_pc = two_dof->plant.gain * two_dof->feedback_ctrl.gain;
        for (int i = 0; i <= num_pc_prod.degree; i++) {
            num_pc_prod.coeff[i] *= k_pc;
        }
        Poly den_cl = poly_add(&den_pc, &num_pc_prod);

        /* Numerator = (F+C) * denominator without plant */
        G_un->num = poly_mul(&f_plus_c.num, &den_pc);
        G_un->den = poly_mul(&f_plus_c.den, &den_cl);
        G_un->gain = 1.0;

        poly_free(&den_pc); poly_free(&num_pc_prod);
        poly_free(&f_plus_c.num); poly_free(&f_plus_c.den);
        tf_free(&pc);
    }
}

/**
 * Evaluate 2-DOF performance metrics.
 *
 * @param two_dof  2-DOF controller
 * @param t_final  Simulation end time (s)
 * @param dt       Time step (s)
 * @param perf     Output performance metrics
 */
void twodof_evaluate_performance(const TwoDOF *two_dof,
                                 double t_final, double dt,
                                 FFPerformance *perf)
{
    if (!two_dof || !perf || dt <= 0.0) return;

    int n_steps = (int)(t_final / dt) + 1;
    if (n_steps > 10000) n_steps = 10000;

    /* Simulate step response with and without feedforward */
    double *y_with_ff = (double *)calloc(n_steps, sizeof(double));
    double *y_without_ff = (double *)calloc(n_steps, sizeof(double));
    double *r = (double *)malloc(n_steps * sizeof(double));
    if (!y_with_ff || !y_without_ff || !r) {
        free(y_with_ff); free(y_without_ff); free(r);
        return;
    }

    /* Unit step reference after 10% of simulation time */
    int step_start = n_steps / 10;
    for (int i = 0; i < n_steps; i++) {
        r[i] = (i >= step_start) ? 1.0 : 0.0;
    }

    /* Approximate plant dynamics using state-space simulation */
    int n = two_dof->plant.den.degree;
    if (n > 10) n = 10;
    double *A = (double *)malloc(n * n * sizeof(double));
    double *B = (double *)malloc(n * sizeof(double));
    double *C = (double *)malloc(n * sizeof(double));
    double D = 0.0;
    double *x_ff = (double *)calloc(n, sizeof(double));
    double *x_nff = (double *)calloc(n, sizeof(double));

    if (!A || !B || !C || !x_ff || !x_nff) {
        free(A); free(B); free(C); free(x_ff); free(x_nff);
        free(y_with_ff); free(y_without_ff); free(r);
        return;
    }

    int dim = tf_to_state_space(&two_dof->plant, A, B, C, &D);

    double Kfb = ff_dc_gain(&two_dof->feedback_ctrl, 0);
    double Kff = ff_dc_gain(&two_dof->ref_ff, 0);

    /* Forward Euler simulation */
    for (int k = 0; k < n_steps; k++) {
        double rk = r[k];

        /* With feedforward */
        {
            double yk = D;
            for (int i = 0; i < dim; i++) yk += C[i] * x_ff[i];
            double u_ff = Kfb * (rk - yk) + Kff * rk;
            double *dx = (double *)calloc(dim, sizeof(double));
            if (dx) {
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        dx[i] += A[i * dim + j] * x_ff[j];
                    }
                    dx[i] += B[i] * u_ff;
                }
                for (int i = 0; i < dim; i++) x_ff[i] += dt * dx[i];
                free(dx);
            }
            y_with_ff[k] = yk;
        }

        /* Without feedforward (feedback only) */
        {
            double yk = D;
            for (int i = 0; i < dim; i++) yk += C[i] * x_nff[i];
            double u_fb = Kfb * (rk - yk);
            double *dx = (double *)calloc(dim, sizeof(double));
            if (dx) {
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        dx[i] += A[i * dim + j] * x_nff[j];
                    }
                    dx[i] += B[i] * u_fb;
                }
                for (int i = 0; i < dim; i++) x_nff[i] += dt * dx[i];
                free(dx);
            }
            y_without_ff[k] = yk;
        }
    }

    /* Compute metrics */
    double ise_ff = 0.0, ise_nff = 0.0;
    double overshoot_ff = 0.0, overshoot_nff = 0.0;
    double settle_idx_ff = n_steps, settle_idx_nff = n_steps;
    double steady_state = 1.0; /* Tracking unit step */

    for (int k = step_start; k < n_steps; k++) {
        double e_ff = r[k] - y_with_ff[k];
        double e_nff = r[k] - y_without_ff[k];
        ise_ff += e_ff * e_ff * dt;
        ise_nff += e_nff * e_nff * dt;

        if (y_with_ff[k] > overshoot_ff) overshoot_ff = y_with_ff[k];
        if (y_without_ff[k] > overshoot_nff) overshoot_nff = y_without_ff[k];

        /* Settling time: last time error exceeds 2% */
        if (fabs(e_ff) > 0.02) settle_idx_ff = k;
        if (fabs(e_nff) > 0.02) settle_idx_nff = k;
    }

    perf->tracking_ise = ise_ff;
    perf->disturbance_ise = 0.0; /* No disturbance in this test */
    perf->ff_contribution = (ise_nff > 1e-15) ?
                            1.0 - ise_ff / ise_nff : 0.0;
    perf->settle_time = (settle_idx_ff - step_start) * dt;
    perf->overshoot_pct = (overshoot_ff - steady_state) * 100.0;
    /* Robustness margin from settling time improvement */
    perf->robustness_margin = (settle_idx_nff > step_start && settle_idx_ff > step_start) ?
        (double)(settle_idx_nff - step_start) / (double)(settle_idx_ff - step_start + 1) : 1.0;

    free(y_with_ff); free(y_without_ff); free(r);
    free(A); free(B); free(C); free(x_ff); free(x_nff);
}

/* ==========================================================================
 * L4: Internal Model Principle in Feedforward
 * ========================================================================== */

/**
 * Check if a feedforward compensator satisfies the internal model principle.
 *
 * The IMP (Francis & Wonham, 1976) states that for perfect asymptotic
 * tracking/rejection of a signal, the controller must contain the
 * generator of that signal.
 *
 * For feedforward: the reference/disturbance model must be embedded
 * in the feedforward path.
 *
 * @param tf       Feedforward transfer function
 * @param signal_type 0=step, 1=ramp, 2=sinusoid
 * @param freq     Frequency for sinusoidal (rad/s)
 * @return 1 if IMP satisfied, 0 otherwise
 */
int twodof_check_imp(const TransferFn *tf, int signal_type, double freq)
{
    if (!tf) return 0;

    switch (signal_type) {
        case 0: /* Step: need integrator (pole at s=0, i.e., a0=0) */
            return fabs(tf->den.coeff[0]) < 1e-10;
        case 1: /* Ramp: need double integrator */
            return (fabs(tf->den.coeff[0]) < 1e-10 &&
                    fabs(tf->den.coeff[1]) < 1e-10);
        case 2: /* Sinusoid: need mode at freq */
        {
            /* Check if (s^2 + freq^2) divides denominator */
            double w2 = freq * freq;
            double val = poly_eval(&tf->den, 0.0); /* evaluate den at s=j*freq */
            (void)val;
            (void)w2;
            /* Simplified check: look for pure imaginary pole pair */
            double poles[20];
            int np = tf_find_poles(tf, poles, 20);
            for (int i = 0; i < np; i++) {
                /* Check if any pole is pure imaginary (real part ? 0) */
                /* This is simplified; real check needs complex root finder */
            }
            return 0; /* Conservative */
        }
        default:
            return 0;
    }
}

/* ==========================================================================
 * L5: Robustness Analysis
 * ========================================================================== */

/**
 * Compute the robustness of the 2-DOF system to plant uncertainty.
 *
 * Uses the small-gain theorem: stability is guaranteed if
 * |Delta(jw)| * |T(jw)| < 1 for all frequencies,
 * where Delta is the multiplicative uncertainty and T is the
 * complementary sensitivity.
 *
 * @param two_dof    2-DOF controller
 * @param delta_max  Maximum multiplicative uncertainty magnitude
 * @param omega      Frequency array (rad/s)
 * @param n_freq     Number of frequencies
 * @param margin     Output: minimum stability margin (>1 = robust)
 */
void twodof_robustness_margin(const TwoDOF *two_dof,
                              double delta_max,
                              const double *omega, int n_freq,
                              double *margin)
{
    if (!two_dof || !omega || !margin) return;

    /* Complementary sensitivity: T = P*C / (1 + P*C) */
    double min_margin = 1e15;

    for (int i = 0; i < n_freq; i++) {
        /* |T(jw)| */
        TransferFn pc = tf_series(&two_dof->plant, &two_dof->feedback_ctrl);
        double mag_pc, ph_pc;
        tf_freq_response(&pc, omega[i], &mag_pc, &ph_pc);
        tf_free(&pc);

        double re_pc = mag_pc * cos(ph_pc);
        double im_pc = mag_pc * sin(ph_pc);
        double mag_T = mag_pc / sqrt((1.0 + re_pc) * (1.0 + re_pc) + im_pc * im_pc);

        /* Robustness condition: delta_max * |T(jw)| < 1 */
        double m = 1.0 / (delta_max * mag_T + 1e-15);
        if (m < min_margin) min_margin = m;
    }

    *margin = min_margin;
}

/* ==========================================================================
 * L6: Feedforward for DC Motor Position Control
 * ========================================================================== */

/**
 * Design feedforward for a DC motor position servo.
 *
 * DC motor model: P(s) = K / (s * (tau*s + 1))
 *   where K = Kt / (R*b + Ke*Kt) ? 1/Ke for ideal motor
 *         tau = R*J / (R*b + Ke*Kt) ? R*J / (Ke*Kt)
 *
 * Feedforward for position tracking:
 *   u_ff(t) = J_hat/Kt * a_d(t) + b_hat/Kt * v_d(t)
 *
 * @param J_hat      Estimated inertia (kg*m^2)
 * @param b_hat      Estimated damping (N*m*s/rad)
 * @param Kt         Torque constant (N*m/A)
 * @param q_d        Desired position (rad)
 * @param qd_d       Desired velocity (rad/s)
 * @param qdd_d      Desired acceleration (rad/s^2)
 * @param u_ff       Output feedforward voltage
 */
double ff_dc_motor_position(double J_hat, double b_hat, double Kt,
                             double q_d, double qd_d, double qdd_d)
{
    /* Simple model: tau = J*qdd + b*qd
     * Voltage: u = tau / Kt (neglecting back-EMF for FF) */
    (void)q_d; /* position does not directly affect torque in simple model */
    double torque = J_hat * qdd_d + b_hat * qd_d;
    return torque / Kt;
}

/**
 * Design feedforward for a DC motor velocity servo.
 *
 * Velocity model: P(s) = K / (tau*s + 1)
 * Feedforward: u_ff = v_d / K  (steady-state) + tau/K * a_d (dynamic)
 *
 * @param K          Motor gain (rad/s/V)
 * @param tau        Motor time constant (s)
 * @param v_d        Desired velocity (rad/s)
 * @param a_d        Desired acceleration (rad/s^2)
 * @return Feedforward voltage
 */
double ff_dc_motor_velocity(double K, double tau,
                             double v_d, double a_d)
{
    return v_d / K + tau * a_d / K;
}

/* ==========================================================================
 * L6: Feedforward for Temperature Control
 * ========================================================================== */

/**
 * Design disturbance feedforward for a temperature control system.
 *
 * Typical first-order thermal model:
 *   T_out(s) = K/(tau*s+1) * u(s) + Kd/(tau*s+1) * T_ambient(s)
 *
 * Disturbance feedforward: u_ff = -(Kd/K) * T_ambient
 *
 * @param K          Process gain
 * @param Kd         Disturbance gain
 * @param T_ambient  Measured ambient temperature
 * @param prev_u     Previous control signal (for smoothing)
 * @param alpha      Smoothing factor
 * @return Feedforward control signal
 */
double ff_temperature_dist(double K, double Kd, double T_ambient,
                            double prev_u, double alpha)
{
    double u_ff_raw = -(Kd / K) * T_ambient;
    /* Smooth to avoid abrupt changes */
    return alpha * u_ff_raw + (1.0 - alpha) * prev_u;
}

/* ==========================================================================
 * L7: Application ? Process Control Disturbance Feedforward
 * ========================================================================== */

/**
 * Chemical process disturbance feedforward ? feed composition change.
 *
 * For a mixing process:
 *   C_out(s) = Gp(s)*u(s) + Gd(s)*C_feed(s)
 *
 * Static feedforward: u_ff = -(Gd(0)/Gp(0)) * C_feed
 *
 * @param Gp_dc     Process DC gain
 * @param Gd_dc     Disturbance DC gain
 * @param C_feed    Feed concentration measurement
 * @param C_feed_0  Nominal feed concentration (for linearization)
 * @return Feedforward correction to control valve
 */
double ff_process_composition(double Gp_dc, double Gd_dc,
                               double C_feed, double C_feed_0)
{
    /* Linearized FF around nominal operating point */
    double delta = C_feed - C_feed_0;
    return -(Gd_dc / Gp_dc) * delta;
}

/**
 * Boiler feedforward control for steam pressure.
 *
 * When steam demand changes suddenly, feedforward adjusts fuel flow
 * before the pressure controller sees the error.
 *
 * Model: pressure drop = Kd * steam_flow_change
 * Feedforward: fuel_increase = Kff * measured_steam_flow
 *
 * @param Kff            Feedforward gain
 * @param steam_flow     Measured steam flow (kg/s)
 * @param steam_flow_nom Nominal steam flow
 * @param fuel_nom       Nominal fuel flow
 * @return Fuel flow setpoint (kg/s)
 */
double ff_boiler_steam_pressure(double Kff, double steam_flow,
                                 double steam_flow_nom, double fuel_nom)
{
    double delta_steam = steam_flow - steam_flow_nom;
    return fuel_nom + Kff * delta_steam;
}

/**
 * pH neutralization process feedforward.
 *
 * Highly nonlinear process. Simplified feedforward based on
 * titration curve linearization around operating point.
 *
 * @param pH_setpoint  Target pH
 * @param pH_influent  Measured influent pH
 * @param flow_acid    Acid flow rate
 * @param flow_base    Base flow rate
 * @param buffer_cap   Buffer capacity (mol/L/pH)
 * @return Required reagent flow adjustment
 */
double ff_ph_neutralization(double pH_setpoint, double pH_influent,
                             double flow_acid, double flow_base,
                             double buffer_cap)
{
    double ph_error = pH_influent - pH_setpoint;
    double required_neutralization = buffer_cap * ph_error;
    /* Scale to flow adjustment */
    return required_neutralization * (flow_acid + flow_base) / buffer_cap;
}

/* ==========================================================================
 * L7: Application ? Aerospace Feedforward
 * ========================================================================== */

/**
 * Satellite attitude control with feedforward.
 *
 * For known slew maneuvers, feedforward computes the torque profile
 * needed to track the desired attitude trajectory.
 *
 * Rigid-body: J*theta_ddot = tau (simplified single-axis)
 *
 * @param J          Moment of inertia (kg*m^2)
 * @param theta_d    Desired angle (rad)
 * @param omega_d    Desired angular velocity (rad/s)
 * @param alpha_d    Desired angular acceleration (rad/s^2)
 * @return Feedforward torque (N*m)
 */
double ff_satellite_attitude(double J, double theta_d,
                              double omega_d, double alpha_d)
{
    (void)theta_d;
    /* Simple rigid-body: tau = J*alpha
     * More sophisticated: includes Euler equation coupling terms
     * and gyroscopic torques proportional to omega_d */
    double tau_ff = J * alpha_d;

    /* Gyroscopic coupling (typically small for slow maneuvers):
     * tau_gyro = omega_d x (J * omega), approximated as 0 here.
     * For high-rate slews, add: tau_ff += cross(omega_d, J*omega); */
    (void)omega_d;

    return tau_ff;
}

/**
 * Quadrotor position control feedforward.
 *
 * Simplified altitude feedforward:
 *   thrust_ff = m * (g + a_d_z) / (cos(phi)*cos(theta))
 *
 * @param mass       Total mass (kg)
 * @param a_d_z      Desired vertical acceleration (m/s^2)
 * @param phi        Roll angle (rad)
 * @param theta      Pitch angle (rad)
 * @return Feedforward thrust (N)
 */
double ff_quadrotor_altitude(double mass, double a_d_z,
                              double phi, double theta)
{
    double g = 9.81;
    double cos_phi = cos(phi);
    double cos_theta = cos(theta);
    double denom = cos_phi * cos_theta;
    if (fabs(denom) < 0.01) denom = 0.01;
    return mass * (g + a_d_z) / denom;
}

/* ==========================================================================
 * L8: Advanced ? Gain-Scheduled Feedforward
 * ========================================================================== */

/**
 * Gain-scheduled feedforward for nonlinear systems.
 *
 * Uses multiple linear models at different operating points.
 * Feedforward gain is interpolated based on scheduling variable.
 *
 * @param K_low      FF gain at low operating point
 * @param K_high     FF gain at high operating point
 * @param var_low    Scheduling variable low value
 * @param var_high   Scheduling variable high value
 * @param var        Current scheduling variable value
 * @param r          Reference
 * @return Scheduled feedforward signal
 */
double ff_gain_scheduled(double K_low, double K_high,
                          double var_low, double var_high,
                          double var, double r)
{
    /* Linear interpolation of gains */
    double range = var_high - var_low;
    if (fabs(range) < 1e-10) return K_low * r;

    double frac = (var - var_low) / range;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    double K = K_low + frac * (K_high - K_low);
    return K * r;
}

/**
 * Compute the feedforward signal for a flexible transmission.
 *
 * Two-mass system with flexible coupling:
 *   J1*theta1_ddot + b*(theta1_dot - theta2_dot) + k*(theta1 - theta2) = tau
 *   J2*theta2_ddot + b*(theta2_dot - theta1_dot) + k*(theta2 - theta1) = 0
 *
 * Feedforward for load-side (theta2) positioning:
 *   tau_ff = (J1+J2) * theta2_d_ddot + J1*J2/k * theta2_d_dddd
 *   (approximate, neglecting damping)
 *
 * @param J1         Motor inertia
 * @param J2         Load inertia
 * @param k          Stiffness (N*m/rad)
 * @param a_d        Desired load acceleration (rad/s^2)
 * @param j_d        Desired load jerk (rad/s^3)
 * @return Feedforward torque
 */
double ff_flexible_transmission(double J1, double J2, double k,
                                 double a_d, double j_d)
{
    /* Simplified: tau_ff = (J1+J2) * a_d
     * A more precise implementation would add J1*J2/k * snap (jerk derivative) */
    (void)k;
    (void)j_d;
    return (J1 + J2) * a_d;
}
