/**
 * @file lag_simulation.c
 * @brief Time-domain simulation of lag-compensated systems
 *
 * Simulates step, ramp, and general input responses for
 * systems with lag compensators in the forward path.
 *
 * L5: Numerical ODE simulation (Euler, RK4), step/ramp response
 * L6: DC motor speed control, temperature control, position servo
 *
 * Textbook: Ogata Ch. 5; Franklin/Powell/Emami-Naeini Ch. 3
 */

#include "lag_compensator.h"
#include "lag_types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==========================================================================
 * L5: State-space model for simulation
 * ========================================================================== */

/**
 * Plant state-space representation for first-order plant:
 *   G(s) = K / (tau*s + 1)
 *
 * State equation:  dx/dt = (-1/tau)*x + (K/tau)*u
 * Output equation: y = x
 *
 * Combined with lag compensator G_c(s) = Kc*(T*s+1)/(beta*T*s+1):
 *   Compensator state: dz/dt = (-1/(beta*T))*z + u_c
 *   Compensator output: u = (Kc*(beta-1)/(beta^2*T))*z + (Kc/beta)*u_c
 *
 * Wait — let me derive the compensator state-space realization properly.
 *
 * G_c(s) = U(s)/E(s) = Kc*(T*s+1)/(beta*T*s+1)
 *
 * Observability canonical form (controller form):
 *   dz/dt = -1/(beta*T) * z + e
 *   u = (Kc*(beta-1)/(beta^2*T)) * z + (Kc/beta) * e
 *
 * Let me re-derive:
 *   G_c(s) = Kc*(Ts+1)/(beta*Ts+1) = Kc/beta + (Kc*(beta-1)/beta)/(beta*Ts+1)
 *
 * So: u = (Kc/beta)*e + u_z
 *     where u_z comes from: (beta*Ts+1)*u_z = (Kc*(beta-1)/beta)*e
 *     dz/dt = -1/(beta*T)*z + e
 *     u_z = (Kc*(beta-1)/(beta*T))*z  — wait, let me use proper state-space.
 *
 * Controllable canonical form for G_c(s) = (b1*s + b0)/(a1*s + a0):
 *   b1=Kc*T, b0=Kc, a1=beta*T, a0=1
 *   dz/dt = -a0/a1 * z + e = -1/(beta*T)*z + e
 *   u = (b0 - b1*a0/a1)*z + (b1/a1)*e
 *     = (Kc - Kc*T*1/(beta*T))*z + (Kc*T/(beta*T))*e
 *     = Kc*(1 - 1/beta)*z + (Kc/beta)*e
 *     = Kc*(beta-1)/beta * z + Kc/beta * e
 */

/**
 * Simulate step response of unity-feedback system with lag compensator
 * and first-order plant.
 *
 * G(s) = K / (tau*s + 1)
 * G_c(s) = Kc*(T*s+1)/(beta*T*s+1)
 *
 * Closed-loop: T_cl(s) = G_c*G / (1 + G_c*G)
 *
 * State equations (2 states: plant x, compensator z):
 *   dx/dt = -1/tau * x + (K/tau) * u
 *   dz/dt = -1/(beta*T) * z + e
 *   e = r - x   (unity feedback)
 *   u = Kc*(beta-1)/beta * z + Kc/beta * e
 *
 * Combined:
 *   dx/dt = -1/tau*x + (K/tau)*[Kc*(beta-1)/beta*z + Kc/beta*(r-x)]
 *         = -(1/tau + K*Kc/(beta*tau))*x + (K*Kc*(beta-1)/(beta*tau))*z
 *           + (K*Kc/(beta*tau))*r
 *
 *   dz/dt = -1/(beta*T)*z + (r - x)
 *         = -x - 1/(beta*T)*z + r
 *
 * Complexity: O(n_points) for simulation, O(1) per step
 *
 * @param lag         lag compensator parameters
 * @param K_plant     plant DC gain
 * @param tau_plant   plant time constant
 * @param t_end       simulation end time
 * @param n_points    number of time points
 * @param[out] resp   step response data (caller frees time/output/error/control_signal)
 * @return 0 on success
 */
int lag_simulate_step_first_order(const LagCompensator *lag,
                                   double K_plant, double tau_plant,
                                   double t_end, int n_points,
                                   LagStepResponse *resp) {
    if (!lag || !resp || n_points < 2 || tau_plant <= 0) return -1;

    /* Allocation */
    resp->num_points = n_points;
    resp->time = (double*)malloc((size_t)n_points * sizeof(double));
    resp->output = (double*)malloc((size_t)n_points * sizeof(double));
    resp->error = (double*)malloc((size_t)n_points * sizeof(double));
    resp->control_signal = (double*)malloc((size_t)n_points * sizeof(double));

    if (!resp->time || !resp->output || !resp->error || !resp->control_signal) {
        free(resp->time); free(resp->output);
        free(resp->error); free(resp->control_signal);
        return -2;
    }

    /* Time step */
    double dt = t_end / (n_points - 1);

    /* Precompute coefficients */
    double a11 = -(1.0/tau_plant + K_plant * lag->Kc / (lag->beta * tau_plant));
    double a12 = K_plant * lag->Kc * (lag->beta - 1.0) / (lag->beta * tau_plant);
    double b1 = K_plant * lag->Kc / (lag->beta * tau_plant);

    double a21 = -1.0;
    double a22 = -1.0 / (lag->beta * lag->T);
    double b2 = 1.0;

    /* Initial conditions: zero */
    double x = 0.0;  /* plant state */
    double z = 0.0;  /* compensator state */
    double r = 1.0;  /* step input magnitude */

    /* Performance tracking */
    double y_ss_est = K_plant * lag->Kc / (1.0 + K_plant * lag->Kc);
    resp->final_value = y_ss_est;
    resp->rise_time = -1.0;
    resp->settling_time = -1.0;
    resp->overshoot_pct = 0.0;
    resp->peak_time = -1.0;
    double y_max = 0.0;
    int rise_done = 0, settled = 0;

    /* Euler integration */
    for (int i = 0; i < n_points; i++) {
        double t = i * dt;
        resp->time[i] = t;

        double y = x;  /* output = state for first-order plant */
        double e = r - y;
        double u = lag->Kc * (lag->beta - 1.0) / lag->beta * z +
                   lag->Kc / lag->beta * e;
        (void)u;  /* control signal computed for completeness; not stored in ramp response */

        resp->output[i] = y;
        resp->error[i] = e;
        resp->control_signal[i] = u;

        /* Track peak */
        if (y > y_max) {
            y_max = y;
            resp->peak_time = t;
        }

        /* Rise time: 10% to 90% */
        if (!rise_done && y >= 0.9 * y_ss_est) {
            /* Find 10% crossing earlier */
            for (int j = i - 1; j >= 0; j--) {
                if (resp->output[j] <= 0.1 * y_ss_est) {
                    double t10 = resp->time[j] +
                        (0.1*y_ss_est - resp->output[j]) /
                        (resp->output[j+1] - resp->output[j]) * dt;
                    double t90 = t - (y - 0.9*y_ss_est) /
                                  (y - resp->output[i-1]) * dt;
                    resp->rise_time = t90 - t10;
                    rise_done = 1;
                    break;
                }
            }
        }

        /* Settling time: within 2% of final value */
        if (!settled && fabs(y - y_ss_est) <= 0.02 * fabs(y_ss_est)) {
            resp->settling_time = t;
        }
        if (settled == 0 && fabs(y - y_ss_est) > 0.02 * fabs(y_ss_est)) {
            resp->settling_time = -1.0;  /* reset if exits band */
        }

        /* Euler update */
        double dx = a11 * x + a12 * z + b1 * r;
        double dz = a21 * x + a22 * z + b2 * r;

        x += dx * dt;
        z += dz * dt;
    }

    /* Overshoot */
    if (y_ss_est > 1e-10 && y_max > y_ss_est) {
        resp->overshoot_pct = (y_max - y_ss_est) / y_ss_est * 100.0;
    }

    /* Steady-state error */
    resp->steady_state_error = r - y_ss_est;

    return 0;
}

/**
 * Simulate ramp response of unity-feedback system.
 *
 * Input: r(t) = t (unit ramp starting at t=0)
 *
 * The steady-state velocity error for a type-0 system is infinite.
 * For a type-1 system, e_ss = 1/Kv.
 *
 * The lag compensator increases Kv by factor Kc, reducing e_ss by Kc.
 *
 * Complexity: O(n_points)
 */
int lag_simulate_ramp_first_order(const LagCompensator *lag,
                                   double K_plant, double tau_plant,
                                   double t_end, int n_points,
                                   LagRampResponse *resp) {
    if (!lag || !resp || n_points < 2 || tau_plant <= 0) return -1;

    resp->num_points = n_points;
    resp->time = (double*)malloc((size_t)n_points * sizeof(double));
    resp->output = (double*)malloc((size_t)n_points * sizeof(double));
    resp->error = (double*)malloc((size_t)n_points * sizeof(double));

    if (!resp->time || !resp->output || !resp->error) {
        free(resp->time); free(resp->output); free(resp->error);
        return -2;
    }

    double dt = t_end / (n_points - 1);

    /* Same coefficients as step simulation */
    double a11 = -(1.0/tau_plant + K_plant * lag->Kc / (lag->beta * tau_plant));
    double a12 = K_plant * lag->Kc * (lag->beta - 1.0) / (lag->beta * tau_plant);
    double b1_val = K_plant * lag->Kc / (lag->beta * tau_plant);

    double a21 = -1.0;
    double a22 = -1.0 / (lag->beta * lag->T);

    double x = 0.0, z = 0.0;

    /* For type-0 system, velocity error -> infinity.
     * Lag compensator cannot make a type-0 track a ramp with finite error.
     * But it can track with constant lag: e(t) -> tau/(K*Kc) as t -> inf.
     *
     * For a first-order plant: e_ss_ramp = 1/Kv = tau/(K*Kc) */
    resp->velocity_error = tau_plant / (K_plant * lag->Kc);

    for (int i = 0; i < n_points; i++) {
        double t = i * dt;
        double r = t;  /* ramp input */
        resp->time[i] = t;

        double y = x;
        double e = r - y;
        double u = lag->Kc * (lag->beta - 1.0) / lag->beta * z +
                   lag->Kc / lag->beta * e;
        (void)u;  /* control signal computed for completeness; not stored in ramp response */

        resp->output[i] = y;
        resp->error[i] = e;

        /* Euler update with time-varying r */
        double dx = a11 * x + a12 * z + b1_val * r;
        double dz = a21 * x + a22 * z + r;

        x += dx * dt;
        z += dz * dt;
    }

    return 0;
}

/**
 * Simulate step response using RK4 for higher accuracy.
 *
 * 4th-order Runge-Kutta method for the 2-state system.
 * Complexity: O(n_points), 4x evaluations per step vs Euler.
 */
int lag_simulate_step_rk4(const LagCompensator *lag,
                           double K_plant, double tau_plant,
                           double t_end, int n_points,
                           LagStepResponse *resp) {
    if (!lag || !resp || n_points < 2 || tau_plant <= 0) return -1;

    resp->num_points = n_points;
    resp->time = (double*)malloc((size_t)n_points * sizeof(double));
    resp->output = (double*)malloc((size_t)n_points * sizeof(double));
    resp->error = (double*)malloc((size_t)n_points * sizeof(double));
    resp->control_signal = (double*)malloc((size_t)n_points * sizeof(double));

    if (!resp->time || !resp->output || !resp->error || !resp->control_signal) {
        free(resp->time); free(resp->output);
        free(resp->error); free(resp->control_signal);
        return -2;
    }

    double dt = t_end / (n_points - 1);
    double r = 1.0;  /* step input */

    /* System coefficients — same as Euler */
    double b1_val = K_plant * lag->Kc / (lag->beta * tau_plant);
    double a22 = -1.0 / (lag->beta * lag->T);

    double x = 0.0, z = 0.0;

    /* Derivative function for 2-state system:
     * dx/dt = f1(x,z) = -(1/tau + K*Kc/(beta*tau))*x
     *                  + K*Kc*(beta-1)/(beta*tau)*z
     *                  + K*Kc/(beta*tau)*r
     * dz/dt = f2(x,z) = -x - 1/(beta*T)*z + r
     */
    double a11_val = -(1.0/tau_plant + K_plant * lag->Kc / (lag->beta * tau_plant));
    double a12_val = K_plant * lag->Kc * (lag->beta - 1.0) / (lag->beta * tau_plant);

    double y_ss_est = K_plant * lag->Kc / (1.0 + K_plant * lag->Kc);
    resp->final_value = y_ss_est;
    double y_max = 0.0;
    resp->peak_time = -1.0;
    resp->rise_time = -1.0;
    resp->settling_time = -1.0;

    for (int i = 0; i < n_points; i++) {
        double t = i * dt;
        resp->time[i] = t;

        double y = x;
        double e = r - y;
        double u = lag->Kc * (lag->beta - 1.0) / lag->beta * z +
                   lag->Kc / lag->beta * e;
        (void)u;  /* control signal computed for completeness; not stored in ramp response */

        resp->output[i] = y;
        resp->error[i] = e;
        resp->control_signal[i] = u;

        /* Track max */
        if (y > y_max) { y_max = y; resp->peak_time = t; }

        /* Settling time */
        if (resp->settling_time < 0 && fabs(y - y_ss_est) <= 0.02 * fabs(y_ss_est))
            resp->settling_time = t;

        /* RK4 step */
        double k1x = (a11_val * x + a12_val * z + b1_val * r);
        double k1z = (-x + a22 * z + r);
        double k2x = (a11_val * (x + 0.5*dt*k1x) + a12_val * (z + 0.5*dt*k1z) + b1_val * r);
        double k2z = (-(x + 0.5*dt*k1x) + a22 * (z + 0.5*dt*k1z) + r);
        double k3x = (a11_val * (x + 0.5*dt*k2x) + a12_val * (z + 0.5*dt*k2z) + b1_val * r);
        double k3z = (-(x + 0.5*dt*k2x) + a22 * (z + 0.5*dt*k2z) + r);
        double k4x = (a11_val * (x + dt*k3x) + a12_val * (z + dt*k3z) + b1_val * r);
        double k4z = (-(x + dt*k3x) + a22 * (z + dt*k3z) + r);

        x += (dt/6.0) * (k1x + 2*k2x + 2*k3x + k4x);
        z += (dt/6.0) * (k1z + 2*k2z + 2*k3z + k4z);
    }

    /* Overshoot */
    if (y_ss_est > 1e-10 && y_max > y_ss_est) {
        resp->overshoot_pct = (y_max - y_ss_est) / y_ss_est * 100.0;
    }

    /* Rise time — find 10% and 90% crossings */
    double t10 = -1, t90 = -1;
    for (int i = 0; i < n_points; i++) {
        if (t10 < 0 && resp->output[i] >= 0.1 * y_ss_est) {
            if (i > 0) t10 = resp->time[i-1] +
                (0.1*y_ss_est - resp->output[i-1]) /
                (resp->output[i] - resp->output[i-1]) * dt;
            else t10 = resp->time[i];
        }
        if (t90 < 0 && resp->output[i] >= 0.9 * y_ss_est) {
            if (i > 0) t90 = resp->time[i-1] +
                (0.9*y_ss_est - resp->output[i-1]) /
                (resp->output[i] - resp->output[i-1]) * dt;
            else t90 = resp->time[i];
            break;
        }
    }
    if (t10 >= 0 && t90 >= 0) resp->rise_time = t90 - t10;

    resp->steady_state_error = r - y_ss_est;

    return 0;
}

/**
 * Simulate step response for a general second-order plant
 * with lag compensator.
 *
 * G(s) = K * w_n^2 / (s^2 + 2*zeta*w_n*s + w_n^2)
 *
 * State-space (controllable canonical form):
 *   [dx1]   [    0         1    ] [x1]   [0]
 *   [dx2] = [-w_n^2   -2*zeta*w_n] [x2] + [K*w_n^2] u
 *   y = [1  0] [x1; x2]
 *
 * Combined with compensator (3rd state z):
 *   Total system is 3rd order.
 *
 * Complexity: O(n_points)
 */
int lag_simulate_step_second_order(const LagCompensator *lag,
                                    double K_plant, double zeta, double wn,
                                    double t_end, int n_points,
                                    LagStepResponse *resp) {
    if (!lag || !resp || n_points < 2 || wn <= 0) return -1;

    resp->num_points = n_points;
    resp->time = (double*)malloc((size_t)n_points * sizeof(double));
    resp->output = (double*)malloc((size_t)n_points * sizeof(double));
    resp->error = (double*)malloc((size_t)n_points * sizeof(double));
    resp->control_signal = (double*)malloc((size_t)n_points * sizeof(double));

    if (!resp->time || !resp->output || !resp->error || !resp->control_signal) {
        free(resp->time); free(resp->output);
        free(resp->error); free(resp->control_signal);
        return -2;
    }

    double dt = t_end / (n_points - 1);
    double r = 1.0;

    /* Plant state-space matrices */
    double A_plant[4] = {0.0, 1.0, -wn*wn, -2.0*zeta*wn};
    double B_plant[2] = {0.0, K_plant * wn * wn};
    /* double C_plant[2] = {1.0, 0.0}; */

    /* Compensator */
    double Kc = lag->Kc;
    double beta = lag->beta;
    double T = lag->T;
    double a_c = -1.0 / (beta * T);
    double b_c = Kc * (beta - 1.0) / beta;
    double d_c = Kc / beta;

    /* State vector: [x1, x2, z] */
    double x1 = 0.0, x2 = 0.0, xz = 0.0;

    /* Estimate final value */
    double dc_gain = K_plant * Kc;
    double y_ss = dc_gain / (1.0 + dc_gain);
    resp->final_value = y_ss;

    double y_max = 0.0;
    resp->peak_time = -1.0;
    resp->settling_time = -1.0;
    resp->rise_time = -1.0;

    for (int i = 0; i < n_points; i++) {
        double t = i * dt;
        resp->time[i] = t;

        double y = x1;
        double e = r - y;
        double u = b_c * xz + d_c * e;

        resp->output[i] = y;
        resp->error[i] = e;
        resp->control_signal[i] = u;

        if (y > y_max) { y_max = y; resp->peak_time = t; }
        if (resp->settling_time < 0 && fabs(y - y_ss) <= 0.02 * fabs(y_ss))
            resp->settling_time = t;

        /* Euler integration */
        double dx1 = A_plant[0]*x1 + A_plant[1]*x2 + B_plant[0]*u;
        double dx2 = A_plant[2]*x1 + A_plant[3]*x2 + B_plant[1]*u;
        double dxz = a_c * xz + e;

        x1 += dx1 * dt;
        x2 += dx2 * dt;
        xz += dxz * dt;
    }

    if (y_ss > 1e-10 && y_max > y_ss)
        resp->overshoot_pct = (y_max - y_ss) / y_ss * 100.0;

    resp->steady_state_error = r - y_ss;

    return 0;
}

/**
 * Free memory allocated for a step response structure.
 * Complexity: O(1)
 */
void lag_step_response_free(LagStepResponse *resp) {
    if (resp) {
        free(resp->time);
        free(resp->output);
        free(resp->error);
        free(resp->control_signal);
        resp->time = NULL;
        resp->output = NULL;
        resp->error = NULL;
        resp->control_signal = NULL;
        resp->num_points = 0;
    }
}

/**
 * Free memory allocated for a ramp response structure.
 * Complexity: O(1)
 */
void lag_ramp_response_free(LagRampResponse *resp) {
    if (resp) {
        free(resp->time);
        free(resp->output);
        free(resp->error);
        resp->time = NULL;
        resp->output = NULL;
        resp->error = NULL;
        resp->num_points = 0;
    }
}