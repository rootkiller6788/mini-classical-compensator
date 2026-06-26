/**
 * @file feedforward_adaptive.c
 * @brief Adaptive and learning feedforward: LMS, ILC, nonlinear FF.
 *
 * L5-L8: LMS adaptive filtering for feedforward, iterative learning
 * control (ILC), computed-torque nonlinear feedforward.
 *
 * L7 Applications: DC motor servo, 2-link robot arm, pendulum.
 */

#include "feedforward_adaptive.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L5: LMS Adaptive Filter
 * ========================================================================== */

void lms_init(LMSFilter *filter, int n_taps, double mu, double leakage)
{
    if (!filter || n_taps <= 0) return;

    filter->n_taps = n_taps;
    filter->mu = mu;
    filter->leakage = leakage;
    filter->norm_factor = 1.0;

    filter->weights = (double *)calloc(n_taps, sizeof(double));
    filter->buffer = (double *)calloc(n_taps, sizeof(double));
    filter->buf_idx = 0;

    /* Initialize first weight to small value to break symmetry */
    if (filter->weights) filter->weights[0] = 0.001;
}

double lms_process(LMSFilter *filter, double x, double desired)
{
    if (!filter || !filter->weights || !filter->buffer) return 0.0;

    /* Insert input into circular buffer */
    filter->buffer[filter->buf_idx] = x;

    /* Compute filter output */
    double y = 0.0;
    for (int k = 0; k < filter->n_taps; k++) {
        int idx = (filter->buf_idx - k + filter->n_taps) % filter->n_taps;
        y += filter->weights[k] * filter->buffer[idx];
    }

    /* LMS weight update: w[n+1] = (1 - mu*leakage)*w[n] + mu*e[n]*x[n] */
    double error = desired - y;
    for (int k = 0; k < filter->n_taps; k++) {
        int idx = (filter->buf_idx - k + filter->n_taps) % filter->n_taps;
        filter->weights[k] = (1.0 - filter->mu * filter->leakage) * filter->weights[k]
                             + filter->mu * error * filter->buffer[idx];
    }

    /* Advance buffer pointer */
    filter->buf_idx = (filter->buf_idx + 1) % filter->n_taps;

    return y;
}

double lms_process_normalized(LMSFilter *filter, double x,
                               double desired, double epsilon)
{
    if (!filter || !filter->weights || !filter->buffer) return 0.0;

    /* Insert input */
    filter->buffer[filter->buf_idx] = x;

    /* Compute output */
    double y = 0.0;
    for (int k = 0; k < filter->n_taps; k++) {
        int idx = (filter->buf_idx - k + filter->n_taps) % filter->n_taps;
        y += filter->weights[k] * filter->buffer[idx];
    }

    /* Compute input power for normalization */
    double power = epsilon;
    for (int k = 0; k < filter->n_taps; k++) {
        int idx = (filter->buf_idx - k + filter->n_taps) % filter->n_taps;
        power += filter->buffer[idx] * filter->buffer[idx];
    }

    /* Normalized LMS: mu_norm = mu / (epsilon + ||x||^2) */
    double mu_norm = filter->mu / power;

    /* Weight update */
    double error = desired - y;
    for (int k = 0; k < filter->n_taps; k++) {
        int idx = (filter->buf_idx - k + filter->n_taps) % filter->n_taps;
        filter->weights[k] = (1.0 - mu_norm * filter->leakage) * filter->weights[k]
                             + mu_norm * error * filter->buffer[idx];
    }

    filter->buf_idx = (filter->buf_idx + 1) % filter->n_taps;
    return y;
}

void lms_get_weights(const LMSFilter *filter, double *weights_out)
{
    if (filter && filter->weights && weights_out) {
        memcpy(weights_out, filter->weights, filter->n_taps * sizeof(double));
    }
}

void lms_free(LMSFilter *filter)
{
    if (filter) {
        free(filter->weights);
        free(filter->buffer);
        filter->weights = NULL;
        filter->buffer = NULL;
        filter->n_taps = 0;
    }
}

/* ==========================================================================
 * L8: Iterative Learning Control
 * ========================================================================== */

void ilc_init(ILCController *ilc, int n_samples, double dt,
              const TransferFn *learn_l, const TransferFn *learn_q)
{
    if (!ilc || n_samples <= 0) return;

    ilc->n_samples = n_samples;
    ilc->dt = dt;
    ilc->trial_num = 0;
    ilc->convergence_rate = 1.0;

    if (learn_l) ilc->learning_filter = *learn_l;
    if (learn_q) ilc->q_filter = *learn_q;

    ilc->trial_data = (double *)calloc(n_samples, sizeof(double));
    ilc->error_data = (double *)calloc(n_samples, sizeof(double));
}

void ilc_update(ILCController *ilc,
                const double *u_current,
                const double *error,
                double *u_next)
{
    if (!ilc || !u_current || !error || !u_next) return;

    /* ILC update law: u_{k+1} = Q * (u_k + L * e_k)
     *
     * In discrete time, Q and L are filters applied to the signals.
     * For simplicity, implement as:
     *   u_next[i] = u_current[i] + learning_gain * error[i]
     * with a first-order Q-filter applied. */

    /* Extract learning gain from L filter (use DC gain as scalar gain) */
    double L_gain = ff_dc_gain(&ilc->learning_filter, 1);

    /* Q-filter smoothing: first-order low-pass on error */
    double q_state = 0.0;
    double q_alpha = ilc->dt / (ilc->dt + 1.0 / (2.0 * M_PI * 10.0)); /* ~10 Hz */

    for (int i = 0; i < ilc->n_samples; i++) {
        /* Apply learning gain */
        double correction = L_gain * error[i];

        /* Apply Q-filter (low-pass smoothing) */
        q_state = q_alpha * correction + (1.0 - q_alpha) * q_state;
        double filtered_correction = q_state;

        u_next[i] = u_current[i] + filtered_correction;
    }

    ilc->trial_num++;
}

void ilc_convergence_check(const ILCController *ilc,
                           const TransferFn *plant,
                           const double *omega, int n_freq,
                           double *margin)
{
    if (!ilc || !plant || !omega || !margin) return;

    /* Convergence condition: |Q(jw) * (1 - L(jw)*P(jw))| < 1
     * Compute max value across frequencies. */
    double max_val = 0.0;

    for (int i = 0; i < n_freq; i++) {
        double mag_q, ph_q, mag_l, ph_l, mag_p, ph_p;
        tf_freq_response(&ilc->q_filter, omega[i], &mag_q, &ph_q);
        tf_freq_response(&ilc->learning_filter, omega[i], &mag_l, &ph_l);
        tf_freq_response(plant, omega[i], &mag_p, &ph_p);

        /* 1 - L*P */
        double mag_lp = mag_l * mag_p;
        double ph_lp = ph_l + ph_p;
        double re_lp = mag_lp * cos(ph_lp);
        double im_lp = mag_lp * sin(ph_lp);
        double mag_one_minus_lp = sqrt((1.0 - re_lp) * (1.0 - re_lp) + im_lp * im_lp);

        double val = mag_q * mag_one_minus_lp;
        if (val > max_val) max_val = val;
    }

    *margin = max_val;
}

void ilc_free(ILCController *ilc)
{
    if (ilc) {
        free(ilc->trial_data);
        free(ilc->error_data);
        ilc->trial_data = NULL;
        ilc->error_data = NULL;
    }
}

/* ==========================================================================
 * L8: Nonlinear Feedforward (Computed Torque)
 * ========================================================================== */

void nl_ff_computed_torque(const NonlinearFFModel *model,
                           const double *q, const double *qd,
                           const double *qdd, double *tau_ff)
{
    if (!model || !q || !qd || !qdd || !tau_ff) return;

    int n = model->n_dof;

    /* M(q) matrix */
    double *M = (double *)malloc(n * n * sizeof(double));
    double *C = (double *)malloc(n * sizeof(double));
    double *G = (double *)malloc(n * sizeof(double));
    double *F = (double *)malloc(n * sizeof(double));

    if (!M || !C || !G || !F) {
        free(M); free(C); free(G); free(F);
        return;
    }

    /* Evaluate model functions */
    if (model->inertia) model->inertia(q, M, n, model->params);
    if (model->coriolis) model->coriolis(q, qd, C, n, model->params);
    if (model->gravity) model->gravity(q, G, n, model->params);
    if (model->friction) model->friction(qd, F, n, model->params);

    /* tau_ff = M(q)*qdd + C(q,qd) + G(q) + F(qd) */
    for (int i = 0; i < n; i++) {
        /* M*qdd term */
        double Mqdd = 0.0;
        for (int j = 0; j < n; j++) {
            Mqdd += M[i * n + j] * qdd[j];
        }
        tau_ff[i] = Mqdd + C[i] + G[i] + F[i];
    }

    free(M); free(C); free(G); free(F);
}

double nl_ff_pendulum(double J, double b, double mgl,
                       double q_d, double qd_d, double qdd_d)
{
    /* Dynamics: J*qdd + b*qd + mgl*sin(q) = tau
     * Feedforward: tau_ff = J*qdd_d + b*qd_d + mgl*sin(q_d) */
    return J * qdd_d + b * qd_d + mgl * sin(q_d);
}

void nl_ff_2link_arm(double m1, double m2, double a1, double a2,
                     const double q[2], const double qd[2],
                     const double qdd[2], double tau[2])
{
    /* Standard two-link manipulator dynamics.
     * References: Craig "Introduction to Robotics" (2005),
     *             Spong "Robot Modeling and Control" (2006).
     *
     * Simplified with point masses at link ends.
     *
     * M(q) = [M11 M12; M21 M22]
     * M11 = (m1+m2)*a1^2 + m2*a2^2 + 2*m2*a1*a2*cos(q2)
     * M12 = M21 = m2*a2^2 + m2*a1*a2*cos(q2)
     * M22 = m2*a2^2
     *
     * C(q,qd) = [-m2*a1*a2*sin(q2)*(2*qd1*qd2 + qd2^2);
     *             m2*a1*a2*sin(q2)*qd1^2]
     *
     * G(q) = [(m1+m2)*g*a1*cos(q1) + m2*g*a2*cos(q1+q2);
     *          m2*g*a2*cos(q1+q2)]
     */

    double g = 9.81; /* gravitational acceleration */
    double a1sq = a1 * a1;
    double a2sq = a2 * a2;
    double a1a2 = a1 * a2;
    double cos_q2 = cos(q[1]);
    double sin_q2 = sin(q[1]);
    double cos_q1 = cos(q[0]);
    double cos_q12 = cos(q[0] + q[1]);

    /* Inertia matrix */
    double M11 = (m1 + m2) * a1sq + m2 * a2sq + 2.0 * m2 * a1a2 * cos_q2;
    double M12 = m2 * a2sq + m2 * a1a2 * cos_q2;
    double M21 = M12;
    double M22 = m2 * a2sq;

    /* Coriolis/centripetal */
    double C1 = -m2 * a1a2 * sin_q2 * (2.0 * qd[0] * qd[1] + qd[1] * qd[1]);
    double C2 = m2 * a1a2 * sin_q2 * qd[0] * qd[0];

    /* Gravity */
    double G1 = (m1 + m2) * g * a1 * cos_q1 + m2 * g * a2 * cos_q12;
    double G2 = m2 * g * a2 * cos_q12;

    /* tau = M*qdd + C + G */
    tau[0] = M11 * qdd[0] + M12 * qdd[1] + C1 + G1;
    tau[1] = M21 * qdd[0] + M22 * qdd[1] + C2 + G2;
}

double ff_servo_motion(double Kv_ff, double Ka_ff, double Ks_ff,
                        double v_d, double a_d)
{
    /* Standard industrial motion feedforward:
     * u_ff = Kv_ff * v_d + Ka_ff * a_d + Ks_ff * sign(v_d)
     *
     * This is the typical structure in commercial servo drives
     * (e.g., Siemens SINAMICS, Bosch Rexroth, Yaskawa Sigma). */

    double u = Kv_ff * v_d + Ka_ff * a_d;
    if (v_d > 0.0) u += Ks_ff;
    else if (v_d < 0.0) u -= Ks_ff;
    return u;
}
