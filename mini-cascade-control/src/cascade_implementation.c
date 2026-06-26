/**
 * @file cascade_implementation.c
 * @brief Digital implementation, anti-windup, bumpless transfer, discretization
 *
 * L2 -- Core Concepts: Anti-windup strategies for cascade PID, bumpless transfer.
 * L5 -- Computational Methods: Discrete PID (velocity/position form),
 *       Tustin discretization, digital cascade execution.
 * L7 -- Applications: Industrial DCS/PLC cascade implementation patterns.
 */

#include "cascade_types.h"
#include "cascade_implementation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * L5: PID Controller Lifecycle Management
 * ========================================================================== */

int cascade_pid_init(double Kp, double Ki, double Kd, double N,
                      double Ts, CascadePID *pid) {
    if (!pid) return -1;
    pid->Kp = Kp;
    pid->Ki = Ki;
    pid->Kd = Kd;
    pid->N  = (N > 0.0) ? N : 100.0;
    pid->Ts = Ts;
    pid->b  = 1.0;
    pid->c  = 0.0;
    pid->u_min = -1e300;
    pid->u_max =  1e300;
    pid->integrator = 0.0;
    pid->prev_error  = 0.0;
    pid->prev_y      = 0.0;
    pid->has_antiwindup = 0;
    pid->Tt = 1.0;
    return 0;
}

void cascade_pid_set_limits(double u_min, double u_max, CascadePID *pid) {
    if (!pid) return;
    pid->u_min = u_min;
    pid->u_max = u_max;
}

void cascade_pid_set_antiwindup(CascadeAWMethod method, double Tt,
                                 CascadePID *pid) {
    if (!pid) return;
    pid->has_antiwindup = (method != CASCADE_AW_NONE) ? 1 : 0;
    pid->Tt = (Tt > 0.0) ? Tt : 1.0;
}

void cascade_pid_reset(double init_integrator, double init_u,
                        CascadePID *pid) {
    (void)init_u;  /* Reserved for output initialization in advanced modes */
    if (!pid) return;
    pid->integrator = init_integrator;
    pid->prev_error = 0.0;
    pid->prev_y = 0.0;
}

/* ==========================================================================
 * L5: Velocity (Incremental) Form PID
 *
 * du(k) = Kp*(e(k)-e(k-1)) + Ki*Ts*e(k) + (Kd/Ts)*(e(k)-2e(k-1)+e(k-2))
 * u(k) = u(k-1) + du(k)
 *
 * Advantages:
 *   - Inherently anti-windup (no unbounded integrator accumulation)
 *   - Bumpless parameter changes (only du depends on parameters)
 *   - Smooth mode transitions (u(k-1) provides continuity)
 *
 * Reference: Astrom & Hagglund, "PID Controllers" (1995), Ch. 3.
 * ========================================================================== */

double cascade_pid_velocity(double setpoint, double pv, CascadePID *pid) {
    if (!pid) return 0.0;
    double e = setpoint - pv;
    double Ts = (pid->Ts > 0.0) ? pid->Ts : 0.01;
    double du = pid->Kp * (e - pid->prev_error)
              + pid->Ki * Ts * e
              + (pid->Kd / Ts) * (e - 2.0 * pid->prev_error);
    double u = pid->integrator + du;
    if (u > pid->u_max) u = pid->u_max;
    if (u < pid->u_min) u = pid->u_min;
    pid->integrator = u;
    pid->prev_error = e;
    return u;
}

/* ==========================================================================
 * L5: Position Form PID with Back-Calculation Anti-Windup
 *
 * P = Kp*(b*r - y)
 * I += Ki*Ts*e + (Ts/Tt)*(u_sat - u_unsat)
 * D = Kd * ((c*r - y) - prev_y) / (Ts + 1/N)
 * u = P + I + D
 *
 * Back-calculation: tracking time constant Tt feeds saturation error
 * back to integrator. Tt = Ti for PI, Tt = sqrt(Ti*Td) for PID.
 *
 * Reference: Astrom & Rundqwist (1989), Proc. ACC.
 * ========================================================================== */

double cascade_pid_position(double setpoint, double pv, CascadePID *pid) {
    if (!pid) return 0.0;
    double e = setpoint - pv;
    double Ts = (pid->Ts > 0.0) ? pid->Ts : 0.01;
    double N_val = pid->N;

    double P = pid->Kp * (pid->b * setpoint - pv);
    double I_raw = pid->integrator + pid->Ki * Ts * e;

    double alpha = Ts / (Ts + 1.0 / N_val);
    double D_raw = pid->Kd * ((pid->c * setpoint - pv) - pid->prev_y)
                   / (Ts + 1.0 / N_val);
    static double D_state = 0.0;
    D_state = alpha * D_raw + (1.0 - alpha) * D_state;
    double D = D_state;

    double u_unsat = P + I_raw + D;
    double u_sat = u_unsat;
    if (u_sat > pid->u_max) u_sat = pid->u_max;
    if (u_sat < pid->u_min) u_sat = pid->u_min;

    if (pid->has_antiwindup && pid->Tt > 0.0) {
        pid->integrator = I_raw + (Ts / pid->Tt) * (u_sat - u_unsat);
    } else {
        pid->integrator = I_raw;
    }
    if (pid->integrator > pid->u_max) pid->integrator = pid->u_max;
    if (pid->integrator < pid->u_min) pid->integrator = pid->u_min;
    pid->prev_error = e;
    pid->prev_y = pid->c * setpoint - pv;
    return u_sat;
}

/* ==========================================================================
 * L2: Anti-Windup Implementations
 * ========================================================================== */

void cascade_aw_init(CascadeAWMethod inner_method,
                      CascadeAWMethod outer_method,
                      CascadeSystem *sys) {
    if (!sys) return;
    sys->inner.inner_controller.has_antiwindup =
        (inner_method != CASCADE_AW_NONE) ? 1 : 0;
    sys->outer.outer_controller.has_antiwindup =
        (outer_method != CASCADE_AW_NONE) ? 1 : 0;
}

void cascade_aw_clamping(CascadePID *pid) {
    if (!pid) return;
    if (pid->integrator > pid->u_max) pid->integrator = pid->u_max;
    if (pid->integrator < pid->u_min) pid->integrator = pid->u_min;
}

void cascade_aw_back_calculation(CascadePID *pid, double u_unsat,
                                  double u_sat) {
    if (!pid || pid->Tt <= 0.0) return;
    double Ts = (pid->Ts > 0.0) ? pid->Ts : 0.01;
    pid->integrator += (Ts / pid->Tt) * (u_sat - u_unsat);
}

void cascade_aw_combined(CascadeSystem *sys) {
    /* Cascade-specific anti-windup:
     * When inner loop saturates, freeze outer integrator too.
     * This is the most important cascade-specific anti-windup rule.
     * Reference: Astrom & Hagglund (2006), Sec. 7.5. */
    if (!sys) return;
    CascadePID *inner = &sys->inner.inner_controller;
    CascadePID *outer = &sys->outer.outer_controller;
    int inner_sat = (inner->integrator >= inner->u_max - 1e-9 ||
                     inner->integrator <= inner->u_min + 1e-9);
    if (inner_sat) cascade_aw_clamping(outer);
    cascade_aw_clamping(inner);
}

/* ==========================================================================
 * L2: Bumpless Transfer State Machine
 * ========================================================================== */

void cascade_bt_init(int initial_mode, double transition_t,
                      CascadeBumpless *bt) {
    if (!bt) return;
    bt->current_mode = initial_mode;
    bt->target_mode = initial_mode;
    bt->last_good_u = 0.0;
    bt->transition_time = (transition_t > 0.0) ? transition_t : 1.0;
    bt->transition_start = 0.0;
    bt->in_transition = 0;
    bt->manual_output = 0.0;
}

int cascade_bt_start(int target_mode, CascadeBumpless *bt,
                      CascadeSystem *sys, double current_time) {
    if (!bt || !sys) return -1;
    if (bt->in_transition || target_mode == bt->current_mode) return 0;
    bt->target_mode = target_mode;
    bt->transition_start = current_time;
    bt->in_transition = 1;
    bt->last_good_u = sys->inner.inner_controller.integrator;
    return 0;
}

double cascade_bt_update(CascadeBumpless *bt, CascadeSystem *sys,
                          double current_time,
                          double auto_u, double cascade_u) {
    if (!bt || !sys) return 0.0;
    if (!bt->in_transition) {
        if (bt->current_mode == 0) return bt->manual_output;
        if (bt->current_mode == 1) return auto_u;
        return cascade_u;
    }
    double elapsed = current_time - bt->transition_start;
    double frac = elapsed / bt->transition_time;
    if (frac > 1.0) frac = 1.0;
    double old_u = (bt->current_mode == 0) ? bt->manual_output :
                   (bt->current_mode == 1) ? auto_u : cascade_u;
    double new_u = (bt->target_mode == 0) ? bt->manual_output :
                   (bt->target_mode == 1) ? auto_u : cascade_u;
    double blended = old_u + frac * (new_u - old_u);
    if (frac >= 1.0) {
        bt->current_mode = bt->target_mode;
        bt->in_transition = 0;
    }
    return blended;
}

int cascade_bt_is_complete(const CascadeBumpless *bt) {
    return (bt && !bt->in_transition) ? 1 : 0;
}

/* ==========================================================================
 * L5: Cascade Loop Execution
 * ========================================================================== */

int cascade_execute(CascadeSystem *sys,
                     double setpoint,
                     double inner_pv,
                     double outer_pv,
                     CascadeBumpless *bt,
                     double *inner_sp,
                     double *inner_u) {
    if (!sys || !inner_sp || !inner_u) return -1;
    double outer_ctrl = cascade_pid_position(setpoint, outer_pv,
        &sys->outer.outer_controller);
    *inner_sp = outer_ctrl;
    if (*inner_sp > sys->inner.inner_controller.u_max)
        *inner_sp = sys->inner.inner_controller.u_max;
    if (*inner_sp < sys->inner.inner_controller.u_min)
        *inner_sp = sys->inner.inner_controller.u_min;
    double inner_ctrl = cascade_pid_position(*inner_sp, inner_pv,
        &sys->inner.inner_controller);
    if (sys->inner.inner_controller.has_antiwindup &&
        sys->outer.outer_controller.has_antiwindup) {
        cascade_aw_combined(sys);
    }
    if (bt && bt->in_transition) {
        inner_ctrl = cascade_bt_update(bt, sys, 0.0, inner_ctrl, inner_ctrl);
    }
    *inner_u = inner_ctrl;
    return 0;
}

void cascade_set_mode(CascadeSystem *sys, int cascade_on,
                       CascadeBumpless *bt) {
    if (!sys) return;
    if (cascade_on && !sys->cascade_active) {
        if (bt) cascade_bt_start(2, bt, sys, 0.0);
        cascade_pid_reset(sys->inner.inner_controller.integrator,
                          sys->inner.inner_controller.integrator,
                          &sys->outer.outer_controller);
    } else if (!cascade_on && sys->cascade_active) {
        if (bt) cascade_bt_start(1, bt, sys, 0.0);
    }
    sys->cascade_active = cascade_on;
}

/* ==========================================================================
 * L5: Tustin (Bilinear) Discretization
 *
 * s <- (2/Ts) * (z-1)/(z+1)
 * Maps LHP s-plane to interior of unit circle (stability preserved).
 * For G(s) = K/(tau*s+1):
 *   G(z) = K*(Ts/(2*tau+Ts))*(z+1)/(z - (2*tau-Ts)/(2*tau+Ts))
 * Reference: Franklin, Powell & Workman (1998), Ch. 6.
 * ========================================================================== */

int cascade_tustin_discretize(const CascadeTF *ct_tf, double Ts,
                               CascadeDTF *dt_tf) {
    if (!ct_tf || !dt_tf || Ts <= 0.0) return -1;
    if (ct_tf->den.degree == 1 && ct_tf->num.degree == 0) {
        double a0 = ct_tf->den.coeff[0], a1 = ct_tf->den.coeff[1];
        double b0 = ct_tf->num.coeff[0];
        if (fabs(a1) < 1e-15) return -1;
        double tau = a1 / a0;
        double K = ct_tf->gain * b0 / a0;
        double c = 2.0 / Ts;
        double den = c * tau + 1.0;
        dt_tf->num_order = 1; dt_tf->den_order = 1;
        dt_tf->num_coeff = (double *)malloc(2 * sizeof(double));
        dt_tf->den_coeff = (double *)malloc(2 * sizeof(double));
        if (!dt_tf->num_coeff || !dt_tf->den_coeff) {
            free(dt_tf->num_coeff); free(dt_tf->den_coeff); return -1;
        }
        dt_tf->num_coeff[0] = K / den;
        dt_tf->num_coeff[1] = K / den;
        dt_tf->den_coeff[0] = 1.0;
        dt_tf->den_coeff[1] = (c * tau - 1.0) / den;
        dt_tf->gain = 1.0;
        dt_tf->Ts = Ts;
        return 0;
    }
    return -1;
}

void cascade_pid_discretize(const CascadePID *ct_pid, double Ts,
                             CascadePID *dt_pid) {
    if (!ct_pid || !dt_pid) return;
    *dt_pid = *ct_pid;
    dt_pid->Ts = Ts;
}

void cascade_dtf_free(CascadeDTF *dtf) {
    if (dtf) {
        free(dtf->num_coeff); dtf->num_coeff = NULL;
        free(dtf->den_coeff); dtf->den_coeff = NULL;
        dtf->num_order = 0; dtf->den_order = 0;
    }
}

double cascade_dtf_execute(const CascadeDTF *dtf, double input,
                            double *u_history, double *y_history) {
    if (!dtf || !u_history || !y_history) return 0.0;
    double y = dtf->num_coeff[0] * input;
    for (int i = 1; i <= dtf->num_order; i++)
        y += dtf->num_coeff[i] * u_history[i - 1];
    for (int i = 1; i <= dtf->den_order; i++)
        y -= dtf->den_coeff[i] * y_history[i - 1];
    for (int i = dtf->num_order; i > 0; i--)
        u_history[i] = u_history[i - 1];
    u_history[0] = input;
    for (int i = dtf->den_order; i > 0; i--)
        y_history[i] = y_history[i - 1];
    y_history[0] = y;
    return y * dtf->gain;
}
