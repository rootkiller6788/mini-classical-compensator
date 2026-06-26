/**
 * pid_tuning_core.c — Core PID Controller Implementation
 *
 * Implements the foundational PID controller runtime:
 * initialization, parameter conversion, update equation,
 * performance evaluation, frequency response analysis.
 *
 * Knowledge: L1 (PID forms), L2 (feedback), L3 (transfer functions),
 *             L5 (PID update algorithm).
 *
 * Reference: Åström & Hägglund, "PID Controllers" (1995), Ch. 3.
 */

#include "pid_tuning.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L5: Initialize PID Controller
 * ────────────────────────────────────────────── */

void pid_init(pid_controller_t *pid, double Kp, double Ki, double Kd, double Ts)
{
    if (!pid) return;

    memset(pid, 0, sizeof(*pid));

    pid->params.Kp = Kp;
    pid->params.Ki = Ki;
    pid->params.Kd = Kd;
    pid->params.Ti = (Ki > 1e-12) ? (Kp / Ki) : 0.0;
    pid->params.Td = (Kp > 1e-12) ? (Kd / Kp) : 0.0;
    pid->params.N  = 10.0;       /* Default derivative filter factor */
    pid->params.b  = 1.0;        /* Default: full P on error */
    pid->params.c  = 1.0;        /* Default: full D on error */
    pid->params.Tt = 1.0;        /* Default tracking time constant */

    pid->form    = PID_FORM_PARALLEL;
    pid->action  = PID_DIRECT_ACTING;
    pid->aw_mode = AW_MODE_NONE;

    pid->I      = 0.0;
    pid->e_prev = 0.0;
    pid->y_prev = 0.0;
    pid->u_sat  = 0.0;
    pid->u_unsat= 0.0;
    pid->u_min  = -1.0;
    pid->u_max  =  1.0;
    pid->Ts     = Ts > 0 ? Ts : 0.01;
    pid->step_count = 0;
}

/* ──────────────────────────────────────────────
 * L3: Parameter Form Conversions
 * ────────────────────────────────────────────── */

void pid_ideal_to_parallel(double Kp, double Ti, double Td,
                           double *Ki, double *Kd_out)
{
    if (Ki) {
        if (Ti > 1e-12) {
            *Ki = Kp / Ti;
        } else {
            *Ki = 0.0;
        }
    }
    if (Kd_out) {
        *Kd_out = Kp * Td;
    }
}

void pid_parallel_to_ideal(double Kp, double Ki, double Kd,
                           double *Ti, double *Td)
{
    if (Ti) {
        if (Ki > 1e-12) {
            *Ti = Kp / Ki;
        } else {
            *Ti = 0.0;
        }
    }
    if (Td) {
        if (Kp > 1e-12) {
            *Td = Kd / Kp;
        } else {
            *Td = 0.0;
        }
    }
}

void pid_series_to_parallel(double Kp_s, double Ti_s, double Td_s,
                            double *Kp_p, double *Ki_p, double *Kd_p)
{
    /**
     * Series (interacting) PID:
     *   C(s) = Kp' * (1 + 1/(Ti'*s)) * (1 + Td'*s)
     *        = Kp' * (1 + Td'/Ti' + 1/(Ti'*s) + Td'*s)
     *
     * Parallel (non-interacting) PID:
     *   C(s) = Kp + Ki/s + Kd*s
     *
     * Matching coefficients:
     *   Kp = Kp' * (1 + Td'/Ti')
     *   Ki = Kp' / Ti'
     *   Kd = Kp' * Td'
     */
    if (Ti_s < 1e-12) {
        *Kp_p = Kp_s;
        *Ki_p = 0.0;
        *Kd_p = Kp_s * Td_s;
        return;
    }

    *Kp_p = Kp_s * (1.0 + Td_s / Ti_s);
    *Ki_p = Kp_s / Ti_s;
    *Kd_p = Kp_s * Td_s;
}

/* ──────────────────────────────────────────────
 * L2: Core PID Update — Parallel Form with AW
 * ────────────────────────────────────────────── */

double pid_update(pid_controller_t *pid, double setpoint, double measurement)
{
    if (!pid) return 0.0;

    /* Apply action direction */
    double sp = setpoint;
    double y  = measurement;
    double e  = pid->action * (sp - y);

    /* Proportional term (with setpoint weighting) */
    double P = pid->params.Kp * (pid->params.b * sp - y);

    /* --- Integral term --- */
    double Ki_Ts = pid->params.Ki * pid->Ts;
    /* Back-calculation anti-windup */
    if (pid->aw_mode == AW_MODE_BACK_CALC && pid->params.Tt > 1e-12) {
        double tracking_correction = (pid->Ts / pid->params.Tt)
                                     * (pid->u_sat - pid->u_unsat);
        pid->I += Ki_Ts * e + tracking_correction;
    } else if (pid->aw_mode == AW_MODE_CLAMPING) {
        /* Conditional integration: only integrate if not saturated
           or if error drives output away from saturation */
        int sat_high = (pid->u_sat >= pid->u_max) && (e > 0);
        int sat_low  = (pid->u_sat <= pid->u_min) && (e < 0);
        if (!sat_high && !sat_low) {
            pid->I += Ki_Ts * e;
        }
    } else if (pid->aw_mode == AW_MODE_VELOCITY) {
        /* Velocity form: handled separately; fall through to standard */
        pid->I += Ki_Ts * e;
    } else {
        /* No anti-windup or unknown mode: standard integration */
        pid->I += Ki_Ts * e;
    }

    double I_term = pid->I;

    /* --- Derivative term (with filtering and measurement-only option) --- */
    double D_term = 0.0;
    double Td = pid->params.Td;
    double N  = pid->params.N;

    if (Td > 1e-12 && N > 1e-12) {
        /* Filter coefficient: α = Td / (Td + N*Ts) */
        double alpha = Td / (Td + N * pid->Ts);

        if (pid->form == PID_FORM_DERIVATIVE_ON_PV) {
            /* Derivative on measurement: D = Kd * d(PV)/dt */
            double dy = y - pid->y_prev;
            double raw_D = -pid->params.Kd * dy / pid->Ts;
            pid->D_filt = alpha * pid->D_filt + (1.0 - alpha) * raw_D;
            D_term = pid->D_filt;
        } else {
            /* Derivative on error */
            double de = e - pid->e_prev;
            double raw_D = pid->params.Kd * de / pid->Ts;
            pid->D_filt = alpha * pid->D_filt + (1.0 - alpha) * raw_D;
            D_term = pid->D_filt;
        }
    }

    /* --- Summation --- */
    double u_unsat = P + I_term + D_term;

    /* --- Saturation --- */
    double u = u_unsat;
    if (u > pid->u_max)  u = pid->u_max;
    if (u < pid->u_min)  u = pid->u_min;

    /* Update controller state for next iteration */
    pid->e_prev  = e;
    pid->y_prev  = y;
    pid->u_unsat = u_unsat;
    pid->u_sat   = u;
    pid->step_count++;

    return u;
}

/* ──────────────────────────────────────────────
 * L5: Reset PID State (Bumpless Transfer Prep)
 * ────────────────────────────────────────────── */

void pid_reset(pid_controller_t *pid)
{
    if (!pid) return;
    pid->I      = 0.0;
    pid->e_prev = 0.0;
    pid->y_prev = 0.0;
    pid->u_sat  = 0.0;
    pid->u_unsat= 0.0;
    pid->step_count = 0;
}

/* ──────────────────────────────────────────────
 * L3: PID Frequency Response
 * ────────────────────────────────────────────── */

void pid_freq_response(const pid_params_t *params, double omega,
                       double *mag, double *phase)
{
    /**
     * PID transfer function (parallel form with derivative filter):
     *   C(jω) = Kp + Ki/(jω) + Kd*jω / (1 + jω*Td/N)
     *
     * Separate into real and imaginary parts:
     *   Re: Kp + Kd*ω²*(Td/N) / (1 + ω²*(Td/N)²)
     *   Im: -Ki/ω + Kd*ω / (1 + ω²*(Td/N)²)
     */
    double Kp = params->Kp;
    double Ki = params->Ki;
    double Kd = params->Kd;
    double Td = params->Td;
    double N  = params->N;

    double omega_sq = omega * omega;
    double TdN = (Td > 1e-12 && N > 1e-12) ? (Td / N) : 0.0;
    double denom = 1.0 + omega_sq * TdN * TdN;

    double Re = Kp;
    double Im = 0.0;

    if (omega > 1e-12) {
        Im = -Ki / omega;
    }

    /* Derivative contribution */
    if (Kd > 1e-12 && denom > 1e-12) {
        Re += Kd * omega_sq * TdN / denom;
        Im += Kd * omega / denom;
    }

    *mag   = sqrt(Re * Re + Im * Im);
    *phase = atan2(Im, Re);
}

/* ──────────────────────────────────────────────
 * L3: PID Transfer Function Coefficients
 * ────────────────────────────────────────────── */

void pid_transfer_function_coeffs(const pid_params_t *p,
                                  double *a2, double *a1, double *a0)
{
    /**
     * PID parallel form with derivative filter:
     *   C(s) = Kp + Ki/s + Kd*s / (1 + s*Td/N)
     *         = [Kp*(1 + s*Td/N)*s + Ki*(1 + s*Td/N) + Kd*s²] / [s*(1 + s*Td/N)]
     *
     * Numerator: (Kd + Kp*Td/N)*s² + (Kp + Ki*Td/N)*s + Ki
     *
     * So:
     *   a2 = Kd + Kp * Td / N
     *   a1 = Kp + Ki * Td / N
     *   a0 = Ki
     */
    double Td = p->Td;
    double N  = p->N;
    double TdN = (Td > 1e-12 && N > 1e-12) ? (Td / N) : 0.0;

    *a2 = p->Kd + p->Kp * TdN;
    *a1 = p->Kp + p->Ki * TdN;
    *a0 = p->Ki;
}

/* ──────────────────────────────────────────────
 * L2: Performance Metric Evaluation
 * ────────────────────────────────────────────── */

void pid_eval_performance(const double *t, const double *y,
                          double sp, int N, pid_perf_t *perf)
{
    if (!t || !y || !perf || N < 2) return;
    memset(perf, 0, sizeof(*perf));

    perf->IAE  = 0.0;
    perf->ISE  = 0.0;
    perf->ITAE = 0.0;
    perf->ITSE = 0.0;

    double y_max = y[0];
    double y_min = y[0];
    int oscillation_sign_changes = 0;

    for (int i = 0; i < N; i++) {
        double e_i = sp - y[i];
        double abs_e = fabs(e_i);
        double t_i = t[i];

        perf->IAE  += abs_e;
        perf->ISE  += e_i * e_i;
        perf->ITAE += t_i * abs_e;
        perf->ITSE += t_i * e_i * e_i;

        if (y[i] > y_max) y_max = y[i];
        if (y[i] < y_min) y_min = y[i];

        /* Detect sign changes for oscillation counting */
        if (i > 0) {
            double e_curr = sp - y[i];
            double e_prev = sp - y[i-1];
            if (e_curr * e_prev < 0) oscillation_sign_changes++;
        }
    }

    if (sp > 1e-12) {
        perf->overshoot_pct = (y_max > sp) ? 100.0 * (y_max - sp) / sp : 0.0;
    } else {
        perf->overshoot_pct = (y_max > 0.01) ? 100.0 * (y_max / 0.01) : 0.0;
    }

    /* Settling time: find when output stays within ±2% of final value */
    double y_final = y[N-1];
    double band = 0.02 * fabs(sp);
    if (band < 1e-9) band = 0.02;
    perf->settling_time = t[N-1]; /* default = end of data */
    for (int i = N - 1; i >= 0; i--) {
        if (fabs(y[i] - y_final) > band) {
            perf->settling_time = (i + 1 < N) ? t[i + 1] : t[i];
            break;
        }
    }

    /* Rise time: 10% to 90% of final value */
    double y10 = 0.1 * y_final;
    double y90 = 0.9 * y_final;
    double t10 = t[0], t90 = t[N-1];
    int found10 = 0, found90 = 0;
    for (int i = 0; i < N; i++) {
        if (!found10 && y[i] >= y10) { t10 = t[i]; found10 = 1; }
        if (!found90 && y[i] >= y90) { t90 = t[i]; found90 = 1; }
    }
    perf->rise_time = t90 - t10;

    /* Steady-state error */
    perf->steady_state_error = sp - y_final;
    perf->oscillation_count = oscillation_sign_changes;

    /* Normalize integrals by time span for fair comparison */
    if (t[N-1] > 1e-12) {
        perf->IAE  /= t[N-1];
        perf->ITAE /= (t[N-1] * t[N-1]);
    }
}

/* ──────────────────────────────────────────────
 * L5: Configure Limits and Anti-Windup
 * ────────────────────────────────────────────── */

void pid_set_limits(pid_controller_t *pid, antiwindup_mode_t mode,
                    double u_min, double u_max)
{
    if (!pid) return;
    pid->aw_mode = mode;
    pid->u_min   = u_min;
    pid->u_max   = u_max;
}

void pid_set_setpoint_weights(pid_controller_t *pid, double b, double c)
{
    if (!pid) return;
    pid->params.b = b;
    pid->params.c = c;
}

/* ──────────────────────────────────────────────
 * L5: Apply Tuning Result to Controller
 * ────────────────────────────────────────────── */

void pid_apply_tuning(pid_controller_t *pid, const pid_params_t *p,
                      pid_form_t form, double Ts)
{
    if (!pid || !p) return;

    pid->form = form;
    pid->params.Kp = p->Kp;
    pid->params.Ti = p->Ti;
    pid->params.Td = p->Td;
    pid->params.N  = p->N;
    pid->params.b  = p->b;
    pid->params.c  = p->c;
    pid->params.Tt = p->Tt;

    /* Recompute parallel gains based on form */
    if (form == PID_FORM_PARALLEL || form == PID_FORM_DERIVATIVE_ON_PV) {
        pid->params.Ki = p->Ki;
        pid->params.Kd = p->Kd;
    } else {
        pid_ideal_to_parallel(p->Kp, p->Ti, p->Td,
                              &pid->params.Ki, &pid->params.Kd);
    }

    pid->Ts = Ts;
    pid_reset(pid);
}

/* ──────────────────────────────────────────────
 * L3: Form Conversion Between PID Representations
 * ────────────────────────────────────────────── */

void pid_convert_form(const pid_params_t *src, pid_form_t src_form,
                      pid_params_t *dst, pid_form_t dst_form)
{
    if (!src || !dst) return;
    memset(dst, 0, sizeof(*dst));

    /* First, ensure src is in parallel representation */
    double Kp, Ki, Kd;
    if (src_form == PID_FORM_IDEAL) {
        Kp = src->Kp;
        pid_ideal_to_parallel(Kp, src->Ti, src->Td, &Ki, &Kd);
    } else if (src_form == PID_FORM_SERIES) {
        pid_series_to_parallel(src->Kp, src->Ti, src->Td, &Kp, &Ki, &Kd);
    } else {
        /* Already parallel or custom */
        Kp = src->Kp;
        Ki = src->Ki;
        Kd = src->Kd;
    }

    /* Now convert to destination form */
    if (dst_form == PID_FORM_IDEAL) {
        dst->Kp = Kp;
        pid_parallel_to_ideal(Kp, Ki, Kd, &dst->Ti, &dst->Td);
    } else if (dst_form == PID_FORM_SERIES) {
        /* Series is approximate — only valid if Td << Ti */
        dst->Kp = Kp;
        if (Ki > 1e-12) {
            dst->Ti = Kp / Ki;
            dst->Td = (Kd > 1e-12) ? (Kd / Kp) : 0.0;
        } else {
            dst->Ti = 0.0;
            dst->Td = 0.0;
        }
    } else {
        dst->Kp = Kp;
        dst->Ki = Ki;
        dst->Kd = Kd;
        pid_parallel_to_ideal(Kp, Ki, Kd, &dst->Ti, &dst->Td);
    }

    /* Copy ancillary parameters */
    dst->N  = src->N;
    dst->b  = src->b;
    dst->c  = src->c;
    dst->Tt = src->Tt;
}

/* ──────────────────────────────────────────────
 * L5: Filter and Sampling Utilities
 * ────────────────────────────────────────────── */

double pid_derive_filter_N(double Td, double cutoff)
{
    /**
     * The filter term is: 1 / (1 + s*Td/N)
     * Its cutoff frequency: ωc = N / Td
     * Therefore: N = ωc * Td
     */
    if (Td <= 1e-12) return 10.0;
    double N = cutoff * Td;
    if (N < 2.0)  N = 2.0;   /* Minimum reasonable N */
    if (N > 50.0) N = 50.0;  /* Maximum reasonable N */
    return N;
}

double pid_recommend_sampling(double omega_bw)
{
    /**
     * Nyquist criterion requires Ts ≤ π/ω_max.
     * For PID, the closed-loop bandwidth is approximately ω_bw.
     * Rule of thumb: 4-10 samples per time constant of fastest mode.
     *   Ts ≤ 1 / (10 * ω_bw)
     * More conservative: Ts ≤ 0.1 / ω_bw
     */
    if (omega_bw <= 1e-12) return 0.01;
    return 0.1 / omega_bw;
}
