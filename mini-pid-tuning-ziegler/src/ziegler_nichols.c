/**
 * ziegler_nichols.c — Ziegler-Nichols PID Tuning Implementation
 *
 * Implements both the step response (1942) and frequency response
 * (ultimate sensitivity) methods of Ziegler and Nichols.
 *
 * Knowledge: L1 (ZN method definitions), L2 (process reaction curve,
 * ultimate gain concept), L4 (empirical rule derivation),
 * L5 (algorithmic tuning), L6 (FOPDT identification).
 *
 * Reference:
 *   Ziegler & Nichols (1942), Trans. ASME, 64, 759-768.
 *   Åström & Hägglund (2004), J. Process Control, 14, 635-650.
 */

#include "ziegler_nichols.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L4: Z-N Tuning Rule Tables
 * ────────────────────────────────────────────── */

/**
 * Original Ziegler-Nichols step response rule coefficients.
 *
 * These give quarter-amplitude damping (ζ ≈ 0.21, decay ratio ≈ 0.25).
 *
 * Controller | Kp_factor (× T/(K*L)) | Ti_factor (× L) | Td_factor (× L)
 * ----------|------------------------|----------------|---------------
 * P          | 1.0                   | —               | —
 * PI         | 0.9                   | 3.333          | —
 * PID        | 1.2                   | 2.0            | 0.5
 *
 * Note: In the literature, PI Ti is sometimes quoted as 3.33*L.
 *       Here we use 3.0*L as a common variant (Åström compilation).
 */

static const zn_rule_entry_t ZN_STEP_RULES[4] = {
    /* ZN_CONTROLLER_P    */  { 1.00, 0.0,  0.0  },
    /* ZN_CONTROLLER_PI   */  { 0.90, 3.333, 0.0  },
    /* ZN_CONTROLLER_PD   */  { 0.80, 0.0,  0.33 }, /* Not in original; Åström extension */
    /* ZN_CONTROLLER_PID  */  { 1.20, 2.0,  0.5  }
};

/**
 * Original Ziegler-Nichols frequency response rule coefficients.
 *
 * Controller | Kp_factor (× Ku) | Ti_factor (× Pu) | Td_factor (× Pu)
 * ----------|------------------|----------------|---------------
 * P          | 0.50             | —               | —
 * PI         | 0.45             | 0.833 (Pu/1.2) | —
 * PID        | 0.60             | 0.50 (Pu/2.0)  | 0.125 (Pu/8)
 */

static const zn_freq_rule_entry_t ZN_FREQ_RULES[4] = {
    /* ZN_CONTROLLER_P    */  { 0.50, 0.0,   0.0   },
    /* ZN_CONTROLLER_PI   */  { 0.45, 0.833, 0.0   },
    /* ZN_CONTROLLER_PD   */  { 0.55, 0.0,   0.15  },
    /* ZN_CONTROLLER_PID  */  { 0.60, 0.50,  0.125 }
};

/* Chien-Hrones-Reswick (1952) rule coefficients */
static const zn_rule_entry_t CHR_RULES_SP[4] = {
    /* Setpoint tracking, 0% overshoot */
    { 0.30, 1.0,  0.0  },  /* P  */
    { 0.35, 4.0,  0.0  },  /* PI */
    { 0.60, 0.0,  0.33 },  /* PD */
    { 0.95, 2.4,  0.42 }   /* PID */
};

static const zn_rule_entry_t CHR_RULES_DIST[4] = {
    /* Disturbance rejection, 20% overshoot */
    { 0.70, 1.0,  0.0  },  /* P  */
    { 0.60, 4.0,  0.0  },  /* PI */
    { 1.00, 0.0,  0.50 },  /* PD */
    { 1.20, 2.0,  0.42 }   /* PID */
};

/* ──────────────────────────────────────────────
 * L5: Find Inflection Point in Step Response
 * ────────────────────────────────────────────── */

static int find_inflection_point(const double *y, int N, double Ts,
                                 int *idx_inf, double *max_slope)
{
    if (N < 5) return -1;

    double best_slope = 0.0;
    int    best_idx   = 0;

    /* Compute finite differences; find maximum slope */
    for (int i = 1; i < N; i++) {
        double slope = (y[i] - y[i-1]) / Ts;
        if (slope > best_slope) {
            best_slope = slope;
            best_idx   = i;
        }
    }

    if (best_slope <= 1e-12) return -1;

    *idx_inf   = best_idx;
    *max_slope = best_slope;
    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Z-N FOPDT Identification (Tangent Method)
 * ────────────────────────────────────────────── */

int zn_identify_fopdt(const step_response_data_t *data, fopdt_model_t *model)
{
    if (!data || !model) return -1;
    if (data->N < 10)   return -1;
    if (data->step_mag < 1e-12) return -1;

    int    idx_inf;
    double max_slope;
    if (find_inflection_point(data->output, data->N, data->Ts,
                              &idx_inf, &max_slope) != 0) {
        return -1;
    }

    double y_inf = data->output[idx_inf];
    double t_inf = data->time[idx_inf];

    /* Tangent line: y = y_inf + max_slope * (t - t_inf) */
    /* Intersection with baseline (y = y[0]): L = t_inf - y_inf/max_slope */
    double y0 = data->output[0];
    double L = t_inf - (y_inf - y0) / max_slope;
    if (L < 0.0) L = 0.0;

    /* Steady-state value */
    double y_ss = data->output[data->N - 1];

    /* T: time for tangent to reach y_ss from baseline */
    double T = (y_ss - y0) / max_slope - L;
    if (T <= 0.0) T = data->Ts; /* Minimum time constant */

    /* Static gain */
    double K = (y_ss - y0) / data->step_mag;

    model->K = K;
    model->T = T;
    model->L = L;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Z-N Step Response Tuning
 * ────────────────────────────────────────────── */

int zn_step_tune(const fopdt_model_t *model, zn_controller_type_t type,
                 pid_params_t *params)
{
    if (!model || !params) return -1;
    if (model->T <= 1e-12 || model->K < 1e-12) return -1;
    if (model->L < 0.0) return -1;

    int idx = (int)type;
    if (idx < 0 || idx > 3) return -1;

    const zn_rule_entry_t *rule = &ZN_STEP_RULES[idx];

    double K = model->K;
    double T = model->T;
    double L = model->L;

    /* a = K * L / T */
    if (L < 1e-12) {
        /* For processes with negligible dead time, use conservative defaults */
        params->Kp = 0.5 * T / K;
        params->Ti = T;
        params->Td = 0.0;
        params->Kd = 0.0;
    } else {
        double a = K * L / T;
        if (a < 1e-12) {
            /* Avoid division by zero */
            params->Kp = 1.0;
            params->Ti = 1.0;
            params->Td = 0.0;
            params->Kd = 0.0;
        } else {
            params->Kp = rule->Kp_factor / a;
            params->Ti = rule->Ti_factor * L;
            params->Td = rule->Td_factor * L;
        }
    }

    /* Compute parallel-form gains */
    pid_ideal_to_parallel(params->Kp, params->Ti, params->Td,
                          &params->Ki, &params->Kd);

    /* Default filter and weighting */
    params->N = 10.0;
    params->b = 1.0;
    params->c = 0.0;
    params->Tt = sqrt(params->Ti * params->Td);
    if (params->Tt < 1e-12) params->Tt = params->Ti;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Modified Z-N Step Response Tuning
 * ────────────────────────────────────────────── */

int zn_step_tune_modified(const fopdt_model_t *model, zn_controller_type_t type,
                          zn_modified_t variant, double Ms,
                          pid_params_t *params)
{
    if (zn_step_tune(model, type, params) != 0) return -1;

    /**
     * The original Z-N rules target ζ ≈ 0.21 (decay ratio ≈ 0.25).
     * Modified rules adjust Kp to achieve different damping.
     *
     * Approximate relationship between damping and Z-N Kp multiplier:
     *   ζ = 0.21 → multiplier = 1.0  (standard)
     *   ζ = 0.50 → multiplier = 0.5  (conservative)
     *   ζ = 0.10 → multiplier = 1.5  (aggressive)
     *   ζ = 1.00 → multiplier = 0.2  (no overshoot)
     *   ζ = 0.35 → multiplier = 0.75 (some overshoot, ~10-20%)
     */

    double Kp_mult;
    switch (variant) {
        case ZN_MOD_CONSERVATIVE:
            Kp_mult = 0.50; break;
        case ZN_MOD_AGGRESSIVE:
            Kp_mult = 1.50; break;
        case ZN_MOD_NO_OVERSHOOT:
            Kp_mult = 0.20; break;
        case ZN_MOD_SOME_OVERSHOOT:
            Kp_mult = 0.75; break;
        case ZN_MOD_STANDARD:
        default:
            Kp_mult = 1.00; break;
    }

    /* Adjust gains proportionally */
    params->Kp *= Kp_mult;
    params->Ki *= Kp_mult;
    params->Kd *= Kp_mult;

    /* Tune Ti/Td for Ms if requested */
    if (Ms > 1.2 && Ms < 2.5 && variant != ZN_MOD_STANDARD) {
        /* Approximate relationship: higher Ms → faster Ti to compensate */
        double Ms_factor = (Ms - 1.2) / 0.8; /* Map to [0, 1] */
        params->Ti *= (1.0 - 0.3 * Ms_factor);
        if (params->Ti < model->L) params->Ti = model->L;
    }

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Z-N Step Response Tuning for SOPDT
 * ────────────────────────────────────────────── */

int zn_step_tune_sopdt(const sopdt_model_t *model, zn_controller_type_t type,
                       pid_params_t *params)
{
    if (!model || !params) return -1;

    /**
     * SOPDT → FOPDT approximation:
     * If T1 and T2 are both real:
     *   T_eff = max(T1, T2) + min(T1, T2)/2  (Skogestad half-rule)
     *   L_eff = L + min(T1, T2)/2
     *
     * If oscillatory (zeta, wn):
     *   T_eff = 2*zeta/wn  (dominant time constant)
     *   L_eff = L (dead time unchanged, but dynamics are different)
     */

    fopdt_model_t fopdt;
    double T1 = model->T1;
    double T2 = model->T2;

    if (model->use_osc) {
        double zeta = model->zeta;
        double wn   = model->wn;
        if (zeta > 1.0) {
            /* Overdamped oscillatory → real poles */
            double tau1 = 0.0, tau2 = 0.0;
            double disc = sqrt(zeta * zeta - 1.0);
            tau1 = 1.0 / (wn * (zeta - disc));
            tau2 = 1.0 / (wn * (zeta + disc));
            T1 = tau1;
            T2 = tau2;
        } else {
            /* Underdamped: use dominant time constant approximation */
            T1 = 2.0 * zeta / wn;
            T2 = 0.1 * T1; /* residual */
        }
    }

    /* Apply half-rule */
    double T_dom, T_sub;
    if (T1 >= T2) { T_dom = T1; T_sub = T2; }
    else          { T_dom = T2; T_sub = T1; }

    fopdt.K = model->K;
    fopdt.T = T_dom + 0.5 * T_sub;
    fopdt.L = model->L + 0.5 * T_sub;

    return zn_step_tune(&fopdt, type, params);
}

/* ──────────────────────────────────────────────
 * L5: Find Ultimate Gain (Frequency Domain Analysis)
 * ────────────────────────────────────────────── */

int zn_find_ultimate_gain(const fopdt_model_t *fopdt,
                          ultimate_gain_result_t *result)
{
    if (!fopdt || !result) return -1;
    if (fopdt->T <= 1e-12) return -1;

    /**
     * The phase crossover frequency ω_pc solves:
     *   ∠G(jω_pc) = -ω_pc * L - arctan(ω_pc * T) = -π
     *
     * This is a transcendental equation. We solve numerically:
     *   f(ω) = ω * L + arctan(ω * T) - π = 0
     *
     * At ω = ω_pc:
     *   Ku = 1 / |G(jω_pc)| = sqrt(1 + ω_pc²*T²) / K
     */

    double K = fopdt->K;
    double T = fopdt->T;
    double L = fopdt->L;

    if (K < 1e-12) return -1;

    /* For processes with no dead time, the phase asymptotically approaches
       -π/2 (for FOPDT) or -π (for SOPDT). FOPDT without L can't cross -π. */
    if (L < 1e-12) {
        /* No dead time → no phase crossover → Ku is infinite (no oscillation) */
        result->Ku = 1e6;
        result->Pu = 1e6;
        result->converged = 0;
        return 0;
    }

    /* Bisection search for ω_pc in [ω_low, ω_high] */
    double omega_low  = 0.01 / (L + T);
    double omega_high = M_PI / L;  /* At this ω, phase ≤ -π */
    double omega_mid;

    for (int iter = 0; iter < 60; iter++) {
        omega_mid = 0.5 * (omega_low + omega_high);
        double phase = -omega_mid * L - atan(omega_mid * T);
        if (phase > -M_PI) {
            omega_low = omega_mid;
        } else {
            omega_high = omega_mid;
        }
        if (omega_high - omega_low < 1e-9) break;
    }

    double omega_pc = omega_mid;
    double mag_G = K / sqrt(1.0 + omega_pc * omega_pc * T * T);

    result->Ku = 1.0 / mag_G;
    result->Pu = (omega_pc > 1e-12) ? (2.0 * M_PI / omega_pc) : 1e6;
    result->converged = 1;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Z-N Frequency Response Tuning
 * ────────────────────────────────────────────── */

int zn_freq_tune(double Ku, double Pu, zn_controller_type_t type,
                 pid_params_t *params)
{
    if (!params) return -1;
    if (Ku <= 1e-12 || Pu <= 1e-12) return -1;

    int idx = (int)type;
    if (idx < 0 || idx > 3) return -1;

    const zn_freq_rule_entry_t *rule = &ZN_FREQ_RULES[idx];

    params->Kp = rule->Kp_factor * Ku;
    params->Ti = rule->Ti_factor * Pu;
    params->Td = rule->Td_factor * Pu;

    pid_ideal_to_parallel(params->Kp, params->Ti, params->Td,
                          &params->Ki, &params->Kd);

    params->N  = 10.0;
    params->b  = 1.0;
    params->c  = 0.0;
    params->Tt = sqrt(params->Ti * params->Td);
    if (params->Tt < 1e-12) params->Tt = params->Ti;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Custom Z-N Frequency Response Tuning
 * ────────────────────────────────────────────── */

int zn_freq_tune_custom(double Ku, double Pu, zn_controller_type_t type,
                        double alpha, double beta, double gamma,
                        pid_params_t *params)
{
    if (!params) return -1;
    if (Ku <= 1e-12 || Pu <= 1e-12) return -1;

    int idx = (int)type;
    if (idx < 0 || idx > 3) return -1;

    params->Kp = alpha * Ku;
    params->Ti = beta  * Pu;
    params->Td = gamma * Pu;

    pid_ideal_to_parallel(params->Kp, params->Ti, params->Td,
                          &params->Ki, &params->Kd);

    params->N  = 10.0;
    params->b  = 1.0;
    params->c  = 0.0;
    params->Tt = (params->Ti > 1e-12 && params->Td > 1e-12)
                 ? sqrt(params->Ti * params->Td) : params->Ti;

    return 0;
}

/* ──────────────────────────────────────────────
 * L4: Margin Verification for Given Z-N Tuning
 * ────────────────────────────────────────────── */

int zn_verify_margins(const fopdt_model_t *fopdt, const pid_params_t *params,
                      double *Gm, double *Pm)
{
    if (!fopdt || !params || !Gm || !Pm) return -1;

    /**
     * Compute open-loop transfer function: L(s) = G(s) * C(s)
     * Find ω_pc: ∠L(jω_pc) = -π → Gm = 1 / |L(jω_pc)|
     * Find ω_gc: |L(jω_gc)| = 1 → Pm = π + ∠L(jω_gc)
     */

    double K = fopdt->K;
    double T = fopdt->T;
    double L = fopdt->L;
    double Kp = params->Kp;
    double Ti = params->Ti;
    double Td = params->Td;

    if (K < 1e-12 || Kp < 1e-12) return -1;

    /* Phase crossover: ∠L(ω) = -ωL - π/2 at low freq (integrator)
       + arctan(...) from PI/PID. For PI only, simplified.
       Solve: -ωL - arctan(ωT) + arctan(ωTi) - π/2 = -π  (for PI)
       Solve: -ωL - arctan(ωT) + arctan(ωTi) + arctan(ωTd) - π/2 = -π (for PID)
    */

    /* Numerical bisection for ω_pc */
    double w_lo = 0.001 / (L + T + 1.0);
    double w_hi = M_PI / (L + 0.001);
    double omega_pc = 0.0;
    double mag_L_pc = 0.0;

    for (int iter = 0; iter < 50; iter++) {
        double w = 0.5 * (w_lo + w_hi);
        double phase_G = -w * L - atan(w * T);
        double mag_G   = K / sqrt(1.0 + w * w * T * T);

        double phase_C;
        if (Ti < 1e-12 && Td < 1e-12) {
            /* P only */
            phase_C = 0.0;
        } else if (Td < 1e-12) {
            /* PI */
            phase_C = -M_PI / 2.0 + atan(w * Ti);
        } else {
            /* PID */
            phase_C = -M_PI / 2.0 + atan(w * Ti) + atan(w * Td);
        }

        double phase_L = phase_G + phase_C;
        if (phase_L > -M_PI) {
            w_lo = w;
        } else {
            w_hi = w;
            omega_pc = w;
            /* Magnitude of PI/PID at ω:
               |C| = Kp * sqrt(1 + (ωTi - 1/(ωTi)...)
               Simplified for PI: |C| = Kp * sqrt(1 + 1/(ω²Ti²))
               Simplified for PID: |C| = Kp * sqrt(1 + (ωTd - 1/(ωTi))²)
               */
            double mag_C;
            if (Ti < 1e-12) {
                mag_C = Kp;
            } else if (Td < 1e-12) {
                mag_C = Kp * sqrt(1.0 + 1.0 / (w * w * Ti * Ti));
            } else {
                double x = w * Td - 1.0 / (w * Ti);
                mag_C = Kp * sqrt(1.0 + x * x);
            }
            mag_L_pc = mag_G * mag_C;
        }

        if (w_hi - w_lo < 1e-8) break;
    }

    if (omega_pc > 1e-12) {
        *Gm = 1.0 / mag_L_pc;
    } else {
        *Gm = 1e6; /* No phase crossover → infinite gain margin */
    }

    /* Gain crossover: solve |L(ω)| = 1 using similar bisection */
    /* For simplicity, estimate from known relationships */
    w_lo = 0.001 / (L + T);
    w_hi = 10.0 / L; /* if L > 0, else wide range */
    double omega_gc = 0.0;
    double phase_at_gc = 0.0;

    for (int iter = 0; iter < 50; iter++) {
        double w = 0.5 * (w_lo + w_hi);
        double mag_G = K / sqrt(1.0 + w * w * T * T);

        double mag_C;
        if (Ti < 1e-12) {
            mag_C = Kp;
        } else if (Td < 1e-12) {
            mag_C = Kp * sqrt(1.0 + 1.0 / (w * w * Ti * Ti));
        } else {
            double x = w * Td - 1.0 / (w * Ti);
            mag_C = Kp * sqrt(1.0 + x * x);
        }

        double mag_L = mag_G * mag_C;
        if (mag_L > 1.0) {
            w_lo = w;
        } else {
            w_hi = w;
            omega_gc = w;
            /* Phase at this frequency */
            double phase_G = -w * L - atan(w * T);
            double phase_C = (Ti < 1e-12) ? 0.0 :
                             -M_PI/2.0 + atan(w * Ti) +
                             ((Td > 1e-12) ? atan(w * Td) : 0.0);
            phase_at_gc = phase_G + phase_C;
        }

        if (w_hi - w_lo < 1e-8) break;
    }

    *Pm = M_PI + phase_at_gc;

    (void)omega_gc; /* Computed, reported indirectly through Pm */
    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Sørensen's Method for Ultimate Gain
 * ────────────────────────────────────────────── */

int zn_sorensen_ultimate_gain(double relay_amplitude, double hysteresis_width,
                              double oscillation_ampl, double oscillation_period,
                              ultimate_gain_result_t *result)
{
    if (!result) return -1;
    if (oscillation_ampl < 1e-12 || oscillation_period < 1e-12) return -1;
    if (relay_amplitude < 1e-12) return -1;

    /**
     * Sørensen (1958) modification: relay with hysteresis ε.
     *
     * For ideal relay (ε = 0):
     *   Ku = 4*d / (π * a)
     *
     * For relay with hysteresis ε:
     *   N(a) = (4*d / (π*a)) * sqrt(1 - (ε/a)²) - j*(4*d*ε) / (π*a²)
     *   |N(a)| = 4*d / (π * a)  (magnitude unchanged for ideal portion)
     *
     * Simplified: Ku ≈ 4*d / (π * a) still holds for small ε.
     * Pu = P_osc (oscillation period, approximately).
     */

    double d = relay_amplitude;
    double a = oscillation_ampl;
    double e = hysteresis_width;

    if (e > 0.0 && e < a) {
        /* Corrected describing function magnitude: */
        double ratio = e / a;
        double mag_N = (4.0 * d) / (M_PI * a) * sqrt(1.0 - ratio * ratio);
        result->Ku = 1.0 / mag_N;
    } else {
        result->Ku = (4.0 * d) / (M_PI * a);
    }

    result->Pu = oscillation_period;
    result->converged = 1;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Chien-Hrones-Reswick (1952) Tuning
 * ────────────────────────────────────────────── */

int zn_chien_hrones_reswick(const fopdt_model_t *model,
                            zn_controller_type_t type, int mode,
                            pid_params_t *params)
{
    if (!model || !params) return -1;

    if (model->T <= 1e-12 || model->K < 1e-12 || model->L < 0) return -1;

    int idx = (int)type;
    if (idx < 0 || idx > 3) return -1;

    const zn_rule_entry_t *rule = (mode == 0) ? &CHR_RULES_SP[idx]
                                               : &CHR_RULES_DIST[idx];

    double K = model->K;
    double T = model->T;
    double L = model->L;

    /**
     * CHR rules use the same form as Z-N but with different coefficients:
     *   Kp = Kp_factor * T / (K * L)
     *   Ti = Ti_factor * L
     *   Td = Td_factor * L
     */

    if (L < 1e-12) {
        /* Use minimum effective dead time */
        L = 0.01 * T;
    }

    params->Kp = rule->Kp_factor * T / (K * L);
    params->Ti = rule->Ti_factor * L;
    params->Td = rule->Td_factor * L;

    pid_ideal_to_parallel(params->Kp, params->Ti, params->Td,
                          &params->Ki, &params->Kd);

    params->N  = 10.0;
    params->b  = 1.0;
    params->c  = 0.0;
    params->Tt = sqrt(params->Ti * params->Td);
    if (params->Tt < 1e-12) params->Tt = params->Ti;

    return 0;
}
