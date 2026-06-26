/**
 * cohen_coon.c — Cohen-Coon PID Tuning Method Implementation
 *
 * Implements the Cohen-Coon (1953) tuning method, which provides
 * PID parameters for FOPDT processes based on closed-loop pole
 * placement with dead-time ratio compensation.
 *
 * Knowledge: L1 (CC definition), L2 (decay ratio), L4 (pole placement),
 *             L5 (algorithmic computation), L6 (FOPDT with dead-time).
 *
 * Reference:
 *   Cohen, G.H. & Coon, G.A. (1953), Trans. ASME, 75, 827-834.
 */

#include "cohen_coon.h"
#include "ziegler_nichols.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L5: Cohen-Coon Tuning Formulas (Main)
 * ────────────────────────────────────────────── */

int cohen_coon_tune(const fopdt_model_t *model, cc_controller_type_t type,
                    pid_params_t *params)
{
    if (!model || !params) return -1;
    if (model->K <= 1e-12 || model->T <= 1e-12 || model->L <= 1e-12) return -1;

    double K = model->K;
    double T = model->T;
    double L = model->L;

    /* Normalized dead time ratio: μ = L / T */
    double mu = L / T;

    /**
     * Cohen-Coon formulas:
     *
     * P:
     *   Kp = (T/(K*L)) * (1 + 0.333*μ)
     *
     * PI:
     *   Kp = (T/(K*L)) * (0.9 + 0.083*μ)
     *   Ti = L * (3.33 + 0.3*μ) / (1 + 2.2*μ)
     *
     * PD:
     *   Kp = (T/(K*L)) * (1.24 + 0.12*μ)
     *   Td = L * (0.27 - 0.036*μ) / (1 + 0.87*μ)
     *
     * PID:
     *   Kp = (T/(K*L)) * (1.35 + 0.27*μ)
     *   Ti = L * (2.5 + 0.46*μ) / (1 + 0.61*μ)
     *   Td = L * (0.37) / (1 + 0.19*μ)
     */

    double Kp = 0.0, Ti = 0.0, Td = 0.0;

    switch (type) {
        case CC_CONTROLLER_P:
            Kp = (T / (K * L)) * (1.0 + 0.333 * mu);
            Ti = 0.0;
            Td = 0.0;
            break;

        case CC_CONTROLLER_PI:
            Kp = (T / (K * L)) * (0.9 + 0.083 * mu);
            Ti = L * (3.33 + 0.3 * mu) / (1.0 + 2.2 * mu);
            Td = 0.0;
            break;

        case CC_CONTROLLER_PD:
            Kp = (T / (K * L)) * (1.24 + 0.12 * mu);
            Ti = 0.0;
            Td = L * (0.27 - 0.036 * mu) / (1.0 + 0.87 * mu);
            if (Td < 0.0) Td = 0.0;
            break;

        case CC_CONTROLLER_PID:
            Kp = (T / (K * L)) * (1.35 + 0.27 * mu);
            Ti = L * (2.5 + 0.46 * mu) / (1.0 + 0.61 * mu);
            Td = L * (0.37) / (1.0 + 0.19 * mu);
            break;

        default:
            return -1;
    }

    params->Kp = Kp;
    params->Ti = Ti;
    params->Td = Td;

    /* Convert to parallel gains */
    pid_ideal_to_parallel(Kp, Ti, Td, &params->Ki, &params->Kd);

    params->N  = 10.0;
    params->b  = 1.0;
    params->c  = 0.0;
    params->Tt = (Ti > 1e-12 && Td > 1e-12) ? sqrt(Ti * Td) : Ti;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Cohen-Coon for SOPDT
 * ────────────────────────────────────────────── */

int cohen_coon_tune_sopdt(const sopdt_model_t *sopdt,
                          cc_controller_type_t type,
                          pid_params_t *params)
{
    if (!sopdt || !params) return -1;

    fopdt_model_t fopdt;
    double T1 = sopdt->T1;
    double T2 = sopdt->T2;

    if (sopdt->use_osc) {
        /* Oscillatory → convert to real poles */
        double zeta = sopdt->zeta;
        double wn   = sopdt->wn;
        if (zeta <= 0) return -1;
        if (zeta > 1.0) {
            double disc = sqrt(zeta * zeta - 1.0);
            T1 = 1.0 / (wn * (zeta - disc));
            T2 = 1.0 / (wn * (zeta + disc));
        } else {
            T1 = 2.0 * zeta / wn;
            T2 = T1 * 0.1;
        }
    }

    /* Skogestad half-rule reduction */
    if (T1 >= T2) {
        fopdt.K = sopdt->K;
        fopdt.T = T1 + 0.5 * T2;
        fopdt.L = sopdt->L + 0.5 * T2;
    } else {
        fopdt.K = sopdt->K;
        fopdt.T = T2 + 0.5 * T1;
        fopdt.L = sopdt->L + 0.5 * T1;
    }

    if (fopdt.L < 1e-12) {
        fopdt.L = fopdt.T * 0.01;
    }

    return cohen_coon_tune(&fopdt, type, params);
}

/* ──────────────────────────────────────────────
 * L5: Cohen-Coon with Custom Decay Ratio
 * ────────────────────────────────────────────── */

int cohen_coon_decay_tune(const fopdt_model_t *model,
                          cc_controller_type_t type,
                          double decay_ratio,
                          pid_params_t *params)
{
    if (cohen_coon_tune(model, type, params) != 0) return -1;

    /**
     * Adjust gain to achieve different decay ratio.
     * Decay ratio DR = amplitude ratio between successive peaks.
     *
     * For a second-order dominant closed loop:
     *   DR = exp(-2πζ / sqrt(1-ζ²))
     *
     * Standard CC: DR = 0.25 → ζ ≈ 0.215
     *
     * Kp_adjusted = Kp * (ζ_desired / ζ_standard)²
     * (since Kp inversely proportional to ζ² for dominant pole)
     */

    if (decay_ratio < 0.01) decay_ratio = 0.01;
    if (decay_ratio > 0.90) decay_ratio = 0.90;

    /* Compute damping ratio from decay ratio */
    double log_dr = log(decay_ratio);
    double zeta_desired = fabs(log_dr) / sqrt(4.0 * M_PI * M_PI + log_dr * log_dr);

    /* Standard ζ ≈ 0.215 from DR = 0.25 */
    double zeta_standard = 0.215;
    double gain_adj = (zeta_desired / zeta_standard) *
                      (zeta_desired / zeta_standard);

    params->Kp *= gain_adj;
    params->Ki *= gain_adj;
    params->Kd *= gain_adj;

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Predict Closed-Loop Response
 * ────────────────────────────────────────────── */

int cohen_coon_predict_response(const fopdt_model_t *model,
                                const pid_params_t *params,
                                double *decay, double *overshoot,
                                double *omega_osc)
{
    if (!model || !params) return -1;

    /**
     * Approximate the closed loop as second-order dominant:
     *   T_cl(s) ≈ ωₙ² / (s² + 2ζωₙs + ωₙ²)
     *
     * From PID + FOPDT, the effective ζ and ωₙ can be approximated
     * using the characteristic equation:
     *   1 + C(s) * G(s) = 0
     *
     * Using Padé approximation for dead time: exp(-Ls) ≈ (1-Ls/2)/(1+Ls/2)
     *
     * Then extract dominant poles to estimate ζ, ωₙ.
     */

    double K = model->K, T = model->T, L = model->L;
    double Kp = params->Kp, Ti = params->Ti, Td = params->Td;

    if (Kp < 1e-12) return -1;

    /* Effective loop gain and time constants */
    double K_loop = Kp * K;
    double tau_I = (Ti > 1e-12) ? Ti : 1e6;
    double L_eff = L;

    /* Approximate second-order dominant poles from
       characteristic: 1 + K_loop * (1 + 1/(s*Ti) + s*Td) * 1/(1+sT) * e^{-sL} ≈ 0
       Using simplified formula (valid for PI, L not too large): */
    double a = T * tau_I;
    double b = tau_I + K_loop * (tau_I * (Td + L_eff) + 0.5*L_eff*L_eff);
    double c = K_loop * (tau_I + L_eff);

    /* For PID, include derivative effect */
    if (Td > 1e-12) {
        a = T * tau_I + K_loop * Td * tau_I;
        b = tau_I + K_loop * (tau_I * L_eff + Td);
        c = K_loop * tau_I;
    }

    if (a < 1e-12) a = 1e-12;

    double zeta = b / (2.0 * sqrt(a * c));
    double wn   = sqrt(c / a);

    if (zeta < 0.0) zeta = 0.01;
    if (zeta > 10.0) zeta = 10.0;

    /* Decay ratio from damping */
    double zeta_eff = (zeta < 1.0) ? zeta : 1.0;
    *decay = exp(-2.0 * M_PI * zeta_eff / sqrt(1.0 - zeta_eff * zeta_eff));
    if (*decay > 1.0) *decay = 1.0;

    /* Overshoot from damping */
    *overshoot = 100.0 * exp(-M_PI * zeta_eff / sqrt(1.0 - zeta_eff * zeta_eff));
    if (*overshoot > 100.0) *overshoot = 100.0;

    /* Oscillation frequency */
    *omega_osc = wn * sqrt(1.0 - zeta_eff * zeta_eff);

    return 0;
}

/* ──────────────────────────────────────────────
 * L5: Cohen-Coon for Integrating Processes
 * ────────────────────────────────────────────── */

int cohen_coon_tune_ipdt(const ipdt_model_t *ipdt,
                         cc_controller_type_t type,
                         pid_params_t *params)
{
    if (!ipdt || !params) return -1;
    if (ipdt->Kv < 1e-12 || ipdt->L < 1e-12) return -1;

    double Kv = ipdt->Kv;
    double L  = ipdt->L;

    /**
     * For integrating processes G(s) = Kv*exp(-L*s)/s:
     *
     * P:    Kp = 1.0 / (Kv * L)   (with offset; rarely used alone)
     * PI:   Kp = 0.9 / (Kv * L);  Ti = 3.33 * L
     * PID:  Kp = 1.2 / (Kv * L);  Ti = 2.0 * L;  Td = 0.5 * L
     */

    switch (type) {
        case CC_CONTROLLER_P:
            params->Kp = 1.0 / (Kv * L);
            params->Ti = 0.0;
            params->Td = 0.0;
            break;
        case CC_CONTROLLER_PI:
            params->Kp = 0.9 / (Kv * L);
            params->Ti = 3.33 * L;
            params->Td = 0.0;
            break;
        case CC_CONTROLLER_PID:
            params->Kp = 1.2 / (Kv * L);
            params->Ti = 2.0 * L;
            params->Td = 0.5 * L;
            break;
        default:
            return -1;
    }

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
 * L5: Cohen-Coon for n-equal-lag Processes
 * ────────────────────────────────────────────── */

int cohen_coon_tune_nlag(double K, double T, double L, int n,
                         cc_controller_type_t type, pid_params_t *params)
{
    if (n < 1) return -1;

    /**
     * High-order process: G(s) = K*exp(-L*s) / (Ts+1)^n
     *
     * Reduce to FOPDT using Skogestad half-rule:
     *   T_eff = T + (n-1)*T/2 = T * (1 + (n-1)/2)
     *   L_eff = L + (n-1)*T/2
     */

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T * (1.0 + (n - 1) * 0.5);
    fopdt.L = L + (n - 1) * T * 0.5;

    return cohen_coon_tune(&fopdt, type, params);
}

/* ──────────────────────────────────────────────
 * L4: Robustness Analysis
 * ────────────────────────────────────────────── */

int cohen_coon_robustness(const fopdt_model_t *nom, const fopdt_model_t *pert,
                          const pid_params_t *params,
                          double *Gm_deg, double *Pm_deg)
{
    if (!nom || !pert || !params || !Gm_deg || !Pm_deg) return -1;

    /**
     * Compute gain and phase margins for both nominal and perturbed models.
     * The degradation quantifies how much the margins shrink due to
     * model mismatch.
     */

    double Gm_nom_lin, Pm_nom_rad;
    double Gm_pert_lin, Pm_pert_rad;

    if (zn_verify_margins(nom, params, &Gm_nom_lin, &Pm_nom_rad) != 0) return -1;
    if (zn_verify_margins(pert, params, &Gm_pert_lin, &Pm_pert_rad) != 0) return -1;

    /* Convert to dB and degrees for reporting */
    double Gm_nom_dB = 20.0 * log10(Gm_nom_lin > 1e-12 ? Gm_nom_lin : 1e-12);
    double Gm_pert_dB = 20.0 * log10(Gm_pert_lin > 1e-12 ? Gm_pert_lin : 1e-12);

    *Gm_deg = Gm_nom_dB - Gm_pert_dB;  /* positive = degradation */
    *Pm_deg = (Pm_nom_rad - Pm_pert_rad) * 180.0 / M_PI;

    return 0;
}
