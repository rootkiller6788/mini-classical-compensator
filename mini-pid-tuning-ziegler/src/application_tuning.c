/**
 * application_tuning.c — Application-Specific PID Tuning Guidelines
 *
 * Implements tuning heuristics and rules of thumb for specific
 * industrial process applications.
 *
 * Knowledge: L6 (canonical process types), L7 (industrial applications:
 *   temperature, flow, pressure, level, pH, cascade, HVAC, DC motor),
 *   L8 (specialized tuning for difficult processes).
 *
 * Reference:
 *   Seborg, Edgar & Mellichamp (2004), "Process Dynamics and Control".
 *   Shinskey, F.G. (1996), "Process Control Systems", McGraw-Hill.
 *   McMillan, G.K. (2005), "Tuning and Control Loop Performance", ISA.
 */

#include "application_tuning.h"
#include "advanced_tuning.h"
#include "ziegler_nichols.h"
#include "cohen_coon.h"
#include "imc_tuning.h"
#include "gain_margin_tuning.h"
#include "fopdt_model.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────
 * L6/L7: Temperature Control Loop Tuning
 * ────────────────────────────────────────────── */

/**
 * Temperature control loops are typically:
 *   - Lag-dominant (T >> L, τ < 0.2)
 *   - Slow (minutes to hours)
 *   - Multi-capacity (high-order, can be reduced to FOPDT)
 *   - Often require PID (derivative beneficial)
 *
 * Typical FOPDT ranges:
 *   K: 0.5-5 °C/%output
 *   T: 1-60 min
 *   L: 0.1-5 min
 *
 * Recommended method: IMC/SIMC with PI or PID.
 *
 * @param K   Static gain.
 * @param T   Time constant [minutes].
 * @param L   Dead time [minutes].
 * @param params [out] Recommended PID parameters.
 * @return    0 on success.
 */
int app_tune_temperature(double K, double T, double L, pid_params_t *params)
{
    if (!params || K < 1e-12 || T < 1e-12) return -1;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T;
    fopdt.L = (L > 0) ? L : T * 0.01;

    /* Temperature loops: use PID if T >> L, otherwise PI */
    int use_pid = (T > 8.0 * fopdt.L) ? 1 : 0;
    double tau_c = fopdt.L * 2.0; /* Conservative for temperature */
    return imc_simc_tune(&fopdt, tau_c, use_pid, params);
}

/* ──────────────────────────────────────────────
 * L7: Flow Control Loop Tuning
 * ────────────────────────────────────────────── */

/**
 * Flow control loops are typically:
 *   - Very fast (seconds)
 *   - Near-integrating or small time constant
 *   - Noisy measurement (requires filtering)
 *   - Fast inner loops in cascade (PI only)
 *
 * Typical FOPDT ranges:
 *   K: 0.5-2 (flow %/%output)
 *   T: 0.1-10 sec
 *   L: 0.05-1 sec
 *
 * Recommended method: PI only (derivative amplifies noise).
 * Gain margin based tuning preferred.
 *
 * @param K   Static gain.
 * @param T   Time constant [seconds].
 * @param L   Dead time [seconds].
 * @param params [out] PI parameters.
 * @return    0 on success.
 */
int app_tune_flow(double K, double T, double L, pid_params_t *params)
{
    if (!params || K < 1e-12) return -1;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = (T > 0) ? T : 1.0;    /* Flow is very fast, default T = 1s */
    fopdt.L = (L > 0) ? L : 0.1;

    /* Flow loops: PI only, aggressive gain */
    double Am = 3.0;                /* Gain margin = 3 */
    double Pm = 45.0 * M_PI / 180.0; /* Phase margin = 45° */
    return gm_pm_tune_pi(&fopdt, Am, Pm, params);
}

/* ──────────────────────────────────────────────
 * L7: Pressure Control Loop Tuning
 * ────────────────────────────────────────────── */

/**
 * Pressure control loops are typically:
 *   - Moderate speed (seconds to minutes)
 *   - Often self-regulating (FOPDT)
 *   - Can be noisy if gas
 *
 * Typical FOPDT ranges:
 *   K: 0.5-5
 *   T: 1-60 sec
 *   L: 0.1-10 sec
 *
 * Recommended: Z-N or Cohen-Coon for moderate dead time.
 *
 * @param K   Static gain.
 * @param T   Time constant [seconds].
 * @param L   Dead time [seconds].
 * @param params [out] PID parameters.
 * @return    0 on success.
 */
int app_tune_pressure(double K, double T, double L, pid_params_t *params)
{
    if (!params || K < 1e-12 || T < 1e-12) return -1;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T;
    fopdt.L = (L > 0) ? L : T * 0.05;

    /* Pressure: PID if moderate dead time */
    if (fopdt_dead_time_ratio(&fopdt) > 0.3) {
        return cohen_coon_tune(&fopdt, CC_CONTROLLER_PID, params);
    } else {
        return zn_step_tune(&fopdt, ZN_CONTROLLER_PID, params);
    }
}

/* ──────────────────────────────────────────────
 * L7: Level Control Loop Tuning
 * ────────────────────────────────────────────── */

/**
 * Level control loops are typically:
 *   - Integrating process (no self-regulation)
 *   - Surge tanks: loose control (P-only acceptable)
 *   - Tight level: PI required
 *
 * Model: G(s) = Kv * exp(-L*s) / s
 *   Kv = 1 / (cross-sectional area)
 *
 * For surge tanks, P-only with Kp = 0.5-2 is typical.
 * For tight level control, PI with moderate gain.
 *
 * @param Kv  Velocity gain [1/time].
 * @param L   Dead time [time].
 * @param surge_tank 1 = surge (P-only), 0 = tight (PI).
 * @param params [out] P or PI parameters.
 * @return    0 on success.
 */
int app_tune_level(double Kv, double L, int surge_tank, pid_params_t *params)
{
    if (!params || Kv < 1e-12) return -1;
    if (L < 1e-12) L = 0.1;

    memset(params, 0, sizeof(*params));

    if (surge_tank) {
        /* P-only: accept offset, avoid interaction with other loops */
        params->Kp = 1.0 / (Kv * L);
        params->Ki = 0.0;
        params->Kd = 0.0;
    } else {
        /* PI for tight level: use IMC for integrating */
        ipdt_model_t ipdt;
        ipdt.Kv = Kv;
        ipdt.L  = L;
        return cohen_coon_tune_ipdt(&ipdt, CC_CONTROLLER_PI, params);
    }

    params->N  = 10.0;
    params->b  = 1.0;
    params->c  = 0.0;
    params->Tt = params->Ti;

    return 0;
}

/* ──────────────────────────────────────────────
 * L7: pH Control Loop Tuning
 * ────────────────────────────────────────────── */

/**
 * pH control is extremely nonlinear (S-shaped titration curve).
 * Gain varies by 100-1000x across the operating range.
 *
 * Strategy:
 *   - Use gain scheduling based on pH setpoint
 *   - Conservative tuning near neutral (pH 7, highest gain)
 *   - Aggressive tuning far from neutral (low gain)
 *
 * @param pH_sp      pH setpoint (0-14).
 * @param valve_size Actuator range [ml/s].
 * @param max_gain   Maximum process gain at steepest slope.
 * @param params     [out] PI parameters.
 * @return    0 on success.
 */
int app_tune_ph(double pH_sp, double valve_size, double max_gain,
                pid_params_t *params)
{
    if (!params) return -1;

    (void)valve_size; /* Documented: actuator capacity, used for output scaling */

    /**
     * The gain at a given pH depends on the titration curve.
     * Near pH 7, the gain is maximum (buffer strength is minimal).
     *
     *   K_eff = max_gain * f(pH_sp)
     *   f(pH) ≈ exp(-(pH-7)² / (2σ²))  with σ ≈ 2
     *
     * Use IMC-based tuning with λ = 5*L for robustness.
     */

    double sigma = 2.0;
    double dpH = pH_sp - 7.0;
    double f_pH = exp(-dpH * dpH / (2.0 * sigma * sigma));
    double K_eff = max_gain * f_pH;
    if (K_eff < max_gain * 0.01) K_eff = max_gain * 0.01;

    /* Assume typical mixing delay */
    double T_eff = 5.0;  /* minutes, typical mixing lag */
    double L_eff = 1.0;  /* minutes, typical transport delay */

    fopdt_model_t fopdt;
    fopdt.K = K_eff;
    fopdt.T = T_eff;
    fopdt.L = L_eff;

    /* Robust IMC for pH: λ = 5*L */
    double lambda = 5.0 * L_eff;
    return imc_simc_tune(&fopdt, lambda, 0, params);
}

/* ──────────────────────────────────────────────
 * L7: DC Motor Speed Control Tuning
 * ────────────────────────────────────────────── */

/**
 * DC motor speed control:
 *   G(s) = K / (τ_m * s + 1)(τ_e * s + 1)
 *
 * where τ_m = mechanical time constant, τ_e = electrical time constant.
 * Typically τ_m >> τ_e, so SOPDT → FOPDT reduction is valid.
 *
 * @param K_motor  Motor gain [rpm/V].
 * @param tau_m    Mechanical time constant [sec].
 * @param tau_e    Electrical time constant [sec].
 * @param params   [out] PID parameters.
 * @return         0 on success.
 */
int app_tune_dc_motor(double K_motor, double tau_m, double tau_e,
                      pid_params_t *params)
{
    if (!params || K_motor < 1e-12) return -1;
    if (tau_m < 1e-12) tau_m = 1.0;
    if (tau_e < 1e-12) tau_e = 0.1;

    /* Reduce to FOPDT: T ≈ τ_m + τ_e/2, L ≈ τ_e/2 */
    fopdt_model_t fopdt;
    fopdt.K = K_motor;
    fopdt.T = tau_m + 0.5 * tau_e;
    fopdt.L = 0.5 * tau_e;

    /* Motor control: IMC with moderate speed */
    double tau_c = tau_e; /* desired closed-loop = electrical time scale */
    return imc_simc_tune(&fopdt, tau_c, 1, params);
}

/* ──────────────────────────────────────────────
 * L7: HVAC Temperature Control
 * ────────────────────────────────────────────── */

/**
 * HVAC zone temperature control:
 *   - Very slow (10-60 minute time constant)
 *   - Significant dead time (air transport)
 *   - External disturbances (solar, occupancy)
 *
 * @param K      Zone gain [°C/%damper].
 * @param T      Zone time constant [minutes].
 * @param L      Transport delay [minutes].
 * @param params [out] PI parameters (PID rarely needed in HVAC).
 * @return       0 on success.
 */
int app_tune_hvac(double K, double T, double L, pid_params_t *params)
{
    if (!params || K < 1e-12 || T < 1e-12) return -1;
    if (L < 1e-12) L = T * 0.05;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T;
    fopdt.L = L;

    /* HVAC: conservative IMC, PI only */
    double tau_c = 3.0 * L; /* Very conservative */
    return imc_simc_tune(&fopdt, tau_c, 0, params);
}

/* ──────────────────────────────────────────────
 * L7: Chemical Reactor Temperature Control
 * ────────────────────────────────────────────── */

/**
 * Exothermic CSTR temperature control:
 *   - Highly nonlinear (Arrhenius kinetics)
 *   - Can be unstable (thermal runaway)
 *   - Requires conservative tuning with gain margin ≥ 3
 *
 * @param K      Reactor gain [°C/%cooling].
 * @param T      Dominant time constant [min].
 * @param L      Effective dead time [min].
 * @param params [out] PID parameters.
 * @return       0 on success.
 */
int app_tune_chemical_reactor(double K, double T, double L,
                              pid_params_t *params)
{
    if (!params || K < 1e-12 || T < 1e-12) return -1;
    if (L < 1e-12) L = T * 0.1;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T;
    fopdt.L = L;

    /* Chemical reactor: need gain margin ≥ 3, phase margin ≥ 60° */
    double Am = 3.5;
    double Pm = 60.0 * M_PI / 180.0;
    return gm_pm_tune_pid(&fopdt, Am, Pm, 1, params);
}

/* ──────────────────────────────────────────────
 * L7: Paper Machine Headbox Control
 * ────────────────────────────────────────────── */

/**
 * Paper machine headbox pressure control:
 *   - FOPDT with moderate dynamics
 *   - Coupled with level control
 *   - PI control typical, moderate tuning
 *
 * @param K      Process gain.
 * @param T      Time constant [sec].
 * @param L      Dead time [sec].
 * @param params [out] PI parameters.
 * @return       0 on success.
 */
int app_tune_paper_headbox(double K, double T, double L,
                           pid_params_t *params)
{
    if (!params || K < 1e-12 || T < 1e-12) return -1;
    if (L < 1e-12) L = T * 0.1;

    fopdt_model_t fopdt;
    fopdt.K = K;
    fopdt.T = T;
    fopdt.L = L;

    /* Moderate tuning: Cohen-Coon PI */
    return cohen_coon_tune(&fopdt, CC_CONTROLLER_PI, params);
}

/* ──────────────────────────────────────────────
 * L8: Time-Varying Process Advisory
 * ────────────────────────────────────────────── */

/**
 * For processes where parameters vary significantly, this function
 * provides a recommendation on whether to use adaptive PID,
 * gain scheduling, or robust fixed-parameter tuning.
 *
 * @param K_variation   Fractional variation in gain (σ_K / K_mean).
 * @param L_variation   Fractional variation in dead time.
 * @param T_variation   Fractional variation in time constant.
 * @return              0 = fixed tuning OK,
 *                      1 = gain scheduling recommended,
 *                      2 = adaptive tuning recommended.
 */
int app_recommend_adaptation(double K_variation, double L_variation,
                             double T_variation)
{
    double total_variation = K_variation + L_variation + T_variation;

    if (total_variation < 0.3) {
        return 0; /* Fixed robust tuning is sufficient */
    } else if (total_variation < 0.8) {
        return 1; /* Gain scheduling recommended */
    } else {
        return 2; /* Full adaptive control recommended */
    }
}

/**
 * Compute recommended tuning method for a FOPDT process based
 * on the normalized dead time ratio and process type.
 *
 * Decision logic (from Seborg et al. 2004, McMillan 2005):
 *   τ < 0.1 (lag-dominant):  IMC/SIMC — most modern
 *   0.1 ≤ τ < 0.3:            Z-N step response — classic
 *   0.3 ≤ τ < 0.5:            Cohen-Coon — handles delay better
 *   τ ≥ 0.5:                   Gain margin based — delay-dominant
 *
 * @param model FOPDT model.
 * @return      Recommended tuning method enum.
 */
pid_tune_method_t app_recommend_method(const fopdt_model_t *model)
{
    if (!model) return TUNE_METHOD_ZN_STEP;

    double tau = fopdt_dead_time_ratio(model);

    if (tau < 0.1) {
        return TUNE_METHOD_IMC;
    } else if (tau < 0.3) {
        return TUNE_METHOD_ZN_STEP;
    } else if (tau < 0.5) {
        return TUNE_METHOD_COHEN_COON;
    } else {
        return TUNE_METHOD_GM_PM;
    }
}

/**
 * Comprehensive application-driven auto-tuning: given minimal
 * process information (K, T, L and application type), selects
 * the best tuning method and returns PID parameters.
 *
 * @param app_type  0=temperature, 1=flow, 2=pressure, 3=level, 4=pH.
 * @param K         Process gain.
 * @param T         Time constant.
 * @param L         Dead time.
 * @param params    [out] PID parameters.
 * @return          0 on success.
 */
int app_auto_tune(int app_type, double K, double T, double L,
                  pid_params_t *params)
{
    if (!params) return -1;

    switch (app_type) {
        case 0: return app_tune_temperature(K, T, L, params);
        case 1: return app_tune_flow(K, T, L, params);
        case 2: return app_tune_pressure(K, T, L, params);
        case 3: return app_tune_level(K, L, 0, params);
        case 4: return app_tune_ph(7.0, 100.0, K, params);
        case 5: return app_tune_dc_motor(K, T, L, params);
        case 6: return app_tune_hvac(K, T, L, params);
        case 7: return app_tune_chemical_reactor(K, T, L, params);
        case 8: return app_tune_paper_headbox(K, T, L, params);
        default:
            /* Default: use recommendation engine */
            {
                fopdt_model_t fopdt;
                fopdt.K = K;
                fopdt.T = (T > 1e-12) ? T : 1.0;
                fopdt.L = (L > 0) ? L : fopdt.T * 0.05;
                pid_tune_method_t method = app_recommend_method(&fopdt);
                switch (method) {
                    case TUNE_METHOD_IMC:
                        return imc_simc_tune(&fopdt, fopdt.L, 1, params);
                    case TUNE_METHOD_COHEN_COON:
                        return cohen_coon_tune(&fopdt, CC_CONTROLLER_PID, params);
                    case TUNE_METHOD_GM_PM: {
                        double Am = 3.0;
                        double Pm = 45.0 * M_PI / 180.0;
                        return gm_pm_tune_pid(&fopdt, Am, Pm, 0, params);
                    }
                    default:
                        return zn_step_tune(&fopdt, ZN_CONTROLLER_PID, params);
                }
            }
    }
}

/**
 * Estimate FOPDT model parameters from simple process tests
 * (bump test / process reaction curve).
 *
 * Given a few summary statistics from a bump test:
 *   y_initial, y_final, y_max_slope, time_to_max_slope
 * estimate K, T, L for quick-and-dirty tuning.
 *
 * @param y_initial      Initial output.
 * @param y_final        Final output (steady state).
 * @param max_slope      Maximum slope during transient.
 * @param t_at_max_slope Time at which max slope occurs.
 * @param step_mag       Input step magnitude.
 * @param model          [out] Estimated FOPDT model.
 * @return               0 on success.
 */
int app_estimate_fopdt(double y_initial, double y_final,
                       double max_slope, double t_at_max_slope,
                       double step_mag, fopdt_model_t *model)
{
    if (!model || step_mag < 1e-12 || max_slope < 1e-12) return -1;

    double dy = y_final - y_initial;

    model->K = dy / step_mag;
    model->L = t_at_max_slope - (y_final - y_initial) / (2.0 * max_slope);
    if (model->L < 0.0) model->L = 0.0;
    model->T = dy / max_slope - model->L;
    if (model->T < 1e-12) model->T = 1e-6;

    return 0;
}
