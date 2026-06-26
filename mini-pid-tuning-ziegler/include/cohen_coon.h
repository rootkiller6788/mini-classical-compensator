#ifndef COHEN_COON_H
#define COHEN_COON_H

/**
 * cohen_coon.h — Cohen-Coon PID Tuning Method (1953)
 *
 * Knowledge coverage:
 *   L1: Cohen-Coon tuning definition — closed-loop pole placement for FOPDT
 *   L2: Quarter-decay ratio criterion (alternative to amplitude ratio)
 *   L4: Cohen-Coon tuning rules and their derivation from pole-placement
 *   L5: Algorithmic implementation and numerical validation
 *   L6: Application to FOPDT processes with significant dead time
 *
 * Reference:
 *   Cohen, G.H. & Coon, G.A. (1953)
 *   "Theoretical Considerations of Retarded Control",
 *   Trans. ASME, 75, 827-834.
 *
 * Key insight: Cohen-Coon explicitly considers the dead-time-to-time-constant
 * ratio (L/T) in the tuning formulas, making them more accurate than Z-N for
 * processes with significant dead time (L/T > 0.5).
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: Cohen-Coon Method Definition
 * ────────────────────────────────────────────── */

/**
 * Cohen-Coon controller types (same categories as Z-N).
 */
typedef enum {
    CC_CONTROLLER_P  = 0,
    CC_CONTROLLER_PI = 1,
    CC_CONTROLLER_PD = 2,
    CC_CONTROLLER_PID= 3
} cc_controller_type_t;

/**
 * Cohen-Coon tuning equations (L4).
 *
 * Based on FOPDT model G(s) = K*exp(-L*s)/(T*s+1), define:
 *   μ = L / T   (normalized dead time ratio)
 *
 * The Cohen-Coon rules are:
 *
 * P:
 *   Kp = (T/KL) * (1 + 0.333*μ)       = (1/(K*μ)) * (1 + 0.333*μ)
 *
 * PI:
 *   Kp = (T/KL) * (0.9 + 0.083*μ)     = (0.9/(K*μ)) * (1 + 0.092*μ)
 *   Ti = L * (3.33 + 0.3*μ) / (1 + 2.2*μ)
 *
 * PID:
 *   Kp = (T/KL) * (1.35 + 0.27*μ)     = (1.35/(K*μ)) * (1 + 0.200*μ)
 *   Ti = L * (2.50 + 0.46*μ) / (1 + 0.61*μ)
 *   Td = L * (0.37) / (1 + 0.19*μ)
 */

/* ──────────────────────────────────────────────
 * L5: Cohen-Coon Tuning Functions
 * ────────────────────────────────────────────── */

/**
 * Compute Cohen-Coon PID tuning parameters for a FOPDT model.
 *
 * Complexity: O(1).
 *
 * @param model  Process model (K, T, L).
 * @param type   Controller type (P/PI/PD/PID).
 * @param params [out] Computed PID parameters (IDEAL form).
 * @return       0 on success, -1 if model invalid (K≤0, T≤0, L≤0).
 */
int cohen_coon_tune(const fopdt_model_t *model, cc_controller_type_t type,
                    pid_params_t *params);

/**
 * Compute Cohen-Coon tuning for SOPDT model.
 *
 * Uses the dominant-pole approximation: transforms SOPDT → FOPDT
 * then applies Cohen-Coon rules. The secondary pole is accounted for
 * by adjusting the effective dead time.
 *
 * @param sopdt  Second-order-plus-dead-time model.
 * @param type   Controller type.
 * @param params [out] PID parameters.
 * @return       0 on success.
 */
int cohen_coon_tune_sopdt(const sopdt_model_t *sopdt,
                          cc_controller_type_t type,
                          pid_params_t *params);

/**
 * Compute Cohen-Coon with a user-specified decay ratio.
 *
 * Original Cohen-Coon is designed for 1/4 decay ratio. This function
 * adjusts the gain to achieve a different decay ratio target.
 *   Decay ratio = (amplitude of 2nd peak)/(amplitude of 1st peak)
 *
 * @param model       FOPDT model.
 * @param type        Controller type.
 * @param decay_ratio Desired decay ratio (0.0 = no overshoot, 0.25 = standard).
 * @param params      [out] PID parameters.
 * @return            0 on success.
 */
int cohen_coon_decay_tune(const fopdt_model_t *model,
                          cc_controller_type_t type,
                          double decay_ratio,
                          pid_params_t *params);

/**
 * Cohen-Coon evaluation: predict closed-loop decay ratio and overshoot
 * for a given FOPDT + PID tuning analytically.
 *
 * Uses the dominant pole approximation of the closed-loop.
 * Complexity: O(1).
 *
 * @param model    FOPDT model.
 * @param params   PID parameters.
 * @param decay    [out] Predicted decay ratio.
 * @param overshoot [out] Predicted % overshoot.
 * @param omega_osc [out] Predicted oscillation frequency.
 * @return         0 on success.
 */
int cohen_coon_predict_response(const fopdt_model_t *model,
                                const pid_params_t *params,
                                double *decay, double *overshoot,
                                double *omega_osc);

/**
 * Compute Cohen-Coon for integrating processes (IPDT).
 *
 * For G(s) = Kv*exp(-L*s)/s, Cohen and Coon provided modified rules:
 *
 * P:   Kp = 1.0 / (Kv*L) * (1 - 0.06)
 * PI:  Kp = 0.9 / (Kv*L);    Ti = 3.33*L
 * PID: Kp = 1.2 / (Kv*L);    Ti = 2.0*L;  Td = 0.5*L
 *
 * @param ipdt   Integrating process model.
 * @param type   Controller type.
 * @param params [out] PID parameters.
 * @return       0 on success.
 */
int cohen_coon_tune_ipdt(const ipdt_model_t *ipdt,
                         cc_controller_type_t type,
                         pid_params_t *params);

/**
 * High-order Cohen-Coon: for processes with multiple lags.
 *
 * When the process is high-order, approximating as FOPDT can be poor.
 * This function uses a higher-order extension:
 *
 *   G(s) = K * exp(-L*s) / (T*s+1)^n
 *
 * The effective FOPDT parameters are derived from the n-lag process
 * using the Skogestad half-rule approximation.
 *
 * @param K        Static gain.
 * @param T        Individual lag time constant.
 * @param L        Transport delay.
 * @param n        Number of equal lags (n ≥ 1).
 * @param type     Controller type.
 * @param params   [out] PID parameters.
 * @return         0 on success.
 */
int cohen_coon_tune_nlag(double K, double T, double L, int n,
                         cc_controller_type_t type, pid_params_t *params);

/**
 * Perturbation analysis: how sensitive is Cohen-Coon tuning to model errors?
 *
 * Given nominal and perturbed FOPDT models, evaluate the change in
 * closed-loop stability margins.
 *
 * @param nom      Nominal FOPDT model.
 * @param pert     Perturbed (actual) FOPDT model.
 * @param params   PID tuning (computed from nominal).
 * @param Gm_deg   [out] Gain margin degradation [dB].
 * @param Pm_deg   [out] Phase margin degradation [degrees].
 * @return         0 on success.
 */
int cohen_coon_robustness(const fopdt_model_t *nom, const fopdt_model_t *pert,
                          const pid_params_t *params,
                          double *Gm_deg, double *Pm_deg);

#ifdef __cplusplus
}
#endif

#endif /* COHEN_COON_H */
