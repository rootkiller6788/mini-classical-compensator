#ifndef IMC_TUNING_H
#define IMC_TUNING_H

/**
 * imc_tuning.h — Internal Model Control (IMC) Based PID Tuning
 *
 * Knowledge coverage:
 *   L1: IMC structure — model-based controller design paradigm
 *   L2: IMC filter design — robustness vs. performance tradeoff via λ
 *   L4: IMC-PID equivalence — analytical derivation for common process models
 *   L5: IMC-based PID computation for FOPDT, SOPDT, integrating processes
 *   L6: Application to chemical process control (Morari & Zafiriou, 1989)
 *   L8: Advanced: IMC with Smith predictor for dominant dead-time
 *
 * Reference:
 *   Garcia, C.E. & Morari, M. (1982)
 *   "Internal Model Control. 1. A Unifying Review and Some New Results",
 *   Ind. Eng. Chem. Process Des. Dev., 21, 308-323.
 *
 *   Rivera, D.E., Morari, M. & Skogestad, S. (1986)
 *   "Internal Model Control. 4. PID Controller Design",
 *   Ind. Eng. Chem. Process Des. Dev., 25, 252-265.
 *
 *   Skogestad, S. (2003)
 *   "Simple Analytic Rules for Model Reduction and PID Controller Tuning",
 *   Journal of Process Control, 13, 291-309.  [SIMC rules]
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: IMC Structure and Filter Definitions
 * ────────────────────────────────────────────── */

/**
 * IMC filter order: determines the roll-off rate in the IMC controller.
 *
 * 1st-order:  F(s) = 1 / (λ*s + 1)
 * 2nd-order:  F(s) = 1 / (λ*s + 1)²
 */
typedef enum {
    IMC_FILTER_FIRST  = 1,
    IMC_FILTER_SECOND = 2
} imc_filter_order_t;

/**
 * IMC tuning rule set (Rivera-Morari-Skogestad, 1986).
 *
 * Different process models yield different IMC-PID equivalences.
 */
typedef enum {
    IMC_MODEL_FOPDT = 0,         /* First-order + dead time */
    IMC_MODEL_SOPDT_OVERDAMPED,  /* Second-order overdamped */
    IMC_MODEL_SOPDT_UNDERDAMPED, /* Second-order underdamped */
    IMC_MODEL_INTEGRATING,       /* Integrating process */
    IMC_MODEL_INTEGRATING_LAG,   /* Integrating + lag */
    IMC_MODEL_UNSTABLE_FOPDT     /* Unstable first-order + dead time */
} imc_model_type_t;

/**
 * IMC-PID tuning parameters container.
 *
 * Stores the IMC-specific interpretation of PID gains
 * with explicit λ (closed-loop time constant) parameter.
 *
 * IMC design steps:
 *   1. Factor plant model: G(s) = G_+(s) * G_-(s)
 *      where G_-(s) is invertible (stable, minimum-phase)
 *            G_+(s) contains non-invertible parts (dead time, RHP zeros)
 *   2. Design IMC controller: Q(s) = G_-(s)⁻¹ * F(s)
 *      where F(s) is the IMC filter.
 *   3. Convert to standard feedback form:
 *        C(s) = Q(s) / (1 - G(s)*Q(s))
 *   4. Match to PID form: C_PID(s) = Kp + Ki/s + Kd*s
 */
typedef struct {
    double lambda;  /* IMC filter time constant (closed-loop speed) */
    double Kp;      /* Equivalent PID proportional gain */
    double Ki;      /* Equivalent PID integral gain */
    double Kd;      /* Equivalent PID derivative gain */
    imc_filter_order_t filter_order;
    imc_model_type_t   model_type;
} imc_tuning_result_t;

/* ──────────────────────────────────────────────
 * L4: IMC-PID Equivalence Derivations
 * ────────────────────────────────────────────── */

/**
 * IMC-PID tuning for FOPDT model using 1st-order Padé approximation
 * of the dead time term:  exp(-L*s) ≈ (1 - L*s/2) / (1 + L*s/2)
 *
 * After algebraic manipulation, the equivalent PID parameters are:
 *
 *   Kp = (T + L/2) / (K * (λ + L/2))
 *   Ti = T + L/2
 *   Td = (T * L) / (2*T + L)
 *
 * The IMC filter time constant λ determines closed-loop speed:
 *   λ small → fast response, less robust
 *   λ large → slow response, more robust
 *   λ ≥ 1.7*L recommended for stability (Morari & Zafiriou, 1989)
 *
 * @param model   FOPDT model.
 * @param lambda  Desired closed-loop time constant.
 * @param filter_order Filter order (1st or 2nd).
 * @param result  [out] IMC-PID tuning result.
 * @return        0 on success.
 */
int imc_tune_fopdt(const fopdt_model_t *model, double lambda,
                   imc_filter_order_t filter_order,
                   imc_tuning_result_t *result);

/**
 * IMC-PID tuning for overdamped SOPDT model.
 *
 * Process: G(s) = K * exp(-L*s) / ((τ₁*s + 1)*(τ₂*s + 1))
 *
 * Equivalent PID (Rivera-Morari-Skogestad, 1986):
 *   Kp = (τ₁ + τ₂) / (K*(λ + L))
 *   Ti = τ₁ + τ₂
 *   Td = τ₁*τ₂ / (τ₁ + τ₂)
 *
 * @param sopdt   SOPDT model (use_t1t2 mode).
 * @param lambda  Closed-loop time constant.
 * @param filter_order Filter order.
 * @param result  [out] IMC-PID tuning result.
 * @return        0 on success.
 */
int imc_tune_sopdt(const sopdt_model_t *sopdt, double lambda,
                   imc_filter_order_t filter_order,
                   imc_tuning_result_t *result);

/**
 * IMC-PID tuning for integrating process: G(s) = Kv * exp(-L*s) / s
 *
 * Equivalent PID:
 *   Kp = (2*λ + L) / (Kv * (λ + L)²)
 *   Ti = 2*λ + L
 *   Td = (λ² + λ*L) / (2*λ + L)
 *
 * @param ipdt    Integrating process model.
 * @param lambda  Closed-loop time constant.
 * @param filter_order Filter order.
 * @param result  [out] IMC-PID tuning result.
 * @return        0 on success.
 */
int imc_tune_integrating(const ipdt_model_t *ipdt, double lambda,
                         imc_filter_order_t filter_order,
                         imc_tuning_result_t *result);

/* ──────────────────────────────────────────────
 * L5: IMC Filter Design and Robustness
 * ────────────────────────────────────────────── */

/**
 * Skogestad IMC (SIMC) tuning rules for FOPDT.
 *
 * The SIMC rules (Skogestad, 2003) are a simplified version that
 * chooses λ = L (the dead time) as default, providing a good balance
 * between performance and robustness for most processes.
 *
 * PI controller (recommended for most cases):
 *   Kp = T / (K * (L + τ_c))    where τ_c = L (default)
 *   Ti = min(T, 4*(L + τ_c))
 *
 * PID controller (if derivative is beneficial, T > 8*L):
 *   Kp = T / (K * (L + τ_c))
 *   Ti = min(T, 4*(L + τ_c))
 *   Td = T / 8    (simplified)
 *
 * @param model     FOPDT model.
 * @param tau_c     Desired closed-loop time constant (τ_c).
 *                  Set τ_c = L for default robust tuning.
 * @param use_pid   1 = PID, 0 = PI only.
 * @param params    [out] PID parameters (IDEAL form).
 * @return          0 on success.
 */
int imc_simc_tune(const fopdt_model_t *model, double tau_c, int use_pid,
                  pid_params_t *params);

/**
 * IMC filter design: compute λ from desired robustness (Ms).
 *
 * The IMC filter λ directly controls robustness.
 * An approximate relationship (Morari & Zafiriou):
 *   Ms ≈ 1 + L / λ
 *
 * This function inverts the relationship:
 *   λ ≈ L / (Ms - 1)
 *
 * Typical Ms values: 1.2 (robust), 1.4 (balanced), 2.0 (aggressive).
 *
 * @param L   Dead time.
 * @param Ms  Desired maximum sensitivity (1.2 to 2.0).
 * @return    Recommended λ.
 */
double imc_lambda_from_Ms(double L, double Ms);

/**
 * IMC robustness analysis: given λ and model uncertainty,
 * compute the guaranteed stability margin.
 *
 * For multiplicative uncertainty:
 *   |G_actual(jω) - G_model(jω)| / |G_model(jω)| ≤ l_m(ω)
 *
 * Robust stability condition:
 *   |F(jω) * T(jω)| * l_m(ω) < 1  for all ω
 *
 * @param model      Nominal FOPDT model.
 * @param lambda     IMC filter time constant.
 * @param delta_K    Fractional gain uncertainty (0 ≤ δ < 1).
 * @param delta_L    Fractional dead-time uncertainty.
 * @param delta_T    Fractional time-constant uncertainty.
 * @return           Upper bound of complementary sensitivity (worst-case).
 */
double imc_robustness_analysis(const fopdt_model_t *model, double lambda,
                               double delta_K, double delta_L, double delta_T);

/**
 * IMC-based PID for unstable FOPDT process.
 *
 * Process: G(s) = K * exp(-L*s) / (T*s - 1)
 *
 * For unstable processes, the IMC filter must be at least 2nd order
 * to ensure properness of the equivalent PID controller.
 *
 * @param model   FOPDT with negative T (unstable pole).
 * @param lambda  Filter time constant (λ > L/2 required).
 * @param filter_order Must be IMC_FILTER_SECOND for unstable.
 * @param result  [out] IMC-PID tuning result.
 * @return        0 on success.
 */
int imc_tune_unstable(const fopdt_model_t *model, double lambda,
                      imc_filter_order_t filter_order,
                      imc_tuning_result_t *result);

/**
 * Convert IMC tuning result to standard PID parameters.
 *
 * Handles the conversion from IMC structure to PID controller form,
 * including the derivative filter selection.
 *
 * @param result  IMC tuning result.
 * @param params  [out] Standard PID parameters.
 * @param N       [out] Recommended derivative filter constant.
 * @return        0 on success.
 */
int imc_to_pid_params(const imc_tuning_result_t *result,
                      pid_params_t *params, double *N);

#ifdef __cplusplus
}
#endif

#endif /* IMC_TUNING_H */
