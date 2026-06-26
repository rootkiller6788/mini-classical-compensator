#ifndef FOPDT_MODEL_H
#define FOPDT_MODEL_H

/**
 * fopdt_model.h — First-Order Plus Dead Time (FOPDT) Model Identification
 *
 * Knowledge coverage:
 *   L1: FOPDT model definition — the lingua franca of process control
 *   L2: Model identification concept — from plant data to parametric model
 *   L3: Laplace domain: G(s) = K*exp(-L*s)/(T*s+1); time domain step response
 *   L4: Identification theory — bias/variance tradeoff, identifiability
 *   L5: Computational methods: tangent, two-point, area, least-squares
 *   L6: Standard identification on canonical processes
 *
 * Reference:
 *   Seborg, D.E., Edgar, T.F. & Mellichamp, D.A. (2004)
 *   "Process Dynamics and Control", 2nd ed., Wiley, Ch.7.
 *   Åström, K.J. & Hägglund, T. (1995) Ch.2 "Process Models".
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: FOPDT Model Structures (forward declared in pid_tuning.h)
 * ────────────────────────────────────────────── */

/**
 * FOPDT identification method enumeration.
 *
 * Each method has different statistical properties and computational
 * complexity. The choice depends on data quality and application.
 */
typedef enum {
    FOPDT_ID_GRAPHICAL = 0,    /* Tangent line at inflection point (Z-N classic) */
    FOPDT_ID_TWO_POINT,        /* Two-point method (28.3% and 63.2% levels) */
    FOPDT_ID_AREA,             /* Area method (integrals of step response) */
    FOPDT_ID_LEAST_SQUARES,    /* Least-squares parameter estimation */
    FOPDT_ID_PREDICTION_ERROR, /* Prediction error minimization (PEM) */
    FOPDT_ID_WEIGHTED_LS,      /* Weighted least squares for noise filtering */
    FOPDT_ID_SKOGESTAD_HALF_RULE, /* Skogestad's half-rule for model reduction */
    FOPDT_ID_SUNDARESAN_KRISHNASWAMY /* Sundaresan-Krishnaswamy (1977) */
} fopdt_id_method_t;

/**
 * L1: FOPDT identification result with quality metrics.
 */
typedef struct {
    fopdt_model_t model;       /* Identified model parameters */
    double r_squared;          /* Coefficient of determination R² */
    double rmse;               /* Root Mean Square Error */
    double fit_pct;            /* NRMSE fit percent (100 = perfect) */
    int    iterations;         /* Number of iterations (for LS/PEM) */
    int    converged;          /* 1 = converged */
    fopdt_id_method_t method;  /* Method used for identification */
} fopdt_id_result_t;

/* ──────────────────────────────────────────────
 * L5: FOPDT Model Identification Algorithms
 * ────────────────────────────────────────────── */

/**
 * **Tangent (graphical) method** — the classical Ziegler-Nichols approach.
 *
 * Algorithm:
 *   1. Find inflection point (maximum slope) in step response y(t).
 *   2. Draw tangent line through inflection point.
 *   3. L = intersection of tangent with baseline (time of crossing).
 *   4. T = (y_ss - intersection_y) / (max slope).
 *   5. K = Δy_ss / Δu (steady-state gain).
 *
 * Complexity: O(N) data scan. N = number of data points.
 * Strengths: Simple, intuitive. Weaknesses: Sensitive to noise, single-point.
 *
 * @param data   Step response data.
 * @param result [out] Identified model with quality metrics.
 * @return       0 on success, -1 if data insufficient.
 */
int fopdt_identify_graphical(const step_response_data_t *data,
                             fopdt_id_result_t *result);

/**
 * **Two-point method** — uses times at 28.3% and 63.2% of final value.
 *
 * Theoretical basis: For a pure FOPDT step response:
 *   y(t) = K * [1 - exp(-(t-L)/T)]  (for t ≥ L)
 *   At y = 0.283*K:  t28 = L + T/3
 *   At y = 0.632*K:  t63 = L + T
 *
 * Solving:
 *   T = 1.5 * (t63 - t28)
 *   L = t63 - T
 *
 * Advantage: More robust to measurement noise than tangent method.
 *
 * @param data   Step response data.
 * @param result [out] Identified model.
 * @return       0 on success.
 */
int fopdt_identify_two_point(const step_response_data_t *data,
                             fopdt_id_result_t *result);

/**
 * **Area method** — uses integrals of the step response.
 *
 * Define the normalized step response: ȳ(t) = y(t) / y_ss.
 * Compute:
 *   A0 = ∫₀^∞ [1 - ȳ(t)] dt = T + L
 *   A1 = ∫₀^A0 ȳ(t) dt  = T / e  (for FOPDT)
 *
 * From which:
 *   T = e * A1   (where e = exp(1))
 *   L = A0 - T
 *
 * This method averages over the entire response, making it robust to
 * noise but sensitive to steady-state estimation error.
 *
 * Reference: Åström & Hägglund (1995).
 *
 * @param data   Step response data.
 * @param result [out] Identified model.
 * @return       0 on success.
 */
int fopdt_identify_area(const step_response_data_t *data,
                        fopdt_id_result_t *result);

/**
 * **Least-squares parameter estimation** for FOPDT model.
 *
 * Minimizes:  J(θ) = Σ [y_measured(t_i) - y_model(t_i; θ)]²
 * where θ = [K, T, L].
 *
 * Uses Levenberg-Marquardt algorithm for nonlinear optimization.
 * The dead time L is treated as a discrete search over possible values.
 *
 * Complexity: O(N * I * M) where I = iterations, M = L-search resolution.
 *
 * @param data         Step response data.
 * @param L_min        Minimum dead time to search.
 * @param L_max        Maximum dead time to search.
 * @param L_steps      Number of L grid points.
 * @param max_iters    Maximum Levenberg-Marquardt iterations.
 * @param result       [out] Identified model.
 * @return             0 on success.
 */
int fopdt_identify_ls(const step_response_data_t *data,
                      double L_min, double L_max, int L_steps,
                      int max_iters, fopdt_id_result_t *result);

/**
 * Skogestad's half-rule for approximating higher-order dynamics as FOPDT.
 *
 * If G(s) = K * exp(-θ*s) / ∏(τ_i*s + 1):
 *   T_eff = τ₁ + τ₂/2  (dominant time constant + half of second)
 *   L_eff = θ + τ₂/2 + Σ τ_i   (i ≥ 3)
 *
 * This is widely used in process control for model reduction
 * before applying PID tuning rules.
 *
 * @param K       Static gain.
 * @param tau     Array of time constants [length n_tau], sorted descending.
 * @param n_tau   Number of time constants.
 * @param theta   Original transport delay.
 * @param model   [out] Reduced FOPDT model.
 * @return        0 on success.
 */
int fopdt_skogestad_half_rule(double K, const double *tau, int n_tau,
                              double theta, fopdt_model_t *model);

/**
 * Sundaresan-Krishnaswamy (1977) identification method.
 *
 * Uses two points on the step response to identify FOPDT:
 *   y(t₁) at 35.3% and y(t₂) at 85.3% of final value.
 *
 * Formula:
 *   T = 0.68 * (t85 - t35)
 *   L = 1.3 * t35 - 0.29 * t85
 *
 * @param data   Step response data.
 * @param result [out] Identified model.
 * @return       0 on success.
 */
int fopdt_identify_sundaresan(const step_response_data_t *data,
                              fopdt_id_result_t *result);

/**
 * Prediction Error Minimization (PEM) identification of FOPDT.
 *
 * Uses one-step-ahead prediction to refine model parameters iteratively.
 * More robust to colored noise than ordinary least squares.
 *
 * @param data       Step response data.
 * @param initial    Initial parameter guess (from fast method).
 * @param max_iters  Maximum iterations.
 * @param result     [out] Refined model.
 * @return           0 on success.
 */
int fopdt_identify_pem(const step_response_data_t *data,
                       const fopdt_model_t *initial, int max_iters,
                       fopdt_id_result_t *result);

/* ──────────────────────────────────────────────
 * L5: FOPDT Simulation and Analysis
 * ────────────────────────────────────────────── */

/**
 * Simulate FOPDT step response: y(t) = K * [1 - exp(-(t-L)/T)] for t ≥ L.
 *
 * Uses analytic solution — no numerical integration needed.
 * Complexity: O(N).
 *
 * @param model   FOPDT model.
 * @param t       Time vector [length N].
 * @param y       [out] Output response [length N].
 * @param u_step  Input step magnitude.
 * @param N       Number of points.
 */
void fopdt_simulate_step(const fopdt_model_t *model, const double *t,
                         double *y, double u_step, int N);

/**
 * Simulate FOPDT response to arbitrary input using Euler integration.
 *
 * State equation: T * dy/dt + y = K * u(t-L)
 *
 * Complexity: O(N).
 *
 * @param model   FOPDT model.
 * @param t       Time vector.
 * @param u       Input signal (pre-sampled to t grid).
 * @param y0      Initial output.
 * @param y       [out] Output response.
 * @param N       Number of points.
 */
void fopdt_simulate_arbitrary(const fopdt_model_t *model, const double *t,
                              const double *u, double y0, double *y, int N);

/**
 * Compute FOPDT frequency response: G(jω) = K * exp(-jωL) / (1 + jωT).
 *
 * Magnitude: |G| = K / sqrt(1 + ω²T²)
 * Phase: ∠G = -ωL - arctan(ωT)
 *
 * @param model  FOPDT model.
 * @param omega  Angular frequency [rad/time].
 * @param mag    [out] Magnitude.
 * @param phase  [out] Phase [radians].
 */
void fopdt_freq_response(const fopdt_model_t *model, double omega,
                         double *mag, double *phase);

/**
 * Compute the normalized dead time ratio τ = L / (L + T).
 *
 * Interpretation:
 *   τ < 0.1: lag-dominant (easy to control)
 *   0.1 ≤ τ < 0.3: balanced
 *   0.3 ≤ τ < 0.7: delay-dominant (harder to control)
 *   τ ≥ 0.7: very delay-dominant (needs Smith predictor or advanced control)
 *
 * @param model FOPDT model.
 * @return      Normalized dead time ratio τ.
 */
double fopdt_dead_time_ratio(const fopdt_model_t *model);

/**
 * Compute the controllability index κ = (1/K) * max(L, T).
 *
 * This index correlates with how much control effort is needed.
 * Smaller κ means easier to control.
 *
 * @param model FOPDT model.
 * @return      Controllability index.
 */
double fopdt_controllability_index(const fopdt_model_t *model);

/**
 * Compare two FOPDT models and return a distance metric.
 *
 * Distance d = sqrt( (ΔK/K)² + (ΔT/T)² + (ΔL/L)² )
 *
 * @param a   First model.
 * @param b   Second model.
 * @return    Normalized distance between models.
 */
double fopdt_distance(const fopdt_model_t *a, const fopdt_model_t *b);

#ifdef __cplusplus
}
#endif

#endif /* FOPDT_MODEL_H */
