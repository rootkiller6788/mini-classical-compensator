#ifndef GAIN_MARGIN_TUNING_H
#define GAIN_MARGIN_TUNING_H

/**
 * gain_margin_tuning.h — Gain & Phase Margin Based PID Tuning
 *
 * Knowledge coverage:
 *   L1: Gain margin (Gm) and phase margin (Pm) definitions
 *   L2: Stability margin concept — distance to Nyquist critical point -1
 *   L4: Nyquist stability criterion applied to PID design
 *   L5: Frequency-domain PID tuning by margin specification
 *   L6: Canonical design: Gm = 2-5 (6-14 dB), Pm = 30-60°
 *   L8: Robust PID design with guaranteed margins
 *
 * Reference:
 *   Åström, K.J. & Hägglund, T. (1995) Ch.5 "Design Methods".
 *   Ho, W.K., Hang, C.C. & Cao, L.S. (1995)
 *   "Tuning of PID Controllers Based on Gain and Phase Margin Specifications",
 *   Automatica, 31(3), 497-502.
 *   Ho, W.K., Lim, K.W. & Xu, W. (1998)
 *   "Optimal Gain and Phase Margin Tuning for PID Controllers", Automatica.
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: Margin Specification Definition
 * ────────────────────────────────────────────── */

/**
 * Gain and phase margin specification structure.
 *
 * Am = gain margin (linear, not dB): Am > 1 means stable.
 *   Am_dB = 20 * log10(Am)
 * Pm = phase margin [radians]. Typically 30° to 60°.
 *
 * Classical recommendations:
 *   Servo (setpoint tracking):  Pm ≥ 45°, Am ≥ 2 (6 dB)
 *   Regulatory (disturbance):   Pm ≥ 30°, Am ≥ 2.5 (8 dB)
 *   Robust (model uncertainty): Pm ≥ 60°, Am ≥ 3 (9.5 dB)
 */
typedef struct {
    double Am;     /* Gain margin (linear: Am > 1 for stability) */
    double Pm;     /* Phase margin [radians] */
} margin_spec_t;

/**
 * Frequency design points for margin-based tuning.
 *
 * ω_gc = gain crossover frequency (where |L(jω)| = 1)
 * ω_pc = phase crossover frequency (where ∠L(jω) = -π)
 *
 * At design:
 *   |C(jω_pc)*G(jω_pc)| = 1/Am
 *   ∠C(jω_gc)*G(jω_gc) = -π + Pm
 */
typedef struct {
    double omega_gc; /* Gain crossover frequency [rad/time] */
    double omega_pc; /* Phase crossover frequency [rad/time] */
} margin_design_freqs_t;

/* ──────────────────────────────────────────────
 * L4: Margin Computation (Nyquist Stability)
 * ────────────────────────────────────────────── */

/**
 * Compute gain and phase margins for a given FOPDT process + PID tuning.
 *
 * Finds ω_pc such that ∠G(jω)*C(jω) = -π,
 * then Am = 1 / |G(jω_pc)*C(jω_pc)|.
 * Finds ω_gc such that |G(jω)*C(jω)| = 1,
 * then Pm = π + ∠G(jω_gc)*C(jω_gc).
 *
 * Uses numerical search (bisection/Newton) for frequency crossing points.
 * Complexity: O(log(ω_max/ω_min)) per search.
 *
 * @param fopdt   Process model.
 * @param params  PID parameters.
 * @param margins [out] Computed gain and phase margins.
 * @return        0 on success, -1 if margins cannot be determined.
 */
int gm_pm_compute_margins(const fopdt_model_t *fopdt,
                          const pid_params_t *params,
                          margin_spec_t *margins);

/* ──────────────────────────────────────────────
 * L5: Margin-Specified PID Tuning
 * ────────────────────────────────────────────── */

/**
 * PI tuning by gain margin and phase margin specification.
 *
 * Given FOPDT: G(s) = K * exp(-L*s) / (T*s + 1)
 * Given desired Am, Pm, this computes Kp and Ti analytically.
 *
 * Approach (Ho-Hang-Cao, 1995):
 *   From the phase condition:
 *     ω_gc * L + arctan(ω_gc * T) = π/2 - Pm + arctan(1/(ω_gc*Ti))
 *   From the gain condition at ω_pc:
 *     Am = (1/K) * sqrt((1 + ω_pc²*T²) * ω_pc²*Ti²) / (Kp*sqrt(1 + ω_pc²*Ti²))
 *
 * If the equations are solvable, this yields a unique (Kp, Ti) for PI.
 *
 * @param fopdt   Process model.
 * @param Am      Desired gain margin (linear, > 1).
 * @param Pm      Desired phase margin [radians].
 * @param params  [out] PI parameters (Kp, Ki, 0 for D).
 * @return        0 on success, -1 if no feasible solution.
 */
int gm_pm_tune_pi(const fopdt_model_t *fopdt, double Am, double Pm,
                  pid_params_t *params);

/**
 * PID tuning by gain margin and phase margin specification.
 *
 * Extends the PI design above by also specifying the derivative gain.
 * The derivative adds phase lead at ω_gc, allowing higher gain.
 *
 * For PID, the design equations are:
 *   ∠C(jω_gc) + ∠G(jω_gc) = -π + Pm   (phase condition)
 *   |C(jω_pc) * G(jω_pc)| = 1/Am       (gain condition)
 *   ∂(∠L(jω))/∂ω |_{ω=ω_gc} = 0       (flat phase — iso-damping)
 *
 * The iso-damping condition ensures robustness to gain variations.
 *
 * Reference: Chen, Yang, & Huang (2003).
 *
 * @param fopdt   Process model.
 * @param Am      Desired gain margin (linear).
 * @param Pm      Desired phase margin [radians].
 * @param use_flat_phase 1 = enforce flat-phase (iso-damping) condition.
 * @param params  [out] PID parameters.
 * @return        0 on success.
 */
int gm_pm_tune_pid(const fopdt_model_t *fopdt, double Am, double Pm,
                   int use_flat_phase, pid_params_t *params);

/**
 * PI tuning for integrating processes by margin specification.
 *
 * For G(s) = Kv * exp(-L*s) / s:
 *
 *   Kp = (ω_gc / Kv) * sqrt( (1 + (ω_gc*Ti)⁻²)⁻¹ )
 *   Ti = (tan(Pm - π + arctan(ω_gc*L) + π/2)) / ω_gc
 *
 * @param ipdt    Integrating process model.
 * @param Am      Desired gain margin.
 * @param Pm      Desired phase margin [radians].
 * @param params  [out] PI parameters.
 * @return        0 on success.
 */
int gm_pm_tune_ipdt(const ipdt_model_t *ipdt, double Am, double Pm,
                    pid_params_t *params);

/**
 * Time-domain margin verification: given PID tuning, simulate step
 * response and verify that the actual gain/phase margins match.
 *
 * This bridges frequency-domain specifications to time-domain performance.
 *
 * @param fopdt       Process model.
 * @param params      PID parameters.
 * @param sim_time    Simulation duration.
 * @param n_steps     Number of simulation steps.
 * @param margins_in  [in] Target margins.
 * @param margins_out [out] Actual margins from frequency analysis.
 * @param perf        [out] Time-domain performance from step response.
 * @return            0 on success.
 */
int gm_pm_verify_tuning(const fopdt_model_t *fopdt,
                        const pid_params_t *params,
                        double sim_time, int n_steps,
                        const margin_spec_t *margins_in,
                        margin_spec_t *margins_out,
                        pid_perf_t *perf);

/**
 * Compute the maximum sensitivity (Ms) from PID + FOPDT.
 *
 * Ms = max_{ω} |S(jω)| = max_{ω} |1 / (1 + G(jω)*C(jω))|
 *
 * Ms is the inverse of the shortest distance from the Nyquist curve
 * to the critical point -1. Typical ranges:
 *   Ms = 1.2-1.4  Conservative
 *   Ms = 1.4-1.6  Normal
 *   Ms = 1.6-2.0  Aggressive
 *   Ms > 2.0      Potentially unstable
 *
 * @param fopdt  Process model.
 * @param params PID parameters.
 * @param Ms     [out] Maximum sensitivity.
 * @param omega_Ms [out] Frequency at which Ms occurs.
 * @return       0 on success.
 */
int gm_pm_max_sensitivity(const fopdt_model_t *fopdt,
                          const pid_params_t *params,
                          double *Ms, double *omega_Ms);

/**
 * Sensitivity-constrained tuning: design PID to achieve a specific
 * maximum sensitivity Ms.
 *
 * This is a frequency-loop-shaping approach: iterate on Kp, Ti until
 * the Nyquist curve touches a circle of radius 1/Ms around -1.
 *
 * @param fopdt    Process model.
 * @param Ms_target Desired maximum sensitivity (1.2 ≤ Ms ≤ 2.0).
 * @param params   [out] PID parameters.
 * @return         0 on success.
 */
int gm_pm_tune_by_Ms(const fopdt_model_t *fopdt, double Ms_target,
                     pid_params_t *params);

#ifdef __cplusplus
}
#endif

#endif /* GAIN_MARGIN_TUNING_H */
