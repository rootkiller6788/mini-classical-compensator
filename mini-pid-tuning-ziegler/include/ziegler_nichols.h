#ifndef ZIEGLER_NICHOLS_H
#define ZIEGLER_NICHOLS_H

/**
 * ziegler_nichols.h — Ziegler-Nichols PID Tuning Methods
 *
 * Knowledge coverage:
 *   L1: Z-N step response method definition (1942)
 *   L1: Z-N frequency response (ultimate sensitivity) method (1942)
 *   L2: Process reaction curve concept — tangent line at inflection
 *   L2: Ultimate gain/period as stability margin boundary
 *   L4: Ziegler-Nichols tuning rules (original empirical correlations)
 *   L5: Algorithmic implementation of both Z-N methods
 *   L6: Canonical FOPDT/SOPDT process identification
 *
 * Reference:
 *   Ziegler, J.G. & Nichols, N.B. (1942)
 *   "Optimum Settings for Automatic Controllers", Trans. ASME, 64, 759-768.
 *
 *   Åström, K.J. & Hägglund, T. (2004)
 *   "Revisiting the Ziegler-Nichols Step Response Method for PID Control",
 *   Journal of Process Control, 14, 635-650.
 */

#include "pid_tuning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────
 * L1: Ziegler-Nichols Tuning Rule Set Definitions
 * ────────────────────────────────────────────── */

/**
 * Z-N tuning rule table entries for step response method.
 *
 * The original Z-N rules relate PID parameters to the process reaction
 * curve parameters: K (static gain), L (apparent dead time), T (time
 * constant), from the tangent line method on the step response.
 *
 * For FOPDT model:  G(s) = K * exp(-L*s) / (T*s + 1)
 *
 * Parameter 'a' = K * L / T  (normalized slope at inflection)
 */
typedef enum {
    ZN_CONTROLLER_P  = 0,  /* Proportional only */
    ZN_CONTROLLER_PI = 1,  /* Proportional + Integral */
    ZN_CONTROLLER_PD = 2,  /* Proportional + Derivative */
    ZN_CONTROLLER_PID= 3   /* Three-term PID */
} zn_controller_type_t;

/**
 * L4: Z-N step response tuning rule table.
 *
 * These coefficients were determined empirically to give
 * quarter-amplitude damping ratio (ζ ≈ 0.21) for a wide
 * range of processes.
 *
 *   Controller | Kp         | Ti      | Td
 *   -----------|------------|---------|------
 *   P          | 1/a        | —       | —
 *   PI         | 0.9/a      | 3.0*L   | —
 *   PID        | 1.2/a      | 2.0*L   | 0.5*L
 *
 * where a = K * L / T.
 */
typedef struct {
    double Kp_factor;  /* Multiplier for 1/a or T/(K*L) */
    double Ti_factor;  /* Multiplier for L (integral time) */
    double Td_factor;  /* Multiplier for L (derivative time) */
} zn_rule_entry_t;

/**
 * L4: Z-N frequency response (ultimate sensitivity) rule table.
 *
 * Based on Ku (ultimate gain) and Pu (ultimate period).
 *
 *   Controller | Kp         | Ti        | Td
 *   -----------|------------|-----------|------
 *   P          | 0.50*Ku    | —         | —
 *   PI         | 0.45*Ku    | Pu/1.2    | —
 *   PID        | 0.60*Ku    | Pu/2.0    | Pu/8.0
 *
 * These give gain margin ≈ 2 and phase margin ≈ 30-45° for typical
 * industrial processes.
 */
typedef struct {
    double Kp_factor;  /* Multiplier for Ku */
    double Ti_factor;  /* Multiplier for Pu */
    double Td_factor;  /* Multiplier for Pu */
} zn_freq_rule_entry_t;

/**
 * Modified Ziegler-Nichols tuning options (for overshoot control).
 *
 * Åström & Hägglund (2004) documented that the original Z-N rules
 * produce ζ ≈ 0.21 (≈45% overshoot), which is often too aggressive.
 * Modified rules provide tuning for different damping targets.
 */
typedef enum {
    ZN_MOD_STANDARD = 0,    /* Original 1942 (ζ ≈ 0.21) */
    ZN_MOD_CONSERVATIVE,    /* Reduced gain (ζ ≈ 0.5) */
    ZN_MOD_AGGRESSIVE,      /* Increased gain (ζ ≈ 0.1) */
    ZN_MOD_NO_OVERSHOOT,    /* Critically damped target */
    ZN_MOD_SOME_OVERSHOOT   /* 10-20% overshoot target */
} zn_modified_t;

/* ──────────────────────────────────────────────
 * L5: Z-N Step Response Method Functions
 * ────────────────────────────────────────────── */

/**
 * Identify FOPDT model parameters from step response data using the
 * classical tangent (slope) method of Ziegler-Nichols.
 *
 * Algorithm (L5):
 *   1. Find the inflection point (maximum slope) of the step response.
 *   2. Draw tangent line at inflection point.
 *   3. L = intercept of tangent line with initial value (time axis).
 *   4. T = (steady state - intercept) / slope = time for tangent to
 *          reach final value.
 *   5. K = Δoutput / Δinput (steady-state gain).
 *
 * Complexity: O(N) data scan. N = number of data points.
 *
 * @param data  Step response experiment.
 * @param model [out] Identified FOPDT model parameters.
 * @return      0 on success, -1 if data insufficient or no inflection found.
 */
int zn_identify_fopdt(const step_response_data_t *data, fopdt_model_t *model);

/**
 * Compute Z-N step response PID tuning parameters.
 *
 * Given an identified FOPDT model, this computes the PID (or P/PI/PD)
 * parameters using the original Ziegler-Nichols step response rules.
 *
 * @param model  FOPDT process model.
 * @param type   Controller type (P/PI/PD/PID).
 * @param params [out] Computed PID parameters (IDEAL form).
 * @return       0 on success, -1 if model is invalid (e.g., T ≤ 0).
 */
int zn_step_tune(const fopdt_model_t *model, zn_controller_type_t type,
                 pid_params_t *params);

/**
 * Compute modified Z-N step response tuning with reduced overshoot.
 *
 * Reference: Åström & Hägglund (2004) revision of Z-N rules.
 *
 * @param model    FOPDT process model.
 * @param type     Controller type.
 * @param variant  Overshoot target variant.
 * @param Ms       Desired maximum sensitivity (1.2 ≤ Ms ≤ 2.0).
 * @param params   [out] Computed PID parameters.
 * @return         0 on success.
 */
int zn_step_tune_modified(const fopdt_model_t *model, zn_controller_type_t type,
                          zn_modified_t variant, double Ms,
                          pid_params_t *params);

/**
 * Compute Z-N step response tuning for SOPDT processes.
 *
 * Transforms SOPDT → approximate FOPDT before applying rules.
 * The approximation preserves the dominant time constant and
 * absorbs the secondary time constant into dead time.
 *
 * @param model  SOPDT process model.
 * @param type   Controller type.
 * @param params [out] Computed PID parameters.
 * @return       0 on success.
 */
int zn_step_tune_sopdt(const sopdt_model_t *model, zn_controller_type_t type,
                       pid_params_t *params);

/* ──────────────────────────────────────────────
 * L5: Z-N Frequency Response Method Functions
 * ────────────────────────────────────────────── */

/**
 * Simulate or analyze a relay feedback experiment to find ultimate gain
 * and ultimate period (Ku, Pu).
 *
 * In an actual plant, this would use an Åström-Hägglund relay experiment.
 * Here we provide the algorithmic computation: given a process transfer
 * function, analytically determine Ku and Pu.
 *
 * For FOPDT: solve |G(jω)| = 1/Ku, ∠G(jω) = -π at ω = 2π/Pu.
 *
 * Complexity: O(iterations) for numerical frequency search.
 *
 * @param fopdt   FOPDT model of the plant.
 * @param result  [out] Ku and Pu.
 * @return        0 on success, -1 if no crossing found.
 */
int zn_find_ultimate_gain(const fopdt_model_t *fopdt,
                          ultimate_gain_result_t *result);

/**
 * Compute Z-N frequency response PID tuning parameters.
 *
 * Given Ku (ultimate gain) and Pu (ultimate period), apply the
 * Ziegler-Nichols frequency response tuning rules.
 *
 * @param Ku      Ultimate gain.
 * @param Pu      Ultimate period.
 * @param type    Controller type (P/PI/PD/PID).
 * @param params  [out] Computed PID parameters (IDEAL form).
 * @return        0 on success, -1 if Ku ≤ 0 or Pu ≤ 0.
 */
int zn_freq_tune(double Ku, double Pu, zn_controller_type_t type,
                 pid_params_t *params);

/**
 * Modified Z-N frequency response tuning with settable parameters.
 *
 * Allows adjusting the Z-N multipliers. The standard factors are:
 *   α:Kp_factor = 0.6, β:Ti_factor = 0.5, γ:Td_factor = 0.125 (PID)
 *
 * @param Ku        Ultimate gain.
 * @param Pu        Ultimate period.
 * @param type      Controller type.
 * @param alpha     Kp factor (e.g., 0.6 for standard, 0.3 for conservative).
 * @param beta      Ti factor (e.g., 0.5 for standard).
 * @param gamma     Td factor (e.g., 0.125 for standard).
 * @param params    [out] PID parameters.
 * @return          0 on success.
 */
int zn_freq_tune_custom(double Ku, double Pu, zn_controller_type_t type,
                        double alpha, double beta, double gamma,
                        pid_params_t *params);

/**
 * Compute the gain margin and phase margin of a PID-controlled FOPDT loop
 * analytically for given PID tuning.
 *
 * This verifies the robustness of a given Z-N tuning.
 *   GM = 1 / |G(jω_pc) * C(jω_pc)|  where ∠GC = -π
 *   PM = π + ∠G(jω_gc) * C(jω_gc)   where |GC| = 1
 *
 * @param fopdt  Process model.
 * @param params PID parameters.
 * @param Gm     [out] Gain margin [linear, not dB].
 * @param Pm     [out] Phase margin [radians].
 * @return       0 on success.
 */
int zn_verify_margins(const fopdt_model_t *fopdt, const pid_params_t *params,
                      double *Gm, double *Pm);

/**
 * Sørensen's method (1958): alternative ultimate gain identification
 * for processes where sustained oscillation cannot be reached safely.
 *
 * Uses a relay with hysteresis width ε to generate a stable oscillation
 * with period P_osc. Then:
 *   Ku = 4*d / (π * a)
 *   Pu = P_osc
 * where d = relay amplitude, a = oscillation amplitude at process output.
 *
 * @param relay_amplitude   Relay output amplitude d.
 * @param hysteresis_width  Relay hysteresis ε.
 * @param oscillation_ampl  Measured oscillation amplitude a.
 * @param oscillation_period Measured oscillation period P_osc.
 * @param result [out] Ku and Pu.
 * @return 0 on success.
 */
int zn_sorensen_ultimate_gain(double relay_amplitude, double hysteresis_width,
                              double oscillation_ampl, double oscillation_period,
                              ultimate_gain_result_t *result);

/**
 * Chien-Hrones-Reswick (1952) modification to Z-N step response method.
 *
 * Provides two sets of rules: one for setpoint tracking (0% overshoot)
 * and one for disturbance rejection (20% overshoot).
 *
 * Reference: Chien, Hrones, Reswick (1952), Trans. ASME.
 *
 * @param model  FOPDT model.
 * @param type   Controller type.
 * @param mode   0 = setpoint (0% OS), 1 = disturbance (20% OS).
 * @param params [out] PID parameters.
 * @return       0 on success.
 */
int zn_chien_hrones_reswick(const fopdt_model_t *model,
                            zn_controller_type_t type, int mode,
                            pid_params_t *params);

#ifdef __cplusplus
}
#endif

#endif /* ZIEGLER_NICHOLS_H */
