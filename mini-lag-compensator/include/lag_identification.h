/**
 * @file lag_identification.h
 * @brief System identification for lag compensator design
 *
 * Methods to identify plant characteristics needed for compensator design:
 * DC gain, dominant time constant, gain crossover frequency, error constants.
 *
 * L5: System identification algorithms (step response, frequency response)
 * L6: DC motor parameter estimation from step response data
 * L7: Industrial process identification (temperature, position servo)
 *
 * Course alignment:
 *   MIT 6.302 — system identification for controller design
 *   Stanford ENGR105 — plant modeling from response data
 *   Berkeley ME132 — experimental system identification
 *   Caltech CDS 110 — data-driven model extraction
 *   ETH 151-0591 — Identifikation aus Sprungantwort
 *   Cambridge 3F2 — system identification methods
 *   Tsinghua 自动控制原理 — 系统辨识
 *
 * Textbook: Ljung, "System Identification" (1999)
 *           Ogata, Ch. 10 — system identification for control
 */

#ifndef LAG_IDENTIFICATION_H
#define LAG_IDENTIFICATION_H

#include "lag_types.h"

/* ==========================================================================
 * L5: Step response identification
 * ========================================================================== */

/**
 * Identify first-order plant parameters from step response data.
 *
 * For a first-order system: G(s) = K / (tau*s + 1)
 *
 * From step response y(t) = K * (1 - exp(-t/tau)):
 *   - K = final_value / step_magnitude (DC gain)
 *   - tau = time when y(t) reaches 63.2% of final value
 *   - Or: tau = time_constant from slope at origin
 *
 * Algorithm:
 *   1. Find steady-state final value (average of last 10% of data)
 *   2. Compute K = y_ss / step_amplitude
 *   3. Find time t_63 where y(t_63) = 0.632 * y_ss
 *   4. tau = t_63 (for first-order system)
 *
 * Complexity: O(n_points)
 *
 * @param response  step response trajectory data
 * @param step_amp  amplitude of the input step
 * @param[out] K    identified DC gain
 * @param[out] tau  identified time constant (seconds)
 * @return 0 on success, negative on insufficient data
 */
int lag_identify_first_order(const LagStepResponse *response,
                              double step_amp,
                              double *K, double *tau);

/**
 * Identify second-order plant parameters from step response.
 *
 * For G(s) = K * w_n^2 / (s^2 + 2*zeta*w_n*s + w_n^2)
 *
 * From step response overshoot and period:
 *   zeta = sqrt(ln^2(M_p) / (pi^2 + ln^2(M_p)))
 *   w_n = pi / (T_p * sqrt(1 - zeta^2))
 *   K = y_ss / step_amp
 *
 * where M_p = overshoot ratio, T_p = peak time.
 *
 * Algorithm (logarithmic decrement method):
 *   1. Compute K from steady-state value
 *   2. Measure overshoot M_p = (y_peak - y_ss)/y_ss
 *   3. Compute damping zeta from M_p formula
 *   4. Measure peak time T_p
 *   5. Compute natural frequency w_n
 *
 * Complexity: O(n_points)
 *
 * @param response  step response trajectory
 * @param step_amp  input step amplitude
 * @param[out] K    identified DC gain
 * @param[out] zeta identified damping ratio
 * @param[out] wn   identified natural frequency (rad/s)
 * @return 0 on success
 */
int lag_identify_second_order(const LagStepResponse *response,
                               double step_amp,
                               double *K, double *zeta, double *wn);

/**
 * Identify DC gain (steady-state gain) from step response.
 *
 * Simply computes: K = (final_value - initial_value) / step_amplitude
 *
 * Complexity: O(1) (uses pre-computed final_value)
 *
 * @param response  step response data
 * @param step_amp  input step amplitude
 * @return identified DC gain, or 0 on error
 */
double lag_identify_dc_gain(const LagStepResponse *response,
                             double step_amp);

/**
 * Identify the dominant time constant using the 63.2% method.
 *
 * For a predominantly first-order response, the time at which
 * the output reaches 63.2% of its final value equals the time constant.
 *
 * Complexity: O(n_points) for linear search
 *
 * @param response  step response data
 * @return dominant time constant in seconds
 */
double lag_identify_time_constant(const LagStepResponse *response);

/**
 * Estimate the time delay (dead time) from step response.
 *
 * Dead time is the interval between step application and
 * the first significant response (e.g., output exceeds 5% of final).
 *
 * This is critical for processes with transport delay
 * (chemical plants, long pipelines, thermal systems).
 *
 * Complexity: O(n_points)
 */
double lag_identify_dead_time(const LagStepResponse *response,
                               double threshold_fraction);

/* ==========================================================================
 * L5: Frequency response identification
 * ========================================================================== */

/**
 * Identify DC gain and dominant pole from low-frequency Bode data.
 *
 * From the low-frequency magnitude asymptote:
 *   |G(jw)| ~= K  for w << 1/tau (flat region)
 *
 * The corner frequency is where magnitude drops by 3 dB.
 *
 * Complexity: O(n_points)
 */
int lag_identify_from_bode(const LagBodeData *bode,
                            double *K, double *corner_freq);

/**
 * Compute the error constants Kp, Kv, Ka from Bode data.
 *
 * Uses low-frequency asymptote extrapolation:
 *   Kp = |G(0)| asymptote
 *   Kv = w * |G(jw)| for w -> 0 (requires -20 dB/dec initial slope)
 *   Ka = w^2 * |G(jw)| for w -> 0 (requires -40 dB/dec initial slope)
 *
 * Complexity: O(1) using lowest-frequency data
 */
void lag_identify_error_constants(const LagBodeData *bode,
                                   double *Kp, double *Kv, double *Ka);

/* ==========================================================================
 * L6: DC motor identification
 * ========================================================================== */

/**
 * Identify DC motor parameters from step response data.
 *
 * For a DC motor, the speed transfer function is:
 *   Omega(s)/V(s) = Kt / ((Js+B)*(Ls+R) + Kt*Kb)
 *
 * Which approximates to first-order when L is small:
 *   Omega(s)/V(s) = Km / (tau_m*s + 1)
 *   where Km = Kt/(B*R + Kt*Kb), tau_m = J*R/(B*R + Kt*Kb)
 *
 * This function identifies Km and tau_m from speed step response.
 *
 * Complexity: O(n_points)
 *
 * @param response       speed step response data
 * @param voltage_step   applied voltage step (V)
 * @param[out] params    identified motor parameters
 * @return 0 on success
 */
int lag_identify_dc_motor(const LagStepResponse *response,
                           double voltage_step,
                           LagDCMotorParams *params);

/**
 * Estimate DC motor electrical time constant from current response.
 *
 * tau_e = L/R can be identified from the current rise time
 * when the rotor is locked (back-EMF = 0).
 *
 * Complexity: O(n_points)
 */
double lag_identify_electrical_tc(const LagStepResponse *current_response,
                                   const LagDCMotorParams *params);

/**
 * Estimate DC motor mechanical time constant from speed response.
 *
 * tau_m is the time for speed to reach 63.2% of final speed
 * under constant voltage input.
 *
 * Complexity: O(n_points)
 */
double lag_identify_mechanical_tc(const LagStepResponse *speed_response);

/* ==========================================================================
 * L7: Industrial process identification
 * ========================================================================== */

/**
 * Identify a FOPDT (First-Order Plus Dead Time) model from step response.
 *
 * G(s) = K * exp(-theta*s) / (tau*s + 1)
 *
 * This is the most common industrial process model, used extensively
 * in chemical plants, refineries, and thermal processes.
 *
 * Uses the two-point method (Smith, 1957):
 *   Find times t1 and t2 where y = 28.3% and 63.2% of final value
 *   tau = 1.5 * (t2 - t1)
 *   theta = t2 - tau
 *   K = y_ss / step_amp
 *
 * Complexity: O(n_points)
 */
int lag_identify_fopdt(const LagStepResponse *response,
                        double step_amp,
                        double *K, double *tau, double *theta);

/**
 * Identify process model using area method (robust to noise).
 *
 * For noisy step response data, the area method provides
 * more robust parameter estimates than the tangent method.
 *
 * Complexity: O(n_points)
 */
int lag_identify_area_method(const LagStepResponse *response,
                              double step_amp,
                              double *K, double *tau, double *theta);

#endif /* LAG_IDENTIFICATION_H */