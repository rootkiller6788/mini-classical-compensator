/**
 * @file lead_analysis.h
 * @brief Lead Compensator Analysis Tools
 *
 * Time-domain and frequency-domain analysis.
 * L4: Sensitivity functions, Nyquist, Routh-Hurwitz
 * L5: Step/ramp response, performance metrics
 * L8: Robustness metrics, noise-PM tradeoff
 *
 * Reference: Ogata Ch.7-8, Skogestad & Postlethwaite Ch.2
 * MIT 6.302, Stanford ENGR207B, Cambridge 4F2
 */

#ifndef LEAD_ANALYSIS_H
#define LEAD_ANALYSIS_H

#include "lead_compensator.h"

#ifdef __cplusplus
extern "C" {
#endif

/* L5 - Step Response Analysis */
void lead_step_response(const lead_compensator_t *comp,
                         const lead_system_t *sys,
                         double t_final, int num_steps,
                         double *t, double *y);
void lead_step_metrics(const double *t, const double *y, int num_steps,
                        lead_performance_t *perf);
void lead_ramp_response(const lead_compensator_t *comp,
                         const lead_system_t *sys,
                         double t_final, int num_steps,
                         double *t, double *y);

/* L4 - Sensitivity Functions */
double lead_sensitivity(const lead_compensator_t *comp,
                         const lead_system_t *sys, double omega);
double lead_complementary_sensitivity(const lead_compensator_t *comp,
                                       const lead_system_t *sys, double omega);
double lead_loop_gain(const lead_compensator_t *comp,
                       const lead_system_t *sys, double omega);
double lead_peak_sensitivity(const lead_compensator_t *comp,
                              const lead_system_t *sys,
                              double w_min, double w_max, int num_points);
double lead_peak_complementary_sensitivity(const lead_compensator_t *comp,
                                            const lead_system_t *sys,
                                            double w_min, double w_max,
                                            int num_points);

/* L2 - Disturbance Rejection */
double lead_output_disturbance_gain(const lead_compensator_t *comp,
                                     const lead_system_t *sys, double omega);
double lead_input_disturbance_gain(const lead_compensator_t *comp,
                                    const lead_system_t *sys, double omega);

/* L2 - Noise and Actuator Analysis */
double lead_noise_amplification(const lead_compensator_t *comp,
                                 const lead_system_t *sys, double omega);
double lead_control_effort(const lead_compensator_t *comp,
                            const lead_system_t *sys, double omega);
double lead_max_control_effort(const lead_compensator_t *comp,
                                const lead_system_t *sys,
                                double w_min, double w_max, int num_points);

/* L2 - Bandwidth Analysis */
double lead_closed_loop_bandwidth(const lead_compensator_t *comp,
                                   const lead_system_t *sys);
double lead_open_loop_crossover(const lead_compensator_t *comp,
                                 const lead_system_t *sys);

/* L8 - Advanced Robustness Metrics */
double lead_modulus_margin(const lead_compensator_t *comp,
                            const lead_system_t *sys,
                            double w_min, double w_max, int num_points);
double lead_delay_margin(const lead_compensator_t *comp,
                          const lead_system_t *sys);
double lead_noise_pm_tradeoff(const lead_compensator_t *comp,
                               const lead_system_t *sys,
                               const lead_compensator_t *uncompensated);
bool lead_excessive_hf_gain(const lead_compensator_t *comp, double threshold_db);

/* L4 - Closed-loop Pole Analysis */
void lead_closed_loop_poles(const lead_compensator_t *comp,
                             const lead_system_t *sys,
                             lead_complex_t *poles, int *num_poles);
int lead_count_unstable_cl_poles(const lead_compensator_t *comp,
                                  const lead_system_t *sys);

#ifdef __cplusplus
}
#endif

#endif /* LEAD_ANALYSIS_H */
