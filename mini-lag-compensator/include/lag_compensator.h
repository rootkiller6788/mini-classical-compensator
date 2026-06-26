/**
 * @file lag_compensator.h
 * @brief Core lag compensator: definition, creation, evaluation, analysis
 *
 * The lag compensator (phase-lag controller) improves steady-state accuracy
 * while attempting to preserve transient response and stability margins.
 *
 * Transfer function (beta form):
 *   G_c(s) = K_c * (T*s + 1) / (beta*T*s + 1)    with beta > 1
 *
 * Knowledge map:
 *   L1 Definitions:   LagCompensator struct, parameter forms
 *   L2 Core Concepts:  DC gain, corner frequencies, lag principle
 *   L3 Math Struct:    rational TF in s-domain, complex evaluation
 *   L4 Theorems:       steady-state error improvement, phase lag bound
 *   L5 Algorithms:     creation, evaluation, validation
 *
 * Course: MIT 6.302 / Stanford ENGR105 / Berkeley ME132/232 / Caltech CDS 110
 *         ETH 151-0591 / Cambridge 3F2 / Tsinghua
 *
 * Textbook: Ogata Ch.7; Franklin/Powell/Emami-Naeini Ch.6
 */

#ifndef LAG_COMPENSATOR_H
#define LAG_COMPENSATOR_H

#include "lag_types.h"

typedef struct {
    double Kc;
    double T;
    double beta;
    double zero;
    double pole;
    double dc_gain;
    double hf_gain;
    double max_phase_lag_rad;
    double max_lag_freq;
} LagCompensator;

LagCompensator lag_create(double Kc, double T, double beta);
LagCompensator lag_create_from_pole_zero(double Kc, double zero, double pole);
LagCompensator lag_create_from_corners(double Kc, double omega_z, double omega_p);
LagCompensator lag_create_identity(void);
double lag_get_zero(const LagCompensator *lag);
double lag_get_pole(const LagCompensator *lag);
double lag_get_dc_gain(const LagCompensator *lag);
double lag_get_hf_gain(const LagCompensator *lag);
double lag_get_beta(const LagCompensator *lag);
double lag_get_time_constant(const LagCompensator *lag);
double lag_get_zero_pole_ratio(const LagCompensator *lag);
double lag_ess_improvement(const LagCompensator *lag, LagESSType ess_type);
LagComplex lag_eval_s(const LagCompensator *lag, LagComplex s);
double lag_eval_magnitude(const LagCompensator *lag, double omega);
double lag_eval_phase(const LagCompensator *lag, double omega);
LagFreqPoint lag_eval_frequency(const LagCompensator *lag, double omega);
void lag_get_corner_frequencies(const LagCompensator *lag, double *omega_low, double *omega_high);
double lag_get_max_lag_frequency(const LagCompensator *lag);
double lag_get_max_phase_lag(const LagCompensator *lag);
double lag_low_freq_asymptote_db(const LagCompensator *lag);
double lag_high_freq_asymptote_db(const LagCompensator *lag);
LagTransferFunction lag_to_transfer_function(const LagCompensator *lag);
LagCompensator lag_from_first_order_tf(double b1, double b0, double a1, double a0);
int lag_is_stable(const LagCompensator *lag);
int lag_is_minimum_phase(const LagCompensator *lag);
int lag_validate(const LagCompensator *lag);
const char* lag_validate_error_string(int error_code);
void lag_to_string(const LagCompensator *lag, char *buf, int bufsz);
int lag_equals(const LagCompensator *a, const LagCompensator *b, double tol);

#endif /* LAG_COMPENSATOR_H */
