#ifndef LEAD_DESIGN_H
#define LEAD_DESIGN_H

/** @file lead_design.h
 * Lead Compensator: C(s) = Kc * (1+T*s)/(1+alpha*T*s), 0 < alpha < 1
 * Max phase: phi_max = arcsin((1-alpha)/(1+alpha))
 * Freq at max: omega_m = 1/(T*sqrt(alpha))
 * Ref: Ogata Ch7, Nise Ch9, Dorf & Bishop Ch10, Franklin Ch6
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double Kc; double T; double alpha; } LeadCompensator;

typedef struct {
    double phase_margin_desired;
    double phase_margin_current;
    double gain_crossover;
    double gain_crossover_desired;
    double phase_margin_extra;
    double max_alpha;
} LeadDesignSpec;

typedef struct {
    double omega;
    double magnitude;
    double magnitude_db;
    double phase_deg;
    double real;
    double imag;
} LeadFreqPoint;

/* L4 - Fundamental Theorems */
double lead_max_phase_deg(double alpha);
double lead_omega_max(double T, double alpha);
double lead_alpha_from_phase(double phase_lead_deg);
int lead_design_from_spec(const LeadDesignSpec *spec, LeadCompensator *comp);
int lead_design_direct(double pm_current, double pm_desired, double wc, double margin_deg, LeadCompensator *comp);

/* L5 - Frequency Analysis */
LeadFreqPoint lead_evaluate(const LeadCompensator *comp, double omega);
double lead_magnitude_db(const LeadCompensator *comp, double omega);
double lead_phase_deg(const LeadCompensator *comp, double omega);
double lead_dc_gain(const LeadCompensator *comp);
double lead_hf_gain(const LeadCompensator *comp);
void lead_pole_zero(const LeadCompensator *comp, double *zero, double *pole);
int lead_is_valid(const LeadCompensator *comp);

/* L5 - Time Domain */
void lead_step_response(const LeadCompensator *comp, double t_final, int num_points, double *t_out, double *y_out);
void lead_impulse_response(const LeadCompensator *comp, double t_final, int num_points, double *t_out, double *y_out);

/* L5 - Discretization */
void lead_discretize_tustin(const LeadCompensator *comp, double Ts, double b0_b1[2], double *a1);
void lead_discretize_zoh(const LeadCompensator *comp, double Ts, double b0_b1[2], double *a1);
double lead_apply_digital(const double b0_b1[2], double a1, double input, double *u_prev, double *y_prev);

/* L5 - Bode */
void lead_bode_data(const LeadCompensator *comp, double omega_min, double omega_max, int num_points, double *omega_out, double *mag_db_out, double *phase_out);

#ifdef __cplusplus
}
#endif
#endif