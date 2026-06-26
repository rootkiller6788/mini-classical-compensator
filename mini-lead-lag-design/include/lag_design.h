#ifndef LAG_DESIGN_H
#define LAG_DESIGN_H

/** @file lag_design.h
 * Lag Compensator: C(s) = Kc*(1+T*s)/(1+beta*T*s), beta > 1
 * Reduces steady-state error. Attenuation at HF: 1/beta.
 * Ref: Ogata Ch7, Nise Ch9, Dorf & Bishop Ch10
 */
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double Kc; double T; double beta; } LagCompensator;

typedef struct {
    double ess_desired;        /* Desired steady-state error ratio */
    double ess_current;        /* Current steady-state error ratio */
    double gain_crossover;     /* Current gain crossover (rad/s) */
    double dc_gain_required;   /* Required DC gain increase factor */
} LagDesignSpec;

typedef struct {
    double omega; double magnitude; double magnitude_db;
    double phase_deg; double real; double imag;
} LagFreqPoint;

/* L4 - Fundamental Theorems */
double lag_attenuation_db(double beta);
double lag_max_phase_lag_deg(double beta);
double lag_beta_from_attenuation(double atten_db);
double lag_beta_from_dc_gain(double gain_ratio);
int lag_design_from_spec(const LagDesignSpec *spec, LagCompensator *comp);
int lag_design_direct(double ess_current, double ess_desired,
                      double wc, double ratio, LagCompensator *comp);

/* L5 - Analysis */
LagFreqPoint lag_evaluate(const LagCompensator *comp, double omega);
double lag_magnitude_db(const LagCompensator *comp, double omega);
double lag_phase_deg(const LagCompensator *comp, double omega);
double lag_dc_gain(const LagCompensator *comp);
double lag_hf_gain(const LagCompensator *comp);
void lag_pole_zero(const LagCompensator *comp, double *zero, double *pole);
int lag_is_valid(const LagCompensator *comp);

/* L5 - Time Domain */
void lag_step_response(const LagCompensator *comp, double t_final,
                       int num_points, double *t_out, double *y_out);
void lag_impulse_response(const LagCompensator *comp, double t_final,
                          int num_points, double *t_out, double *y_out);

/* L5 - Discretization */
void lag_discretize_tustin(const LagCompensator *comp, double Ts,
                           double b0_b1[2], double *a1);
void lag_discretize_zoh(const LagCompensator *comp, double Ts,
                        double b0_b1[2], double *a1);
double lag_apply_digital(const double b0_b1[2], double a1,
                         double input, double *u_prev, double *y_prev);

/* L5 - Bode */
void lag_bode_data(const LagCompensator *comp,
                   double omega_min, double omega_max, int num_points,
                   double *omega_out, double *mag_db_out, double *phase_out);

#ifdef __cplusplus
}
#endif
#endif
