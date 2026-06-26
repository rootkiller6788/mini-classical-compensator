#ifndef LEAD_LAG_DESIGN_H
#define LEAD_LAG_DESIGN_H

/** @file lead_lag_design.h
 * Lead-Lag Compensator: C(s) = Kc * (T1*s+1)/(alpha*T1*s+1) * (T2*s+1)/(beta*T2*s+1)
 * Lead part: 0 < alpha < 1 (phase boost), Lag part: beta > 1 (DC gain boost).
 * Combined: improves both transient and steady-state response.
 * Ref: Ogata Ch7, Nise Ch9, Dorf & Bishop Ch10, Astrom Ch11
 */
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double Kc;        /* Overall gain */
    double T_lead;    /* Lead time constant */
    double alpha;     /* Lead attenuation, 0<alpha<1 */
    double T_lag;     /* Lag time constant */
    double beta;      /* Lag ratio, beta>1 */
} LeadLagCompensator;

typedef struct {
    double phase_margin_desired;
    double phase_margin_current;
    double gain_crossover;
    double ess_current;
    double ess_desired;
    double dc_gain_required;
    double phase_margin_extra;
} LeadLagDesignSpec;

typedef struct {
    double omega; double magnitude; double magnitude_db;
    double phase_deg; double real; double imag;
} LeadLagFreqPoint;

/* L4 - Design */
int lead_lag_design_from_spec(const LeadLagDesignSpec *spec, LeadLagCompensator *comp);
int lead_lag_design_direct(double pm_cur, double pm_des,
                           double wc, double ess_cur, double ess_des,
                           double margin, LeadLagCompensator *comp);

/* L5 - Analysis */
LeadLagFreqPoint lead_lag_evaluate(const LeadLagCompensator *comp, double omega);
double lead_lag_magnitude_db(const LeadLagCompensator *comp, double omega);
double lead_lag_phase_deg(const LeadLagCompensator *comp, double omega);
double lead_lag_dc_gain(const LeadLagCompensator *comp);
void lead_lag_pole_zero(const LeadLagCompensator *comp,
                        double *z_lead, double *p_lead,
                        double *z_lag, double *p_lag);
int lead_lag_is_valid(const LeadLagCompensator *comp);

/* L5 - Time Domain */
void lead_lag_step_response(const LeadLagCompensator *comp, double t_final,
                            int num_points, double *t_out, double *y_out);

/* L5 - Bode */
void lead_lag_bode_data(const LeadLagCompensator *comp,
                        double wmin, double wmax, int n,
                        double *wo, double *md, double *ph);

/* L5 - Discretization */
void lead_lag_discretize_tustin(const LeadLagCompensator *comp, double Ts,
                                double b_lead[2], double *a1_lead,
                                double b_lag[2], double *a1_lag);
double lead_lag_apply_digital(const double b_lead[2], double a1_lead,
                              const double b_lag[2], double a1_lag,
                              double input,
                              double *u_lead_prev, double *y_lead_prev,
                              double *u_lag_prev, double *y_lag_prev);

/* L7 - Closed-loop analysis */
double lead_lag_closed_loop_pm(const LeadLagCompensator *comp,
                               double plant_dc, double plant_pole,
                               double wc);
double lead_lag_bandwidth_estimate(const LeadLagCompensator *comp);

#ifdef __cplusplus
}
#endif
#endif
