/**
 * @file lead_compensator.h
 * @brief Lead Compensator - Core Definitions and Data Structures
 *
 * Lead compensator transfer function:
 *   C(s) = K * (T*s + 1) / (alpha*T*s + 1),   0 < alpha < 1
 *   C(s) = K_c * (s + z_c) / (s + p_c),  |z_c| < |p_c|
 *
 * z_c = 1/T, p_c = 1/(alpha*T), K_c = K/alpha
 *
 * Reference: Ogata Ch.7, Dorf&Bishop Ch.10, Franklin&Powell Ch.6
 * MIT 6.302 �� Stanford ENGR105 �� Berkeley ME132 �� Caltech CDS 110
 * ETH 151-0591 �� Cambridge 3F2 �� Georgia Tech ECE 6550
 * Purdue ECE 602 �� Tsinghua Auto. Control
 */

#ifndef LEAD_COMPENSATOR_H
#define LEAD_COMPENSATOR_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 - Fundamental Constants
 * ========================================================================= */
#define LEAD_MAX_FREQ_POINTS    1024
#define LEAD_MAX_ITERATIONS     200
#define LEAD_MAX_ORDER          32
#define LEAD_TOLERANCE          1e-10
#define LEAD_FREQ_MIN_DEFAULT   1e-3
#define LEAD_FREQ_MAX_DEFAULT   1e6
#define LEAD_PHASE_MARGIN_EPS   5.0
#define LEAD_ALPHA_MIN          0.05
#define LEAD_MAX_SINGLE_STAGE   65.0

/* =========================================================================
 * L1 - Core Enumerations
 * ========================================================================= */
typedef enum {
    LEAD_TYPE_PASSIVE          = 0,
    LEAD_TYPE_ACTIVE           = 1,
    LEAD_TYPE_DIGITAL          = 2,
    LEAD_TYPE_ANALYTIC         = 3
} lead_type_t;

typedef enum {
    LEAD_METHOD_FREQUENCY      = 0,
    LEAD_METHOD_ROOT_LOCUS     = 1,
    LEAD_METHOD_ANALYTIC       = 2,
    LEAD_METHOD_OPTIMIZATION   = 3
} lead_design_method_t;

typedef enum {
    LEAD_PM_STABLE             = 0,
    LEAD_PM_MARGINAL           = 1,
    LEAD_PM_UNSTABLE           = 2,
    LEAD_PM_UNDEFINED          = 3
} lead_pm_status_t;

typedef enum {
    LEAD_STABLE                = 0,
    LEAD_MARGINALLY_STABLE     = 1,
    LEAD_UNSTABLE              = 2
} lead_stability_t;

/* =========================================================================
 * L1 - Core Data Structures
 * ========================================================================= */

/** Complex number (L3 - Complex analysis) */
typedef struct { double re; double im; } lead_complex_t;

/** Rational transfer function G(s) = num(s)/den(s) * gain */
typedef struct {
    double num[LEAD_MAX_ORDER + 1];
    double den[LEAD_MAX_ORDER + 1];
    int    num_order;
    int    den_order;
    double gain;
} lead_tf_t;

/** Lead compensator: C(s)=K_c*(s+z_c)/(s+p_c), |z_c|<|p_c|, alpha=z_c/p_c */
typedef struct {
    double K_c;
    double z_c;
    double p_c;
    double T;
    double alpha;
    double dc_gain;
    lead_type_t type;
} lead_compensator_t;

/** Single frequency response data point */
typedef struct {
    double frequency;
    double magnitude;
    double magnitude_db;
    double phase_rad;
    double phase_deg;
    double real_part;
    double imag_part;
} lead_freq_point_t;

/** Bode plot data */
typedef struct {
    lead_freq_point_t points[LEAD_MAX_FREQ_POINTS];
    int    num_points;
    double freq_min;
    double freq_max;
    bool   log_spacing;
} lead_bode_data_t;

/** Design specifications */
typedef struct {
    double phase_margin_desired;
    double gain_margin_desired;
    double bandwidth_desired;
    double steady_state_error;
    double velocity_error_const;
    double damping_ratio_desired;
    double natural_freq_desired;
    double max_overshoot_pct;
    double settling_time;
    bool   use_frequency_domain;
    lead_design_method_t method;
} lead_specs_t;

/** Design result */
typedef struct {
    lead_compensator_t compensator;
    lead_specs_t       specs;
    double achieved_phase_margin;
    double achieved_gain_margin;
    double crossover_freq_old;
    double crossover_freq_new;
    double phase_lead_added;
    double alpha_actual;
    double T_actual;
    int    iterations;
    bool   converged;
} lead_design_result_t;

/** Closed-loop performance metrics */
typedef struct {
    double rise_time;
    double settling_time;
    double peak_time;
    double percent_overshoot;
    double steady_state_error;
    double damping_ratio;
    double natural_freq;
    double bandwidth;
} lead_performance_t;

/** Plant model for compensation */
typedef struct {
    lead_tf_t  tf;
    bool       has_state_space;
    double     A[LEAD_MAX_ORDER][LEAD_MAX_ORDER];
    double     B[LEAD_MAX_ORDER];
    double     C[LEAD_MAX_ORDER];
    double     D;
    int        order;
    char       name[64];
} lead_system_t;

/* =========================================================================
 * L1-L2 - Core API
 * ========================================================================= */

/** Initialize from pole-zero: C(s)=K_c*(s+z_c)/(s+p_c) */
bool lead_init(lead_compensator_t *comp, double K_c, double z_c,
               double p_c, lead_type_t type);

/** Initialize from K-T-alpha: C(s)=K*(T*s+1)/(alpha*T*s+1) */
bool lead_init_from_KT_alpha(lead_compensator_t *comp, double K, double T,
                              double alpha, lead_type_t type);

void lead_copy(lead_compensator_t *dest, const lead_compensator_t *src);
bool lead_validate(const lead_compensator_t *comp);
void lead_to_transfer_function(const lead_compensator_t *comp, lead_tf_t *tf);
lead_complex_t lead_evaluate(const lead_compensator_t *comp, lead_complex_t s);
lead_complex_t lead_evaluate_jw(const lead_compensator_t *comp, double omega);

/** Phase: phi(w) = atan(wT) - atan(alpha*wT) */
double lead_phase_at(const lead_compensator_t *comp, double omega);

/** Magnitude: |C(jw)| = K*sqrt(1+(wT)^2)/sqrt(1+(alpha*wT)^2) */
double lead_magnitude_at(const lead_compensator_t *comp, double omega);

/* =========================================================================
 * L2 - Phase Lead Analysis
 * ========================================================================= */

/** L4 Theorem: phi_max = asin((1-alpha)/(1+alpha)) */
double lead_phi_max_from_alpha(double alpha);

/** L4 Inverse: alpha = (1-sin(phi_m))/(1+sin(phi_m)) */
double lead_alpha_from_phi_max(double phi_m);

/** L4: omega_m = 1/(T*sqrt(alpha)) */
double lead_omega_max(const lead_compensator_t *comp);

/** L4: |C(j*omega_m)| = K/sqrt(alpha) */
double lead_magnitude_at_omega_max(const lead_compensator_t *comp);

double lead_compensated_phase(const lead_compensator_t *comp,
                               const lead_system_t *system, double omega);
double lead_compensated_magnitude(const lead_compensator_t *comp,
                                   const lead_system_t *system, double omega);
double lead_dc_gain(const lead_compensator_t *comp);
double lead_hf_gain(const lead_compensator_t *comp);
double lead_pm_to_damping(double phase_margin_deg);
double lead_damping_to_pm(double zeta);

/* L3 - Additional zeta/wn/PO/ts conversions */
double lead_overshoot_from_zeta(double zeta);
double lead_zeta_from_overshoot(double po_pct);
double lead_settling_time_2pct(double zeta, double wn);
double lead_wn_from_settling(double zeta, double ts);
double lead_peak_time(double zeta, double wn);
double lead_rise_time(double zeta, double wn);
double lead_damped_frequency(double zeta, double wn);
double lead_bandwidth_from_zeta_wn(double zeta, double wn);

/* L2 - Frequency response and analysis */
void lead_freq_response(const lead_compensator_t *comp, double omega,
                         double *mag_out, double *phase_out);
double lead_efficiency(const lead_compensator_t *comp, double omega);
double lead_phase_bandwidth(const lead_compensator_t *comp, double threshold_ratio);

/* Utility */
void lead_print(const lead_compensator_t *comp);

/* L8 - Advanced: multi-stage and discrete-time (from lead_advanced.c) */
bool lead_design_two_stage(const lead_system_t *system,
                            const lead_specs_t *specs,
                            lead_compensator_t *stage1,
                            lead_compensator_t *stage2);
double lead_cascade_phase(const lead_compensator_t *c1,
                           const lead_compensator_t *c2,
                           double omega);
double lead_cascade_magnitude(const lead_compensator_t *c1,
                               const lead_compensator_t *c2,
                               double omega);
void lead_cascade_transfer_function(const lead_compensator_t *c1,
                                     const lead_compensator_t *c2,
                                     lead_tf_t *tf);
void lead_to_discrete_tustin(const lead_compensator_t *comp, double Ts,
                              double *b0, double *b1, double *a1);
double lead_discrete_update(double b0, double b1, double a1,
                             double x_new, double *x_prev, double *y_prev);
bool lead_design_robust(const lead_system_t *system,
                         const lead_specs_t *specs,
                         lead_design_result_t *result);
double lead_gm_degradation(const lead_compensator_t *comp,
                            const lead_system_t *sys);
bool lead_robust_stability_check(const lead_compensator_t *comp,
                                  const lead_system_t *sys,
                                  double uncertainty_weight);
double lead_half_phase_frequency(const lead_compensator_t *comp, bool upper);

#ifdef __cplusplus
}
#endif

#endif /* LEAD_COMPENSATOR_H */
