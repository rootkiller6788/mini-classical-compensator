/**
 * @file cascade_types.h
 * @brief Core type definitions for cascade control systems
 *
 * L1 --- Definitions: Cascade control structure with primary (outer/master)
 * and secondary (inner/slave) loops. The inner loop provides fast disturbance
 * rejection; the outer loop handles the main process variable setpoint tracking.
 *
 * Cascade control decomposes a complex control problem into two nested loops:
 *   - Outer loop: slow, handles primary variable (product quality)
 *   - Inner loop: fast, handles secondary variable (actuator/subsystem)
 *
 * Fundamental principle: the inner loop must be 3-10x faster than the outer
 * loop for effective cascade operation (bandwidth separation).
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems -- cascade compensation
 *   Stanford ENGR105 Feedback Control -- nested loops
 *   Berkeley ME232 Advanced Control -- cascade architecture
 *   Caltech CDS 110 Introduction to Control -- cascade design
 *   ETH 151-0591 Control I -- Kaskadenregelung
 *   Cambridge 3F2 Systems & Control -- cascade systems
 *   Tsinghua Automatic Control -- cascade control
 *   Georgia Tech ECE 6550 -- nonlinear cascade
 *   Purdue ME 575 -- industrial cascade control
 *
 * Textbook: Seborg, Edgar, Mellichamp & Doyle,
 *           "Process Dynamics and Control" (2017), Ch. 13
 *           Astrom & Hagglund, "Advanced PID Control" (2006), Ch. 7
 *           Marlin, "Process Control" (2000), Ch. 14
 */

#ifndef CASCADE_TYPES_H
#define CASCADE_TYPES_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Polynomial and Transfer Function Types
 * ========================================================================== */

/** Real-coefficient polynomial: p(s) = a0 + a1*s + ... + an*s^n
 *  Coefficients stored as [a0, a1, ..., an] of length degree+1. */
typedef struct {
    int      degree;
    double  *coeff;
} CascadePoly;

/** Continuous-time rational transfer function: G(s) = K * num(s) / den(s) */
typedef struct {
    CascadePoly  num;
    CascadePoly  den;
    double       gain;
    double       delay;
    int          has_delay;
} CascadeTF;

/** Discrete-time transfer function: G(z) = K * num(z) / den(z) */
typedef struct {
    int      num_order;
    int      den_order;
    double  *num_coeff;
    double  *den_coeff;
    double   gain;
    double   Ts;
} CascadeDTF;

/* ==========================================================================
 * L1: PID Controller Types
 * ========================================================================== */

/** PID controller in parallel (ISA) form:
 *  C(s) = Kp + Ki/s + Kd*s/(N*s+1) */
typedef struct {
    double   Kp;
    double   Ki;
    double   Kd;
    double   N;
    double   b;
    double   c;
    double   Ts;
    double   u_min;
    double   u_max;
    double   integrator;
    double   prev_error;
    double   prev_y;
    int      has_antiwindup;
    double   Tt;
} CascadePID;

/* ==========================================================================
 * L1: Cascade Control Structure
 * ========================================================================== */

/** Inner (secondary/slave) loop specification */
typedef struct {
    CascadeTF    inner_process;
    CascadePID   inner_controller;
    CascadeTF    inner_cl;
    double       inner_bandwidth;
    double       inner_rise_time;
    double       inner_settle_time;
    double       inner_overshoot;
    int          inner_is_stable;
    const char  *inner_var_name;
} CascadeInner;

/** Outer (primary/master) loop specification */
typedef struct {
    CascadeTF    outer_process;
    CascadePID   outer_controller;
    CascadeTF    equivalent_plant;
    CascadeTF    outer_cl;
    double       outer_bandwidth;
    double       outer_rise_time;
    double       outer_settle_time;
    double       outer_overshoot;
    int          outer_is_stable;
    const char  *outer_var_name;
} CascadeOuter;

/** Complete cascade control system */
typedef struct {
    CascadeInner  inner;
    CascadeOuter  outer;
    double        bandwidth_ratio;
    double        decoupling_factor;
    int           cascade_active;
    int           bumpless_ready;
    const char   *system_name;
} CascadeSystem;

/* ==========================================================================
 * L2: Anti-Windup and Bumpless Transfer Types
 * ========================================================================== */

typedef enum {
    CASCADE_AW_NONE,
    CASCADE_AW_CLAMPING,
    CASCADE_AW_BACK_CALC,
    CASCADE_AW_INCREMENTAL,
    CASCADE_AW_COMBINED
} CascadeAWMethod;

typedef struct {
    CascadeAWMethod  method;
    double           sat_upper;
    double           sat_lower;
    double           Tt;
    int              is_saturated;
    double           saturated_value;
} CascadeAWState;

/** Bumpless transfer state for mode switching */
typedef struct {
    int      current_mode;
    int      target_mode;
    double   last_good_u;
    double   transition_time;
    double   transition_start;
    int      in_transition;
    double   manual_output;
} CascadeBumpless;

/* ==========================================================================
 * L1: Disturbance and Measurement Types
 * ========================================================================== */

typedef struct {
    CascadeTF    dist_inner_model;
    CascadeTF    dist_outer_model;
    double       dist_magnitude;
    int          dist_is_measurable;
    int          dist_is_stationary;
    const char  *dist_name;
} CascadeDisturbance;

typedef struct {
    CascadeTF    sensor_tf;
    double       noise_std;
    double       bias;
    double       calibration_err;
    double       sample_time;
} CascadeSensor;

/* ==========================================================================
 * L1: Design Specifications
 * ========================================================================== */

typedef struct {
    double   inner_pm_target;
    double   outer_pm_target;
    double   inner_gm_target;
    double   outer_gm_target;
    double   inner_bw_min;
    double   outer_bw_max;
    double   bw_ratio_min;
    double   bw_ratio_optimal;
    double   inner_settle_max;
    double   outer_settle_max;
    double   outer_overshoot_max;
    double   outer_ess_max;
    double   inner_noise_sens_max;
    int      use_freq_domain;
    int      inner_controller_type;
    int      outer_controller_type;
} CascadeDesignSpec;

/* ==========================================================================
 * L5: Frequency Response Types
 * ========================================================================== */

typedef struct {
    double   omega;
    double   magnitude;
    double   magnitude_db;
    double   phase_rad;
    double   phase_deg;
    double   real_part;
    double   imag_part;
} CascadeFreqPoint;

typedef struct {
    CascadeFreqPoint *points;
    int               num_points;
    double            freq_min;
    double            freq_max;
    double            dc_gain_db;
    double            gain_crossover;
    double            phase_crossover;
    double            phase_margin;
    double            gain_margin;
    double            bandwidth;
} CascadeFreqResponse;

/* ==========================================================================
 * L2: Performance Metrics
 * ========================================================================== */

typedef struct {
    double   outer_rise_time;
    double   outer_settle_time;
    double   outer_overshoot;
    double   outer_peak_time;
    double   outer_ess;
    double   outer_decay_ratio;
    double   inner_rise_time;
    double   inner_settle_time;
    double   inner_overshoot;
    double   inner_ess;
    double   dist_rejection_ratio;
    double   cascade_efficiency;
    double   robustness_margin;
} CascadePerformance;

/* ==========================================================================
 * L6: Canonical System Models
 * ========================================================================== */

typedef struct {
    double   R, L, Kb, Kt, J, B, gear_ratio;
} CascadeDCMotor;

typedef struct {
    double   V_r, V_j, rho, Cp, UA, dH, k0, Ea, R_gas;
    double   T_in, T_j_in, F_in, F_j;
} CascadeReactor;

typedef struct {
    double   pipe_length, pipe_diameter, roughness;
    double   fluid_density, fluid_viscosity, valve_Cv;
    double   pump_head, pump_curve_slope;
    double   tank_level, outlet_pressure;
} CascadeFlowPressure;

typedef struct {
    double   tank_area, max_level, min_level;
    double   inflow_nominal, outflow_valve_Cv;
    double   pump_max_flow, level_setpoint;
} CascadeLevelTank;

/* ==========================================================================
 * L8: Advanced Cascade Types
 * ========================================================================== */

typedef struct {
    CascadeTF   inner_imc_filter;
    CascadeTF   outer_imc_filter;
    double      inner_lambda;
    double      outer_lambda;
    int         inner_filter_order;
    int         outer_filter_order;
} CascadeIMCParams;

typedef struct {
    int      num_schedules;
    double  *schedule_var;
    double  *schedule_Kc_i;
    double  *schedule_Ti_i;
    double  *schedule_Kc_o;
    double  *schedule_Ti_o;
    int      current_schedule;
    int      interpolation;
} CascadeGainSchedule;

typedef struct {
    CascadeTF   inner_predictor_model;
    CascadeTF   outer_predictor_model;
    double      inner_delay;
    double      outer_delay;
    int         inner_sp_active;
    int         outer_sp_active;
} CascadeSmithPredictor;

#endif /* CASCADE_TYPES_H */
