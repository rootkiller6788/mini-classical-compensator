/**
 * @file lag_types.h
 * @brief Shared type definitions for the lag compensator module
 *
 * L1 --- Definitions: Complex frequency response, controller structures,
 * design specifications, and time-domain response data types.
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems ? compensator design
 *   Stanford ENGR105 Feedback Control ? lag/lead design
 *   Berkeley ME232 Advanced Control ? frequency domain methods
 *   Caltech CDS 110 Introduction to Control ? classical compensators
 *   ETH 151-0591 Control I ? Korrekturglieder (correction elements)
 *   Cambridge 3F2 Systems & Control ? compensator synthesis
 *   Tsinghua ?????? ? ??????
 *
 * Textbook: Ogata, "Modern Control Engineering" (2010), Ch. 7
 *           Franklin, Powell, Emami-Naeini, "Feedback Control of Dynamic Systems"
 */

#ifndef LAG_TYPES_H
#define LAG_TYPES_H

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stddef.h>

/* --------------------------------------------------------------------------
 * L1: Complex number type for frequency-domain analysis
 *     z = re + i * im, representing points in the s-plane
 * -------------------------------------------------------------------------- */

/** Complex number in Cartesian representation: z = re + i*im */
typedef struct {
    double re;  /**< real part sigma (damping in s-plane) */
    double im;  /**< imaginary part omega (frequency in s-plane) */
} LagComplex;

/* --------------------------------------------------------------------------
 * L1: Frequency response data point
 *     H(jw) = M(w) * e^{j*phi(w)}
 * -------------------------------------------------------------------------- */

/** Single frequency response data point (Bode/Nyquist) */
typedef struct {
    double omega;        /**< frequency in rad/s */
    double magnitude;    /**< |G(jw)| */
    double magnitude_db; /**< 20*log10(|G(jw)|) in dB */
    double phase_rad;    /**< angle of G(jw) in radians */
    double phase_deg;    /**< angle of G(jw) in degrees */
    double real_part;    /**< Re{G(jw)} for Nyquist plot */
    double imag_part;    /**< Im{G(jw)} for Nyquist plot */
} LagFreqPoint;

/* --------------------------------------------------------------------------
 * L1: Design specification for lag compensator
 *     Captures the engineering requirements that drive the design
 * -------------------------------------------------------------------------- */

/** Steady-state error specification type */
typedef enum {
    LAG_ESS_STEP,       /**< step input -> position error constant Kp */
    LAG_ESS_RAMP,       /**< ramp input -> velocity error constant Kv */
    LAG_ESS_PARABOLIC   /**< parabolic input -> acceleration error constant Ka */
} LagESSType;

/** Design specification for lag compensator synthesis */
typedef struct {
    double phase_margin_target;      /**< desired phase margin in degrees */
    double gain_margin_target;       /**< desired gain margin in dB */
    double ess_target;               /**< target max steady-state error */
    LagESSType ess_type;             /**< type of steady-state error spec */
    double ess_improvement_factor;   /**< factor beta = error_old / error_new >= 1 */
    double gain_crossover_current;   /**< current uncompensated w_gc in rad/s */
    double dc_gain_required;         /**< required DC gain for ESS spec */
    double settling_time_target;     /**< desired settling time (seconds) */
    double overshoot_target;         /**< desired max overshoot (fraction 0-1) */
    double bandwidth_min;            /**< minimum bandwidth in rad/s */
    double bandwidth_max;            /**< maximum bandwidth in rad/s */
    int    safety_margin_decades;    /**< # decades below w_gc for corner freq */
} LagDesignSpec;

/* --------------------------------------------------------------------------
 * L1: Time-domain simulation results
 * -------------------------------------------------------------------------- */

/** Step response trajectory */
typedef struct {
    double *time;           /**< time vector (seconds) */
    double *output;         /**< output response y(t) */
    double *error;          /**< error e(t) = r(t) - y(t) */
    double *control_signal; /**< control signal u(t) from compensator */
    int     num_points;     /**< number of data points */
    double  final_value;    /**< steady-state final value */
    double  settling_time;  /**< 2% settling time */
    double  rise_time;      /**< 10%-90% rise time */
    double  overshoot_pct;  /**< percent overshoot */
    double  peak_time;      /**< time to first peak */
    double  steady_state_error; /**< steady-state error */
} LagStepResponse;

/** Ramp response trajectory */
typedef struct {
    double *time;           /**< time vector */
    double *output;         /**< output response */
    double *error;          /**< tracking error */
    int     num_points;     /**< number of data points */
    double  velocity_error; /**< steady-state velocity error */
} LagRampResponse;

/* --------------------------------------------------------------------------
 * L2: Transfer function representation
 *     General rational transfer function G(s) = N(s)/D(s)
 * -------------------------------------------------------------------------- */

/** Polynomial coefficient representation */
typedef struct {
    int      order;    /**< polynomial order */
    double  *coeff;    /**< coefficients [a0, a1, ..., an] for sum ai*s^i */
} LagPolynomial;

/** Rational transfer function G(s) = num(s)/den(s) */
typedef struct {
    LagPolynomial numerator;
    LagPolynomial denominator;
    double        dc_gain;     /**< G(0) = num(0)/den(0) */
} LagTransferFunction;

/* --------------------------------------------------------------------------
 * L2: Closed-loop system representation
 * -------------------------------------------------------------------------- */

/** Unity-feedback closed-loop system */
typedef struct {
    LagTransferFunction plant;           /**< plant G(s) */
    LagTransferFunction compensator;     /**< compensator G_c(s) */
    LagTransferFunction open_loop;       /**< G_c(s) * G(s) */
    LagTransferFunction closed_loop;     /**< G_c*G / (1 + G_c*G) */
    LagTransferFunction sensitivity;     /**< S = 1/(1+G_c*G) */
    LagTransferFunction comp_sensitivity;/**< T = G_c*G/(1+G_c*G) */
} LagClosedLoop;

/* --------------------------------------------------------------------------
 * L5: Bode plot data container
 * -------------------------------------------------------------------------- */

/** Complete Bode plot data (magnitude + phase vs frequency) */
typedef struct {
    LagFreqPoint *points;       /**< array of frequency response points */
    int           num_points;   /**< number of points */
    double        freq_min;     /**< minimum frequency (rad/s) */
    double        freq_max;     /**< maximum frequency (rad/s) */
    double        dc_gain_db;   /**< DC gain in dB */
    double        gain_crossover; /**< w where |G| = 1 (0 dB) */
    double        phase_crossover;/**< w where angle G = -180 deg */
    double        phase_margin;  /**< phase margin in degrees */
    double        gain_margin;   /**< gain margin in dB */
} LagBodeData;

/** Nyquist plot data container */
typedef struct {
    LagFreqPoint *points;       /**< Nyquist trajectory points */
    int           num_points;   /**< number of points */
    int           encirclements;/**< encirclements of -1 point */
    int           is_stable;    /**< Nyquist stability verdict */
} LagNyquistData;

/* --------------------------------------------------------------------------
 * L5: Optimization result
 * -------------------------------------------------------------------------- */

/** Numerical optimization result */
typedef struct {
    double optimal_beta;     /**< optimal beta parameter */
    double optimal_T;        /**< optimal time constant T */
    double optimal_Kc;       /**< optimal gain Kc */
    double objective_value;  /**< final objective function value */
    int    iterations;       /**< iterations to converge */
    int    converged;        /**< 1 if converged, 0 otherwise */
} LagOptResult;

/* --------------------------------------------------------------------------
 * L7: DC motor system parameters
 * -------------------------------------------------------------------------- */

/** DC motor model parameters for compensator design application */
typedef struct {
    double R;       /**< armature resistance (Ohm) */
    double L;       /**< armature inductance (H) */
    double Kb;      /**< back-EMF constant (V*s/rad) */
    double Kt;      /**< torque constant (N*m/A) */
    double J;       /**< rotor inertia (kg*m^2) */
    double B;       /**< viscous friction (N*m*s/rad) */
    double rated_speed;  /**< rated speed (rad/s) */
    double rated_torque; /**< rated torque (N*m) */
    double supply_voltage; /**< supply voltage (V) */
} LagDCMotorParams;

#endif /* LAG_TYPES_H */
