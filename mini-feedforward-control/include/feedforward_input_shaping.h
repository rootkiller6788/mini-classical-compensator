/**
 * @file feedforward_input_shaping.h
 * @brief Input shaping techniques for vibration-free motion control.
 *
 * L5 --- Computational methods: ZV, ZVD, EI, multi-mode shapers.
 * Input shaping convolves a sequence of impulses with the reference
 * command to cancel residual vibration in flexible systems.
 *
 * Core principle (Singer & Seering, 1990):
 *   For a second-order underdamped mode with natural frequency wn
 *   and damping zeta, a sequence of N impulses with amplitudes A_i
 *   and times t_i achieves zero residual vibration if:
 *     sum(A_i * exp(zeta*wn*t_i) * cos(wd*t_i)) = 0
 *     sum(A_i * exp(zeta*wn*t_i) * sin(wd*t_i)) = 0
 *   where wd = wn * sqrt(1 - zeta^2).
 *
 * Course alignment:
 *   MIT 6.302 ? flexible structure control
 *   Georgia Tech AE 6530 ? spacecraft vibration suppression
 *   Berkeley ME232 ? motion control with input shaping
 *   Cambridge 4F2 ? vibration control
 */

#ifndef FEEDFORWARD_INPUT_SHAPING_H
#define FEEDFORWARD_INPUT_SHAPING_H

#include "feedforward_core.h"

/* ==========================================================================
 * L1: Input Shaper Definitions
 * ========================================================================== */

/** Single impulse in a shaper sequence. */
typedef struct {
    double amplitude;   /**< normalized amplitude (sum = 1.0) */
    double time;        /**< application time (s) */
} Impulse;

/** Input shaper: sequence of impulses that cancel vibration. */
typedef struct {
    Impulse *impulses;  /**< array of impulses */
    int      n_imp;     /**< number of impulses */
    double   total_time;/**< total shaper duration (s) */
    double   wn;         /**< natural frequency designed for (rad/s) */
    double   zeta;       /**< damping ratio designed for */
} InputShaper;

/** Type of input shaper. */
typedef enum {
    SHAPER_ZV,          /**< Zero-Vibration: 2 impulses */
    SHAPER_ZVD,         /**< Zero-Vibration-Derivative: 3 impulses */
    SHAPER_EI,          /**< Extra-Insensitive: 3 impulses */
    SHAPER_TWO_MODE,    /**< Two-mode shaper: 4-5 impulses */
    SHAPER_NEGATIVE,    /**< Negative impulse shaper (faster) */
    SHAPER_SPECIFIED    /**< Custom specified impulse sequence */
} ShaperType;

/* ==========================================================================
 * L5: Input Shaper Design
 * ========================================================================== */

/**
 * Design a Zero-Vibration (ZV) shaper.
 *
 * The simplest shaper: 2 impulses that cancel a single mode exactly.
 *   A1 = 1/(1+K),  A2 = K/(1+K)
 *   t1 = 0,       t2 = pi / wd
 *   where K = exp(-zeta*pi / sqrt(1-zeta^2)), wd = wn*sqrt(1-zeta^2)
 *
 * Theorem (Smith 1958, Singer & Seering 1990):
 *   The ZV shaper produces zero residual vibration if the system
 *   natural frequency and damping match the design values exactly.
 *
 * @param wn     Natural frequency (rad/s)
 * @param zeta   Damping ratio
 * @param shaper Output shaper (caller frees with shaper_free)
 */
void shaper_design_zv(double wn, double zeta, InputShaper *shaper);

/**
 * Design a Zero-Vibration-Derivative (ZVD) shaper.
 *
 * 3 impulses. More robust to frequency modeling errors than ZV.
 * Adds derivative constraint: dV/dw = 0 at design frequency.
 *   A1 = 1/(1+2K+K^2), A2 = 2K/(1+2K+K^2), A3 = K^2/(1+2K+K^2)
 *   t1 = 0, t2 = pi/wd, t3 = 2*pi/wd
 *
 * @param wn     Natural frequency (rad/s)
 * @param zeta   Damping ratio
 * @param shaper Output shaper
 */
void shaper_design_zvd(double wn, double zeta, InputShaper *shaper);

/**
 * Design an Extra-Insensitive (EI) shaper.
 *
 * 3 impulses. Designed to keep vibration below a specified threshold V_tol
 * over a wider frequency range than ZVD.
 *
 * @param wn     Natural frequency (rad/s)
 * @param zeta   Damping ratio
 * @param v_tol  Vibration tolerance (typically 0.05 for 5 percent)
 * @param shaper Output shaper
 */
void shaper_design_ei(double wn, double zeta, double v_tol,
                      InputShaper *shaper);

/**
 * Design a two-mode shaper for systems with two dominant vibration modes.
 *
 * Uses convolution of two single-mode shapers.
 *
 * @param wn1     First mode natural frequency (rad/s)
 * @param zeta1   First mode damping
 * @param wn2     Second mode natural frequency (rad/s)
 * @param zeta2   Second mode damping
 * @param type    Shaper type for each mode (ZV or ZVD)
 * @param shaper  Output combined shaper
 */
void shaper_design_two_mode(double wn1, double zeta1,
                            double wn2, double zeta2,
                            ShaperType type, InputShaper *shaper);

/**
 * Design a negative-input shaper for faster response.
 *
 * Uses negative impulses to shorten shaper duration.
 * Trade-off: possible actuator saturation, excitation of unmodeled modes.
 *
 * @param wn     Natural frequency (rad/s)
 * @param zeta   Damping ratio
 * @param shaper Output shaper
 */
void shaper_design_negative_zv(double wn, double zeta, InputShaper *shaper);

/* ==========================================================================
 * L5: Shaper Application and Analysis
 * ========================================================================== */

/**
 * Convolve input shaper with a reference trajectory.
 *
 * y_shaped(t) = sum_i A_i * r(t - t_i)
 *
 * @param shaper    Input shaper
 * @param ref_in    Input reference array
 * @param n_in      Number of input points
 * @param dt        Time step (s)
 * @param ref_out   Output shaped reference (caller allocates n_in + extra)
 * @param n_out     Output number of points
 */
void shaper_apply(const InputShaper *shaper,
                  const double *ref_in, int n_in, double dt,
                  double *ref_out, int *n_out);

/**
 * Compute residual vibration percentage for a given shaper
 * at a frequency ratio r = w / wn.
 *
 * V(r) = exp(-zeta*wn*t_last) * sqrt(C^2 + S^2)
 * where C = sum(A_i*exp(zeta*wn*t_i)*cos(wd*t_i))
 *       S = sum(A_i*exp(zeta*wn*t_i)*sin(wd*t_i))
 *
 * @param shaper  Input shaper
 * @param r       Frequency ratio w / wn
 * @param zeta    Actual damping ratio
 * @return Vibration amplitude percentage (0..1)
 */
double shaper_residual_vibration(const InputShaper *shaper,
                                 double r, double zeta);

/**
 * Compute the sensitivity curve (vibration vs frequency ratio).
 *
 * @param shaper  Input shaper
 * @param zeta    Damping ratio
 * @param r_min   Minimum frequency ratio
 * @param r_max   Maximum frequency ratio
 * @param npts    Number of evaluation points
 * @param r_out   Output frequency ratios
 * @param v_out   Output vibration percentages
 */
void shaper_sensitivity_curve(const InputShaper *shaper,
                              double zeta, double r_min, double r_max,
                              int npts, double *r_out, double *v_out);

/**
 * Compute shaper robustness: frequency range where V < V_tol.
 *
 * @param shaper   Input shaper
 * @param zeta     Damping ratio
 * @param v_tol    Vibration tolerance
 * @param w_low    Output lower frequency bound
 * @param w_high   Output upper frequency bound
 * @return Width of insensitive frequency range (normalized)
 */
double shaper_insensitivity_range(const InputShaper *shaper,
                                  double zeta, double v_tol,
                                  double *w_low, double *w_high);

/**
 * Compute the time-optimal shaper given constraints.
 *
 * Uses convex optimization approach to find minimum-duration
 * impulse sequence satisfying vibration and actuator constraints.
 *
 * @param wn         Natural frequency (rad/s)
 * @param zeta       Damping ratio
 * @param max_amps   Maximum number of impulses allowed
 * @param v_tol      Vibration tolerance
 * @param amp_limit  Per-impulse amplitude limit
 * @param shaper     Output shaper
 */
void shaper_optimal_min_time(double wn, double zeta,
                             int max_amps, double v_tol,
                             double amp_limit, InputShaper *shaper);

/**
 * Free input shaper memory.
 */
void shaper_free(InputShaper *shaper);

/**
 * Print shaper impulse sequence to string.
 */
void shaper_to_string(const InputShaper *shaper, char *buf, int bufsz);

#endif /* FEEDFORWARD_INPUT_SHAPING_H */
