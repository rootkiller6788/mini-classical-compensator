/**
 * @file feedforward_filter.h
 * @brief Prefilter and filter design for feedforward signal conditioning.
 *
 * L3 --- Mathematical structures: filter design (Butterworth, Bessel),
 * frequency-domain prefilter specifications.
 *
 * L5 --- Computational methods: FIR/IIR prefilter design,
 * low-pass filters for properness enforcement, notch filters.
 *
 * Course alignment:
 *   MIT 6.302 ? prefilter design
 *   Georgia Tech ECE 6550 ? signal processing for control
 *   Berkeley EE221A ? filter theory
 *   Cambridge 3F2 ? signal conditioning
 */

#ifndef FEEDFORWARD_FILTER_H
#define FEEDFORWARD_FILTER_H

#include "feedforward_core.h"

/* ==========================================================================
 * L1: Filter Types
 * ========================================================================== */

/** FIR (Finite Impulse Response) filter. */
typedef struct {
    double *coeff;      /**< filter coefficients b[0..n_taps-1] */
    int     n_taps;     /**< number of taps (filter order + 1) */
    double *buffer;     /**< delay-line buffer */
    int     buf_idx;    /**< current position in circular buffer */
} FIRFilter;

/** IIR (Infinite Impulse Response) filter (Direct Form II Transposed). */
typedef struct {
    double *b;          /**< numerator coefficients b[0..order] */
    double *a;          /**< denominator coefficients a[1..order] (a[0]=1) */
    int     order;      /**< filter order */
    double *state;      /**< filter states */
} IIRFilter;

/** Notch filter: removes a narrow frequency band.
 *  H(s) = (s^2 + 2*zeta_z*wn*s + wn^2) / (s^2 + 2*zeta_p*wn*s + wn^2)
 */
typedef struct {
    double wn;          /**< notch frequency (rad/s) */
    double zeta_z;      /**< zero damping */
    double zeta_p;      /**< pole damping (zeta_p >> zeta_z for deep notch) */
    double depth;       /**< notch depth at center frequency */
} NotchFilter;

/** Low-pass Butterworth filter design specification. */
typedef struct {
    int     order;      /**< filter order */
    double  cutoff;     /**< cutoff frequency (rad/s) */
    double  num[10];    /**< numerator coefficients */
    double  den[10];    /**< denominator coefficients */
    int     den_order;  /**< actual denominator order */
} ButterworthLP;

/* ==========================================================================
 * L5: FIR Filter Design and Operation
 * ========================================================================== */

/**
 * Initialize an FIR filter with given coefficients.
 *
 * @param fir     Filter to initialize
 * @param coeff   Coefficient array
 * @param n_taps  Number of taps
 */
void fir_init(FIRFilter *fir, const double *coeff, int n_taps);

/**
 * Process one sample through FIR filter.
 *
 * y[n] = sum_{k=0}^{N-1} b[k] * x[n-k]
 *
 * @param fir  FIR filter
 * @param x    Input sample
 * @return Filtered output
 */
double fir_process(FIRFilter *fir, double x);

/**
 * Reset FIR filter state (clear delay line).
 */
void fir_reset(FIRFilter *fir);

/**
 * Free FIR filter memory.
 */
void fir_free(FIRFilter *fir);

/**
 * Design a moving average FIR filter.
 *
 * y[n] = (1/N) * sum_{k=0}^{N-1} x[n-k]
 *
 * Used for smoothing reference trajectories.
 *
 * @param window_size Number of samples to average
 * @param fir         Output FIR filter
 */
void fir_design_moving_average(int window_size, FIRFilter *fir);

/**
 * Design a low-pass FIR filter using the window method.
 *
 * Ideal low-pass: h[n] = sin(wc*(n-M/2)) / (pi*(n-M/2))
 * Windowed with Hamming window.
 *
 * @param cutoff    Cutoff frequency (normalized, 0..0.5)
 * @param n_taps    Number of taps (odd recommended)
 * @param fir       Output FIR filter
 */
void fir_design_lowpass(double cutoff, int n_taps, FIRFilter *fir);

/* ==========================================================================
 * L5: IIR Filter Design and Operation
 * ========================================================================== */

/**
 * Initialize an IIR filter.
 *
 * @param iir    Filter to initialize
 * @param b      Numerator coefficients (order+1 elements)
 * @param a      Denominator coefficients (order elements, a[0]=1 implied)
 * @param order  Filter order
 */
void iir_init(IIRFilter *iir, const double *b, const double *a, int order);

/**
 * Process one sample through IIR filter (Direct Form II Transposed).
 *
 * @param iir  IIR filter
 * @param x    Input sample
 * @return Filtered output
 */
double iir_process(IIRFilter *iir, double x);

/**
 * Reset IIR filter state.
 */
void iir_reset(IIRFilter *iir);

/**
 * Free IIR filter memory.
 */
void iir_free(IIRFilter *iir);

/* ==========================================================================
 * L5: Butterworth Low-Pass Filter
 * ========================================================================== */

/**
 * Design a Butterworth low-pass filter.
 *
 * Maximally flat magnitude response in passband:
 *   |H(jw)|^2 = 1 / (1 + (w/wc)^{2n})
 *
 * @param order   Filter order (1-8)
 * @param cutoff  Cutoff frequency (rad/s)
 * @param lp      Output filter specification
 */
void butterworth_lp_design(int order, double cutoff, ButterworthLP *lp);

/**
 * Get analog Butterworth pole locations.
 *
 * p_k = wc * exp(j*pi*(0.5 + (2k+1)/(2n))) for k = 0..n-1
 *
 * @param order   Filter order
 * @param cutoff  Cutoff frequency (rad/s)
 * @param poles   Output complex poles (real/imag pairs, 2 per complex pair)
 * @param n_poles Output number of poles stored
 */
void butterworth_poles(int order, double cutoff,
                       double *poles, int *n_poles);

/* ==========================================================================
 * L5: Notch Filter Design
 * ========================================================================== */

/**
 * Design a notch filter.
 *
 * H(s) = (s^2 + wn^2) / (s^2 + 2*zeta_p*wn*s + wn^2)
 *
 * At w = wn: |H(j*wn)| = depth
 *
 * @param wn     Notch frequency (rad/s)
 * @param depth  Notch depth (0 = infinite, 0.1 = moderate)
 * @param width  Notch width parameter (larger = wider notch)
 * @param nf     Output notch filter specification
 */
void notch_design(double wn, double depth, double width, NotchFilter *nf);

/**
 * Biquad (second-order section) filter for discrete-time notch.
 *
 * Designed via bilinear transform from analog prototype.
 *
 * @param wn     Notch frequency (rad/s)
 * @param zeta   Damping ratio
 * @param Ts     Sampling period (s)
 * @param b      Output numerator [b0, b1, b2]
 * @param a      Output denominator [a0, a1, a2]
 */
void notch_biquad_design(double wn, double zeta, double Ts,
                         double b[3], double a[3]);

/* ==========================================================================
 * L3: Prefilter for Reference Smoothing
 * ========================================================================== */

/**
 * Design a velocity-limiting prefilter.
 *
 * Limits the rate of change of reference signal:
 *   dr/dt < v_max
 *
 * @param v_max     Maximum velocity (units/s)
 * @param dt        Sample time (s)
 * @param r_in      Input reference array
 * @param n         Number of samples
 * @param r_out     Output filtered reference
 */
void prefilter_velocity_limit(double v_max, double dt,
                              const double *r_in, int n,
                              double *r_out);

/**
 * Design an acceleration-limiting prefilter (S-curve).
 *
 * Limits both velocity and acceleration of reference:
 *   dr/dt < v_max, d^2r/dt^2 < a_max
 *
 * @param v_max   Maximum velocity
 * @param a_max   Maximum acceleration
 * @param dt      Sample time
 * @param target  Target position
 * @param r_out   Output S-curve trajectory
 * @param n_out   Output number of trajectory points
 */
void prefilter_scurve(double v_max, double a_max, double dt,
                      double target, double *r_out, int *n_out);

/**
 * Exponential smoothing prefilter.
 *
 * y[n] = alpha * r[n] + (1 - alpha) * y[n-1]
 *
 * @param alpha   Smoothing factor (0..1, higher = faster response)
 * @param r_in    Input reference
 * @param n       Number of samples
 * @param r_out   Output smoothed reference
 */
void prefilter_exp_smooth(double alpha, const double *r_in,
                          int n, double *r_out);

#endif /* FEEDFORWARD_FILTER_H */
