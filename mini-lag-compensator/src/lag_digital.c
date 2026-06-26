/**
 * @file lag_digital.c
 * @brief Digital implementation of lag compensators
 *
 * Converts continuous-time lag compensators to discrete-time
 * using bilinear (Tustin) transform, zero-order hold, and
 * matched pole-zero methods.
 *
 * L5: Bilinear transform, discrete-time difference equations
 * L8: Anti-windup for digital lag compensators
 *
 * The discrete lag compensator implements:
 *   u[k] = -a1*u[k-1] + b0*e[k] + b1*e[k-1]
 *
 * where coefficients are derived from the continuous-time
 * compensator via the chosen discretization method.
 *
 * Textbook: Ogata, "Discrete-Time Control Systems" (1995)
 *           Franklin/Powell/Workman, "Digital Control of Dynamic Systems"
 */

#include "lag_compensator.h"
#include "lag_types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ==========================================================================
 * L1: Digital lag compensator data structure
 * ========================================================================== */

/**
 * Digital (discrete-time) lag compensator.
 *
 * Implements the difference equation:
 *   u[k] = b0*e[k] + b1*e[k-1] - a1*u[k-1]
 *
 * where e[k] is the error input and u[k] is the control output.
 *
 * Transfer function in z-domain:
 *   D(z) = (b0 + b1*z^{-1}) / (1 + a1*z^{-1})
 */
typedef struct {
    double b0;         /**< feedforward coefficient for e[k] */
    double b1;         /**< feedforward coefficient for e[k-1] */
    double a1;         /**< feedback coefficient for u[k-1] (note: D(z) numerator is b0+b1*z^-1, denominator is 1+a1*z^-1) */
    double e_prev;     /**< previous error e[k-1] for state memory */
    double u_prev;     /**< previous output u[k-1] for state memory */
    double Ts;         /**< sampling period (seconds) */
} LagDigital;

/* ==========================================================================
 * L5: Bilinear (Tustin) transform discretization
 * ========================================================================== */

/**
 * Convert continuous lag compensator to discrete using bilinear transform.
 *
 * The bilinear transform maps: s = (2/Ts) * (z-1)/(z+1)
 *
 * G_c(s) = Kc * (T*s + 1) / (beta*T*s + 1)
 *
 * Substituting s:
 *   D(z) = Kc * (T*(2/Ts)*(z-1)/(z+1) + 1) /
 *               (beta*T*(2/Ts)*(z-1)/(z+1) + 1)
 *
 * Let alpha = 2*T/Ts:
 *   D(z) = Kc * (alpha*(z-1)/(z+1) + 1) /
 *               (beta*alpha*(z-1)/(z+1) + 1)
 *
 * Multiply numerator and denominator by (z+1):
 *   D(z) = Kc * (alpha*(z-1) + (z+1)) /
 *               (beta*alpha*(z-1) + (z+1))
 *        = Kc * ((alpha+1)*z + (1-alpha)) /
 *               ((beta*alpha+1)*z + (1-beta*alpha))
 *
 * Divide through by the denominator's z coefficient:
 *   D(z) = Kc * ((alpha+1)*z + (1-alpha)) / (beta*alpha+1) /
 *               (z + (1-beta*alpha)/(beta*alpha+1))
 *
 * So: b0 = Kc*(alpha+1)/(beta*alpha+1)
 *     b1 = Kc*(1-alpha)/(beta*alpha+1)
 *     a1 = (1-beta*alpha)/(beta*alpha+1)
 *
 * The difference equation is:
 *   u[k] + a1*u[k-1] = b0*e[k] + b1*e[k-1]
 *   u[k] = b0*e[k] + b1*e[k-1] - a1*u[k-1]
 *
 * Theorem (Bilinear transform properties):
 *   - Maps LHP s-plane to unit circle in z-plane (stability preserved)
 *   - Maps jw axis to unit circle (frequency warping occurs)
 *   - Pre-warping may be needed for critical frequencies
 *
 * Complexity: O(1)
 *
 * @param lag  continuous-time lag compensator
 * @param Ts   sampling period (seconds)
 * @param[out] dig  initialized digital compensator
 * @return 0 on success, negative on error
 */
int lag_to_digital_bilinear(const LagCompensator *lag, double Ts,
                              LagDigital *dig) {
    if (!lag || !dig || Ts <= 0) return -1;

    double alpha = 2.0 * lag->T / Ts;
    double beta_alpha = lag->beta * alpha;

    double denom = beta_alpha + 1.0;
    if (fabs(denom) < 1e-15) return -2;

    dig->b0 = lag->Kc * (alpha + 1.0) / denom;
    dig->b1 = lag->Kc * (1.0 - alpha) / denom;
    dig->a1 = (1.0 - beta_alpha) / denom;

    dig->e_prev = 0.0;
    dig->u_prev = 0.0;
    dig->Ts = Ts;

    return 0;
}

/* ==========================================================================
 * L5: Zero-order hold (ZOH) discretization
 * ========================================================================== */

/**
 * Convert to discrete using zero-order hold method.
 *
 * For a continuous TF G_c(s), the ZOH equivalent is:
 *   D(z) = (1 - z^{-1}) * Z{L^{-1}[G_c(s)/s]}
 *
 * For lag compensator G_c(s) = Kc*(Ts+1)/(beta*Ts+1):
 *   G_c(s)/s = Kc*(Ts+1)/(s*(beta*Ts+1))
 *            = Kc/s + Kc*(T-beta*T)/(beta*Ts+1)
 *            = Kc/s + Kc*(1-beta)*T/(beta*Ts+1)
 *
 * Let p_c = -1/(beta*T) (pole), z_c = -1/T (zero).
 *   G_c(s)/s = Kc/s + Kc*(p_c - z_c)/(p_c) * 1/(s - p_c)
 *   Wait, let me use partial fractions properly:
 *
 *   G_c(s) = Kc*(T*s+1)/(beta*T*s+1) = Kc/beta + Kc*(beta-1)/beta * 1/(beta*T*s+1)
 *   G_c(s)/s = Kc/(beta*s) + Kc*(beta-1)/(beta*s*(beta*T*s+1))
 *
 *   Let a = beta*T. Then:
 *   G_c(s)/s = Kc/(beta*s) + Kc*(beta-1)/beta * [1/s - a/(a*s+1)]/a
 *   = Kc/(beta*s) + Kc*(beta-1)/(a*beta)/s - Kc*(beta-1)/beta/(a*s+1)
 *   = Kc/beta * [1/s + (beta-1)/(a*s)] - Kc*(beta-1)/beta * 1/(a*s+1)
 *
 *   Hmm, this is getting messy. Let me use a more direct approach.
 *
 * The ZOH discrete equivalent for a first-order system:
 *   G(s) = K * (s+z) / (s+p)
 *
 * Using step-invariant transform:
 *   D(z) = K * [1 - (z-p_d)/(z-1) * (1 - exp(-p*Ts))/p_d]  ...
 *
 * Actually, let me implement the general formula:
 *
 * For G_c(s) = Kc*(Ts+1)/(beta*Ts+1):
 *   Let z_c = 1/T (zero in numerator, positive), p_c = 1/(beta*T) (pole, positive).
 *   G_c(s) = Kc * (s + z_c) / (s + p_c)
 *
 * Step-invariant (ZOH) transform:
 *   D(z) = Kc * [1 - (z-1)*A/z * ...]
 *
 * The formula for ZOH of (s+a)/(s+b) is:
 *   D(z) = [z*(a/b) + (1 - a/b - exp(-b*Ts))] / [z - exp(-b*Ts)]
 *
 * Equivalent to:
 *   b0 = a/b
 *   b1 = 1 - a/b - exp(-b*Ts)
 *   a1 = -exp(-b*Ts)
 *
 * For our case: a = z_c = 1/T, b = p_c = 1/(beta*T)
 *   a/b = (1/T) / (1/(beta*T)) = beta
 *
 *   b0 = beta
 *   b1 = 1 - beta - exp(-Ts/(beta*T))
 *   a1 = -exp(-Ts/(beta*T))
 *
 * But this gives DC gain = (b0+b1)/(1+a1) = (beta+1-beta-exp(-p_c*Ts))/(1-exp(-p_c*Ts))
 * = (1-exp(-p_c*Ts))/(1-exp(-p_c*Ts)) = 1. This is NOT Kc!
 *
 * We need to scale by Kc: D'(z) = Kc * D(z).
 * But then DC gain = Kc, which matches the continuous DC gain.
 *
 * So for the ZOH equivalent:
 *   b0 = Kc * beta
 *   b1 = Kc * (1 - beta - exp(-Ts/(beta*T)))
 *   a1 = -exp(-Ts/(beta*T))
 *
 * The difference equation: u[k] + a1*u[k-1] = b0*e[k] + b1*e[k-1]
 *
 * Complexity: O(1)
 */
int lag_to_digital_zoh(const LagCompensator *lag, double Ts,
                        LagDigital *dig) {
    if (!lag || !dig || Ts <= 0) return -1;

    double p_c = 1.0 / (lag->beta * lag->T);  /* pole magnitude (positive) */
    double exp_term = exp(-p_c * Ts);

    dig->b0 = lag->Kc * lag->beta;
    dig->b1 = lag->Kc * (1.0 - lag->beta - exp_term);
    dig->a1 = -exp_term;

    /* Also represented as: u[k] = b0*e[k] + b1*e[k-1] - a1*u[k-1]
     * Since a1 is negative, -a1 is positive: u[k] = b0*e[k] + b1*e[k-1] + |a1|*u[k-1] */

    dig->e_prev = 0.0;
    dig->u_prev = 0.0;
    dig->Ts = Ts;

    return 0;
}

/* ==========================================================================
 * L5: Matched pole-zero discretization
 * ========================================================================== */

/**
 * Convert using matched pole-zero method.
 *
 * Maps poles and zeros via: s -> exp(s*Ts)
 *   zero at s = -1/T         -> z_zero = exp(-Ts/T)
 *   pole at s = -1/(beta*T)  -> z_pole = exp(-Ts/(beta*T))
 *
 * D(z) = K_d * (z - z_zero) / (z - z_pole)
 *      = K_d * (1 - z_zero*z^{-1}) / (1 - z_pole*z^{-1})
 *
 * DC gain matching: D(1) = G_c(0) = Kc
 *   K_d * (1 - z_zero) / (1 - z_pole) = Kc
 *   K_d = Kc * (1 - z_pole) / (1 - z_zero)
 *
 * So: b0 = K_d
 *     b1 = -K_d * z_zero
 *     a1 = -z_pole
 *
 * Complexity: O(1)
 */
int lag_to_digital_matched(const LagCompensator *lag, double Ts,
                            LagDigital *dig) {
    if (!lag || !dig || Ts <= 0) return -1;

    double z_zero = exp(-Ts / lag->T);
    double z_pole = exp(-Ts / (lag->beta * lag->T));

    double K_d = lag->Kc * (1.0 - z_pole) / (1.0 - z_zero);

    dig->b0 = K_d;
    dig->b1 = -K_d * z_zero;
    dig->a1 = -z_pole;

    dig->e_prev = 0.0;
    dig->u_prev = 0.0;
    dig->Ts = Ts;

    return 0;
}

/* ==========================================================================
 * L5: Digital compensator step update
 * ========================================================================== */

/**
 * Execute one step of the digital lag compensator.
 *
 * Implements: u[k] = b0*e[k] + b1*e[k-1] - a1*u[k-1]
 *
 * Note: our a1 is stored such that the denominator is 1 + a1*z^{-1}.
 * So the difference equation is: u[k] + a1*u[k-1] = b0*e[k] + b1*e[k-1]
 * Therefore: u[k] = b0*e[k] + b1*e[k-1] - a1*u[k-1]
 *
 * Complexity: O(1)
 *
 * @param dig  digital compensator state
 * @param e_k  current error input
 * @return current control output u[k]
 */
double lag_digital_step(LagDigital *dig, double e_k) {
    if (!dig) return 0.0;

    double u_k = dig->b0 * e_k + dig->b1 * dig->e_prev - dig->a1 * dig->u_prev;

    /* Anti-windup: clamp output to reasonable range */
    if (u_k > 1e6)  u_k = 1e6;
    if (u_k < -1e6) u_k = -1e6;

    /* Update state memory */
    dig->e_prev = e_k;
    dig->u_prev = u_k;

    return u_k;
}

/**
 * Reset the digital compensator state to zero.
 * Complexity: O(1)
 */
void lag_digital_reset(LagDigital *dig) {
    if (dig) {
        dig->e_prev = 0.0;
        dig->u_prev = 0.0;
    }
}

/* ==========================================================================
 * L5: Sample rate selection
 * ========================================================================== */

/**
 * Recommend a sampling period for digital implementation.
 *
 * Rule of thumb: sample at 10-30 times the bandwidth.
 * For a lag compensator, the highest meaningful frequency is
 * approximately 1/T (the zero corner frequency).
 *
 * Ts_recommended = 1 / (30 * w_zero) = T / 30
 *
 * This ensures at least 30 samples per time constant of the fastest dynamics.
 *
 * Complexity: O(1)
 */
double lag_recommend_sample_period(const LagCompensator *lag) {
    if (!lag || lag->T <= 0) return 0.01;
    return lag->T / 30.0;
}

/* ==========================================================================
 * L5: Frequency response pre-warping for bilinear transform
 * ========================================================================== */

/**
 * Compute the pre-warped continuous-time frequency for bilinear transform.
 *
 * The bilinear transform maps continuous frequency w_c to
 * discrete frequency w_d via:
 *   w_d = (2/Ts) * arctan(w_c * Ts / 2)
 *
 * Pre-warping: given desired discrete frequency w_d, find w_c:
 *   w_c = (2/Ts) * tan(w_d * Ts / 2)
 *
 * This is critical when designing digital compensators that
 * must achieve specific phase characteristics at a target frequency.
 *
 * Complexity: O(1)
 */
double lag_pre_warp_frequency(double w_desired, double Ts) {
    if (Ts <= 0) return w_desired;
    return (2.0 / Ts) * tan(w_desired * Ts / 2.0);
}

/* ==========================================================================
 * L5: Continuous-to-discrete conversion helper
 * ========================================================================== */

/**
 * Convert continuous lag compensator to digital using the specified method.
 *
 * Method 0: Bilinear (Tustin) — frequency warping, stability preserved
 * Method 1: Zero-Order Hold — step-invariant, exact for step inputs
 * Method 2: Matched Pole-Zero — exact DC gain, natural mapping
 *
 * @param lag     continuous compensator
 * @param Ts      sampling period
 * @param method  0=bilinear, 1=ZOH, 2=matched
 * @param[out] dig  digital compensator
 * @return 0 on success
 */
int lag_to_digital(const LagCompensator *lag, double Ts, int method,
                    LagDigital *dig) {
    switch (method) {
        case 0: return lag_to_digital_bilinear(lag, Ts, dig);
        case 1: return lag_to_digital_zoh(lag, Ts, dig);
        case 2: return lag_to_digital_matched(lag, Ts, dig);
        default: return -1;
    }
}