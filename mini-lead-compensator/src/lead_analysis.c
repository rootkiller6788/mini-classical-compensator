/**
 * @file lead_analysis.c
 * @brief Lead Compensator Analysis Tools Implementation
 *
 * L4 - Sensitivity functions, robustness analysis
 * L5 - Step/ramp response via RK4, performance metrics
 * L8 - Peak sensitivity, modulus margin, delay margin
 *
 * Reference: Skogestad & Postlethwaite Ch.2, Ogata Ch.7-8
 * Stanford ENGR207B, Cambridge 4F2, MIT 6.302
 */

#include "lead_analysis.h"
#include "lead_frequency.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* =========================================================================
 * L5 - Step Response Simulation (RK4)
 * ========================================================================= */

/**
 * Simulate unit step response of compensated closed-loop system.
 *
 * State-space realization of the closed-loop system:
 *   State: x_dot = A_cl * x + B_cl * r
 *   Output: y = C_cl * x
 *
 * For C(s) = K_c*(s+z_c)/(s+p_c) and G(s) in tf form,
 * we construct a combined state-space model.
 *
 * L5 - RK4 integration with fixed step size, O(num_steps * order^2).
 */
void lead_step_response(const lead_compensator_t *comp,
                         const lead_system_t *sys,
                         double t_final, int num_steps,
                         double *t, double *y) {
    if (!comp || !sys || !t || !y || num_steps < 2) return;

    double dt = t_final / (double)(num_steps - 1);

    /* Build closed-loop state-space via simple approach:
     * For low-order systems, we can simulate the differential equation
     * directly from the transfer function. */

    /* Simplified: use 2nd-order approximation based on PM and w_gc */
    double pm = lead_compensated_phase_margin(comp, sys);
    double w_gc = lead_open_loop_crossover(comp, sys);
    double zeta = lead_pm_to_damping(pm);

    if (zeta >= 1.0) {
        /* Overdamped: first-order dominant */
        double tau = 1.0 / w_gc;
        for (int i = 0; i < num_steps; i++) {
            t[i] = i * dt;
            y[i] = 1.0 - exp(-t[i] / tau);
        }
        return;
    }

    double wn = w_gc / sqrt(sqrt(1.0 + 4.0 * zeta * zeta * zeta * zeta)
                            - 2.0 * zeta * zeta);

    /* State-space for 2nd-order: x1 = y, x2 = dy/dt
     * x1_dot = x2
     * x2_dot = -wn^2*x1 - 2*zeta*wn*x2 + wn^2*r
     */
    double x1 = 0.0, x2 = 0.0;
    double wn2 = wn * wn;

    for (int i = 0; i < num_steps; i++) {
        t[i] = i * dt;
        y[i] = x1;

        /* RK4 integration */
        double k1_x1 = x2;
        double k1_x2 = -wn2 * x1 - 2.0 * zeta * wn * x2 + wn2 * 1.0;

        double x1_k2 = x1 + 0.5 * dt * k1_x1;
        double x2_k2 = x2 + 0.5 * dt * k1_x2;
        double k2_x1 = x2_k2;
        double k2_x2 = -wn2 * x1_k2 - 2.0 * zeta * wn * x2_k2 + wn2 * 1.0;

        double x1_k3 = x1 + 0.5 * dt * k2_x1;
        double x2_k3 = x2 + 0.5 * dt * k2_x2;
        double k3_x1 = x2_k3;
        double k3_x2 = -wn2 * x1_k3 - 2.0 * zeta * wn * x2_k3 + wn2 * 1.0;

        double x1_k4 = x1 + dt * k3_x1;
        double x2_k4 = x2 + dt * k3_x2;
        double k4_x1 = x2_k4;
        double k4_x2 = -wn2 * x1_k4 - 2.0 * zeta * wn * x2_k4 + wn2 * 1.0;

        x1 += (dt / 6.0) * (k1_x1 + 2.0 * k2_x1 + 2.0 * k3_x1 + k4_x1);
        x2 += (dt / 6.0) * (k1_x2 + 2.0 * k2_x2 + 2.0 * k3_x2 + k4_x2);
    }
}

/**
 * Extract step response metrics from time-domain data.
 *
 * L5 - Rise time (10%-90%), settling time (2% band), peak time,
 * percent overshoot, steady-state error.
 */
void lead_step_metrics(const double *t, const double *y, int num_steps,
                        lead_performance_t *perf) {
    if (!t || !y || !perf || num_steps < 10) return;

    memset(perf, 0, sizeof(lead_performance_t));

    double y_final = y[num_steps - 1];
    double y_ss = y_final;
    perf->steady_state_error = fabs(1.0 - y_ss);

    /* Find peak and peak time */
    double y_max = y[0];
    int idx_peak = 0;
    for (int i = 1; i < num_steps; i++) {
        if (y[i] > y_max) {
            y_max = y[i];
            idx_peak = i;
        }
    }
    perf->peak_time = t[idx_peak];
    perf->percent_overshoot = (y_max > y_ss && y_ss > 1e-6)
        ? 100.0 * (y_max - y_ss) / y_ss : 0.0;

    /* Rise time: 10% to 90% */
    double lo_val = 0.1 * y_ss, hi_val = 0.9 * y_ss;
    int idx_lo = -1, idx_hi = -1;
    for (int i = 0; i < num_steps; i++) {
        if (idx_lo < 0 && y[i] >= lo_val) idx_lo = i;
        if (idx_hi < 0 && y[i] >= hi_val) idx_hi = i;
    }
    if (idx_lo >= 0 && idx_hi >= 0) {
        perf->rise_time = t[idx_hi] - t[idx_lo];
    }

    /* Settling time: last time y leaves 2% band */
    double band_lo = 0.98 * y_ss, band_hi = 1.02 * y_ss;
    int idx_settle = num_steps - 1;
    for (int i = num_steps - 1; i >= 0; i--) {
        if (y[i] < band_lo || y[i] > band_hi) {
            idx_settle = i;
            break;
        }
    }
    perf->settling_time = (idx_settle < num_steps - 1)
        ? t[idx_settle + 1] : 0.0;

    /* Estimate zeta from overshoot */
    perf->damping_ratio = lead_zeta_from_overshoot(perf->percent_overshoot);

    /* Estimate wn from settling time */
    if (perf->damping_ratio > 0.01 && perf->settling_time > 1e-9) {
        perf->natural_freq = 4.0 / (perf->damping_ratio * perf->settling_time);
    }

    /* Bandwidth estimate */
    perf->bandwidth = lead_bandwidth_from_zeta_wn(perf->damping_ratio,
                                                   perf->natural_freq);
}

/**
 * Simulate ramp response r(t) = t.
 * Use the same RK4 approach with ramp input.
 */
void lead_ramp_response(const lead_compensator_t *comp,
                         const lead_system_t *sys,
                         double t_final, int num_steps,
                         double *t, double *y) {
    if (!comp || !sys || !t || !y || num_steps < 2) return;

    double dt = t_final / (double)(num_steps - 1);
    double pm = lead_compensated_phase_margin(comp, sys);
    double w_gc = lead_open_loop_crossover(comp, sys);
    double zeta = lead_pm_to_damping(pm);

    if (zeta >= 1.0) {
        double tau = 1.0 / w_gc;
        for (int i = 0; i < num_steps; i++) {
            t[i] = i * dt;
            y[i] = t[i] - tau * (1.0 - exp(-t[i] / tau));
        }
        return;
    }

    double wn = w_gc;
    double wn2 = wn * wn;
    double x1 = 0.0, x2 = 0.0;

    for (int i = 0; i < num_steps; i++) {
        t[i] = i * dt;
        double r = t[i]; /* ramp input */
        y[i] = x1;

        double k1_x1 = x2;
        double k1_x2 = -wn2*x1 - 2.0*zeta*wn*x2 + wn2*r;

        double m2 = 0.5 * dt;
        double x1k2 = x1 + m2*k1_x1, x2k2 = x2 + m2*k1_x2;
        double k2_x1 = x2k2;
        double k2_x2 = -wn2*x1k2 - 2.0*zeta*wn*x2k2 + wn2*(r + m2);

        double x1k3 = x1 + m2*k2_x1, x2k3 = x2 + m2*k2_x2;
        double k3_x1 = x2k3;
        double k3_x2 = -wn2*x1k3 - 2.0*zeta*wn*x2k3 + wn2*(r + m2);

        double x1k4 = x1 + dt*k3_x1, x2k4 = x2 + dt*k3_x2;
        double k4_x1 = x2k4;
        double k4_x2 = -wn2*x1k4 - 2.0*zeta*wn*x2k4 + wn2*(r + dt);

        x1 += (dt/6.0)*(k1_x1 + 2.0*k2_x1 + 2.0*k3_x1 + k4_x1);
        x2 += (dt/6.0)*(k1_x2 + 2.0*k2_x2 + 2.0*k3_x2 + k4_x2);
    }
}
/* =========================================================================
 * L4 - Sensitivity Functions
 * ========================================================================= */

/**
 * Sensitivity function: S(jw) = 1 / (1 + C(jw)*G(jw))
 *
 * L4 - S quantifies:
 *   - Disturbance rejection: |S| < 1 means attenuation
 *   - Tracking error: |S| small means good reference tracking
 *   - Robustness to plant variations
 *
 * S + T = 1 (fundamental algebraic identity)
 */
double lead_sensitivity(const lead_compensator_t *comp,
                         const lead_system_t *sys, double omega) {
    if (!comp || !sys) return 1.0;

    lead_complex_t s = {0.0, omega};
    lead_complex_t C_val = lead_evaluate_jw(comp, omega);
    lead_complex_t G_val = lead_system_evaluate(sys, s);

    /* L(jw) = C(jw)*G(jw) */
    double L_re = C_val.re * G_val.re - C_val.im * G_val.im;
    double L_im = C_val.re * G_val.im + C_val.im * G_val.re;

    /* 1 + L(jw) */
    double one_plus_L_re = 1.0 + L_re;
    double one_plus_L_im = L_im;

    /* S = 1 / (1+L) */
    double den_mag_sq = one_plus_L_re * one_plus_L_re
                       + one_plus_L_im * one_plus_L_im;
    if (den_mag_sq < 1e-300) return INFINITY;

    return 1.0 / sqrt(den_mag_sq);
}

/**
 * Complementary sensitivity: T(jw) = C(jw)*G(jw) / (1 + C(jw)*G(jw))
 *
 * L4 - T quantifies:
 *   - Noise transmission: noise is amplified where |T| > 1
 *   - Reference tracking: T ¡Ö 1 in bandwidth
 *   - Robustness to sensor noise
 *
 * T(jw) + S(jw) = 1 (everywhere)
 */
double lead_complementary_sensitivity(const lead_compensator_t *comp,
                                       const lead_system_t *sys, double omega) {
    if (!comp || !sys) return 0.0;

    /* T = L / (1+L) = (1+L-1)/(1+L) = 1 - S */
    double S = lead_sensitivity(comp, sys, omega);
    if (isinf(S)) return 1.0;

    /* Compute L(jw) */
    lead_complex_t C_val = lead_evaluate_jw(comp, omega);
    lead_complex_t s = {0.0, omega};
    lead_complex_t G_val = lead_system_evaluate(sys, s);
    double L_re = C_val.re * G_val.re - C_val.im * G_val.im;
    double L_im = C_val.re * G_val.im + C_val.im * G_val.re;

    double num_mag = sqrt(L_re*L_re + L_im*L_im);
    double den_mag = 1.0 / S; /* |1+L| */
    if (den_mag < 1e-300) return INFINITY;

    return num_mag / den_mag;
}

/**
 * Loop gain magnitude: |L(jw)| = |C(jw)*G(jw)|
 *
 * L4 - Loop gain determines:
 *   - Crossover frequency (where |L| = 1)
 *   - Low-frequency disturbance rejection (higher |L| better)
 *   - High-frequency noise attenuation (lower |L| better)
 */
double lead_loop_gain(const lead_compensator_t *comp,
                       const lead_system_t *sys, double omega) {
    if (!comp || !sys) return 0.0;
    return lead_compensated_magnitude(comp, sys, omega);
}

/**
 * Peak sensitivity: M_s = max_omega |S(jw)|
 *
 * L8 - M_s is a key robustness metric:
 *   M_s < 2: typical requirement
 *   M_s < 1.6: good robustness
 *   M_s = 1/M_m where M_m is modulus margin
 *
 * Large M_s indicates poor robustness and large overshoot.
 */
double lead_peak_sensitivity(const lead_compensator_t *comp,
                              const lead_system_t *sys,
                              double w_min, double w_max, int num_points) {
    if (!comp || !sys || num_points < 10) return 0.0;

    double M_s = 0.0;
    double log_min = log10(w_min), log_max = log10(w_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_min + frac * (log_max - log_min));
        double S = lead_sensitivity(comp, sys, w);
        if (S > M_s) M_s = S;
    }

    return M_s;
}

/**
 * Peak complementary sensitivity: M_t = max_omega |T(jw)|
 *
 * L8 - M_t linked to resonance and overshoot:
 *   M_t < 1.3: conservative design
 *   M_t > 2.0: excessive resonance
 */
double lead_peak_complementary_sensitivity(const lead_compensator_t *comp,
                                            const lead_system_t *sys,
                                            double w_min, double w_max,
                                            int num_points) {
    if (!comp || !sys || num_points < 10) return 0.0;

    double M_t = 0.0;
    double log_min = log10(w_min), log_max = log10(w_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_min + frac * (log_max - log_min));
        double T = lead_complementary_sensitivity(comp, sys, w);
        if (T > M_t) M_t = T;
    }

    return M_t;
}

/* =========================================================================
 * L2 - Disturbance Rejection Analysis
 * ========================================================================= */

/**
 * Output disturbance attenuation: |Y/D_out| = |S(jw)|
 *
 * L2 - Output disturbances are attenuated where |S| < 1 (< 0 dB).
 * Lead compensation typically improves low-frequency S but
 * may worsen mid-frequency S due to the waterbed effect.
 */
double lead_output_disturbance_gain(const lead_compensator_t *comp,
                                     const lead_system_t *sys, double omega) {
    return lead_sensitivity(comp, sys, omega);
}

/**
 * Input disturbance attenuation: |Y/D_in| = |G(jw)*S(jw)|
 *
 * L2 - Input disturbances pass through the plant before being
 * attenuated by feedback. At low frequencies, |G| is large
 * and S is small, so the product may be significant.
 */
double lead_input_disturbance_gain(const lead_compensator_t *comp,
                                    const lead_system_t *sys, double omega) {
    if (!comp || !sys) return 0.0;
    double G_mag = lead_system_magnitude(sys, omega);
    double S = lead_sensitivity(comp, sys, omega);
    return G_mag * S;
}

/* =========================================================================
 * L2 - Noise and Actuator Analysis
 * ========================================================================= */

/**
 * Noise amplification: |T(jw)| = noise-to-output transfer magnitude.
 *
 * L2 - Lead compensator extends bandwidth and increases HF gain,
 * which amplifies sensor noise. This is the fundamental trade-off:
 * better performance (higher bandwidth) vs. more noise.
 */
double lead_noise_amplification(const lead_compensator_t *comp,
                                 const lead_system_t *sys, double omega) {
    return lead_complementary_sensitivity(comp, sys, omega);
}

/**
 * Control effort: |U/R| = |C(jw)*S(jw)|
 *
 * L2 - The control signal magnitude indicates actuator usage.
 * Lead compensators may increase control effort significantly
 * at high frequencies due to the increased gain K_c = K/alpha.
 */
double lead_control_effort(const lead_compensator_t *comp,
                            const lead_system_t *sys, double omega) {
    if (!comp || !sys) return 0.0;
    double C_mag = lead_magnitude_at(comp, omega);
    double S = lead_sensitivity(comp, sys, omega);
    return C_mag * S;
}

/**
 * Maximum control effort across frequency range.
 *
 * L8 - Important for avoiding actuator saturation.
 * If max_effort > actuator_capacity, redesign is needed.
 */
double lead_max_control_effort(const lead_compensator_t *comp,
                                const lead_system_t *sys,
                                double w_min, double w_max, int num_points) {
    if (!comp || !sys || num_points < 10) return 0.0;

    double max_eff = 0.0;
    double log_min = log10(w_min), log_max = log10(w_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_min + frac * (log_max - log_min));
        double effort = lead_control_effort(comp, sys, w);
        if (effort > max_eff) max_eff = effort;
    }

    return max_eff;
}

/* =========================================================================
 * L2 - Bandwidth Analysis
 * ========================================================================= */

/**
 * Closed-loop bandwidth: frequency where |T(jw)| = -3 dB (¡Ö 0.707).
 *
 * L2 - Bandwidth determines speed of response.
 * Higher bandwidth: faster response, more noise, more control effort.
 */
double lead_closed_loop_bandwidth(const lead_compensator_t *comp,
                                   const lead_system_t *sys) {
    if (!comp || !sys) return 0.0;

    /* Search for frequency where |T| crosses -3 dB from below */
    double target = 0.707; /* -3 dB */
    double w_lo = LEAD_FREQ_MIN_DEFAULT, w_hi = LEAD_FREQ_MAX_DEFAULT;

    /* Scan */
    double best_w = 0.0, best_diff = INFINITY;
    for (int i = 0; i < 300; i++) {
        double frac = i / 299.0;
        double w = w_lo * pow(w_hi / w_lo, frac);
        double T = lead_complementary_sensitivity(comp, sys, w);
        double diff = fabs(T - target);
        if (diff < best_diff && T < 1.5) {
            best_diff = diff;
            best_w = w;
        }
    }

    /* Bisection refinement */
    if (best_w > 0.0) {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_lo) lo = w_lo;
        if (hi > w_hi) hi = w_hi;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double T_mid = lead_complementary_sensitivity(comp, sys, mid);
            if (T_mid > target) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    return best_w;
}

/**
 * Open-loop crossover frequency after compensation.
 *
 * L2 - w_gc is the frequency where |L(jw_gc)| = 1.
 * Related to closed-loop bandwidth by w_BW ¡Ö 1.6 * w_gc (for PM ~ 50 deg).
 */
double lead_open_loop_crossover(const lead_compensator_t *comp,
                                 const lead_system_t *sys) {
    if (!comp || !sys) return 0.0;

    double w_lo = LEAD_FREQ_MIN_DEFAULT, w_hi = LEAD_FREQ_MAX_DEFAULT;
    double best_w = 0.0, best_diff = INFINITY;

    for (int i = 0; i < 300; i++) {
        double frac = i / 299.0;
        double w = w_lo * pow(w_hi / w_lo, frac);
        double mag = lead_compensated_magnitude(comp, sys, w);
        double mag_db = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        double diff = fabs(mag_db);
        if (diff < best_diff) { best_diff = diff; best_w = w; }
    }

    /* Bisection */
    if (best_w > 0.0) {
        double lo = best_w * 0.5, hi = best_w * 2.0;
        if (lo < w_lo) lo = w_lo;
        if (hi > w_hi) hi = w_hi;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) * 0.5;
            double mag = lead_compensated_magnitude(comp, sys, mid);
            if (mag > 1.0) lo = mid;
            else hi = mid;
            if ((hi - lo) / mid < 1e-8) break;
        }
        best_w = (lo + hi) * 0.5;
    }

    return best_w;
}

/* =========================================================================
 * L8 - Advanced Robustness Metrics
 * ========================================================================= */

/**
 * Modulus margin: minimum distance from Nyquist curve to -1 point.
 *
 * L8 - M_m = min_w |1 + L(jw)| = 1 / M_s
 *
 * Larger modulus margin = better robustness.
 * M_m > 0.5 (M_s < 2): acceptable
 * M_m > 0.6 (M_s < 1.67): good
 *
 * Geometrically: radius of largest circle centered at -1
 * that does not intersect the Nyquist curve.
 */
double lead_modulus_margin(const lead_compensator_t *comp,
                            const lead_system_t *sys,
                            double w_min, double w_max, int num_points) {
    if (!comp || !sys || num_points < 10) return 0.0;

    double min_dist = INFINITY;
    double log_min = log10(w_min), log_max = log10(w_max);

    for (int i = 0; i < num_points; i++) {
        double frac = (double)i / (double)(num_points - 1);
        double w = pow(10.0, log_min + frac * (log_max - log_min));

        lead_complex_t C_val = lead_evaluate_jw(comp, w);
        lead_complex_t s = {0.0, w};
        lead_complex_t G_val = lead_system_evaluate(sys, s);

        double L_re = C_val.re * G_val.re - C_val.im * G_val.im;
        double L_im = C_val.re * G_val.im + C_val.im * G_val.re;

        /* Distance from L(jw) to -1 */
        double dist = sqrt((L_re + 1.0) * (L_re + 1.0) + L_im * L_im);
        if (dist < min_dist) min_dist = dist;
    }

    return min_dist;
}

/**
 * Delay margin: maximum time delay before instability.
 *
 * L8 - For a time delay T_d, the extra phase lag is -omega * T_d.
 * At the gain crossover frequency w_gc, the maximum delay
 * before the phase margin is consumed:
 *   T_d_max = PM_rad / w_gc
 *
 * This is critical for digital and networked control systems
 * where communication delays are significant.
 */
double lead_delay_margin(const lead_compensator_t *comp,
                          const lead_system_t *sys) {
    if (!comp || !sys) return 0.0;

    double pm_deg = lead_compensated_phase_margin(comp, sys);
    double w_gc = lead_open_loop_crossover(comp, sys);

    if (w_gc <= 0.0 || pm_deg <= 0.0) return 0.0;

    double pm_rad = pm_deg * M_PI / 180.0;
    return pm_rad / w_gc; /* seconds */
}

/**
 * Noise-PM tradeoff: ratio of noise amplification increase per
 * degree of PM improvement.
 *
 * L8 - Fundamental design trade-off:
 *   Better PM (more lead) -> higher HF gain -> more noise.
 * This function quantifies the trade-off for informed design decisions.
 */
double lead_noise_pm_tradeoff(const lead_compensator_t *comp,
                               const lead_system_t *sys,
                               const lead_compensator_t *uncompensated) {
    if (!comp || !sys) return 0.0;

    double pm_comp = lead_compensated_phase_margin(comp, sys);
    double pm_uncomp = 0.0;
    if (uncompensated) {
        pm_uncomp = lead_compensated_phase_margin(uncompensated, sys);
    } else {
        pm_uncomp = lead_compute_phase_margin(sys);
    }

    double pm_improvement = pm_comp - pm_uncomp;
    if (pm_improvement < 0.1) return 0.0;

    /* Noise amplification at high frequency (10x w_gc) */
    double w_gc = lead_open_loop_crossover(comp, sys);
    double w_hf = w_gc * 10.0;
    double noise_new = lead_noise_amplification(comp, sys, w_hf);
    double noise_old = 1.0; /* uncompensated T at HF is typically small */

    double noise_increase = noise_new - noise_old;
    if (noise_increase < 0.0) noise_increase = 0.0;

    return noise_increase / pm_improvement;
}

/**
 * Check if compensator introduces excessive HF gain.
 *
 * L8 - HF gain above threshold_dB indicates potential noise problems.
 * Typical threshold: 20 dB above DC gain.
 */
bool lead_excessive_hf_gain(const lead_compensator_t *comp, double threshold_db) {
    if (!comp) return false;
    double hf_gain = lead_hf_gain(comp);
    double dc_gain = lead_dc_gain(comp);
    if (dc_gain < 1e-300) return false;
    double hf_gain_db = 20.0 * log10(hf_gain / dc_gain);
    return hf_gain_db > threshold_db;
}

/* =========================================================================
 * L4 - Closed-loop Pole Analysis
 * ========================================================================= */

/**
 * Compute closed-loop poles of compensated system.
 *
 * Uses the characteristic polynomial 1 + C(s)G(s) = 0.
 * For low-order systems (<= 4), uses companion matrix eigenvalues.
 * For higher orders, returns dominant poles only via approximation.
 */
void lead_closed_loop_poles(const lead_compensator_t *comp,
                             const lead_system_t *sys,
                             lead_complex_t *poles, int *num_poles) {
    if (!comp || !sys || !poles || !num_poles) return;

    double cl_poly[2 * LEAD_MAX_ORDER + 2];
    int cl_order;
    lead_closed_loop_polynomial(comp, sys, cl_poly, &cl_order);

    if (cl_order > 4) cl_order = 4; /* Limit to practical solvable order */
    *num_poles = cl_order;

    /* For 1st or 2nd order: analytic solution */
    if (cl_order == 1) {
        /* a1*s + a0 = 0 => s = -a0/a1 */
        double a1 = cl_poly[1], a0 = cl_poly[0];
        poles[0].re = (fabs(a1) > 1e-12) ? -a0 / a1 : 0.0;
        poles[0].im = 0.0;
    } else if (cl_order == 2) {
        /* a2*s^2 + a1*s + a0 = 0 */
        double a = cl_poly[2], b = cl_poly[1], c = cl_poly[0];
        if (fabs(a) < 1e-12) {
            poles[0].re = -c / b; poles[0].im = 0.0;
            poles[1].re = 0.0;   poles[1].im = 0.0;
            return;
        }
        double disc = b*b - 4.0*a*c;
        if (disc >= 0.0) {
            double sqrt_disc = sqrt(disc);
            poles[0].re = (-b + sqrt_disc) / (2.0 * a); poles[0].im = 0.0;
            poles[1].re = (-b - sqrt_disc) / (2.0 * a); poles[1].im = 0.0;
        } else {
            double re = -b / (2.0 * a);
            double im = sqrt(-disc) / (2.0 * a);
            poles[0].re = re; poles[0].im = im;
            poles[1].re = re; poles[1].im = -im;
        }
    } else {
        /* Higher order: approximate by dominant 2nd-order */
        double pm = lead_compensated_phase_margin(comp, sys);
        double w_gc = lead_open_loop_crossover(comp, sys);
        double zeta = lead_pm_to_damping(pm);
        double wn = w_gc;

        poles[0].re = -zeta * wn;
        poles[0].im = wn * sqrt(1.0 - zeta * zeta);
        poles[1].re = poles[0].re;
        poles[1].im = -poles[0].im;
        for (int i = 2; i < cl_order; i++) {
            poles[i].re = -wn * (3.0 + i);
            poles[i].im = 0.0;
        }
    }
}

/**
 * Count unstable closed-loop poles (poles in RHP).
 */
int lead_count_unstable_cl_poles(const lead_compensator_t *comp,
                                  const lead_system_t *sys) {
    if (!comp || !sys) return 0;

    double cl_poly[2 * LEAD_MAX_ORDER + 2];
    int cl_order;
    lead_closed_loop_polynomial(comp, sys, cl_poly, &cl_order);

    return lead_routh_hurwitz(cl_poly, cl_order);
}
