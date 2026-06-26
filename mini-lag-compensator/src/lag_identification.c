/**
 * @file lag_identification.c
 * @brief System identification from step/frequency response data
 *
 * L5: Step response identification, frequency response identification
 * L6: DC motor parameter estimation
 * L7: Industrial process FOPDT identification
 *
 * Textbook: Ljung, "System Identification" (1999)
 *           Ogata Ch. 10; Seborg/Edgar/Mellichamp, "Process Dynamics and Control"
 */

#include "lag_identification.h"
#include "lag_compensator.h"
#include "lag_frequency.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ==========================================================================
 * L5: First-order identification from step response
 * ========================================================================== */

int lag_identify_first_order(const LagStepResponse *response,
                              double step_amp,
                              double *K, double *tau) {
    /* Identify K and tau from step response of a first-order system.
     *
     * y(t) = K * step_amp * (1 - exp(-t/tau))
     *
     * Method:
     *   1. K = y_ss / step_amp
     *   2. tau = time at which y(t) = 0.632 * y_ss
     */
    if (!response || !K || !tau || step_amp == 0.0) return -1;
    if (response->num_points < 2) return -2;

    /* Use final_value if available, otherwise last data point */
    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0) {
        y_ss = response->output[response->num_points - 1];
    }
    if (fabs(y_ss) < 1e-10) return -3;

    /* DC gain */
    *K = y_ss / step_amp;

    /* Find time constant: time to reach 63.2% of final value */
    double target = 0.632 * y_ss;
    *tau = -1.0;

    for (int i = 1; i < response->num_points; i++) {
        if (response->output[i] >= target) {
            /* Linear interpolation */
            double y0 = response->output[i-1];
            double y1 = response->output[i];
            double t0 = response->time[i-1];
            double t1 = response->time[i];
            if (fabs(y1 - y0) > 1e-15) {
                *tau = t0 + (target - y0) * (t1 - t0) / (y1 - y0);
            } else {
                *tau = t0;
            }
            break;
        }
    }

    if (*tau < 0) {
        /* Fallback: use initial slope method */
        /* tau = K * step_amp / (dy/dt at t=0) */
        double dy = response->output[1] - response->output[0];
        double dt = response->time[1] - response->time[0];
        if (fabs(dy) > 1e-15 && dt > 0) {
            *tau = y_ss * dt / dy;
        } else {
            *tau = 1.0;  /* default */
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: Second-order identification
 * ========================================================================== */

int lag_identify_second_order(const LagStepResponse *response,
                               double step_amp,
                               double *K, double *zeta, double *wn) {
    /* Identify K, zeta, wn from step response of second-order system.
     *
     * From overshoot M_p and peak time T_p:
     *   zeta = sqrt(ln^2(M_p) / (pi^2 + ln^2(M_p)))
     *   w_n = pi / (T_p * sqrt(1 - zeta^2))
     *   K = y_ss / step_amp
     */
    if (!response || !K || !zeta || !wn || step_amp == 0.0) return -1;
    if (response->num_points < 2) return -2;

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0) {
        y_ss = response->output[response->num_points - 1];
    }

    *K = y_ss / step_amp;

    /* Find peak value and time */
    double y_max = response->output[0];
    double t_peak = response->time[0];
    for (int i = 1; i < response->num_points; i++) {
        if (response->output[i] > y_max) {
            y_max = response->output[i];
            t_peak = response->time[i];
        }
    }

    /* Compute overshoot */
    double M_p = (y_max - y_ss) / y_ss;
    if (M_p < 0.001) {
        /* Essentially no overshoot => heavily damped */
        *zeta = 1.0;
        *wn = 4.0 / (response->settling_time > 0 ? response->settling_time : 1.0);
    } else {
        /* Damping ratio from overshoot */
        double ln_Mp = log(M_p);
        *zeta = sqrt(ln_Mp * ln_Mp / (M_PI * M_PI + ln_Mp * ln_Mp));
        if (*zeta >= 1.0) *zeta = 0.999;

        /* Natural frequency from peak time */
        if (t_peak > 0) {
            *wn = M_PI / (t_peak * sqrt(1.0 - (*zeta) * (*zeta)));
        } else {
            *wn = 1.0;
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: DC gain and time constant identification
 * ========================================================================== */

double lag_identify_dc_gain(const LagStepResponse *response,
                             double step_amp) {
    if (!response || step_amp == 0.0) return 0.0;

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0) {
        /* Use average of last 10% of data */
        int start = (int)(response->num_points * 0.9);
        if (start < 0) start = 0;
        double sum = 0.0;
        int count = 0;
        for (int i = start; i < response->num_points; i++) {
            sum += response->output[i];
            count++;
        }
        y_ss = (count > 0) ? sum / count : response->output[response->num_points - 1];
    }

    return y_ss / step_amp;
}

double lag_identify_time_constant(const LagStepResponse *response) {
    /* 63.2% method for dominant time constant */
    if (!response || response->num_points < 2) return -1.0;

    double y_ss = response->final_value;
    if (y_ss <= 0) {
        /* Use last value */
        y_ss = response->output[response->num_points - 1];
    }
    if (fabs(y_ss) < 1e-10) return -1.0;

    double target = 0.632 * fabs(y_ss);

    for (int i = 1; i < response->num_points; i++) {
        double val = (y_ss > 0) ? response->output[i]
                                : -response->output[i];
        if (val >= target) {
            double v0 = (y_ss > 0) ? response->output[i-1]
                                    : -response->output[i-1];
            double v1 = val;
            double t0 = response->time[i-1];
            double t1 = response->time[i];
            if (fabs(v1 - v0) > 1e-15) {
                return t0 + (target - v0) * (t1 - t0) / (v1 - v0);
            }
            return t0;
        }
    }
    return -1.0;
}

double lag_identify_dead_time(const LagStepResponse *response,
                               double threshold_fraction) {
    /* Identify dead time: time until output exceeds threshold of final */
    if (!response || response->num_points < 2) return -1.0;
    if (threshold_fraction <= 0.0) threshold_fraction = 0.05;

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0)
        y_ss = response->output[response->num_points - 1];
    if (fabs(y_ss) < 1e-10) return 0.0;

    double threshold = threshold_fraction * fabs(y_ss);

    for (int i = 0; i < response->num_points; i++) {
        double val = (y_ss > 0) ? response->output[i]
                                : -response->output[i];
        if (val >= threshold) {
            if (i == 0) return response->time[0];
            /* Interpolate */
            double v0 = (y_ss > 0) ? response->output[i-1]
                                    : -response->output[i-1];
            double v1 = val;
            double t0 = response->time[i-1];
            double t1 = response->time[i];
            if (fabs(v1 - v0) > 1e-15)
                return t0 + (threshold - v0) * (t1 - t0) / (v1 - v0);
            return t0;
        }
    }
    return response->time[response->num_points - 1];
}

/* ==========================================================================
 * L5: Frequency response identification
 * ========================================================================== */

int lag_identify_from_bode(const LagBodeData *bode,
                            double *K, double *corner_freq) {
    if (!bode || !K || !corner_freq || bode->num_points < 2) return -1;

    /* DC gain from lowest frequency point */
    *K = bode->points[0].magnitude;

    /* Find corner frequency: where magnitude drops 3 dB from DC */
    double dc_db = bode->points[0].magnitude_db;
    double target_db = dc_db - 3.0;

    *corner_freq = -1.0;
    for (int i = 1; i < bode->num_points; i++) {
        double db1 = bode->points[i-1].magnitude_db;
        double db2 = bode->points[i].magnitude_db;
        if ((db1 - target_db) * (db2 - target_db) <= 0.0) {
            /* Interpolate */
            double t = (target_db - db1) / (db2 - db1);
            double log_w1 = log10(bode->points[i-1].omega);
            double log_w2 = log10(bode->points[i].omega);
            double log_w = log_w1 + t * (log_w2 - log_w1);
            *corner_freq = pow(10.0, log_w);
            break;
        }
    }

    return (*corner_freq > 0) ? 0 : -2;
}

void lag_identify_error_constants(const LagBodeData *bode,
                                   double *Kp, double *Kv, double *Ka) {
    if (!bode || bode->num_points < 1) {
        if (Kp) *Kp = 0.0;
        if (Kv) *Kv = 0.0;
        if (Ka) *Ka = 0.0;
        return;
    }

    double mag_lf = bode->points[0].magnitude;
    double w_lf = bode->points[0].omega;

    if (Kp) *Kp = mag_lf;
    if (Kv) *Kv = w_lf * mag_lf;
    if (Ka) *Ka = w_lf * w_lf * mag_lf;
}

/* ==========================================================================
 * L6: DC motor identification
 * ========================================================================== */

int lag_identify_dc_motor(const LagStepResponse *response,
                           double voltage_step,
                           LagDCMotorParams *params) {
    if (!response || !params || voltage_step == 0.0) return -1;
    if (response->num_points < 2) return -2;

    /* From speed step response to voltage step:
     *   w_ss = Km * V_step
     *   Km = w_ss / V_step
     *
     *   tau_m = identified time constant
     *
     * From Km and known parameters, we can estimate Kt, J, B.
     * But we typically know R, L, Kb, Kt from motor datasheet.
     * Here we estimate Km and tau_m.
     */

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0)
        y_ss = response->output[response->num_points - 1];

    double Km = y_ss / voltage_step;
    double tau_m = lag_identify_time_constant(response);
    if (tau_m < 0) tau_m = 0.1;  /* default */

    /* Fill in estimated parameters.
     * We cannot uniquely determine all parameters from step response alone
     * without additional measurements. We set the known constants
     * and estimate J based on tau_m. */
    /* Km = Kt/(B*R + Kt*Kb), tau_m = J*R/(B*R + Kt*Kb)
     * If we know Kt=Kb (SI units), R, and B:
     *   J = tau_m * (B*R + Kt*Kb) / R
     *   = tau_m * Kt / (Km * R)  [since Km = Kt/(B*R+Kt*Kb)]
     */

    double Kt = params->Kt;
    double R = params->R;
    if (Kt <= 0) Kt = 0.05;
    if (R <= 0) R = 2.0;

    double B = 0.0;
    if (Km > 0) {
        double denom = Kt / Km;
        B = (denom - Kt * params->Kb) / R;
        if (B < 0) B = 0.00001;
    }

    params->J = tau_m * Kt / (Km * R);
    params->B = B;

    return 0;
}

double lag_identify_electrical_tc(const LagStepResponse *current_response,
                                   const LagDCMotorParams *params) {
    /* Electrical time constant = L/R */
    if (!params || params->R <= 0) return -1.0;
    double tau_from_response = lag_identify_time_constant(current_response);
    if (tau_from_response > 0) return tau_from_response;
    return params->L / params->R;
}

double lag_identify_mechanical_tc(const LagStepResponse *speed_response) {
    return lag_identify_time_constant(speed_response);
}

/* ==========================================================================
 * L7: FOPDT identification
 * ========================================================================== */

int lag_identify_fopdt(const LagStepResponse *response,
                        double step_amp,
                        double *K, double *tau, double *theta) {
    /* Smith's two-point method (1957) for FOPDT identification.
     *
     * Find times t1 (28.3%) and t2 (63.2%):
     *   tau = 1.5 * (t2 - t1)
     *   theta = t2 - tau
     *   K = y_ss / step_amp
     */
    if (!response || !K || !tau || !theta || step_amp == 0.0) return -1;
    if (response->num_points < 3) return -2;

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0)
        y_ss = response->output[response->num_points - 1];
    if (fabs(y_ss) < 1e-10) return -3;

    *K = y_ss / step_amp;

    /* Find t1 (28.3%) and t2 (63.2%) */
    double t1 = -1.0, t2 = -1.0;
    double target1 = 0.283 * fabs(y_ss);
    double target2 = 0.632 * fabs(y_ss);

    for (int i = 1; i < response->num_points; i++) {
        double v0 = (y_ss > 0) ? response->output[i-1]
                                : -response->output[i-1];
        double v1 = (y_ss > 0) ? response->output[i]
                                : -response->output[i];

        if (t1 < 0 && v0 <= target1 && v1 >= target1) {
            if (fabs(v1 - v0) > 1e-15)
                t1 = response->time[i-1] +
                     (target1 - v0) * (response->time[i] - response->time[i-1])
                     / (v1 - v0);
            else
                t1 = response->time[i-1];
        }

        if (t2 < 0 && v0 <= target2 && v1 >= target2) {
            if (fabs(v1 - v0) > 1e-15)
                t2 = response->time[i-1] +
                     (target2 - v0) * (response->time[i] - response->time[i-1])
                     / (v1 - v0);
            else
                t2 = response->time[i-1];
        }

        if (t1 >= 0 && t2 >= 0) break;
    }

    if (t1 < 0 || t2 < 0) {
        /* Fallback: use 63.2% method for first-order (theta=0) */
        *tau = lag_identify_time_constant(response);
        *theta = 0.0;
    } else {
        *tau = 1.5 * (t2 - t1);
        *theta = t2 - *tau;
        if (*theta < 0) *theta = 0.0;
    }

    return 0;
}

int lag_identify_area_method(const LagStepResponse *response,
                              double step_amp,
                              double *K, double *tau, double *theta) {
    /* Area method for noisy step response data.
     *
     * More robust than the tangent or two-point methods.
     *
     * A0 = integral_0^inf (y_ss - y(t)) dt
     * tau = A0 * e / (K * step_amp) = A0 / y_ss
     * theta = (integral up to inflection) - tau
     *
     * Simplified: use the time-delay integral method.
     */
    if (!response || !K || !tau || !theta || step_amp == 0.0) return -1;
    if (response->num_points < 3) return -2;

    double y_ss = response->final_value;
    if (y_ss <= 0 && response->num_points > 0)
        y_ss = response->output[response->num_points - 1];
    if (fabs(y_ss) < 1e-10) return -3;

    *K = y_ss / step_amp;

    /* Compute area A0 = integral of (y_ss - y(t)) dt using trapezoidal rule */
    double area = 0.0;
    for (int i = 1; i < response->num_points; i++) {
        double e0 = y_ss - response->output[i-1];
        double e1 = y_ss - response->output[i];
        if (e0 < 0) e0 = 0.0;
        if (e1 < 0) e1 = 0.0;
        double dt = response->time[i] - response->time[i-1];
        area += 0.5 * (e0 + e1) * dt;
    }

    /* tau + theta = A0 / y_ss */
    double sum_tau_theta = area / y_ss;

    /* For a FOPDT: area / K_step = tau + theta
     * We also know: dy/dt_max * (tau + theta) ~= K_step  (approximately)
     *
     * Without the maximum slope, we use an alternative: estimate tau
     * from the time to reach 63.2% (already includes theta). */
    double t632 = lag_identify_time_constant(response);
    if (t632 < 0) t632 = sum_tau_theta;

    /* t632 = theta + tau (for FOPDT, 63.2% of the first-order part)
     * sum_tau_theta = theta + tau
     * So they should be approximately equal. */
    *tau = sum_tau_theta * 0.7;  /* rough split */
    *theta = sum_tau_theta * 0.3;

    if (*theta < 0) *theta = 0.0;
    if (*tau < 0.01) *tau = 0.01;

    return 0;
}