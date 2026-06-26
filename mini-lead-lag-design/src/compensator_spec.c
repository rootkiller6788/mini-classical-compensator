/**
 * @file compensator_spec.c
 * Design Specification Translation and Performance Metrics.
 * Converts between time-domain (PO, ts, tr) and frequency-domain
 * (PM, GM, BW) specs using 2nd-order approximations.
 * L1-L5 coverage. Ref: Ogata Ch10, Nise Ch10, Franklin Ch6.
 */
#include "compensator_spec.h"
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L2/L4: PO <-> Damping Ratio ---- */

double po_to_damping(double po) {
    /* PO = 100*exp(-pi*zeta/sqrt(1-zeta^2))
     * Solve: zeta = -ln(PO/100)/sqrt(pi^2 + ln^2(PO/100))
     * Valid for 0 < PO < 100. PO=0 -> zeta=1. PO=100 -> zeta=0.
     */
    if (po <= 0.0) return 1.0;
    if (po >= 100.0) return 0.0;
    double ln_po = log(po / 100.0);
    return -ln_po / sqrt(M_PI * M_PI + ln_po * ln_po);
}

double damping_to_po(double zeta) {
    /* PO = 100 * exp(-pi*zeta/sqrt(1-zeta^2))
     * Valid for 0 <= zeta < 1.
     */
    if (zeta >= 1.0) return 0.0;
    if (zeta <= 0.0) return 100.0;
    return 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
}

/* ---- L2/L4: Damping <-> Phase Margin ---- */

double damping_to_pm(double zeta) {
    /* Exact formula for 2nd-order unity feedback:
     * PM = atan(2*zeta / sqrt(sqrt(1+4*zeta^4) - 2*zeta^2))
     * Engineering approximation: PM ~ 100*zeta (for zeta < 0.7)
     */
    if (zeta <= 0.0) return 0.0;
    if (zeta >= 1.0) return 90.0;
    double z4 = zeta * zeta * zeta * zeta;
    double inner = sqrt(1.0 + 4.0 * z4) - 2.0 * zeta * zeta;
    if (inner <= 0.0) return 90.0;
    double term = sqrt(inner);
    return atan2(2.0 * zeta, term) * 180.0 / M_PI;
}

double pm_to_damping(double pm) {
    /* Inverse of damping_to_pm().
     * Uses Newton refinement for accuracy.
     * Engineering rule: zeta ~ PM/100.
     */
    if (pm <= 0.0) return 0.0;
    if (pm >= 90.0) return 1.0;
    double zeta = pm / 100.0;  /* Initial estimate */
    /* Newton iteration: zeta_{n+1} = zeta_n - (PM(zeta_n) - PM_target) / PM'(zeta_n) */
    for (int i = 0; i < 4; i++) {
        double pm_est = damping_to_pm(zeta);
        double pm_plus  = damping_to_pm(zeta + 0.0005);
        double pm_minus = damping_to_pm(zeta - 0.0005);
        double deriv = (pm_plus - pm_minus) / 0.001;
        if (fabs(deriv) > 1e-12)
            zeta -= (pm_est - pm) / deriv;
        if (zeta < 0.001) zeta = 0.001;
        if (zeta > 0.999) zeta = 0.999;
    }
    return zeta;
}

/* ---- L4: Bandwidth <-> Settling Time ---- */

double settling_time_to_bandwidth(double ts, double zeta) {
    /* 2nd-order: ts = 4/(zeta*wn). BW ~ wn.
     * Conservative: BW = 4/(zeta*ts)
     */
    if (ts <= 0.0) return INFINITY;
    if (zeta <= 0.0) zeta = 0.1;
    if (zeta > 1.0) zeta = 1.0;
    return 4.0 / (zeta * ts);
}

double bandwidth_to_settling_time(double bw, double zeta) {
    if (bw <= 0.0) return INFINITY;
    if (zeta <= 0.0) zeta = 0.1;
    if (zeta > 1.0) zeta = 1.0;
    return 4.0 / (zeta * bw);
}

double wc_to_bandwidth(double wc, double pm) {
    /* Closed-loop BW relates to open-loop wc:
     * PM ~ 45 deg => BW ~ wc
     * PM ~ 60 deg => BW ~ 1.2*wc
     * PM ~ 30 deg => BW ~ 0.8*wc (with resonant peak)
     */
    if (pm <= 0.0) return wc;
    double factor = 1.0 + (pm - 45.0) / 100.0;
    if (factor < 0.5) factor = 0.5;
    if (factor > 2.0) factor = 2.0;
    return wc * factor;
}

double settling_time_to_wc(double ts, double zeta) {
    /* Conservative design: wc ~ 1.5*wn, wn ~ 4/(zeta*ts)
     * So wc ~ 6/(zeta*ts)
     */
    if (ts <= 0.0) return INFINITY;
    if (zeta <= 0.0) zeta = 0.1;
    return 6.0 / (zeta * ts);
}

double rise_time_to_bandwidth(double tr) {
    /* Empirical: BW * tr ~ 2.2 for many practical systems */
    if (tr <= 0.0) return INFINITY;
    return 2.2 / tr;
}

/* ---- L4: Steady-State Error <-> DC Gain ---- */

double ess_to_dc_gain(double ess) {
    /* Unity feedback, type-0: ess_step = 1/(1+Kp)
     * Kp = 1/ess - 1
     * For type-1: ess_ramp = 1/Kv
     */
    if (ess <= 0.0) return INFINITY;
    if (ess >= 1.0) return 0.0;
    return 1.0 / ess - 1.0;
}

double ess_ratio_to_gain(double ess_cur, double ess_des) {
    /* Gain increase factor: K_new/K_old = ess_cur/ess_des
     * Example: ess_cur=0.2 (20%), ess_des=0.02 (2%)
     *   => need 10x DC gain increase.
     */
    if (ess_des <= 0.0) return INFINITY;
    if (ess_cur <= 0.0) return 1.0;
    return ess_cur / ess_des;
}

/* ---- L4: Performance Estimation from PM/wc ---- */

double cl_damping_from_pm(double pm) {
    return pm_to_damping(pm);
}

double expected_po_from_pm(double pm) {
    double zeta = pm_to_damping(pm);
    return damping_to_po(zeta);
}

double expected_ts_from_pm_wc(double pm, double wc) {
    /* ts = 4/(zeta*wn). Conservative: wn ~ wc */
    double zeta = pm_to_damping(pm);
    if (zeta <= 0.0) zeta = 0.1;
    return 4.0 / (zeta * wc);
}

/* ---- L4: Spec Validation ---- */

int time_spec_is_achievable(const TimeDomainSpec *s) {
    if (!s) return 0;
    if (s->percent_overshoot < 0.0 || s->percent_overshoot > 100.0) return 0;
    if (s->settling_time <= 0.0 || s->rise_time <= 0.0 || s->peak_time <= 0.0) return 0;
    if (s->rise_time >= s->settling_time) return 0;
    if (s->peak_time >= s->settling_time) return 0;
    double zeta = po_to_damping(s->percent_overshoot);
    if (zeta < 0.0 || zeta > 2.0) return 0;
    return 1;
}

int freq_spec_is_achievable(const FrequencyDomainSpec *s) {
    if (!s) return 0;
    if (s->phase_margin <= 0.0 || s->phase_margin > 90.0) return 0;
    if (s->gain_margin < 3.0) return 0; /* <3dB is too marginal */
    if (s->gain_crossover <= 0.0 || s->bandwidth <= 0.0) return 0;
    return 1;
}

/* ---- L5: Specification Conversion Functions ---- */

int time_to_freq_spec(const TimeDomainSpec *td, FrequencyDomainSpec *fd) {
    /* Convert time-domain spec to frequency-domain requirements.
     * Uses 2nd-order approximations as a starting point.
     */
    if (!td || !fd) return -1;
    double zeta = po_to_damping(td->percent_overshoot);
    fd->phase_margin = damping_to_pm(zeta);
    fd->bandwidth = 4.0 / (zeta * td->settling_time);
    fd->gain_crossover = fd->bandwidth * 0.8;
    fd->phase_crossover = fd->gain_crossover * 2.5;
    fd->gain_margin = 10.0; /* Conservative default */
    if (zeta < 0.7 && zeta > 0.0) {
        double denom = 2.0 * zeta * sqrt(1.0 - zeta * zeta);
        fd->resonance_peak = (denom > 0) ? 20.0 * log10(1.0 / denom) : 0.0;
    } else {
        fd->resonance_peak = 0.0;
    }
    fd->resonance_freq = fd->gain_crossover * sqrt(1.0 - 2.0 * zeta * zeta);
    if (fd->resonance_freq < 0.0) fd->resonance_freq = 0.0;
    return 0;
}

int freq_to_time_spec(const FrequencyDomainSpec *fd, TimeDomainSpec *td) {
    /* Convert frequency-domain spec to time-domain estimates. */
    if (!fd || !td) return -1;
    double zeta = pm_to_damping(fd->phase_margin);
    td->damping_ratio = zeta;
    td->percent_overshoot = damping_to_po(zeta);
    td->settling_time = 4.0 / (zeta * fd->bandwidth);
    td->rise_time = 2.2 / fd->bandwidth;
    if (zeta < 1.0) {
        td->peak_time = M_PI / (fd->bandwidth * sqrt(1.0 - zeta * zeta));
    } else {
        td->peak_time = td->settling_time;
    }
    td->steady_state_error = -1.0; /* Unknown from freq spec alone */
    return 0;
}

int merge_specs(const TimeDomainSpec *td, const FrequencyDomainSpec *fd,
                CompensatorSpec *merged) {
    if (!merged) return -1;
    if (td) {
        merged->time = *td;
        merged->use_time_domain = 1;
    } else {
        merged->use_time_domain = 0;
    }
    if (fd) {
        merged->freq = *fd;
        merged->use_freq_domain = 1;
    } else {
        merged->use_freq_domain = 0;
    }
    /* If both provided, prefer frequency-domain as primary for design */
    return 0;
}
