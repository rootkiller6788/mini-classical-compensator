#ifndef COMPENSATOR_SPEC_H
#define COMPENSATOR_SPEC_H

/** @file compensator_spec.h
 * Design Specifications for Compensator Synthesis.
 * Translates between time-domain (PO, ts, tr, ess) and frequency-domain
 * (PM, GM, wc, BW) specifications. Provides performance metrics.
 * Ref: Ogata Ch10, Nise Ch10, Dorf & Bishop Ch10, Franklin Ch6
 */
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ---- L1: Specification Structs ---- */

/** Time-domain performance specification */
typedef struct {
    double percent_overshoot;   /* PO in percent (e.g., 15.0 = 15%%) */
    double settling_time;       /* ts (seconds, 2%% criterion) */
    double rise_time;           /* tr (seconds, 10-90%%) */
    double peak_time;           /* tp (seconds) */
    double steady_state_error;  /* ess (absolute or ratio) */
    double damping_ratio;       /* zeta, derived from PO */
} TimeDomainSpec;

/** Frequency-domain performance specification */
typedef struct {
    double phase_margin;        /* PM (degrees) */
    double gain_margin;         /* GM (dB) */
    double gain_crossover;      /* wc (rad/s) */
    double phase_crossover;     /* w180 (rad/s) */
    double bandwidth;           /* closed-loop BW (rad/s) */
    double resonance_peak;      /* Mr (dB), closed-loop */
    double resonance_freq;      /* wr (rad/s) */
} FrequencyDomainSpec;

/** Combined compensator design target */
typedef struct {
    TimeDomainSpec time;
    FrequencyDomainSpec freq;
    int use_time_domain;        /* 1 if time-domain is primary */
    int use_freq_domain;        /* 1 if freq-domain is primary */
} CompensatorSpec;

/* ---- L2/L4: Specification Translation Theorems ---- */

/** PO to damping ratio: zeta = -ln(PO/100)/sqrt(pi^2 + ln^2(PO/100)) */
double po_to_damping(double percent_overshoot);

/** Damping ratio to PO: PO = 100*exp(-pi*zeta/sqrt(1-zeta^2)) */
double damping_to_po(double zeta);

/** Damping ratio to phase margin (approx): PM ~ 100*zeta (for 2nd order) */
double damping_to_pm(double zeta);

/** Phase margin to damping ratio: zeta ~ PM/100 */
double pm_to_damping(double pm);

/** Settling time to bandwidth: BW ~ 4/(zeta*ts) * sqrt(...) */
double settling_time_to_bandwidth(double ts, double zeta);

/** Bandwidth to settling time */
double bandwidth_to_settling_time(double bw, double zeta);

/** Gain crossover to bandwidth: wc ~ BW for moderate PM */
double wc_to_bandwidth(double wc, double pm);

/** Settling time to gain crossover (2nd order approx) */
double settling_time_to_wc(double ts, double zeta);

/** Rise time to bandwidth: BW ~ 2.2/tr (approximate) */
double rise_time_to_bandwidth(double tr);

/** Steady-state error to required DC gain: Kdc = 1/ess - 1 */
double ess_to_dc_gain(double ess);

/** Required gain increase from ess ratio */
double ess_ratio_to_gain(double ess_current, double ess_desired);

/* ---- L4: Performance Metrics ---- */

/** Compute closed-loop damping ratio from phase margin */
double cl_damping_from_pm(double pm);

/** Compute expected percent overshoot from PM (2nd order approx) */
double expected_po_from_pm(double pm);

/** Compute expected settling time from PM and crossover */
double expected_ts_from_pm_wc(double pm, double wc);

/** Validate that a time-domain spec is achievable */
int time_spec_is_achievable(const TimeDomainSpec *spec);

/** Validate frequency-domain spec */
int freq_spec_is_achievable(const FrequencyDomainSpec *spec);

/* ---- L5: Specification Conversion ---- */

/** Convert time-domain spec to frequency-domain requirements */
int time_to_freq_spec(const TimeDomainSpec *td, FrequencyDomainSpec *fd);

/** Convert frequency-domain spec to time-domain estimates */
int freq_to_time_spec(const FrequencyDomainSpec *fd, TimeDomainSpec *td);

/** Merge time and frequency specs, resolving conflicts */
int merge_specs(const TimeDomainSpec *td, const FrequencyDomainSpec *fd,
                CompensatorSpec *merged);

#ifdef __cplusplus
}
#endif
#endif
