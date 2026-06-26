#ifndef COMPENSATOR_FREQ_DESIGN_H
#define COMPENSATOR_FREQ_DESIGN_H

/** @file compensator_freq_design.h
 * Frequency-Domain Compensator Design.
 * Bode-based design methods for lead, lag, and lead-lag compensators.
 * Includes: asymptotic Bode plots, phase margin targeting,
 * gain crossover placement, and closed-loop bandwidth shaping.
 * Ref: Ogata Ch7, Nise Ch9, Franklin Ch6, Astrom Ch11
 */
#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include "compensator_spec.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ---- L1: Data Structures ---- */

/** Bode plot data for a transfer function */
typedef struct {
    int num_points;
    double *omega;       /* Frequency array [rad/s] */
    double *mag_db;      /* Magnitude [dB] */
    double *phase_deg;   /* Phase [degrees] */
    double omega_min;
    double omega_max;
} BodeData;

/** First-order plant model: G(s) = Kdc / (tau*s + 1) */
typedef struct {
    double Kdc;  /* DC gain */
    double tau;  /* Time constant (s) */
} FirstOrderPlant;

/** Second-order plant: G(s) = Kdc*wn^2 / (s^2 + 2*zeta*wn*s + wn^2) */
typedef struct {
    double Kdc;   /* DC gain */
    double wn;    /* Natural frequency (rad/s) */
    double zeta;  /* Damping ratio */
} SecondOrderPlant;

/** General plant described by frequency response points */
typedef struct {
    int num_points;
    double *omega;
    double *mag_db;
    double *phase_deg;
} PlantFreqData;

/* ---- L4: Frequency-Domain Design Methods ---- */

/** Design lead compensator to achieve target PM at specified crossover */
int freq_design_lead_pm(double pm_current, double pm_target,
                        double wc_target, double margin_deg,
                        LeadCompensator *comp);

/** Design lead compensator from plant Bode data */
int freq_design_lead_from_bode(const PlantFreqData *plant,
                               double pm_target, double margin_deg,
                               LeadCompensator *comp);

/** Design lag compensator for DC gain improvement */
int freq_design_lag_gain(double gain_ratio, double wc,
                         LagCompensator *comp);

/** Design lag compensator from plant Bode data */
int freq_design_lag_from_bode(const PlantFreqData *plant,
                              double gain_ratio,
                              LagCompensator *comp);

/** Design lead-lag compensator from combined specs */
int freq_design_lead_lag(double pm_current, double pm_target,
                         double wc, double ess_ratio,
                         double margin_deg,
                         LeadLagCompensator *comp);

/* ---- L5: Bode Plot Utilities ---- */

/** Allocate and initialize BodeData */
BodeData *bode_data_create(int num_points);
void bode_data_free(BodeData *bd);

/** Compute Bode data for first-order plant */
void bode_first_order(const FirstOrderPlant *plant,
                      double wmin, double wmax, BodeData *bd);

/** Compute Bode data for second-order plant */
void bode_second_order(const SecondOrderPlant *plant,
                       double wmin, double wmax, BodeData *bd);

/** Find gain crossover frequency from Bode data (magnitude crosses 0 dB) */
double bode_find_gain_crossover(const BodeData *bd);

/** Find phase crossover frequency (phase crosses -180 deg) */
double bode_find_phase_crossover(const BodeData *bd);

/** Compute phase margin from Bode data */
double bode_phase_margin(const BodeData *bd);

/** Compute gain margin from Bode data */
double bode_gain_margin(const BodeData *bd);

/** Compute magnitude at a specific frequency via log interpolation */
double bode_mag_at_freq(const BodeData *bd, double omega);

/** Compute phase at a specific frequency via log interpolation */
double bode_phase_at_freq(const BodeData *bd, double omega);

/** Compute open-loop Bode data: C(s) * G(s) */
void bode_open_loop(const BodeData *plant, const BodeData *comp,
                    BodeData *result);

/** Compute closed-loop Bode data: C*G / (1 + C*G) */
void bode_closed_loop(const BodeData *loop, BodeData *result);

/** Compute sensitivity function: 1 / (1 + C*G) */
void bode_sensitivity(const BodeData *loop, BodeData *result);

/** Compute complementary sensitivity: C*G / (1 + C*G) */
void bode_complementary_sensitivity(const BodeData *loop, BodeData *result);

/* ---- L5: Plant Model Fitting ---- */

/** Fit first-order model to Bode data (low-frequency asymptote) */
int fit_first_order(const PlantFreqData *data, FirstOrderPlant *plant);

/** Fit second-order model to Bode data */
int fit_second_order(const PlantFreqData *data, SecondOrderPlant *plant);

#ifdef __cplusplus
}
#endif
#endif
