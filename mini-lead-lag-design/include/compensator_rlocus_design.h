#ifndef COMPENSATOR_RLOCUS_DESIGN_H
#define COMPENSATOR_RLOCUS_DESIGN_H
/** @file compensator_rlocus_design.h
 * Root-Locus Based Compensator Design.
 * Lead/lag design via angle and magnitude conditions on the s-plane.
 * Dominant pole placement, pole-zero cancellation,
 * and analytical compensator synthesis from root locus.
 * Ref: Ogata Ch7, Nise Ch8-9, Dorf Ch10 */
#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double real;       /* Desired closed-loop pole real part */
    double imag;       /* Desired closed-loop pole imag part */
    double zeta;       /* Damping ratio */
    double wn;         /* Natural frequency (rad/s) */
} DominantPoleSpec;

typedef struct {
    int num_poles;
    double poles[20];  /* Open-loop poles (real values or complex conj) */
    int num_zeros;
    double zeros[20];  /* Open-loop zeros */
    double gain;       /* Open-loop gain K */
} OpenLoopPZ;

/* L4: Angle and Magnitude Conditions */
double rl_angle_condition(const OpenLoopPZ *ol, double s_real, double s_imag);
double rl_magnitude_condition(const OpenLoopPZ *ol, double s_real, double s_imag);
int rl_is_on_locus(const OpenLoopPZ *ol, double s_real, double s_imag, double tol);

/* L4: Lead Design via Root Locus */
int rl_design_lead(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                   LeadCompensator *comp);
int rl_lead_angle_deficit(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                          double *angle_deficit_deg);
int rl_lead_pick_zero(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                      double *zero_loc);

/* L4: Lag Design via Root Locus */
int rl_design_lag(const OpenLoopPZ *ol, double ess_ratio,
                  LagCompensator *comp);

/* L5: Root Locus Utilities */
double rl_angle_to_point(double px, double py, double sx, double sy);
double rl_distance(double px, double py, double sx, double sy);
void rl_compute_locus(const OpenLoopPZ *ol, double K_min, double K_max,
                      int nK, int n_branches, double **locus_real, double **locus_imag);
int rl_find_gain_for_pole(const OpenLoopPZ *ol, double s_real, double s_imag, double *K);
double rl_breakaway_point(const OpenLoopPZ *ol, double range_min, double range_max);
int rl_angle_of_departure(const OpenLoopPZ *ol, int pole_idx, double *angle_deg);
int rl_angle_of_arrival(const OpenLoopPZ *ol, int zero_idx, double *angle_deg);

#ifdef __cplusplus
}
#endif
#endif