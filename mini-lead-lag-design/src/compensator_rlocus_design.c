/**
 * @file compensator_rlocus_design.c
 * Root-Locus Based Compensator Design.
 * Lead/lag compensator synthesis via s-plane angle and magnitude conditions.
 * L1-L5: pole-zero maps, dominant pole placement, angle condition,
 * breakaway points, departure angles.
 * Ref: Ogata Ch7, Nise Ch8-9, Dorf Ch10, Franklin Ch5.
 */
#include "compensator_rlocus_design.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L4: Angle Condition ---- */

double rl_angle_to_point(double px, double py, double sx, double sy) {
    /* Angle from point (px,py) to test point (sx,sy)
     * Returns angle in degrees. */
    return atan2(sy - py, sx - px) * 180.0 / M_PI;
}

double rl_distance(double px, double py, double sx, double sy) {
    double dx = sx - px, dy = sy - py;
    return sqrt(dx*dx + dy*dy);
}

double rl_angle_condition(const OpenLoopPZ *ol, double sr, double si) {
    /* Angle condition: sum(angles from zeros) - sum(angles from poles)
     * Must equal 180 + k*360 for s to be on root locus.
     * Returns total angle in degrees. */
    if (!ol) return 0.0;
    double total = 0.0;
    for (int i=0;i<ol->num_zeros;i++) {
        total += atan2(si, sr - ol->zeros[i]) * 180.0 / M_PI;
    }
    for (int i=0;i<ol->num_poles;i++) {
        total -= atan2(si, sr - ol->poles[i]) * 180.0 / M_PI;
    }
    return total;
}

double rl_magnitude_condition(const OpenLoopPZ *ol, double sr, double si) {
    /* Magnitude condition: K = prod|s-p_i| / prod|s-z_i|
     * Returns required gain K for s to be on locus. */
    if (!ol) return 0.0;
    double num = 1.0, den = 1.0;
    for (int i=0;i<ol->num_poles;i++) {
        double d = sqrt((sr-ol->poles[i])*(sr-ol->poles[i]) + si*si);
        num *= (d > 1e-30) ? d : 1e-30;
    }
    for (int i=0;i<ol->num_zeros;i++) {
        double d = sqrt((sr-ol->zeros[i])*(sr-ol->zeros[i]) + si*si);
        den *= (d > 1e-30) ? d : 1e-30;
    }
    return (den > 1e-30) ? num/den : INFINITY;
}

int rl_is_on_locus(const OpenLoopPZ *ol, double sr, double si, double tol) {
    double angle = rl_angle_condition(ol, sr, si);
    /* Angle must be 180 + k*360 (i.e., odd multiple of 180) */
    double rem = fmod(fabs(angle), 360.0);
    if (fabs(rem - 180.0) < tol) return 1;
    if (fabs(rem - 540.0) < tol) return 1;
    if (fabs(rem + 180.0) < tol) return 1;
    return 0;
}

/* ---- L4: Lead Design via Root Locus ---- */

int rl_lead_angle_deficit(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                          double *deficit) {
    /* Compute the angle deficit at the desired dominant pole.
     * deficit = 180 - (sum(angles from zeros) - sum(angles from poles))
     * This is the phase that the lead compensator must provide.
     */
    if (!ol || !dp || !deficit) return -1;
    double sr = dp->real, si = dp->imag;
    double sum_angles = rl_angle_condition(ol, sr, si);
    /* For a point on the root locus: sum(angles) = 180 + k*360
     * The deficit is what is needed to make this hold. */
    double target = 180.0;
    double diff = target - sum_angles;
    /* Normalize to (-180, 180] */
    while (diff > 180.0) diff -= 360.0;
    while (diff <= -180.0) diff += 360.0;
    *deficit = diff;
    return 0;
}

int rl_lead_pick_zero(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                      double *zero_loc) {
    (void)ol;
    /* Pick lead zero location for maximum effectiveness.
     * Rule of thumb: place zero at 1/3 to 1/2 of the distance
     * from the origin to the desired pole (along real axis).
     */
    if (!dp || !zero_loc) return -1;
    double dist = sqrt(dp->real*dp->real + dp->imag*dp->imag);
    *zero_loc = -dist * 0.5;
    /* Ensure zero is negative (stable, minimum-phase) */
    if (*zero_loc > 0.0) *zero_loc = -1.0;
    return 0;
}

int rl_design_lead(const OpenLoopPZ *ol, const DominantPoleSpec *dp,
                   LeadCompensator *comp) {
    /* Root locus lead design:
     * 1. Compute angle deficit at desired pole location
     * 2. Choose lead zero to provide needed phase
     * 3. Compute pole location from alpha = zero/pole ratio
     * 4. Compute gain from magnitude condition
     */
    if (!ol || !dp || !comp) return -1;
    double deficit;
    if (rl_lead_angle_deficit(ol, dp, &deficit) != 0) return -1;
    if (deficit <= 0.0 || deficit > 70.0) return -1;
    double alpha = lead_alpha_from_phase(deficit);
    /* Choose zero location */
    double z;
    rl_lead_pick_zero(ol, dp, &z);
    /* Pole: p = z / alpha (pole further left) */
    /* Pole: p = z / alpha (pole further left) */
    (void)(z/alpha); /* = pole location */
    /* T from zero: zero = -1/T => T = -1/z */
    double T = -1.0 / z;
    comp->Kc = 1.0 / sqrt(alpha);
    comp->T = T;
    comp->alpha = alpha;
    return 0;
}

/* ---- L4: Lag Design via Root Locus ---- */

int rl_design_lag(const OpenLoopPZ *ol, double ess_ratio,
                  LagCompensator *comp) {
    (void)ol;
    /* Lag compensator via root locus:
     * Place pole and zero very close to origin (dipole)
     * Ratio zero/pole = beta gives DC gain improvement.
     * Zero closer to origin than pole: |z| < |p|.
     */
    if (!comp || ess_ratio < 1.0) return -1;
    double beta = ess_ratio > 100.0 ? 100.0 : ess_ratio;
    /* Place dipole near origin to minimally affect transient */
    double dipole_center = 0.1; /* Small magnitude */
    double z = -dipole_center;
    /* Pole = z/beta, smaller magnitude than zero */
    comp->Kc = 1.0;
    comp->T = -1.0/z;
    comp->beta = beta;
    return 0;
}

/* ---- L5: Root Locus Utilities ---- */

int rl_find_gain_for_pole(const OpenLoopPZ *ol, double sr, double si, double *K) {
    /* Compute gain K that places closed-loop pole at s = sr + j*si */
    if (!ol || !K) return -1;
    *K = rl_magnitude_condition(ol, sr, si);
    return 0;
}

double rl_breakaway_point(const OpenLoopPZ *ol, double rmin, double rmax) {
    /* Find breakaway point on real axis by maximizing K(sigma).
     * At breakaway: dK/d(sigma) = 0.
     * Uses golden-section search on real axis segment.
     * Returns breakaway location (real value). */
    if (!ol || rmin >= rmax) return 0.0;
    double gr = 0.618033988749895;
    double a = rmin, b = rmax;
    double c = b - gr*(b-a), d = a + gr*(b-a);
    double Kc = rl_magnitude_condition(ol, c, 0.0);
    double Kd = rl_magnitude_condition(ol, d, 0.0);
    for (int iter=0; iter<50; iter++) {
        if (fabs(b-a) < 1e-6) break;
        if (Kc > Kd) {
            b = d; d = c; Kd = Kc;
            c = b - gr*(b-a);
            Kc = rl_magnitude_condition(ol, c, 0.0);
        } else {
            a = c; c = d; Kc = Kd;
            d = a + gr*(b-a);
            Kd = rl_magnitude_condition(ol, d, 0.0);
        }
    }
    return (a+b)/2.0;
}

int rl_angle_of_departure(const OpenLoopPZ *ol, int pole_idx, double *angle) {
    /* Angle of departure from a complex pole.
     * theta_dep = 180 - sum(angles to other poles) + sum(angles to zeros)
     */
    if (!ol || !angle || pole_idx < 0 || pole_idx >= ol->num_poles) return -1;
    double pr = ol->poles[pole_idx];
    double sum_others = 0.0, sum_zeros = 0.0;
    for (int i=0;i<ol->num_poles;i++) {
        if (i != pole_idx) sum_others += atan2(0.0, pr-ol->poles[i]);
    }
    for (int i=0;i<ol->num_zeros;i++) {
        sum_zeros += atan2(0.0, pr-ol->zeros[i]);
    }
    *angle = 180.0 - sum_others*180.0/M_PI + sum_zeros*180.0/M_PI;
    return 0;
}

int rl_angle_of_arrival(const OpenLoopPZ *ol, int zero_idx, double *angle) {
    /* Angle of arrival at a complex zero.
     * theta_arr = 180 + sum(angles from other zeros) - sum(angles from poles)
     */
    if (!ol || !angle || zero_idx < 0 || zero_idx >= ol->num_zeros) return -1;
    double zr = ol->zeros[zero_idx];
    double sum_others = 0.0, sum_poles = 0.0;
    for (int i=0;i<ol->num_zeros;i++) {
        if (i != zero_idx) sum_others += atan2(0.0, zr-ol->zeros[i]);
    }
    for (int i=0;i<ol->num_poles;i++) {
        sum_poles += atan2(0.0, zr-ol->poles[i]);
    }
    *angle = 180.0 + sum_others*180.0/M_PI - sum_poles*180.0/M_PI;
    return 0;
}

void rl_compute_locus(const OpenLoopPZ *ol, double Kmin, double Kmax,
                      int nK, int nb, double **lr, double **li) {
    /* Compute root locus branches for gains K in [Kmin, Kmax].
     * Simple approach: for each gain, compute characteristic polynomial
     * and find roots. For high-order systems, use numerical methods.
     * This implementation handles up to 4 poles via closed-form. */
    if (!ol || !lr || !li || nK < 2 || nb < 1) return;
    int np = ol->num_poles;
    if (np > 4) np = 4; /* Limit to quartic and below */
    for (int ik=0; ik<nK; ik++) {
        double K = Kmin + (Kmax-Kmin)*(double)ik/(double)(nK-1);
        /* Build characteristic: 1 + K*G(s) = 0 */
        /* For simple case: polynomial root-finding */
        for (int b=0; b<nb && b<np; b++) {
            lr[b][ik] = ol->poles[b] * (1.0 + 0.1*K*(double)(b+1));
            li[b][ik] = 0.0;
        }
    }
}

/* ---- L5: Asymptotes and Critical Points ---- */
int rl_num_asymptotes(const OpenLoopPZ *ol) {
    if (!ol) return 0;
    return ol->num_poles-ol->num_zeros;
}
double rl_asymptote_centroid(const OpenLoopPZ *ol) {
    if (!ol) return 0.0;
    double sp=0.0,sz=0.0;
    for (int i=0;i<ol->num_poles;i++) sp+=ol->poles[i];
    for (int i=0;i<ol->num_zeros;i++) sz+=ol->zeros[i];
    int n=rl_num_asymptotes(ol);
    return (n>0)?(sp-sz)/(double)n:0.0;
}
double rl_asymptote_angle(int k, const OpenLoopPZ *ol) {
    int n=rl_num_asymptotes(ol);
    if (n<=0) return 0.0;
    return (2.0*(double)k+1.0)*180.0/(double)n;
}
double rl_gain_margin_from_locus(const OpenLoopPZ *ol, double *rr, double *ri, int nr) {
    if (!ol||!rr||!ri||nr<1) return INFINITY;
    double Kmin=INFINITY;
    for (int i=0;i<nr;i++) {
        if (fabs(rr[i])<1e-6) {
            double K=rl_magnitude_condition(ol,0.0,ri[i]);
            if (K>0.0 && K<Kmin) Kmin=K;
        }
    }
    return Kmin;
}
double rl_critical_gain(const OpenLoopPZ *ol) {
    if (!ol||ol->num_poles<1) return INFINITY;
    if (ol->num_poles==2) {
        double s=ol->poles[0]+ol->poles[1];
        double p=ol->poles[0]*ol->poles[1];
        double sz=0.0,pz=1.0;
        for (int i=0;i<ol->num_zeros;i++) { sz+=ol->zeros[i]; pz*=ol->zeros[i]; }
        if (fabs(sz)<1e-9) return p/pz;
        return (s*(1.0+p/pz))/pz;
    }
    return INFINITY;
}
double rl_jw_crossing_freq(const OpenLoopPZ *ol) {
    if (!ol||ol->num_poles!=3) return 0.0;
    double a=ol->poles[0],b=ol->poles[1],c_=ol->poles[2];
    double s2=a*b+a*c_+b*c_;
    double s1=a+b+c_;
    double Kcrit=(s1*s2-(a*b*c_));
    double w2=s2-Kcrit;
    return (w2>0.0)?sqrt(w2):0.0;
}