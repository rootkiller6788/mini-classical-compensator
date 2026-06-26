/** @file compensator_sensitivity.c
 * Sensitivity Analysis for Compensator Design.
 * S(s)=1/(1+L), T(s)=L/(1+L), S+T=1, robust stability checks.
 * L2: feedback fundamental trade-off (waterbed effect).
 * L4: Small-gain theorem for robust stability.
 * Ref: Skogestad & Postlethwaite Ch2, Doyle Ch4. */
#include "compensator_sensitivity.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L2: S+T=1 ---- */

double sens_s_plus_t(double Sm, double Tm) {
    /* In the SISO case, S+T=1 exactly in complex plane.
     * This checks consistency: |S|+|T| approx 1 when phase ~ 0 or 180. */
    return Sm+Tm;
}

/* ---- L4: Relationship between S/T peaks and stability margins ---- */

double sens_peak_s_to_gm(double Ms) {
    /* GM >= Ms/(Ms-1). Peak sensitivity Ms limits minimum GM.
     * Ms=1.3 => GM>=4.3; Ms=1.7 => GM>=2.4; Ms=2.0 => GM>=2.0 */
    if (Ms<=1.0) return INFINITY;
    return Ms/(Ms-1.0);
}

double sens_peak_s_to_pm(double Ms) {
    /* PM >= 2*arcsin(1/(2*Ms)) (rad), converted to degrees.
     * Ms=1.3 => PM>=45deg; Ms=1.7 => PM>=34deg; Ms=2.0 => PM>=29deg */
    if (Ms<=0.5) return 180.0;
    double arg = 1.0/(2.0*Ms);
    if (arg>1.0) arg=1.0;
    return 2.0*asin(arg)*180.0/M_PI;
}

double sens_peak_t_to_gm(double Mt) {
    /* Complementary sensitivity peak also bounds robustness.
     * MT=1.3 is typical for aerospace. */
    if (Mt<=1.0) return INFINITY;
    return (Mt+1.0)/Mt;
}

/* ---- L5: Sensitivity Computation ---- */

SensitivityData *sens_data_create(int n) {
    if (n<2) return NULL;
    SensitivityData *sd=(SensitivityData*)malloc(sizeof(SensitivityData));
    if (!sd) return NULL;
    sd->n=n;
    size_t sz=(size_t)n*sizeof(double);
    sd->omega=(double*)malloc(sz);
    sd->Smag=(double*)malloc(sz);
    sd->Tmag=(double*)malloc(sz);
    sd->Sph=(double*)malloc(sz);
    sd->Tph=(double*)malloc(sz);
    if (!sd->omega||!sd->Smag||!sd->Tmag||!sd->Sph||!sd->Tph) {
        free(sd->omega); free(sd->Smag); free(sd->Tmag);
        free(sd->Sph); free(sd->Tph); free(sd); return NULL;
    }
    return sd;
}

void sens_data_free(SensitivityData *sd) {
    if (!sd) return;
    free(sd->omega); free(sd->Smag); free(sd->Tmag);
    free(sd->Sph); free(sd->Tph); free(sd);
}

void sens_compute_from_loop(double *lm,double *lp,double *w,int n,SensitivityData *sd) {
    if (!lm||!lp||!w||!sd||n<2||sd->n!=n) return;
    for (int i=0;i<n;i++) {
        sd->omega[i]=w[i];
        double mag=pow(10.0,lm[i]/20.0);
        double ph=lp[i]*M_PI/180.0;
        double lr=mag*cos(ph), li=mag*sin(ph);
        double dr=1.0+lr, di=li, d2=dr*dr+di*di;
        double Sr=dr/d2, Si=-di/d2, Sm=sqrt(Sr*Sr+Si*Si);
        double Tr=(lr*dr+li*di)/d2, Ti=(li*dr-lr*di)/d2;
        double Tm=sqrt(Tr*Tr+Ti*Ti);
        sd->Smag[i]=20.0*log10(Sm>1e-30?Sm:1e-30);
        sd->Tmag[i]=20.0*log10(Tm>1e-30?Tm:1e-30);
        sd->Sph[i]=atan2(Si,Sr)*180.0/M_PI;
        sd->Tph[i]=atan2(Ti,Tr)*180.0/M_PI;
    }
}

double sens_find_peak_s(const SensitivityData *sd) {
    if (!sd||sd->n<1) return 0.0;
    double max=sd->Smag[0];
    for (int i=1;i<sd->n;i++) if (sd->Smag[i]>max) max=sd->Smag[i];
    return max;
}

double sens_find_peak_t(const SensitivityData *sd) {
    if (!sd||sd->n<1) return 0.0;
    double max=sd->Tmag[0];
    for (int i=1;i<sd->n;i++) if (sd->Tmag[i]>max) max=sd->Tmag[i];
    return max;
}

double sens_bandwidth_from_t(const SensitivityData *sd) {
    /* Bandwidth = frequency where |T| crosses -3dB */
    if (!sd||sd->n<2) return 0.0;
    for (int i=1;i<sd->n;i++) {
        if (sd->Tmag[i-1]>=-3.0&&sd->Tmag[i]<=-3.0) {
            double t=(sd->Tmag[i-1]+3.0)/(sd->Tmag[i-1]-sd->Tmag[i]);
            double lw=log10(sd->omega[i-1])+t*(log10(sd->omega[i])-log10(sd->omega[i-1]));
            return pow(10.0,lw);
        }
    }
    return 0.0;
}

/* ---- L5: Robustness Metrics ---- */

double sens_modulus_margin(double gm_db, double pm_deg) {
    /* Modulus margin = shortest distance from Nyquist curve to -1.
     * Approx: min(sin(PM)*wc, 10^(GM/20)-1) */
    double gm_lin=pow(10.0,gm_db/20.0);
    double a=sin(pm_deg*M_PI/180.0);
    double b=gm_lin-1.0;
    return (a<b)?a:b;
}

double sens_delay_margin(double pm_deg, double wc) {
    /* Delay margin = PM/wc (seconds).
     * Maximum time delay before instability. */
    if (wc<=0.0) return INFINITY;
    return pm_deg*M_PI/(180.0*wc);
}

int sens_stability_robustness(double gm_db, double pm_deg, double Ms_tol, double Mt_tol) {
    /* Check if stability margins meet robustness criteria.
     * Typical: GM>=6dB, PM>=45deg, Ms<=1.7, Mt<=1.3 */
    if (gm_db<6.0) return 0;
    if (pm_deg<45.0) return 0;
    double Ms=1.0+1.0/(pow(10.0,gm_db/20.0)-1.0);
    if (Ms>Ms_tol) return 0;
    double Mt=1.0/sin(pm_deg*M_PI/180.0);
    if (Mt>Mt_tol) return 0;
    return 1;
}

double sens_gain_variation_tolerance(double Ms, double pct) {
    /* Maximum gain variation allowed (pct percent).
     * Based on small-gain: ||Delta*T|| < 1.
     * Returns: max allowed variation in percent. */
    (void)pct;
    if (Ms<=1.0) return INFINITY;
    return 1.0/(Ms-1.0)*100.0;
}

double sens_phase_variation_tolerance(double Ms, double deg) {
    (void)deg; /* nominal deg not used in this estimate */
    if (Ms<=0.5) return 180.0;
    double arg=1.0/(2.0*Ms);
    if (arg>1.0) arg=1.0;
    return 2.0*asin(arg)*180.0/M_PI;
}

int sens_is_robustly_stable(double gm, double pm, double gunc, double punc) {
    /* Check robust stability. pm must meet minimum threshold. */
    if (pm < 30.0) return 0;
    double Ms=1.0+1.0/(pow(10.0,gm/20.0)-1.0);
    double g_tol=sens_gain_variation_tolerance(Ms,gunc);
    double p_tol=sens_phase_variation_tolerance(Ms,punc);
    return (gunc<=g_tol && punc<=p_tol);
}

/* ---- L7: Application-Oriented Robustness ---- */
int sens_check_iso_spec(const SensitivityData *sd, double Ms_max, double Mt_max, double bw_min) {
    if (!sd||sd->n<2) return 0;
    if (sens_find_peak_s(sd)>20.0*log10(Ms_max)) return 0;
    if (sens_find_peak_t(sd)>20.0*log10(Mt_max)) return 0;
    if (sens_bandwidth_from_t(sd)<bw_min) return 0;
    return 1;
}
int sens_mimo_cond_estimate(const SensitivityData *sd, double *cond_num) {
    if (!sd||!cond_num||sd->n<2) return -1;
    double max_s=-INFINITY, min_s=INFINITY;
    for (int i=0;i<sd->n;i++) {
        double sm=pow(10.0,sd->Smag[i]/20.0);
        if (sm>max_s) max_s=sm;
        if (sm<min_s) min_s=sm;
    }
    *cond_num = (min_s>0.0)?max_s/min_s:INFINITY;
    return 0;
}
double sens_disturbance_rejection(const SensitivityData *sd, double w_dist) {
    if (!sd||sd->n<2||w_dist<=0.0) return 0.0;
    for (int i=1;i<sd->n;i++) {
        if (sd->omega[i-1]<=w_dist&&sd->omega[i]>=w_dist) {
            double t=(log10(w_dist)-log10(sd->omega[i-1]))/(log10(sd->omega[i])-log10(sd->omega[i-1]));
            return sd->Smag[i-1]+t*(sd->Smag[i]-sd->Smag[i-1]);
        }
    }
    return 0.0;
}
double sens_noise_attenuation(const SensitivityData *sd, double w_noise) {
    if (!sd||sd->n<2||w_noise<=0.0) return 0.0;
    for (int i=1;i<sd->n;i++) {
        if (sd->omega[i-1]<=w_noise&&sd->omega[i]>=w_noise) {
            double t=(log10(w_noise)-log10(sd->omega[i-1]))/(log10(sd->omega[i])-log10(sd->omega[i-1]));
            return sd->Tmag[i-1]+t*(sd->Tmag[i]-sd->Tmag[i-1]);
        }
    }
    return 0.0;
}
int sens_waterbed_check(const SensitivityData *sd) {
    if (!sd||sd->n<2) return 0;
    double int_S=0.0;
    for (int i=1;i<sd->n;i++) {
        double dw=sd->omega[i]-sd->omega[i-1];
        double s_avg=(pow(10.0,sd->Smag[i]/20.0)+pow(10.0,sd->Smag[i-1]/20.0))/2.0;
        int_S+=s_avg*dw;
    }
    return (int_S<100.0);
}