/**
 * @file lag_design.c - Lag Compensator Implementation
 * C(s) = Kc*(1+T*s)/(1+beta*T*s), beta > 1
 * L1-L5: structs, theorems, frequency/time analysis, discretization
 * Ref: Ogata Ch7, Nise Ch9, Dorf & Bishop Ch10
 */
#include "lag_design.h"
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

/* ---- L4: Fundamental Formulas ---- */

double lag_attenuation_db(double beta) {
    /* HF attenuation: 20*log10(1/beta) = -20*log10(beta) dB
     * For beta=10, attenuation = -20 dB.
     */
    if (beta <= 1.0) return 0.0;
    return -20.0 * log10(beta);
}

double lag_max_phase_lag_deg(double beta) {
    /* Maximum phase lag of lag compensator:
     * phi_min = -arcsin((beta-1)/(beta+1))
     * Occurs at omega_m = 1/(T*sqrt(beta))
     */
    if (beta <= 1.0) return 0.0;
    double s = (beta - 1.0)/(beta + 1.0);
    return -RAD2DEG(asin(s));
}

double lag_beta_from_attenuation(double atten_db) {
    /* beta = 10^(|atten_db|/20) */
    if (atten_db >= 0.0) return 1.0;
    return pow(10.0, -atten_db/20.0);
}

double lag_beta_from_dc_gain(double gain_ratio) {
    /* DC gain of lag compensator = Kc*beta (if Kc=1: just beta)
     * beta = gain_ratio / Kc
     */
    if (gain_ratio <= 1.0) return 1.0;
    return gain_ratio;
}

/* ---- L4: Lag Design Algorithm ---- */

int lag_design_from_spec(const LagDesignSpec *spec, LagCompensator *comp) {
    /* Ogata lag design:
     * 1. Determine required DC gain increase: ratio = ess_current/ess_desired
     * 2. beta = ratio (or from dc_gain_required)
     * 3. Place 1/T one decade below gain crossover: T = 10/wc
     * 4. Kc = 1 (unity DC gain with lag) or adjust as needed
     */
    if (!spec || !comp) return -1;
    double beta;
    if (spec->dc_gain_required > 1.0) {
        beta = spec->dc_gain_required;
    } else if (spec->ess_current > 0 && spec->ess_desired > 0) {
        beta = spec->ess_current / spec->ess_desired;
    } else {
        beta = 10.0; /* default 10:1 improvement */
    }
    if (beta < 1.0) beta = 1.0;
    if (beta > 100.0) beta = 100.0; /* practical limit */
    double wc = spec->gain_crossover;
    if (wc <= 0.0) wc = 1.0;
    /* Place zero at wc/10 so phase lag at wc is < 6 deg */
    double T = 10.0 / wc;
    comp->Kc = 1.0;
    comp->T = T;
    comp->beta = beta;
    return 0;
}

int lag_design_direct(double ess_cur, double ess_des, double wc,
                      double ratio, LagCompensator *comp) {
    if (!comp) return -1;
    LagDesignSpec s = {ess_des, ess_cur, wc, ratio};
    return lag_design_from_spec(&s, comp);
}

/* ---- L5: Frequency Response ---- */

LagFreqPoint lag_evaluate(const LagCompensator *comp, double omega) {
    /* C(jw) = Kc*(1+jwT)/(1+jw*beta*T) */
    LagFreqPoint r; r.omega = omega;
    if (!comp || !lag_is_valid(comp)) {
        r.magnitude=0; r.magnitude_db=-INFINITY;
        r.phase_deg=0; r.real=0; r.imag=0; return r;
    }
    double wT=omega*comp->T, wbT=omega*comp->beta*comp->T;
    double nr=comp->Kc, ni=comp->Kc*wT;
    double dr=1.0, di=wbT, d2=dr*dr+di*di;
    if (d2<1e-30) {
        r.magnitude=INFINITY; r.magnitude_db=INFINITY;
        r.phase_deg=0; r.real=INFINITY; r.imag=0; return r;
    }
    r.real=(nr*dr+ni*di)/d2; r.imag=(ni*dr-nr*di)/d2;
    double m2=r.real*r.real+r.imag*r.imag;
    r.magnitude=sqrt(m2);
    r.magnitude_db=(m2<1e-30)?-INFINITY:20.0*log10(r.magnitude);
    r.phase_deg=RAD2DEG(atan2(r.imag,r.real));
    return r;
}

double lag_magnitude_db(const LagCompensator *c, double w)
{ return lag_evaluate(c,w).magnitude_db; }
double lag_phase_deg(const LagCompensator *c, double w)
{ return lag_evaluate(c,w).phase_deg; }
double lag_dc_gain(const LagCompensator *c)
{ return c?c->Kc*c->beta:0.0; }
double lag_hf_gain(const LagCompensator *c)
{ return c?c->Kc:0.0; }

void lag_pole_zero(const LagCompensator *c, double *z, double *p) {
    /* zero: s=-1/T (closer to origin - - - smaller magnitude)
     * pole: s=-1/(beta*T) (farther from origin since beta>1) */
    if (!c||!lag_is_valid(c)) { if(z)*z=0; if(p)*p=0; return; }
    if(z)*z = -1.0/c->T;
    if(p)*p = -1.0/(c->beta*c->T);
}

int lag_is_valid(const LagCompensator *c) {
    return c && c->Kc>0 && c->T>0 && c->beta>=1.0;
}

/* ---- L5: Time Domain ---- */

void lag_step_response(const LagCompensator *c, double tf,
                       int n, double *t, double *y) {
    /* y(t) = Kc*beta + Kc*(1-beta)*exp(-t/(beta*T))
     * y(0)=Kc, y(inf)=Kc*beta */
    if (!c||!t||!y||n<2) return;
    double dt=tf/(n-1), tau=c->beta*c->T;
    double amp=c->Kc*(1.0-c->beta);
    for (int i=0;i<n;i++) {
        double ti=i*dt; t[i]=ti;
        y[i]=(ti<0)?0.0:c->Kc*c->beta+amp*exp(-ti/tau);
    }
}

void lag_impulse_response(const LagCompensator *c, double tf,
                          int n, double *t, double *y) {
    /* y_delta(t) = Kc*(beta-1)/(beta*beta*T)*exp(-t/(beta*T)), t>0 */
    if (!c||!t||!y||n<2) return;
    double dt=tf/(n-1), tau=c->beta*c->T;
    double coeff=c->Kc*(c->beta-1.0)/(c->beta*c->beta*c->T);
    for (int i=0;i<n;i++) {
        double ti=i*dt; t[i]=ti;
        y[i]=(ti<=0)?0.0:coeff*exp(-ti/tau);
    }
}

/* ---- L5: Discretization ---- */

void lag_discretize_tustin(const LagCompensator *c, double Ts,
                           double b[2], double *a1) {
    if (!c||!b||!a1||Ts<=0) { if(b){b[0]=b[1]=0;} if(a1)*a1=0; return; }
    double k1=2.0*c->T/Ts, k2=2.0*c->beta*c->T/Ts;
    double den=1.0+k2;
    if (fabs(den)<1e-15) { b[0]=c->Kc/c->beta; b[1]=0; *a1=0; return; }
    b[0]=c->Kc*(1.0+k1)/den;
    b[1]=c->Kc*(1.0-k1)/den;
    *a1=(1.0-k2)/den;
}

void lag_discretize_zoh(const LagCompensator *c, double Ts,
                        double b[2], double *a1) {
    if (!c||!b||!a1||Ts<=0) { if(b){b[0]=b[1]=0;} if(a1)*a1=0; return; }
    double tau=c->beta*c->T, e=exp(-Ts/tau);
    *a1 = -e;
    double g0=c->Kc, ginf=c->Kc*c->beta;
    double corr=(ginf-g0)*(tau/Ts)*(1.0-e);
    b[0]=g0+corr;
    b[1]=-g0*e+corr*e;
}

double lag_apply_digital(const double b[2], double a1,
                         double u, double *up, double *yp) {
    if (!b||!up||!yp) return 0.0;
    double y=-a1*(*yp)+b[0]*u+b[1]*(*up);
    *up=u; *yp=y; return y;
}

/* ---- L5: Bode Data ---- */

void lag_bode_data(const LagCompensator *c,
                   double wmin, double wmax, int n,
                   double *wo, double *md, double *ph) {
    if (!c||!wo||!md||!ph||n<2) return;
    if (wmin<=0) wmin=0.001;
    if (wmax<=wmin) wmax=wmin*1000.0;
    double lm=log10(wmin), lx=log10(wmax), dl=(lx-lm)/(n-1);
    for (int i=0;i<n;i++) {
        double w=pow(10.0,lm+i*dl); wo[i]=w;
        LagFreqPoint fp=lag_evaluate(c,w);
        md[i]=fp.magnitude_db; ph[i]=fp.phase_deg;
    }
}
