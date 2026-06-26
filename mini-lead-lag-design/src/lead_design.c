/**
 * @file lead_design.c - Lead Compensator Full Implementation
 * C(s) = Kc*(1+T*s)/(1+alpha*T*s), 0<alpha<1
 * L1-L5 coverage: structs, theorems, algorithms, discretization
 */
#include "lead_design.h"
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

/* ---- L4: Fundamental Formulas ---- */

double lead_max_phase_deg(double alpha) {
    /* phi_max = arcsin((1-alpha)/(1+alpha)) */
    if (alpha <= 0.0 || alpha >= 1.0) return 0.0;
    double s = (1.0 - alpha)/(1.0 + alpha);
    if (s >= 1.0) return 90.0;
    if (s <= -1.0) return -90.0;
    return RAD2DEG(asin(s));
}

double lead_omega_max(double T, double alpha) {
    /* omega_m = 1/(T*sqrt(alpha)) */
    if (T <= 0.0 || alpha <= 0.0 || alpha >= 1.0) return 0.0;
    return 1.0 / (T * sqrt(alpha));
}

double lead_alpha_from_phase(double phase_lead_deg) {
    /* alpha = (1-sin(phi))/(1+sin(phi)) */
    if (phase_lead_deg <= 0.0) return 0.99;
    if (phase_lead_deg >= 90.0) return 0.01;
    double sin_phi = sin(DEG2RAD(phase_lead_deg));
    double alpha = (1.0 - sin_phi)/(1.0 + sin_phi);
    if (alpha < 0.01) alpha = 0.01;
    if (alpha > 0.99) alpha = 0.99;
    return alpha;
}

/* ---- L4: Ogata Design Method ---- */

int lead_design_from_spec(const LeadDesignSpec *spec, LeadCompensator *comp) {
    if (!spec || !comp) return -1;
    if (spec->phase_margin_desired <= 0.0) return -1;
    double phi_lead = spec->phase_margin_desired - spec->phase_margin_current
                      + spec->phase_margin_extra;
    if (phi_lead > 70.0) phi_lead = 70.0;
    if (phi_lead <= 0.0) {
        comp->Kc = 1.0; comp->T = 1.0; comp->alpha = 0.99;
        return 0;
    }
    double alpha = lead_alpha_from_phase(phi_lead);
    if (spec->max_alpha > 0.0 && alpha < spec->max_alpha)
        alpha = spec->max_alpha;
    double wm = spec->gain_crossover_desired > 0.0
                ? spec->gain_crossover_desired : spec->gain_crossover;
    if (wm <= 0.0) return -1;
    comp->T = 1.0/(wm * sqrt(alpha));
    comp->Kc = 1.0/sqrt(alpha);
    comp->alpha = alpha;
    return 0;
}

int lead_design_direct(double pm_cur, double pm_des, double wc,
                       double margin, LeadCompensator *comp) {
    if (!comp) return -1;
    LeadDesignSpec s = {pm_des, pm_cur, wc, wc, margin, 0.05};
    return lead_design_from_spec(&s, comp);
}

/* ---- L5: Frequency Response ---- */

LeadFreqPoint lead_evaluate(const LeadCompensator *comp, double omega) {
    LeadFreqPoint r; r.omega = omega;
    if (!comp || !lead_is_valid(comp)) {
        r.magnitude=0; r.magnitude_db=-INFINITY;
        r.phase_deg=0; r.real=0; r.imag=0; return r;
    }
    double wT = omega*comp->T, waT = omega*comp->alpha*comp->T;
    double nr = comp->Kc, ni = comp->Kc*wT;
    double dr = 1.0, di = waT, d2 = dr*dr + di*di;
    if (d2 < 1e-30) {
        r.magnitude=INFINITY; r.magnitude_db=INFINITY;
        r.phase_deg=0; r.real=INFINITY; r.imag=0; return r;
    }
    r.real = (nr*dr + ni*di)/d2;
    r.imag = (ni*dr - nr*di)/d2;
    double m2 = r.real*r.real + r.imag*r.imag;
    r.magnitude = sqrt(m2);
    r.magnitude_db = (m2<1e-30) ? -INFINITY : 20.0*log10(r.magnitude);
    r.phase_deg = RAD2DEG(atan2(r.imag, r.real));
    return r;
}

double lead_magnitude_db(const LeadCompensator *c, double w)
{ return lead_evaluate(c,w).magnitude_db; }
double lead_phase_deg(const LeadCompensator *c, double w)
{ return lead_evaluate(c,w).phase_deg; }
double lead_dc_gain(const LeadCompensator *c)
{ return c?c->Kc:0.0; }
double lead_hf_gain(const LeadCompensator *c)
{ return (c&&c->alpha>0)?c->Kc/c->alpha:INFINITY; }

void lead_pole_zero(const LeadCompensator *c, double *z, double *p) {
    if (!c||!lead_is_valid(c)) { if(z)*z=0; if(p)*p=0; return; }
    if(z)*z = -1.0/c->T;
    if(p)*p = -1.0/(c->alpha*c->T);
}

int lead_is_valid(const LeadCompensator *c) {
    return c && c->Kc>0 && c->T>0 && c->alpha>0 && c->alpha<1.0;
}

/* ---- L5: Time Domain ---- */

void lead_step_response(const LeadCompensator *c, double tf,
                        int n, double *t, double *y) {
    /* y(t) = Kc + Kc*(1/alpha-1)*exp(-t/(alpha*T)) */
    if (!c||!t||!y||n<2) return;
    double dt = tf/(n-1), tau=c->alpha*c->T;
    double amp = c->Kc*(1.0/c->alpha - 1.0);
    for (int i=0;i<n;i++) {
        double ti = i*dt; t[i]=ti;
        y[i] = (ti<0)?0.0 : c->Kc + amp*exp(-ti/tau);
    }
}

void lead_impulse_response(const LeadCompensator *c, double tf,
                           int n, double *t, double *y) {
    /* y_delta(t) = Kc*(1-alpha)/(alpha^2*T)*exp(-t/(alpha*T)), t>0 */
    if (!c||!t||!y||n<2) return;
    double dt = tf/(n-1), tau=c->alpha*c->T;
    double coeff = c->Kc*(1.0-c->alpha)/(c->alpha*c->alpha*c->T);
    for (int i=0;i<n;i++) {
        double ti = i*dt; t[i]=ti;
        y[i] = (ti<=0)?0.0 : coeff*exp(-ti/tau);
    }
}

/* ---- L5: Discretization ---- */

void lead_discretize_tustin(const LeadCompensator *c, double Ts,
                            double b[2], double *a1) {
    if (!c||!b||!a1||Ts<=0) { if(b){b[0]=b[1]=0;} if(a1)*a1=0; return; }
    double k1=2.0*c->T/Ts, k2=2.0*c->alpha*c->T/Ts, den=1.0+k2;
    if (fabs(den)<1e-15) { b[0]=c->Kc/c->alpha; b[1]=0; *a1=0; return; }
    b[0] = c->Kc*(1.0+k1)/den;
    b[1] = c->Kc*(1.0-k1)/den;
    *a1 = (1.0-k2)/den;
}

void lead_discretize_zoh(const LeadCompensator *c, double Ts,
                         double b[2], double *a1) {
    if (!c||!b||!a1||Ts<=0) { if(b){b[0]=b[1]=0;} if(a1)*a1=0; return; }
    double tau=c->alpha*c->T, e=exp(-Ts/tau);
    *a1 = -e;
    double g0 = c->Kc/c->alpha;
    double corr = c->Kc*(1.0-c->alpha)/c->alpha * (tau/Ts)*(1.0-e);
    b[0] = g0 - corr;
    b[1] = -g0*e + corr*e;
}

double lead_apply_digital(const double b[2], double a1,
                          double u, double *up, double *yp) {
    if (!b||!up||!yp) return 0.0;
    double y = -a1*(*yp) + b[0]*u + b[1]*(*up);
    *up = u; *yp = y;
    return y;
}

/* ---- L5: Bode Data ---- */

void lead_bode_data(const LeadCompensator *c,
                    double wmin, double wmax, int n,
                    double *wo, double *md, double *ph) {
    if (!c||!wo||!md||!ph||n<2) return;
    if (wmin<=0) wmin=0.001;
    if (wmax<=wmin) wmax=wmin*1000.0;
    double lmin=log10(wmin), lmax=log10(wmax);
    double dl=(lmax-lmin)/(n-1);
    for (int i=0;i<n;i++) {
        double w=pow(10.0,lmin+i*dl); wo[i]=w;
        LeadFreqPoint fp=lead_evaluate(c,w);
        md[i]=fp.magnitude_db; ph[i]=fp.phase_deg;
    }
}
