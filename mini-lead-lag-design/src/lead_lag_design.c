/**
 * @file lead_lag_design.c
 * Lead-Lag Compensator - Combined Implementation
 * C(s) = Kc * (T1*s+1)/(alpha*T1*s+1) * (T2*s+1)/(beta*T2*s+1)
 * alpha<1 for phase lead, beta>1 for DC gain boost.
 * L1-L5 coverage. Ref: Ogata Ch7, Nise Ch9, Dorf Ch10.
 */
#include "lead_lag_design.h"
#include "lead_design.h"
#include "lag_design.h"
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

int lead_lag_design_from_spec(const LeadLagDesignSpec *s, LeadLagCompensator *c) {
    if (!s || !c) return -1;
    double phi_lead = s->phase_margin_desired - s->phase_margin_current + s->phase_margin_extra;
    if (phi_lead > 70.0) phi_lead = 70.0;
    if (phi_lead < 0.0) phi_lead = 0.0;
    double alpha = 0.99;
    if (phi_lead > 0.0) {
        double sp = sin(DEG2RAD(phi_lead));
        alpha = (1.0 - sp)/(1.0 + sp);
        if (alpha < 0.05) alpha = 0.05;
        if (alpha > 0.99) alpha = 0.99;
    }
    double wc = s->gain_crossover;
    if (wc <= 0.0) wc = 1.0;
    double T_lead = 1.0/(wc * sqrt(alpha));
    double beta = s->dc_gain_required > 1.0 ? s->dc_gain_required :
                  (s->ess_current > 0 && s->ess_desired > 0 ?
                   s->ess_current/s->ess_desired : 10.0);
    if (beta < 1.0) beta = 1.0;
    if (beta > 100.0) beta = 100.0;
    double wm = 1.0/(T_lead * sqrt(alpha));
    double T_lag = 10.0 / wm;
    c->Kc = 1.0/sqrt(alpha);
    c->T_lead = T_lead;
    c->alpha = alpha;
    c->T_lag = T_lag;
    c->beta = beta;
    return 0;
}

int lead_lag_design_direct(double pmc, double pmd, double wc,
                           double essc, double essd, double m,
                           LeadLagCompensator *c) {
    if (!c) return -1;
    LeadLagDesignSpec s = {pmd, pmc, wc, essc, essd, 0.0, m};
    return lead_lag_design_from_spec(&s, c);
}

LeadLagFreqPoint lead_lag_evaluate(const LeadLagCompensator *c, double w) {
    LeadLagFreqPoint r; r.omega = w;
    if (!c || !lead_lag_is_valid(c)) {
        r.magnitude=0; r.magnitude_db=-INFINITY; r.phase_deg=0; r.real=0; r.imag=0;
        return r;
    }
    double wT1=w*c->T_lead, waT1=w*c->alpha*c->T_lead;
    double n1r=c->Kc, n1i=c->Kc*wT1;
    double d1r=1.0, d1i=waT1, d1m2=d1r*d1r+d1i*d1i;
    double lr=(n1r*d1r+n1i*d1i)/d1m2, li=(n1i*d1r-n1r*d1i)/d1m2;
    double wT2=w*c->T_lag, wbT2=w*c->beta*c->T_lag;
    double n2r=1.0, n2i=wT2;
    double d2r=1.0, d2i=wbT2, d2m2=d2r*d2r+d2i*d2i;
    double gr=(n2r*d2r+n2i*d2i)/d2m2, gi=(n2i*d2r-n2r*d2i)/d2m2;
    r.real = lr*gr - li*gi;
    r.imag = lr*gi + li*gr;
    double m2=r.real*r.real+r.imag*r.imag;
    r.magnitude=sqrt(m2);
    r.magnitude_db=(m2<1e-30)?-INFINITY:20.0*log10(r.magnitude);
    r.phase_deg=RAD2DEG(atan2(r.imag,r.real));
    return r;
}

double lead_lag_magnitude_db(const LeadLagCompensator *c, double w)
{ return lead_lag_evaluate(c,w).magnitude_db; }
double lead_lag_phase_deg(const LeadLagCompensator *c, double w)
{ return lead_lag_evaluate(c,w).phase_deg; }
double lead_lag_dc_gain(const LeadLagCompensator *c)
{ return c ? c->Kc : 0.0; }

void lead_lag_pole_zero(const LeadLagCompensator *c,
                        double *zl, double *pl, double *zg, double *pg) {
    if (!c||!lead_lag_is_valid(c)) {
        if(zl)*zl=0;
        if(pl)*pl=0;
        if(zg)*zg=0;
        if(pg)*pg=0;
        return;
    }
    if(zl)*zl = -1.0/c->T_lead;
    if(pl)*pl = -1.0/(c->alpha*c->T_lead);
    if(zg)*zg = -1.0/c->T_lag;
    if(pg)*pg = -1.0/(c->beta*c->T_lag);
}

int lead_lag_is_valid(const LeadLagCompensator *c) {
    return c && c->Kc>0 && c->T_lead>0 && c->alpha>0 && c->alpha<1.0
           && c->T_lag>0 && c->beta>=1.0;
}

void lead_lag_step_response(const LeadLagCompensator *c, double tf,
                            int n, double *t, double *y) {
    if (!c||!t||!y||n<2) return;
    double dt=tf/(n-1);
    double tau1=c->alpha*c->T_lead, tau2=c->beta*c->T_lag;
    double A1=-1.0/tau1, B1=c->Kc*(1.0-c->alpha)/(c->alpha*c->alpha*c->T_lead);
    double C1=1.0, D1=c->Kc/c->alpha;
    double A2=-1.0/tau2, B2=(c->beta-1.0)/(c->beta*c->beta*c->T_lag);
    double C2=1.0, D2=1.0/c->beta;
    double x1=0, x2=0;
    for (int i=0;i<n;i++) {
        double ti=i*dt; t[i]=ti;
        double u=(ti<0)?0.0:1.0;
        double yl=C1*x1+D1*u;
        x1 += (A1*x1+B1*u)*dt;
        double yg=C2*x2+D2*yl;
        x2 += (A2*x2+B2*yl)*dt;
        y[i]=yg;
    }
}

void lead_lag_bode_data(const LeadLagCompensator *c,
                        double wmin, double wmax, int n,
                        double *wo, double *md, double *ph) {
    if (!c||!wo||!md||!ph||n<2) return;
    if (wmin<=0) wmin=0.001;
    if (wmax<=wmin) wmax=wmin*1000.0;
    double lm=log10(wmin), lx=log10(wmax), dl=(lx-lm)/(n-1);
    for (int i=0;i<n;i++) {
        double w=pow(10.0,lm+i*dl); wo[i]=w;
        LeadLagFreqPoint fp=lead_lag_evaluate(c,w);
        md[i]=fp.magnitude_db; ph[i]=fp.phase_deg;
    }
}

void lead_lag_discretize_tustin(const LeadLagCompensator *c, double Ts,
                                double bl[2], double *a1l,
                                double bg[2], double *a1g) {
    if (!c||!bl||!a1l||!bg||!a1g||Ts<=0) {
        if(bl){bl[0]=bl[1]=0;} if(a1l)*a1l=0;
        if(bg){bg[0]=bg[1]=0;} if(a1g)*a1g=0; return;
    }
    double k1=2.0*c->T_lead/Ts, k2=2.0*c->alpha*c->T_lead/Ts;
    double den1=1.0+k2;
    bl[0]=c->Kc*(1.0+k1)/den1; bl[1]=c->Kc*(1.0-k1)/den1;
    *a1l=(1.0-k2)/den1;
    double k3=2.0*c->T_lag/Ts, k4=2.0*c->beta*c->T_lag/Ts;
    double den2=1.0+k4;
    bg[0]=(1.0+k3)/den2; bg[1]=(1.0-k3)/den2;
    *a1g=(1.0-k4)/den2;
}

double lead_lag_apply_digital(const double bl[2], double a1l,
                              const double bg[2], double a1g,
                              double u,
                              double *ulp, double *ylp,
                              double *ugp, double *ygp) {
    if (!bl||!bg||!ulp||!ylp||!ugp||!ygp) return 0.0;
    double yl = -a1l*(*ylp) + bl[0]*u + bl[1]*(*ulp);
    *ulp=u; *ylp=yl;
    double y = -a1g*(*ygp) + bg[0]*yl + bg[1]*(*ugp);
    *ugp=yl; *ygp=y;
    return y;
}

double lead_lag_closed_loop_pm(const LeadLagCompensator *c,
                               double pdc, double pp, double wc) {
    (void)pdc; /* DC gain affects magnitude but not phase at crossover */
    LeadLagFreqPoint fp=lead_lag_evaluate(c,wc);
    double pp_deg = -RAD2DEG(atan2(wc,pp));
    return 180.0 + fp.phase_deg + pp_deg;
}

double lead_lag_bandwidth_estimate(const LeadLagCompensator *c) {
    if (!c||!lead_lag_is_valid(c)) return 0.0;
    return 1.0/(c->T_lead*sqrt(c->alpha));
}
