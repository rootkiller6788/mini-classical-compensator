/**
 * @file compensator_freq_design.c
 * Frequency-Domain Compensator Design using Bode Methods.
 * L1-L5: BodeData, stability margins, lead/lag/lag-lead design,
 * loop shaping, sensitivity functions, model fitting.
 * Ref: Ogata Ch7, Nise Ch9, Franklin Ch6, Astrom Ch11.
 */
#include "compensator_freq_design.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L4: Lead Design from PM Spec ---- */

int freq_design_lead_pm(double pm_cur, double pm_tgt,
                        double wc_tgt, double margin,
                        LeadCompensator *comp) {
    if (!comp || pm_tgt <= 0.0 || wc_tgt <= 0.0) return -1;
    double phi = pm_tgt - pm_cur + margin;
    if (phi <= 0.0) { comp->Kc=1.0; comp->T=1.0; comp->alpha=0.99; return 0; }
    if (phi > 65.0) phi = 65.0;
    double sp = sin(phi * M_PI / 180.0);
    double alpha = (1.0 - sp)/(1.0 + sp);
    if (alpha < 0.05) alpha = 0.05;
    comp->T = 1.0/(wc_tgt * sqrt(alpha));
    comp->Kc = 1.0/sqrt(alpha);
    comp->alpha = alpha;
    return 0;
}

int freq_design_lead_from_bode(const PlantFreqData *p,
                               double pm_tgt, double margin,
                               LeadCompensator *comp) {
    if (!p || !comp || p->num_points < 2) return -1;
    double wc = 0.0, pm_cur = 0.0;
    for (int i=1; i<p->num_points; i++) {
        if ((p->mag_db[i-1]>=0 && p->mag_db[i]<=0)||(p->mag_db[i-1]<=0 && p->mag_db[i]>=0)) {
            double t = fabs(p->mag_db[i-1])/(fabs(p->mag_db[i-1])+fabs(p->mag_db[i]));
            double lw = log10(p->omega[i-1])+t*(log10(p->omega[i])-log10(p->omega[i-1]));
            wc = pow(10.0, lw); break;
        }
    }
    if (wc <= 0.0) return -1;
    for (int i=1; i<p->num_points; i++) {
        if (p->omega[i-1]<=wc && p->omega[i]>=wc) {
            double t = (log10(wc)-log10(p->omega[i-1]))/(log10(p->omega[i])-log10(p->omega[i-1]));
            pm_cur = 180.0 + p->phase_deg[i-1] + t*(p->phase_deg[i]-p->phase_deg[i-1]);
            break;
        }
    }
    return freq_design_lead_pm(pm_cur, pm_tgt, wc, margin, comp);
}

/* ---- L4: Lag Design ---- */

int freq_design_lag_gain(double ratio, double wc, LagCompensator *comp) {
    if (!comp || ratio < 1.0 || wc <= 0.0) return -1;
    double beta = ratio > 100.0 ? 100.0 : ratio;
    comp->Kc = 1.0;
    comp->T = 10.0/wc;
    comp->beta = beta;
    return 0;
}

int freq_design_lag_from_bode(const PlantFreqData *p,
                              double ratio, LagCompensator *comp) {
    if (!p || !comp || p->num_points < 2) return -1;
    double wc = 1.0;
    for (int i=1; i<p->num_points; i++) {
        if ((p->mag_db[i-1]>=0 && p->mag_db[i]<=0)||(p->mag_db[i-1]<=0 && p->mag_db[i]>=0)) {
            double t = fabs(p->mag_db[i-1])/(fabs(p->mag_db[i-1])+fabs(p->mag_db[i]));
            double lw = log10(p->omega[i-1])+t*(log10(p->omega[i])-log10(p->omega[i-1]));
            wc = pow(10.0, lw); break;
        }
    }
    return freq_design_lag_gain(ratio, wc, comp);
}

/* ---- L4: Lead-Lag Combined ---- */

int freq_design_lead_lag(double pmc, double pmt, double wc,
                         double ess_ratio, double margin,
                         LeadLagCompensator *comp) {
    if (!comp || wc <= 0.0) return -1;
    LeadCompensator lead;
    freq_design_lead_pm(pmc, pmt, wc, margin, &lead);
    double beta = ess_ratio > 1.0 ? ess_ratio : 10.0;
    if (beta > 100.0) beta = 100.0;
    double wm = 1.0/(lead.T * sqrt(lead.alpha));
    comp->Kc = lead.Kc;
    comp->T_lead = lead.T;
    comp->alpha = lead.alpha;
    comp->T_lag = 10.0/wm;
    comp->beta = beta;
    return 0;
}

/* ---- L5: BodeData Allocation ---- */

BodeData *bode_data_create(int n) {
    if (n < 2) return NULL;
    BodeData *bd = (BodeData*)malloc(sizeof(BodeData));
    if (!bd) return NULL;
    bd->num_points = n;
    bd->omega = (double*)calloc((size_t)n, sizeof(double));
    bd->mag_db = (double*)calloc((size_t)n, sizeof(double));
    bd->phase_deg = (double*)calloc((size_t)n, sizeof(double));
    if (!bd->omega||!bd->mag_db||!bd->phase_deg) {
        free(bd->omega); free(bd->mag_db); free(bd->phase_deg); free(bd);
        return NULL;
    }
    return bd;
}

void bode_data_free(BodeData *bd) {
    if (!bd) return;
    free(bd->omega); free(bd->mag_db); free(bd->phase_deg);
    free(bd);
}

/* ---- L5: Bode for Standard Plants ---- */

void bode_first_order(const FirstOrderPlant *p, double wmin, double wmax, BodeData *bd) {
    if (!p||!bd||bd->num_points<2||wmin<=0.0) return;
    bd->omega_min=wmin; bd->omega_max=wmax;
    double lm=log10(wmin), lx=log10(wmax), dl=(lx-lm)/(bd->num_points-1);
    for (int i=0;i<bd->num_points;i++) {
        double w=pow(10.0,lm+i*dl); bd->omega[i]=w;
        double mag=p->Kdc/sqrt(1.0+w*w*p->tau*p->tau);
        bd->mag_db[i]=20.0*log10(mag>1e-30?mag:1e-30);
        bd->phase_deg[i]=-atan2(w*p->tau,1.0)*180.0/M_PI;
    }
}

void bode_second_order(const SecondOrderPlant *p, double wmin, double wmax, BodeData *bd) {
    if (!p||!bd||bd->num_points<2||wmin<=0.0) return;
    bd->omega_min=wmin; bd->omega_max=wmax;
    double lm=log10(wmin), lx=log10(wmax), dl=(lx-lm)/(bd->num_points-1);
    double wn2=p->wn*p->wn;
    for (int i=0;i<bd->num_points;i++) {
        double w=pow(10.0,lm+i*dl); bd->omega[i]=w;
        double w2=w*w;
        double dr=wn2-w2, di=2.0*p->zeta*p->wn*w;
        double dm=sqrt(dr*dr+di*di);
        double mag=p->Kdc*wn2/(dm>1e-30?dm:1e-30);
        bd->mag_db[i]=20.0*log10(mag>1e-30?mag:1e-30);
        bd->phase_deg[i]=-atan2(di,dr)*180.0/M_PI;
    }
}
/* ---- L5: Stability Margin Computation ---- */

double bode_find_gain_crossover(const BodeData *bd) {
    if (!bd||bd->num_points<2) return 0.0;
    for (int i=1;i<bd->num_points;i++) {
        double m0=bd->mag_db[i-1], m1=bd->mag_db[i];
        if ((m0>=0.0&&m1<=0.0)||(m0<=0.0&&m1>=0.0)) {
            double t=fabs(m0)/(fabs(m0)+fabs(m1));
            double lw=log10(bd->omega[i-1])+t*(log10(bd->omega[i])-log10(bd->omega[i-1]));
            return pow(10.0,lw);
        }
    }
    return 0.0;
}

double bode_find_phase_crossover(const BodeData *bd) {
    if (!bd||bd->num_points<2) return 0.0;
    for (int i=1;i<bd->num_points;i++) {
        double p0=bd->phase_deg[i-1], p1=bd->phase_deg[i];
        if ((p0>=-180.0&&p1<=-180.0)||(p0<=-180.0&&p1>=-180.0)) {
            double t=fabs(p0+180.0)/(fabs(p0+180.0)+fabs(p1+180.0));
            double lw=log10(bd->omega[i-1])+t*(log10(bd->omega[i])-log10(bd->omega[i-1]));
            return pow(10.0,lw);
        }
    }
    return 0.0;
}

double bode_phase_margin(const BodeData *bd) {
    if (!bd||bd->num_points<2) return 0.0;
    double wc=bode_find_gain_crossover(bd);
    if (wc<=0.0) return 90.0;
    double ph=bode_phase_at_freq(bd,wc);
    return 180.0+ph;
}

double bode_gain_margin(const BodeData *bd) {
    if (!bd||bd->num_points<2) return INFINITY;
    double w180=bode_find_phase_crossover(bd);
    if (w180<=0.0) return INFINITY;
    double mag=bode_mag_at_freq(bd,w180);
    return -mag;
}

double bode_mag_at_freq(const BodeData *bd, double w) {
    if (!bd||bd->num_points<2||w<=0.0) return -INFINITY;
    if (w<=bd->omega[0]) return bd->mag_db[0];
    if (w>=bd->omega[bd->num_points-1]) return bd->mag_db[bd->num_points-1];
    for (int i=1;i<bd->num_points;i++) {
        if (bd->omega[i-1]<=w&&bd->omega[i]>=w) {
            double lo=log10(bd->omega[i-1]),hi=log10(bd->omega[i]);
            double lw=log10(w),t=(lw-lo)/(hi-lo);
            return bd->mag_db[i-1]+t*(bd->mag_db[i]-bd->mag_db[i-1]);
        }
    }
    return bd->mag_db[bd->num_points-1];
}

double bode_phase_at_freq(const BodeData *bd, double w) {
    if (!bd||bd->num_points<2||w<=0.0) return 0.0;
    if (w<=bd->omega[0]) return bd->phase_deg[0];
    if (w>=bd->omega[bd->num_points-1]) return bd->phase_deg[bd->num_points-1];
    for (int i=1;i<bd->num_points;i++) {
        if (bd->omega[i-1]<=w&&bd->omega[i]>=w) {
            double lo=log10(bd->omega[i-1]),hi=log10(bd->omega[i]);
            double lw=log10(w),t=(lw-lo)/(hi-lo);
            return bd->phase_deg[i-1]+t*(bd->phase_deg[i]-bd->phase_deg[i-1]);
        }
    }
    return bd->phase_deg[bd->num_points-1];
}

/* ---- L5: Loop Transfer Functions ---- */

void bode_open_loop(const BodeData *p, const BodeData *c, BodeData *r) {
    if (!p||!c||!r||r->num_points<2) return;
    for (int i=0;i<r->num_points;i++) {
        double w=r->omega[i];
        r->mag_db[i]=bode_mag_at_freq(c,w)+bode_mag_at_freq(p,w);
        double ph=bode_phase_at_freq(c,w)+bode_phase_at_freq(p,w);
        while (ph>180.0) ph-=360.0;
        while (ph<-180.0) ph+=360.0;
        r->phase_deg[i]=ph;
    }
}

void bode_closed_loop(const BodeData *loop, BodeData *r) {
    if (!loop||!r||r->num_points<2) return;
    for (int i=0;i<r->num_points;i++) {
        double w=r->omega[i];
        double m_db=bode_mag_at_freq(loop,w);
        double ph_deg=bode_phase_at_freq(loop,w);
        double mag=pow(10.0,m_db/20.0);
        double ph=ph_deg*M_PI/180.0;
        double lr=mag*cos(ph),li=mag*sin(ph);
        double dr=1.0+lr,di=li,d2=dr*dr+di*di;
        double tr=(lr*dr+li*di)/d2,ti=(li*dr-lr*di)/d2;
        double tm=sqrt(tr*tr+ti*ti);
        r->mag_db[i]=20.0*log10(tm>1e-30?tm:1e-30);
        r->phase_deg[i]=atan2(ti,tr)*180.0/M_PI;
    }
}

void bode_sensitivity(const BodeData *loop, BodeData *r) {
    if (!loop||!r||r->num_points<2) return;
    for (int i=0;i<r->num_points;i++) {
        double m_db=bode_mag_at_freq(loop,r->omega[i]);
        double ph_rad=bode_phase_at_freq(loop,r->omega[i])*M_PI/180.0;
        double mag=pow(10.0,m_db/20.0);
        double lr=mag*cos(ph_rad),li=mag*sin(ph_rad);
        double dr=1.0+lr,di=li,d2=dr*dr+di*di;
        double sr=dr/d2,si=-di/d2,sm=sqrt(sr*sr+si*si);
        r->mag_db[i]=20.0*log10(sm>1e-30?sm:1e-30);
        r->phase_deg[i]=atan2(si,sr)*180.0/M_PI;
    }
}

void bode_complementary_sensitivity(const BodeData *loop, BodeData *r) {
    bode_closed_loop(loop,r);
}
/* ---- L5: Plant Model Fitting ---- */

int fit_first_order(const PlantFreqData *d, FirstOrderPlant *p) {
    if (!d||!p||d->num_points<2) return -1;
    p->Kdc = pow(10.0, d->mag_db[0]/20.0);
    double dc = d->mag_db[0], corner = dc-3.0, wc = d->omega[d->num_points-1];
    for (int i=1;i<d->num_points;i++) {
        if (d->mag_db[i-1]>=corner && d->mag_db[i]<=corner) {
            double t=(d->mag_db[i-1]-corner)/(d->mag_db[i-1]-d->mag_db[i]);
            double lw=log10(d->omega[i-1])+t*(log10(d->omega[i])-log10(d->omega[i-1]));
            wc=pow(10.0,lw); break;
        }
    }
    p->tau = 1.0/wc;
    return 0;
}

int fit_second_order(const PlantFreqData *d, SecondOrderPlant *p) {
    if (!d||!p||d->num_points<2) return -1;
    p->Kdc = pow(10.0, d->mag_db[0]/20.0);
    double max_mag=d->mag_db[0], wr=d->omega[0];
    for (int i=1;i<d->num_points;i++) {
        if (d->mag_db[i]>max_mag) { max_mag=d->mag_db[i]; wr=d->omega[i]; }
    }
    double Mr=max_mag-d->mag_db[0];
    if (Mr>0.0) {
        double Mr_lin=pow(10.0,Mr/20.0);
        double z=1.0/(2.0*Mr_lin);
        if (z<0.05) z=0.05;
        if (z>1.0) z=1.0;
        p->zeta=z; p->wn=wr/sqrt(1.0-2.0*z*z);
    } else {
        p->zeta=0.7;
        p->wn=d->omega[d->num_points-1];
        double dc=d->mag_db[0];
        double c3=dc-3.0;
        for (int i=1;i<d->num_points;i++) {
            if (d->mag_db[i-1]>=c3&&d->mag_db[i]<=c3) {
                double t=(d->mag_db[i-1]-c3)/(d->mag_db[i-1]-d->mag_db[i]);
                double lw=log10(d->omega[i-1])+t*(log10(d->omega[i])-log10(d->omega[i-1]));
                p->wn=pow(10.0,lw);
                break;
            }
        }
    }
    return 0;
}

/* ---- L7: Design Verification ---- */
int freq_verify_lead_design(const LeadCompensator *c, const PlantFreqData *p) {
    if (!c||!p||p->num_points<2||!lead_is_valid(c)) return 0;
    for (int i=0;i<p->num_points;i++) {
        LeadFreqPoint cf=lead_evaluate(c,p->omega[i]);
        double loop_mag=p->mag_db[i]+cf.magnitude_db;
        double loop_ph=p->phase_deg[i]+cf.phase_deg;
        if (loop_mag>=-1.0 && loop_mag<=1.0) {
            double pm=180.0+loop_ph;
            if (pm<30.0) return 0;
        }
    }
    return 1;
}
int freq_verify_lag_design(const LagCompensator *c, const PlantFreqData *p) {
    if (!c||!p||p->num_points<2||!lag_is_valid(c)) return 0;
    double dc_gain_orig=pow(10.0,p->mag_db[0]/20.0);
    double dc_gain_new=dc_gain_orig*c->Kc*c->beta;
    return (dc_gain_new>=dc_gain_orig)?1:0;
}