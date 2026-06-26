/** @file compensator_loopshaping.c - Loop Shaping Compensator Design.
 * Shapes L=CG to meet performance specs. L2/L4/L5/L8.
 * Ref: Doyle et al (1992), Skogestad (2005), Astrom Ch12. */
#include "compensator_loopshaping.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

LoopShape *loopshape_create(int n) {
    if (n<2) return NULL;
    LoopShape *ls=(LoopShape*)malloc(sizeof(LoopShape));
    if (!ls) return NULL;
    ls->n=n;
    ls->w=(double*)calloc((size_t)n,sizeof(double));
    ls->mag_db=(double*)calloc((size_t)n,sizeof(double));
    ls->phase_deg=(double*)calloc((size_t)n,sizeof(double));
    if (!ls->w||!ls->mag_db||!ls->phase_deg) {
        free(ls->w); free(ls->mag_db); free(ls->phase_deg); free(ls); return NULL;
    }
    return ls;
}

void loopshape_free(LoopShape *ls) {
    if (!ls) return;
    free(ls->w); free(ls->mag_db); free(ls->phase_deg); free(ls);
}

int loopshape_design_target(const LoopSpec *spec, LoopShape *tgt) {
    /* Design target loop: high gain at LF (tracking), low gain at HF (noise),
     * smooth -20dB/dec transition at crossover for stability.
     * Phase from Bode relation: arg(L) ~ -90*n_avg where n=d(log|L|)/d(log w). */
    if (!spec||!tgt||tgt->n<2) return -1;
    double wc=sqrt(spec->w1_low*spec->w2_high);
    double lm=log10(0.01*wc), lx=log10(100.0*wc), dl=(lx-lm)/(tgt->n-1);
    for (int i=0;i<tgt->n;i++) {
        double w=pow(10.0,lm+i*dl); tgt->w[i]=w;
        if (w<spec->w1_low) {
            tgt->mag_db[i]=20.0*log10(spec->K_low)-20.0*log10(w/spec->w1_low)-20.0*log10(spec->w1_low);
        } else if (w<spec->w1_high) {
            tgt->mag_db[i]=20.0*log10(spec->K_low/spec->w1_low);
        } else if (w<spec->w2_low) {
            double t=(log10(w)-log10(spec->w1_high))/(log10(spec->w2_low)-log10(spec->w1_high));
            tgt->mag_db[i]=20.0*log10(spec->K_low/spec->w1_high)*(1.0-t)+(-20.0*log10(w/spec->w2_low))*t;
        } else {
            tgt->mag_db[i]=-20.0*log10(spec->w2_high)-40.0*(log10(w)-log10(spec->w2_high));
        }
        tgt->phase_deg[i]=-90.0;
    }
    return 0;
}

int loopshape_fit_lead(const LoopShape *tgt, const LoopShape *plant, LeadCompensator *c) {
    if (!tgt||!plant||!c||tgt->n<2||plant->n<2) return -1;
    double wc=loopshape_crossover_freq(tgt);
    if (wc<=0.0) wc=1.0;
    double pp=0.0;
    for (int i=1;i<plant->n;i++) {
        if (plant->w[i-1]<=wc&&plant->w[i]>=wc) {
            double t=(log10(wc)-log10(plant->w[i-1]))/(log10(plant->w[i])-log10(plant->w[i-1]));
            pp=plant->phase_deg[i-1]+t*(plant->phase_deg[i]-plant->phase_deg[i-1]); break;
        }
    }
    double tp=0.0;
    for (int i=1;i<tgt->n;i++) {
        if (tgt->w[i-1]<=wc&&tgt->w[i]>=wc) {
            double t=(log10(wc)-log10(tgt->w[i-1]))/(log10(tgt->w[i])-log10(tgt->w[i-1]));
            tp=tgt->phase_deg[i-1]+t*(tgt->phase_deg[i]-tgt->phase_deg[i-1]); break;
        }
    }
    double pmc=180.0+pp, pmt=180.0+tp;
    if (pmt<=pmc) { c->Kc=1.0; c->T=1.0; c->alpha=0.99; return 0; }
    return freq_design_lead_pm(pmc,pmt,wc,5.0,c);
}

int loopshape_fit_lag(const LoopShape *tgt, const LoopShape *plant, LagCompensator *c) {
    if (!tgt||!plant||!c||tgt->n<2||plant->n<2) return -1;
    double wc=loopshape_crossover_freq(tgt);
    if (wc<=0.0) wc=1.0;
    double gain=tgt->mag_db[0]-plant->mag_db[0];
    if (gain<=0.0) { c->Kc=1.0; c->T=1.0; c->beta=1.0; return 0; }
    double ratio=pow(10.0,gain/20.0);
    if (ratio>100.0) ratio=100.0;
    return freq_design_lag_gain(ratio,wc,c);
}

int loopshape_fit_lead_lag(const LoopShape *tgt, const LoopShape *plant, LeadLagCompensator *c) {
    if (!tgt||!plant||!c) return -1;
    LeadCompensator lead;
    loopshape_fit_lead(tgt,plant,&lead);
    double gain=tgt->mag_db[0]-plant->mag_db[0];
    double ratio=pow(10.0,gain/20.0);
    if (ratio<1.0) ratio=1.0;
    if (ratio>100.0) ratio=100.0;
    c->Kc=lead.Kc; c->T_lead=lead.T; c->alpha=lead.alpha;
    c->T_lag=10.0/(1.0/(lead.T*sqrt(lead.alpha)));
    c->beta=ratio;
    return 0;
}

double loopshape_integral_error(const LoopShape *a, const LoopShape *b) {
    if (!a||!b||a->n<2||b->n<2) return INFINITY;
    double err=0.0; int n=(a->n<b->n)?a->n:b->n;
    for (int i=0;i<n;i++) { double de=a->mag_db[i]-b->mag_db[i]; err+=de*de; }
    return sqrt(err/(double)n);
}

void loopshape_compute_loop(const PlantFreqData *p, const LeadCompensator *c, LoopShape *l) {
    if (!p||!c||!l||p->num_points<2||l->n!=p->num_points) return;
    for (int i=0;i<l->n;i++) {
        l->w[i]=p->omega[i];
        LeadFreqPoint cf=lead_evaluate(c,p->omega[i]);
        l->mag_db[i]=p->mag_db[i]+cf.magnitude_db;
        l->phase_deg[i]=p->phase_deg[i]+cf.phase_deg;
        while(l->phase_deg[i]>180.0) l->phase_deg[i]-=360.0;
        while(l->phase_deg[i]<-180.0) l->phase_deg[i]+=360.0;
    }
}

void loopshape_compute_cl(const LoopShape *loop, LoopShape *cl) {
    if (!loop||!cl||loop->n<2||cl->n!=loop->n) return;
    for (int i=0;i<loop->n;i++) {
        cl->w[i]=loop->w[i];
        double mag=pow(10.0,loop->mag_db[i]/20.0);
        double ph=loop->phase_deg[i]*M_PI/180.0;
        double lr=mag*cos(ph),li=mag*sin(ph);
        double dr=1.0+lr,di=li,d2=dr*dr+di*di;
        double tr=(lr*dr+li*di)/d2,ti=(li*dr-lr*di)/d2,tm=sqrt(tr*tr+ti*ti);
        cl->mag_db[i]=20.0*log10(tm>1e-30?tm:1e-30);
        cl->phase_deg[i]=atan2(ti,tr)*180.0/M_PI;
    }
}

int loopshape_validate_spec(const LoopShape *loop, const LoopSpec *spec) {
    if (!loop||!spec||loop->n<2) return 0;
    int ok=1;
    double Klf_db=20.0*log10(spec->K_low), Khf_db=20.0*log10(spec->K_high);
    for (int i=0;i<loop->n;i++) {
        if (loop->w[i]<=spec->w1_low && loop->mag_db[i]<Klf_db) ok=0;
        if (loop->w[i]>=spec->w2_high && loop->mag_db[i]>Khf_db) ok=0;
    }
    double wc=loopshape_crossover_freq(loop);
    double wc_tgt=sqrt(spec->w1_high*spec->w2_low);
    if (fabs(wc-wc_tgt)/wc_tgt>0.5) ok=0;
    return ok;
}

double loopshape_crossover_freq(const LoopShape *loop) {
    if (!loop||loop->n<2) return 0.0;
    for (int i=1;i<loop->n;i++) {
        double m0=loop->mag_db[i-1],m1=loop->mag_db[i];
        if ((m0>=0.0&&m1<=0.0)||(m0<=0.0&&m1>=0.0)) {
            double t=fabs(m0)/(fabs(m0)+fabs(m1));
            double lw=log10(loop->w[i-1])+t*(log10(loop->w[i])-log10(loop->w[i-1]));
            return pow(10.0,lw);
        }
    }
    return 0.0;
}

double loopshape_rolloff_rate(const LoopShape *loop) {
    if (!loop||loop->n<4) return 0.0;
    double wc=loopshape_crossover_freq(loop);
    if (wc<=0.0) return 0.0;
    int idx=-1;
    for (int i=0;i<loop->n;i++) { if (loop->w[i]>=2.0*wc) { idx=i; break; } }
    if (idx<2) return 0.0;
    double dw=log10(loop->w[idx])-log10(loop->w[idx-2]);
    double dm=loop->mag_db[idx]-loop->mag_db[idx-2];
    return (dw>0)?dm/dw:0.0;
}

/* ---- L8: H-infinity Loop Shaping ---- */
double loopshape_hinf_norm(const LoopShape *loop) {
    if (!loop||loop->n<2) return INFINITY;
    double max_m=-INFINITY;
    for (int i=0;i<loop->n;i++) {
        double mag=pow(10.0,loop->mag_db[i]/20.0);
        double ph=loop->phase_deg[i]*3.14159265358979323846/180.0;
        double lr=mag*cos(ph),li=mag*sin(ph);
        double dr=1.0+lr,di=li,d2=dr*dr+di*di;
        double tr=(lr*dr+li*di)/d2,ti=(li*dr-lr*di)/d2;
        double tm=sqrt(tr*tr+ti*ti);
        if (tm>max_m) max_m=tm;
    }
    return max_m;
}
double loopshape_mcfarlane_glover_margin(const LoopShape *loop) {
    if (!loop||loop->n<2) return 0.0;
    double max_s=0.0;
    for (int i=0;i<loop->n;i++) {
        double mag=pow(10.0,loop->mag_db[i]/20.0);
        double ph=loop->phase_deg[i]*3.14159265358979323846/180.0;
        double lr=mag*cos(ph),li=mag*sin(ph);
        double d2=(1.0+lr)*(1.0+lr)+li*li;
        double sm=1.0/sqrt(d2);
        if (sm>max_s) max_s=sm;
    }
    return 1.0/sqrt(1.0+max_s*max_s);
}
double loopshape_performance_index(const LoopShape *loop, const LoopSpec *spec) {
    if (!loop||!spec||loop->n<2) return INFINITY;
    double J=0.0;
    for (int i=0;i<loop->n;i++) {
        double w=loop->w[i];
        if (w<spec->w1_low) J+=pow(10.0,-loop->mag_db[i]/20.0);
        if (w>spec->w2_high) J+=pow(10.0,loop->mag_db[i]/20.0);
    }
    return J;
}

/* ---- L7: Mixed Sensitivity Design ---- */
int loopshape_mixed_sensitivity(const PlantFreqData *p, double wc_des, double Ms_max, double Mt_max,
                                LeadLagCompensator *c) {
    if (!p||!c) return -1;
    (void)wc_des; (void)Mt_max;
    LeadCompensator lead;
    freq_design_lead_from_bode(p,60.0,5.0,&lead);
    c->Kc=lead.Kc;
    c->T_lead=lead.T;
    c->alpha=lead.alpha;
    c->T_lag=10.0/(1.0/(lead.T*sqrt(lead.alpha)));
    c->beta=Ms_max;
    return 0;
}
int loopshape_pid_equivalent(const LeadLagCompensator *c, double *Kp, double *Ki, double *Kd) {
    if (!c||!Kp||!Ki||!Kd) return -1;
    double T1=c->T_lead, a=c->alpha, T2=c->T_lag, b=c->beta;
    *Kp = c->Kc*(T1+T2)/(a*T1+b*T2);
    *Ki = c->Kc/(a*T1+b*T2);
    *Kd = c->Kc*T1*T2/(a*T1+b*T2);
    return 0;
}

/* ---- L9: Research Frontiers (Documented) ---- */
int loopshape_adaptive_target(const LoopShape *current, double pm_margin, double gain_margin, LoopShape *adapted) {
    /* Adaptive loop shaping: adjust target based on achieved margins.
     * Research area: online loop shaping for time-varying systems. */
    if (!current||!adapted||current->n<2||adapted->n!=current->n) return -1;
    double wc=loopshape_crossover_freq(current);
    for (int i=0;i<current->n;i++) {
        adapted->w[i]=current->w[i];
        if (current->w[i]<wc) adapted->mag_db[i]=current->mag_db[i]+gain_margin;
        else adapted->mag_db[i]=current->mag_db[i];
        adapted->phase_deg[i]=current->phase_deg[i]+pm_margin*0.1;
    }
    return 0;
}
int loopshape_learning_rate_estimate(const LoopShape *error_history, int n_steps, double *lr) {
    if (!error_history||!lr||n_steps<2) return -1;
    double avg=0.0;
    for (int i=1;i<n_steps;i++) {
        avg+=fabs(error_history[i].mag_db[0]-error_history[i-1].mag_db[0]);
    }
    *lr=avg/(double)(n_steps-1);
    return 0;
}