/** @file compensator_analytical.c
 * Analytical Compensator Design: Direct Synthesis, IMC, SIMC, Lambda.
 * L4: Internal model principle, pole-zero cancellation,
 *     half-rule reduction, SIMC tuning rules (Skogestad 2003).
 * L5: FOPDT/SOPDT model fitting, lambda selection.
 * Ref: Seborg Ch12, Skogestad JPC 2003, Astrom Ch10. */
#include "compensator_analytical.h"
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L4: Direct Synthesis ---- */

int ds_design_pi(const FOPDTModel *p, double lambda, double *Kp, double *Ki) {
    /* DS PI for FOPDT: G(s)=K*exp(-theta*s)/(tau*s+1)
     * Desired: T(s)=exp(-theta*s)/(lambda*s+1)
     * C(s)=T/(G*(1-T)), PI approx after Pade:
     * Kc=tau/(K*(lambda+theta)), Ti=tau */
    if (!p||!Kp||!Ki||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    *Kp = p->tau/(p->Kdc*(lambda+p->theta));
    *Ki = *Kp/p->tau;
    return 0;
}

int ds_design_pid(const FOPDTModel *p, double lambda, double *Kp, double *Ki, double *Kd) {
    /* DS PID for FOPDT: Kc=tau/(K*(lambda+theta/2)), Ti=tau, Td=theta/2 */
    if (!p||!Kp||!Ki||!Kd||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    *Kp = p->tau/(p->Kdc*(lambda+0.5*p->theta));
    *Ki = *Kp/p->tau;
    *Kd = *Kp*0.5*p->theta;
    return 0;
}

int ds_design_lead_lag(const SOPDTModel *p, double lambda, LeadLagCompensator *c) {
    if (!p||!c||lambda<=0.0) return -1;
    c->Kc = p->tau1*p->tau2/(p->Kdc*lambda*lambda);
    c->T_lead = p->tau1; c->alpha = 0.1;
    c->T_lag = lambda; c->beta = 10.0;
    return 0;
}

/* ---- L4: IMC-Based Design ---- */

int imc_design_pi(const FOPDTModel *p, double lambda, double *Kp, double *Ki) {
    /* IMC PI: G=K*exp(-theta*s)/(tau*s+1), factor as invertible*non-invertible.
     * C=Q/(1-G*Q), Q=inv(G_minus)*F, F=1/(lambda*s+1)
     * PI: Kc=tau/(K*(lambda+theta)), Ti=tau */
    if (!p||!Kp||!Ki||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    double tau_c = lambda + p->theta;
    *Kp = p->tau/(p->Kdc*tau_c);
    *Ki = *Kp/p->tau;
    return 0;
}

int imc_design_pid(const FOPDTModel *p, double lambda, double *Kp, double *Ki, double *Kd) {
    /* IMC PID: Kc=(tau+theta/2)/(K*(lambda+theta/2)),
     * Ti=tau+theta/2, Td=tau*theta/(2*tau+theta) */
    if (!p||!Kp||!Ki||!Kd||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    double tau_c = lambda + 0.5*p->theta;
    *Kp = (p->tau+0.5*p->theta)/(p->Kdc*tau_c);
    double Ti = p->tau + 0.5*p->theta;
    *Ki = *Kp/Ti;
    *Kd = *Kp*(p->tau*p->theta)/(2.0*p->tau+p->theta);
    return 0;
}

int imc_lead_lag(const FOPDTModel *p, double lambda, LeadLagCompensator *c) {
    if (!p||!c||lambda<=0.0) return -1;
    c->Kc = p->tau/(p->Kdc*lambda);
    c->T_lead = p->tau; c->alpha = 0.1;
    c->T_lag = lambda; c->beta = p->tau/(0.1*lambda);
    return 0;
}

/* ---- L4: SIMC Tuning (Skogestad 2003) ---- */

int simc_pi(const FOPDTModel *p, double *Kp, double *Ki) {
    /* SIMC PI rules (Skogestad, J. Process Control, 2003):
     * Kc = tau/(K*(theta+tau_c)) where tau_c = theta (tight)
     * Ti = min(tau, 4*tau_c) */
    if (!p||!Kp||!Ki) return -1;
    if (p->Kdc<=0.0) return -1;
    double tau_c = p->theta;
    if (tau_c<=0.0) tau_c = p->tau*0.1;
    *Kp = p->tau/(p->Kdc*(p->theta+tau_c));
    double Ti = (p->tau < 4.0*tau_c) ? p->tau : 4.0*tau_c;
    *Ki = *Kp/Ti;
    return 0;
}

int simc_pid(const FOPDTModel *p, double *Kp, double *Ki, double *Kd) {
    /* SIMC PID: cascade form. For tau>8*theta, use PID.
     * Kc = tau/(K*(theta+tau_c)), Ti=tau, Td=theta/2 */
    if (!p||!Kp||!Ki||!Kd) return -1;
    if (p->Kdc<=0.0) return -1;
    double tau_c = p->theta;
    if (tau_c<=0.0) tau_c = p->tau*0.1;
    *Kp = p->tau/(p->Kdc*(p->theta+tau_c));
    *Ki = *Kp/p->tau;
    *Kd = *Kp*0.5*p->theta;
    return 0;
}

/* ---- L4: Lambda Tuning ---- */

int lambda_tune_pi(const FOPDTModel *p, double lambda, double *Kp, double *Ki) {
    /* Lambda tuning: Kc = tau/(K*(lambda+theta)), Ti = tau */
    if (!p||!Kp||!Ki||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    *Kp = p->tau/(p->Kdc*(lambda+p->theta));
    *Ki = *Kp/p->tau;
    return 0;
}

int lambda_tune_pid(const FOPDTModel *p, double lambda, double *Kp, double *Ki, double *Kd) {
    if (!p||!Kp||!Ki||!Kd||lambda<=0.0) return -1;
    if (p->Kdc<=0.0||p->tau<=0.0) return -1;
    *Kp = p->tau/(p->Kdc*(lambda+0.5*p->theta));
    *Ki = *Kp/p->tau;
    *Kd = *Kp*0.5*p->theta;
    return 0;
}

/* ---- L5: Model Reduction ---- */

int reduce_to_fopdt(double Kdc, double tau1, double tau2, double theta, FOPDTModel *r) {
    /* Skogestad half-rule: combine two time constants into one.
     * tau = tau1 + tau2/2, effective theta = theta + tau2/2 */
    if (!r) return -1;
    r->Kdc = Kdc;
    r->tau = tau1 + 0.5*tau2;
    r->theta = theta + 0.5*tau2;
    return 0;
}

int reduce_skogestad_half(const FOPDTModel *d, FOPDTModel *r) {
    /* Apply half-rule to convert detailed model to simplified FOPDT */
    if (!d||!r) return -1;
    *r = *d;
    return 0;
}

double compute_lambda_optimal(const FOPDTModel *p) {
    /* Optimal lambda for IMC/SIMC:
     * Tight control: lambda = theta (maximum bandwidth)
     * Robust control: lambda = 2*theta
     * Conservative: lambda = 3*theta */
    if (!p||p->theta<=0.0) return p?p->tau:1.0;
    return p->theta;
}
/* ---- L5: PID Tuning Correlations (Ziegler-Nichols, Cohen-Coon, AMIGO) ---- */
int zn_pi(const FOPDTModel *p, double *Kp, double *Ki) {
    /* Ziegler-Nichols open-loop PI: Kc=0.9*tau/(K*theta), Ti=3.33*theta */
    if (!p||!Kp||!Ki) return -1;
    if (p->theta<=0.0||p->Kdc<=0.0) return -1;
    *Kp = 0.9*p->tau/(p->Kdc*p->theta);
    *Ki = *Kp/(3.33*p->theta);
    return 0;
}
int zn_pid(const FOPDTModel *p, double *Kp, double *Ki, double *Kd) {
    /* ZN open-loop PID: Kc=1.2*tau/(K*theta), Ti=2*theta, Td=0.5*theta */
    if (!p||!Kp||!Ki||!Kd) return -1;
    if (p->theta<=0.0||p->Kdc<=0.0) return -1;
    *Kp = 1.2*p->tau/(p->Kdc*p->theta);
    *Ki = *Kp/(2.0*p->theta);
    *Kd = *Kp*0.5*p->theta;
    return 0;
}
int cohen_coon_pi(const FOPDTModel *p, double *Kp, double *Ki) {
    /* Cohen-Coon PI: Kc=(1/K)*(tau/theta)*(0.9+theta/(12*tau)), Ti=theta*(30+3*theta/tau)/(9+20*theta/tau) */
    if (!p||!Kp||!Ki||p->theta<=0.0||p->Kdc<=0.0) return -1;
    double r=p->theta/p->tau;
    *Kp = (p->tau/(p->Kdc*p->theta))*(0.9+r/12.0);
    double Ti = p->theta*(30.0+3.0*r)/(9.0+20.0*r);
    *Ki = *Kp/Ti;
    return 0;
}
int amigo_pi(const FOPDTModel *p, double *Kp, double *Ki) {
    /* AMIGO (Astrom-Hagglund) PI: Kc=(0.15/K)+(0.35*tau/(K*theta)), Ti=0.35*theta+13.6*theta*tau/(tau+25.2*theta) */
    if (!p||!Kp||!Ki||p->Kdc<=0.0||p->tau<=0.0) return -1;
    *Kp = 0.15/p->Kdc + 0.35*p->tau/(p->Kdc*(p->theta>0.001?p->theta:0.001));
    double Ti = 0.35*(p->theta>0.001?p->theta:0.001) + 13.6*(p->theta>0.001?p->theta:0.001)*p->tau/(p->tau+25.2*(p->theta>0.001?p->theta:0.001));
    *Ki = *Kp/Ti;
    return 0;
}
int chien_hrones_reswick_pi(const FOPDTModel *p, int mode, double *Kp, double *Ki) {
    /* CHR PI tuning: mode=0 (0% overshoot), mode=1 (20% overshoot) */
    if (!p||!Kp||!Ki||p->theta<=0.0||p->Kdc<=0.0) return -1;
    double r=p->theta/p->tau;
    if (mode==0) {
        *Kp = 0.35/(p->Kdc*r);
        *Ki = *Kp/(1.17*p->tau);
    } else {
        *Kp = 0.6/(p->Kdc*r);
        *Ki = *Kp/(p->tau);
    }
    return 0;
}
/* ---- L5: Compensator Selection Guide ---- */
int compensator_select_type(double pm_cur, double pm_des, double ess_cur, double ess_des) {
    /* Returns: 0=no compensation needed, 1=lead only, 2=lag only, 3=lead-lag */
    double pm_deficit = pm_des - pm_cur;
    double ess_ratio = ess_cur/ess_des;
    if (pm_deficit <= 5.0 && ess_ratio <= 1.1) return 0;
    if (pm_deficit > 5.0 && ess_ratio <= 1.1) return 1;
    if (pm_deficit <= 5.0 && ess_ratio > 1.1) return 2;
    return 3;
}
int compensator_compare_performance(const LeadCompensator *lead, const LagCompensator *lag,
                                    const LeadLagCompensator *lead_lag,
                                    double *pm_lead, double *pm_lag, double *pm_ll) {
    if (!pm_lead||!pm_lag||!pm_ll) return -1;
    *pm_lead = lead?lead_max_phase_deg(lead->alpha):0.0;
    *pm_lag = lag?lag_attenuation_db(lag->beta):0.0;
    *pm_ll = lead_lag?lead_max_phase_deg(lead_lag->alpha):0.0;
    return 0;
}
int compensator_order(const LeadLagCompensator *c) {
    if (!c) return 0;
    if (c->alpha>0.99 && c->beta<1.01) return 0;
    if (c->alpha<0.99 && c->beta<1.01) return 1;
    if (c->alpha>0.99 && c->beta>1.01) return 1;
    return 2;
}