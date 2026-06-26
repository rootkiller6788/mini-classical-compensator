/** @file compensator_digital.c - Digital Compensator Implementation.
 * Tustin/ZOH/matched PZ discretization, fixed-point, anti-windup.
 * Ref: Ogata Ch13, Franklin Ch8, Astrom Ch13 */
#include "compensator_digital.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- L5: Tustin Discretization ---- */

int digital_lead_tustin(const LeadCompensator *c, double Ts, DigitalFilter *df) {
    if (!c||!df||Ts<=0||!lead_is_valid(c)) return -1;
    double k1=2.0*c->T/Ts, k2=2.0*c->alpha*c->T/Ts;
    double den=1.0+k2;
    if (fabs(den)<1e-15) return -1;
    df->b0=c->Kc*(1.0+k1)/den;  df->b1=c->Kc*(1.0-k1)/den;
    df->b2=0.0;  df->a1=(1.0-k2)/den;  df->a2=0.0;
    df->order=1;
    memset(df->u_prev,0,2*sizeof(double));
    memset(df->y_prev,0,2*sizeof(double));
    return 0;
}

int digital_lag_tustin(const LagCompensator *c, double Ts, DigitalFilter *df) {
    if (!c||!df||Ts<=0||!lag_is_valid(c)) return -1;
    double k1=2.0*c->T/Ts, k2=2.0*c->beta*c->T/Ts;
    double den=1.0+k2;
    if (fabs(den)<1e-15) return -1;
    df->b0=c->Kc*(1.0+k1)/den;  df->b1=c->Kc*(1.0-k1)/den;
    df->b2=0.0;  df->a1=(1.0-k2)/den;  df->a2=0.0;
    df->order=1;
    memset(df->u_prev,0,2*sizeof(double));
    memset(df->y_prev,0,2*sizeof(double));
    return 0;
}

int digital_lead_lag_tustin(const LeadLagCompensator *c, double Ts,
                             DigitalFilter *ldf, DigitalFilter *gdf) {
    if (!c||!ldf||!gdf||Ts<=0||!lead_lag_is_valid(c)) return -1;
    LeadCompensator lc = {c->Kc, c->T_lead, c->alpha};
    LagCompensator gc = {1.0, c->T_lag, c->beta};
    if (digital_lead_tustin(&lc,Ts,ldf)!=0) return -1;
    if (digital_lag_tustin(&gc,Ts,gdf)!=0) return -1;
    return 0;
}

int digital_matched_pz(const LeadCompensator *c, double Ts, DigitalFilter *df) {
    /* Match pole-zero: z_z = exp(-Ts/T), z_p = exp(-Ts/(alpha*T))
     * DC gain matching: b0+b1 = (1+a1)*Kc */
    if (!c||!df||Ts<=0||!lead_is_valid(c)) return -1;
    double zz=exp(-Ts/c->T), zp=exp(-Ts/(c->alpha*c->T));
    df->a1=-zp; df->a2=0.0; df->order=1;
    double target_dc=c->Kc*(1.0+df->a1);
    df->b1=-zz*(target_dc-0.0)/(zz-zp+1e-30);
    df->b0=target_dc-df->b1;
    df->b2=0.0;
    memset(df->u_prev,0,2*sizeof(double));
    memset(df->y_prev,0,2*sizeof(double));
    return 0;
}

int digital_zoh_equivalent(const LeadCompensator *c, double Ts, DigitalFilter *df) {
    if (!c||!df||Ts<=0||!lead_is_valid(c)) return -1;
    double tau=c->alpha*c->T, e=exp(-Ts/tau);
    df->a1=-e; df->a2=0.0; df->order=1;
    double g0=c->Kc/c->alpha, corr=c->Kc*(1.0-c->alpha)/c->alpha*(tau/Ts)*(1.0-e);
    df->b0=g0-corr; df->b1=-g0*e+corr*e; df->b2=0.0;
    memset(df->u_prev,0,2*sizeof(double));
    memset(df->y_prev,0,2*sizeof(double));
    return 0;
}

/* ---- L5: Digital Filter Operations ---- */

double digital_apply(DigitalFilter *df, double u) {
    if (!df) return 0.0;
    double y;
    if (df->order==1) {
        y=-df->a1*df->y_prev[0]+df->b0*u+df->b1*df->u_prev[0];
        df->u_prev[0]=u; df->y_prev[0]=y;
    } else {
        y=-df->a1*df->y_prev[0]-df->a2*df->y_prev[1]
          +df->b0*u+df->b1*df->u_prev[0]+df->b2*df->u_prev[1];
        df->u_prev[1]=df->u_prev[0]; df->u_prev[0]=u;
        df->y_prev[1]=df->y_prev[0]; df->y_prev[0]=y;
    }
    return y;
}

void digital_reset(DigitalFilter *df) {
    if (!df) return;
    memset(df->u_prev,0,2*sizeof(double));
    memset(df->y_prev,0,2*sizeof(double));
}

int digital_is_stable(const DigitalFilter *df) {
    /* Check poles inside unit circle: |a1|<1, |a2|<1, |a1|-a2<1 */
    if (!df) return 0;
    if (df->order==1) return fabs(df->a1)<1.0;
    return fabs(df->a2)<1.0 && fabs(df->a1)<1.0+df->a2 && fabs(df->a1)+df->a2<1.0;
}

double digital_dc_gain(const DigitalFilter *df) {
    if (!df) return 0.0;
    double num=df->b0+df->b1+df->b2;
    double den=1.0+df->a1+df->a2;
    return (fabs(den)>1e-15)?num/den:0.0;
}

double digital_min_sampling_rate(const LeadCompensator *c) {
    /* fs >= 20*f_max where f_max is the compensator corner frequency.
     * Corner at 1/(alpha*T) rad/s. fs_min = 20/(2*pi*alpha*T) Hz. */
    if (!c||!lead_is_valid(c)) return 0.0;
    double fmax = 1.0/(2.0*M_PI*c->alpha*c->T);
    return 20.0*fmax;
}

double digital_nyquist_freq(const LeadCompensator *c) {
    if (!c||!lead_is_valid(c)) return 0.0;
    return 1.0/(2.0*M_PI*c->alpha*c->T);
}

int digital_check_aliasing(const LeadCompensator *c, double Ts) {
    if (!c||Ts<=0) return 0;
    double fn=0.5/Ts;
    double fmax=1.0/(2.0*M_PI*c->alpha*c->T);
    return fn>2.0*fmax;
}

double digital_apply_aw(DigitalFilter *df, double u, double umin, double umax) {
    double y=digital_apply(df,u);
    if (y>umax) {
        y=umax;
        df->y_prev[0]=y;
    } else if (y<umin) {
        y=umin;
        df->y_prev[0]=y;
    }
    return y;
}

int digital_to_fixed_point(const DigitalFilter *df, int fb, int *b0f, int *b1f, int *a1f) {
    if (!df||!b0f||!b1f||!a1f||fb<1||fb>30) return -1;
    double scale=(double)(1<<fb);
    *b0f=(int)(df->b0*scale);
    *b1f=(int)(df->b1*scale);
    *a1f=(int)(df->a1*scale);
    return 0;
}

double digital_from_fixed_point(int x, int fb) {
    return (double)x/(double)(1<<fb);
}
/* ---- PID Digital Implementation ---- */
int digital_pid_init(DigitalPID *pid, double Kp, double Ki, double Kd, double Tf, double Ts) {
    if (!pid||Ts<=0.0) return -1;
    pid->Kp=Kp; pid->Ki=Ki; pid->Kd=Kd; pid->Tf=Tf; pid->Ts=Ts;
    pid->u_max=1e6; pid->u_min=-1e6; pid->i_accum=0.0; pid->prev_error=0.0;
    return 0;
}
double digital_pid_apply(DigitalPID *pid, double sp, double mv) {
    if (!pid) return 0.0;
    double e=sp-mv, P=pid->Kp*e;
    pid->i_accum+=pid->Ki*e*pid->Ts;
    if (pid->i_accum>pid->u_max) pid->i_accum=pid->u_max;
    if (pid->i_accum<pid->u_min) pid->i_accum=pid->u_min;
    double D=pid->Kd*(e-pid->prev_error)/pid->Ts;
    D=D*pid->Tf/(pid->Tf+pid->Ts);
    pid->prev_error=e;
    double u=P+pid->i_accum+D;
    if (u>pid->u_max) u=pid->u_max;
    if (u<pid->u_min) u=pid->u_min;
    return u;
}
void digital_pid_reset(DigitalPID *pid) {
    if (!pid) return; pid->i_accum=0.0; pid->prev_error=0.0;
}
int digital_pid_set_limits(DigitalPID *pid, double umin, double umax) {
    if (!pid||umin>=umax) return -1;
    pid->u_min=umin; pid->u_max=umax; return 0;
}
double digital_filter_bw(const DigitalFilter *df, double Ts) {
    if (!df||Ts<=0.0) return 0.0;
    double zp=-df->a1;
    if (fabs(zp)>=1.0) return 0.0;
    return -log(fabs(zp))/Ts;
}
int digital_filter_freq_resp(const DigitalFilter *df, double Ts, double w, double *mag, double *ph) {
    if (!df||!mag||!ph||Ts<=0.0) return -1;
    double wT=w*Ts, c=cos(wT), s=sin(wT);
    double nr=df->b0+df->b1*c+df->b2*cos(2.0*wT);
    double ni=-df->b1*s-df->b2*sin(2.0*wT);
    double dr=1.0+df->a1*c+df->a2*cos(2.0*wT);
    double di=-df->a1*s-df->a2*sin(2.0*wT);
    double d2=dr*dr+di*di;
    double re=(nr*dr+ni*di)/d2, im=(ni*dr-nr*di)/d2;
    double m=sqrt(re*re+im*im);
    *mag=20.0*log10(m>1e-30?m:1e-30);
    *ph=atan2(im,re)*180.0/3.14159265358979323846;
    return 0;
}
int digital_rate_limiter(DigitalFilter *df, double input, double max_rate, double *output) {
    if (!df||!output) return -1;
    double raw=digital_apply(df,input);
    double diff=raw-df->y_prev[0];
    if (diff>max_rate) raw=df->y_prev[0]+max_rate;
    if (diff<-max_rate) raw=df->y_prev[0]-max_rate;
    *output=raw;
    return 0;
}