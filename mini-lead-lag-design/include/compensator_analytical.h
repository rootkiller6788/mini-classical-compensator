#ifndef COMPENSATOR_ANALYTICAL_H
#define COMPENSATOR_ANALYTICAL_H
#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { double Kdc,tau,theta; } FOPDTModel;
typedef struct { double Kdc,tau1,tau2,theta; } SOPDTModel;
int ds_design_pi(const FOPDTModel *p,double lambda,double *Kp,double *Ki);
int ds_design_pid(const FOPDTModel *p,double lambda,double *Kp,double *Ki,double *Kd);
int ds_design_lead_lag(const SOPDTModel *p,double lambda,LeadLagCompensator *c);
int imc_design_pi(const FOPDTModel *p,double lambda,double *Kp,double *Ki);
int imc_design_pid(const FOPDTModel *p,double lambda,double *Kp,double *Ki,double *Kd);
int imc_lead_lag(const FOPDTModel *p,double lambda,LeadLagCompensator *c);
int simc_pi(const FOPDTModel *p,double *Kp,double *Ki);
int simc_pid(const FOPDTModel *p,double *Kp,double *Ki,double *Kd);
int lambda_tune_pi(const FOPDTModel *p,double lambda,double *Kp,double *Ki);
int lambda_tune_pid(const FOPDTModel *p,double lambda,double *Kp,double *Ki,double *Kd);
int reduce_to_fopdt(double Kdc,double tau1,double tau2,double theta,FOPDTModel *r);
int reduce_skogestad_half(const FOPDTModel *d,FOPDTModel *r);
double compute_lambda_optimal(const FOPDTModel *p);
#ifdef __cplusplus
}
#endif
#endif