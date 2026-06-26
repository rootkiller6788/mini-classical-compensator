#ifndef COMPENSATOR_SENSITIVITY_H
#define COMPENSATOR_SENSITIVITY_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { int n; double *omega,*Smag,*Tmag,*Sph,*Tph; } SensitivityData;
double sens_s_plus_t(double Sm, double Tm);
double sens_peak_s_to_gm(double Ms);
double sens_peak_s_to_pm(double Ms);
double sens_peak_t_to_gm(double Mt);
SensitivityData *sens_data_create(int n);
void sens_data_free(SensitivityData *sd);
void sens_compute_from_loop(double *lm,double *lp,double *w,int n,SensitivityData *sd);
double sens_find_peak_s(const SensitivityData *sd);
double sens_find_peak_t(const SensitivityData *sd);
double sens_bandwidth_from_t(const SensitivityData *sd);
double sens_modulus_margin(double gm_db,double pm_deg);
double sens_delay_margin(double pm_deg,double wc);
int sens_stability_robustness(double gm,double pm,double Ms_tol,double Mt_tol);
double sens_gain_variation_tolerance(double Ms,double pct);
double sens_phase_variation_tolerance(double Ms,double deg);
int sens_is_robustly_stable(double gm,double pm,double gunc,double punc);
#ifdef __cplusplus
}
#endif
#endif