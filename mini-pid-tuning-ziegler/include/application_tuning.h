#ifndef APPLICATION_TUNING_H
#define APPLICATION_TUNING_H
#include "pid_tuning.h"

int app_tune_temperature(double K, double T, double L, pid_params_t *params);
int app_tune_flow(double K, double T, double L, pid_params_t *params);
int app_tune_pressure(double K, double T, double L, pid_params_t *params);
int app_tune_level(double Kv, double L, int surge_tank, pid_params_t *params);
int app_tune_ph(double pH_sp, double valve_size, double max_gain, pid_params_t *params);
int app_tune_dc_motor(double K_motor, double tau_m, double tau_e, pid_params_t *params);
int app_tune_hvac(double K, double T, double L, pid_params_t *params);
int app_tune_chemical_reactor(double K, double T, double L, pid_params_t *params);
int app_tune_paper_headbox(double K, double T, double L, pid_params_t *params);
int app_recommend_adaptation(double K_var, double L_var, double T_var);
pid_tune_method_t app_recommend_method(const fopdt_model_t *model);
int app_auto_tune(int app_type, double K, double T, double L, pid_params_t *params);
int app_estimate_fopdt(double y_initial, double y_final, double max_slope,
                       double t_at_max_slope, double step_mag, fopdt_model_t *model);

#endif
