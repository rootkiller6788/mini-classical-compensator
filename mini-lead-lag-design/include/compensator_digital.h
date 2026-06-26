#ifndef COMPENSATOR_DIGITAL_H
#define COMPENSATOR_DIGITAL_H
#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { double b0,b1,b2,a1,a2; double u_prev[2],y_prev[2]; int order; } DigitalFilter;
typedef struct { double Kp,Ki,Kd,Tf,Ts,u_max,u_min,i_accum,prev_error; } DigitalPID;
int digital_lead_tustin(const LeadCompensator *c, double Ts, DigitalFilter *df);
int digital_lag_tustin(const LagCompensator *c, double Ts, DigitalFilter *df);
int digital_lead_lag_tustin(const LeadLagCompensator *c, double Ts, DigitalFilter *ldf, DigitalFilter *gdf);
int digital_matched_pz(const LeadCompensator *c, double Ts, DigitalFilter *df);
int digital_zoh_equivalent(const LeadCompensator *c, double Ts, DigitalFilter *df);
double digital_apply(DigitalFilter *df, double input);
void digital_reset(DigitalFilter *df);
int digital_is_stable(const DigitalFilter *df);
double digital_dc_gain(const DigitalFilter *df);
double digital_min_sampling_rate(const LeadCompensator *c);
double digital_nyquist_freq(const LeadCompensator *c);
int digital_check_aliasing(const LeadCompensator *c, double Ts);
double digital_apply_aw(DigitalFilter *df, double input, double u_min, double u_max);
int digital_to_fixed_point(const DigitalFilter *df, int frac_bits, int *b0f, int *b1f, int *a1f);
double digital_from_fixed_point(int input_fix, int frac_bits);
#ifdef __cplusplus
}
#endif
#endif