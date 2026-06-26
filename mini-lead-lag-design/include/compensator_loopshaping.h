#ifndef COMPENSATOR_LOOPSHAPING_H
#define COMPENSATOR_LOOPSHAPING_H
#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include "compensator_freq_design.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { double *w; double *mag_db; double *phase_deg; int n; } LoopShape;
typedef struct { double w1_low; double w1_high; double w2_low; double w2_high; double K_low; double K_high; } LoopSpec;
LoopShape *loopshape_create(int n);
void loopshape_free(LoopShape *ls);
int loopshape_design_target(const LoopSpec *spec, LoopShape *target);
int loopshape_fit_lead(const LoopShape *target, const LoopShape *plant, LeadCompensator *c);
int loopshape_fit_lag(const LoopShape *target, const LoopShape *plant, LagCompensator *c);
int loopshape_fit_lead_lag(const LoopShape *target, const LoopShape *plant, LeadLagCompensator *c);
double loopshape_integral_error(const LoopShape *a, const LoopShape *b);
void loopshape_compute_loop(const PlantFreqData *p, const LeadCompensator *c, LoopShape *l);
void loopshape_compute_cl(const LoopShape *loop, LoopShape *cl);
int loopshape_validate_spec(const LoopShape *loop, const LoopSpec *spec);
double loopshape_crossover_freq(const LoopShape *loop);
double loopshape_rolloff_rate(const LoopShape *loop);
#ifdef __cplusplus
}
#endif
#endif