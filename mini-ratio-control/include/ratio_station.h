#ifndef RATIO_STATION_H
#define RATIO_STATION_H

#include "ratio_control.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file ratio_station.h
 * @brief Ratio Station — the computational block that produces
 *        the slave setpoint from the master process variable.
 *
 * In ISA-S5.1 notation, the Ratio Station is the FY (ratio relay)
 * function block. It implements:
 *   SP_slave = K_r * PV_master + bias
 *
 * This header covers:
 *   L1: Ratio station state structures
 *   L2: Ratio station types and limits
 *   L3: Linear fractional transforms for ratio biasing
 *   L4: Scaling laws (span and zero shifts are ratio-preserving)
 *   L5: Dynamic ratio compensation with lead-lag
 *
 * MIT 6.302 §8: Feedforward & Ratio Control
 * Cambridge 3F2: Ratio Control Architectures
 */

typedef enum {
    RATIO_FORMULA_LINEAR = 0,
    RATIO_FORMULA_LINEAR_BIAS,
    RATIO_FORMULA_SQUARE_ROOT,
    RATIO_FORMULA_POLY_2ND,
    RATIO_FORMULA_POLY_3RD,
    RATIO_FORMULA_EXPONENTIAL,
    RATIO_FORMULA_LOGARITHMIC,
    RATIO_FORMULA_COUNT
} ratio_formula_t;

typedef enum {
    RATIO_COMP_DISABLED = 0,
    RATIO_COMP_LEAD_LAG,
    RATIO_COMP_DEADTIME,
    RATIO_COMP_LAG_ONLY,
    RATIO_COMP_COUNT
} ratio_compensation_t;

typedef struct {
    double              numerator_coeffs[4];
    double              denominator_coeffs[4];
    int                 num_order;
    double              state[3];
    double              output_prev;
} ratio_dynamic_compensator_t;

typedef struct {
    ratio_formula_t     formula_type;
    ratio_compensation_t compensation;
    double              poly_coeffs[3];
    double              lead_time;
    double              lag_time;
    double              dead_time;
    double              square_root_cutoff;
    ratio_dynamic_compensator_t dyn_comp;
    double              ratio_bias_internal;
    double              ratio_bias_external;
    int                 bias_tracking;
    double              raw_master_value;
    double              conditioned_master_value;
    double              raw_ratio;
    double              computed_ratio;
    double              compensated_ratio;
    double              final_slave_setpoint;
    double              slave_setpoint_max;
    double              slave_setpoint_min;
    double              setpoint_rate_limit;
    double              prev_setpoint;
    int                 initialized;
} ratio_station_t;

int  ratio_station_init(ratio_station_t *rs,
                        ratio_formula_t formula,
                        ratio_compensation_t comp);
int  ratio_station_set_formula(ratio_station_t *rs,
                               ratio_formula_t formula,
                               const double *coeffs, int n_coeffs);
int  ratio_station_set_dynamics(ratio_station_t *rs,
                                double lead, double lag, double deadtime);
int  ratio_station_set_limits(ratio_station_t *rs,
                              double sp_min, double sp_max, double rate_limit);
int  ratio_station_compute(ratio_station_t *rs,
                           double master_pv, double external_bias, double dt);
double ratio_station_get_setpoint(const ratio_station_t *rs);
double ratio_station_get_computed_ratio(const ratio_station_t *rs);

const char *ratio_formula_name(ratio_formula_t formula);
const char *ratio_compensation_name(ratio_compensation_t comp);

double ratio_apply_formula(ratio_formula_t formula,
                           double master, const double *coeffs);
int ratio_dynamic_compensator_init(ratio_dynamic_compensator_t *comp,
                                   int order,
                                   const double *num, const double *den);
double ratio_dynamic_compensator_step(ratio_dynamic_compensator_t *comp,
                                      double input, double dt);
void ratio_dynamic_compensator_reset(ratio_dynamic_compensator_t *comp);

double ratio_square_root_extract(double value, double cutoff);
double ratio_lead_lag(double input, double lead, double lag,
                      double *state, double dt);

#ifdef __cplusplus
}
#endif

#endif /* RATIO_STATION_H */
