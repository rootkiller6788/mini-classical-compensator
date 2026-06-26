#ifndef CROSS_LIMITING_H
#define CROSS_LIMITING_H

#include "ratio_control.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cross_limiting.h
 * @brief Cross-Limiting Control — prevents dangerous process conditions
 *        by constraining one loop's setpoint based on another's available capacity.
 *
 * Classic application: Combustion control where fuel flow is limited
 * by available air (prevents fuel-rich / smoking), and air flow is
 * limited by fuel demand (prevents excess air / thermal inefficiency).
 *
 * Algorithm (High/Low Select):
 *   SP_fuel_raw = demand
 *   SP_fuel_hi  = air_flow / K_fa          (fuel cannot exceed air/ratio)
 *   SP_fuel     = MIN(SP_fuel_raw, SP_fuel_hi)
 *
 *   SP_air_raw  = demand * K_fa            (stoichiometric air)
 *   SP_air_lo   = fuel_flow * K_fa         (air must at least serve fuel)
 *   SP_air      = MAX(SP_air_raw, SP_air_lo)
 *
 * L1: Cross-limiting signal structures
 * L2: High/Low select architecture
 * L3: Boolean selection logic as algebraic min/max
 * L4: Fuel-rich / excess-air guard theorems
 * L5: Dynamic cross-limiting with rate-of-change limits
 * L6: Boiler combustion control
 * L7: Power plant burner management, refinery fired heaters
 * L8: Advanced cross-limiting with adaptive ratio
 *
 * MIT 6.302: Combustion cross-limiting example
 * Purdue ME 575: Industrial boiler control
 * Cambridge 4F2: Safety-critical ratio control
 */

typedef enum {
    CROSS_LIMIT_INACTIVE = 0,
    CROSS_LIMIT_FUEL_LIMITED_BY_AIR,
    CROSS_LIMIT_AIR_LIMITED_BY_FUEL,
    CROSS_LIMIT_DEMAND_LIMITED,
    CROSS_LIMIT_FULL_CROSS,
    CROSS_LIMIT_COUNT
} cross_limit_status_t;

typedef struct {
    double          fuel_demand_raw;
    double          air_demand_raw;
    double          fuel_flow_actual;
    double          air_flow_actual;
    double          fuel_to_air_ratio;
    double          fuel_to_air_ratio_min;
    double          fuel_to_air_ratio_max;
    double          excess_air_target;
    double          excess_air_actual;
    double          fuel_sp;
    double          air_sp;
    double          fuel_sp_prev;
    double          air_sp_prev;
    double          fuel_hi_limit;
    double          air_lo_limit;
    double          fuel_rate_limit_up;
    double          fuel_rate_limit_down;
    double          air_rate_limit_up;
    double          air_rate_limit_down;
    int             fuel_high_selected;
    int             air_low_selected;
    int             cross_limiting_active;
    unsigned long   cycle_count;
    cross_limit_status_t status;
} cross_limit_t;

typedef struct {
    const char     *name;
    double          lo_selector_inputs[4];
    double          hi_selector_inputs[4];
    size_t          num_lo_inputs;
    size_t          num_hi_inputs;
    double          lo_selected;
    double          hi_selected;
    size_t          lo_winner_index;
    size_t          hi_winner_index;
    double          median_selected;
} signal_selector_t;

int  cross_limit_init(cross_limit_t *cl,
                      double fuel_air_ratio, double excess_air_pct,
                      double ratio_min, double ratio_max);
int  cross_limit_set_demands(cross_limit_t *cl,
                             double fuel_demand, double air_demand);
int  cross_limit_set_actuals(cross_limit_t *cl,
                             double fuel_flow, double air_flow);
int  cross_limit_execute(cross_limit_t *cl, double dt);
int  cross_limit_get_fuel_sp(const cross_limit_t *cl, double *sp);
int  cross_limit_get_air_sp(const cross_limit_t *cl, double *sp);
cross_limit_status_t cross_limit_get_status(const cross_limit_t *cl);
double cross_limit_get_excess_air(const cross_limit_t *cl);
int  cross_limit_set_rate_limits(cross_limit_t *cl,
                                 double fuel_up, double fuel_down,
                                 double air_up, double air_down);
int  cross_limit_reset(cross_limit_t *cl);
int  cross_limit_override(cross_limit_t *cl,
                          double fuel_sp_forced, double air_sp_forced);

int  signal_selector_init(signal_selector_t *ss, const char *name);
int  signal_selector_add_lo(signal_selector_t *ss, double value);
int  signal_selector_add_hi(signal_selector_t *ss, double value);
int  signal_selector_execute(signal_selector_t *ss);
double signal_selector_get_lo(const signal_selector_t *ss);
double signal_selector_get_hi(const signal_selector_t *ss);
double signal_selector_get_median(const signal_selector_t *ss);
size_t signal_selector_get_lo_index(const signal_selector_t *ss);
size_t signal_selector_get_hi_index(const signal_selector_t *ss);
void   signal_selector_clear(signal_selector_t *ss);

const char *cross_limit_status_name(cross_limit_status_t status);

double hi_select(double a, double b);
double lo_select(double a, double b);
double mid_select(double a, double b, double c);
double hi_select_n(const double *values, size_t n, size_t *winning_idx);
double lo_select_n(const double *values, size_t n, size_t *winning_idx);
double median_select_n(const double *values, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* CROSS_LIMITING_H */
