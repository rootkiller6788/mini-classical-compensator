#include "cross_limiting.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @file cross_limiting.c
 * @brief Cross-Limiting Control Implementation
 *
 * Implements the high/low-select logic that prevents dangerous
 * process conditions by constraining one loop's setpoint based
 * on another's available capacity.
 *
 * Classic application — Combustion Control:
 *
 * Fuel-limited-by-air (prevents fuel-rich / smoking):
 *   SP_fuel = MIN(fuel_demand, air_flow / K_fa)
 *
 * Air-limited-by-fuel (prevents excess air / heat loss):
 *   SP_air  = MAX(air_demand, fuel_flow * K_fa)
 *
 * where K_fa = fuel-to-air ratio (stoichiometric ratio).
 *
 * L1: Cross-limit signal path definitions
 * L2: High/low select architecture and logic
 * L3: Boolean min/max algebra
 * L4: Combustion safety theorems (air-rich guaranteed)
 * L5: Dynamic cross-limiting with rate-of-change limits
 * L6: Boiler drum-level cross-limiting
 * L7: Power plant burner management (Fukushima lessons)
 * L8: Adaptive cross-limiting ratio
 *
 * Reference: Shinskey, "Process Control Systems" Ch.8
 * MIT 6.302: Combustion cross-limiting
 * Cambridge 4F2: Safety-critical ratio control
 */

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* -------------------------------------------------------------------
 * L2: Cross-Limiting Core
 * ------------------------------------------------------------------- */

int cross_limit_init(cross_limit_t *cl,
                     double fuel_air_ratio, double excess_air_pct,
                     double ratio_min, double ratio_max)
{
    if (!cl) return -1;
    if (fuel_air_ratio <= 0.0) return -2;
    if (ratio_min > ratio_max) return -3;

    memset(cl, 0, sizeof(*cl));

    cl->fuel_to_air_ratio     = fuel_air_ratio;
    cl->fuel_to_air_ratio_min = ratio_min;
    cl->fuel_to_air_ratio_max = ratio_max;
    cl->excess_air_target     = excess_air_pct;
    cl->excess_air_actual     = 0.0;
    cl->fuel_sp               = 0.0;
    cl->air_sp                = 0.0;
    cl->fuel_sp_prev          = 0.0;
    cl->air_sp_prev           = 0.0;
    cl->fuel_hi_limit         = INFINITY;
    cl->air_lo_limit          = 0.0;
    cl->fuel_high_selected    = 0;
    cl->air_low_selected      = 0;
    cl->cross_limiting_active = 0;
    cl->status                = CROSS_LIMIT_INACTIVE;

    /* Default: no rate limits */
    cl->fuel_rate_limit_up   = INFINITY;
    cl->fuel_rate_limit_down = INFINITY;
    cl->air_rate_limit_up    = INFINITY;
    cl->air_rate_limit_down  = INFINITY;

    return 0;
}

int cross_limit_set_demands(cross_limit_t *cl,
                            double fuel_demand, double air_demand)
{
    if (!cl) return -1;

    cl->fuel_demand_raw = fuel_demand;
    cl->air_demand_raw  = air_demand;
    return 0;
}

int cross_limit_set_actuals(cross_limit_t *cl,
                            double fuel_flow, double air_flow)
{
    if (!cl) return -1;
    if (fuel_flow < 0.0 || air_flow < 0.0) return -2;

    cl->fuel_flow_actual = fuel_flow;
    cl->air_flow_actual  = air_flow;

    /* Compute actual excess air:
     *   excess_air% = (actual_air - stoichiometric_air) / stoichiometric_air * 100
     *   stoichiometric_air = fuel_flow * ratio
     */
    double stoich_air = fuel_flow * cl->fuel_to_air_ratio;
    if (stoich_air > 1e-12) {
        cl->excess_air_actual = (air_flow - stoich_air) / stoich_air * 100.0;
    } else {
        cl->excess_air_actual = INFINITY;
    }

    return 0;
}

int cross_limit_execute(cross_limit_t *cl, double dt)
{
    if (!cl) return -1;
    if (dt <= 0.0) return -2;

    double fuel_demand = cl->fuel_demand_raw;
    double air_demand  = cl->air_demand_raw;
    double actual_fuel = cl->fuel_flow_actual;
    double actual_air  = cl->air_flow_actual;
    double K_fa = cl->fuel_to_air_ratio;

    /* -----------------------------------------------------------------
     * Cross-Limiting Algorithm:
     *
     * Fuel path:  SP_fuel = MIN( fuel_demand,  (actual_air  - bias) / K_fa )
     * Air  path:  SP_air  = MAX( air_demand,   (actual_fuel + bias) * K_fa )
     *
     * This ensures:
     *   - Fuel never exceeds what available air can burn
     *   - Air never falls below what current fuel requires
     *
     * Result: Always air-rich (fuel-lean) during transients = safe.
     * ----------------------------------------------------------------- */

    /* Fuel High Limit: maximum fuel that available air can support */
    double fuel_hi_limit;
    if (K_fa > 1e-12) {
        fuel_hi_limit = actual_air / K_fa;
    } else {
        fuel_hi_limit = INFINITY;
    }

    /* Add excess air margin to fuel limit */
    if (cl->excess_air_target > 0.0) {
        fuel_hi_limit *= (100.0 / (100.0 + cl->excess_air_target));
    }

    cl->fuel_hi_limit = fuel_hi_limit;

    /* Air Low Limit: minimum air needed for current fuel */
    double air_lo_limit = actual_fuel * K_fa;
    if (cl->excess_air_target > 0.0) {
        air_lo_limit *= (100.0 + cl->excess_air_target) / 100.0;
    }
    cl->air_lo_limit = air_lo_limit;

    /* --- Fuel setpoint with cross-limiting --- */
    double fuel_sp_raw = fuel_demand;
    cl->fuel_high_selected = 0;

    if (fuel_demand > fuel_hi_limit) {
        fuel_sp_raw = fuel_hi_limit;
        cl->fuel_high_selected = 1;
        cl->cross_limiting_active = 1;
    }

    /* Rate limit fuel SP */
    if (1) {
        double delta = fuel_sp_raw - cl->fuel_sp_prev;
        double max_up   = cl->fuel_rate_limit_up   * dt;
        double max_down = cl->fuel_rate_limit_down * dt;
        if (delta > max_up)   delta = max_up;
        if (delta < -max_down) delta = -max_down;
        fuel_sp_raw = cl->fuel_sp_prev + delta;
    }

    cl->fuel_sp_prev = cl->fuel_sp;
    cl->fuel_sp      = fuel_sp_raw;

    /* --- Air setpoint with cross-limiting --- */
    double air_sp_raw = air_demand;
    cl->air_low_selected = 0;

    if (air_demand < air_lo_limit) {
        air_sp_raw = air_lo_limit;
        cl->air_low_selected = 1;
        cl->cross_limiting_active = 1;
    }

    /* Rate limit air SP */
    if (1) {
        double delta = air_sp_raw - cl->air_sp_prev;
        double max_up   = cl->air_rate_limit_up   * dt;
        double max_down = cl->air_rate_limit_down * dt;
        if (delta > max_up)   delta = max_up;
        if (delta < -max_down) delta = -max_down;
        air_sp_raw = cl->air_sp_prev + delta;
    }

    cl->air_sp_prev = cl->air_sp;
    cl->air_sp      = air_sp_raw;

    /* Determine status */
    if (cl->fuel_high_selected && cl->air_low_selected) {
        cl->status = CROSS_LIMIT_FULL_CROSS;
    } else if (cl->fuel_high_selected) {
        cl->status = CROSS_LIMIT_FUEL_LIMITED_BY_AIR;
    } else if (cl->air_low_selected) {
        cl->status = CROSS_LIMIT_AIR_LIMITED_BY_FUEL;
    } else if (!cl->cross_limiting_active) {
        cl->status = CROSS_LIMIT_DEMAND_LIMITED;
    }

    cl->cycle_count++;
    return 0;
}

int cross_limit_get_fuel_sp(const cross_limit_t *cl, double *sp)
{
    if (!cl || !sp) return -1;
    *sp = cl->fuel_sp;
    return 0;
}

int cross_limit_get_air_sp(const cross_limit_t *cl, double *sp)
{
    if (!cl || !sp) return -1;
    *sp = cl->air_sp;
    return 0;
}

cross_limit_status_t cross_limit_get_status(const cross_limit_t *cl)
{
    if (!cl) return CROSS_LIMIT_COUNT;
    return cl->status;
}

double cross_limit_get_excess_air(const cross_limit_t *cl)
{
    if (!cl) return NAN;
    return cl->excess_air_actual;
}

int cross_limit_set_rate_limits(cross_limit_t *cl,
                                double fuel_up, double fuel_down,
                                double air_up, double air_down)
{
    if (!cl) return -1;
    if (fuel_up < 0.0 || fuel_down < 0.0 || air_up < 0.0 || air_down < 0.0)
        return -2;

    cl->fuel_rate_limit_up   = fuel_up;
    cl->fuel_rate_limit_down = fuel_down;
    cl->air_rate_limit_up    = air_up;
    cl->air_rate_limit_down  = air_down;
    return 0;
}

int cross_limit_reset(cross_limit_t *cl)
{
    if (!cl) return -1;

    cl->fuel_sp = 0.0;
    cl->air_sp  = 0.0;
    cl->fuel_sp_prev = 0.0;
    cl->air_sp_prev  = 0.0;
    cl->cross_limiting_active = 0;
    cl->fuel_high_selected = 0;
    cl->air_low_selected   = 0;
    cl->status = CROSS_LIMIT_INACTIVE;
    cl->cycle_count = 0;
    return 0;
}

int cross_limit_override(cross_limit_t *cl,
                         double fuel_sp_forced, double air_sp_forced)
{
    if (!cl) return -1;

    cl->fuel_sp = fuel_sp_forced;
    cl->air_sp  = air_sp_forced;
    cl->fuel_sp_prev = fuel_sp_forced;
    cl->air_sp_prev  = air_sp_forced;
    cl->cross_limiting_active = 0;
    cl->status = CROSS_LIMIT_INACTIVE;
    return 0;
}

/* -------------------------------------------------------------------
 * L2: Signal Selector (High/Low/Median)
 * ------------------------------------------------------------------- */

int signal_selector_init(signal_selector_t *ss, const char *name)
{
    if (!ss) return -1;
    memset(ss, 0, sizeof(*ss));
    ss->name = name ? name : "";
    return 0;
}

int signal_selector_add_lo(signal_selector_t *ss, double value)
{
    if (!ss) return -1;
    if (ss->num_lo_inputs >= 4) return -2;
    ss->lo_selector_inputs[ss->num_lo_inputs++] = value;
    return 0;
}

int signal_selector_add_hi(signal_selector_t *ss, double value)
{
    if (!ss) return -1;
    if (ss->num_hi_inputs >= 4) return -2;
    ss->hi_selector_inputs[ss->num_hi_inputs++] = value;
    return 0;
}

int signal_selector_execute(signal_selector_t *ss)
{
    if (!ss) return -1;

    /* Low select */
    ss->lo_selected = INFINITY;
    ss->lo_winner_index = 0;
    for (size_t i = 0; i < ss->num_lo_inputs; i++) {
        if (ss->lo_selector_inputs[i] < ss->lo_selected) {
            ss->lo_selected = ss->lo_selector_inputs[i];
            ss->lo_winner_index = i;
        }
    }

    /* High select */
    ss->hi_selected = -INFINITY;
    ss->hi_winner_index = 0;
    for (size_t i = 0; i < ss->num_hi_inputs; i++) {
        if (ss->hi_selector_inputs[i] > ss->hi_selected) {
            ss->hi_selected = ss->hi_selector_inputs[i];
            ss->hi_winner_index = i;
        }
    }

    /* Median: combine all inputs from both selectors */
    double all_vals[8];
    size_t num_all = 0;
    for (size_t i = 0; i < ss->num_lo_inputs; i++)
        all_vals[num_all++] = ss->lo_selector_inputs[i];
    for (size_t i = 0; i < ss->num_hi_inputs; i++)
        all_vals[num_all++] = ss->hi_selector_inputs[i];

    if (num_all > 0) {
        ss->median_selected = median_select_n(all_vals, num_all);
    }

    return 0;
}

double signal_selector_get_lo(const signal_selector_t *ss)
{
    if (!ss) return NAN;
    return ss->lo_selected;
}

double signal_selector_get_hi(const signal_selector_t *ss)
{
    if (!ss) return NAN;
    return ss->hi_selected;
}

double signal_selector_get_median(const signal_selector_t *ss)
{
    if (!ss) return NAN;
    return ss->median_selected;
}

size_t signal_selector_get_lo_index(const signal_selector_t *ss)
{
    if (!ss) return 0;
    return ss->lo_winner_index;
}

size_t signal_selector_get_hi_index(const signal_selector_t *ss)
{
    if (!ss) return 0;
    return ss->hi_winner_index;
}

void signal_selector_clear(signal_selector_t *ss)
{
    if (!ss) return;
    ss->num_lo_inputs = 0;
    ss->num_hi_inputs = 0;
    ss->lo_selected   = 0.0;
    ss->hi_selected   = 0.0;
    ss->median_selected = 0.0;
}

/* -------------------------------------------------------------------
 * L4: Selector Functions — Algebraic min/max
 *
 * Theorem: MIN(a,b) = (a + b - |a - b|) / 2
 * Theorem: MAX(a,b) = (a + b + |a - b|) / 2
 * ------------------------------------------------------------------- */

double hi_select(double a, double b) { return (a > b) ? a : b; }
double lo_select(double a, double b) { return (a < b) ? a : b; }

double mid_select(double a, double b, double c)
{
    if ((a >= b && a <= c) || (a >= c && a <= b)) return a;
    if ((b >= a && b <= c) || (b >= c && b <= a)) return b;
    return c;
}

double hi_select_n(const double *values, size_t n, size_t *winning_idx)
{
    if (!values || n == 0) return NAN;
    double best = values[0];
    size_t idx  = 0;
    for (size_t i = 1; i < n; i++) {
        if (values[i] > best) { best = values[i]; idx = i; }
    }
    if (winning_idx) *winning_idx = idx;
    return best;
}

double lo_select_n(const double *values, size_t n, size_t *winning_idx)
{
    if (!values || n == 0) return NAN;
    double best = values[0];
    size_t idx  = 0;
    for (size_t i = 1; i < n; i++) {
        if (values[i] < best) { best = values[i]; idx = i; }
    }
    if (winning_idx) *winning_idx = idx;
    return best;
}

double median_select_n(const double *values, size_t n)
{
    if (!values || n == 0) return NAN;

    /* Copy and sort */
    double *sorted = (double *)malloc(n * sizeof(double));
    if (!sorted) return NAN;

    for (size_t i = 0; i < n; i++) sorted[i] = values[i];

    /* Simple insertion sort for small n */
    for (size_t i = 1; i < n; i++) {
        double key = sorted[i];
        size_t j = i;
        while (j > 0 && sorted[j - 1] > key) {
            sorted[j] = sorted[j - 1];
            j--;
        }
        sorted[j] = key;
    }

    double result;
    if (n % 2 == 1) {
        result = sorted[n / 2];
    } else {
        result = 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    }

    free(sorted);
    return result;
}

/* -------------------------------------------------------------------
 * L2: Utility
 * ------------------------------------------------------------------- */

const char *cross_limit_status_name(cross_limit_status_t status)
{
    static const char *names[] = {
        "Inactive", "Fuel Limited by Air", "Air Limited by Fuel",
        "Demand Limited", "Full Cross-Limiting"
    };
    if (status < 0 || status >= CROSS_LIMIT_COUNT) return "Unknown";
    return names[status];
}
