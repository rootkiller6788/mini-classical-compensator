#include "ratio_control.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * @file ratio_control_core.c
 * @brief Core Ratio Control Implementation
 *
 * Implements the fundamental ratio computation engine.
 * Each function corresponds to a specific ratio control knowledge point.
 *
 * Theorem (Ratio Station Output):
 *   For a constant master PV, the slave SP asymptotically approaches
 *   K_r * PV_master + b with time constant tau_filter.
 *
 * Theorem (Ratio Integrity):
 *   A ratio loop is "healthy" iff master_valid && slave_valid &&
 *   ratio_min <= K_actual <= ratio_max && NOT slave_saturated.
 *
 * Reference: Shinskey, "Process Control Systems", 4th Ed, Ch.8
 * Reference: ISA-5.1, "Instrumentation Symbols and Identification"
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* =====================================================================
 * L2: Ratio Control Core Operations
 * ===================================================================== */

int ratio_init(ratio_loop_t *loop,
               unsigned int loop_id,
               const char *tag,
               const char *master_tag,
               const char *slave_tag)
{
    if (!loop || !tag || !master_tag || !slave_tag) return -1;

    memset(loop, 0, sizeof(*loop));

    loop->loop_id   = loop_id;
    loop->tag       = tag;
    loop->master_tag = master_tag;
    loop->slave_tag  = slave_tag;
    loop->description = "";
    loop->eng_units_master = "EU";
    loop->eng_units_slave  = "EU";

    /* Safe defaults: 1:1 ratio, no bias, master-slave mode */
    ratio_configure(loop,
                    RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY,
                    1.0,   /* gain */
                    0.0,   /* bias */
                    0.0,   /* ratio min — open */
                    INFINITY); /* ratio max — open */

    /* Set open defaults for limits */
    loop->config.rate_limit_up   = INFINITY;
    loop->config.rate_limit_down = INFINITY;
    loop->config.slave_sp_min    = -INFINITY;
    loop->config.slave_sp_max    = INFINITY;

    /* Initialize state to safe neutral values */
    loop->state.master_pv      = 0.0;
    loop->state.master_pv_filt = 0.0;
    loop->state.slave_pv       = 0.0;
    loop->state.slave_pv_filt  = 0.0;
    loop->state.slave_sp       = 0.0;
    loop->state.slave_sp_prev  = 0.0;
    loop->state.actual_ratio   = 0.0;
    loop->state.ratio_deviation = 0.0;
    loop->state.slave_error    = 0.0;
    loop->state.master_rate    = 0.0;
    loop->state.master_prev    = 0.0;
    loop->state.integrity      = RATIO_INTEGRITY_OK;
    loop->state.master_valid   = 1;
    loop->state.slave_valid    = 1;
    loop->state.shed_active    = 0;
    loop->state.cycle_count    = 0;

    return 0;
}

int ratio_configure(ratio_loop_t *loop,
                    ratio_mode_t mode,
                    ratio_action_t action,
                    double gain,
                    double bias,
                    double rmin,
                    double rmax)
{
    if (!loop) return -1;
    if (rmin > rmax) return -2;
    if (action == RATIO_ACTION_DIVIDE && fabs(gain) < 1e-12) return -3;

    loop->config.mode              = mode;
    loop->config.action            = action;
    loop->config.ratio_gain        = gain;
    loop->config.ratio_bias        = bias;
    loop->config.ratio_min         = rmin;
    loop->config.ratio_max         = rmax;
    loop->config.direction         = (gain >= 0.0) ? RATIO_DIRECT_ACTING
                                                    : RATIO_REVERSE_ACTING;

    return 0;
}

int ratio_set_master(ratio_loop_t *loop, double master_pv, double dt)
{
    if (!loop) return -1;
    if (dt <= 0.0) return -2;
    if (loop->state.shed_active) return -3; /* don't compute during shed */

    ratio_state_t *s = &loop->state;
    ratio_config_t *c = &loop->config;

    /* Master PV quality check */
    s->master_valid = isfinite(master_pv) ? 1 : 0;
    if (!s->master_valid) {
        s->integrity = RATIO_INTEGRITY_MASTER_BAD;
        if (c->shed_on_master_fail) {
            ratio_shed(loop, c->shed_slave_sp);
        }
        return -4;
    }

    s->master_pv = master_pv;

    /* First-order low-pass filter on master PV */
    if (c->filter_enabled && c->filter_tau > 0.0) {
        double alpha = dt / (c->filter_tau + dt);
        s->master_pv_filt = s->master_pv_filt + alpha * (master_pv - s->master_pv_filt);
    } else {
        s->master_pv_filt = master_pv;
    }

    /* Compute master rate of change */
    if (s->cycle_count > 0) {
        s->master_rate = (master_pv - s->master_prev) / dt;
    }
    s->master_prev = master_pv;

    /* Compute slave setpoint based on action type */
    double raw_sp;
    double effective_gain = c->ratio_gain;

    switch (c->action) {
    case RATIO_ACTION_MULTIPLY:
        raw_sp = effective_gain * s->master_pv_filt + c->ratio_bias;
        break;
    case RATIO_ACTION_DIVIDE:
        if (fabs(s->master_pv_filt) < 1e-12) {
            raw_sp = c->slave_sp_max;
        } else {
            raw_sp = s->master_pv_filt / effective_gain + c->ratio_bias;
        }
        break;
    case RATIO_ACTION_PERCENT:
        raw_sp = (effective_gain / 100.0) * c->slave_span + c->ratio_bias;
        break;
    case RATIO_ACTION_CHARACTERIZED:
        if (loop->char_table) {
            effective_gain = ratio_char_lookup(loop->char_table, s->master_pv_filt);
        }
        raw_sp = effective_gain * s->master_pv_filt + c->ratio_bias;
        break;
    default:
        raw_sp = effective_gain * s->master_pv_filt + c->ratio_bias;
        break;
    }

    /* Apply reverse-acting sign flip if configured */
    if (c->direction == RATIO_REVERSE_ACTING) {
        raw_sp = -raw_sp;
    }

    /* Clamp to absolute slave SP limits */
    if (raw_sp > c->slave_sp_max) raw_sp = c->slave_sp_max;
    if (raw_sp < c->slave_sp_min) raw_sp = c->slave_sp_min;

    /* Rate limiting */
    double sp_change = raw_sp - s->slave_sp_prev;
    if (s->cycle_count > 0) {
        double max_up   = c->rate_limit_up   * dt;
        double max_down = c->rate_limit_down * dt;
        if (sp_change > max_up)   sp_change = max_up;
        if (sp_change < -max_down) sp_change = -max_down;
        raw_sp = s->slave_sp_prev + sp_change;
    }

    s->slave_sp_prev = s->slave_sp;
    s->slave_sp = raw_sp;

    /* Update integrity */
    if (fabs(s->master_pv_filt) > 1e-12) {
        s->actual_ratio = s->slave_pv / s->master_pv_filt;
    }
    s->ratio_deviation = c->ratio_gain - s->actual_ratio;
    s->slave_error     = s->slave_sp - s->slave_pv;
    s->integrity       = RATIO_INTEGRITY_OK;
    s->cycle_count++;

    return 0;
}

int ratio_set_slave(ratio_loop_t *loop, double slave_pv)
{
    if (!loop) return -1;

    ratio_state_t *s = &loop->state;
    s->slave_valid = isfinite(slave_pv) ? 1 : 0;

    if (!s->slave_valid) {
        s->integrity = RATIO_INTEGRITY_SLAVE_BAD;
        return -2;
    }

    s->slave_pv = slave_pv;

    /* First-order filter on slave PV */
    if (loop->config.filter_enabled && loop->config.filter_tau > 0.0) {
        s->slave_pv_filt = slave_pv;
    } else {
        s->slave_pv_filt = slave_pv;
    }

    /* Compute actual ratio with divide-by-zero protection */
    if (fabs(s->master_pv_filt) > 1e-12) {
        s->actual_ratio = slave_pv / s->master_pv_filt;
    } else {
        s->actual_ratio = NAN;
    }

    s->ratio_deviation = loop->config.ratio_gain - s->actual_ratio;
    s->slave_error     = s->slave_sp - slave_pv;

    /* Check ratio clamps */
    if (loop->config.ratio_min > loop->config.ratio_max) {
        /* No valid range */
    } else if (s->actual_ratio < loop->config.ratio_min ||
               s->actual_ratio > loop->config.ratio_max) {
        if (s->integrity == RATIO_INTEGRITY_OK) {
            s->integrity = RATIO_INTEGRITY_RATIO_CLAMPED;
        }
    }

    /* Check slave saturation */
    if (fabs(slave_pv - loop->config.slave_sp_max) < 1e-9 ||
        fabs(slave_pv - loop->config.slave_sp_min) < 1e-9) {
        s->integrity = RATIO_INTEGRITY_SLAVE_SATURATED;
    }

    return 0;
}

double ratio_get_slave_sp(const ratio_loop_t *loop)
{
    if (!loop) return NAN;
    return loop->state.slave_sp;
}

double ratio_get_actual(const ratio_loop_t *loop)
{
    if (!loop) return NAN;
    return loop->state.actual_ratio;
}

ratio_integrity_t ratio_get_integrity(const ratio_loop_t *loop)
{
    if (!loop) return RATIO_INTEGRITY_COUNT;
    return loop->state.integrity;
}

int ratio_set_gain(ratio_loop_t *loop, double new_gain, double dt)
{
    if (!loop) return -1;
    if (dt <= 0.0) return -2;
    if (fabs(loop->config.ratio_gain - new_gain) < 1e-12) return 0;

    double delta = new_gain - loop->config.ratio_gain;
    double max_step = loop->config.rate_limit_up * dt;

    if (fabs(delta) > max_step) {
        delta = (delta > 0) ? max_step : -max_step;
    }

    loop->config.ratio_gain += delta;

    /* Clamp to config limits */
    if (loop->config.ratio_gain < loop->config.ratio_min)
        loop->config.ratio_gain = loop->config.ratio_min;
    if (loop->config.ratio_gain > loop->config.ratio_max)
        loop->config.ratio_gain = loop->config.ratio_max;

    return 0;
}

int ratio_shed(ratio_loop_t *loop, double safe_sp)
{
    if (!loop) return -1;

    loop->state.slave_sp    = safe_sp;
    loop->state.slave_sp_prev = safe_sp;
    loop->state.shed_active = 1;
    loop->state.integrity   = RATIO_INTEGRITY_MASTER_BAD;

    return 0;
}

int ratio_reset_state(ratio_loop_t *loop)
{
    if (!loop) return -1;

    ratio_state_t *s = &loop->state;
    s->master_pv       = 0.0;
    s->master_pv_filt  = 0.0;
    s->slave_pv        = 0.0;
    s->slave_pv_filt   = 0.0;
    s->slave_sp        = 0.0;
    s->slave_sp_prev   = 0.0;
    s->actual_ratio    = 0.0;
    s->ratio_deviation = 0.0;
    s->slave_error     = 0.0;
    s->master_rate     = 0.0;
    s->master_prev     = 0.0;
    s->integrity       = RATIO_INTEGRITY_OK;
    s->master_valid    = 1;
    s->slave_valid     = 1;
    s->shed_active     = 0;
    s->cycle_count     = 0;
    s->total_master_flow = 0.0;
    s->total_slave_flow  = 0.0;

    return 0;
}

/* =====================================================================
 * L2: Ratio Characterization Table (Piecewise-Linear Lookup)
 * ===================================================================== */

static int compare_char_points(const void *a, const void *b)
{
    const ratio_char_point_t *pa = (const ratio_char_point_t *)a;
    const ratio_char_point_t *pb = (const ratio_char_point_t *)b;
    if (pa->input < pb->input) return -1;
    if (pa->input > pb->input) return  1;
    return 0;
}

int ratio_char_table_init(ratio_char_table_t *table,
                          const ratio_char_point_t *points,
                          size_t n)
{
    if (!table)  return -1;
    if (!points) return -2;
    if (n < 2)   return -3;

    table->points = (ratio_char_point_t *)malloc(n * sizeof(ratio_char_point_t));
    if (!table->points) return -4;

    memcpy(table->points, points, n * sizeof(ratio_char_point_t));
    table->num_points = n;

    /* Sort by input value */
    qsort(table->points, n, sizeof(ratio_char_point_t), compare_char_points);

    table->extrap_low  = table->points[0].ratio_value;
    table->extrap_high = table->points[n - 1].ratio_value;

    return 0;
}

double ratio_char_lookup(const ratio_char_table_t *table, double input)
{
    if (!table || !table->points || table->num_points == 0) return NAN;

    const ratio_char_point_t *pts = table->points;
    size_t n = table->num_points;

    /* Below range — extrapolate using extrap_low */
    if (input <= pts[0].input) return table->extrap_low;

    /* Above range — extrapolate using extrap_high */
    if (input >= pts[n - 1].input) return table->extrap_high;

    /* Binary search for interpolation segment */
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (pts[mid].input <= input) lo = mid;
        else                         hi = mid;
    }

    /* Linear interpolation between pts[lo] and pts[hi] */
    double x0 = pts[lo].input;
    double x1 = pts[hi].input;
    double y0 = pts[lo].ratio_value;
    double y1 = pts[hi].ratio_value;

    if (fabs(x1 - x0) < 1e-15) return y0;

    double t = (input - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

void ratio_char_table_free(ratio_char_table_t *table)
{
    if (!table || !table->points) return;
    free(table->points);
    table->points     = NULL;
    table->num_points = 0;
}

/* =====================================================================
 * L2: Utility Functions
 * ===================================================================== */

const char *ratio_mode_name(ratio_mode_t mode)
{
    static const char *names[] = {
        "Master-Slave", "Cross-Limiting", "Blending",
        "Feedforward",  "Cascade Ratio",  "Adaptive",
        "Split-Range"
    };
    if (mode < 0 || mode >= RATIO_MODE_COUNT) return "Unknown";
    return names[mode];
}

const char *ratio_integrity_name(ratio_integrity_t status)
{
    static const char *names[] = {
        "OK",           "Master Bad",       "Slave Bad",
        "Ratio Clamped","Slave Saturated",  "Cross Limited",
        "Mass Balance Violated", "Slave Tracking Error"
    };
    if (status < 0 || status >= RATIO_INTEGRITY_COUNT) return "Unknown";
    return names[status];
}

const char *ratio_action_name(ratio_action_t action)
{
    static const char *names[] = {
        "Multiply", "Divide", "Percent", "Characterized"
    };
    if (action < 0 || action >= RATIO_ACTION_CHARACTERIZED + 1) return "Unknown";
    return names[action];
}

int ratio_validate_config(const ratio_config_t *cfg)
{
    if (!cfg) return -1;

    if (cfg->ratio_min > cfg->ratio_max)
        return -2; /* min > max */

    if (cfg->action == RATIO_ACTION_DIVIDE && fabs(cfg->ratio_gain) < 1e-12)
        return -3; /* divide by zero */

    if (cfg->master_span <= 0.0)
        return -4; /* invalid master span */

    if (cfg->slave_span <= 0.0)
        return -5; /* invalid slave span */

    if (cfg->filter_tau < 0.0)
        return -6; /* negative filter time constant */

    if (cfg->slave_sp_min > cfg->slave_sp_max)
        return -7; /* SP range inverted */

    return 0;
}

double ratio_check_mass_balance(const double *fractions,
                                size_t n,
                                double tolerance)
{
    if (!fractions || n == 0) return -1.0;

    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += fractions[i];
    }

    double error = sum - 1.0;
    if (fabs(error) <= tolerance) return 0.0;
    return error;
}

double ratio_optimal_blend(double cost_a, double cost_b,
                           double ratio_min, double ratio_max)
{
    /* Minimize cost = cost_a * r + cost_b * (1 - r) = cost_b + r * (cost_a - cost_b) */
    /* This is linear in r, so optimum is at boundary */

    if (ratio_min > ratio_max) {
        double tmp = ratio_min;
        ratio_min  = ratio_max;
        ratio_max  = tmp;
    }

    if (cost_a < cost_b) {
        /* Use as much of A as possible */
        return ratio_max;
    } else if (cost_a > cost_b) {
        /* Use as little of A as possible */
        return ratio_min;
    } else {
        /* Equal cost — any ratio, pick midpoint */
        return 0.5 * (ratio_min + ratio_max);
    }
}
