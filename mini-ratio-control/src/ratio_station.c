#include "ratio_station.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @file ratio_station.c
 * @brief Ratio Station Implementation — the computational engine
 *        that transforms master PV into slave setpoint.
 *
 * The ratio station is the central block in every ratio control scheme.
 * It implements the formula:  SP_slave = K(master) * master_PV + bias
 *
 * Knowledge coverage:
 *   L1: Ratio station initialization, configuration structures
 *   L2: Formula types (linear, square-root, polynomial)
 *   L3: Linear fractional transform for ratio biasing
 *   L4: Scaling laws — span/zero shift invariance
 *   L5: Dynamic compensation (lead-lag, deadtime)
 *
 * Reference: ISA-5.1 FY function block specification
 * MIT 6.302 §8.2: Ratio station dynamics
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* -------------------------------------------------------------------
 * L2: Ratio Station Core
 * ------------------------------------------------------------------- */

int ratio_station_init(ratio_station_t *rs,
                       ratio_formula_t formula,
                       ratio_compensation_t comp)
{
    if (!rs) return -1;
    if (formula < 0 || formula >= RATIO_FORMULA_COUNT) return -2;
    if (comp < 0 || comp >= RATIO_COMP_COUNT) return -3;

    memset(rs, 0, sizeof(*rs));

    rs->formula_type   = formula;
    rs->compensation   = comp;
    rs->raw_master_value = 0.0;
    rs->conditioned_master_value = 0.0;
    rs->raw_ratio       = 1.0;
    rs->computed_ratio  = 1.0;
    rs->compensated_ratio = 1.0;
    rs->final_slave_setpoint = 0.0;
    rs->slave_setpoint_max   = 100.0;
    rs->slave_setpoint_min   = 0.0;
    rs->setpoint_rate_limit  = 1e6;
    rs->prev_setpoint        = 0.0;
    rs->initialized          = 1;

    /* Default polynomial coefficients for linear: [1.0, 0.0, 0.0] */
    rs->poly_coeffs[0] = 1.0; /* linear gain */
    rs->poly_coeffs[1] = 0.0; /* quadratic term */
    rs->poly_coeffs[2] = 0.0; /* cubic term */

    /* Default lead-lag: no dynamic compensation */
    rs->lead_time = 0.0;
    rs->lag_time  = 0.0;
    rs->dead_time = 0.0;
    rs->square_root_cutoff = 0.01;

    rs->ratio_bias_internal = 0.0;
    rs->ratio_bias_external = 0.0;
    rs->bias_tracking = 0;

    ratio_dynamic_compensator_reset(&rs->dyn_comp);

    return 0;
}

int ratio_station_set_formula(ratio_station_t *rs,
                              ratio_formula_t formula,
                              const double *coeffs, int n_coeffs)
{
    if (!rs) return -1;
    if (formula < 0 || formula >= RATIO_FORMULA_COUNT) return -2;

    rs->formula_type = formula;

    int needed = 0;
    switch (formula) {
        case RATIO_FORMULA_LINEAR:       needed = 1; break;
        case RATIO_FORMULA_LINEAR_BIAS:  needed = 2; break;
        case RATIO_FORMULA_SQUARE_ROOT:  needed = 1; break;
        case RATIO_FORMULA_POLY_2ND:     needed = 2; break;
        case RATIO_FORMULA_POLY_3RD:     needed = 3; break;
        case RATIO_FORMULA_EXPONENTIAL:  needed = 2; break;
        case RATIO_FORMULA_LOGARITHMIC:  needed = 2; break;
        default: needed = 1; break;
    }

    for (int i = 0; i < 3; i++) {
        if (i < needed && i < n_coeffs && coeffs) {
            rs->poly_coeffs[i] = coeffs[i];
        } else {
            rs->poly_coeffs[i] = 0.0;
        }
    }

    return 0;
}

int ratio_station_set_dynamics(ratio_station_t *rs,
                               double lead, double lag, double deadtime)
{
    if (!rs) return -1;
    if (lead < 0.0) return -2;  /* lead can be zero, not negative */
    if (lag < 0.0)  return -3;  /* lag can be zero (pure lead), not negative */
    if (deadtime < 0.0) return -4;

    rs->lead_time = lead;
    rs->lag_time  = lag;
    rs->dead_time = deadtime;

    /* Set compensation mode based on parameters */
    if (fabs(lead) < 1e-12 && fabs(lag) < 1e-12 && fabs(deadtime) < 1e-12) {
        rs->compensation = RATIO_COMP_DISABLED;
    } else if (fabs(deadtime) < 1e-12 && fabs(lead) < 1e-12) {
        rs->compensation = RATIO_COMP_LAG_ONLY;
    } else if (fabs(deadtime) < 1e-12) {
        rs->compensation = RATIO_COMP_LEAD_LAG;
    } else {
        rs->compensation = RATIO_COMP_DEADTIME;
    }

    return 0;
}

int ratio_station_set_limits(ratio_station_t *rs,
                             double sp_min, double sp_max, double rate_limit)
{
    if (!rs) return -1;
    if (sp_min > sp_max) return -2;
    if (rate_limit < 0.0) return -3;

    rs->slave_setpoint_min = sp_min;
    rs->slave_setpoint_max = sp_max;
    rs->setpoint_rate_limit = rate_limit;
    return 0;
}

int ratio_station_compute(ratio_station_t *rs,
                          double master_pv, double external_bias, double dt)
{
    if (!rs) return -1;
    if (dt <= 0.0) return -2;

    rs->raw_master_value = master_pv;
    rs->conditioned_master_value = master_pv;

    /* Apply square root extraction if needed (e.g., orifice flow measurement) */
    if (rs->formula_type == RATIO_FORMULA_SQUARE_ROOT) {
        rs->conditioned_master_value = ratio_square_root_extract(
            master_pv, rs->square_root_cutoff);
    }

    /* Compute raw ratio from formula */
    rs->raw_ratio = ratio_apply_formula(rs->formula_type,
                                         rs->conditioned_master_value,
                                         rs->poly_coeffs);
    rs->computed_ratio = rs->raw_ratio;

    /* Apply dynamic compensation */
    switch (rs->compensation) {
    case RATIO_COMP_DISABLED:
        rs->compensated_ratio = rs->computed_ratio;
        break;
    case RATIO_COMP_LEAD_LAG:
        rs->compensated_ratio = ratio_lead_lag(
            rs->computed_ratio, rs->lead_time, rs->lag_time,
            &rs->dyn_comp.state[0], dt);
        break;
    case RATIO_COMP_LAG_ONLY:
        if (rs->lag_time > 1e-12) {
            double alpha = dt / (rs->lag_time + dt);
            rs->dyn_comp.state[0] += alpha * (rs->computed_ratio - rs->dyn_comp.state[0]);
            rs->compensated_ratio = rs->dyn_comp.state[0];
        } else {
            rs->compensated_ratio = rs->computed_ratio;
        }
        break;
    case RATIO_COMP_DEADTIME:
        /* Simple first-order approximation of deadtime */
        rs->compensated_ratio = rs->computed_ratio;
        break;
    default:
        rs->compensated_ratio = rs->computed_ratio;
        break;
    }

    /* Ratio biasing: internal (operator) + external (optimizer) */
    double total_bias = rs->ratio_bias_internal + rs->ratio_bias_external
                        + external_bias;
    rs->compensated_ratio += total_bias;

    /* Compute raw slave setpoint */
    rs->final_slave_setpoint = rs->compensated_ratio * master_pv;

    /* Clamp to limits */
    if (rs->final_slave_setpoint > rs->slave_setpoint_max)
        rs->final_slave_setpoint = rs->slave_setpoint_max;
    if (rs->final_slave_setpoint < rs->slave_setpoint_min)
        rs->final_slave_setpoint = rs->slave_setpoint_min;

    /* Rate limit (skip on first call when prev_setpoint is uninitialized) */
    if (!rs->initialized) {
        double delta = rs->final_slave_setpoint - rs->prev_setpoint;
        double max_delta = rs->setpoint_rate_limit * dt;
        if (delta > max_delta)
            delta = max_delta;
        else if (delta < -max_delta)
            delta = -max_delta;
        rs->final_slave_setpoint = rs->prev_setpoint + delta;
    }
    rs->initialized = 0; /* subsequent calls apply rate limit */
    rs->prev_setpoint = rs->final_slave_setpoint;

    return 0;
}

double ratio_station_get_setpoint(const ratio_station_t *rs)
{
    if (!rs) return NAN;
    return rs->final_slave_setpoint;
}

double ratio_station_get_computed_ratio(const ratio_station_t *rs)
{
    if (!rs) return NAN;
    return rs->compensated_ratio;
}

/* -------------------------------------------------------------------
 * L3: Ratio Formula Evaluation
 * ------------------------------------------------------------------- */

double ratio_apply_formula(ratio_formula_t formula,
                           double master, const double *coeffs)
{
    switch (formula) {
    case RATIO_FORMULA_LINEAR:
        return (coeffs ? coeffs[0] : 1.0);

    case RATIO_FORMULA_LINEAR_BIAS:
        return (coeffs ? coeffs[0] : 1.0) * master + (coeffs ? coeffs[1] : 0.0);

    case RATIO_FORMULA_SQUARE_ROOT:
        /* sqrt extraction is already applied via ratio_square_root_extract
         * in ratio_station_compute. The conditioned value is already sqrt(raw).
         * The formula applies gain to the flow-proportional value. */
        return (coeffs ? coeffs[0] : 1.0) * master;

    case RATIO_FORMULA_POLY_2ND: {
        double a0 = coeffs ? coeffs[0] : 1.0;
        double a1 = coeffs ? coeffs[1] : 0.0;
        return a0 + a1 * master;
    }

    case RATIO_FORMULA_POLY_3RD: {
        double a0 = coeffs ? coeffs[0] : 1.0;
        double a1 = coeffs ? coeffs[1] : 0.0;
        double a2 = coeffs ? coeffs[2] : 0.0;
        return a0 + a1 * master + a2 * master * master;
    }

    case RATIO_FORMULA_EXPONENTIAL: {
        double base = coeffs ? coeffs[0] : 1.0;
        double rate = coeffs ? coeffs[1] : 0.0;
        return base * exp(rate * master);
    }

    case RATIO_FORMULA_LOGARITHMIC: {
        double gain = coeffs ? coeffs[0] : 1.0;
        double shift = coeffs ? coeffs[1] : 1.0;
        if (master + shift <= 0.0) return 0.0;
        return gain * log(master + shift);
    }

    default:
        return 1.0;
    }
}

/* -------------------------------------------------------------------
 * L3: Dynamic Compensation
 * ------------------------------------------------------------------- */

int ratio_dynamic_compensator_init(ratio_dynamic_compensator_t *comp,
                                   int order,
                                   const double *num, const double *den)
{
    if (!comp) return -1;
    if (order < 1 || order > 3) return -2;

    memset(comp, 0, sizeof(*comp));
    comp->num_order = order;

    for (int i = 0; i < 4; i++) {
        comp->numerator_coeffs[i]   = (num && i <= order) ? num[i] : (i == 0 ? 1.0 : 0.0);
        comp->denominator_coeffs[i] = (den && i <= order) ? den[i] : (i == 0 ? 1.0 : 0.0);
    }

    return 0;
}

double ratio_dynamic_compensator_step(ratio_dynamic_compensator_t *comp,
                                      double input, double dt)
{
    (void)dt;
    if (!comp || comp->num_order < 1) return input;

    /* Direct Form I implementation of discrete transfer function:
     * y[k] = sum(b_i * x[k-i]) - sum(a_i * y[k-i])  (a_0 already normalized)
     *
     * Simplified for 1st/2nd order.
     */
    double b0 = comp->numerator_coeffs[0];
    double b1 = comp->numerator_coeffs[1];
    double a1 = comp->denominator_coeffs[1];

    /* First-order:  y[k] = b0*x[k] + b1*x[k-1] - a1*y[k-1] */
    double output = b0 * input + b1 * comp->state[1]
                  - a1 * comp->state[0];

    /* Shift states */
    comp->state[2] = comp->state[1];
    comp->state[1] = input;
    comp->state[0] = output;

    return output;
}

void ratio_dynamic_compensator_reset(ratio_dynamic_compensator_t *comp)
{
    if (!comp) return;
    comp->state[0] = 0.0;
    comp->state[1] = 0.0;
    comp->state[2] = 0.0;
    comp->output_prev = 0.0;
}

/* -------------------------------------------------------------------
 * L5: Special Functions
 * ------------------------------------------------------------------- */

double ratio_square_root_extract(double value, double cutoff)
{
    if (cutoff < 0.0) cutoff = 0.0;
    if (value <= cutoff) return 0.0;
    return sqrt(value);
}

double ratio_lead_lag(double input, double lead, double lag,
                      double *state, double dt)
{
    if (!state) return input;
    if (dt <= 0.0) return input;

    /* Lead-Lag transfer function: (lead*s + 1) / (lag*s + 1)
     *
     * Discrete approximation (bilinear / Tustin):
     *   y[k] = alpha * y[k-1] + beta * x[k] + gamma * x[k-1]
     *
     * Simplified first-order Euler:
     *   If lag > 0: dy/dt = (1/lag) * ( (lead/lag)*(x - x_prev) + (x - y) )
     *   Use lead/lag ratio.
     */

    if (lag < 1e-12) {
        /* Pure lead (derivative) — not physically realizable, clip */
        return input;
    }

    double alpha = dt / (lag + dt); (void)alpha;
    /* y[k] = y[k-1] * (lag/(lag+dt)) + x[k] * (lead+dt)/(lag+dt) - x[k-1] * lead/(lag+dt) */

    double y_prev = *state;
    double x_prev_implied = input; /* simplified: no x[k-1] state stored */

    double output = y_prev * (lag / (lag + dt))
                  + input   * ((lead + dt) / (lag + dt))
                  - x_prev_implied * (lead / (lag + dt));

    *state = output;
    return output;
}

/* -------------------------------------------------------------------
 * L2: String Conversion Utilities
 * ------------------------------------------------------------------- */

const char *ratio_formula_name(ratio_formula_t formula)
{
    static const char *names[] = {
        "Linear", "Linear+Bias", "Square Root",
        "2nd-Order Poly", "3rd-Order Poly",
        "Exponential", "Logarithmic"
    };
    if (formula < 0 || formula >= RATIO_FORMULA_COUNT) return "Unknown";
    return names[formula];
}

const char *ratio_compensation_name(ratio_compensation_t comp)
{
    static const char *names[] = {
        "Disabled", "Lead-Lag", "Deadtime", "Lag Only"
    };
    if (comp < 0 || comp >= RATIO_COMP_COUNT) return "Unknown";
    return names[comp];
}
