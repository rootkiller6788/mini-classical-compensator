#include "blending_control.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * @file blending_control.c
 * @brief Multi-Component Blending Control Implementation
 *
 * Implements mass/volume-balance-preserving blending of N streams.
 * Each component is maintained at a prescribed fraction of total flow.
 *
 * Core mass balance:  sum(f_i * Q_total) = Q_total
 *                     sum(f_i) = 1.0  (fraction constraint)
 * where f_i = Q_i / Q_total, Q_total = sum(Q_i)
 *
 * L4: Mass conservation — Lagrange multipliers for blend optimization
 * L5: Linear programming for least-cost blend
 * L6: Gasoline blending (octane, RVP, sulfur constraints)
 * L7: Refinery blending operations (Boeing fuel spec, ASTM standards)
 *
 * Reference: Shinskey, "Process Control Systems" Ch.12
 * Purdue ME 575: Industrial Blending
 * MIT 6.302: Feedforward for multi-component
 */

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* -------------------------------------------------------------------
 * L2: Blend System Lifecycle
 * ------------------------------------------------------------------- */

int blend_system_init(blend_system_t *bs, blend_basis_t basis,
                      blend_strategy_t strategy, size_t max_components)
{
    if (!bs) return -1;
    if (max_components < 2) return -2;

    memset(bs, 0, sizeof(*bs));

    bs->basis        = basis;
    bs->strategy     = strategy;
    bs->max_components = max_components;
    bs->num_components = 0;

    bs->components = (blend_component_t *)calloc(max_components,
                                                  sizeof(blend_component_t));
    if (!bs->components) return -3;

    bs->total_flow_sp          = 100.0;
    bs->total_flow_actual      = 0.0;
    bs->fraction_tolerance     = 1e-4;
    bs->mass_balance_tolerance = 1e-4;
    bs->min_total_flow         = 0.0;
    bs->max_total_flow         = 1000.0;
    bs->blend_quality_sp       = 0.0;
    bs->blend_cost_per_unit    = 0.0;

    return 0;
}

int blend_add_component(blend_system_t *bs, const char *name,
                        double fraction_sp, double density,
                        double cost, double quality,
                        double min_frac, double max_frac,
                        double min_flow, double max_flow)
{
    if (!bs || !bs->components) return -1;
    if (bs->num_components >= bs->max_components) return -2;
    if (!name) return -3;
    if (fraction_sp < 0.0 || fraction_sp > 1.0) return -4;

    blend_component_t *bc = &bs->components[bs->num_components];
    memset(bc, 0, sizeof(*bc));

    bc->name              = name;
    bc->fraction_setpoint = fraction_sp;
    bc->fraction_actual   = 0.0;
    bc->flow_rate         = 0.0;
    bc->density           = density;
    bc->cost_per_unit     = cost;
    bc->quality_value     = quality;
    bc->min_fraction      = min_frac;
    bc->max_fraction      = max_frac;
    bc->min_flow          = min_flow;
    bc->max_flow          = max_flow;
    bc->active            = 1;

    bs->num_components++;
    return 0;
}

int blend_remove_component(blend_system_t *bs, size_t index)
{
    if (!bs) return -1;
    if (index >= bs->num_components) return -2;

    /* Shift remaining components down */
    for (size_t i = index; i < bs->num_components - 1; i++) {
        bs->components[i] = bs->components[i + 1];
    }
    bs->num_components--;
    return 0;
}

void blend_system_free(blend_system_t *bs)
{
    if (!bs) return;
    if (bs->components) {
        free(bs->components);
        bs->components = NULL;
    }
    bs->num_components = 0;
    bs->max_components = 0;
}

/* -------------------------------------------------------------------
 * L2: Blend Computation Functions
 * ------------------------------------------------------------------- */

int blend_set_total_flow(blend_system_t *bs, double total_flow_sp)
{
    if (!bs) return -1;
    if (total_flow_sp < bs->min_total_flow || total_flow_sp > bs->max_total_flow)
        return -2;

    bs->total_flow_sp = total_flow_sp;
    return blend_compute_individual_flows(bs);
}

int blend_compute_individual_flows(blend_system_t *bs)
{
    if (!bs || !bs->components) return -1;
    if (bs->num_components == 0) return -2;

    /* Q_i = f_i * Q_total */
    for (size_t i = 0; i < bs->num_components; i++) {
        blend_component_t *bc = &bs->components[i];
        if (!bc->active) {
            bc->flow_rate = 0.0;
            continue;
        }

        double target_flow = bc->fraction_setpoint * bs->total_flow_sp;

        /* Clamp to individual flow limits */
        if (target_flow < bc->min_flow) target_flow = bc->min_flow;
        if (target_flow > bc->max_flow) target_flow = bc->max_flow;

        bc->flow_rate = target_flow;
    }

    return 0;
}

int blend_check_feasibility(const blend_system_t *bs)
{
    if (!bs) return -1;

    /* Check 1: fractions sum to ~1.0 */
    double sum = 0.0;
    for (size_t i = 0; i < bs->num_components; i++) {
        if (bs->components[i].active) {
            sum += bs->components[i].fraction_setpoint;
        }
    }
    if (fabs(sum - 1.0) > bs->fraction_tolerance) return -2;

    /* Check 2: individual fractions within bounds */
    for (size_t i = 0; i < bs->num_components; i++) {
        const blend_component_t *bc = &bs->components[i];
        if (!bc->active) continue;
        if (bc->fraction_setpoint < bc->min_fraction ||
            bc->fraction_setpoint > bc->max_fraction) {
            return -3;
        }
    }

    /* Check 3: total flow within bounds */
    if (bs->total_flow_sp < bs->min_total_flow ||
        bs->total_flow_sp > bs->max_total_flow) {
        return -4;
    }

    /* Check 4: each component's computed flow within bounds */
    for (size_t i = 0; i < bs->num_components; i++) {
        const blend_component_t *bc = &bs->components[i];
        if (!bc->active) continue;
        if (bc->flow_rate < bc->min_flow || bc->flow_rate > bc->max_flow) {
            return -5;
        }
    }

    return 0;
}

int blend_update_actual(blend_system_t *bs,
                        const double *actual_flows, size_t n)
{
    if (!bs || !actual_flows) return -1;
    if (n != bs->num_components) return -2;

    double total = 0.0;
    for (size_t i = 0; i < n; i++) {
        bs->components[i].flow_rate = actual_flows[i];
        total += actual_flows[i];
    }

    bs->total_flow_actual = total;

    /* Update actual fractions: f_i = Q_i / Q_total */
    if (total > 1e-12) {
        for (size_t i = 0; i < n; i++) {
            bs->components[i].fraction_actual = actual_flows[i] / total;
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            bs->components[i].fraction_actual = 0.0;
        }
    }

    return 0;
}

int blend_normalize_fractions(blend_system_t *bs)
{
    if (!bs) return -1;
    if (bs->num_components == 0) return -2;

    /* Compute sum of active fractions */
    double sum = 0.0;
    int n_active = 0;
    for (size_t i = 0; i < bs->num_components; i++) {
        if (bs->components[i].active) {
            sum += bs->components[i].fraction_setpoint;
            n_active++;
        }
    }

    if (n_active == 0) return -3;
    if (fabs(sum) < 1e-12) return -4;

    /* Normalize so active fractions sum to 1.0 */
    double scale = 1.0 / sum;
    for (size_t i = 0; i < bs->num_components; i++) {
        if (bs->components[i].active) {
            bs->components[i].fraction_setpoint *= scale;
        }
    }

    return 0;
}

/* -------------------------------------------------------------------
 * L2: Blend Quality and Cost
 * ------------------------------------------------------------------- */

double blend_compute_quality(const blend_system_t *bs)
{
    if (!bs || bs->num_components == 0) return NAN;

    /* Quality = sum(f_i * quality_i)  (linear blending rule) */
    double quality = 0.0;
    for (size_t i = 0; i < bs->num_components; i++) {
        if (bs->components[i].active) {
            quality += bs->components[i].fraction_actual
                     * bs->components[i].quality_value;
        }
    }
    return quality;
}

double blend_compute_cost(const blend_system_t *bs)
{
    if (!bs || bs->num_components == 0) return NAN;

    /* Cost = sum(Q_i * cost_i) / Q_total = sum(f_i * cost_i) */
    double weighted_cost = 0.0;
    for (size_t i = 0; i < bs->num_components; i++) {
        if (bs->components[i].active) {
            weighted_cost += bs->components[i].fraction_actual
                           * bs->components[i].cost_per_unit;
        }
    }
    return weighted_cost;
}

double blend_mass_balance_error(const blend_system_t *bs)
{
    if (!bs || bs->num_components == 0) return -1.0;

    double sum = 0.0;
    for (size_t i = 0; i < bs->num_components; i++) {
        sum += bs->components[i].fraction_actual;
    }
    return sum - 1.0;
}

/* -------------------------------------------------------------------
 * L5: Blend Optimization via Linear Programming
 * ------------------------------------------------------------------- */

int blend_optimize_cost(blend_system_t *bs,
                        const double *component_costs, size_t n)
{
    if (!bs || !component_costs) return -1;
    if (n != bs->num_components) return -2;

    /* Simple LP for 2-component case:
     *   min  sum(c_i * f_i)
     *   s.t. sum(f_i) = 1.0
     *        f_i_min <= f_i <= f_i_max
     *
     * For general N, use simplex or search.
     * Here we implement a greedy approach: allocate more to cheaper components.
     */

    /* Build index array sorted by cost */
    size_t *idx = (size_t *)malloc(n * sizeof(size_t));
    if (!idx) return -3;

    for (size_t i = 0; i < n; i++) idx[i] = i;

    /* Simple bubble sort on cost (n is small in practice) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (component_costs[idx[i]] > component_costs[idx[j]]) {
                size_t tmp = idx[i];
                idx[i] = idx[j];
                idx[j] = tmp;
            }
        }
    }

    /* Allocate fractions: fill cheapest first to max, then next... */
    double remaining = 1.0;
    for (size_t i = 0; i < n; i++) {
        size_t ci = idx[i];
        double min_f = bs->components[ci].min_fraction;
        double max_f = bs->components[ci].max_fraction;

        /* Subtract minimum from remaining for all */
        if (i == 0) {
            /* First pass: compute sum of minimums */
            double sum_min = 0.0;
            for (size_t j = 0; j < n; j++) {
                sum_min += bs->components[j].min_fraction;
            }
            if (sum_min > 1.0) {
                free(idx);
                return -4; /* infeasible */
            }
            remaining = 1.0 - sum_min;
            /* Set all to minimum first */
            for (size_t j = 0; j < n; j++) {
                bs->components[j].fraction_setpoint = bs->components[j].min_fraction;
            }
        }

        /* Allocate remaining to cheapest components */
        double available = max_f - min_f;
        if (available > remaining) available = remaining;
        bs->components[ci].fraction_setpoint = min_f + available;
        remaining -= available;
        if (remaining <= 1e-10) break;
    }

    free(idx);
    return 0;
}

/* -------------------------------------------------------------------
 * L5: General Linear Programming Solver for Blending
 * ------------------------------------------------------------------- */

double blend_linear_program_solve(const double *costs,
                                  const double *quality_coeffs,
                                  const double target_quality,
                                  const double *min_fracs,
                                  const double *max_fracs,
                                  double *optimal_fractions,
                                  size_t n)
{
    if (!costs || !min_fracs || !max_fracs || !optimal_fractions) return NAN;
    if (n == 0) return NAN;
    if (n > 10) return NAN; /* practical limit for enumeration */

    /* For small n, enumerate vertices of the feasible polytope.
     * Each fraction is at either min or max at vertices (corner points).
     * With quality constraint, we use a 2^n enumeration.
     */

    double best_cost = INFINITY;
    double best_fracs[10] = {0};

    unsigned long long max_combos = 1ULL << n;
    for (unsigned long long combo = 0; combo < max_combos; combo++) {
        double fracs[10];
        double sum = 0.0;
        int feasible = 1;

        for (size_t i = 0; i < n; i++) {
            int at_max = (combo >> i) & 1;
            fracs[i] = at_max ? max_fracs[i] : min_fracs[i];
            sum += fracs[i];
        }

        /* Must sum to 1.0 */
        if (fabs(sum - 1.0) > 1e-6) feasible = 0;

        /* Quality constraint (if specified) */
        if (feasible && quality_coeffs) {
            double q = 0.0;
            for (size_t i = 0; i < n; i++) {
                q += fracs[i] * quality_coeffs[i];
            }
            if (q < target_quality - 1e-6) feasible = 0;
        }

        if (feasible) {
            double cost = 0.0;
            for (size_t i = 0; i < n; i++) {
                cost += fracs[i] * costs[i];
            }
            if (cost < best_cost) {
                best_cost = cost;
                for (size_t i = 0; i < n; i++) {
                    best_fracs[i] = fracs[i];
                }
            }
        }
    }

    if (best_cost < INFINITY) {
        for (size_t i = 0; i < n; i++) {
            optimal_fractions[i] = best_fracs[i];
        }
        return best_cost;
    }

    return NAN; /* no feasible solution */
}

/* -------------------------------------------------------------------
 * L6: Recipe Management
 * ------------------------------------------------------------------- */

int blend_load_recipe(blend_system_t *bs, const blend_recipe_t *recipe)
{
    if (!bs || !recipe) return -1;
    if (!recipe->fractions || recipe->num_components == 0) return -2;

    for (size_t i = 0; i < recipe->num_components && i < bs->num_components; i++) {
        bs->components[i].fraction_setpoint = recipe->fractions[i];
    }
    bs->blend_quality_sp = recipe->target_quality;
    bs->blend_cost_per_unit = recipe->target_cost;

    return blend_normalize_fractions(bs);
}

int blend_save_recipe(const blend_system_t *bs, blend_recipe_t *recipe)
{
    if (!bs || !recipe) return -1;
    if (!recipe->fractions) return -2;

    recipe->num_components = bs->num_components;
    recipe->target_quality = bs->blend_quality_sp;
    recipe->target_cost    = bs->blend_cost_per_unit;

    for (size_t i = 0; i < bs->num_components; i++) {
        recipe->fractions[i] = bs->components[i].fraction_setpoint;
    }
    return 0;
}

/* -------------------------------------------------------------------
 * L2: Utility
 * ------------------------------------------------------------------- */

const char *blend_basis_name(blend_basis_t basis)
{
    static const char *names[] = {"Volume", "Mass", "Mole"};
    if (basis < 0 || basis >= BLEND_MODE_COUNT) return "Unknown";
    return names[basis];
}

const char *blend_strategy_name(blend_strategy_t strategy)
{
    static const char *names[] = {"Parallel", "Sequential", "In-Line", "Tank"};
    if (strategy < 0 || strategy >= BLEND_STRATEGY_COUNT) return "Unknown";
    return names[strategy];
}
