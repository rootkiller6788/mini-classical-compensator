#ifndef BLENDING_CONTROL_H
#define BLENDING_CONTROL_H

#include "ratio_control.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file blending_control.h
 * @brief Multi-Component Blending Control — mixing N streams
 *        in prescribed proportions while respecting mass balance.
 *
 * Core equation: sum(component_i) = total_flow
 *                component_i / total_flow = setpoint_fraction_i
 *
 * L1: Blending system definitions (components, fractions, tanks)
 * L2: Blend ratio control concepts (split-range, sequential)
 * L3: Linear programming for blend optimization
 * L4: Mass balance conservation (Dalton's law of mixtures)
 * L5: Blend recipe management and setpoint calculation
 * L6: Gasoline blending, chemical mixing tank
 * L7: Refinery blending, wastewater pH neutralization
 *
 * MIT 6.302: Feedforward blending
 * Purdue ME 575: Industrial blending control
 */

typedef enum {
    BLEND_MODE_VOLUME = 0,
    BLEND_MODE_MASS,
    BLEND_MODE_MOLE,
    BLEND_MODE_COUNT
} blend_basis_t;

typedef enum {
    BLEND_STRATEGY_PARALLEL = 0,
    BLEND_STRATEGY_SEQUENTIAL,
    BLEND_STRATEGY_INLINE,
    BLEND_STRATEGY_TANK,
    BLEND_STRATEGY_COUNT
} blend_strategy_t;

typedef struct {
    const char *name;
    double      fraction_setpoint;
    double      fraction_actual;
    double      flow_rate;
    double      density;
    double      cost_per_unit;
    double      quality_value;
    double      min_fraction;
    double      max_fraction;
    double      min_flow;
    double      max_flow;
    int         active;
} blend_component_t;

typedef struct {
    blend_basis_t       basis;
    blend_strategy_t    strategy;
    blend_component_t  *components;
    size_t              num_components;
    size_t              max_components;
    double              total_flow_sp;
    double              total_flow_actual;
    double              blend_quality_sp;
    double              blend_quality_actual;
    double              blend_cost_per_unit;
    double              fraction_tolerance;
    double              mass_balance_tolerance;
    double              min_total_flow;
    double              max_total_flow;
    ratio_loop_t       *ratio_loops;
    size_t              num_ratio_loops;
} blend_system_t;

typedef struct {
    const char  *recipe_name;
    double      *fractions;
    size_t       num_components;
    double       target_quality;
    double       target_cost;
} blend_recipe_t;

int  blend_system_init(blend_system_t *bs, blend_basis_t basis,
                       blend_strategy_t strategy, size_t max_components);
int  blend_add_component(blend_system_t *bs, const char *name,
                         double fraction_sp, double density,
                         double cost, double quality,
                         double min_frac, double max_frac,
                         double min_flow, double max_flow);
int  blend_remove_component(blend_system_t *bs, size_t index);
int  blend_set_total_flow(blend_system_t *bs, double total_flow_sp);
int  blend_compute_individual_flows(blend_system_t *bs);
int  blend_check_feasibility(const blend_system_t *bs);
int  blend_update_actual(blend_system_t *bs,
                         const double *actual_flows, size_t n);
int  blend_normalize_fractions(blend_system_t *bs);
double blend_compute_quality(const blend_system_t *bs);
double blend_compute_cost(const blend_system_t *bs);
double blend_mass_balance_error(const blend_system_t *bs);
int  blend_optimize_cost(blend_system_t *bs,
                         const double *component_costs, size_t n);
int  blend_load_recipe(blend_system_t *bs, const blend_recipe_t *recipe);
int  blend_save_recipe(const blend_system_t *bs, blend_recipe_t *recipe);
void blend_system_free(blend_system_t *bs);

const char *blend_basis_name(blend_basis_t basis);
const char *blend_strategy_name(blend_strategy_t strategy);
double blend_linear_program_solve(const double *costs,
                                  const double *quality_coeffs,
                                  const double target_quality,
                                  const double *min_fracs,
                                  const double *max_fracs,
                                  double *optimal_fractions,
                                  size_t n);

#ifdef __cplusplus
}
#endif

#endif /* BLENDING_CONTROL_H */
