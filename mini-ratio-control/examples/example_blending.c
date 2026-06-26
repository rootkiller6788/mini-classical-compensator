#include <stdio.h>
#include <math.h>
#include "../include/ratio_control.h"
#include "../include/blending_control.h"

/**
 * @file example_blending.c
 * @brief End-to-End Gasoline Blending Ratio Control Example
 *
 * This example demonstrates multi-component blending ratio control
 * for a gasoline blending terminal. Three components (reformate,
 * alkylate, butane) are blended to meet octane and RVP specs
 * while minimizing cost.
 *
 * L6: Gasoline blending — canonical multi-component ratio problem
 * L7: Refinery blending application (ISO 8217 fuel standards)
 *
 * Knowledge points demonstrated:
 *   - Ratio station computation for each component
 *   - Mass balance verification
 *   - Cost optimization via LP
 *   - Quality constraint checking
 *   - Total flow control with individual ratio loops
 */

static void print_header(const char *title) {
    printf("\n============================================================\n");
    printf("  %s\n", title);
    printf("============================================================\n\n");
}

static void print_blend_state(const blend_system_t *bs) {
    printf("Blend System State:\n");
    printf("  Total Flow SP:   %8.1f  Actual: %8.1f\n",
           bs->total_flow_sp, bs->total_flow_actual);
    printf("  Blend Quality:   %8.1f  Cost/unit: $%6.2f\n",
           blend_compute_quality(bs), blend_compute_cost(bs));
    printf("  Mass Balance Err: %+.6f\n\n", blend_mass_balance_error(bs));

    printf("  %-15s %8s %8s %8s %8s\n",
           "Component", "Frac SP", "Frac Act", "Flow", "Quality");
    printf("  %-15s %8s %8s %8s %8s\n",
           "---------", "-------", "-------", "------", "-------");
    for (size_t i = 0; i < bs->num_components; i++) {
        printf("  %-15s %7.3f  %7.3f  %8.1f %8.1f\n",
               bs->components[i].name,
               bs->components[i].fraction_setpoint,
               bs->components[i].fraction_actual,
               bs->components[i].flow_rate,
               bs->components[i].quality_value);
    }
    printf("\n");
}

int main(void) {
    print_header("Gasoline Blending Ratio Control Example (L6/L7)");

    /* Step 1: Initialize blend system */
    blend_system_t bs;
    int rc = blend_system_init(&bs, BLEND_MODE_VOLUME,
                                BLEND_STRATEGY_INLINE, 5);
    if (rc != 0) {
        printf("ERROR: blend_system_init failed: %d\n", rc);
        return 1;
    }
    printf("[Step 1] Blend system initialized. Max components: %zu\n\n",
           bs.max_components);

    /* Step 2: Define components (reformate, alkylate, butane) */
    printf("[Step 2] Adding blend components...\n");

    /* Reformate: high octane, moderate RVP, expensive */
    blend_add_component(&bs, "Reformate", 0.50,
                        0.78, 3.20, 98.0,     /* density, $/gal, (R+M)/2 */
                        0.30, 0.60, 0.0, 5000.0);

    /* Alkylate: very high octane, low RVP, moderate cost */
    blend_add_component(&bs, "Alkylate", 0.35,
                        0.72, 2.80, 96.0,
                        0.20, 0.50, 0.0, 5000.0);

    /* Butane: high RVP, low cost, octane booster ceiling */
    blend_add_component(&bs, "n-Butane", 0.15,
                        0.58, 1.50, 93.0,
                        0.05, 0.20, 0.0, 3000.0);

    printf("  Added %zu components\n\n", bs.num_components);

    /* Step 3: Normalize fractions to sum to 1.0 */
    printf("[Step 3] Normalizing fractions (mass balance check)...\n");
    double sum_before = 0.0;
    for (size_t i = 0; i < bs.num_components; i++) {
        sum_before += bs.components[i].fraction_setpoint;
    }
    printf("  Sum before normalization: %.3f\n", sum_before);

    blend_normalize_fractions(&bs);

    double sum_after = 0.0;
    for (size_t i = 0; i < bs.num_components; i++) {
        sum_after += bs.components[i].fraction_setpoint;
    }
    printf("  Sum after normalization:  %.3f\n\n", sum_after);

    /* Step 4: Set total flow target */
    printf("[Step 4] Setting total flow target...\n");
    blend_set_total_flow(&bs, 1000.0);
    printf("  Total flow SP: %.1f barrels/day\n", bs.total_flow_sp);

    /* Compute individual flows */
    blend_compute_individual_flows(&bs);
    for (size_t i = 0; i < bs.num_components; i++) {
        printf("  %-15s target flow: %8.1f bbl/day\n",
               bs.components[i].name, bs.components[i].flow_rate);
    }
    printf("\n");

    /* Step 5: Simulate actual flows with measurement noise */
    printf("[Step 5] Simulating actual blend process...\n");
    double actual_flows[] = {498.0, 352.0, 150.0}; /* ~1% noise */
    blend_update_actual(&bs, actual_flows, 3);
    print_blend_state(&bs);

    /* Step 6: Cost optimization */
    printf("[Step 6] Running blend cost optimization...\n");
    double costs[] = {3.20, 2.80, 1.50};
    double quality_coeffs[] = {98.0, 96.0, 93.0};
    double min_fracs[] = {0.30, 0.20, 0.05};
    double max_fracs[] = {0.60, 0.50, 0.20};

    double optimal_fracs[3];
    double optimal_cost = blend_linear_program_solve(
        costs, quality_coeffs, 91.0,  /* target octane */
        min_fracs, max_fracs, optimal_fracs, 3);

    if (!isnan(optimal_cost)) {
        printf("  Optimal cost: $%.4f per unit\n", optimal_cost);
        printf("  Optimal fractions:\n");
        const char *names[] = {"Reformate", "Alkylate", "n-Butane"};
        for (int i = 0; i < 3; i++) {
            printf("    %-15s: %.4f\n", names[i], optimal_fracs[i]);
        }
    } else {
        printf("  No feasible solution found.\n");
    }
    printf("\n");

    /* Step 7: Use the built-in optimizer */
    printf("[Step 7] Built-in blend optimizer (greedy)...\n");
    blend_optimize_cost(&bs, costs, 3);
    printf("  Optimized fractions:\n");
    const char *comp_names[] = {"Reformate", "Alkylate", "n-Butane"};
    for (size_t i = 0; i < bs.num_components; i++) {
        printf("    %-15s: %.4f\n",
               comp_names[i], bs.components[i].fraction_setpoint);
    }
    printf("\n");

    /* Step 8: Feasibility check */
    printf("[Step 8] Feasibility check...\n");
    rc = blend_check_feasibility(&bs);
    printf("  Feasibility: %s (rc=%d)\n\n",
           rc == 0 ? "FEASIBLE ✓" : "INFEASIBLE ✗", rc);

    /* Step 9: Recipe save/load */
    printf("[Step 9] Recipe management...\n");
    double recipe_fracs[5] = {0, 0, 0, 0, 0};
    blend_recipe_t recipe = {"Summer Blend 2020", recipe_fracs, 0, 0.0, 0.0};
    blend_save_recipe(&bs, &recipe);
    printf("  Saved recipe: \"%s\"\n", recipe.recipe_name);
    for (size_t i = 0; i < recipe.num_components; i++) {
        printf("    %-15s: %.4f\n", comp_names[i], recipe_fracs[i]);
    }

    /* Change fractions, then reload */
    bs.components[0].fraction_setpoint = 1.0;
    blend_load_recipe(&bs, &recipe);
    printf("  Reloaded recipe — fractions restored.\n\n");

    /* Cleanup */
    blend_system_free(&bs);
    printf("[Done] Blend system freed.\n");

    printf("\n============================================================\n");
    printf("  Gasoline Blending Example — COMPLETE\n");
    printf("============================================================\n\n");
    return 0;
}
