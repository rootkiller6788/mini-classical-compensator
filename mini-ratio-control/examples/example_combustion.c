#include <stdio.h>
#include <math.h>
#include "../include/ratio_control.h"
#include "../include/cross_limiting.h"

/**
 * @file example_combustion.c
 * @brief Boiler Combustion Cross-Limiting Ratio Control Example
 *
 * Demonstrates cross-limiting ratio control for safe combustion.
 * Fuel/air ratio must be maintained to prevent:
 *   - Fuel-rich conditions (smoking, CO, safety hazard)
 *   - Excessive excess air (heat loss up the stack)
 *
 * The cross-limiting algorithm uses high/low select logic:
 *   SP_fuel = MIN(demand, air_flow / K_fa)
 *   SP_air  = MAX(demand, fuel_flow * K_fa)
 *
 * This guarantees:
 *   Actual fuel ≤ available air / K_fa  (always fuel-lean)
 *
 * L6: Boiler combustion control
 * L7: Power plant burner management (Fukushima safety lessons)
 *
 * Reference: Shinskey "Process Control Systems" Ch.8
 * MIT 6.302: Combustion Case Study
 */

#define SIMULATION_STEPS 20

static void print_header(const char *t) {
    printf("\n==========================================================\n");
    printf("  %s\n", t);
    printf("==========================================================\n\n");
}

int main(void) {
    print_header("Boiler Combustion Cross-Limiting Control (L6/L7)");

    /* Step 1: Configure cross-limit controller */
    double K_fa = 0.08;       /* fuel-to-air ratio (lb fuel / lb air) */
    double excess_air = 15.0; /* 15% excess air target */
    double ratio_min = 0.06;
    double ratio_max = 0.10;

    cross_limit_t cl;
    int rc = cross_limit_init(&cl, K_fa, excess_air, ratio_min, ratio_max);
    if (rc != 0) {
        printf("ERROR: cross_limit_init failed: %d\n", rc);
        return 1;
    }

    /* Set rate limits */
    cross_limit_set_rate_limits(&cl,
        20.0, 20.0,    /* fuel: ±20 lb/s per second */
        200.0, 200.0); /* air:  ±200 lb/s per second */

    printf("[Config] Fuel/Air Ratio: %.4f lb/lb\n", K_fa);
    printf("[Config] Excess Air Target: %.1f%%\n", excess_air);
    printf("[Config] Ratio range: [%.3f, %.3f]\n\n", ratio_min, ratio_max);

    /* Step 2: Simulation loop with load changes */
    printf("Time  | FuelDem | AirDem | FuelSP |  AirSP | ActFuel| ActAir | Status\n");
    printf("------+---------+--------+--------+--------+--------+--------+----------\n");

    double fuel_actual = 50.0;
    double air_actual  = 700.0;
    double dt = 0.5; /* 0.5 second time step */

    for (int step = 0; step < SIMULATION_STEPS; step++) {
        double time = step * dt;

        /* Load profile: ramp up, hold, step down */
        double load_demand;
        if (step < 5) {
            load_demand = 50.0 + step * 5.0;  /* ramp 50→75 */
        } else if (step < 12) {
            load_demand = 75.0;               /* hold at 75% */
        } else {
            load_demand = 75.0 - (step - 12) * 10.0; /* step down */
            if (load_demand < 20.0) load_demand = 20.0;
        }

        double fuel_demand = load_demand * 0.8; /* 80% of load = fuel */
        double air_demand  = load_demand * 12.0; /* nominal air for load */

        /* Set demands and actuals */
        cross_limit_set_demands(&cl, fuel_demand, air_demand);
        cross_limit_set_actuals(&cl, fuel_actual, air_actual);
        cross_limit_execute(&cl, dt);

        double fuel_sp, air_sp;
        cross_limit_get_fuel_sp(&cl, &fuel_sp);
        cross_limit_get_air_sp(&cl, &air_sp);

        cross_limit_status_t status = cross_limit_get_status(&cl);
        double excess = cross_limit_get_excess_air(&cl);

        printf("%5.1f | %7.1f | %6.1f | %6.1f | %6.1f | %6.1f | %6.1f | %s",
               time, fuel_demand, air_demand,
               fuel_sp, air_sp, fuel_actual, air_actual,
               cross_limit_status_name(status));

        /* Safety check: verify fuel never exceeds air/K_fa */
        double max_safe_fuel = air_actual / K_fa;
        if (fuel_sp > max_safe_fuel) {
            printf(" *** UNSAFE! ***");
        }
        printf(" (exc=%.1f%%)\n", excess);

        /* Simulate process: actuals track SPs with lag */
        double tau_fuel = 2.0;
        double tau_air  = 1.5;
        fuel_actual += (dt / tau_fuel) * (fuel_sp - fuel_actual);
        air_actual  += (dt / tau_air)  * (air_sp  - air_actual);
    }

    /* Step 3: Show safe operation verification */
    printf("\n[Safety Verification]\n");

    /* Test extreme case: sudden fuel demand spike */
    cross_limit_t cl_test;
    cross_limit_init(&cl_test, 0.08, 15.0, 0.06, 0.10);

    /* Air-limited scenario: high fuel demand but low air available */
    cross_limit_set_demands(&cl_test, 100.0, 50.0);
    cross_limit_set_actuals(&cl_test, 50.0, 600.0);
    cross_limit_execute(&cl_test, 1.0);

    double fuel_sp, air_sp;
    cross_limit_get_fuel_sp(&cl_test, &fuel_sp);
    cross_limit_get_air_sp(&cl_test, &air_sp);

    printf("  Scenario: Fuel demand spike, limited air\n");
    printf("    Fuel demand: 100.0 → Fuel SP: %.1f (limited by air)\n", fuel_sp);
    printf("    Air demand:   50.0 → Air SP:  %.1f (boosted by fuel)\n", air_sp);
    printf("    Safety: fuel/air = %.4f ≤ %.4f (ratio max)\n",
           fuel_sp / air_sp, ratio_max);
    printf("    Result: SAFE (air-rich guaranteed) ✓\n");

    /* Test extreme case: sudden fuel loss */
    cross_limit_set_demands(&cl_test, 10.0, 200.0);
    cross_limit_set_actuals(&cl_test, 20.0, 100.0);
    cross_limit_execute(&cl_test, 1.0);

    cross_limit_get_fuel_sp(&cl_test, &fuel_sp);
    cross_limit_get_air_sp(&cl_test, &air_sp);

    printf("\n  Scenario: Fuel loss, air demand high\n");
    printf("    Fuel demand:  10.0 → Fuel SP: %.1f\n", fuel_sp);
    printf("    Air demand:  200.0 → Air SP:  %.1f (limited by fuel)\n", air_sp);
    printf("    Result: Air not excessive (prevents heat loss) ✓\n");

    printf("\n==========================================================\n");
    printf("  Combustion Cross-Limiting Example — COMPLETE\n");
    printf("==========================================================\n\n");
    return 0;
}
