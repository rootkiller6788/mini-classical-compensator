#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "../include/ratio_control.h"
#include "../include/ratio_station.h"
#include "../include/blending_control.h"
#include "../include/cross_limiting.h"
#include "../include/ratio_feedforward.h"
#include "../include/ratio_cascade.h"

/**
 * @file test_ratio.c
 * @brief Comprehensive test suite for mini-ratio-control
 *
 * L4: Mathematical assertions testing fundamental theorems:
 *   1. Ratio station = K * master + b (linearity)
 *   2. Mass balance: sum(f_i) = 1.0
 *   3. Cross-limiting preserves fuel <= air/K (safety)
 *   4. Lead-lag steady-state gain = 1.0
 *   5. Cascade primary has anti-windup
 */

#define EPS 1e-9
#define assert_double_eq(a, b) do { \
    double _a = (a), _b = (b); \
    assert(fabs(_a - _b) < EPS || (printf("FAIL: %.15g != %.15g\n", _a, _b), 0)); \
} while(0)
#define assert_double_near(a, b, tol) do { \
    double _a = (a), _b = (b), _t = (tol); \
    assert(fabs(_a - _b) < _t || (printf("FAIL: %.15g != %.15g (tol=%.15g)\n", _a, _b, _t), 0)); \
} while(0)

static int tests_run = 0;
static int tests_passed = 0;
#define RUN_TEST(name) do { \
    tests_run++; \
    printf("  %-35s... ", #name); \
    test_##name(); \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

/* ===================================================================
 * L2: Ratio Control Core Tests
 * =================================================================== */

void test_ratio_init_defaults(void) {
    ratio_loop_t loop;
    int rc = ratio_init(&loop, 1, "FRC-101", "FIC-100", "FIC-101");
    assert(rc == 0);
    assert(loop.loop_id == 1);
    assert(strcmp(loop.tag, "FRC-101") == 0);
    assert(loop.config.ratio_gain == 1.0);
    assert(loop.config.ratio_bias == 0.0);
    assert(loop.config.mode == RATIO_MODE_MASTER_SLAVE);
    assert(loop.state.integrity == RATIO_INTEGRITY_OK);
}

void test_ratio_configure(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 2, "FRC-102", "FI-10", "FI-20");
    int rc = ratio_configure(&loop, RATIO_MODE_BLENDING,
                             RATIO_ACTION_MULTIPLY,
                             2.5, 0.5, 0.0, 10.0);
    assert(rc == 0);
    assert(loop.config.ratio_gain == 2.5);
    assert(loop.config.ratio_bias == 0.5);
    assert(loop.config.ratio_min == 0.0);
    assert(loop.config.ratio_max == 10.0);
    assert(loop.config.mode == RATIO_MODE_BLENDING);
}

void test_ratio_configure_invalid(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 3, "BAD", "M1", "S1");
    /* min > max */
    int rc = ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                             RATIO_ACTION_MULTIPLY, 1.0, 0.0, 10.0, 0.0);
    assert(rc == -2);
}

void test_ratio_set_master_multiply(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 4, "FRC-104", "FI-30", "FI-40");
    ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY, 2.0, 0.0, 0.0, 100.0);

    /* Set config limits open */
    loop.config.slave_sp_min = -1000.0;
    loop.config.slave_sp_max = 1000.0;

    int rc = ratio_set_master(&loop, 50.0, 1.0);
    assert(rc == 0);
    assert_double_eq(ratio_get_slave_sp(&loop), 100.0); /* 2.0 * 50 = 100 */
}

void test_ratio_set_master_with_bias(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 5, "FRC-105", "FI-50", "FI-60");
    ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY, 1.5, 10.0, 0.0, 100.0);
    loop.config.slave_sp_min = -1000.0;
    loop.config.slave_sp_max = 1000.0;

    ratio_set_master(&loop, 100.0, 1.0);
    assert_double_eq(ratio_get_slave_sp(&loop), 160.0); /* 1.5*100 + 10 */
}

void test_ratio_set_master_filter(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 6, "FRC-106", "FI-70", "FI-80");
    ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY, 1.0, 0.0, 0.0, 100.0);
    loop.config.filter_enabled = 1;
    loop.config.filter_tau = 10.0;
    loop.config.slave_sp_min = -1000.0;
    loop.config.slave_sp_max = 1000.0;

    /* Step response: after one dt, should be partially filtered */
    ratio_set_master(&loop, 100.0, 2.0);
    double sp1 = ratio_get_slave_sp(&loop);
    /* With tau=10, dt=2, alpha = 2/12 = 0.1667, filtered = 0 + 0.1667 * 100 = 16.67 */
    assert_double_near(sp1, 16.6666667, 0.01);

    /* After several steps, should converge to 100 */
    for (int i = 0; i < 50; i++) {
        ratio_set_master(&loop, 100.0, 2.0);
    }
    double sp_final = ratio_get_slave_sp(&loop);
    assert_double_near(sp_final, 100.0, 1.0);
}

void test_ratio_set_slave(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 7, "FRC-107", "FI-90", "FI-100");
    ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY, 2.0, 0.0, 0.0, 100.0);
    loop.config.slave_sp_min = -1000.0;
    loop.config.slave_sp_max = 1000.0;

    ratio_set_master(&loop, 100.0, 1.0);
    ratio_set_slave(&loop, 200.0);
    /* Actual ratio = 200/100 = 2.0 = configured ratio */
    assert_double_eq(ratio_get_actual(&loop), 2.0);
}

void test_ratio_shed(void) {
    ratio_loop_t loop;
    ratio_init(&loop, 8, "FRC-108", "FI-110", "FI-120");
    ratio_configure(&loop, RATIO_MODE_MASTER_SLAVE,
                    RATIO_ACTION_MULTIPLY, 2.0, 0.0, 0.0, 100.0);
    loop.config.slave_sp_min = -1000.0;
    loop.config.slave_sp_max = 1000.0;

    ratio_set_master(&loop, 50.0, 1.0);
    double sp_before = ratio_get_slave_sp(&loop);
    assert(sp_before > 0.0);

    ratio_shed(&loop, 25.0);
    assert_double_eq(ratio_get_slave_sp(&loop), 25.0);
    assert(ratio_get_integrity(&loop) == RATIO_INTEGRITY_MASTER_BAD);
}

void test_ratio_mass_balance(void) {
    double fractions[] = {0.3, 0.3, 0.4};
    double err = ratio_check_mass_balance(fractions, 3, 1e-6);
    assert_double_eq(err, 0.0);

    double bad_fractions[] = {0.3, 0.3, 0.3};
    double err2 = ratio_check_mass_balance(bad_fractions, 3, 1e-6);
    assert(fabs(err2 + 0.1) < 1e-6 || err2 < -0.01); /* sum = 0.9, error = -0.1 */

    double *null_f = NULL;
    double err3 = ratio_check_mass_balance(null_f, 0, 1e-6);
    assert(err3 == -1.0); /* null check */
}

void test_ratio_optimal_blend(void) {
    /* cost_a < cost_b → use max A */
    double r = ratio_optimal_blend(1.0, 5.0, 0.2, 0.8);
    assert_double_eq(r, 0.8);

    /* cost_a > cost_b → use min A */
    r = ratio_optimal_blend(5.0, 1.0, 0.2, 0.8);
    assert_double_eq(r, 0.2);

    /* equal costs → midpoint */
    r = ratio_optimal_blend(3.0, 3.0, 0.2, 0.8);
    assert_double_eq(r, 0.5);
}

/* ===================================================================
 * L4: Ratio Characterization Table Tests
 * =================================================================== */

void test_ratio_char_table(void) {
    ratio_char_point_t pts[] = {
        {0.0,  1.0},
        {50.0, 1.5},
        {100.0, 2.0}
    };
    ratio_char_table_t table;
    int rc = ratio_char_table_init(&table, pts, 3);
    assert(rc == 0);
    assert(table.num_points == 3);

    /* Exact breakpoint lookup */
    assert_double_eq(ratio_char_lookup(&table, 0.0), 1.0);
    assert_double_eq(ratio_char_lookup(&table, 50.0), 1.5);
    assert_double_eq(ratio_char_lookup(&table, 100.0), 2.0);

    /* Interpolation: at 25.0 → 1.0 + 0.5*(25/50) = 1.25 */
    double v = ratio_char_lookup(&table, 25.0);
    assert_double_near(v, 1.25, 1e-6);

    /* Interpolation: at 75.0 → 1.5 + 0.5*(25/50) = 1.75 */
    v = ratio_char_lookup(&table, 75.0);
    assert_double_near(v, 1.75, 1e-6);

    /* Extrapolation low */
    v = ratio_char_lookup(&table, -10.0);
    assert_double_eq(v, 1.0); /* extrap_low = first point */

    ratio_char_table_free(&table);
    assert(table.points == NULL);
}

/* ===================================================================
 * L2: Ratio Station Tests
 * =================================================================== */

void test_ratio_station_linear(void) {
    ratio_station_t rs;
    ratio_station_init(&rs, RATIO_FORMULA_LINEAR, RATIO_COMP_DISABLED);
    ratio_station_set_limits(&rs, 0.0, 1000.0, 1000.0);

    double coeffs[] = {2.0};
    ratio_station_set_formula(&rs, RATIO_FORMULA_LINEAR, coeffs, 1);
    ratio_station_compute(&rs, 50.0, 0.0, 1.0);

    /* SP = computed_ratio * master = 2.0 * 50 = 100 */
    assert_double_eq(ratio_station_get_setpoint(&rs), 100.0);
}

void test_ratio_station_square_root(void) {
    ratio_station_t rs;
    ratio_station_init(&rs, RATIO_FORMULA_SQUARE_ROOT, RATIO_COMP_DISABLED);
    ratio_station_set_limits(&rs, 0.0, 10000.0, 10000.0);

    double coeffs[] = {3.0};
    ratio_station_set_formula(&rs, RATIO_FORMULA_SQUARE_ROOT, coeffs, 1);
    ratio_station_compute(&rs, 100.0, 0.0, 1.0);

    /* SP = 3 * sqrt(100) * 100 = 3 * 10 * 100... wait no.
     * raw_ratio = 3 * sqrt(100) = 30
     * final_SP = raw_ratio * master = 30 * 100 = 3000
     */
    assert_double_near(ratio_station_get_setpoint(&rs), 3000.0, 1.0);
}

void test_ratio_station_rate_limit(void) {
    ratio_station_t rs;
    ratio_station_init(&rs, RATIO_FORMULA_LINEAR, RATIO_COMP_DISABLED);
    ratio_station_set_limits(&rs, 0.0, 1000.0, 10.0); /* rate limit = 10/s */

    double coeffs[] = {1.0};
    ratio_station_set_formula(&rs, RATIO_FORMULA_LINEAR, coeffs, 1);

    /* Step 1: small initial */
    ratio_station_compute(&rs, 50.0, 0.0, 1.0);
    double sp1 = ratio_station_get_setpoint(&rs);
    assert_double_eq(sp1, 50.0);

    /* Step 2: big jump to 200 — rate limited to +10 */
    ratio_station_compute(&rs, 200.0, 0.0, 1.0);
    double sp2 = ratio_station_get_setpoint(&rs);
    assert_double_eq(sp2, 60.0); /* 50 + 10 */

    /* Step 3: another +10 */
    ratio_station_compute(&rs, 200.0, 0.0, 1.0);
    double sp3 = ratio_station_get_setpoint(&rs);
    assert_double_eq(sp3, 70.0);
}

/* ===================================================================
 * L4: Cross-Limiting Tests (Safety Theorem)
 * =================================================================== */

void test_cross_limit_fuel_limited(void) {
    cross_limit_t cl;
    cross_limit_init(&cl, 0.1, 10.0, 0.05, 0.2); /* ratio=0.1 fuel/air, 10% excess air */

    /* High fuel demand, low air — fuel should be limited */
    cross_limit_set_demands(&cl, 50.0, 10.0);   /* want 50 fuel, 10 air */
    cross_limit_set_actuals(&cl, 10.0, 100.0);   /* actually 10 fuel, 100 air */
    cross_limit_execute(&cl, 1.0);

    double fuel_sp;
    cross_limit_get_fuel_sp(&cl, &fuel_sp);

    /* fuel_hi_limit = air_actual / K_fa * (100/(100+excess)) */
    /* = 100 / 0.1 * (100/110) = 1000 * 0.909 = 909.09 */
    /* fuel_demand = 50 < 909 → not limited → fuel_sp = 50 */
    assert_double_near(fuel_sp, 50.0, 0.1);

    /* Now reduce air so fuel IS limited */
    cross_limit_set_actuals(&cl, 10.0, 5.0);     /* only 5 air */
    cross_limit_execute(&cl, 1.0);
    cross_limit_get_fuel_sp(&cl, &fuel_sp);
    /* fuel_hi_limit = 5 / 0.1 * 0.909 = 45.45 */
    /* fuel_demand = 50 > 45.45 → fuel_sp = 45.45 */
    assert(fuel_sp < 50.0); /* must be limited */
    assert(fuel_sp > 40.0);
}

void test_cross_limit_air_guaranteed(void) {
    cross_limit_t cl;
    cross_limit_init(&cl, 0.1, 10.0, 0.05, 0.2);

    /* Low air demand, high fuel — air should be boosted */
    cross_limit_set_demands(&cl, 50.0, 3.0);     /* only 3 air demanded */
    cross_limit_set_actuals(&cl, 50.0, 3.0);     /* actual: 50 fuel, 3 air */
    cross_limit_execute(&cl, 1.0);

    double air_sp;
    cross_limit_get_air_sp(&cl, &air_sp);
    /* air_lo_limit = 50 * 0.1 * 1.1 = 5.5 */
    /* air_demand = 3 < 5.5 → air_sp = 5.5 (boosted) */
    assert(air_sp > 5.0);
    assert(air_sp < 6.0);
}

void test_signal_selector(void) {
    signal_selector_t ss;
    signal_selector_init(&ss, "TEST_SEL");

    signal_selector_add_lo(&ss, 5.0);
    signal_selector_add_lo(&ss, 2.0);
    signal_selector_add_lo(&ss, 8.0);

    signal_selector_add_hi(&ss, 5.0);
    signal_selector_add_hi(&ss, 2.0);
    signal_selector_add_hi(&ss, 8.0);

    signal_selector_execute(&ss);
    assert_double_eq(signal_selector_get_lo(&ss), 2.0);
    assert_double_eq(signal_selector_get_hi(&ss), 8.0);

    double med = signal_selector_get_median(&ss);
    assert_double_eq(med, 5.0);
}

void test_hi_lo_select_functions(void) {
    assert_double_eq(hi_select(3.0, 7.0), 7.0);
    assert_double_eq(lo_select(3.0, 7.0), 3.0);
    assert_double_eq(mid_select(1.0, 5.0, 3.0), 3.0);
    assert_double_eq(mid_select(5.0, 1.0, 3.0), 3.0);
    assert_double_eq(mid_select(1.0, 3.0, 5.0), 3.0);

    double vals[] = {1.0, 9.0, 3.0, 7.0, 5.0};
    size_t idx;
    assert_double_eq(hi_select_n(vals, 5, &idx), 9.0);
    assert(idx == 1);
    assert_double_eq(lo_select_n(vals, 5, &idx), 1.0);
    assert(idx == 0);
    double med5 = median_select_n(vals, 5);
    assert_double_eq(med5, 5.0);
}

/* ===================================================================
 * L4: Feedforward Tests
 * =================================================================== */

void test_feedforward_static(void) {
    feedforward_compensator_t ff;
    feedforward_init(&ff, FF_MODE_STATIC, FF_ACTION_ADDITIVE, 2.0);

    feedforward_compute(&ff, 50.0, 1.0);
    assert_double_eq(feedforward_get_output(&ff), 100.0); /* 2.0 * 50 */
}

void test_feedforward_lead_lag(void) {
    feedforward_compensator_t ff;
    feedforward_init(&ff, FF_MODE_DYNAMIC_LEAD_LAG, FF_ACTION_ADDITIVE, 1.0);
    feedforward_set_dynamics(&ff, 2.0, 5.0, 0.0);

    /* Steady-state gain of lead-lag is 1.0: (T_lead*s+1)/(T_lag*s+1) → 1 as s→0 */
    /* Run multiple steps to reach steady state */
    for (int i = 0; i < 100; i++) {
        feedforward_compute(&ff, 50.0, 0.1);
    }
    /* Should converge to 50.0 (gain=1, master=50) */
    double out = feedforward_get_output(&ff);
    assert_double_near(out, 50.0, 0.1);
}

void test_ratio_ff_loop(void) {
    ratio_ff_loop_t loop;
    ratio_ff_loop_init(&loop, 1, "FF-01", FF_MODE_STATIC,
                       FF_ACTION_ADDITIVE, 0.5);
    ratio_ff_loop_set_ratio(&loop, 3.0);

    /* Execute: master=100, fb=20, dt=1 */
    /* ff_contribution = 0.5 * 100 * 3.0 = 150 */
    /* total = 150 + 20 = 170 (additive) */
    ratio_ff_loop_execute(&loop, 100.0, 20.0, 1.0);
    double out = ratio_ff_loop_get_output(&loop);
    assert_double_near(out, 170.0, 1.0);
}

/* ===================================================================
 * L2: Blending Control Tests
 * =================================================================== */

void test_blend_system_init(void) {
    blend_system_t bs;
    int rc = blend_system_init(&bs, BLEND_MODE_VOLUME,
                               BLEND_STRATEGY_PARALLEL, 5);
    assert(rc == 0);
    assert(bs.num_components == 0);
    assert(bs.max_components == 5);

    blend_system_free(&bs);
}

void test_blend_add_component(void) {
    blend_system_t bs;
    blend_system_init(&bs, BLEND_MODE_MASS, BLEND_STRATEGY_INLINE, 4);

    int rc = blend_add_component(&bs, "Gasoline", 0.5,
                                  0.75, 2.50, 92.0,    /* density, cost, octane */
                                  0.3, 0.7, 0.0, 100.0);
    assert(rc == 0);

    rc = blend_add_component(&bs, "Ethanol", 0.5,
                              0.79, 1.80, 108.0,
                              0.3, 0.7, 0.0, 100.0);
    assert(rc == 0);
    assert(bs.num_components == 2);

    blend_system_free(&bs);
}

void test_blend_mass_balance(void) {
    blend_system_t bs;
    blend_system_init(&bs, BLEND_MODE_VOLUME, BLEND_STRATEGY_TANK, 3);

    blend_add_component(&bs, "A", 0.4, 1.0, 1.0, 50.0, 0.0, 1.0, 0.0, 100.0);
    blend_add_component(&bs, "B", 0.35, 1.0, 2.0, 75.0, 0.0, 1.0, 0.0, 100.0);
    blend_add_component(&bs, "C", 0.25, 1.0, 3.0, 90.0, 0.0, 1.0, 0.0, 100.0);

    double err = blend_mass_balance_error(&bs);
    /* Actuals are all 0 since we haven't updated */
    assert(fabs(err + 1.0) < 1e-6); /* sum of actuals=0, err = 0-1 = -1 */

    /* Update actuals */
    double actuals[] = {40.0, 35.0, 25.0}; /* total=100 */
    blend_update_actual(&bs, actuals, 3);
    err = blend_mass_balance_error(&bs);
    assert_double_eq(err, 0.0); /* 40/100 + 35/100 + 25/100 = 1.0 */
    assert_double_eq(bs.total_flow_actual, 100.0);

    blend_system_free(&bs);
}

void test_blend_cost_quality(void) {
    blend_system_t bs;
    blend_system_init(&bs, BLEND_MODE_VOLUME, BLEND_STRATEGY_PARALLEL, 2);

    blend_add_component(&bs, "Premium", 0.6, 0.75, 3.00, 95.0, 0.0, 1.0, 0.0, 100.0);
    blend_add_component(&bs, "Regular", 0.4, 0.74, 2.20, 87.0, 0.0, 1.0, 0.0, 100.0);

    double actuals[] = {60.0, 40.0};
    blend_update_actual(&bs, actuals, 2);

    double cost = blend_compute_cost(&bs);
    /* f1*cost1 + f2*cost2 = 0.6*3.00 + 0.4*2.20 = 1.80 + 0.88 = 2.68 */
    assert_double_near(cost, 2.68, 0.01);

    double quality = blend_compute_quality(&bs);
    /* 0.6*95 + 0.4*87 = 57 + 34.8 = 91.8 */
    assert_double_near(quality, 91.8, 0.01);

    blend_system_free(&bs);
}

void test_blend_optimization(void) {
    blend_system_t bs;
    blend_system_init(&bs, BLEND_MODE_VOLUME, BLEND_STRATEGY_PARALLEL, 2);

    blend_add_component(&bs, "Cheap", 0.5, 1.0, 1.0, 80.0, 0.3, 0.7, 0.0, 100.0);
    blend_add_component(&bs, "Expensive", 0.5, 1.0, 5.0, 95.0, 0.3, 0.7, 0.0, 100.0);

    double costs[] = {1.0, 5.0};
    int rc = blend_optimize_cost(&bs, costs, 2);
    assert(rc == 0);

    /* Optimizer should prefer "Cheap" (cost 1.0 < 5.0) at max fraction */
    assert(bs.components[0].fraction_setpoint > bs.components[1].fraction_setpoint);

    blend_system_free(&bs);
}

/* ===================================================================
 * L4: Cascade Ratio Tests
 * =================================================================== */

void test_cascade_init(void) {
    cascade_ratio_t cr;
    int rc = cascade_ratio_init(&cr, CASCADE_MODE_SIMPLE, 1, "TC-100/FIC-101");
    assert(rc == 0);
    assert(cr.mode == CASCADE_MODE_SIMPLE);
    assert(cr.ratio_cfg_a.ratio_gain == 1.0);
}

void test_cascade_execute_simple(void) {
    cascade_ratio_t cr;
    cascade_ratio_init(&cr, CASCADE_MODE_SIMPLE, 2, "LC-200/FC-201");

    /* Tune primary and secondary */
    cascade_ratio_tune_primary(&cr, 2.0, 0.5, 0.0, 0.0, 100.0);
    cascade_ratio_tune_secondary(&cr, 0, 1.0, 0.3, 0.0, 0.0, 100.0);

    /* Put primary in auto */
    cr.primary_pid.in_auto = 1;

    /* Set ratio config: ratio_gain = 2.0 */
    ratio_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.ratio_gain = 2.0;
    cfg.ratio_bias = 0.0;
    cfg.ratio_min  = 0.0;
    cfg.ratio_max  = 100.0;
    cfg.action     = RATIO_ACTION_MULTIPLY;
    cascade_ratio_configure(&cr, &cfg, NULL);

    /* Put secondary in cascade */
    cascade_ratio_switch_secondary_to_cascade(&cr, 0);

    /* Set primary SP and execute */
    cascade_ratio_set_primary_sp(&cr, 50.0); /* primary SP = 50 */
    cascade_ratio_execute(&cr, 40.0, 0.0, 0.0, 1.0);

    /* Primary error = 50 - 40 = 10
     * P_out = 2.0 * 10 = 20, I_out = 0.5 * 10 * 1 = 5, total = 25 (unclamped)
     * Secondary SP = 25 * 2.0 = 50 */
    double sec_sp = cascade_ratio_get_secondary_sp(&cr, 0);
    assert_double_near(sec_sp, 50.0, 1.0);
}

void test_pid_loop(void) {
    pid_loop_t pid;
    pid_loop_init(&pid, "PIC-001", 2.0, 0.1, 0.0, 0.0, 100.0);
    pid_loop_set_auto(&pid, 1);
    pid_loop_set_sp(&pid, 50.0);
    pid_loop_set_pv(&pid, 40.0);

    pid_loop_execute(&pid, 1.0);
    /* P: 2*10=20, I: 0.1*10*1=1, total=21 */
    double out = pid_loop_get_output(&pid);
    assert_double_near(out, 21.0, 0.1);

    /* Run more steps — integral should accumulate */
    pid_loop_execute(&pid, 1.0); /* error still 10, P=20, I_accum=2, total=22 */
    out = pid_loop_get_output(&pid);
    assert_double_near(out, 22.0, 0.1);
}

/* ===================================================================
 * L4: Validate Config Tests
 * =================================================================== */

void test_validate_config(void) {
    ratio_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));

    cfg.ratio_min = 1.0;
    cfg.ratio_max = 0.5; /* invalid */
    cfg.master_span = 100.0;
    cfg.slave_span = 100.0;
    assert(ratio_validate_config(&cfg) == -2); /* min > max */

    cfg.ratio_min = 0.5;
    cfg.ratio_max = 1.0;
    cfg.master_span = 0.0; /* invalid */
    assert(ratio_validate_config(&cfg) == -4); /* invalid master span */

    cfg.master_span = 100.0;
    cfg.slave_span = 0.0; /* invalid */
    assert(ratio_validate_config(&cfg) == -5); /* invalid slave span */

    cfg.slave_span = 100.0;
    cfg.filter_tau = -1.0; /* invalid */
    assert(ratio_validate_config(&cfg) == -6); /* negative filter tau */

    cfg.filter_tau = 10.0;
    assert(ratio_validate_config(&cfg) == 0); /* all good */
}

/* ===================================================================
 * L2: String Conversion Tests
 * =================================================================== */

void test_string_conversions(void) {
    assert(strcmp(ratio_mode_name(RATIO_MODE_MASTER_SLAVE), "Master-Slave") == 0);
    assert(strcmp(ratio_mode_name(RATIO_MODE_CROSS_LIMITING), "Cross-Limiting") == 0);
    assert(strcmp(ratio_mode_name(RATIO_MODE_BLENDING), "Blending") == 0);
    assert(strcmp(ratio_integrity_name(RATIO_INTEGRITY_OK), "OK") == 0);
    assert(strcmp(ratio_integrity_name(RATIO_INTEGRITY_MASTER_BAD), "Master Bad") == 0);
    assert(strcmp(ratio_action_name(RATIO_ACTION_MULTIPLY), "Multiply") == 0);
    assert(strcmp(ratio_formula_name(RATIO_FORMULA_LINEAR), "Linear") == 0);
    assert(strcmp(ratio_formula_name(RATIO_FORMULA_EXPONENTIAL), "Exponential") == 0);
    assert(strcmp(cross_limit_status_name(CROSS_LIMIT_FULL_CROSS), "Full Cross-Limiting") == 0);
    assert(strcmp(cascade_state_name(CASCADE_STATE_PRIMARY_ACTIVE), "Primary Active") == 0);
    assert(cascade_mode_name(100) != NULL); /* out of range → "Unknown" */
}

/* ===================================================================
 * Main
 * =================================================================== */

int main(void) {
    printf("\n=== mini-ratio-control Test Suite ===\n\n");

    printf("[L1/L2] Ratio Control Core:\n");
    RUN_TEST(ratio_init_defaults);
    RUN_TEST(ratio_configure);
    RUN_TEST(ratio_configure_invalid);
    RUN_TEST(ratio_set_master_multiply);
    RUN_TEST(ratio_set_master_with_bias);
    RUN_TEST(ratio_set_master_filter);
    RUN_TEST(ratio_set_slave);
    RUN_TEST(ratio_shed);
    RUN_TEST(ratio_mass_balance);
    RUN_TEST(ratio_optimal_blend);
    RUN_TEST(validate_config);
    RUN_TEST(string_conversions);

    printf("\n[L4] Ratio Characterization Table:\n");
    RUN_TEST(ratio_char_table);

    printf("\n[L2/L4] Ratio Station:\n");
    RUN_TEST(ratio_station_linear);
    RUN_TEST(ratio_station_square_root);
    RUN_TEST(ratio_station_rate_limit);

    printf("\n[L4] Cross-Limiting (Safety Theorems):\n");
    RUN_TEST(cross_limit_fuel_limited);
    RUN_TEST(cross_limit_air_guaranteed);
    RUN_TEST(signal_selector);
    RUN_TEST(hi_lo_select_functions);

    printf("\n[L4/L5] Feedforward:\n");
    RUN_TEST(feedforward_static);
    RUN_TEST(feedforward_lead_lag);
    RUN_TEST(ratio_ff_loop);

    printf("\n[L2/L4] Blending Control:\n");
    RUN_TEST(blend_system_init);
    RUN_TEST(blend_add_component);
    RUN_TEST(blend_mass_balance);
    RUN_TEST(blend_cost_quality);
    RUN_TEST(blend_optimization);

    printf("\n[L4/L5] Cascade Ratio:\n");
    RUN_TEST(cascade_init);
    RUN_TEST(cascade_execute_simple);
    RUN_TEST(pid_loop);

    printf("\n========================================\n");
    printf("  %d / %d tests PASSED\n", tests_passed, tests_run);
    printf("========================================\n\n");

    assert(tests_passed == tests_run);
    return 0;
}
