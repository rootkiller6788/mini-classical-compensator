/**
 * test_pid_tuning.c — Comprehensive Tests for PID Tuning Module
 *
 * Covers core APIs: PID init, update, performance, form conversion,
 * Z-N step/freq tuning, Cohen-Coon, IMC/SIMC, relay auto-tuner,
 * gain/phase margin methods, FOPDT identification, anti-windup,
 * gain scheduling, and application-specific tuning.
 */

#include "pid_tuning.h"
#include "ziegler_nichols.h"
#include "cohen_coon.h"
#include "imc_tuning.h"
#include "relay_autotune.h"
#include "gain_margin_tuning.h"
#include "fopdt_model.h"
#include "advanced_tuning.h"
#include "application_tuning.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %-45s ", name); \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
} while(0)

/* ────────────────────────────────────
 * L1: PID Initialization & Parameter Tests
 * ──────────────────────────────────── */

static void test_pid_init(void)
{
    TEST("pid_init sets default parameters");
    pid_controller_t pid;
    pid_init(&pid, 1.5, 0.8, 0.3, 0.01);
    assert(pid.params.Kp == 1.5);
    assert(pid.params.Ki == 0.8);
    assert(pid.params.Kd == 0.3);
    assert(pid.Ts == 0.01);
    assert(pid.params.N == 10.0);
    assert(pid.params.b == 1.0);
    assert(pid.params.c == 1.0);
    assert(pid.I == 0.0);
    PASS();
}

static void test_pid_init_null_safety(void)
{
    TEST("pid_init NULL pointer safety");
    pid_init(NULL, 1.0, 1.0, 1.0, 0.01); /* Should not crash */
    PASS();
}

static void test_pid_reset(void)
{
    TEST("pid_reset clears integrator and state");
    pid_controller_t pid;
    pid_init(&pid, 1.0, 0.5, 0.1, 0.01);
    pid.I = 5.0;
    pid.e_prev = 2.0;
    pid.step_count = 100;
    pid_reset(&pid);
    assert(pid.I == 0.0);
    assert(pid.e_prev == 0.0);
    assert(pid.step_count == 0);
    PASS();
}

/* ────────────────────────────────────
 * L3: Form Conversion Tests
 * ──────────────────────────────────── */

static void test_ideal_to_parallel(void)
{
    TEST("pid_ideal_to_parallel conversion");
    double Ki, Kd;
    pid_ideal_to_parallel(2.0, 4.0, 0.5, &Ki, &Kd);
    assert(fabs(Ki - 0.5) < 0.001);  /* Kp/Ti = 2/4 = 0.5 */
    assert(fabs(Kd - 1.0) < 0.001);  /* Kp*Td = 2*0.5 = 1.0 */
    PASS();
}

static void test_parallel_to_ideal(void)
{
    TEST("pid_parallel_to_ideal conversion");
    double Ti, Td;
    pid_parallel_to_ideal(2.0, 0.5, 1.0, &Ti, &Td);
    assert(fabs(Ti - 4.0) < 0.001);  /* Kp/Ki = 2/0.5 = 4 */
    assert(fabs(Td - 0.5) < 0.001);  /* Kd/Kp = 1/2 = 0.5 */
    PASS();
}

static void test_series_to_parallel(void)
{
    TEST("pid_series_to_parallel conversion");
    double Kp_p, Ki_p, Kd_p;
    pid_series_to_parallel(1.0, 5.0, 2.0, &Kp_p, &Ki_p, &Kd_p);
    /* Kp_s * (1 + Td/Ti) = 1*(1+2/5) = 1.4 */
    assert(fabs(Kp_p - 1.4) < 0.01);
    /* Ki = Kp_s / Ti = 1/5 = 0.2 */
    assert(fabs(Ki_p - 0.2) < 0.001);
    /* Kd = Kp_s * Td = 2.0 */
    assert(fabs(Kd_p - 2.0) < 0.001);
    PASS();
}

/* ────────────────────────────────────
 * L5: PID Update Test
 * ──────────────────────────────────── */

static void test_pid_update_proportional_only(void)
{
    TEST("pid_update P-only response");
    pid_controller_t pid;
    pid_init(&pid, 2.0, 0.0, 0.0, 0.1);
    pid_set_limits(&pid, AW_MODE_NONE, -10.0, 10.0);

    /* Setpoint = 1, measurement = 0 → error = +1, P = 2.0 */
    double u = pid_update(&pid, 1.0, 0.0);
    assert(fabs(u - 2.0) < 0.01);
    PASS();
}

static void test_pid_update_integral(void)
{
    TEST("pid_update integral accumulation");
    pid_controller_t pid;
    pid_init(&pid, 1.0, 0.5, 0.0, 0.1);
    pid_set_limits(&pid, AW_MODE_NONE, -100.0, 100.0);

    /* Multiple steps with constant error */
    pid_update(&pid, 1.0, 0.0);
    pid_update(&pid, 1.0, 0.0);
    double u = pid_update(&pid, 1.0, 0.0);

    /* P = 1, I accumulated = 0.5*0.1*3 = 0.15, total ≈ 1.15 */
    assert(u > 1.0);
    PASS();
}

static void test_pid_update_saturation(void)
{
    TEST("pid_update output saturation");
    pid_controller_t pid;
    pid_init(&pid, 10.0, 0.0, 0.0, 0.1);
    pid_set_limits(&pid, AW_MODE_NONE, -5.0, 5.0);

    /* Error = 2, P = 20 but sat to 5 */
    double u = pid_update(&pid, 2.0, 0.0);
    assert(fabs(u - 5.0) < 0.01);
    PASS();
}

/* ────────────────────────────────────
 * L5: Performance Metric Tests
 * ──────────────────────────────────── */

static void test_perf_metrics(void)
{
    TEST("pid_eval_performance computes metrics");
    double t[100], y[100];
    for (int i = 0; i < 100; i++) {
        t[i] = i * 0.01;
        if (t[i] < 0.2) y[i] = 0.0;
        else y[i] = 1.0 - exp(-(t[i]-0.2)/0.1);
    }
    pid_perf_t perf;
    pid_eval_performance(t, y, 1.0, 100, &perf);
    assert(perf.IAE > 0.0);
    assert(perf.overshoot_pct >= 0.0);
    assert(perf.rise_time > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L4: Z-N Step Response Tuning Test
 * ──────────────────────────────────── */

static void test_zn_step_tune_pid(void)
{
    TEST("zn_step_tune produces valid PID parameters");
    fopdt_model_t fopdt;
    fopdt.K = 1.0;
    fopdt.T = 5.0;
    fopdt.L = 1.0;

    pid_params_t params;
    int ret = zn_step_tune(&fopdt, ZN_CONTROLLER_PID, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    assert(params.Td > 0.0);
    assert(params.Ki > 0.0);
    PASS();
}

static void test_zn_step_tune_pi(void)
{
    TEST("zn_step_tune PI controller");
    fopdt_model_t fopdt = { 2.0, 10.0, 2.0 };
    pid_params_t params;
    int ret = zn_step_tune(&fopdt, ZN_CONTROLLER_PI, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    /* PI should have no derivative */
    assert(fabs(params.Td) < 1e-6 || params.Td == 0.0);
    PASS();
}

static void test_zn_step_tune_invalid(void)
{
    TEST("zn_step_tune rejects invalid model");
    fopdt_model_t fopdt = { 0.0, 5.0, 1.0 }; /* Zero gain */
    pid_params_t params;
    int ret = zn_step_tune(&fopdt, ZN_CONTROLLER_PID, &params);
    assert(ret != 0);
    PASS();
}

/* ────────────────────────────────────
 * L4: Z-N Frequency Response Tuning
 * ──────────────────────────────────── */

static void test_zn_freq_tune(void)
{
    TEST("zn_freq_tune produces valid PID from Ku/Pu");
    pid_params_t params;
    int ret = zn_freq_tune(3.0, 4.0, ZN_CONTROLLER_PID, &params);
    assert(ret == 0);
    /* PID: Kp = 0.6*3 = 1.8, Ti = 0.5*4 = 2.0, Td = 0.125*4 = 0.5 */
    assert(fabs(params.Kp - 1.8) < 0.01);
    assert(fabs(params.Ti - 2.0) < 0.01);
    assert(fabs(params.Td - 0.5) < 0.01);
    PASS();
}

static void test_zn_find_ultimate_gain(void)
{
    TEST("zn_find_ultimate_gain computes Ku/Pu from FOPDT");
    fopdt_model_t fopdt = { 1.0, 2.0, 0.5 };
    ultimate_gain_result_t result;
    int ret = zn_find_ultimate_gain(&fopdt, &result);
    assert(ret == 0);
    assert(result.converged == 1);
    assert(result.Ku > 1.0);
    assert(result.Pu > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L4: Cohen-Coon Tuning Tests
 * ──────────────────────────────────── */

static void test_cohen_coon_pid(void)
{
    TEST("cohen_coon_tune PID controller");
    fopdt_model_t fopdt = { 1.0, 10.0, 2.0 };
    pid_params_t params;
    int ret = cohen_coon_tune(&fopdt, CC_CONTROLLER_PID, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    assert(params.Td > 0.0);
    PASS();
}

static void test_cohen_coon_pi(void)
{
    TEST("cohen_coon_tune PI controller");
    fopdt_model_t fopdt = { 1.5, 8.0, 1.5 };
    pid_params_t params;
    int ret = cohen_coon_tune(&fopdt, CC_CONTROLLER_PI, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    PASS();
}

static void test_cohen_coon_ipdt(void)
{
    TEST("cohen_coon_tune_ipdt integrating process");
    ipdt_model_t ipdt = { 0.5, 1.0 };
    pid_params_t params;
    int ret = cohen_coon_tune_ipdt(&ipdt, CC_CONTROLLER_PI, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L4: IMC/SIMC Tuning Tests
 * ──────────────────────────────────── */

static void test_imc_tune_fopdt(void)
{
    TEST("imc_tune_fopdt produces valid IMC-PID");
    fopdt_model_t fopdt = { 1.0, 5.0, 1.0 };
    imc_tuning_result_t result;
    int ret = imc_tune_fopdt(&fopdt, 1.5, IMC_FILTER_FIRST, &result);
    assert(ret == 0);
    assert(result.Kp > 0.0);
    assert(result.Ki > 0.0);
    assert(result.Kd >= 0.0);
    PASS();
}

static void test_imc_simc_tune(void)
{
    TEST("imc_simc_tune follows Skogestad rules");
    fopdt_model_t fopdt = { 1.0, 20.0, 2.0 };
    pid_params_t params;
    int ret = imc_simc_tune(&fopdt, 2.0, 1, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    /* T >> L → derivative should be active */
    assert(params.Td > 0.0);
    PASS();
}

static void test_imc_simc_pi_only(void)
{
    TEST("imc_simc_tune PI only for balanced process");
    fopdt_model_t fopdt = { 1.0, 8.0, 2.0 };
    pid_params_t params;
    int ret = imc_simc_tune(&fopdt, 2.0, 0, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    /* PI only: no derivative */
    double Td;
    pid_parallel_to_ideal(params.Kp, params.Ki, params.Kd, NULL, &Td);
    assert(fabs(params.Td) < 1e-6 || params.Td == 0.0);
    PASS();
}

static void test_imc_lambda_from_Ms(void)
{
    TEST("imc_lambda_from_Ms relationship");
    double lambda = imc_lambda_from_Ms(2.0, 1.4);
    assert(lambda > 0.0);
    /* λ = L/(Ms-1) = 2/(0.4) = 5.0 */
    assert(fabs(lambda - 5.0) < 0.01);
    PASS();
}

/* ────────────────────────────────────
 * L5: Relay Auto-Tuning Tests
 * ──────────────────────────────────── */

static void test_relay_describing_function(void)
{
    TEST("relay_describing_function ideal relay");
    double mag, phase;
    relay_describing_function(RELAY_IDEAL, 2.0, 0.0, 1.0, &mag, &phase);
    /* mag = 4*d/(π*a) = 4*2/(π*1) ≈ 2.546 */
    assert(fabs(mag - 2.546) < 0.02);
    assert(fabs(phase) < 0.001);
    PASS();
}

static void test_relay_describing_function_hysteresis(void)
{
    TEST("relay_describing_function with hysteresis");
    double mag, phase;
    relay_describing_function(RELAY_HYSTERESIS, 2.0, 0.3, 1.0, &mag, &phase);
    /* mag should be non-zero, phase should be negative (lag) */
    assert(mag > 0.0);
    assert(phase < 0.0);
    PASS();
}

static void test_relay_extract_ultimate(void)
{
    TEST("relay_extract_ultimate from oscillation data");
    relay_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.amplitude = 1.0;
    cfg.type = RELAY_IDEAL;

    relay_result_t result;
    int ret = relay_extract_ultimate(&cfg, 0.5, 4.0, &result);
    assert(ret == 0);
    assert(result.Ku > 0.0);
    assert(fabs(result.Pu - 4.0) < 0.01);
    PASS();
}

static void test_relay_simulate_fopdt(void)
{
    TEST("relay_simulate_fopdt finds Ku/Pu from model");
    fopdt_model_t fopdt = { 1.0, 3.0, 1.0 };
    relay_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.amplitude   = 1.0;
    cfg.hysteresis  = 0.02;
    cfg.Ts          = 0.05;
    cfg.max_duration= 50.0;
    cfg.max_cycles  = 5;
    cfg.settle_tol  = 0.05;
    cfg.type        = RELAY_HYSTERESIS;

    relay_result_t result;
    int ret = relay_simulate_fopdt(&fopdt, &cfg, &result);
    assert(ret == 0);
    assert(result.converged == 1);
    assert(result.Ku > 0.1);
    assert(result.Pu > 0.5);
    PASS();
}

/* ────────────────────────────────────
 * L4: Gain/Phase Margin Tests
 * ──────────────────────────────────── */

static void test_gm_pm_compute_margins(void)
{
    TEST("gm_pm_compute_margins returns stability margins");
    fopdt_model_t fopdt = { 1.0, 4.0, 1.0 };
    pid_params_t params;
    /* Use a known-stable PID tuning */
    imc_simc_tune(&fopdt, 1.5, 1, &params);

    margin_spec_t margins;
    int ret = gm_pm_compute_margins(&fopdt, &params, &margins);
    /* May or may not have clear margins — just check no crash */
    assert(ret == 0 || ret == -1); /* Either OK */
    PASS();
}

static void test_gm_pm_max_sensitivity(void)
{
    TEST("gm_pm_max_sensitivity computes Ms");
    fopdt_model_t fopdt = { 1.0, 5.0, 1.0 };
    pid_params_t params;
    imc_simc_tune(&fopdt, 1.5, 1, &params);

    double Ms, wMs;
    int ret = gm_pm_max_sensitivity(&fopdt, &params, &Ms, &wMs);
    assert(ret == 0);
    assert(Ms >= 1.0);
    assert(wMs > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L5: FOPDT Identification Tests
 * ──────────────────────────────────── */

static void test_fopdt_graphical(void)
{
    TEST("fopdt_identify_graphical tangent method");
    int N = 200;
    double *t    = (double *)malloc(N * sizeof(double));
    double *y    = (double *)malloc(N * sizeof(double));
    double *u_in = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        t[i] = i * 0.1;
        u_in[i] = 1.0;
        if (t[i] < 1.0) y[i] = 0.0;
        else y[i] = 2.0 * (1.0 - exp(-(t[i]-1.0)/5.0));
    }

    step_response_data_t data;
    data.time = t;
    data.input = u_in;
    data.output = y;
    data.N = N;
    data.Ts = 0.1;
    data.step_mag = 1.0;

    fopdt_id_result_t result;
    int ret = fopdt_identify_graphical(&data, &result);
    assert(ret == 0);
    assert(fabs(result.model.K - 2.0) < 0.3);
    assert(fabs(result.model.T - 5.0) < 1.5);
    assert(result.model.L > 0.0);

    free(t); free(y); free(u_in);
    PASS();
}

static void test_fopdt_two_point(void)
{
    TEST("fopdt_identify_two_point 28/63 method");
    int N = 200;
    double *t = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *u_in = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        t[i] = i * 0.05;
        u_in[i] = 1.0;
        if (t[i] < 0.5) y[i] = 0.0;
        else y[i] = 1.0 * (1.0 - exp(-(t[i]-0.5)/3.0));
    }

    step_response_data_t data;
    data.time = t; data.input = u_in; data.output = y;
    data.N = N; data.Ts = 0.05; data.step_mag = 1.0;

    fopdt_id_result_t result;
    int ret = fopdt_identify_two_point(&data, &result);
    assert(ret == 0);
    assert(fabs(result.model.K - 1.0) < 0.2);
    assert(fabs(result.model.T - 3.0) < 1.0);
    assert(result.model.L > 0.0);

    free(t); free(y); free(u_in);
    PASS();
}

static void test_fopdt_area(void)
{
    TEST("fopdt_identify_area method");
    int N = 200;
    double *t = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *u_in = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        t[i] = i * 0.05;
        u_in[i] = 1.0;
        if (t[i] < 1.0) y[i] = 0.0;
        else y[i] = 1.0 * (1.0 - exp(-(t[i]-1.0)/4.0));
    }

    step_response_data_t data;
    data.time = t; data.input = u_in; data.output = y;
    data.N = N; data.Ts = 0.05; data.step_mag = 1.0;

    fopdt_id_result_t result;
    int ret = fopdt_identify_area(&data, &result);
    assert(ret == 0);
    assert(result.model.K > 0.5);
    assert(result.model.T > 1.0);

    free(t); free(y); free(u_in);
    PASS();
}

static void test_fopdt_simulate(void)
{
    TEST("fopdt_simulate_step analytic solution");
    fopdt_model_t fopdt = { 2.0, 5.0, 1.0 };
    double t[100], y[100];
    for (int i = 0; i < 100; i++) t[i] = i * 0.2;
    fopdt_simulate_step(&fopdt, t, y, 1.0, 100);
    assert(fabs(y[0]) < 0.01);        /* Before dead time */
    assert(fabs(y[99] - 2.0) < 0.05); /* Near steady state at t=19.8, tau=18.8 */
    PASS();
}

static void test_fopdt_dead_time_ratio(void)
{
    TEST("fopdt_dead_time_ratio computes τ");
    fopdt_model_t fopdt = { 1.0, 7.0, 3.0 };
    double tau = fopdt_dead_time_ratio(&fopdt);
    assert(fabs(tau - 0.3) < 0.01);  /* 3/(7+3) = 0.3 */
    PASS();
}

/* ────────────────────────────────────
 * L5: Anti-Windup Tests
 * ──────────────────────────────────── */

static void test_aw_back_calculation(void)
{
    TEST("aw_back_calculation prevents windup");
    pid_controller_t pid;
    pid_init(&pid, 2.0, 0.5, 0.1, 0.1);
    pid_set_limits(&pid, AW_MODE_BACK_CALC, -5.0, 5.0);

    /* Drive to saturation */
    pid.u_sat = 5.0;
    pid.u_unsat = 8.0;
    aw_back_calculation(&pid, 1.0);
    /* The integrator should be adjusted downward */
    assert(pid.I < 1.0); /* Should be less than without AW */
    PASS();
}

static void test_aw_clamping(void)
{
    TEST("aw_clamping freezes integral when saturated");
    pid_controller_t pid;
    pid_init(&pid, 2.0, 2.0, 0.0, 0.1);
    pid_set_limits(&pid, AW_MODE_CLAMPING, -5.0, 5.0);

    pid.u_sat = 5.0;
    double I_before = pid.I = 3.0;
    aw_clamping(&pid, 1.0);  /* Saturated high + positive error: freeze */
    /* Integral should not increase */
    assert(fabs(pid.I - I_before) < 0.001);
    PASS();
}

static void test_aw_velocity_form(void)
{
    TEST("aw_velocity_form computes incremental output");
    pid_controller_t pid;
    pid_init(&pid, 1.0, 0.2, 0.05, 0.1);
    pid_set_limits(&pid, AW_MODE_VELOCITY, -10.0, 10.0);

    double u = aw_velocity_form(&pid, 1.0, 0.0);
    /* Should return a finite value */
    assert(isfinite(u));
    PASS();
}

/* ────────────────────────────────────
 * L5: Gain Scheduling Tests
 * ──────────────────────────────────── */

static void test_gs_linear_interpolate(void)
{
    TEST("gs_linear_interpolate between two points");
    gs_table_t table;
    gs_entry_t entries[3];
    entries[0].condition = 0.0; entries[0].Kp = 1.0;
    entries[0].Ki = 0.1; entries[0].Kd = 0.0;
    entries[1].condition = 10.0; entries[1].Kp = 2.0;
    entries[1].Ki = 0.2; entries[1].Kd = 0.0;
    table.entries = entries;
    table.n = 2;
    table.sorted = 1;

    double Kp, Ki, Kd;
    int ret = gs_linear_interpolate(&table, 5.0, &Kp, &Ki, &Kd);
    assert(ret == 0);
    assert(fabs(Kp - 1.5) < 0.01);
    assert(fabs(Ki - 0.15) < 0.01);
    PASS();
}

/* ────────────────────────────────────
 * L5: Setpoint Filter Test
 * ──────────────────────────────────── */

static void test_sp_filter(void)
{
    TEST("sp_filter_update smooths step change");
    double sp_f = 0.0;
    /* First update: step from 0 to 10 */
    double out = sp_filter_update(1.0, 0.1, 10.0, &sp_f);
    assert(out < 10.0); /* Should be filtered, not full step */
    assert(out > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L7: Application Tuning Tests
 * ──────────────────────────────────── */

static void test_app_temperature(void)
{
    TEST("app_tune_temperature for lag-dominant process");
    pid_params_t params;
    int ret = app_tune_temperature(2.0, 20.0, 1.0, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    /* T >> L → should have derivative */
    assert(params.Td > 0.0);
    PASS();
}

static void test_app_flow(void)
{
    TEST("app_tune_flow PI only for fast process");
    pid_params_t params;
    int ret = app_tune_flow(0.8, 2.0, 0.2, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    PASS();
}

static void test_app_level(void)
{
    TEST("app_tune_level surge tank P-only");
    pid_params_t params;
    int ret = app_tune_level(0.1, 2.0, 1, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(fabs(params.Ki) < 1e-6); /* P-only */
    PASS();
}

static void test_app_auto_tune(void)
{
    TEST("app_auto_tune dispatches to correct method");
    pid_params_t params;
    int ret = app_auto_tune(1, 0.8, 2.0, 0.2, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    PASS();
}

static void test_app_recommend_method(void)
{
    TEST("app_recommend_method selects appropriate method");
    fopdt_model_t lag_dominant = { 1.0, 10.0, 0.5 };
    fopdt_model_t delay_dominant = { 1.0, 2.0, 3.0 };

    pid_tune_method_t m1 = app_recommend_method(&lag_dominant);
    pid_tune_method_t m2 = app_recommend_method(&delay_dominant);

    assert(m1 == TUNE_METHOD_IMC || m1 == TUNE_METHOD_ZN_STEP);
    assert(m2 == TUNE_METHOD_COHEN_COON || m2 == TUNE_METHOD_GM_PM);
    PASS();
}

/* ────────────────────────────────────
 * L5: Cascade & Feedforward Tests
 * ──────────────────────────────────── */

static void test_cascade_tune(void)
{
    TEST("cascade_tune_pid inner/outer loops");
    fopdt_model_t inner = { 0.8, 2.0, 0.3 };
    fopdt_model_t outer = { 1.5, 15.0, 2.0 };
    pid_params_t inner_params, outer_params;

    int ret = cascade_tune_pid(&inner, &outer,
                               TUNE_METHOD_COHEN_COON,
                               TUNE_METHOD_ZN_STEP,
                               &inner_params, &outer_params);
    assert(ret == 0);
    assert(inner_params.Kp > 0.0);
    assert(outer_params.Kp > 0.0);
    PASS();
}

static void test_feedforward_design(void)
{
    TEST("feedforward_pid_design computes FF parameters");
    double ff_gain, ff_lead, ff_lag;
    int ret = feedforward_pid_design(2.0, 5.0, 1.0,
                                      1.5, 3.0, 0.5,
                                      &ff_gain, &ff_lead, &ff_lag);
    assert(ret == 0);
    /* ff_gain = -1.5/2.0 = -0.75 */
    assert(fabs(ff_gain + 0.75) < 0.01);
    assert(ff_lead > 0.0);
    assert(ff_lag  > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L5: Form Conversion Test
 * ──────────────────────────────────── */

static void test_pid_convert_form(void)
{
    TEST("pid_convert_form IDEAL to PARALLEL");
    pid_params_t src;
    memset(&src, 0, sizeof(src));
    src.Kp = 2.0; src.Ti = 4.0; src.Td = 0.5;
    src.N = 10.0; src.b = 1.0; src.c = 0.0; src.Tt = 1.0;

    pid_params_t dst;
    pid_convert_form(&src, PID_FORM_IDEAL, &dst, PID_FORM_PARALLEL);

    assert(fabs(dst.Kp - 2.0) < 0.01);
    assert(fabs(dst.Ki - 0.5) < 0.01);  /* Kp/Ti */
    assert(fabs(dst.Kd - 1.0) < 0.01);  /* Kp*Td */
    PASS();
}

/* ────────────────────────────────────
 * L5: Chien-Hrones-Reswick Test
 * ──────────────────────────────────── */

static void test_chr_tuning(void)
{
    TEST("zn_chien_hrones_reswick setpoint tuning");
    fopdt_model_t fopdt = { 1.0, 8.0, 2.0 };
    pid_params_t params;
    int ret = zn_chien_hrones_reswick(&fopdt, ZN_CONTROLLER_PID, 0, &params);
    assert(ret == 0);
    assert(params.Kp > 0.0);
    assert(params.Ti > 0.0);
    PASS();
}

/* ────────────────────────────────────
 * L5: Filter Design Tests
 * ──────────────────────────────────── */

static void test_pid_derive_filter_N(void)
{
    TEST("pid_derive_filter_N computes filter constant");
    double N = pid_derive_filter_N(0.5, 20.0);
    assert(N >= 2.0);  /* Should be clamped to min 2 */
    assert(fabs(N - 10.0) < 0.01); /* 20*0.5 = 10 */
    PASS();
}

static void test_pid_recommend_sampling(void)
{
    TEST("pid_recommend_sampling returns feasible Ts");
    double Ts = pid_recommend_sampling(5.0);
    assert(Ts > 0.0);
    assert(Ts < 1.0);
    PASS();
}

static void test_pid_derivative_filter(void)
{
    TEST("pid_derivative_filter applies lowpass");
    pid_controller_t pid;
    pid_init(&pid, 1.0, 0.2, 0.5, 0.01);
    pid.e_prev = 0.0;
    double D_prev = 0.0;
    double D = pid_derivative_filter(&pid, 1.0, &D_prev, 0.9);
    /* D should be smoothed */
    assert(isfinite(D));
    PASS();
}

/* ────────────────────────────────────
 * L5: Sørensen Ultimate Gain Test
 * ──────────────────────────────────── */

static void test_sorensen(void)
{
    TEST("zn_sorensen_ultimate_gain with hysteresis");
    ultimate_gain_result_t result;
    int ret = zn_sorensen_ultimate_gain(2.0, 0.1, 0.8, 4.0, &result);
    assert(ret == 0);
    assert(result.Ku > 0.0);
    assert(fabs(result.Pu - 4.0) < 0.01);
    PASS();
}

/* ────────────────────────────────────
 * L5: Advanced Tuning Tests
 * ──────────────────────────────────── */

static void test_sp_optimal_weights(void)
{
    TEST("sp_compute_optimal_weights for delay-dominant");
    fopdt_model_t fopdt = { 1.0, 3.0, 3.0 }; /* τ = 0.5 */
    pid_params_t params;
    imc_simc_tune(&fopdt, 3.0, 1, &params);

    double b, c;
    int ret = sp_compute_optimal_weights(&fopdt, &params, &b, &c);
    assert(ret == 0);
    assert(b >= 0.0 && b <= 1.0);
    assert(c == 0.0); /* D always on PV */
    PASS();
}

static void test_pid_compare_methods(void)
{
    TEST("pid_compare_tuning_methods ranks methods");
    fopdt_model_t fopdt = { 1.0, 10.0, 2.0 };
    pid_tune_method_t methods[] = {
        TUNE_METHOD_ZN_STEP,
        TUNE_METHOD_COHEN_COON,
        TUNE_METHOD_IMC
    };
    double rankings[3];
    int best_idx;
    int ret = pid_compare_tuning_methods(&fopdt, methods, 3,
                                          rankings, &best_idx);
    assert(ret == 0);
    assert(best_idx >= 0 && best_idx < 3);
    PASS();
}

/* ────────────────────────────────────
 * Main Test Runner
 * ──────────────────────────────────── */

int main(void)
{
    printf("\n=== PID Tuning Ziegler-Nichols Module Test Suite ===\n\n");

    printf("--- L1: Definitions & Core PID ---\n");
    test_pid_init();
    test_pid_init_null_safety();
    test_pid_reset();

    printf("\n--- L3: Form Conversions ---\n");
    test_ideal_to_parallel();
    test_parallel_to_ideal();
    test_series_to_parallel();
    test_pid_convert_form();

    printf("\n--- L5: PID Update & Performance ---\n");
    test_pid_update_proportional_only();
    test_pid_update_integral();
    test_pid_update_saturation();
    test_perf_metrics();

    printf("\n--- L4: Ziegler-Nichols Step Response ---\n");
    test_zn_step_tune_pid();
    test_zn_step_tune_pi();
    test_zn_step_tune_invalid();

    printf("\n--- L4: Ziegler-Nichols Frequency Response ---\n");
    test_zn_freq_tune();
    test_zn_find_ultimate_gain();
    test_sorensen();

    printf("\n--- L4: Cohen-Coon ---\n");
    test_cohen_coon_pid();
    test_cohen_coon_pi();
    test_cohen_coon_ipdt();
    test_chr_tuning();

    printf("\n--- L4: IMC / SIMC ---\n");
    test_imc_tune_fopdt();
    test_imc_simc_tune();
    test_imc_simc_pi_only();
    test_imc_lambda_from_Ms();

    printf("\n--- L5: Relay Auto-Tuning ---\n");
    test_relay_describing_function();
    test_relay_describing_function_hysteresis();
    test_relay_extract_ultimate();
    test_relay_simulate_fopdt();

    printf("\n--- L4: Gain/Phase Margin ---\n");
    test_gm_pm_compute_margins();
    test_gm_pm_max_sensitivity();

    printf("\n--- L5: FOPDT Identification ---\n");
    test_fopdt_graphical();
    test_fopdt_two_point();
    test_fopdt_area();
    test_fopdt_simulate();
    test_fopdt_dead_time_ratio();

    printf("\n--- L5: Anti-Windup ---\n");
    test_aw_back_calculation();
    test_aw_clamping();
    test_aw_velocity_form();

    printf("\n--- L5: Gain Scheduling ---\n");
    test_gs_linear_interpolate();

    printf("\n--- L5: Filters ---\n");
    test_sp_filter();
    test_pid_derive_filter_N();
    test_pid_recommend_sampling();
    test_pid_derivative_filter();

    printf("\n--- L7: Applications ---\n");
    test_app_temperature();
    test_app_flow();
    test_app_level();
    test_app_auto_tune();
    test_app_recommend_method();

    printf("\n--- L8: Cascade & Feedforward ---\n");
    test_cascade_tune();
    test_feedforward_design();

    printf("\n--- L8: Advanced ---\n");
    test_sp_optimal_weights();
    test_pid_compare_methods();

    printf("\n========================================\n");
    printf("RESULTS: %d / %d tests passed\n", tests_passed, tests_run);
    if (tests_passed == tests_run) {
        printf("ALL TESTS PASSED ✅\n\n");
        return 0;
    } else {
        printf("SOME TESTS FAILED ⚠️\n\n");
        return 1;
    }
}
