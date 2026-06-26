/**
 * @file test_cascade.c
 * @brief Comprehensive test suite for cascade control module
 *
 * Tests cover:
 *   - Polynomial operations (create, eval, mul, add)
 *   - Transfer function creation and frequency response
 *   - PID controller conversion and frequency response
 *   - Loop closure (inner, equivalent plant, overall)
 *   - Routh-Hurwitz stability criterion
 *   - Step response computation
 *   - Cascade design (sequential, tuning rules)
 *   - Digital implementation (velocity/position PID, discretization)
 *   - Anti-windup and bumpless transfer
 *   - Performance metrics
 *   - Canonical system models
 *   - Internal stability verification
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)
#define CHECK_CLOSE(a, b, tol, msg) \
    do { if (fabs((a)-(b)) < (tol)) PASS(); else { \
        printf("FAIL: %s (%.6f vs %.6f)\n", msg, a, b); tests_failed++; \
    } } while(0)

/* ==========================================================================
 * L3: Polynomial Operations Tests
 * ========================================================================== */

static void test_polynomial_create(void) {
    TEST("poly_create");
    double c[] = {1.0, 2.0, 3.0};
    CascadePoly p = cascade_poly_create(c, 2);
    CHECK(p.degree == 2 && p.coeff && p.coeff[0] == 1.0 &&
          p.coeff[1] == 2.0 && p.coeff[2] == 3.0, "coefficient mismatch");
    cascade_poly_free(&p);
}

static void test_polynomial_eval(void) {
    TEST("poly_eval_complex");
    /* p(s) = 1 + 2s + s^2, evaluate at s = 1 + j0: p(1) = 1+2+1 = 4 */
    double c[] = {1.0, 2.0, 1.0};
    CascadePoly p = cascade_poly_create(c, 2);
    double re, im;
    cascade_poly_eval_complex(&p, 1.0, 0.0, &re, &im);
    CHECK_CLOSE(re, 4.0, 1e-10, "p(1) should be 4");
    CHECK_CLOSE(im, 0.0, 1e-10, "imag should be 0");

    /* Evaluate at s = j: p(j) = 1 + 2j - 1 = j*2 */
    cascade_poly_eval_complex(&p, 0.0, 1.0, &re, &im);
    CHECK_CLOSE(re, 0.0, 1e-10, "Re{p(j)} should be 0");
    CHECK_CLOSE(im, 2.0, 1e-10, "Im{p(j)} should be 2");
    cascade_poly_free(&p);
}

static void test_polynomial_mul(void) {
    TEST("poly_mul");
    /* (1 + s) * (1 + s) = 1 + 2s + s^2 */
    double c1[] = {1.0, 1.0}, c2[] = {1.0, 1.0};
    CascadePoly p1 = cascade_poly_create(c1, 1);
    CascadePoly p2 = cascade_poly_create(c2, 1);
    CascadePoly r = cascade_poly_mul(&p1, &p2);
    CHECK(r.degree == 2 && r.coeff[0] == 1.0 &&
          r.coeff[1] == 2.0 && r.coeff[2] == 1.0, "multiplication error");
    cascade_poly_free(&p1); cascade_poly_free(&p2);
    cascade_poly_free(&r);
}

static void test_polynomial_add(void) {
    TEST("poly_add");
    double c1[] = {1.0, 2.0}, c2[] = {3.0, 4.0, 1.0};
    CascadePoly p1 = cascade_poly_create(c1, 1);
    CascadePoly p2 = cascade_poly_create(c2, 2);
    CascadePoly r = cascade_poly_add(&p1, &p2);
    /* (1+2s) + (3+4s+s^2) = 4 + 6s + s^2 */
    CHECK(r.degree == 2 && r.coeff[0] == 4.0 &&
          r.coeff[1] == 6.0 && r.coeff[2] == 1.0, "addition error");
    cascade_poly_free(&p1); cascade_poly_free(&p2);
    cascade_poly_free(&r);
}

/* ==========================================================================
 * L1: Transfer Function Tests
 * ========================================================================== */

static void test_tf_create_and_freq_resp(void) {
    TEST("tf_freq_response");
    /* G(s) = 1 / (s + 1): first-order lowpass
     * At w=1: |G| = 1/sqrt(2) ~ 0.707, phase = -45 deg */
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF tf = cascade_tf_create(num, 0, den, 1, 1.0);
    double mag, phase;
    cascade_tf_freq_response(&tf, 1.0, &mag, &phase);
    CHECK_CLOSE(mag, 0.70710678, 0.01, "|G(j1)| should be ~0.707");
    CHECK_CLOSE(phase, -M_PI/4.0, 0.01, "phase should be -45 deg");
    cascade_tf_free(&tf);
}

static void test_tf_close_loop(void) {
    TEST("tf_close_loop");
    /* L(s) = 1/s, unity feedback: G_cl = 1/(s+1) */
    double numL[] = {1.0}, denL[] = {0.0, 1.0};
    CascadeTF L = cascade_tf_create(numL, 0, denL, 1, 1.0);
    CascadeTF cl;
    cascade_tf_close_loop(&L, &cl);
    /* Expected: G_cl = 1/(s+1) */
    CHECK(cl.den.coeff[0] == 1.0 && cl.den.coeff[1] == 1.0, "denominator wrong");
    cascade_tf_free(&L); cascade_tf_free(&cl);
}

/* ==========================================================================
 * L1: PID Tests
 * ========================================================================== */

static void test_pid_to_tf(void) {
    TEST("pid_to_tf_PI");
    CascadePID pid;
    cascade_pid_init(2.0, 5.0, 0.0, 10.0, 0.0, &pid);
    CascadeTF tf;
    cascade_pid_to_tf(&pid, &tf);
    /* PI: C(s) = 2 + 5/s = (2s + 5)/s */
    CHECK(tf.num.degree == 1 && tf.den.degree == 1, "PI degree wrong");
    CHECK_CLOSE(tf.num.coeff[0], 5.0, 1e-10, "PI num coeff[0]");
    CHECK_CLOSE(tf.num.coeff[1], 2.0, 1e-10, "PI num coeff[1]");
    cascade_tf_free(&tf);
}

static void test_pid_freq_response(void) {
    TEST("pid_freq_response");
    CascadePID pid;
    cascade_pid_init(1.0, 2.0, 0.0, 10.0, 0.0, &pid);
    double mag, phase;
    /* At w=1: C(j1) = 1 + 2/j = 1 - j*2, |C| = sqrt(5) ~ 2.236 */
    cascade_pid_freq_response(&pid, 1.0, &mag, &phase);
    CHECK_CLOSE(mag, sqrt(5.0), 0.01, "|C(j1)| should be sqrt(5)");
}

/* ==========================================================================
 * L4: Stability Tests
 * ========================================================================== */

static void test_routh_hurwitz_stable(void) {
    TEST("routh_hurwitz_stable");
    /* s^2 + 2*s + 1 = (s+1)^2: stable */
    double den[] = {1.0, 2.0, 1.0};
    int n_rhp;
    cascade_routh_hurwitz(den, 2, &n_rhp);
    CHECK(n_rhp == 0, "stable poly should have 0 RHP roots");
}

static void test_routh_hurwitz_unstable(void) {
    TEST("routh_hurwitz_unstable");
    /* s^2 - 1 = (s-1)(s+1): 1 RHP root */
    double den[] = {-1.0, 0.0, 1.0};
    int n_rhp;
    cascade_routh_hurwitz(den, 2, &n_rhp);
    CHECK(n_rhp > 0, "unstable poly should have RHP roots");
}

static void test_is_stable(void) {
    TEST("cascade_is_stable");
    double num[] = {1.0}, den_stable[] = {1.0, 1.0}, den_unstable[] = {-1.0, 1.0};
    CascadeTF stable_tf = cascade_tf_create(num, 0, den_stable, 1, 1.0);
    CascadeTF unstable_tf = cascade_tf_create(num, 0, den_unstable, 1, 1.0);
    CHECK(cascade_is_stable(&stable_tf) == 1, "should be stable");
    CHECK(cascade_is_stable(&unstable_tf) == 0, "should be unstable");
    cascade_tf_free(&stable_tf); cascade_tf_free(&unstable_tf);
}

/* ==========================================================================
 * L2: Cascade Loop Formation Tests
 * ========================================================================== */

static void test_inner_closed_loop(void) {
    TEST("inner_closed_loop");
    /* Gi(s) = 1/(0.2s+1), Ci(s) = PI with Kp=1, Ki=5
     * Verify inner closed loop is stable */
    double num[] = {1.0}, den[] = {1.0, 0.2};
    CascadeTF Gi = cascade_tf_create(num, 0, den, 1, 1.0);
    CascadePID Ci;
    cascade_pid_init(1.0, 5.0, 0.0, 10.0, 0.0, &Ci);
    CascadeTF Gi_cl;
    int ret = cascade_inner_closed_loop(&Gi, &Ci, &Gi_cl);
    CHECK(ret == 0, "inner_closed_loop should succeed");
    int stable = cascade_is_stable(&Gi_cl);
    CHECK(stable == 1, "inner CL should be stable");
    cascade_tf_free(&Gi); cascade_tf_free(&Gi_cl);
}

static void test_equivalent_plant(void) {
    TEST("equivalent_plant");
    /* Gi = 1/(0.2s+1), Go = 1/(2s+1), with PI inner controller
     * Geq should be proper and stable */
    double num_fast[] = {1.0}, den_fast[] = {1.0, 0.2};
    double num_slow[] = {1.0}, den_slow[] = {1.0, 2.0};
    CascadeTF Gi = cascade_tf_create(num_fast, 0, den_fast, 1, 1.0);
    CascadeTF Go = cascade_tf_create(num_slow, 0, den_slow, 1, 1.0);
    CascadePID Ci;
    cascade_pid_init(1.0, 5.0, 0.0, 10.0, 0.0, &Ci);
    CascadeTF Geq;
    int ret = cascade_form_equivalent_plant(&Gi, &Ci, &Go, &Geq);
    CHECK(ret == 0, "form_equivalent_plant should succeed");
    CHECK(Geq.den.degree > 0, "Geq should have non-zero denominator");
    cascade_tf_free(&Gi); cascade_tf_free(&Go);
    cascade_tf_free(&Geq);
}

/* ==========================================================================
 * L5: Cascade Design Tests
 * ========================================================================== */

static void test_sequential_design(void) {
    TEST("sequential_design");
    /* Design cascade for Gi=1/(0.2s+1), Go=1/(2s+1) */
    double nf[] = {1.0}, df[] = {1.0, 0.2};
    double ns[] = {1.0}, ds[] = {1.0, 2.0};
    CascadeTF inner = cascade_tf_create(nf, 0, df, 1, 1.0);
    CascadeTF outer = cascade_tf_create(ns, 0, ds, 1, 1.0);
    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    CascadeSystem sys;
    int ret = cascade_design_sequential(&inner, &outer, &spec, &sys);
    CHECK(ret == 0, "sequential design should succeed");
    CHECK(sys.inner.inner_is_stable == 1, "inner should be stable");
    CHECK(sys.outer.outer_is_stable == 1, "outer should be stable");
    CHECK(sys.bandwidth_ratio >= 1.0, "BW ratio should be >= 1");
    cascade_tf_free(&inner); cascade_tf_free(&outer);
    cascade_system_free(&sys);
}

static void test_skogestad_tuning(void) {
    TEST("skogestad_tuning");
    CascadePID pid;
    int ret = cascade_tune_inner_skogestad(1.0, 2.0, 0.1, 0.5, &pid);
    CHECK(ret == 0, "Skogestad tuning should succeed");
    CHECK(pid.Kp > 0.0, "Kp should be positive");
    CHECK(pid.Ki > 0.0, "Ki should be positive");
}

static void test_direct_synthesis(void) {
    TEST("direct_synthesis");
    CascadePID pid;
    int ret = cascade_tune_outer_direct_synthesis(1.0, 5.0, 0.5, 2.0, &pid);
    CHECK(ret == 0, "DS tuning should succeed");
    CHECK(pid.Kp > 0.0, "Kp should be positive");
}

/* ==========================================================================
 * L5: Digital Implementation Tests
 * ========================================================================== */

static void test_pid_velocity_form(void) {
    TEST("pid_velocity_form");
    CascadePID pid;
    cascade_pid_init(1.0, 0.5, 0.0, 10.0, 0.01, &pid);
    cascade_pid_set_limits(-10.0, 10.0, &pid);
    double u = cascade_pid_velocity(1.0, 0.0, &pid);
    CHECK(fabs(u) < 10.1, "u should be within limits");
    /* Second call: error reduces, du should be smaller in magnitude */
    double u2 = cascade_pid_velocity(1.0, 0.5, &pid);
    CHECK(fabs(u2) <= 10.0, "u2 should be within limits");
}

static void test_pid_position_form(void) {
    TEST("pid_position_form");
    CascadePID pid;
    cascade_pid_init(1.0, 0.5, 0.0, 10.0, 0.01, &pid);
    cascade_pid_set_limits(-10.0, 10.0, &pid);
    cascade_pid_set_antiwindup(CASCADE_AW_BACK_CALC, 1.0, &pid);
    double u = cascade_pid_position(1.0, 0.0, &pid);
    CHECK(fabs(u) < 10.1, "u should be in bounds");
    /* Apply large error to test anti-windup */
    for (int i = 0; i < 100; i++) {
        cascade_pid_position(100.0, 0.0, &pid);
    }
    double u_sat = cascade_pid_position(100.0, 0.0, &pid);
    CHECK(fabs(u_sat) <= 10.0, "saturated output should be clamped");
}

static void test_tustin_discretize(void) {
    TEST("tustin_discretize");
    /* G(s) = 1/(s+1), Ts = 0.1 */
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF ct = cascade_tf_create(num, 0, den, 1, 1.0);
    CascadeDTF dt;
    int ret = cascade_tustin_discretize(&ct, 0.1, &dt);
    CHECK(ret == 0, "Tustin discretize should succeed");
    CHECK(dt.num_order == 1 && dt.den_order == 1, "order preserved");
    /* Check DC gain: G(1) = 1.0 */
    double u_hist[2] = {1.0, 0.0}, y_hist[2] = {0.0, 0.0};
    double y = cascade_dtf_execute(&dt, 1.0, u_hist, y_hist);
    CHECK_CLOSE(y, dt.num_coeff[0] + dt.num_coeff[1], 0.1, "DC gain check");
    cascade_tf_free(&ct); cascade_dtf_free(&dt);
}

/* ==========================================================================
 * L2: Anti-Windup and Bumpless Transfer Tests
 * ========================================================================== */

static void test_aw_clamping(void) {
    TEST("aw_clamping");
    CascadePID pid;
    cascade_pid_init(1.0, 1.0, 0.0, 10.0, 0.01, &pid);
    cascade_pid_set_limits(-5.0, 5.0, &pid);
    pid.integrator = 20.0;  /* Force integrator beyond limit */
    cascade_aw_clamping(&pid);
    CHECK(pid.integrator <= 5.0, "integrator should be clamped");
}

static void test_bumpless_transfer(void) {
    TEST("bumpless_transfer");
    CascadeBumpless bt;
    cascade_bt_init(1, 2.0, &bt);
    CHECK(bt.current_mode == 1, "initial mode should be 1");
    CHECK(bt.in_transition == 0, "should not be in transition");
    CascadeSystem sys;
    cascade_system_init(&sys, "test");
    int ret = cascade_bt_start(2, &bt, &sys, 0.0);
    CHECK(ret == 0, "bt_start should succeed");
    CHECK(bt.in_transition == 1, "should be in transition");
    double u = cascade_bt_update(&bt, &sys, 1.0, 3.0, 5.0);
    CHECK(u >= 3.0 && u <= 5.0, "blended u should be between old and new");
    cascade_system_free(&sys);
}

/* ==========================================================================
 * L5: Performance Metrics Tests
 * ========================================================================== */

static void test_performance_indices(void) {
    TEST("performance_indices");
    double t[] = {0.0, 0.5, 1.0, 1.5, 2.0};
    double e[] = {1.0, 0.6, 0.3, 0.1, 0.0};
    int n = 5;
    double ise = cascade_compute_ise(e, t, n);
    double iae = cascade_compute_iae(e, t, n);
    double itae = cascade_compute_itae(e, t, n);
    CHECK(ise >= 0.0, "ISE should be non-negative");
    CHECK(iae >= 0.0, "IAE should be non-negative");
    CHECK(itae >= 0.0, "ITAE should be non-negative");
    CHECK(iae >= ise, "IAE >= ISE for error < 1");
}

static void test_tv_computation(void) {
    TEST("total_variation");
    double u[] = {0.0, 1.0, 0.5, 2.0, 1.0};
    double tv = cascade_compute_tv(u, 5);
    CHECK_CLOSE(tv, 4.0, 1e-10, "TV should be 4.0");
}

/* ==========================================================================
 * L6: Canonical System Model Tests
 * ========================================================================== */

static void test_dc_motor_model(void) {
    TEST("dc_motor_model");
    CascadeDCMotor motor = {1.0, 0.001, 0.05, 0.05, 0.001, 0.0001, 1.0};
    CascadeTF vel_tf, pos_tf;
    int ret = cascade_create_dc_motor_model(&motor, &vel_tf, &pos_tf);
    CHECK(ret == 0, "DC motor model creation should succeed");
    CHECK(cascade_is_stable(&vel_tf) == 1, "velocity TF should be stable");
    cascade_tf_free(&vel_tf); cascade_tf_free(&pos_tf);
}

static void test_reactor_model(void) {
    TEST("reactor_model");
    CascadeReactor rx = {1.0, 0.2, 1000.0, 4200.0, 500.0, -50000.0,
                         1e10, 80000.0, 8.314, 300.0, 280.0, 0.001, 0.002};
    CascadeTF jacket_tf, reactor_tf;
    int ret = cascade_create_reactor_model(&rx, &jacket_tf, &reactor_tf);
    CHECK(ret == 0, "reactor model creation should succeed");
    CHECK(cascade_is_stable(&jacket_tf) == 1, "jacket TF should be stable");
    cascade_tf_free(&jacket_tf); cascade_tf_free(&reactor_tf);
}

static void test_flow_pressure_model(void) {
    TEST("flow_pressure_model");
    CascadeFlowPressure fp = {100.0, 0.1, 0.0001, 1000.0, 0.001,
                               10.0, 50.0, 0.01, 5.0, 200000.0};
    CascadeTF pressure_tf, flow_tf;
    int ret = cascade_create_flow_pressure_model(&fp, &pressure_tf, &flow_tf);
    CHECK(ret == 0, "flow-pressure model should succeed");
    cascade_tf_free(&pressure_tf); cascade_tf_free(&flow_tf);
}

static void test_level_tank_model(void) {
    TEST("level_tank_model");
    CascadeLevelTank tank = {2.0, 5.0, 1.0, 0.01, 5.0, 0.1, 3.0};
    CascadeTF flow_tf, level_tf;
    int ret = cascade_create_level_tank_model(&tank, &flow_tf, &level_tf);
    CHECK(ret == 0, "level tank model should succeed");
    cascade_tf_free(&flow_tf); cascade_tf_free(&level_tf);
}

/* ==========================================================================
 * L5: Frequency Response Tests
 * ========================================================================== */

static void test_bode_analysis(void) {
    TEST("bode_analysis");
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF tf = cascade_tf_create(num, 0, den, 1, 1.0);
    CascadeFreqResponse resp;
    int ret = cascade_bode_analysis(&tf, 0.1, 10.0, 50, &resp);
    CHECK(ret == 0, "bode analysis should succeed");
    CHECK(resp.num_points == 50, "should have 50 points");
    CHECK(resp.dc_gain_db < 0.1 && resp.dc_gain_db > -0.1, "DC gain should be ~0 dB");
    cascade_freq_response_free(&resp);
    cascade_tf_free(&tf);
}

static void test_stability_margins(void) {
    TEST("stability_margins");
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF tf = cascade_tf_create(num, 0, den, 1, 1.0);
    CascadeFreqResponse resp;
    cascade_bode_analysis(&tf, 0.01, 100.0, 100, &resp);
    double gm, pm;
    int ret = cascade_stability_margins(&resp, &gm, &pm);
    CHECK(ret == 0, "margin computation should succeed");
    /* First-order system: infinite GM, PM = 90 deg at low frequencies */
    CHECK(gm > 0.0, "GM should be positive");
    cascade_freq_response_free(&resp);
    cascade_tf_free(&tf);
}

static void test_closed_loop_bandwidth(void) {
    TEST("closed_loop_bandwidth");
    /* G_cl = 1/(s+1): bandwidth = 1 rad/s */
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF cl = cascade_tf_create(num, 0, den, 1, 1.0);
    double bw;
    int ret = cascade_closed_loop_bandwidth(&cl, &bw);
    CHECK(ret == 0, "bandwidth computation should succeed");
    CHECK_CLOSE(bw, 1.0, 0.1, "bandwidth should be ~1 rad/s");
    cascade_tf_free(&cl);
}

/* ==========================================================================
 * L5: Step Response Tests
 * ========================================================================== */

static void test_step_response_first_order(void) {
    TEST("step_response_1st_order");
    double num[] = {1.0}, den[] = {1.0, 1.0};
    CascadeTF tf = cascade_tf_create(num, 0, den, 1, 1.0);
    double y1 = cascade_step_response(&tf, 1.0, 10);
    CHECK_CLOSE(y1, 1.0 - exp(-1.0), 0.01, "y(1) = 1-e^-1");
    double y5 = cascade_step_response(&tf, 5.0, 10);
    CHECK_CLOSE(y5, 1.0, 0.01, "y(5) should be near steady state");
    cascade_tf_free(&tf);
}

static void test_step_response_second_order(void) {
    TEST("step_response_2nd_order");
    /* G(s) = 1/(s^2 + s + 1): wn=1, zeta=0.5, underdamped */
    double num[] = {1.0}, den[] = {1.0, 1.0, 1.0};
    CascadeTF tf = cascade_tf_create(num, 0, den, 2, 1.0);
    double y_final = cascade_step_response(&tf, 50.0, 10);
    CHECK_CLOSE(y_final, 1.0, 0.05, "steady state should be 1");
    /* At t = pi/wd = pi/sqrt(0.75) ~ 3.63: should be first peak > 1 */
    double t_peak = M_PI / sqrt(0.75);
    double y_peak = cascade_step_response(&tf, t_peak, 10);
    CHECK(y_peak > 1.0, "should overshoot");
    cascade_tf_free(&tf);
}

/* ==========================================================================
 * L5: Optimization Test
 * ========================================================================== */

static void test_optimize_simultaneous(void) {
    TEST("optimize_simultaneous");
    double nf[] = {1.0}, df[] = {1.0, 0.2};
    double ns[] = {1.0}, ds[] = {1.0, 2.0};
    CascadeTF inner = cascade_tf_create(nf, 0, df, 1, 1.0);
    CascadeTF outer = cascade_tf_create(ns, 0, ds, 1, 1.0);
    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    CascadeSystem sys;
    int ret = cascade_optimize_simultaneous(&inner, &outer, &spec, &sys);
    CHECK(ret == 0, "optimization should succeed");
    CHECK(sys.inner.inner_is_stable == 1, "inner should be stable");
    cascade_tf_free(&inner); cascade_tf_free(&outer);
    cascade_system_free(&sys);
}

/* ==========================================================================
 * L4: Internal Stability Test
 * ========================================================================== */

static void test_internal_stability(void) {
    TEST("internal_stability");
    double nf[] = {1.0}, df[] = {1.0, 0.2};
    double ns[] = {1.0}, ds[] = {1.0, 2.0};
    CascadeTF Gi = cascade_tf_create(nf, 0, df, 1, 1.0);
    CascadeTF Go = cascade_tf_create(ns, 0, ds, 1, 1.0);
    CascadePID Ci, Co;
    cascade_pid_init(1.0, 5.0, 0.0, 10.0, 0.0, &Ci);
    cascade_pid_init(0.5, 0.2, 0.0, 10.0, 0.0, &Co);
    int stable = cascade_verify_internal_stability(&Gi, &Ci, &Go, &Co);
    CHECK(stable == 1, "should be internally stable");
    cascade_tf_free(&Gi); cascade_tf_free(&Go);
}

/* ==========================================================================
 * L5: Performance Comparison Test
 * ========================================================================== */

static void test_cascade_vs_single(void) {
    TEST("cascade_vs_single_loop");
    double nf[] = {1.0}, df[] = {1.0, 0.2};
    double ns[] = {1.0}, ds[] = {1.0, 2.0};
    CascadeTF Gi = cascade_tf_create(nf, 0, df, 1, 1.0);
    CascadeTF Go = cascade_tf_create(ns, 0, ds, 1, 1.0);
    CascadePID Ci, Co;
    cascade_pid_init(1.0, 5.0, 0.0, 10.0, 0.0, &Ci);
    cascade_pid_init(0.5, 0.2, 0.0, 10.0, 0.0, &Co);
    double improvement;
    int ret = cascade_compare_single_vs_cascade(&Gi, &Go, &Ci, &Co,
                                                  10.0, 100, &improvement);
    CHECK(ret == 0, "comparison should succeed");
    CHECK(improvement > 0.5, "improvement factor should be reasonable");
    cascade_tf_free(&Gi); cascade_tf_free(&Go);
}

/* ==========================================================================
 * L5: Validate Design Test
 * ========================================================================== */

static void test_validate_design(void) {
    TEST("validate_design");
    double nf[] = {1.0}, df[] = {1.0, 0.2};
    double ns[] = {1.0}, ds[] = {1.0, 2.0};
    CascadeTF inner = cascade_tf_create(nf, 0, df, 1, 1.0);
    CascadeTF outer = cascade_tf_create(ns, 0, ds, 1, 1.0);
    CascadeDesignSpec spec;
    cascade_set_default_spec(&spec);
    CascadeSystem sys;
    cascade_design_sequential(&inner, &outer, &spec, &sys);
    int failures = cascade_validate_design(&sys, &spec);
    CHECK(failures == 0, "design should pass validation");
    cascade_tf_free(&inner); cascade_tf_free(&outer);
    cascade_system_free(&sys);
}

/* ==========================================================================
 * L6: Application Assessment Tests
 * ========================================================================== */

static void test_assess_dc_motor(void) {
    TEST("assess_dc_motor");
    CascadeDCMotor motor = {1.0, 0.001, 0.05, 0.05, 0.001, 0.0001, 1.0};
    CascadePerformance perf;
    int ret = cascade_assess_dc_motor(&motor, &perf);
    CHECK(ret == 0, "DC motor assessment should succeed");
    CHECK(perf.inner_settle_time > 0, "inner settle time should be positive");
    CHECK(perf.dist_rejection_ratio > 0, "dist rejection ratio should be positive");
}

static void test_assess_reactor(void) {
    TEST("assess_reactor");
    CascadeReactor rx = {1.0, 0.2, 1000.0, 4200.0, 500.0, -50000.0,
                         1e10, 80000.0, 8.314, 300.0, 280.0, 0.001, 0.002};
    CascadePerformance perf;
    int ret = cascade_assess_reactor(&rx, &perf);
    CHECK(ret == 0, "reactor assessment should succeed");
    CHECK(perf.outer_settle_time > 0, "outer settle time should be positive");
}

static void test_assess_level_tank(void) {
    TEST("assess_level_tank");
    CascadeLevelTank tank = {2.0, 5.0, 1.0, 0.01, 5.0, 0.1, 3.0};
    CascadePerformance perf;
    int ret = cascade_assess_level_tank(&tank, &perf);
    CHECK(ret == 0, "level tank assessment should succeed");
    CHECK(perf.dist_rejection_ratio >= 3.0, "dist rejection should be >= 3");
}

/* ==========================================================================
 * Main Test Runner
 * ========================================================================== */

int main(void) {
    printf("=== Cascade Control Module Test Suite ===\n\n");

    printf("--- L3: Polynomial Operations ---\n");
    test_polynomial_create();
    test_polynomial_eval();
    test_polynomial_mul();
    test_polynomial_add();

    printf("\n--- L1: Transfer Functions ---\n");
    test_tf_create_and_freq_resp();
    test_tf_close_loop();

    printf("\n--- L1: PID Controllers ---\n");
    test_pid_to_tf();
    test_pid_freq_response();

    printf("\n--- L4: Stability Analysis ---\n");
    test_routh_hurwitz_stable();
    test_routh_hurwitz_unstable();
    test_is_stable();
    test_internal_stability();

    printf("\n--- L2: Cascade Loop Formation ---\n");
    test_inner_closed_loop();
    test_equivalent_plant();

    printf("\n--- L5: Cascade Design ---\n");
    test_sequential_design();
    test_skogestad_tuning();
    test_direct_synthesis();
    test_validate_design();
    test_optimize_simultaneous();

    printf("\n--- L5: Digital Implementation ---\n");
    test_pid_velocity_form();
    test_pid_position_form();
    test_tustin_discretize();

    printf("\n--- L2: Anti-Windup & Bumpless ---\n");
    test_aw_clamping();
    test_bumpless_transfer();

    printf("\n--- L5: Frequency Response ---\n");
    test_bode_analysis();
    test_stability_margins();
    test_closed_loop_bandwidth();

    printf("\n--- L5: Step Response ---\n");
    test_step_response_first_order();
    test_step_response_second_order();

    printf("\n--- L5: Performance Metrics ---\n");
    test_performance_indices();
    test_tv_computation();
    test_cascade_vs_single();

    printf("\n--- L6: Canonical Models ---\n");
    test_dc_motor_model();
    test_reactor_model();
    test_flow_pressure_model();
    test_level_tank_model();

    printf("\n--- L6: Application Assessment ---\n");
    test_assess_dc_motor();
    test_assess_reactor();
    test_assess_level_tank();

    printf("\n========================================\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("========================================\n");

    return (tests_failed > 0) ? 1 : 0;
}
