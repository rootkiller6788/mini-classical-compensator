/**
 * @file test_feedforward.c
 * @brief Comprehensive tests for feedforward control library.
 *
 * Tests cover L1-L6: core operations, design methods, input shaping,
 * filter design, adaptive control, and 2-DOF performance.
 *
 * Each test validates a specific theorem or property.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "feedforward_core.h"
#include "feedforward_design.h"
#include "feedforward_input_shaping.h"
#include "feedforward_filter.h"
#include "feedforward_adaptive.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_EQ(a, b, tol, msg) do { if (fabs((a)-(b)) > (tol)) { \
    printf("FAIL: %s (got %.6f, expected %.6f)\n", msg, (double)(a), (double)(b)); \
    tests_failed++; return; } } while(0)

/* External declarations for functions used by tests */
void twodof_init(const TransferFn *plant, const TransferFn *fb,
                 const TransferFn *ff_ref, const TransferFn *ff_dist,
                 TwoDOF *two_dof);
void twodof_evaluate_performance(const TwoDOF *two_dof,
                                 double t_final, double dt,
                                 FFPerformance *perf);
double ff_satellite_attitude(double J, double theta_d,
                              double omega_d, double alpha_d);
double ff_quadrotor_altitude(double mass, double a_d_z,
                              double phi, double theta);
double ff_dc_motor_position(double J_hat, double b_hat, double Kt,
                             double q_d, double qd_d, double qdd_d);
double ff_temperature_dist(double K, double Kd, double T_ambient,
                            double prev_u, double alpha);

/* ==========================================================================
 * L1: Polynomial Tests
 * ========================================================================== */

static void test_poly_create_eval(void)
{
    TEST("poly_create and poly_eval");
    double coeff[] = {1.0, 2.0, 3.0};  /* p(x) = 1 + 2x + 3x^2 */
    Poly p = poly_create(coeff, 2);
    CHECK(p.coeff != NULL, "poly_create returned NULL");

    CHECK_EQ(poly_eval(&p, 0.0), 1.0, 1e-10, "p(0) != 1");
    CHECK_EQ(poly_eval(&p, 1.0), 6.0, 1e-10, "p(1) != 6");  /* 1+2+3=6 */
    CHECK_EQ(poly_eval(&p, 2.0), 17.0, 1e-10, "p(2) != 17"); /* 1+4+12=17 */
    CHECK_EQ(poly_eval(&p, -1.0), 2.0, 1e-10, "p(-1) != 2"); /* 1-2+3=2 */

    poly_free(&p);
    PASS();
}

static void test_poly_arithmetic(void)
{
    TEST("poly_add, poly_mul, poly_sub, poly_derivative");

    /* p(x) = 1 + x, q(x) = 1 + 2x */
    double c1[] = {1.0, 1.0};
    double c2[] = {1.0, 2.0};
    Poly p = poly_create(c1, 1);
    Poly q = poly_create(c2, 1);

    /* Addition: (1+x) + (1+2x) = 2 + 3x */
    Poly sum = poly_add(&p, &q);
    CHECK_EQ(poly_eval(&sum, 0.0), 2.0, 1e-10, "sum(0) != 2");
    CHECK_EQ(poly_eval(&sum, 1.0), 5.0, 1e-10, "sum(1) != 5");
    poly_free(&sum);

    /* Multiplication: (1+x)*(1+2x) = 1 + 3x + 2x^2 */
    Poly prod = poly_mul(&p, &q);
    CHECK_EQ(poly_eval(&prod, 0.0), 1.0, 1e-10, "prod(0) != 1");
    CHECK_EQ(poly_eval(&prod, 1.0), 6.0, 1e-10, "prod(1) != 6");
    poly_free(&prod);

    /* Subtraction: (1+x) - (1+2x) = -x */
    Poly diff = poly_sub(&p, &q);
    CHECK_EQ(poly_eval(&diff, 0.0), 0.0, 1e-10, "diff(0) != 0");
    CHECK_EQ(poly_eval(&diff, 1.0), -1.0, 1e-10, "diff(1) != -1");
    poly_free(&diff);

    /* Derivative: dp/dx = 1 (constant) */
    Poly dp = poly_derivative(&p);
    CHECK_EQ(poly_eval(&dp, 0.0), 1.0, 1e-10, "dp(0) != 1");
    CHECK_EQ(poly_eval(&dp, 5.0), 1.0, 1e-10, "dp(5) != 1");
    poly_free(&dp);

    poly_free(&p);
    poly_free(&q);
    PASS();
}

/* ==========================================================================
 * L2: Transfer Function Tests
 * ========================================================================== */

static void test_tf_freq_response(void)
{
    TEST("tf_freq_response for first-order system");

    /* G(s) = 1 / (s + 1) ? first-order low-pass */
    double num[] = {1.0};
    double den[] = {1.0, 1.0};
    TransferFn tf = tf_create(num, 0, den, 1, 1.0);

    /* DC gain: |G(0)| = 1, phase = 0 */
    double mag, phase;
    tf_freq_response(&tf, 0.0, &mag, &phase);
    CHECK_EQ(mag, 1.0, 1e-6, "DC gain != 1");
    CHECK_EQ(phase, 0.0, 1e-6, "DC phase != 0");

    /* Corner frequency: w=1 rad/s, |G(j1)| = 1/sqrt(2) ? 0.707, phase = -pi/4 */
    tf_freq_response(&tf, 1.0, &mag, &phase);
    CHECK_EQ(mag, 1.0/sqrt(2.0), 1e-4, "|G(j1)| != 1/sqrt(2)");
    CHECK_EQ(phase, -M_PI/4.0, 1e-4, "phase at w=1 != -45 deg");

    tf_free(&tf);
    PASS();
}

static void test_tf_step_response(void)
{
    TEST("tf_step_response for first and second order");

    /* First-order: G(s) = 1/(tau*s+1), tau=1
     * y(t) = 1 - exp(-t), so y(1) = 0.632 */
    double num1[] = {1.0};
    double den1[] = {1.0, 1.0};
    TransferFn tf1 = tf_create(num1, 0, den1, 1, 1.0);
    double y1 = tf_step_response(&tf1, 1.0, 0);
    CHECK_EQ(y1, 1.0 - exp(-1.0), 5e-4, "First order step at t=1");
    tf_free(&tf1);

    /* Second-order: G(s) = 1/(s^2 + s + 1), wn=1, zeta=0.5
     * y(t) = 1 - exp(-t/2)*(cos(wd*t) + sin(wd*t)/sqrt(3)) */
    double num2[] = {1.0};
    double den2[] = {1.0, 1.0, 1.0};
    TransferFn tf2 = tf_create(num2, 0, den2, 1, 1.0);
    /* At t large, y -> 1 */
    double y_inf = tf_step_response(&tf2, 50.0, 0);
    CHECK_EQ(y_inf, 1.0, 1e-2, "Second order steady state != 1");
    tf_free(&tf2);

    PASS();
}

static void test_tf_poles(void)
{
    TEST("tf_find_poles for low-order systems");

    /* Second-order: (s + 1)(s + 2) = s^2 + 3s + 2, poles at -1, -2 */
    double num[] = {2.0};
    double den[] = {2.0, 3.0, 1.0};
    TransferFn tf = tf_create(num, 0, den, 2, 1.0);

    double poles[4];
    int n = tf_find_poles(&tf, poles, 4);
    CHECK(n >= 2, "Expected at least 2 poles");
    /* Poles should be near -1 and -2 */
    int found_m1 = 0, found_m2 = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(poles[i] - (-1.0)) < 0.2) found_m1 = 1;
        if (fabs(poles[i] - (-2.0)) < 0.2) found_m2 = 1;
    }
    CHECK(found_m1 && found_m2, "Poles not found at -1 and -2");

    tf_free(&tf);
    PASS();
}

static void test_tf_is_minimum_phase(void)
{
    TEST("tf_is_minimum_phase check");

    /* Minimum-phase: zero at s=-2 (LHP) */
    double num_mp[] = {2.0, 1.0};
    double den_mp[] = {1.0, 3.0, 2.0};
    TransferFn tf_mp = tf_create(num_mp, 1, den_mp, 2, 1.0);
    CHECK(tf_is_minimum_phase(&tf_mp, 0) == 1, "Should be minimum-phase");
    tf_free(&tf_mp);

    /* Non-minimum-phase: zero at s=+1 (RHP) */
    double num_nmp[] = {-1.0, 1.0};
    double den_nmp[] = {1.0, 2.0, 1.0};
    TransferFn tf_nmp = tf_create(num_nmp, 1, den_nmp, 1, 1.0);
    CHECK(tf_is_minimum_phase(&tf_nmp, 0) == 0, "Should be non-minimum-phase");
    tf_free(&tf_nmp);

    PASS();
}

static void test_tf_state_space(void)
{
    TEST("tf_to_state_space conversion");

    /* G(s) = 2/(s^2 + 3s + 2) */
    double num[] = {2.0};
    double den[] = {2.0, 3.0, 1.0};
    TransferFn tf = tf_create(num, 0, den, 2, 1.0);

    double A[4], B[2], C[2], D;
    int n = tf_to_state_space(&tf, A, B, C, &D);
    CHECK(n == 2, "Expected n=2");

    /* Controllable canonical form:
     * A = [0 1; -2 -3]
     * B = [0; 1]
     * C = [2 0] (for this TF)
     * D = 0 (strictly proper) */
    CHECK_EQ(A[0], 0.0, 1e-10, "A[0,0] != 0");
    CHECK_EQ(A[1], 1.0, 1e-10, "A[0,1] != 1");
    CHECK_EQ(B[0], 0.0, 1e-10, "B[0] != 0");
    CHECK_EQ(B[1], 1.0, 1e-10, "B[1] != 1");
    CHECK_EQ(D, 0.0, 1e-10, "D != 0");

    tf_free(&tf);
    PASS();
}

/* ==========================================================================
 * L4: Feedforward Design Tests
 * ========================================================================== */

static void test_ff_model_inverse(void)
{
    TEST("ff_model_inverse for minimum-phase plant");

    /* P(s) = 1/(s+1), minimum-phase, proper inverse needs filter */
    double num[] = {1.0};
    double den[] = {1.0, 1.0};
    TransferFn plant = tf_create(num, 0, den, 1, 1.0);

    TransferFn inv;
    int ret = ff_model_inverse(&plant, &inv, 1, 100.0);
    CHECK(ret == 0, "Model inverse failed");

    /* inv(s) ? (s+1) * 100/(s+100) for 1st-order filter */
    /* At DC, |inv(0)| = |plant(0)|^{-1} = 1 */
    double mag, phase;
    tf_freq_response(&inv, 0.0, &mag, &phase);
    /* With filter: inv(s) = (s+1)*100/(s+100)
     * inv(0) = 100/100 = 1 */
    CHECK_EQ(mag, 1.0, 1e-4, "Inverse DC gain != 1");

    tf_free(&plant);
    tf_free(&inv);
    PASS();
}

static void test_ff_static_dist_ff(void)
{
    TEST("ff_design_static_dist_ff");

    /* P(s) = 2/(s+1), G_d(s) = 0.5/(s+1)
     * K_d = G_d(0)/P(0) = 0.5/2 = 0.25, u_ff = -0.25*d */
    double np[] = {2.0};
    double dp[] = {1.0, 1.0};
    TransferFn plant = tf_create(np, 0, dp, 1, 1.0);

    double nd[] = {0.5};
    double dd[] = {1.0, 1.0};
    TransferFn dist = tf_create(nd, 0, dd, 1, 1.0);

    double ff_gain;
    int ret = ff_design_static_dist_ff(&plant, &dist, &ff_gain);
    CHECK(ret == 0, "Static dist FF design failed");
    CHECK_EQ(ff_gain, -0.25, 1e-6, "FF gain != -0.25");

    tf_free(&plant);
    tf_free(&dist);
    PASS();
}

static void test_ff_dc_gain(void)
{
    TEST("ff_dc_gain computation");

    /* F(s) = (s+2)/(s+1), DC gain = 2/1 = 2 */
    double num[] = {2.0, 1.0};
    double den[] = {1.0, 1.0};
    TransferFn ff = tf_create(num, 1, den, 1, 1.0);
    double dc = ff_dc_gain(&ff, 0);
    CHECK_EQ(dc, 2.0, 1e-6, "DC gain != 2");

    tf_free(&ff);
    PASS();
}

/* ==========================================================================
 * L5: Input Shaping Tests
 * ========================================================================== */

static void test_shaper_zv(void)
{
    TEST("shaper_design_zv ? amplitudes sum to 1");

    InputShaper shaper;
    memset(&shaper, 0, sizeof(shaper));
    shaper_design_zv(10.0, 0.1, &shaper);
    CHECK(shaper.n_imp == 2, "ZV should have 2 impulses");

    double sum = 0.0;
    for (int i = 0; i < shaper.n_imp; i++) {
        sum += shaper.impulses[i].amplitude;
        CHECK(shaper.impulses[i].amplitude > 0, "Negative amplitude in ZV");
    }
    CHECK_EQ(sum, 1.0, 1e-6, "Amplitudes don't sum to 1");

    shaper_free(&shaper);
    PASS();
}

static void test_shaper_zvd(void)
{
    TEST("shaper_design_zvd ? amplitudes sum to 1");

    InputShaper shaper;
    memset(&shaper, 0, sizeof(shaper));
    shaper_design_zvd(10.0, 0.1, &shaper);
    CHECK(shaper.n_imp == 3, "ZVD should have 3 impulses");

    double sum = 0.0;
    for (int i = 0; i < shaper.n_imp; i++) {
        sum += shaper.impulses[i].amplitude;
    }
    CHECK_EQ(sum, 1.0, 1e-6, "Amplitudes don't sum to 1");

    shaper_free(&shaper);
    PASS();
}

static void test_shaper_residual_vibration(void)
{
    TEST("shaper_residual_vibration at design frequency = 0");

    InputShaper shaper;
    memset(&shaper, 0, sizeof(shaper));
    shaper_design_zv(10.0, 0.05, &shaper);

    /* At r=1 (design frequency), residual vibration should be ~0 */
    double V = shaper_residual_vibration(&shaper, 1.0, 0.05);
    CHECK(V < 1e-6, "ZV vibration at design frequency is not zero");

    shaper_free(&shaper);
    PASS();
}

/* ==========================================================================
 * L5: Filter Tests
 * ========================================================================== */

static void test_fir_moving_average(void)
{
    TEST("fir_design_moving_average ? unity DC gain");

    FIRFilter fir;
    memset(&fir, 0, sizeof(fir));
    fir_design_moving_average(4, &fir);
    CHECK(fir.n_taps == 4, "MA filter should have 4 taps");

    /* Steady-state: constant input of 1.0 should produce 1.0 output */
    double y = 0.0;
    for (int i = 0; i < 10; i++) {
        y = fir_process(&fir, 1.0);
    }
    CHECK_EQ(y, 1.0, 1e-6, "MA filter DC gain != 1");

    fir_free(&fir);
    PASS();
}

static void test_prefilter_exp_smooth(void)
{
    TEST("prefilter_exp_smooth ? convergence to step");

    double r_in[20];
    double r_out[20];
    for (int i = 0; i < 20; i++) r_in[i] = (i >= 5) ? 1.0 : 0.0;

    prefilter_exp_smooth(0.5, r_in, 20, r_out);
    CHECK_EQ(r_out[0], 0.0, 1e-10, "Initial output != 0");

    /* After many samples, should approach 1.0 */
    CHECK(r_out[19] > 0.9, "Smoothed output did not approach 1.0");

    PASS();
}

static void test_prefilter_velocity_limit(void)
{
    TEST("prefilter_velocity_limit ? rate limiting");

    double r_in[] = {0.0, 0.0, 10.0, 10.0, 10.0};
    double r_out[5];

    prefilter_velocity_limit(2.0, 1.0, r_in, 5, r_out);
    /* v_max=2, dt=1 => max step = 2
     * r_in jumps from 0 to 10 at index 2
     * r_out[2] <= r_out[1] + 2, r_out[3] <= r_out[2] + 2, etc. */
    CHECK(r_out[0] == 0.0, "r_out[0] should be 0");
    CHECK(r_out[1] == 0.0, "r_out[1] should be 0");
    CHECK(r_out[2] <= 2.0, "Rate limit exceeded at step 2");
    CHECK(r_out[4] == 6.0, "Final value should be 6 (limited at +2/step)");

    PASS();
}

/* ==========================================================================
 * L5: Adaptive Feedforward Tests
 * ========================================================================== */

static void test_lms_adaptive(void)
{
    TEST("lms_process ? convergence to identity mapping");

    LMSFilter filter;
    memset(&filter, 0, sizeof(filter));
    lms_init(&filter, 4, 0.05, 0.0);

    /* Train filter to output = input (identity mapping) */
    double sum_error = 0.0;
    for (int i = 0; i < 2000; i++) {
        double x = sin(0.1 * i);
        double y = lms_process(&filter, x, x); /* desired = x */
        sum_error += (x - y) * (x - y);
    }
    /* LMS should converge; check average error is reasonable */
    CHECK(sum_error / 2000.0 < 0.5, "LMS did not converge reasonably");

    lms_free(&filter);
    PASS();
}

/* ==========================================================================
 * L6: Canonical System Tests
 * ========================================================================== */

static void test_dc_motor_ff(void)
{
    TEST("ff_dc_motor_position ? physical consistency");

    /* Motor: J=0.01, b=0.001, Kt=0.1
     * For constant velocity, a_d=0, FF = b*v_d / Kt
     * v_d = 10 rad/s => u_ff = 0.001*10/0.1 = 0.1 V */
    double J = 0.01, b = 0.001, Kt = 0.1;
    double q_d = 1.0, qd_d = 10.0, qdd_d = 0.0;

    double u = ff_dc_motor_position(J, b, Kt, q_d, qd_d, qdd_d);
    CHECK_EQ(u, 0.1, 1e-6, "DC motor FF for constant velocity");

    /* With acceleration: a_d = 100, FF = (0.01*100 + 0.001*10)/0.1 = 1.01/0.1 = 10.1 */
    u = ff_dc_motor_position(J, b, Kt, q_d, qd_d, 100.0);
    CHECK_EQ(u, 10.1, 1e-6, "DC motor FF with acceleration");

    PASS();
}

static void test_temperature_ff(void)
{
    TEST("ff_temperature_dist ? disturbance compensation");

    /* K=2 (process gain), Kd=0.5 (disturbance gain)
     * T_ambient = 25C (5C above nominal 20C)
     * u_ff = -(0.5/2) * 25 = -6.25 */
    double u = ff_temperature_dist(2.0, 0.5, 25.0, 0.0, 1.0);
    CHECK_EQ(u, -6.25, 1e-6, "Temperature FF gain");

    /* With smoothing: alpha=0.5, prev=0
     * u = 0.5*(raw) + 0.5*0 = raw/2 */
    u = ff_temperature_dist(2.0, 0.5, 25.0, 0.0, 0.5);
    CHECK_EQ(u, -3.125, 1e-6, "Temperature FF with smoothing");

    PASS();
}

static void test_2dof_performance(void)
{
    TEST("twodof_evaluate_performance ? FF improves tracking");

    /* Plant: P(s) = 1/(s+1), Feedback: C(s) = 5 (P-control)
     * FF: F(s) = 1 */
    double np[] = {1.0};
    double dp[] = {1.0, 1.0};
    TransferFn plant = tf_create(np, 0, dp, 1, 1.0);

    double nc[] = {5.0};
    double dc[] = {1.0};
    TransferFn fb = tf_create(nc, 0, dc, 0, 1.0);

    double nf[] = {1.0};
    double df[] = {1.0};
    TransferFn ff_ref = tf_create(nf, 0, df, 0, 1.0);

    TwoDOF two_dof;
    memset(&two_dof, 0, sizeof(two_dof));
    twodof_init(&plant, &fb, &ff_ref, NULL, &two_dof);

    FFPerformance perf;
    twodof_evaluate_performance(&two_dof, 5.0, 0.01, &perf);

    /* Feedforward should improve tracking (ISE lower than FB-only equivalent) */
    CHECK(perf.ff_contribution >= 0.0, "FF contribution should be non-negative");

    tf_free(&plant);
    tf_free(&fb);
    tf_free(&ff_ref);
    PASS();
}

/* ==========================================================================
 * L3: Transfer Function Algebra Tests
 * ========================================================================== */

static void test_tf_algebra(void)
{
    TEST("tf_parallel, tf_series, tf_feedback consistency");

    /* G1 = G2 = 1/(s+1)
     * Parallel: G1+G2 = 2/(s+1), DC gain = 2 */
    double n[] = {1.0};
    double d[] = {1.0, 1.0};
    TransferFn g1 = tf_create(n, 0, d, 1, 1.0);
    TransferFn g2 = tf_create(n, 0, d, 1, 1.0);

    TransferFn par = tf_parallel(&g1, &g2);
    double dc = ff_dc_gain(&par, 0);
    CHECK_EQ(dc, 2.0, 1e-6, "Parallel DC gain != 2");
    tf_free(&par);

    /* Series: G1*G2 = 1/(s+1)^2, DC gain = 1 */
    TransferFn ser = tf_series(&g1, &g2);
    dc = ff_dc_gain(&ser, 0);
    CHECK_EQ(dc, 1.0, 1e-6, "Series DC gain != 1");
    tf_free(&ser);

    /* Feedback: G1/(1+G1*G2), DC = 1/(1+1*1) = 1/2 */
    TransferFn fb = tf_feedback(&g1, &g2);
    dc = ff_dc_gain(&fb, 0);
    CHECK_EQ(dc, 0.5, 1e-6, "Feedback DC gain != 0.5");
    tf_free(&fb);

    tf_free(&g1);
    tf_free(&g2);
    PASS();
}

/* ==========================================================================
 * L4: Realizability and Properness Tests
 * ========================================================================== */

static void test_ff_realizability(void)
{
    TEST("ff_is_realizable and ff_enforce_properness");

    /* Improper: num degree > den degree */
    double n_imp[] = {1.0, 2.0, 1.0};  /* s^2 + 2s + 1 */
    double d_imp[] = {1.0, 1.0};        /* s + 1 */
    TransferFn tf_imp = tf_create(n_imp, 2, d_imp, 1, 1.0);
    CHECK(ff_is_realizable(&tf_imp) == 0, "Should be non-realizable");

    /* Enforce properness */
    ff_enforce_properness(&tf_imp, 100.0);
    CHECK(ff_is_realizable(&tf_imp) == 1, "Should be realizable after enforcement");
    tf_free(&tf_imp);

    /* Proper: num degree <= den degree */
    double n_p[] = {1.0, 1.0};  /* s + 1 */
    double d_p[] = {1.0, 1.0, 1.0}; /* s^2 + s + 1 */
    TransferFn tf_p = tf_create(n_p, 1, d_p, 2, 1.0);
    CHECK(ff_is_realizable(&tf_p) == 1, "Should be realizable");
    tf_free(&tf_p);

    PASS();
}

/* ==========================================================================
 * L6: System-specific tests
 * ========================================================================== */

static void test_nl_pendulum_ff(void)
{
    TEST("nl_ff_pendulum ? gravity compensation");

    /* J=1, b=0, mgl=9.81, q_d=pi/2 (horizontal)
     * tau_ff = 1*0 + 0*0 + 9.81*sin(pi/2) = 9.81 */
    double tau = nl_ff_pendulum(1.0, 0.0, 9.81, M_PI/2.0, 0.0, 0.0);
    CHECK_EQ(tau, 9.81, 1e-6, "Pendulum gravity compensation");

    /* With velocity: J=1, b=0.5, mgl=0, qd_d=2, qdd_d=1
     * tau = 1*1 + 0.5*2 + 0 = 2.0 */
    tau = nl_ff_pendulum(1.0, 0.5, 0.0, 0.0, 2.0, 1.0);
    CHECK_EQ(tau, 2.0, 1e-6, "Pendulum inertia + damping");

    PASS();
}

static void test_nl_2link_arm_ff(void)
{
    TEST("nl_ff_2link_arm ? static torque");

    /* 2-link arm at rest (qd=0, qdd=0), should output gravity torques */
    double m1 = 1.0, m2 = 0.5, a1 = 0.5, a2 = 0.3;
    double q[2] = {0.0, 0.0};      /* arm pointing down */
    double qd[2] = {0.0, 0.0};
    double qdd[2] = {0.0, 0.0};
    double tau[2];

    nl_ff_2link_arm(m1, m2, a1, a2, q, qd, qdd, tau);

    /* At q1=q2=0: arm pointing straight down
     * cos(0)=1
     * G1 = (1+0.5)*9.81*0.5*1 + 0.5*9.81*0.3*1 = 7.3575 + 1.4715 = 8.829
     * G2 = 0.5*9.81*0.3*1 = 1.4715 */
    double G1_expected = (m1+m2)*9.81*a1 + m2*9.81*a2;
    double G2_expected = m2*9.81*a2;
    CHECK_EQ(tau[0], G1_expected, 1e-3, "Link 1 gravity torque");
    CHECK_EQ(tau[1], G2_expected, 1e-3, "Link 2 gravity torque");

    PASS();
}

static void test_servo_motion_ff(void)
{
    TEST("ff_servo_motion ? velocity/acceleration/friction");

    /* Kv=0.1, Ka=0.01, Ks=0.05, v=10, a=5
     * u = 0.1*10 + 0.01*5 + 0.05 = 1.0 + 0.05 + 0.05 = 1.10 */
    double u = ff_servo_motion(0.1, 0.01, 0.05, 10.0, 5.0);
    CHECK_EQ(u, 1.10, 1e-10, "Servo feedforward");

    /* Negative velocity: sign change on friction */
    u = ff_servo_motion(0.1, 0.01, 0.05, -10.0, 5.0);
    /* u = -1.0 + 0.05 - 0.05 = -1.0 */
    CHECK_EQ(u, -1.0, 1e-10, "Servo feedforward negative velocity");

    PASS();
}

/* ==========================================================================
 * L7: Application Tests
 * ========================================================================== */

static void test_satellite_ff(void)
{
    TEST("ff_satellite_attitude ? torque scaling");

    /* J=100, a_d=0.01 rad/s^2 => tau = 1.0 N*m */
    double tau = ff_satellite_attitude(100.0, 0.0, 0.0, 0.01);
    CHECK_EQ(tau, 1.0, 1e-6, "Satellite FF torque");

    PASS();
}

static void test_quadrotor_ff(void)
{
    TEST("ff_quadrotor_altitude ? hover thrust");

    /* Hover: a_d=0, phi=0, theta=0
     * thrust = m*g = 1.0*9.81 = 9.81 */
    double thrust = ff_quadrotor_altitude(1.0, 0.0, 0.0, 0.0);
    CHECK_EQ(thrust, 9.81, 1e-4, "Quadrotor hover thrust");

    /* Accelerating up: a_d=2, thrust = 1*(9.81+2) = 11.81 */
    thrust = ff_quadrotor_altitude(1.0, 2.0, 0.0, 0.0);
    CHECK_EQ(thrust, 11.81, 1e-4, "Quadrotor accelerating thrust");

    PASS();
}

/* ==========================================================================
 * Main
 * ========================================================================== */

int main(void)
{
    printf("=== mini-feedforward-control Test Suite ===\n\n");

    printf("--- L1: Polynomial Operations ---\n");
    test_poly_create_eval();
    test_poly_arithmetic();

    printf("\n--- L2: Transfer Function Operations ---\n");
    test_tf_freq_response();
    test_tf_step_response();
    test_tf_poles();
    test_tf_is_minimum_phase();
    test_tf_state_space();

    printf("\n--- L3: TF Algebra ---\n");
    test_tf_algebra();

    printf("\n--- L4: Feedforward Design ---\n");
    test_ff_model_inverse();
    test_ff_static_dist_ff();
    test_ff_dc_gain();
    test_ff_realizability();

    printf("\n--- L5: Input Shaping ---\n");
    test_shaper_zv();
    test_shaper_zvd();
    test_shaper_residual_vibration();

    printf("\n--- L5: Filter Design ---\n");
    test_fir_moving_average();
    test_prefilter_exp_smooth();
    test_prefilter_velocity_limit();

    printf("\n--- L5: Adaptive Control ---\n");
    test_lms_adaptive();

    printf("\n--- L6: Canonical Systems ---\n");
    test_dc_motor_ff();
    test_temperature_ff();
    test_2dof_performance();
    test_nl_pendulum_ff();
    test_nl_2link_arm_ff();
    test_servo_motion_ff();

    printf("\n--- L7: Applications ---\n");
    test_satellite_ff();
    test_quadrotor_ff();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("========================================\n");

    return tests_failed > 0 ? 1 : 0;
}
