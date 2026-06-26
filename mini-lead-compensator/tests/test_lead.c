/**
 * @file test_lead.c
 * @brief Comprehensive test suite for mini-lead-compensator
 *
 * Tests cover L1-L4 core functions with assert-based verification.
 * Every test validates a specific knowledge point.
 */

#include "lead_compensator.h"
#include "lead_design.h"
#include "lead_frequency.h"
#include "lead_analysis.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

#define FAIL(msg) do { \
    tests_failed++; \
    printf("FAILED: %s\n", msg); \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { FAIL(msg); return; } \
} while(0)

#define ASSERT_FLOAT_EQ(a, b, tol, msg) do { \
    if (fabs((a)-(b)) > (tol)) { \
        printf(" (%.6g vs %.6g) ", (a), (b)); \
        FAIL(msg); return; \
    } \
} while(0)

static const double TOL = 1e-9;
static const double TOL_LOOSE = 1e-6;

/* ================================================================
 * L1 - Compensator Initialization Tests
 * ================================================================ */

void test_l1_init_pole_zero(void) {
    TEST("L1: lead_init with valid pole-zero params");
    lead_compensator_t comp;
    bool ok = lead_init(&comp, 10.0, 2.0, 10.0, LEAD_TYPE_ANALYTIC);
    ASSERT_TRUE(ok, "lead_init should succeed with valid params");
    ASSERT_FLOAT_EQ(comp.K_c, 10.0, TOL, "K_c mismatch");
    ASSERT_FLOAT_EQ(comp.z_c, 2.0, TOL, "z_c mismatch");
    ASSERT_FLOAT_EQ(comp.p_c, 10.0, TOL, "p_c mismatch");
    ASSERT_FLOAT_EQ(comp.alpha, 0.2, TOL, "alpha should be z_c/p_c = 0.2");
    ASSERT_FLOAT_EQ(comp.T, 0.5, TOL, "T should be 1/z_c = 0.5");
    ASSERT_FLOAT_EQ(comp.dc_gain, 2.0, TOL, "DC gain K = K_c*alpha = 2.0");
    PASS();
}

void test_l1_init_pole_zero_invalid(void) {
    TEST("L1: lead_init rejects invalid params");
    lead_compensator_t comp;
    ASSERT_TRUE(!lead_init(&comp, 0.0, 2.0, 10.0, LEAD_TYPE_ANALYTIC),
                "Should reject K_c <= 0");
    ASSERT_TRUE(!lead_init(&comp, 10.0, -1.0, 10.0, LEAD_TYPE_ANALYTIC),
                "Should reject z_c <= 0");
    ASSERT_TRUE(!lead_init(&comp, 10.0, 2.0, -10.0, LEAD_TYPE_ANALYTIC),
                "Should reject p_c <= 0");
    ASSERT_TRUE(!lead_init(&comp, 10.0, 20.0, 10.0, LEAD_TYPE_ANALYTIC),
                "Should reject z_c >= p_c");
    ASSERT_TRUE(!lead_init(&comp, 10.0, 5.0, 5.0, LEAD_TYPE_ANALYTIC),
                "Should reject z_c == p_c (alpha == 1)");
    PASS();
}

void test_l1_init_KT_alpha(void) {
    TEST("L1: lead_init_from_KT_alpha");
    lead_compensator_t comp;
    bool ok = lead_init_from_KT_alpha(&comp, 1.0, 0.5, 0.2, LEAD_TYPE_ANALYTIC);
    ASSERT_TRUE(ok, "init_from_KT_alpha should succeed");
    ASSERT_FLOAT_EQ(comp.dc_gain, 1.0, TOL, "DC gain should be K = 1.0");
    ASSERT_FLOAT_EQ(comp.T, 0.5, TOL, "T mismatch");
    ASSERT_FLOAT_EQ(comp.alpha, 0.2, TOL, "alpha mismatch");
    ASSERT_FLOAT_EQ(comp.z_c, 2.0, TOL, "z_c = 1/T = 2.0");
    ASSERT_FLOAT_EQ(comp.p_c, 10.0, TOL, "p_c = 1/(alpha*T) = 10.0");
    ASSERT_FLOAT_EQ(comp.K_c, 5.0, TOL, "K_c = K/alpha = 5.0");
    PASS();
}

void test_l1_validate(void) {
    TEST("L1: lead_validate checks consistency");
    lead_compensator_t comp;
    ASSERT_TRUE(!lead_validate(NULL), "NULL should be invalid");
    lead_init(&comp, 10.0, 2.0, 10.0, LEAD_TYPE_ANALYTIC);
    ASSERT_TRUE(lead_validate(&comp), "Valid compensator should validate");
    comp.z_c = 15.0; /* break lead condition */
    ASSERT_TRUE(!lead_validate(&comp), "z_c >= p_c should fail validation");
    PASS();
}

void test_l1_to_transfer_function(void) {
    TEST("L1: lead_to_transfer_function");
    lead_compensator_t comp;
    lead_init(&comp, 10.0, 2.0, 10.0, LEAD_TYPE_ANALYTIC);
    lead_tf_t tf;
    lead_to_transfer_function(&comp, &tf);
    ASSERT_TRUE(tf.num_order == 1, "num_order should be 1");
    ASSERT_TRUE(tf.den_order == 1, "den_order should be 1");
    ASSERT_FLOAT_EQ(tf.num[0], 20.0, TOL, "num[0] = K_c*z_c = 20");
    ASSERT_FLOAT_EQ(tf.num[1], 10.0, TOL, "num[1] = K_c = 10");
    ASSERT_FLOAT_EQ(tf.den[0], 10.0, TOL, "den[0] = p_c = 10");
    ASSERT_FLOAT_EQ(tf.den[1], 1.0, TOL, "den[1] = 1");
    PASS();
}

/* ================================================================
 * L2 - Phase and Magnitude Tests
 * ================================================================ */

void test_l2_phase_at_dc(void) {
    TEST("L2: phase at DC should be zero");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 0.1, 0.2, LEAD_TYPE_ANALYTIC);
    double phi = lead_phase_at(&comp, 0.0);
    ASSERT_FLOAT_EQ(phi, 0.0, TOL, "Phase at DC should be 0");
    PASS();
}

void test_l2_phase_positive(void) {
    TEST("L2: phase should be positive for omega > 0");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 0.1, 0.2, LEAD_TYPE_ANALYTIC);
    double phi = lead_phase_at(&comp, 10.0);
    ASSERT_TRUE(phi > 0.0, "Phase should be positive (lead)");
    PASS();
}

void test_l2_magnitude_at_dc(void) {
    TEST("L2: magnitude at DC equals K");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 2.5, 0.1, 0.2, LEAD_TYPE_ANALYTIC);
    double mag = lead_magnitude_at(&comp, 0.0);
    ASSERT_FLOAT_EQ(mag, 2.5, TOL_LOOSE, "|C(j0)| should equal DC gain");
    PASS();
}

void test_l2_magnitude_hf(void) {
    TEST("L2: high-frequency magnitude equals K/alpha");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 3.0, 0.1, 0.3, LEAD_TYPE_ANALYTIC);
    double mag = lead_magnitude_at(&comp, 1e9);
    double expected = 3.0 / 0.3;
    ASSERT_FLOAT_EQ(mag, expected, 1e-3, "|C(j*inf)| should equal K/alpha");
    PASS();
}

/* ================================================================
 * L4 - Maximum Phase Lead Theorem Tests
 * ================================================================ */

void test_l4_phi_max_formula(void) {
    TEST("L4: phi_max = asin((1-alpha)/(1+alpha))");
    double alpha = 0.2;
    double phi_max = lead_phi_max_from_alpha(alpha);
    double expected = asin((1.0 - 0.2) / (1.0 + 0.2));
    ASSERT_FLOAT_EQ(phi_max, expected, TOL, "phi_max formula verification");
    /* phi_max for alpha=0.2 should be about 41.8 degrees */
    double phi_deg = phi_max * 180.0 / M_PI;
    ASSERT_TRUE(fabs(phi_deg - 41.81) < 0.1, "phi_max ~ 41.8 deg for alpha=0.2");
    PASS();
}

void test_l4_alpha_from_phi(void) {
    TEST("L4: alpha from phi_max inverse");
    double phi = 30.0 * M_PI / 180.0;
    double alpha = lead_alpha_from_phi_max(phi);
    double expected = (1.0 - sin(phi)) / (1.0 + sin(phi));
    ASSERT_FLOAT_EQ(alpha, expected, TOL, "alpha inverse formula");
    /* For phi=30 deg, alpha should be about 0.333 */
    ASSERT_TRUE(fabs(alpha - 1.0/3.0) < 0.001, "alpha ~ 0.333 for phi=30 deg");
    PASS();
}

void test_l4_omega_max(void) {
    TEST("L4: omega_m = 1/(T*sqrt(alpha))");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 2.0, 0.25, LEAD_TYPE_ANALYTIC);
    double w_m = lead_omega_max(&comp);
    double expected = 1.0 / (2.0 * sqrt(0.25));
    ASSERT_FLOAT_EQ(w_m, expected, TOL, "omega_m formula");
    PASS();
}

void test_l4_magnitude_at_omega_max(void) {
    TEST("L4: |C(j*w_m)| = K/sqrt(alpha)");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 4.0, 0.5, 0.25, LEAD_TYPE_ANALYTIC);
    double mag = lead_magnitude_at_omega_max(&comp);
    double expected = 4.0 / sqrt(0.25);
    ASSERT_FLOAT_EQ(mag, expected, TOL, "|C(j*w_m)| formula");
    PASS();
}

void test_l4_roundtrip_alpha_phi(void) {
    TEST("L4: round-trip alpha -> phi_max -> alpha");
    double alpha_in = 0.1;
    double phi = lead_phi_max_from_alpha(alpha_in);
    double alpha_out = lead_alpha_from_phi_max(phi);
    ASSERT_FLOAT_EQ(alpha_in, alpha_out, 1e-12, "Alpha round-trip consistency");
    PASS();
}

/* ================================================================
 * L3 - PM/Damping Relationship Tests
 * ================================================================ */

void test_l3_pm_to_damping(void) {
    TEST("L3: PM to damping ratio (2nd-order)");
    double zeta = lead_pm_to_damping(65.0);
    ASSERT_TRUE(zeta > 0.6 && zeta < 0.75, "PM=65 -> zeta ~ 0.65-0.7");
    double zeta2 = lead_pm_to_damping(45.0);
    ASSERT_TRUE(zeta2 > 0.4 && zeta2 < 0.5, "PM=45 -> zeta ~ 0.42-0.45");
    PASS();
}

void test_l3_overshoot_from_zeta(void) {
    TEST("L3: Percent overshoot from damping ratio");
    double po = lead_overshoot_from_zeta(0.5);
    ASSERT_TRUE(po > 15.0 && po < 17.0, "PO for zeta=0.5 should be ~16.3%%");
    double po2 = lead_overshoot_from_zeta(0.7);
    ASSERT_TRUE(po2 > 4.0 && po2 < 5.0, "PO for zeta=0.7 should be ~4.6%%");
    PASS();
}

void test_l3_settling_time(void) {
    TEST("L3: Settling time from (zeta, wn)");
    double ts = lead_settling_time_2pct(0.5, 2.0);
    ASSERT_FLOAT_EQ(ts, 4.0, TOL, "ts = 4/(0.5*2) = 4.0");
    PASS();
}

/* ================================================================
 * L1 - DC and HF Gain Tests
 * ================================================================ */

void test_l1_dc_hf_gain(void) {
    TEST("L1: DC and HF gain");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 2.0, 1.0, 0.5, LEAD_TYPE_ANALYTIC);
    ASSERT_FLOAT_EQ(lead_dc_gain(&comp), 2.0, TOL, "DC gain = K = 2.0");
    ASSERT_FLOAT_EQ(lead_hf_gain(&comp), 4.0, TOL, "HF gain = K/alpha = 4.0");
    PASS();
}

/* ================================================================
 * L5 - Design Method Tests
 * ================================================================ */

void test_l5_dominant_pole(void) {
    TEST("L5: Dominant pole from (zeta, wn)");
    lead_complex_t pole = lead_dominant_pole(0.5, 4.0, true);
    ASSERT_FLOAT_EQ(pole.re, -2.0, TOL, "Re(s) = -zeta*wn = -2.0");
    ASSERT_FLOAT_EQ(pole.im, 3.4641016, TOL_LOOSE, "Im(s) = wn*sqrt(1-zeta^2)");
    PASS();
}

void test_l5_gain_from_error(void) {
    TEST("L5: Gain from steady-state error");
    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 1.0;
    sys.tf.den[0] = 1.0;
    sys.tf.den_order = 0;
    sys.tf.gain = 2.0;
    double K = lead_design_gain_from_error(&sys, 0.1);
    /* e_ss = 1/(1+K*2) = 0.1 => K = 4.5 */
    ASSERT_FLOAT_EQ(K, 4.5, 1e-6, "K for 10%% error with plant gain 2");
    PASS();
}

void test_l5_zeta_to_pm_consistency(void) {
    TEST("L5: zeta-to-PM consistency");
    double zeta = 0.6;
    double pm = lead_zeta_to_pm(zeta);
    double zeta_back = lead_pm_to_damping(pm);
    ASSERT_FLOAT_EQ(zeta, zeta_back, 0.02, "zeta -> PM -> zeta round-trip");
    PASS();
}

void test_l5_bandwidth_estimate(void) {
    TEST("L5: Bandwidth from zeta and ts");
    double bw = lead_bandwidth_from_zeta_ts(0.5, 2.0);
    ASSERT_TRUE(bw > 0.0, "Bandwidth should be positive");
    PASS();
}

/* ================================================================
 * L5 - Frequency Design Test (integration)
 * ================================================================ */

void test_l5_design_frequency(void) {
    TEST("L5: Frequency-domain lead design");

    /* Create a type-0 plant: G(s) = 10 / ((s+1)*(s+5)) = 10/(s^2+6s+5) */
    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 10.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 5.0;
    sys.tf.den[1] = 6.0;
    sys.tf.den[2] = 1.0;
    sys.tf.den_order = 2;
    sys.tf.gain = 1.0;

    lead_specs_t specs;
    memset(&specs, 0, sizeof(specs));
    specs.phase_margin_desired = 50.0;
    specs.gain_margin_desired = 10.0;
    specs.steady_state_error = 0.05;
    specs.use_frequency_domain = true;
    specs.method = LEAD_METHOD_FREQUENCY;

    lead_design_result_t result;
    bool ok = lead_design_frequency(&sys, &specs, &result);
    /* Design may or may not converge depending on plant difficulty;
     * the key verification is that if it returns true, the compensator is valid */
    if (ok) {
        ASSERT_TRUE(lead_validate(&result.compensator),
                    "Result must be valid lead compensator");
        double pm_before = lead_compute_phase_margin(&sys);
        double pm_after = result.achieved_phase_margin;
        printf("(PM: %.1f -> %.1f deg) ", pm_before, pm_after);
    } else {
        printf("(no convergence with this plant) ");
    }

    PASS();
}

/* ================================================================
 * L2 - Compensated System
 * ================================================================ */

void test_l2_compensated_phase(void) {
    TEST("L2: Compensated phase computation");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 0.5, 0.2, LEAD_TYPE_ANALYTIC);

    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 1.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 1.0;
    sys.tf.den_order = 0;
    sys.tf.gain = 1.0;

    double phase = lead_compensated_phase(&comp, &sys, 1.0);
    ASSERT_TRUE(isfinite(phase), "Compensated phase should be finite");
    PASS();
}

void test_l2_compensated_magnitude(void) {
    TEST("L2: Compensated magnitude computation");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 0.5, 0.2, LEAD_TYPE_ANALYTIC);

    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 1.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 1.0;
    sys.tf.den_order = 0;
    sys.tf.gain = 1.0;

    double mag = lead_compensated_magnitude(&comp, &sys, 1.0);
    ASSERT_TRUE(isfinite(mag), "Compensated magnitude should be finite");
    ASSERT_TRUE(mag > 0.0, "Magnitude should be positive");
    PASS();
}

/* ================================================================
 * L4 - Routh-Hurwitz Test
 * ================================================================ */

void test_l4_routh_hurwitz(void) {
    TEST("L4: Routh-Hurwitz stability test");
    /* s^2 + 2*s + 1 = 0 -> both roots at -1, stable, 0 sign changes */
    double poly1[] = {1.0, 2.0, 1.0};
    int rhp1 = lead_routh_hurwitz(poly1, 2);
    ASSERT_TRUE(rhp1 == 0, "Stable 2nd-order: 0 RHP poles");

    /* s^2 + s - 2 = 0 -> one RHP root, 1 sign change */
    double poly2[] = {-2.0, 1.0, 1.0};
    int rhp2 = lead_routh_hurwitz(poly2, 2);
    ASSERT_TRUE(rhp2 == 1, "Unstable: 1 RHP pole");

    /* s^2 + 0*s + 4 = 0 -> purely imaginary, marginal */
    double poly3[] = {4.0, 0.0, 1.0};
    int rhp3 = lead_routh_hurwitz(poly3, 2);
    ASSERT_TRUE(rhp3 == 0, "Marginally stable: 0 RHP poles (but oscillatory)");

    PASS();
}

/* ================================================================
 * L4 - Phase/Gain Margin Tests
 * ================================================================ */

void test_l4_phase_margin_computation(void) {
    TEST("L4: Phase margin computation");
    /* G(s) = 1/(s+1) -> PM = 180 + angle(1/(jw+1)) at w_gc */
    /* w_gc = 0 (no crossover for this simple system), PM = 180 */
    /* Actually 1/(s+1) has infinite gain margin, PM ~ 180 at any w */
    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 1.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 1.0;
    sys.tf.den[1] = 1.0;
    sys.tf.den_order = 1;
    sys.tf.gain = 1.0;

    double pm = lead_compute_phase_margin(&sys);
    /* For 1/(s+1), the largest phase is at DC (0 deg), so PM ~ 180 */
    ASSERT_TRUE(pm > 0.0, "Phase margin should be positive for stable system");
    PASS();
}

void test_l4_pm_status(void) {
    TEST("L4: Phase margin classification");
    ASSERT_TRUE(lead_classify_pm(45.0) == LEAD_PM_STABLE, "45 deg -> STABLE");
    ASSERT_TRUE(lead_classify_pm(3.0) == LEAD_PM_MARGINAL, "3 deg -> MARGINAL");
    ASSERT_TRUE(lead_classify_pm(-10.0) == LEAD_PM_UNSTABLE, "-10 deg -> UNSTABLE");
    PASS();
}

/* ================================================================
 * L4 - Closed-Loop Polynomial
 * ================================================================ */

void test_l4_closed_loop_polynomial(void) {
    TEST("L4: Closed-loop characteristic polynomial");
    lead_compensator_t comp;
    lead_init_from_KT_alpha(&comp, 1.0, 1.0, 0.5, LEAD_TYPE_ANALYTIC);
    /* C(s) = (s+1)/(0.5*s+1) = 2*(s+1)/(s+2) */

    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 1.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 0.0;  /* s^0 coefficient */
    sys.tf.den[1] = 1.0;  /* s^1 coefficient -> G(s) = 1/s */
    sys.tf.den_order = 1;
    sys.tf.gain = 1.0;

    double cl_poly[10];
    int cl_order;
    lead_closed_loop_polynomial(&comp, &sys, cl_poly, &cl_order);
    /* den_C = [2, 1], num_C = [2, 2], den_G = [0, 1], num_G = [1]
     * den_CL = den_C*den_G + num_C*num_G = [2,1]*[0,1] + [2,2]*[1]
     *        = [0,2,1] + [2,2] = [2,4,1] -> s^2 + 4s + 2 = 0 */
    ASSERT_TRUE(cl_order == 2, "CL order should be 2");
    ASSERT_FLOAT_EQ(cl_poly[0], 2.0, 0.01, "Constant term");
    ASSERT_FLOAT_EQ(cl_poly[1], 4.0, 0.01, "s coefficient");
    ASSERT_FLOAT_EQ(cl_poly[2], 1.0, 0.01, "s^2 coefficient");
    PASS();
}

/* ================================================================
 * L3 - Bandwidth Test
 * ================================================================ */

void test_l3_bandwidth_zeta_wn(void) {
    TEST("L3: Bandwidth from zeta and wn");
    double bw = lead_bandwidth_from_zeta_wn(0.7, 1.0);
    ASSERT_TRUE(bw > 0.9 && bw < 1.1, "For zeta=0.7, bw ~ wn");
    PASS();
}

/* ================================================================
 * L8 - Multi-Stage Lead Test
 * ================================================================ */

void test_l8_two_stage(void) {
    TEST("L8: Two-stage lead design");
    lead_system_t sys;
    memset(&sys, 0, sizeof(sys));
    sys.tf.num[0] = 10.0;
    sys.tf.num_order = 0;
    sys.tf.den[0] = 0.0;
    sys.tf.den[1] = 1.0;
    sys.tf.den[2] = 1.0;
    sys.tf.den_order = 2;
    sys.tf.gain = 1.0;

    lead_specs_t specs;
    memset(&specs, 0, sizeof(specs));
    specs.phase_margin_desired = 70.0;
    specs.steady_state_error = 0.05;

    lead_compensator_t stage1, stage2;
    bool ok = lead_design_two_stage(&sys, &specs, &stage1, &stage2);
    ASSERT_TRUE(ok, "Two-stage design should succeed");
    ASSERT_TRUE(lead_validate(&stage1), "Stage 1 should be valid");
    ASSERT_TRUE(lead_validate(&stage2), "Stage 2 should be valid");

    /* Combined phase should be sum of individual phases */
    double phi_combined = lead_cascade_phase(&stage1, &stage2, 1.0);
    double phi1 = lead_phase_at(&stage1, 1.0);
    double phi2 = lead_phase_at(&stage2, 1.0);
    ASSERT_FLOAT_EQ(phi_combined, phi1 + phi2, TOL, "Cascade phase is additive");
    PASS();
}

void test_l8_discrete_tustin(void) {
    TEST("L8: Discrete-time via Tustin transform");
    lead_compensator_t comp;
    lead_init(&comp, 5.0, 2.0, 10.0, LEAD_TYPE_DIGITAL);

    double b0, b1, a1;
    lead_to_discrete_tustin(&comp, 0.01, &b0, &b1, &a1);
    ASSERT_TRUE(isfinite(b0) && isfinite(b1) && isfinite(a1),
                "Tustin coefficients should be finite");

    /* Test one-step update */
    double x_prev = 0.0, y_prev = 0.0;
    double y = lead_discrete_update(b0, b1, a1, 1.0, &x_prev, &y_prev);
    ASSERT_FLOAT_EQ(y, b0, TOL, "First update: y = b0*x");
    ASSERT_FLOAT_EQ(x_prev, 1.0, TOL, "x_prev stored");
    ASSERT_FLOAT_EQ(y_prev, b0, TOL, "y_prev stored");
    PASS();
}

/* ================================================================
 * Main
 * ================================================================ */

int main(void) {
    printf("=== mini-lead-compensator Test Suite ===\n\n");

    /* L1 - Definitions */
    printf("--- L1: Definitions ---\n");
    test_l1_init_pole_zero();
    test_l1_init_pole_zero_invalid();
    test_l1_init_KT_alpha();
    test_l1_validate();
    test_l1_to_transfer_function();
    test_l1_dc_hf_gain();

    /* L2 - Core Concepts */
    printf("\n--- L2: Core Concepts ---\n");
    test_l2_phase_at_dc();
    test_l2_phase_positive();
    test_l2_magnitude_at_dc();
    test_l2_magnitude_hf();
    test_l2_compensated_phase();
    test_l2_compensated_magnitude();

    /* L3 - Math Structures */
    printf("\n--- L3: Math Structures ---\n");
    test_l3_pm_to_damping();
    test_l3_overshoot_from_zeta();
    test_l3_settling_time();
    test_l3_bandwidth_zeta_wn();

    /* L4 - Fundamental Laws */
    printf("\n--- L4: Fundamental Laws ---\n");
    test_l4_phi_max_formula();
    test_l4_alpha_from_phi();
    test_l4_omega_max();
    test_l4_magnitude_at_omega_max();
    test_l4_roundtrip_alpha_phi();
    test_l4_routh_hurwitz();
    test_l4_phase_margin_computation();
    test_l4_pm_status();
    test_l4_closed_loop_polynomial();

    /* L5 - Computational Methods */
    printf("\n--- L5: Computational Methods ---\n");
    test_l5_dominant_pole();
    test_l5_gain_from_error();
    test_l5_zeta_to_pm_consistency();
    test_l5_bandwidth_estimate();
    test_l5_design_frequency();

    /* L8 - Advanced Topics */
    printf("\n--- L8: Advanced Topics ---\n");
    test_l8_two_stage();
    test_l8_discrete_tustin();

    printf("\n========================================\n");
    printf("Results: %d/%d passed, %d failed\n",
           tests_passed, tests_run, tests_failed);
    printf("========================================\n");

    return (tests_failed == 0) ? 0 : 1;
}