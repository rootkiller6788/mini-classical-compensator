/**
 * @file test_lag.c
 * @brief Comprehensive tests for lag compensator module
 *
 * Tests all core APIs with mathematical assertions.
 * L4: At least 5 mathematical assertions beyond trivial assert(1).
 */

#include "lag_compensator.h"
#include "lag_design.h"
#include "lag_frequency.h"
#include "lag_identification.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>


/* External declarations for functions without dedicated headers */
extern int lag_design_dc_motor_speed(const LagDCMotorParams *params,
    double speed_accuracy_req, LagCompensator *comp);
extern int lag_simulate_step_first_order(const LagCompensator *lag,
    double K_plant, double tau_plant, double t_end, int n_points,
    LagStepResponse *resp);
extern void lag_step_response_free(LagStepResponse *resp);

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  TEST: %s ... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define EPS 1e-9

/* ==========================================================================
 * L1: Construction and accessors
 * ========================================================================== */

static void test_create(void) {
    TEST("lag_create basic");
    LagCompensator lag = lag_create(10.0, 0.5, 5.0);
    assert(lag_validate(&lag) == 0);
    assert(fabs(lag_get_dc_gain(&lag) - 10.0) < EPS);
    assert(fabs(lag_get_beta(&lag) - 5.0) < EPS);
    assert(fabs(lag_get_zero(&lag) + 2.0) < EPS);   /* zero = -1/T = -2 */
    assert(fabs(lag_get_pole(&lag) + 0.4) < EPS);    /* pole = -1/(beta*T) = -0.4 */
    PASS();
}

static void test_create_pole_zero(void) {
    TEST("lag_create_from_pole_zero");
    /* G_c(s) = 2 * (s+5)/(s+1): zero=-5, pole=-1, Kc_in=2 */
    /* Effective: T=1/5=0.2, beta=5/1=5, Kc_eff=2*5=10 */
    LagCompensator lag = lag_create_from_pole_zero(2.0, -5.0, -1.0);
    assert(lag_validate(&lag) == 0);
    assert(fabs(lag_get_dc_gain(&lag) - 10.0) < EPS);
    assert(fabs(lag_get_beta(&lag) - 5.0) < EPS);
    PASS();
}

static void test_create_corners(void) {
    TEST("lag_create_from_corners");
    LagCompensator lag = lag_create_from_corners(2.0, 5.0, 1.0);
    assert(lag_validate(&lag) == 0);
    assert(fabs(lag_get_beta(&lag) - 5.0) < EPS);
    PASS();
}

static void test_identity(void) {
    TEST("lag_create_identity");
    LagCompensator lag = lag_create_identity();
    assert(fabs(lag_get_dc_gain(&lag) - 1.0) < EPS);
    assert(fabs(lag_get_beta(&lag) - 1.0) < EPS);
    PASS();
}

/* ==========================================================================
 * L2: Steady-state error improvement
 * ========================================================================== */

static void test_ess_improvement(void) {
    TEST("lag_ess_improvement");
    LagCompensator lag = lag_create(5.0, 0.5, 5.0);
    double improvement = lag_ess_improvement(&lag, LAG_ESS_STEP);
    assert(fabs(improvement - 5.0) < EPS);
    /* For ramp and parabolic, the improvement factor is also Kc */
    assert(fabs(lag_ess_improvement(&lag, LAG_ESS_RAMP) - 5.0) < EPS);
    PASS();
}

/* ==========================================================================
 * L3: s-domain evaluation
 * ========================================================================== */

static void test_eval_s_dc(void) {
    TEST("lag_eval_s at s=0 (DC gain)");
    LagCompensator lag = lag_create(4.0, 1.0, 2.0);
    LagComplex s0 = {0.0, 0.0};
    LagComplex result = lag_eval_s(&lag, s0);
    /* G_c(0) = Kc = 4 */
    assert(fabs(result.re - 4.0) < EPS);
    assert(fabs(result.im - 0.0) < EPS);
    PASS();
}

static void test_eval_magnitude(void) {
    TEST("lag_eval_magnitude at various frequencies");
    LagCompensator lag = lag_create(10.0, 1.0, 10.0);

    /* At w=0: magnitude = Kc = 10 */
    double mag0 = lag_eval_magnitude(&lag, 0.0);
    assert(fabs(mag0 - 10.0) < EPS);

    /* At w->inf: magnitude -> Kc/beta = 1 */
    double mag_inf = lag_eval_magnitude(&lag, 1e6);
    assert(fabs(mag_inf - 1.0) < 0.01);  /* ~1% tolerance at high freq */

    /* At w = w_max = 1/(T*sqrt(beta)): known phase lag */
    double w_max = lag_get_max_lag_frequency(&lag);
    double mag_at_wmax = lag_eval_magnitude(&lag, w_max);
    /* |G_c(jw_max)| = Kc * sqrt(1+1/beta) / sqrt(1+beta) = Kc/sqrt(beta) */
    double expected_mag = 10.0 / sqrt(10.0);  /* = 3.1623 */
    assert(fabs(mag_at_wmax - expected_mag) < EPS);
    PASS();
}

static void test_eval_phase(void) {
    TEST("lag_eval_phase — phase is always negative for w>0");
    LagCompensator lag = lag_create(5.0, 0.5, 3.0);
    /* Phase should be negative for positive frequencies (lag) */
    double phi1 = lag_eval_phase(&lag, 0.1);
    double phi2 = lag_eval_phase(&lag, 1.0);
    double phi3 = lag_eval_phase(&lag, 10.0);
    assert(phi1 < 0.0);
    assert(phi2 < 0.0);
    assert(phi3 < 0.0);
    PASS();
}

static void test_frequency_response(void) {
    TEST("lag_eval_frequency full response");
    LagCompensator lag = lag_create(1.0, 1.0, 2.0);
    LagFreqPoint fp = lag_eval_frequency(&lag, 1.0);
    /* Verify magnitude and phase are consistent */
    double mag_check = sqrt(fp.real_part*fp.real_part + fp.imag_part*fp.imag_part);
    assert(fabs(fp.magnitude - mag_check) < EPS);
    assert(fabs(fp.phase_rad - atan2(fp.imag_part, fp.real_part)) < EPS);
    PASS();
}

/* ==========================================================================
 * L4: Phase lag theorem — maximum phase lag
 * ========================================================================== */

static void test_max_phase_lag(void) {
    TEST("Maximum phase lag theorem verification");
    /* Theorem: max |phase| occurs at w_max = 1/(T*sqrt(beta))
     * and phi_max = arcsin((1-beta)/(1+beta)) */
    LagCompensator lag = lag_create(1.0, 1.0, 4.0);
    double w_max = lag_get_max_lag_frequency(&lag);
    double expected_wmax = 1.0 / (1.0 * sqrt(4.0));  /* 0.5 */
    assert(fabs(w_max - expected_wmax) < EPS);

    double phi_max = lag_get_max_phase_lag(&lag);
    double expected_phi = asin((1.0 - 4.0) / (1.0 + 4.0));  /* arcsin(-3/5) */
    /* phi_max is stored as magnitude (positive) */
    assert(fabs(phi_max - fabs(expected_phi)) < EPS);
    PASS();
}

/* ==========================================================================
 * L4: Corner frequencies
 * ========================================================================== */

static void test_corner_frequencies(void) {
    TEST("Corner frequencies: w_low = 1/(beta*T), w_high = 1/T");
    LagCompensator lag = lag_create(1.0, 2.0, 5.0);
    double w_low, w_high;
    lag_get_corner_frequencies(&lag, &w_low, &w_high);
    assert(fabs(w_low - 0.1) < EPS);   /* 1/(5*2) = 0.1 */
    assert(fabs(w_high - 0.5) < EPS);  /* 1/2 = 0.5 */
    PASS();
}

/* ==========================================================================
 * L4: Stability checks
 * ========================================================================== */

static void test_stability(void) {
    TEST("lag compensator is always stable and minimum-phase for valid params");
    LagCompensator lag = lag_create(2.0, 0.3, 3.0);
    assert(lag_is_stable(&lag) == 1);
    assert(lag_is_minimum_phase(&lag) == 1);
    PASS();
}

/* ==========================================================================
 * L5: Validation
 * ========================================================================== */

static void test_validate_valid(void) {
    TEST("lag_validate — valid compensator");
    LagCompensator lag = lag_create(1.0, 1.0, 2.0);
    assert(lag_validate(&lag) == 0);
    PASS();
}

static void test_validate_invalid_Kc(void) {
    TEST("lag_validate — invalid Kc (negative)");
    LagCompensator lag = lag_create(-1.0, 1.0, 2.0);
    assert(lag_validate(&lag) == -1);
    PASS();
}

static void test_validate_invalid_beta(void) {
    TEST("lag_validate — invalid beta (must be > 1)");
    LagCompensator lag = lag_create(1.0, 1.0, 0.5);
    assert(lag_validate(&lag) == -3);
    PASS();
}

static void test_validate_invalid_T(void) {
    TEST("lag_validate — invalid T (negative)");
    LagCompensator lag = lag_create(1.0, -1.0, 2.0);
    assert(lag_validate(&lag) == -2);
    PASS();
}

/* ==========================================================================
 * L5: Transfer function conversion
 * ========================================================================== */

static void test_to_transfer_function(void) {
    TEST("lag_to_transfer_function — verify coefficients");
    LagCompensator lag = lag_create(3.0, 2.0, 4.0);
    LagTransferFunction tf = lag_to_transfer_function(&lag);
    assert(tf.numerator.order == 1);
    assert(tf.denominator.order == 1);
    /* num: [3, 6], den: [1, 8] */
    assert(fabs(tf.numerator.coeff[0] - 3.0) < EPS);
    assert(fabs(tf.numerator.coeff[1] - 6.0) < EPS);
    assert(fabs(tf.denominator.coeff[0] - 1.0) < EPS);
    assert(fabs(tf.denominator.coeff[1] - 8.0) < EPS);
    assert(fabs(tf.dc_gain - 3.0) < EPS);
    free(tf.numerator.coeff);
    free(tf.denominator.coeff);
    PASS();
}

/* ==========================================================================
 * L5: Design from error constants
 * ========================================================================== */

static void test_design_from_error_constants(void) {
    TEST("lag_design_from_error_constants");
    LagCompensator result;
    int ret = lag_design_from_error_constants(1.0, 5.0, 2.0, &result);
    assert(ret == 0);
    assert(lag_validate(&result) == 0);
    /* beta should be approx 5, Kc = beta */
    assert(fabs(lag_get_beta(&result) - 5.0) < 1.0);
    PASS();
}

/* ==========================================================================
 * L5: Design spec defaults
 * ========================================================================== */

static void test_design_spec_default(void) {
    TEST("lag_design_spec_default");
    LagDesignSpec spec = lag_design_spec_default();
    assert(fabs(spec.phase_margin_target - 45.0) < EPS);
    assert(fabs(spec.gain_margin_target - 10.0) < EPS);
    assert(spec.ess_type == LAG_ESS_STEP);
    PASS();
}

/* ==========================================================================
 * L5: Safety margin heuristic
 * ========================================================================== */

static void test_safety_margin(void) {
    TEST("lag_recommended_safety_margin — within bounds [5,15]");
    double s1 = lag_recommended_safety_margin(2.0);
    double s2 = lag_recommended_safety_margin(10.0);
    double s3 = lag_recommended_safety_margin(1.0);
    assert(s1 >= 5.0 && s1 <= 15.0);
    assert(s2 >= 5.0 && s2 <= 15.0);
    assert(s3 >= 5.0 && s3 <= 15.0);
    PASS();
}

/* ==========================================================================
 * L5: Asymptotes
 * ========================================================================== */

static void test_asymptotes(void) {
    TEST("Low and high frequency asymptotes in dB");
    LagCompensator lag = lag_create(10.0, 1.0, 5.0);
    double low_db = lag_low_freq_asymptote_db(&lag);
    double high_db = lag_high_freq_asymptote_db(&lag);
    assert(fabs(low_db - 20.0) < 0.01);   /* 20*log10(10) = 20 dB */
    assert(fabs(high_db - 20.0*log10(2.0)) < 0.01);  /* 20*log10(2) */
    PASS();
}

/* ==========================================================================
 * L5: String and comparison
 * ========================================================================== */

static void test_to_string(void) {
    TEST("lag_to_string — non-null output");
    LagCompensator lag = lag_create(2.0, 0.5, 3.0);
    char buf[256];
    lag_to_string(&lag, buf, sizeof(buf));
    assert(strlen(buf) > 0);
    PASS();
}

static void test_equals(void) {
    TEST("lag_equals — equality and inequality");
    LagCompensator a = lag_create(1.0, 1.0, 2.0);
    LagCompensator b = lag_create(1.0, 1.0, 2.0);
    LagCompensator c = lag_create(1.0, 1.0, 3.0);
    assert(lag_equals(&a, &b, 1e-9) == 1);
    assert(lag_equals(&a, &c, 1e-9) == 0);
    PASS();
}

/* ==========================================================================
 * L5: Bode computation — basic
 * ========================================================================== */

static void test_compute_bode_basic(void) {
    TEST("lag_compute_bode — first-order TF");
    /* Create a simple first-order TF: G(s) = 2/(s+1) */
    LagTransferFunction tf;
    tf.numerator.order = 0;
    tf.numerator.coeff = (double*)malloc(1 * sizeof(double));
    tf.numerator.coeff[0] = 2.0;
    tf.denominator.order = 1;
    tf.denominator.coeff = (double*)malloc(2 * sizeof(double));
    tf.denominator.coeff[0] = 1.0;
    tf.denominator.coeff[1] = 1.0;
    tf.dc_gain = 2.0;

    LagBodeData bode;
    int ret = lag_compute_bode(&tf, 100, &bode);
    assert(ret == 0);
    assert(bode.num_points == 100);
    assert(bode.points != NULL);
    /* DC gain in dB: 20*log10(2) ~= 6.02 */
    assert(fabs(bode.dc_gain_db - 20.0*log10(2.0)) < 1.0);

    free(bode.points);
    free(tf.numerator.coeff);
    free(tf.denominator.coeff);
    PASS();
}

/* ==========================================================================
 * L4: Gain crossover and phase margin (mathematical assertions)
 * ========================================================================== */

static void test_gain_crossover(void) {
    TEST("lag_find_gain_crossover — detects 0dB crossing");
    /* Create Bode data with known crossover at w=2 */
    LagBodeData bode;
    bode.num_points = 3;
    bode.points = (LagFreqPoint*)malloc(3 * sizeof(LagFreqPoint));
    bode.points[0].omega = 0.1;
    bode.points[0].magnitude_db = 20.0;
    bode.points[1].omega = 1.0;
    bode.points[1].magnitude_db = 6.0;
    bode.points[2].omega = 10.0;
    bode.points[2].magnitude_db = -10.0;

    double w_gc = lag_find_gain_crossover(&bode);
    assert(w_gc > 1.0 && w_gc < 10.0);  /* crossover between points */

    free(bode.points);
    PASS();
}

static void test_phase_margin(void) {
    TEST("lag_compute_phase_margin — basic case");
    LagBodeData bode;
    bode.num_points = 3;
    bode.points = (LagFreqPoint*)malloc(3 * sizeof(LagFreqPoint));
    /* Crossover at w=1 with magnitude crossing 0dB */
    bode.points[0].omega = 0.5;
    bode.points[0].magnitude_db = 6.0;
    bode.points[0].phase_deg = -90.0;
    bode.points[1].omega = 1.0;
    bode.points[1].magnitude_db = -6.0;
    bode.points[1].phase_deg = -120.0;
    bode.points[2].omega = 2.0;
    bode.points[2].magnitude_db = -12.0;
    bode.points[2].phase_deg = -135.0;

    double pm = lag_compute_phase_margin(&bode);
    /* PM = 180 + phase_at_gc. Phase at gc ~= -105 (interpolated) */
    assert(pm > 0.0);  /* should be stable */

    free(bode.points);
    PASS();
}

/* ==========================================================================
 * L5: First-order identification
 * ========================================================================== */

static void test_identify_first_order(void) {
    TEST("lag_identify_first_order — from known step response");
    LagStepResponse resp;
    int n = 100;
    resp.num_points = n;
    resp.time = (double*)malloc((size_t)n * sizeof(double));
    resp.output = (double*)malloc((size_t)n * sizeof(double));
    resp.error = NULL;
    resp.control_signal = NULL;
    resp.final_value = 5.0;

    /* Generate first-order step response: y(t) = 5*(1-exp(-t/0.5)) */
    double K_actual = 5.0, tau_actual = 0.5;
    for (int i = 0; i < n; i++) {
        resp.time[i] = i * 0.02;
        resp.output[i] = K_actual * (1.0 - exp(-resp.time[i] / tau_actual));
    }

    double K_id, tau_id;
    int ret = lag_identify_first_order(&resp, 1.0, &K_id, &tau_id);
    assert(ret == 0);
    assert(fabs(K_id - K_actual) < 0.1);    /* 10% tolerance */
    assert(fabs(tau_id - tau_actual) < 0.1);

    free(resp.time);
    free(resp.output);
    PASS();
}

/* ==========================================================================
 * L5: DC gain identification
 * ========================================================================== */

static void test_identify_dc_gain(void) {
    TEST("lag_identify_dc_gain");
    LagStepResponse resp;
    resp.final_value = 3.5;
    resp.num_points = 10;
    double K = lag_identify_dc_gain(&resp, 1.0);
    assert(fabs(K - 3.5) < EPS);
    PASS();
}

/* ==========================================================================
 * L7: DC motor design (application)
 * ========================================================================== */

static void test_dc_motor_design(void) {
    TEST("lag_design_dc_motor_speed");
    LagDCMotorParams params;
    params.R = 1.0;
    params.L = 0.01;
    params.Kb = 0.1;
    params.Kt = 0.1;
    params.J = 0.01;
    params.B = 0.001;
    params.rated_speed = 300.0;
    params.rated_torque = 1.0;
    params.supply_voltage = 24.0;

    LagCompensator comp;
    int ret = lag_design_dc_motor_speed(&params, 0.02, &comp);
    assert(ret == 0);
    assert(lag_validate(&comp) == 0);
    PASS();
}

/* ==========================================================================
 * L5: Simulation — step response
 * ========================================================================== */

static void test_simulate_step(void) {
    TEST("lag_simulate_step_first_order — basic simulation");
    LagCompensator lag = lag_create(5.0, 0.2, 5.0);
    LagStepResponse resp;
    int ret = lag_simulate_step_first_order(&lag, 2.0, 0.5, 5.0, 200, &resp);
    assert(ret == 0);
    assert(resp.num_points == 200);
    assert(resp.time != NULL);
    assert(resp.output != NULL);
    assert(resp.final_value > 0.0);
    /* Steady-state error should be small with lag compensation */
    assert(fabs(resp.steady_state_error) < 0.2);

    lag_step_response_free(&resp);
    PASS();
}

/* ==========================================================================
 * L4: Nyquist stability — encirclement detection
 * ========================================================================== */

static void test_nyquist_stable(void) {
    TEST("lag_compute_nyquist — stable open-loop system");
    LagTransferFunction tf;
    tf.numerator.order = 0;
    tf.numerator.coeff = (double*)malloc(1 * sizeof(double));
    tf.numerator.coeff[0] = 1.0;
    tf.denominator.order = 1;
    tf.denominator.coeff = (double*)malloc(2 * sizeof(double));
    tf.denominator.coeff[0] = 1.0;
    tf.denominator.coeff[1] = 1.0;
    tf.dc_gain = 1.0;

    LagNyquistData nyquist;
    int ret = lag_compute_nyquist(&tf, 100, 0, &nyquist);
    assert(ret == 0);
    /* With P=0 and no encirclements, Z=0 => stable */
    assert(nyquist.is_stable == 1);

    free(nyquist.points);
    free(tf.numerator.coeff);
    free(tf.denominator.coeff);
    PASS();
}

/* ==========================================================================
 * Test runner
 * ========================================================================== */

int main(void) {
    printf("=== mini-lag-compensator Test Suite ===\n\n");

    /* L1: Construction */
    test_create();
    test_create_pole_zero();
    test_create_corners();
    test_identity();

    /* L2: Core concepts */
    test_ess_improvement();

    /* L3: s-domain evaluation */
    test_eval_s_dc();
    test_eval_magnitude();
    test_eval_phase();
    test_frequency_response();

    /* L4: Theorems */
    test_max_phase_lag();
    test_corner_frequencies();
    test_stability();
    test_gain_crossover();
    test_phase_margin();
    test_nyquist_stable();

    /* L5: Algorithms */
    test_validate_valid();
    test_validate_invalid_Kc();
    test_validate_invalid_beta();
    test_validate_invalid_T();
    test_to_transfer_function();
    test_design_from_error_constants();
    test_design_spec_default();
    test_safety_margin();
    test_asymptotes();
    test_to_string();
    test_equals();
    test_compute_bode_basic();
    test_identify_first_order();
    test_identify_dc_gain();

    /* L6-L7: Applications */
    test_dc_motor_design();
    test_simulate_step();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}