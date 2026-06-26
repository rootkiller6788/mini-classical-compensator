/**
 * @file feedforward_design.c
 * @brief Feedforward controller design implementations.
 *
 * L4-L5: Model inverse computation, ZPETC, Diophantine design,
 * properness enforcement, reference/disturbance FF design.
 */

#include "feedforward_design.h"
#include "feedforward_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * L4: Model Inverse Computation
 * ========================================================================== */

int ff_model_inverse(const TransferFn *plant, TransferFn *inv,
                     int filter_order, double filter_bw)
{
    if (!plant || !inv) return -1;

    /* Check minimum-phase condition */
    if (!tf_is_minimum_phase(plant, 0)) {
        /* Non-minimum-phase ? exact inverse would be unstable.
         * Return error; caller should use ZPETC or approximate inverse. */
        return -1;
    }

    /* P^{-1}(s) = den_P(s) / (gain_P * num_P(s))
     * Swap numerator and denominator */
    inv->num = poly_create(plant->den.coeff, plant->den.degree);
    inv->den = poly_create(plant->num.coeff, plant->num.degree);
    inv->gain = 1.0 / plant->gain;

    /* If plant is strictly proper (deg(num) < deg(den)), the inverse
     * is improper. Add low-pass filter to enforce properness. */
    if (plant->den.degree > plant->num.degree) {
        int rel_deg = plant->den.degree - plant->num.degree;
        int needed_order = (filter_order > 0) ? filter_order : rel_deg;

        /* Add (filter_bw/(s + filter_bw))^needed_order */
        for (int i = 0; i < needed_order; i++) {
            double num_lp[] = {filter_bw};
            double den_lp[] = {filter_bw, 1.0};
            Poly num_p = poly_create(num_lp, 0);
            Poly den_p = poly_create(den_lp, 1);
            Poly new_num = poly_mul(&inv->num, &num_p);
            Poly new_den = poly_mul(&inv->den, &den_p);
            poly_free(&inv->num);
            poly_free(&inv->den);
            poly_free(&num_p);
            poly_free(&den_p);
            inv->num = new_num;
            inv->den = new_den;
        }
        /* Note: filter_bw is already embedded in the numerator
         * polynomial via multiplication. DC gain remains 1.0. */
    }

    return 0;
}

int ff_causal_inverse(const TransferFn *plant, TransferFn *inv,
                      int max_order)
{
    if (!plant || !inv) return -1;

    /* For strictly proper plant: B(s)/A(s) with deg(A) > deg(B).
     * We want A(s)/B(s) but it is improper.
     *
     * Use polynomial long division:
     *   A(s) = Q(s)*B(s) + R(s)
     *   A(s)/B(s) = Q(s) + R(s)/B(s)
     *
     * Truncate at max_order. */
    int degA = plant->den.degree;
    int degB = plant->num.degree;

    if (degB == 0 && fabs(plant->num.coeff[0]) < 1e-15) return -1;

    /* Copy polynomials for division */
    double *A = (double *)malloc((degA + 1) * sizeof(double));
    double *B = (double *)malloc((degB + 1) * sizeof(double));
    if (!A || !B) { free(A); free(B); return -1; }

    for (int i = 0; i <= degA; i++) A[i] = plant->den.coeff[i] / plant->gain;
    for (int i = 0; i <= degB; i++) B[i] = plant->num.coeff[i];

    /* Normalize B to monic */
    double B_lead = B[degB];
    if (fabs(B_lead) < 1e-15) { free(A); free(B); return -1; }
    for (int i = 0; i <= degB; i++) B[i] /= B_lead;

    int q_deg = degA - degB;
    if (q_deg < 0) q_deg = 0;
    if (q_deg > max_order) q_deg = max_order;
    double *Q = (double *)calloc(q_deg + 1, sizeof(double));
    double *R = (double *)malloc((degA + 1) * sizeof(double));
    if (!Q || !R) { free(A); free(B); free(Q); free(R); return -1; }
    memcpy(R, A, (degA + 1) * sizeof(double));

    /* Long division */
    for (int k = degA; k >= degB; k--) {
        int q_idx = k - degB;
        if (q_idx > q_deg) continue;
        double factor = R[k] / B[degB];
        Q[q_idx] = factor;
        for (int j = 0; j <= degB; j++) {
            R[k - degB + j] -= factor * B[j];
        }
    }

    /* Construct approximate inverse: Q(s) (causal part) */
    inv->num = poly_create(Q, q_deg);
    double one[] = {1.0};
    inv->den = poly_create(one, 0);
    inv->gain = 1.0;

    free(A); free(B); free(Q); free(R);
    return 0;
}

int ff_zpetc_design(const TransferFn *plant_dt, TransferFn *zp_ff,
                    double Ts)
{
    if (!plant_dt || !zp_ff) return -1;
    (void)Ts;

    /* ZPETC (Tomizuka 1987):
     * Factor P(z) = z^{-d} * B^+(z) * B^-(z) / A(z)
     * where B^+ has cancellable (stable) zeros
     * and B^- has uncancellable (unstable or lightly damped) zeros.
     *
     * F(z) = z^{d} * A(z) * B^-(z^{-1}) / [B^+(z) * (B^-(1))^2]
     *
     * For this implementation, we do a simplified version assuming
     * the plant can be factored. */

    /* Step 1: compute A(z) * B^-(z^{-1})
     * B^-(z) = b0 + b1*z + ... + bm*z^m
     * B^-(z^{-1}) = b0 + b1*z^{-1} + ... + bm*z^{-m}
     * = z^{-m} * (b0*z^m + b1*z^{m-1} + ... + bm)
     * = z^{-m} * B_rev(z) where B_rev reverses coefficients */

    int m = plant_dt->num.degree;
    double *B_rev = (double *)malloc((m + 1) * sizeof(double));
    if (!B_rev) return -1;
    for (int i = 0; i <= m; i++) {
        B_rev[i] = plant_dt->num.coeff[m - i];
    }

    /* Numerator: A(z) * B_rev(z) */
    Poly num_A = poly_create(plant_dt->den.coeff, plant_dt->den.degree);
    Poly num_Brev = poly_create(B_rev, m);
    Poly num_prod = poly_mul(&num_A, &num_Brev);

    /* Denominator: B^+(z) * (B^-(1))^2
     * For simplicity, treat all zeros as in B^+ if stable, else in B^-
     * B^-(1) = sum of B^- coefficients */
    double B_minus_at_1 = 0.0;
    for (int i = 0; i <= m; i++) {
        B_minus_at_1 += plant_dt->num.coeff[i];
    }
    double B_minus_sq = B_minus_at_1 * B_minus_at_1;
    if (fabs(B_minus_sq) < 1e-15) B_minus_sq = 1.0;

    /* Denominator: polynomial from B^+ and (B^-(1))^2 */
    Poly den_Bplus = poly_create(plant_dt->num.coeff, m);
    double scale = B_minus_sq / plant_dt->gain;
    for (int i = 0; i <= den_Bplus.degree; i++) {
        den_Bplus.coeff[i] *= scale;
    }

    zp_ff->num = num_prod;
    zp_ff->den = den_Bplus;
    zp_ff->gain = 1.0;

    free(B_rev);
    return 0;
}

int ff_diophantine_design(const TransferFn *plant,
                          const TransferFn *desired,
                          Poly *ff_num, Poly *ff_den)
{
    if (!plant || !desired || !ff_num || !ff_den) return -1;

    /* Model-following 2-DOF design:
     * Solve A*R + B*S = A_o for R and S (feedback part).
     * Then T = A_o * B_m / B (feedforward numerator).
     *
     * In our simplified version:
     *   ff_num = A * B_m (feedforward numerator)
     *   ff_den = B * A_o (feedforward denominator, simplified)
     *
     * This gives F(s) = T(s)/R(s) for the 2-DOF structure. */

    /* T = plant->den * desired->num (numerator of FF) */
    *ff_num = poly_mul(&plant->den, &desired->num);

    /* R is typically solved from the Diophantine equation.
     * Here we provide a simplified R = plant->num (cancels plant dynamics) */
    *ff_den = poly_create(plant->num.coeff, plant->num.degree);

    return 0;
}

int ff_enforce_properness(TransferFn *tf, double filter_bw)
{
    if (!tf) return -1;

    int rel_deg = tf->num.degree - tf->den.degree;
    if (rel_deg <= 0) return 0; /* Already proper */

    /* Add (rel_deg) poles at s = -filter_bw */
    double num_lp[] = {filter_bw};
    double den_lp[] = {filter_bw, 1.0};

    for (int i = 0; i < rel_deg; i++) {
        Poly np = poly_create(num_lp, 0);
        Poly dp = poly_create(den_lp, 1);
        Poly new_num = poly_mul(&tf->num, &np);
        Poly new_den = poly_mul(&tf->den, &dp);
        poly_free(&tf->num);
        poly_free(&tf->den);
        poly_free(&np);
        poly_free(&dp);
        tf->num = new_num;
        tf->den = new_den;
        tf->gain *= filter_bw;
    }

    return 0;
}

/* ==========================================================================
 * L5: Reference Feedforward Design
 * ========================================================================== */

int ff_design_ref_feedforward(const TransferFn *plant,
                              const TransferFn *fb_ctrl,
                              const TransferFn *desired,
                              TransferFn *ff_ref)
{
    if (!plant || !fb_ctrl || !desired || !ff_ref) return -1;

    /* F(s) = [T_d(s)*(1 + P(s)*C(s)) - P(s)*C(s)] / P(s)
     *     = T_d(s)*(1/(P(s)) + C(s)) - C(s)
     *
     * This design places the closed-loop poles at desired locations
     * while the feedforward shapes the zeros for tracking. */

    /* Compute 1/(P*C) related terms */
    TransferFn pc = tf_series(plant, fb_ctrl);

    /* Denominator: den_P * den_C, Numerator: gain_P*gain_C * num_P * num_C */
    /* 1 + PC */
    Poly den_pc = poly_mul(&plant->den, &fb_ctrl->den);
    Poly num_pc_prod = poly_mul(&plant->num, &fb_ctrl->num);
    double k_pc = plant->gain * fb_ctrl->gain;
    for (int i = 0; i <= num_pc_prod.degree; i++) {
        num_pc_prod.coeff[i] *= k_pc;
    }
    Poly one_plus_pc_den = poly_add(&den_pc, &num_pc_prod);

    /* Numerator: T_d * (1 + PC) = num_T * (den_PC + num_PC_prod)
     * as polynomial; this is a transfer function multiplication */
    Poly num_t = poly_create(desired->num.coeff, desired->num.degree);
    Poly num_t_s = poly_mul(&num_t, &one_plus_pc_den);
    /* Scale */
    for (int i = 0; i <= num_t_s.degree; i++) num_t_s.coeff[i] *= desired->gain;

    /* Subtract P*C to get the FF numerator:
     * num_ff = T_d*(1+PC) - PC = (T_d*(1+PC) numerator)*(den_ff) - PC */
    Poly num_pc = poly_mul(&plant->num, &fb_ctrl->num);
    for (int i = 0; i <= num_pc.degree; i++) num_pc.coeff[i] *= k_pc;

    /* We need common denominator. FF = (T_d*(1+PC) - PC) / P
     * = [T_d*(den_PC+num_PC_prod) - num_PC*den_T] / (den_T * ...) / P
     * This gets complex; simplified implementation: */
    Poly num_ff_temp = poly_sub(&num_t_s, &num_pc);

    /* Divide by P(s): multiply numerator by den_P */
    Poly num_ff_final = poly_mul(&num_ff_temp, &plant->den);

    /* Denominator of FF: den_T * num_P */
    Poly den_ff = poly_mul(&desired->den, &plant->num);

    ff_ref->num = num_ff_final;
    ff_ref->den = den_ff;
    ff_ref->gain = 1.0 / plant->gain;

    /* Cleanup */
    poly_free(&num_t); poly_free(&num_t_s);
    poly_free(&num_pc); poly_free(&num_ff_temp);
    poly_free(&num_pc_prod);
    poly_free(&den_pc); poly_free(&one_plus_pc_den);
    tf_free(&pc);

    return 0;
}

void ff_design_prefilter_1st(double rise_time, Prefilter *prefilter)
{
    /* F(s) = 1 / (tau*s + 1)
     * Rise time (10%-90% for first order): tr = 2.2 * tau
     * So tau = tr / 2.2 */
    double tau = rise_time / 2.2;
    if (tau < 1e-9) tau = 1e-9;

    double num[] = {1.0};
    double den[] = {1.0, tau};

    prefilter->filter = tf_create(num, 0, den, 1, 1.0);
    prefilter->rise_time = rise_time;
    prefilter->bandwidth = 1.0 / tau;
    prefilter->order = 1;
    prefilter->damping = 0.0;
}

void ff_design_prefilter_2nd(double rise_time, double damping,
                             Prefilter *prefilter)
{
    /* F(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
     * For second order: tr ? 1.8 / wn (approx for zeta = 0.7)
     * More generally: wn ? 2.2 / (tr * sqrt(1 - zeta^2)) */
    double zeta = damping;
    if (zeta >= 1.0) zeta = 0.999;
    double wn = 2.2 / (rise_time * sqrt(1.0 - zeta * zeta));
    if (wn < 0.01) wn = 0.01;

    double num[] = {wn * wn};
    double den[] = {wn * wn, 2.0 * zeta * wn, 1.0};

    prefilter->filter = tf_create(num, 0, den, 2, 1.0);
    prefilter->rise_time = rise_time;
    prefilter->bandwidth = wn;
    prefilter->order = 2;
    prefilter->damping = zeta;
}

/* ==========================================================================
 * L5: Disturbance Feedforward Design
 * ========================================================================== */

int ff_design_static_dist_ff(const TransferFn *plant,
                             const TransferFn *dist_model,
                             double *ff_gain)
{
    if (!plant || !dist_model || !ff_gain) return -1;

    /* K_d = G_d(0) / P(0)
     * Evaluate DC gains */
    double Gd_dc = ff_dc_gain(dist_model, 0);
    double P_dc = ff_dc_gain(plant, 0);

    if (fabs(P_dc) < 1e-15) {
        /* Integrating plant ? P(0) = infinity */
        return -1;
    }

    *ff_gain = -Gd_dc / P_dc;
    return 0;
}

int ff_design_dynamic_dist_ff(const TransferFn *plant,
                              const TransferFn *dist_model,
                              TransferFn *comp,
                              double filter_bw)
{
    if (!plant || !dist_model || !comp) return -1;

    /* D_ff(s) = -G_d(s) / P(s) */
    int ret = ff_compute_dist_compensator(plant, dist_model, comp);
    if (ret != 0) return ret;

    /* Enforce properness */
    if (comp->num.degree > comp->den.degree) {
        ff_enforce_properness(comp, filter_bw);
    }

    return 0;
}

int ff_eval_dist_rejection(const TransferFn *plant,
                           const TransferFn *fb_ctrl,
                           const TransferFn *dist_ff,
                           const TransferFn *dist_model,
                           double omega, double *attenuation)
{
    if (!plant || !fb_ctrl || !dist_ff || !dist_model || !attenuation)
        return -1;

    /* Sensitivity: S(s) = 1 / (1 + P(s)*C(s)) */
    TransferFn pc = tf_series(plant, fb_ctrl);
    double mag_pc, phase_pc;
    tf_freq_response(&pc, omega, &mag_pc, &phase_pc);
    tf_free(&pc);

    /* |S(jw)| = 1 / |1 + PC(jw)| */
    /* |1 + PC| = sqrt((1+Re(PC))^2 + Im(PC)^2) */
    double re_pc = mag_pc * cos(phase_pc);
    double im_pc = mag_pc * sin(phase_pc);
    double mag_one_plus_pc = sqrt((1.0 + re_pc) * (1.0 + re_pc) + im_pc * im_pc);
    double mag_S = 1.0 / (mag_one_plus_pc + 1e-15);

    /* With FF: y_ff/d = P*S + P*D_ff*G_d*S (simplified)
     * Without FF: y/d = P*S
     * Attenuation = |P*S + P*D_ff*G_d*S| / |P*S|
     *            = |1 + D_ff*G_d| */
    double mag_plant, ph_plant;
    tf_freq_response(plant, omega, &mag_plant, &ph_plant);
    /* D_ff * G_d: series connection */
    TransferFn dff_gd = tf_series(dist_ff, dist_model);
    double mag_dg, ph_dg;
    tf_freq_response(&dff_gd, omega, &mag_dg, &ph_dg);
    tf_free(&dff_gd);

    double re_dg = mag_dg * cos(ph_dg);
    double im_dg = mag_dg * sin(ph_dg);
    double mag_one_plus_dg = sqrt((1.0 + re_dg) * (1.0 + re_dg) + im_dg * im_dg);

    /* With disturbance feedforward, disturbance effect is scaled by
     * |1 + D_ff(jw)*G_d(jw)| factor. Ideally = 0 (perfect cancellation). */
    (void)mag_S; (void)mag_plant;
    *attenuation = mag_one_plus_dg;

    return 0;
}
