/**
 * @file feedforward_core.c
 * @brief Core feedforward control implementations.
 *
 * Implements: polynomial arithmetic, transfer function operations,
 * feedforward controller initialization, frequency response,
 * step response computation, pole finding, state-space conversion.
 *
 * L1-L3: Core data structures and mathematical operations
 * for feedforward control analysis and synthesis.
 */

#include "feedforward_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L2: Polynomial Operations (Horner, arithmetic, derivative)
 * ========================================================================== */

Poly poly_create(const double *coeff, int degree)
{
    Poly p;
    p.degree = degree;
    p.coeff = (double *)malloc((degree + 1) * sizeof(double));
    if (p.coeff) {
        memcpy(p.coeff, coeff, (degree + 1) * sizeof(double));
    } else {
        p.degree = -1;
    }
    return p;
}

void poly_free(Poly *p)
{
    if (p && p->coeff) {
        free(p->coeff);
        p->coeff = NULL;
        p->degree = -1;
    }
}

double poly_eval(const Poly *p, double x)
{
    if (!p || !p->coeff || p->degree < 0) return 0.0;
    /* Horner's method: O(n), numerically stable */
    double result = p->coeff[p->degree];
    for (int i = p->degree - 1; i >= 0; i--) {
        result = result * x + p->coeff[i];
    }
    return result;
}

Poly poly_mul(const Poly *p, const Poly *q)
{
    int r_deg = p->degree + q->degree;
    double *r_coeff = (double *)calloc(r_deg + 1, sizeof(double));
    if (!r_coeff) {
        Poly empty = {-1, NULL};
        return empty;
    }
    for (int i = 0; i <= p->degree; i++) {
        for (int j = 0; j <= q->degree; j++) {
            r_coeff[i + j] += p->coeff[i] * q->coeff[j];
        }
    }
    /* Trim leading zeros */
    while (r_deg > 0 && fabs(r_coeff[r_deg]) < 1e-15) r_deg--;
    Poly result;
    result.degree = r_deg;
    result.coeff = r_coeff;
    return result;
}

Poly poly_add(const Poly *p, const Poly *q)
{
    int max_deg = (p->degree > q->degree) ? p->degree : q->degree;
    double *r_coeff = (double *)calloc(max_deg + 1, sizeof(double));
    if (!r_coeff) {
        Poly empty = {-1, NULL};
        return empty;
    }
    for (int i = 0; i <= p->degree; i++) r_coeff[i] += p->coeff[i];
    for (int i = 0; i <= q->degree; i++) r_coeff[i] += q->coeff[i];
    /* Trim leading zeros */
    while (max_deg > 0 && fabs(r_coeff[max_deg]) < 1e-15) max_deg--;
    Poly result;
    result.degree = max_deg;
    result.coeff = r_coeff;
    return result;
}

Poly poly_sub(const Poly *p, const Poly *q)
{
    int max_deg = (p->degree > q->degree) ? p->degree : q->degree;
    double *r_coeff = (double *)calloc(max_deg + 1, sizeof(double));
    if (!r_coeff) {
        Poly empty = {-1, NULL};
        return empty;
    }
    for (int i = 0; i <= p->degree; i++) r_coeff[i] = p->coeff[i];
    for (int i = 0; i <= q->degree; i++) r_coeff[i] -= q->coeff[i];
    while (max_deg > 0 && fabs(r_coeff[max_deg]) < 1e-15) max_deg--;
    Poly result;
    result.degree = max_deg;
    result.coeff = r_coeff;
    return result;
}

Poly poly_derivative(const Poly *p)
{
    if (p->degree <= 0) {
        double zero = 0.0;
        return poly_create(&zero, 0);
    }
    int r_deg = p->degree - 1;
    double *r_coeff = (double *)malloc((r_deg + 1) * sizeof(double));
    if (!r_coeff) {
        Poly empty = {-1, NULL};
        return empty;
    }
    for (int i = 0; i <= r_deg; i++) {
        r_coeff[i] = p->coeff[i + 1] * (i + 1);
    }
    Poly result;
    result.degree = r_deg;
    result.coeff = r_coeff;
    return result;
}

/* ==========================================================================
 * L2: Transfer Function Creation and Management
 * ========================================================================== */

TransferFn tf_create(const double *num, int num_deg,
                     const double *den, int den_deg,
                     double gain)
{
    TransferFn tf;
    tf.num = poly_create(num, num_deg);
    tf.den = poly_create(den, den_deg);
    tf.gain = gain;
    return tf;
}

void tf_free(TransferFn *tf)
{
    if (tf) {
        poly_free(&tf->num);
        poly_free(&tf->den);
        tf->gain = 1.0;
    }
}

/* ==========================================================================
 * L2: Feedforward Controller Operations
 * ========================================================================== */

int ff_compute_ideal_filter(const TransferFn *plant,
                            const TransferFn *desired,
                            TransferFn *ff)
{
    if (!plant || !desired || !ff) return -1;
    /* F(s) = T_d(s) / P(s)
     * In polynomial form: num_F / den_F = (num_T * den_P) / (den_T * num_P) */
    Poly num_F = poly_mul(&desired->num, &plant->den);
    Poly den_F = poly_mul(&desired->den, &plant->num);
    double gain_F = desired->gain / plant->gain;

    /* Set output transfer function.
     * Note: if num_F.degree > den_F.degree, the feedforward filter
     * is improper (non-causal). Caller should check with ff_is_realizable()
     * and enforce properness with ff_enforce_properness() if needed. */
    ff->num = num_F;
    ff->den = den_F;
    ff->gain = gain_F;
    return 0;
}

int ff_compute_dist_compensator(const TransferFn *plant,
                                const TransferFn *dist,
                                TransferFn *comp)
{
    if (!plant || !dist || !comp) return -1;
    /* D_ff(s) = -G_d(s) / P(s) = -(num_Gd * den_P) / (den_Gd * num_P) */
    Poly num_comp = poly_mul(&dist->num, &plant->den);
    Poly den_comp = poly_mul(&dist->den, &plant->num);
    double gain_comp = -dist->gain / plant->gain;

    comp->num = num_comp;
    comp->den = den_comp;
    comp->gain = gain_comp;

    return 0;
}

int ff_compute_closed_loop_tf(const TransferFn *plant,
                               const TransferFn *fb_ctrl,
                               const TransferFn *ff_ref,
                               TransferFn *cl_tf)
{
    if (!plant || !fb_ctrl || !ff_ref || !cl_tf) return -1;
    /* y/r = (P*F + P*C) / (1 + P*C) */
    TransferFn pf = tf_series(plant, ff_ref);
    TransferFn pc = tf_series(plant, fb_ctrl);
    TransferFn num_sum = tf_parallel(&pf, &pc);
    /* Denominator: 1 + P*C = (denP*denC + numP*numC) / (denP*denC)
     * Actually we need: 1 + P*C = (den_P*den_C + gain_PC*num_P*num_C)/(den_P*den_C) */
    Poly den_pc = poly_mul(&plant->den, &fb_ctrl->den);
    Poly num_pc_prod = poly_mul(&plant->num, &fb_ctrl->num);
    /* Scale num_pc_prod by gain_PC */
    double k_pc = plant->gain * fb_ctrl->gain;
    for (int i = 0; i <= num_pc_prod.degree; i++) {
        num_pc_prod.coeff[i] *= k_pc;
    }
    Poly den_cl_poly = poly_add(&den_pc, &num_pc_prod);
    /* cl_tf = (P*F + P*C) / (1 + P*C)
     * = (num_sum) / (den_sum) * (den_pc) / (den_cl_poly) simplified... */
    /* For simplicity, set cl_tf to the ratio:
     * num_cl = num_sum * den_pc (simplified)
     * den_cl = den_sum * den_cl_poly
     * This is approximate and may need pole-zero cancellation */
    Poly num_cl = poly_mul(&num_sum.num, &den_pc);
    Poly den_cl = poly_mul(&num_sum.den, &den_cl_poly);

    cl_tf->num = num_cl;
    cl_tf->den = den_cl;
    cl_tf->gain = 1.0;

    /* Cleanup temporary allocations */
    tf_free(&pf);
    tf_free(&pc);
    poly_free(&num_sum.num);
    poly_free(&num_sum.den);
    poly_free(&den_pc);
    poly_free(&num_pc_prod);
    poly_free(&den_cl_poly);

    return 0;
}

int ff_is_realizable(const TransferFn *ff)
{
    if (!ff) return 0;
    return ff->num.degree <= ff->den.degree;
}

double ff_dc_gain(const TransferFn *ff, int is_dt)
{
    if (!ff) return 0.0;
    /* Continuous: evaluate at s=0 (DC)
     * Discrete:   evaluate at z=1 (DC) */
    double eval_point = is_dt ? 1.0 : 0.0;
    double num_val = poly_eval(&ff->num, eval_point);
    double den_val = poly_eval(&ff->den, eval_point);
    if (fabs(den_val) < 1e-15) {
        /* Integrating system: DC gain is infinite in theory */
        return 1e15;
    }
    return ff->gain * num_val / den_val;
}

void ff_init(const TransferFn *plant, FeedforwardCtrl *ff)
{
    if (!plant || !ff) return;
    ff->plant_model = *plant;
    /* Initialize with identity feedforward */
    double one_arr[] = {1.0};
    ff->ff_filter = tf_create(one_arr, 0, one_arr, 0, 1.0);
    ff->dc_gain = 1.0;
    ff->is_causal = 1;
    ff->is_stable = tf_is_minimum_phase(plant, 0);
}

void ff_free(FeedforwardCtrl *ff)
{
    if (ff) {
        tf_free(&ff->ff_filter);
        tf_free(&ff->plant_inverse);
        /* plant_model may alias ? don't double-free */
        ff->dc_gain = 0.0;
        ff->is_causal = 0;
        ff->is_stable = 0;
    }
}

/* ==========================================================================
 * L3: Frequency and Time Response
 * ========================================================================== */

void tf_freq_response(const TransferFn *tf, double omega,
                      double *mag, double *phase)
{
    if (!tf || !mag || !phase) return;

    /* Evaluate G(j*omega) = K * num(jw) / den(jw)
     * Separate real and imaginary parts via Horner */
    double num_re = 0.0, num_im = 0.0;
    double den_re = 0.0, den_im = 0.0;

    /* num(jw) */
    for (int i = 0; i <= tf->num.degree; i++) {
        double c = tf->num.coeff[i];
        /* (jw)^i = j^i * w^i
         * i=0: 1, i=1: jw, i=2: -w^2, i=3: -jw^3, ... */
        double w_pow = pow(omega, i);
        switch (i % 4) {
            case 0: num_re += c * w_pow; break;
            case 1: num_im += c * w_pow; break;
            case 2: num_re -= c * w_pow; break;
            case 3: num_im -= c * w_pow; break;
        }
    }

    /* den(jw) */
    for (int i = 0; i <= tf->den.degree; i++) {
        double c = tf->den.coeff[i];
        double w_pow = pow(omega, i);
        switch (i % 4) {
            case 0: den_re += c * w_pow; break;
            case 1: den_im += c * w_pow; break;
            case 2: den_re -= c * w_pow; break;
            case 3: den_im -= c * w_pow; break;
        }
    }

    /* Complex division */
    double den_mag2 = den_re * den_re + den_im * den_im;
    if (den_mag2 < 1e-30) {
        *mag = 1e15;
        *phase = 0.0;
        return;
    }
    double res_re = tf->gain * (num_re * den_re + num_im * den_im) / den_mag2;
    double res_im = tf->gain * (num_im * den_re - num_re * den_im) / den_mag2;

    *mag = sqrt(res_re * res_re + res_im * res_im);
    *phase = atan2(res_im, res_re);
}

double tf_step_response(const TransferFn *tf, double t, int n_approx)
{
    if (!tf || t < 0.0) return 0.0;
    /* Approximate step response via dominant pole expansion.
     * For simple first/second order, exact formula is used.
     * For higher order, use forward Euler ODE simulation. */

    if (t == 0.0) return 0.0;

    /* For a first-order system: K/(tau*s + 1), y(t) = K*(1 - exp(-t/tau)) */
    if (tf->den.degree == 1 && tf->num.degree == 0) {
        double tau = tf->den.coeff[1] / tf->den.coeff[0];
        double K = tf->gain * tf->num.coeff[0] / tf->den.coeff[0];
        return K * (1.0 - exp(-t / tau));
    }

    /* For a second-order system:
     * K*wn^2/(s^2 + 2*zeta*wn*s + wn^2)
     * y(t) = K*[1 - exp(-zeta*wn*t)*(cos(wd*t) + zeta/sqrt(1-zeta^2)*sin(wd*t))] */
    if (tf->den.degree == 2 && tf->num.degree <= 1) {
        double a0 = tf->den.coeff[0];
        double a1 = tf->den.coeff[1];
        double a2 = tf->den.coeff[2];
        double wn = sqrt(a0 / a2);
        double zeta = a1 / (2.0 * a2 * wn);
        double K = tf->gain * tf->num.coeff[0] / a0;

        if (zeta < 1.0) {
            /* Underdamped */
            double wd = wn * sqrt(1.0 - zeta * zeta);
            double decay = exp(-zeta * wn * t);
            double cos_term = cos(wd * t);
            double sin_term = sin(wd * t);
            double damp_factor = zeta / sqrt(1.0 - zeta * zeta);
            return K * (1.0 - decay * (cos_term + damp_factor * sin_term));
        } else if (fabs(zeta - 1.0) < 1e-6) {
            /* Critically damped */
            return K * (1.0 - exp(-wn * t) * (1.0 + wn * t));
        } else {
            /* Overdamped */
            double s1 = -wn * (zeta + sqrt(zeta * zeta - 1.0));
            double s2 = -wn * (zeta - sqrt(zeta * zeta - 1.0));
            return K * (1.0 - (s2 * exp(s1 * t) - s1 * exp(s2 * t)) / (s2 - s1));
        }
    }

    /* For higher-order systems, use coarse numerical integration */
    (void)n_approx;
    double dt = t / 100.0;
    double y = 0.0;
    /* State-space simulation via forward Euler */
    int n = tf->den.degree;
    if (n > 10) n = 10; /* limit */
    double *x = (double *)calloc(n, sizeof(double));
    if (!x) return 0.0;
    double *A = (double *)malloc(n * n * sizeof(double));
    double *B = (double *)malloc(n * sizeof(double));
    double *C = (double *)malloc(n * sizeof(double));
    double D = 0.0;
    if (A && B && C) {
        int dim = tf_to_state_space(tf, A, B, C, &D);
        if (dim > 0) {
            for (double ti = 0.0; ti < t; ti += dt) {
                /* x' = A*x + B*u, u=1 (step) */
                double *x_dot = (double *)calloc(dim, sizeof(double));
                if (!x_dot) break;
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        x_dot[i] += A[i * dim + j] * x[j];
                    }
                    x_dot[i] += B[i]; /* u = 1 */
                }
                for (int i = 0; i < dim; i++) x[i] += dt * x_dot[i];
                free(x_dot);
            }
            y = D; /* D*u */
            for (int i = 0; i < dim; i++) y += C[i] * x[i];
        }
    }
    free(x); free(A); free(B); free(C);
    return y;
}

int tf_find_poles(const TransferFn *tf, double *poles, int n)
{
    if (!tf || !poles) return 0;

    int deg = tf->den.degree;
    /* For degree 1: p = -a0/a1 */
    if (deg == 1 && n >= 1) {
        poles[0] = -tf->den.coeff[0] / tf->den.coeff[1];
        return 1;
    }
    /* For degree 2: quadratic formula */
    if (deg == 2 && n >= 2) {
        double a = tf->den.coeff[2];
        double b = tf->den.coeff[1];
        double c = tf->den.coeff[0];
        double disc = b * b - 4.0 * a * c;
        if (disc >= 0) {
            double sqrt_disc = sqrt(disc);
            poles[0] = (-b - sqrt_disc) / (2.0 * a);
            poles[1] = (-b + sqrt_disc) / (2.0 * a);
        } else {
            /* Complex conjugate ? return real part */
            poles[0] = -b / (2.0 * a);
            poles[1] = -b / (2.0 * a);
        }
        return 2;
    }
    /* For higher degree: try Newton's method with deflation.
     * Start with initial guesses spread along real axis. */
    int found = 0;
    int max_n = (deg < n) ? deg : n;

    /* Aberth-Ehrlich method simplified for real polynomial with real initial guesses */
    double *roots = (double *)calloc(max_n, sizeof(double));
    double *den_copy = (double *)malloc((deg + 1) * sizeof(double));
    if (!roots || !den_copy) {
        free(roots); free(den_copy);
        return 0;
    }
    memcpy(den_copy, tf->den.coeff, (deg + 1) * sizeof(double));

    for (int iter = 0; iter < 100 && found < max_n; iter++) {
        double x_low = -10.0;
        double x_high = 10.0;
        /* Simple deflation + Newton */
        for (int k = 0; k < max_n; k++) {
            double x = x_low + (x_high - x_low) * (k + 0.5) / max_n;
            for (int ni = 0; ni < 50; ni++) {
                /* Horner evaluation with deflation for previously found roots */
                double p_val = 0.0, dp_val = 0.0;
                double px = 1.0;
                for (int j = 0; j <= deg; j++) {
                    p_val += den_copy[j] * px;
                    if (j < deg) dp_val += den_copy[j + 1] * (j + 1) * px;
                    px *= x;
                }
                /* Deflate: divide by (x - root_i) for previously found */
                for (int r = 0; r < found; r++) {
                    double denom = x - roots[r];
                    if (fabs(denom) < 1e-10) { x += 0.1; denom = 0.1; }
                    p_val /= denom;
                }
                /* Newton step */
                if (fabs(dp_val) < 1e-15) break;
                double dx = -p_val / dp_val;
                x += dx;
                if (fabs(dx) < 1e-10) break;
            }
            roots[found] = x;
            /* Verify it's a real root by substitution */
            double check = 0.0, xp = 1.0;
            for (int j = 0; j <= deg; j++) {
                check += den_copy[j] * xp;
                xp *= x;
            }
            if (fabs(check) < 1e-6) {
                poles[found] = x;
                found++;
            }
        }
        if (iter > 10) break; /* Convergence check */
    }

    free(roots);
    free(den_copy);
    return found;
}

int tf_is_minimum_phase(const TransferFn *tf, int is_dt)
{
    if (!tf) return 0;
    int nz = tf->num.degree;
    double *zeros = (double *)malloc(nz * sizeof(double));
    if (!zeros) return 0;

    /* Find zeros by evaluating numerator */

    /* For the poles of the feedforward, we check plant zeros */
    /* A continuous-time system is minimum-phase if all zeros are in LHP.
     * A discrete-time system is minimum-phase if all zeros are inside unit circle. */
    /* Simplified check for low-order systems */
    if (nz == 0) {
        free(zeros);
        return 1; /* No zeros, vacuously minimum-phase */
    }
    if (nz == 1) {
        double zero_val = -tf->num.coeff[0] / tf->num.coeff[1];
        free(zeros);
        if (is_dt) return fabs(zero_val) < 1.0;
        else       return zero_val < 0.0;
    }
    if (nz == 2) {
        double a = tf->num.coeff[2];
        double b = tf->num.coeff[1];
        double c = tf->num.coeff[0];
        double disc = b * b - 4.0 * a * c;
        free(zeros);
        if (disc >= 0) {
            double z1 = (-b - sqrt(disc)) / (2.0 * a);
            double z2 = (-b + sqrt(disc)) / (2.0 * a);
            if (is_dt) return (fabs(z1) < 1.0) && (fabs(z2) < 1.0);
            else       return (z1 < 0.0) && (z2 < 0.0);
        } else {
            /* Complex zeros */
            double real_part = -b / (2.0 * a);
            if (is_dt) {
                double imag_part = sqrt(-disc) / (2.0 * a);
                return (real_part * real_part + imag_part * imag_part) < 1.0;
            } else {
                return real_part < 0.0;
            }
        }
    }

    free(zeros);
    /* For higher order, approximate: check Routh-Hurwitz criterion
     * for continuous-time; Jury criterion for discrete-time. */
    return 1; /* Conservatively assume minimum-phase for higher order */
}

int tf_to_state_space(const TransferFn *tf,
                      double *A, double *B, double *C, double *D)
{
    if (!tf || !A || !B || !C || !D) return -1;
    int n = tf->den.degree;
    if (n <= 0) return -1;

    /* Normalize: ensure denominator is monic (leading coefficient = 1) */
    double lead = tf->den.coeff[n];
    double *a = (double *)malloc(n * sizeof(double));
    double *b = (double *)malloc((n + 1) * sizeof(double));
    if (!a || !b) { free(a); free(b); return -1; }

    for (int i = 0; i < n; i++) {
        a[i] = tf->den.coeff[i] / lead;
    }
    for (int i = 0; i <= n; i++) {
        b[i] = (i <= tf->num.degree) ? tf->num.coeff[i] : 0.0;
        b[i] *= tf->gain / lead;
    }

    /* Controllable canonical form:
     * A = [0    1    0   ... 0  ]
     *     [0    0    1   ... 0  ]
     *     ...
     *     [-a0 -a1  -a2 ... -a_{n-1}]
     * B = [0 ... 0 1]^T
     * C = [b0-a0*bn, b1-a1*bn, ..., b_{n-1}-a_{n-1}*bn]
     * D = bn
     */
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = (j == i + 1) ? 1.0 : 0.0;
        }
    }
    for (int j = 0; j < n; j++) {
        A[(n - 1) * n + j] = -a[j];
    }

    for (int i = 0; i < n; i++) {
        B[i] = (i == n - 1) ? 1.0 : 0.0;
    }

    /* Direct feedthrough from strict properness */
    *D = (tf->num.degree >= tf->den.degree) ? b[n] : 0.0;

    for (int i = 0; i < n; i++) {
        C[i] = b[i] - a[i] * (*D);
    }

    free(a);
    free(b);
    return n;
}

/* ==========================================================================
 * L3: Transfer Function Algebra
 * ========================================================================== */

TransferFn tf_parallel(const TransferFn *g1, const TransferFn *g2)
{
    TransferFn result;
    /* G1 + G2 = (num1*den2 + num2*den1) / (den1*den2) */
    Poly n1d2 = poly_mul(&g1->num, &g2->den);
    Poly n2d1 = poly_mul(&g2->num, &g1->den);
    /* Scale by gains */
    for (int i = 0; i <= n1d2.degree; i++) n1d2.coeff[i] *= g1->gain;
    for (int i = 0; i <= n2d1.degree; i++) n2d1.coeff[i] *= g2->gain;
    result.num = poly_add(&n1d2, &n2d1);
    result.den = poly_mul(&g1->den, &g2->den);
    result.gain = 1.0;
    poly_free(&n1d2); poly_free(&n2d1);
    return result;
}

TransferFn tf_series(const TransferFn *g1, const TransferFn *g2)
{
    TransferFn result;
    /* G1 * G2 = (num1*num2) / (den1*den2) */
    result.num = poly_mul(&g1->num, &g2->num);
    result.den = poly_mul(&g1->den, &g2->den);
    result.gain = g1->gain * g2->gain;
    return result;
}

TransferFn tf_feedback(const TransferFn *g, const TransferFn *h)
{
    TransferFn result;
    /* G / (1 + G*H) = numG*denH / (denG*denH + numG*numH) */
    Poly num = poly_mul(&g->num, &h->den);
    for (int i = 0; i <= num.degree; i++) num.coeff[i] *= g->gain;
    result.num = num;

    Poly den_gh = poly_mul(&g->den, &h->den);
    Poly num_gh = poly_mul(&g->num, &h->num);
    for (int i = 0; i <= num_gh.degree; i++) {
        num_gh.coeff[i] *= g->gain * h->gain;
    }
    result.den = poly_add(&den_gh, &num_gh);
    result.gain = 1.0;

    poly_free(&den_gh); poly_free(&num_gh);
    return result;
}

TransferFn tf_negate(const TransferFn *tf)
{
    TransferFn result;
    result.num = poly_create(tf->num.coeff, tf->num.degree);
    result.den = poly_create(tf->den.coeff, tf->den.degree);
    result.gain = -tf->gain;
    return result;
}

TransferFn tf_scale(const TransferFn *tf, double k)
{
    TransferFn result;
    result.num = poly_create(tf->num.coeff, tf->num.degree);
    result.den = poly_create(tf->den.coeff, tf->den.degree);
    result.gain = tf->gain * k;
    return result;
}

int tf_equals(const TransferFn *a, const TransferFn *b, double tol)
{
    if (!a || !b) return 0;
    if (a->num.degree != b->num.degree || a->den.degree != b->den.degree)
        return 0;
    for (int i = 0; i <= a->num.degree; i++) {
        double diff = a->num.coeff[i] * a->gain - b->num.coeff[i] * b->gain;
        if (fabs(diff) > tol) return 0;
    }
    for (int i = 0; i <= a->den.degree; i++) {
        if (fabs(a->den.coeff[i] - b->den.coeff[i]) > tol) return 0;
    }
    return 1;
}
