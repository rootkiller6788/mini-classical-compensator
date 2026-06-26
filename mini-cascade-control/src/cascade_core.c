/**
 * @file cascade_core.c
 * @brief Core cascade control operations
 *
 * L1 -- Definitions: Transfer function creation, manipulation, evaluation.
 * L2 -- Core Concepts: Cascade loop formation (inner CL, equivalent plant).
 * L3 -- Mathematical Structures: Polynomial algebra over R[s], TF composition.
 * L4 -- Fundamental Laws: Routh-Hurwitz stability criterion.
 * L5 -- Computational Methods: Frequency response, pole finding, step response.
 * L6 -- Canonical Systems: DC motor, reactor, flow-pressure, level-tank models.
 *
 * Course alignment:
 *   MIT 6.302 -- Feedback loop algebra and stability
 *   Stanford ENGR105 -- Controller implementation
 *   Berkeley ME232 -- Cascade control structures
 *   Caltech CDS 110 -- Stability criteria
 *   Tsinghua -- Transfer function operations and analysis
 */

#include "cascade_types.h"
#include "cascade_analysis.h"
#include "cascade_design.h"
#include "cascade_implementation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==========================================================================
 * L1: Safe Memory Allocation
 * ========================================================================== */

static void *safe_alloc(size_t size) {
    void *ptr = malloc(size);
    if (!ptr && size > 0) {
        fprintf(stderr, "cascade_core: malloc(%zu) failed\n", size);
    }
    return ptr;
}

static void *safe_calloc(size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    if (!ptr && nmemb > 0 && size > 0) {
        fprintf(stderr, "cascade_core: calloc(%zu,%zu) failed\n", nmemb, size);
    }
    return ptr;
}

/* ==========================================================================
 * L3: Polynomial Operations over R[s]
 * ========================================================================== */

CascadePoly cascade_poly_create(const double *coeff, int degree) {
    CascadePoly p;
    p.degree = degree;
    if (degree >= 0 && coeff) {
        p.coeff = (double *)safe_alloc((size_t)(degree + 1) * sizeof(double));
        if (p.coeff) {
            memcpy(p.coeff, coeff, (size_t)(degree + 1) * sizeof(double));
        }
    } else {
        p.coeff = NULL;
    }
    return p;
}

void cascade_poly_free(CascadePoly *p) {
    if (p && p->coeff) {
        free(p->coeff);
        p->coeff = NULL;
        p->degree = 0;
    }
}

void cascade_poly_eval_complex(const CascadePoly *p,
                                double sigma, double omega,
                                double *re, double *im) {
    /* Horner's method for complex arguments.
     * p(s) = a0 + s*(a1 + s*(a2 + ... + s*an)...)
     * Complexity: O(n) exactly n multiplications and n additions.
     * Theorem (Horner, 1819): Optimal evaluation of polynomials. */
    *re = 0.0;
    *im = 0.0;
    if (!p || !p->coeff || p->degree < 0) return;
    for (int i = p->degree; i >= 0; i--) {
        double new_re = (*re) * sigma - (*im) * omega + p->coeff[i];
        double new_im = (*re) * omega + (*im) * sigma;
        *re = new_re;
        *im = new_im;
    }
}

CascadePoly cascade_poly_mul(const CascadePoly *p, const CascadePoly *q) {
    /* Convolution: r_k = sum_{i+j=k} p_i * q_j
     * Makes (R[s], +, *) a commutative ring.
     * Complexity: O(deg(p) * deg(q)) */
    CascadePoly r;
    if (!p || !q || p->degree < 0 || q->degree < 0) {
        r.degree = -1; r.coeff = NULL; return r;
    }
    int r_deg = p->degree + q->degree;
    r.degree = r_deg;
    r.coeff = (double *)safe_calloc((size_t)(r_deg + 1), sizeof(double));
    if (r.coeff) {
        for (int i = 0; i <= p->degree && p->coeff; i++) {
            for (int j = 0; j <= q->degree && q->coeff; j++) {
                r.coeff[i + j] += p->coeff[i] * q->coeff[j];
            }
        }
    }
    return r;
}

CascadePoly cascade_poly_add(const CascadePoly *p, const CascadePoly *q) {
    if (!p || !q) { CascadePoly r = {-1, NULL}; return r; }
    int max_deg = (p->degree > q->degree) ? p->degree : q->degree;
    CascadePoly r;
    r.degree = max_deg;
    r.coeff = (double *)safe_calloc((size_t)(max_deg + 1), sizeof(double));
    if (r.coeff) {
        for (int i = 0; i <= p->degree && p->coeff; i++)
            r.coeff[i] += p->coeff[i];
        for (int i = 0; i <= q->degree && q->coeff; i++)
            r.coeff[i] += q->coeff[i];
    }
    return r;
}

/* ==========================================================================
 * L1: Transfer Function Construction
 * ========================================================================== */

CascadeTF cascade_tf_create(const double *num_coeff, int num_deg,
                              const double *den_coeff, int den_deg,
                              double gain) {
    CascadeTF tf;
    tf.num = cascade_poly_create(num_coeff, num_deg);
    tf.den = cascade_poly_create(den_coeff, den_deg);
    tf.gain = gain;
    tf.delay = 0.0;
    tf.has_delay = 0;
    return tf;
}

void cascade_tf_free(CascadeTF *tf) {
    if (tf) {
        cascade_poly_free(&tf->num);
        cascade_poly_free(&tf->den);
        tf->gain = 1.0;
        tf->delay = 0.0;
        tf->has_delay = 0;
    }
}

/* ==========================================================================
 * L5: Frequency Response Computation
 * ========================================================================== */

int cascade_tf_freq_response(const CascadeTF *tf, double omega,
                              double *mag, double *phase) {
    /* G(jw) = K * num(jw) / den(jw)
     * Complex division: (Nr+j*Ni)/(Dr+j*Di) = (Nr*Dr+Ni*Di + j*(Ni*Dr-Nr*Di))/(Dr^2+Di^2)
     * For delayed systems: G(jw) = G0(jw) * exp(-j*w*theta) */
    if (!tf || !mag || !phase) return -1;
    double num_re, num_im, den_re, den_im;
    cascade_poly_eval_complex(&tf->num, 0.0, omega, &num_re, &num_im);
    cascade_poly_eval_complex(&tf->den, 0.0, omega, &den_re, &den_im);
    double den_mag_sq = den_re * den_re + den_im * den_im;
    if (den_mag_sq < 1e-30) { *mag = 1e300; *phase = 0.0; return -1; }
    double g_re = (num_re * den_re + num_im * den_im) / den_mag_sq;
    double g_im = (num_im * den_re - num_re * den_im) / den_mag_sq;
    double absK = fabs(tf->gain);
    g_re *= tf->gain; g_im *= tf->gain;
    double abs_g = sqrt(g_re * g_re + g_im * g_im);
    double arg_g = atan2(g_im, g_re);
    if (tf->has_delay && tf->delay > 0.0) {
        arg_g -= omega * tf->delay;
        while (arg_g > M_PI) arg_g -= 2.0 * M_PI;
        while (arg_g < -M_PI) arg_g += 2.0 * M_PI;
    }
    *mag = abs_g; *phase = arg_g;
    (void)absK;
    return 0;
}

int cascade_tf_close_loop(const CascadeTF *L, CascadeTF *cl) {
    /* Unity negative feedback: G_cl = L / (1 + L)
     * L = K*N_L/D_L => G_cl = K*N_L/(D_L + K*N_L) */
    if (!L || !cl) return -1;
    CascadePoly num_scaled = cascade_poly_create(L->num.coeff, L->num.degree);
    if (num_scaled.coeff) {
        for (int i = 0; i <= num_scaled.degree; i++)
            num_scaled.coeff[i] *= L->gain;
    }
    CascadePoly den_cl = cascade_poly_add(&L->den, &num_scaled);
    cl->num = num_scaled; cl->den = den_cl;
    cl->gain = 1.0; cl->delay = 0.0; cl->has_delay = 0;
    return 0;
}

/* ==========================================================================
 * L1: PID <-> Transfer Function Conversion
 * ========================================================================== */

int cascade_pid_to_tf(const CascadePID *pid, CascadeTF *tf) {
    /* C(s) = Kp + Ki/s + Kd*s/(N*s+1)
     * = [(Kp*N+Kd)*s^2 + (Kp+Ki*N)*s + Ki] / [N*s^2 + s] */
    if (!pid || !tf) return -1;
    if (fabs(pid->Kd) < 1e-15 && fabs(pid->Ki) < 1e-15) {
        double nc[] = {pid->Kp}, dc[] = {1.0};
        *tf = cascade_tf_create(nc, 0, dc, 0, 1.0);
        return 0;
    }
    if (fabs(pid->Kd) < 1e-15) {
        double nc[] = {pid->Ki, pid->Kp}, dc[] = {0.0, 1.0};
        *tf = cascade_tf_create(nc, 1, dc, 1, 1.0);
        return 0;
    }
    double Nv = (pid->N > 0.0) ? pid->N : 100.0;
    double nc[] = {pid->Ki, pid->Kp + pid->Ki * Nv, pid->Kp * Nv + pid->Kd};
    double dc[] = {0.0, 1.0, Nv};
    *tf = cascade_tf_create(nc, 2, dc, 2, 1.0);
    return 0;
}

int cascade_pid_freq_response(const CascadePID *pid, double omega,
                               double *mag, double *phase) {
    if (!pid || !mag || !phase) return -1;
    double w = (omega < 1e-15) ? 1e-15 : omega;
    double c_re = pid->Kp;
    double c_im = -pid->Ki / w;
    double Nv = (pid->N > 0.0) ? pid->N : 100.0;
    double den = 1.0 + Nv * Nv * w * w;
    c_re += pid->Kd * Nv * w * w / den;
    c_im += pid->Kd * w / den;
    *mag = sqrt(c_re * c_re + c_im * c_im);
    *phase = atan2(c_im, c_re);
    return 0;
}

/* ==========================================================================
 * L2: Cascade Loop Formation
 * ========================================================================== */

int cascade_inner_closed_loop(const CascadeTF *Gi,
                               const CascadePID *Ci,
                               CascadeTF *Gi_cl) {
    /* Gi_cl = Ci*Gi/(1 + Ci*Gi) -- step 1 of sequential closure */
    if (!Gi || !Ci || !Gi_cl) return -1;
    CascadeTF ci_tf;
    if (cascade_pid_to_tf(Ci, &ci_tf) != 0) return -1;
    CascadePoly num_L = cascade_poly_mul(&ci_tf.num, &Gi->num);
    CascadePoly den_L = cascade_poly_mul(&ci_tf.den, &Gi->den);
    double gL = Gi->gain;
    if (num_L.coeff && fabs(gL - 1.0) > 1e-15)
        for (int i = 0; i <= num_L.degree; i++) num_L.coeff[i] *= gL;
    CascadePoly den_cl = cascade_poly_add(&den_L, &num_L);
    *Gi_cl = cascade_tf_create(num_L.coeff, num_L.degree,
                                den_cl.coeff, den_cl.degree, 1.0);
    cascade_poly_free(&num_L); cascade_poly_free(&den_L);
    cascade_poly_free(&den_cl); cascade_tf_free(&ci_tf);
    return 0;
}

int cascade_form_equivalent_plant(const CascadeTF *Gi,
                                   const CascadePID *Ci,
                                   const CascadeTF *Go,
                                   CascadeTF *eq_plant) {
    /* Geq = Gi_cl * Go */
    if (!Gi || !Ci || !Go || !eq_plant) return -1;
    CascadeTF Gi_cl;
    if (cascade_inner_closed_loop(Gi, Ci, &Gi_cl) != 0) return -1;
    CascadePoly num_eq = cascade_poly_mul(&Gi_cl.num, &Go->num);
    CascadePoly den_eq = cascade_poly_mul(&Gi_cl.den, &Go->den);
    double gEq = Go->gain;
    if (num_eq.coeff && fabs(gEq - 1.0) > 1e-15)
        for (int i = 0; i <= num_eq.degree; i++) num_eq.coeff[i] *= gEq;
    eq_plant->num = num_eq; eq_plant->den = den_eq;
    eq_plant->gain = 1.0; eq_plant->delay = 0.0; eq_plant->has_delay = 0;
    cascade_tf_free(&Gi_cl);
    return 0;
}

int cascade_overall_closed_loop(const CascadeTF *Gi_cl,
                                 const CascadePID *Co,
                                 const CascadeTF *Go,
                                 CascadeTF *overall) {
    /* y/r = Co*Gi_cl*Go / (1 + Co*Gi_cl*Go) */
    if (!Gi_cl || !Co || !Go || !overall) return -1;
    CascadeTF co_tf;
    if (cascade_pid_to_tf(Co, &co_tf) != 0) return -1;
    CascadePoly nt = cascade_poly_mul(&co_tf.num, &Gi_cl->num);
    CascadePoly num_L = cascade_poly_mul(&nt, &Go->num);
    CascadePoly dt = cascade_poly_mul(&co_tf.den, &Gi_cl->den);
    CascadePoly den_L = cascade_poly_mul(&dt, &Go->den);
    double gL = Go->gain;
    if (num_L.coeff && fabs(gL - 1.0) > 1e-15)
        for (int i = 0; i <= num_L.degree; i++) num_L.coeff[i] *= gL;
    CascadePoly den_cl = cascade_poly_add(&den_L, &num_L);
    overall->num = num_L; overall->den = den_cl;
    overall->gain = 1.0; overall->delay = 0.0; overall->has_delay = 0;
    cascade_poly_free(&nt); cascade_poly_free(&dt);
    cascade_poly_free(&den_L); cascade_tf_free(&co_tf);
    return 0;
}

/* ==========================================================================
 * L4: Routh-Hurwitz Stability Criterion
 * ========================================================================== */

int cascade_routh_hurwitz(const double *den, int degree, int *n_rhp) {
    /* Theorem (Routh, 1874; Hurwitz, 1895):
     * All roots of polynomial have Re < 0 iff all first-column elements
     * of Routh array have the same sign.
     * Number of RHP roots = number of sign changes in first column.
     * Complexity: O(n^2). Reference: Gantmacher (1959), Vol. II. */
    if (!den || degree <= 0 || !n_rhp) return -1;
    int rows = degree + 1;
    int mc = (degree + 2) / 2;
    double **r = (double **)malloc((size_t)rows * sizeof(double *));
    if (!r) return -1;
    int ok = 1;
    for (int i = 0; i < rows; i++) {
        r[i] = (double *)calloc((size_t)mc, sizeof(double));
        if (!r[i]) ok = 0;
    }
    if (!ok) {
        for (int i = 0; i < rows; i++) free(r[i]);
        free(r); return -1;
    }
    for (int i = 0; i <= degree; i++) {
        int ci = degree - i;
        if (i % 2 == 0) r[0][i / 2] = den[ci];
        else r[1][i / 2] = den[ci];
    }
    int nrhp = 0;
    double ps = (r[0][0] != 0.0) ? r[0][0] : 1e-300;
    int pp = (ps > 0.0) ? 1 : 0;

    /* Check row 1 sign for low-degree polynomials */
    if (rows >= 2 && r[1][0] != 0.0) {
        int cp1 = (r[1][0] > 0.0) ? 1 : 0;
        if (cp1 != pp) { nrhp++; pp = cp1; }
    }

    for (int i = 2; i < rows; i++) {
        double pv = r[i - 1][0];
        if (fabs(pv) < 1e-15) { pv = 1e-10; r[i - 1][0] = pv; }
        for (int j = 0; j < mc - 1; j++) {
            double a = r[i - 2][0], b = r[i - 2][j + 1];
            double c = r[i - 1][0], d = r[i - 1][j + 1];
            r[i][j] = -(a * d - b * c) / pv;
        }
        if (fabs(r[i][0]) < 1e-15) r[i][0] = 0.0;
        if (r[i][0] != 0.0) {
            int cp = (r[i][0] > 0.0) ? 1 : 0;
            if (cp != pp) { nrhp++; pp = cp; }
        }
    }
    *n_rhp = nrhp;
    for (int i = 0; i < rows; i++) free(r[i]);
    free(r);
    return 0;
}

int cascade_is_stable(const CascadeTF *tf) {
    if (!tf || tf->den.degree < 0) return 0;
    int nrhp = 0;
    if (cascade_routh_hurwitz(tf->den.coeff, tf->den.degree, &nrhp) != 0)
        return 0;
    return (nrhp == 0) ? 1 : 0;
}

int cascade_verify_internal_stability(const CascadeTF *Gi,
                                       const CascadePID *Ci,
                                       const CascadeTF *Go,
                                       const CascadePID *Co) {
    if (!Gi || !Ci || !Go || !Co) return 0;
    CascadeTF Gi_cl;
    if (cascade_inner_closed_loop(Gi, Ci, &Gi_cl) != 0) return 0;
    int inner_stable = cascade_is_stable(&Gi_cl);
    if (!inner_stable) { cascade_tf_free(&Gi_cl); return 0; }
    CascadeTF overall;
    if (cascade_overall_closed_loop(&Gi_cl, Co, Go, &overall) != 0) {
        cascade_tf_free(&Gi_cl); return 0;
    }
    int outer_stable = cascade_is_stable(&overall);
    cascade_tf_free(&Gi_cl); cascade_tf_free(&overall);
    if (!cascade_is_stable(Gi)) return 0;
    if (!cascade_is_stable(Go)) return 0;
    return outer_stable ? 1 : 0;
}

/* ==========================================================================
 * L5: Frequency Response Analysis
 * ========================================================================== */

int cascade_bode_analysis(const CascadeTF *tf,
                           double freq_min, double freq_max,
                           int n_points,
                           CascadeFreqResponse *resp) {
    if (!tf || !resp || n_points < 2 || freq_min <= 0 || freq_max <= freq_min)
        return -1;
    resp->num_points = n_points;
    resp->freq_min = freq_min; resp->freq_max = freq_max;
    resp->points = (CascadeFreqPoint *)safe_alloc(
        (size_t)n_points * sizeof(CascadeFreqPoint));
    if (!resp->points) return -1;
    double lr = log(freq_max / freq_min);
    double pmdb = 0.0, po = 0.0;
    int fgc = 0, fpc = 0;
    for (int i = 0; i < n_points; i++) {
        double w = freq_min * exp(lr * (double)i / (double)(n_points - 1));
        double mag, ph;
        cascade_tf_freq_response(tf, w, &mag, &ph);
        double mdb = 20.0 * log10(mag > 1e-300 ? mag : 1e-300);
        double pd = ph * 180.0 / M_PI;
        resp->points[i].omega = w;
        resp->points[i].magnitude = mag;
        resp->points[i].magnitude_db = mdb;
        resp->points[i].phase_rad = ph;
        resp->points[i].phase_deg = pd;
        resp->points[i].real_part = mag * cos(ph);
        resp->points[i].imag_part = mag * sin(ph);
        if (i > 0 && !fgc &&
            ((pmdb > 0.0 && mdb < 0.0) || (pmdb < 0.0 && mdb > 0.0))) {
            double alpha = fabs(pmdb) / (fabs(pmdb) + fabs(mdb) + 1e-15);
            resp->gain_crossover = po + alpha * (w - po);
            fgc = 1;
        }
        if (i > 0 && !fpc && pmdb > 0.0 && pd < -175.0) {
            resp->phase_crossover = po; fpc = 1;
        }
        if (i == 0) resp->dc_gain_db = mdb;
        pmdb = mdb; po = w;
    }
    if (!fgc) resp->gain_crossover = freq_max;
    if (!fpc) resp->phase_crossover = freq_max;
    return 0;
}

int cascade_stability_margins(const CascadeFreqResponse *resp,
                               double *gm, double *pm) {
    if (!resp || !gm || !pm || !resp->points || resp->num_points < 2)
        return -1;
    double bwgc = resp->freq_min, be = 1e300;
    for (int i = 0; i < resp->num_points; i++) {
        double e = fabs(resp->points[i].magnitude_db);
        if (e < be) { be = e; bwgc = resp->points[i].omega; }
    }
    double pmv = 0.0;
    for (int i = 0; i < resp->num_points; i++) {
        if (fabs(resp->points[i].omega - bwgc) < 1e-6) {
            pmv = 180.0 + resp->points[i].phase_deg;
            if (pmv > 180.0) pmv -= 360.0;
            if (pmv < -180.0) pmv += 360.0;
            break;
        }
    }
    double bwpc = resp->freq_min; be = 1e300;
    for (int i = 0; i < resp->num_points; i++) {
        double e = fabs(resp->points[i].phase_deg + 180.0);
        if (e < be) { be = e; bwpc = resp->points[i].omega; }
    }
    double gmv = 20.0;
    for (int i = 0; i < resp->num_points; i++) {
        if (fabs(resp->points[i].omega - bwpc) < 1e-6) {
            gmv = -resp->points[i].magnitude_db; break;
        }
    }
    *pm = pmv; *gm = gmv;
    return 0;
}

double cascade_max_sensitivity(const CascadeFreqResponse *resp) {
    if (!resp || !resp->points || resp->num_points < 1) return 1e300;
    double ms = 1.0;
    for (int i = 0; i < resp->num_points; i++) {
        double re = resp->points[i].real_part;
        double im = resp->points[i].imag_part;
        double sm = 1.0 / sqrt((1.0 + re) * (1.0 + re) + im * im + 1e-15);
        if (sm > ms) ms = sm;
    }
    return ms;
}

void cascade_freq_response_free(CascadeFreqResponse *resp) {
    if (resp) { free(resp->points); resp->points = NULL; resp->num_points = 0; }
}

int cascade_loop_freq_response(const CascadeTF *Gi,
                                const CascadePID *Ci,
                                const CascadeTF *Go,
                                const CascadePID *Co,
                                double omega,
                                double *mag, double *phase) {
    if (!Gi || !Ci || !Go || !Co || !mag || !phase) return -1;
    double gim, gip, cim, cip;
    cascade_tf_freq_response(Gi, omega, &gim, &gip);
    cascade_pid_freq_response(Ci, omega, &cim, &cip);
    double lim = cim * gim, lip = cip + gip;
    double lr = lim * cos(lip), li = lim * sin(lip);
    double dr = 1.0 + lr, di = li;
    double dm = sqrt(dr * dr + di * di);
    double gclm = lim / (dm > 1e-15 ? dm : 1e-15);
    double gclp = lip - atan2(di, dr);
    double com, cop, gom, gop;
    cascade_pid_freq_response(Co, omega, &com, &cop);
    cascade_tf_freq_response(Go, omega, &gom, &gop);
    *mag = com * gclm * gom;
    *phase = cop + gclp + gop;
    while (*phase > M_PI) *phase -= 2.0 * M_PI;
    while (*phase < -M_PI) *phase += 2.0 * M_PI;
    return 0;
}

int cascade_closed_loop_bandwidth(const CascadeTF *cl_tf, double *bw) {
    if (!cl_tf || !bw) return -1;
    double wl = 0.001, wh = 1000.0, ml, mh, ph;
    cascade_tf_freq_response(cl_tf, wl, &ml, &ph);
    cascade_tf_freq_response(cl_tf, wh, &mh, &ph);
    double dcdb = 20.0 * log10(ml > 1e-300 ? ml : 1e-300);
    double tdb = dcdb - 3.0;
    double hdb = 20.0 * log10(mh > 1e-300 ? mh : 1e-300);
    if (hdb > tdb) { *bw = wh; return 0; }
    for (int i = 0; i < 50; i++) {
        double wm = (wl + wh) / 2.0, mm;
        cascade_tf_freq_response(cl_tf, wm, &mm, &ph);
        double mdb = 20.0 * log10(mm > 1e-300 ? mm : 1e-300);
        if (mdb > tdb) wl = wm; else wh = wm;
        if (wh - wl < 1e-6) break;
    }
    *bw = (wl + wh) / 2.0;
    return 0;
}

/* ==========================================================================
 * L5: Time-Domain Analysis
 * ========================================================================== */

double cascade_step_response(const CascadeTF *tf, double t, int n_terms) {
    if (!tf || t < 0.0) return 0.0;
    /* First-order: K/(tau*s+1) => y(t) = K*(1 - exp(-t/tau)) */
    if (tf->den.degree == 1 && tf->num.degree == 0 && !tf->has_delay) {
        double a0 = tf->den.coeff[0], a1 = tf->den.coeff[1];
        double b0 = tf->num.coeff[0];
        if (fabs(a1) < 1e-15) return 0.0;
        double K = tf->gain * b0 / a0;
        double tau = a1 / a0;
        if (tau < 0) tau = -tau;
        return K * (1.0 - exp(-t / (tau > 1e-15 ? tau : 1e-15)));
    }
    /* Second-order: K*wn^2/(s^2+2*zeta*wn*s+wn^2) */
    if (tf->den.degree == 2 && tf->num.degree == 0 && !tf->has_delay) {
        double a0 = tf->den.coeff[0], a1 = tf->den.coeff[1], a2 = tf->den.coeff[2];
        if (fabs(a2) < 1e-15) return 0.0;
        double wn = sqrt(fabs(a0) / a2);
        double zeta = a1 / (2.0 * a2 * wn);
        double K = tf->gain * tf->num.coeff[0] / (a0 > 1e-15 ? a0 : 1e-15);
        if (zeta < 1.0 && wn > 0) {
            double wd = wn * sqrt(1.0 - zeta * zeta);
            double phi = atan2(sqrt(1.0 - zeta * zeta), zeta);
            return K * (1.0 - exp(-zeta * wn * t) *
                        sin(wd * t + phi) / sqrt(1.0 - zeta * zeta));
        } else if (fabs(zeta - 1.0) < 1e-6) {
            return K * (1.0 - (1.0 + wn * t) * exp(-wn * t));
        } else {
            double s1 = -wn * (zeta - sqrt(zeta * zeta - 1.0));
            double s2 = -wn * (zeta + sqrt(zeta * zeta - 1.0));
            return K * (1.0 + (s2 * exp(s1 * t) - s1 * exp(s2 * t)) / (s1 - s2));
        }
    }
    /* General: approximate with dominant time constant */
    (void)n_terms;
    double Kdc = 1.0;
    if (tf->den.coeff && fabs(tf->den.coeff[0]) > 1e-15)
        Kdc = tf->gain * tf->num.coeff[0] / tf->den.coeff[0];
    double tau_d = 1.0;
    if (tf->den.degree >= 1 && fabs(tf->den.coeff[0]) > 1e-15) {
        double sa = 0.0;
        for (int i = 1; i <= tf->den.degree; i++) sa += fabs(tf->den.coeff[i]);
        tau_d = sa / fabs(tf->den.coeff[0]);
        if (tau_d < 1e-6) tau_d = 1.0;
    }
    if (t > 10.0 * tau_d) return Kdc;
    return Kdc * (1.0 - exp(-t / tau_d));
}

int cascade_step_response_vector(const CascadeTF *tf,
                                  const double *t, int n, double *y) {
    if (!tf || !t || !y || n <= 0) return -1;
    for (int i = 0; i < n; i++)
        y[i] = cascade_step_response(tf, t[i], 20);
    return 0;
}

int cascade_find_poles(const CascadeTF *tf, double *poles, int max_n) {
    if (!tf || !poles || max_n <= 0) return -1;
    int n = tf->den.degree;
    if (n <= 0) return 0;
    int count = (n < max_n) ? n : max_n;
    double an = tf->den.coeff[n];
    if (fabs(an) < 1e-15) return -1;
    if (n == 1) {
        poles[0] = -tf->den.coeff[0] / tf->den.coeff[1];
        poles[1] = 0.0;
        return 1;
    }
    if (n == 2) {
        double a = tf->den.coeff[2] / an;
        double b = tf->den.coeff[1] / an;
        double c = tf->den.coeff[0] / an;
        double disc = b * b - 4.0 * a * c;
        if (disc >= 0) {
            double sq = sqrt(disc);
            poles[0] = (-b + sq) / (2.0 * a); poles[1] = 0.0;
            poles[2] = (-b - sq) / (2.0 * a); poles[3] = 0.0;
        } else {
            double re = -b / (2.0 * a);
            double im = sqrt(-disc) / (2.0 * a);
            poles[0] = re; poles[1] = im;
            poles[2] = re; poles[3] = -im;
        }
        return 2;
    }
    for (int i = 0; i < count; i++) {
        if (i < n && fabs(tf->den.coeff[n - 1 - i]) > 1e-15)
            poles[2 * i] = -tf->den.coeff[i] / tf->den.coeff[n - 1 - i];
        else
            poles[2 * i] = -1.0;
        poles[2 * i + 1] = 0.0;
    }
    return count;
}

/* ==========================================================================
 * L1: Cascade System Lifecycle
 * ========================================================================== */

void cascade_set_default_spec(CascadeDesignSpec *spec) {
    if (!spec) return;
    spec->inner_pm_target = 60.0;
    spec->outer_pm_target = 45.0;
    spec->inner_gm_target = 10.0;
    spec->outer_gm_target = 8.0;
    spec->inner_bw_min = 2.0;
    spec->outer_bw_max = 0.5;
    spec->bw_ratio_min = 3.0;
    spec->bw_ratio_optimal = 5.0;
    spec->inner_settle_max = 2.0;
    spec->outer_settle_max = 15.0;
    spec->outer_overshoot_max = 0.10;
    spec->outer_ess_max = 0.02;
    spec->inner_noise_sens_max = 1.8;
    spec->use_freq_domain = 0;
    spec->inner_controller_type = 1;
    spec->outer_controller_type = 1;
}

static void init_default_pid(CascadePID *pid, int is_inner) {
    if (is_inner) {
        pid->Kp = 1.0; pid->Ki = 5.0; pid->Kd = 0.0; pid->N = 10.0;
    } else {
        pid->Kp = 0.5; pid->Ki = 0.2; pid->Kd = 0.0; pid->N = 10.0;
    }
    pid->b = 1.0; pid->c = 0.0; pid->Ts = 0.0;
    pid->u_min = -100.0; pid->u_max = 100.0;
    pid->integrator = 0.0; pid->prev_error = 0.0; pid->prev_y = 0.0;
    pid->has_antiwindup = 1; pid->Tt = 1.0;
}

void cascade_system_init(CascadeSystem *sys, const char *name) {
    if (!sys) return;
    sys->system_name = name;
    sys->cascade_active = 0;
    sys->bumpless_ready = 0;
    sys->bandwidth_ratio = 5.0;
    sys->decoupling_factor = 0.0;

    double n1[] = {1.0}, df[] = {1.0, 0.2}, ds[] = {1.0, 2.0};
    sys->inner.inner_process = cascade_tf_create(n1, 0, df, 1, 1.0);
    init_default_pid(&sys->inner.inner_controller, 1);
    sys->inner.inner_bandwidth = 5.0;
    sys->inner.inner_rise_time = 0.3;
    sys->inner.inner_settle_time = 1.0;
    sys->inner.inner_overshoot = 0.05;
    sys->inner.inner_is_stable = 1;
    sys->inner.inner_var_name = "Inner PV";
    sys->inner.inner_cl = cascade_tf_create(n1, 0, n1, 0, 1.0);

    sys->outer.outer_process = cascade_tf_create(n1, 0, ds, 1, 1.0);
    init_default_pid(&sys->outer.outer_controller, 0);
    sys->outer.outer_bandwidth = 0.8;
    sys->outer.outer_rise_time = 2.0;
    sys->outer.outer_settle_time = 10.0;
    sys->outer.outer_overshoot = 0.10;
    sys->outer.outer_is_stable = 1;
    sys->outer.outer_var_name = "Outer PV";
    sys->outer.equivalent_plant = cascade_tf_create(n1, 0, ds, 1, 1.0);
    sys->outer.outer_cl = cascade_tf_create(n1, 0, ds, 1, 1.0);
}

void cascade_system_free(CascadeSystem *sys) {
    if (!sys) return;
    cascade_tf_free(&sys->inner.inner_process);
    cascade_tf_free(&sys->inner.inner_cl);
    cascade_tf_free(&sys->outer.outer_process);
    cascade_tf_free(&sys->outer.equivalent_plant);
    cascade_tf_free(&sys->outer.outer_cl);
}

/* ==========================================================================
 * L6: Canonical System Models
 * ========================================================================== */

int cascade_create_dc_motor_model(const CascadeDCMotor *motor,
                                   CascadeTF *vel_tf, CascadeTF *pos_tf) {
    /* DC motor velocity TF: w(s)/V(s) = Kt / (J*L*s^2 + (J*R+B*L)*s + (B*R+Kb*Kt))
     * Position TF: theta(s)/w_ref(s) = 1/s
     * References: Krause, Wasynczuk & Sudhoff (2002), "Analysis of Electric Machinery" */
    if (!motor || !vel_tf || !pos_tf) return -1;
    double a0 = motor->B * motor->R + motor->Kb * motor->Kt;
    double a1 = motor->J * motor->R + motor->B * motor->L;
    double a2 = motor->J * motor->L;
    double vn[] = {motor->Kt}, vd[] = {a0, a1, a2};
    double pn[] = {1.0}, pd[] = {0.0, 1.0};
    *vel_tf = cascade_tf_create(vn, 0, vd, 2, 1.0);
    *pos_tf = cascade_tf_create(pn, 0, pd, 1, 1.0);
    return 0;
}

int cascade_create_reactor_model(const CascadeReactor *rx,
                                  CascadeTF *jacket_tf,
                                  CascadeTF *reactor_tf) {
    /* Jacketed CSTR linearized TFs:
     * Inner (jacket): F_j -> T_j, tau_j = V_j/F_j
     * Outer (reactor): T_j -> T, tau_r = V_r/F_in
     * Reference: Luyben (2007), "Chemical Reactor Design and Control" */
    if (!rx || !jacket_tf || !reactor_tf) return -1;
    double Fjs = (rx->F_j > 1e-10) ? rx->F_j : 1e-10;
    double tau_j = rx->V_j / Fjs;
    double K_j = fabs(rx->T_j_in - 300.0) / (rx->rho * rx->Cp * Fjs + 1e-10);
    if (K_j > 10.0) K_j = 1.0;
    double jn[] = {K_j}, jd[] = {1.0, tau_j};
    *jacket_tf = cascade_tf_create(jn, 0, jd, 1, 1.0);
    double Fis = (rx->F_in > 1e-10) ? rx->F_in : 1e-10;
    double tau_r = rx->V_r / Fis;
    double K_r = rx->UA / (rx->UA + rx->rho * rx->Cp * Fis + 1e-10);
    double rn[] = {K_r}, rd[] = {1.0, tau_r};
    *reactor_tf = cascade_tf_create(rn, 0, rd, 1, 1.0);
    return 0;
}

int cascade_create_flow_pressure_model(const CascadeFlowPressure *fp,
                                        CascadeTF *pressure_tf,
                                        CascadeTF *flow_tf) {
    /* Flow-pressure cascade for fluid systems.
     * Inner (pressure): fast acoustic dynamics
     * Outer (flow): quasi-static valve equation */
    if (!fp || !pressure_tf || !flow_tf) return -1;
    double tau_p = fp->pipe_length / 340.0;
    if (tau_p < 0.01) tau_p = 0.1;
    double K_p = fp->pump_head / 10.0;
    if (K_p < 0.01) K_p = 1.0;
    double pn[] = {K_p}, pd[] = {1.0, tau_p};
    *pressure_tf = cascade_tf_create(pn, 0, pd, 1, 1.0);
    double dP = fp->pump_head * fp->fluid_density * 9.81;
    if (dP < 100.0) dP = 100000.0;
    double K_f = fp->valve_Cv / (2.0 * sqrt(dP));
    if (K_f < 1e-10) K_f = 0.01;
    double fn[] = {K_f}, fd[] = {1.0};
    *flow_tf = cascade_tf_create(fn, 0, fd, 0, 1.0);
    return 0;
}

int cascade_create_level_tank_model(const CascadeLevelTank *tank,
                                     CascadeTF *flow_tf,
                                     CascadeTF *level_tf) {
    /* Level-on-flow cascade.
     * Inner (flow): valve -> outflow, first-order
     * Outer (level): 1/(A*s), integrating process */
    if (!tank || !flow_tf || !level_tf) return -1;
    double tau_v = 0.5;
    double K_v = tank->pump_max_flow / 100.0;
    if (K_v < 1e-10) K_v = 0.01;
    double fn[] = {K_v}, fd[] = {1.0, tau_v};
    *flow_tf = cascade_tf_create(fn, 0, fd, 1, 1.0);
    double A = (tank->tank_area > 0.01) ? tank->tank_area : 1.0;
    double Kl = 1.0 / A;
    double ln[] = {Kl}, ld[] = {0.0, 1.0};
    *level_tf = cascade_tf_create(ln, 0, ld, 1, 1.0);
    return 0;
}
