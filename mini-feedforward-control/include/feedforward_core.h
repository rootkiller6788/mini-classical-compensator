/**
 * @file feedforward_core.h
 * @brief Core feedforward control definitions and elementary operations
 *
 * L1 --- Definitions: Feedforward controller, model inverse, 2-DOF structure,
 * prefilter, disturbance feedforward path. Feedforward control is the
 * second degree of freedom in a 2-DOF architecture, operating on reference
 * or measured disturbance signals to improve tracking/rejection beyond
 * what feedback alone can achieve.
 *
 * Key insight: feedback reacts to errors (past), feedforward anticipates
 * them (future). Together they form the complementary basis of modern
 * control design.
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems — feedforward compensation
 *   Stanford ENGR105 — 2-DOF control architecture
 *   Berkeley ME232 — feedforward in motion control
 *   Cambridge 3F2 — disturbance feedforward
 *   ETH 151-0591 — Störgrößenaufschaltung
 *   Tsinghua 自动控制原理 — 前馈控制
 *
 * Textbook: Åström & Hägglund, "Advanced PID Control" (2006)
 *           Goodwin, Graebe & Salgado, "Control System Design" (2001)
 *           Skogestad & Postlethwaite, "Multivariable Feedback Control" (2005)
 */

#ifndef FEEDFORWARD_CORE_H
#define FEEDFORWARD_CORE_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Core Data Types
 * ========================================================================== */

/** Polynomial: p(s) = a0 + a1*s + a2*s^2 + ... + aN*s^N
 *  Represented as array [a0, a1, ..., aN] of length N+1. */
typedef struct {
    int      degree;
    double  *coeff;
} Poly;

/** Transfer function: G(s) = num(s) / den(s) */
typedef struct {
    Poly  num;
    Poly  den;
    double gain;
} TransferFn;

/** Feedforward controller structure.
 *  In 2-DOF: u = C(s)*(r-y) + F(s)*r (reference FF)
 *          or u = C(s)*(r-y) + D(s)*d_m (disturbance FF)
 */
typedef struct {
    TransferFn  ff_filter;
    TransferFn  plant_model;
    TransferFn  plant_inverse;
    double      dc_gain;
    int         is_causal;
    int         is_stable;
} FeedforwardCtrl;

/** Two-Degree-of-Freedom (2-DOF) controller. */
typedef struct {
    TransferFn   feedback_ctrl;
    TransferFn   ref_ff;
    TransferFn   dist_ff;
    TransferFn   plant;
    double        ff_gain_ratio;
    double        coupling_factor;
} TwoDOF;

/** Prefilter / reference shaping filter. */
typedef struct {
    TransferFn  filter;
    double      rise_time;
    double      bandwidth;
    int         order;
    double      damping;
} Prefilter;

/** Disturbance model for feedforward compensation. */
typedef struct {
    TransferFn   dist_model;
    TransferFn   comp_filter;
    double        meas_delay;
    double        comp_gain;
    int           is_measurable;
} DistFFModel;

/** Feedforward performance metrics. */
typedef struct {
    double   tracking_ise;
    double   disturbance_ise;
    double   ff_contribution;
    double   settle_time;
    double   overshoot_pct;
    double   robustness_margin;
} FFPerformance;

/* ==========================================================================
 * L2: Core Feedforward Operations
 * ========================================================================== */

int ff_compute_ideal_filter(const TransferFn *plant,
                            const TransferFn *desired,
                            TransferFn *ff);

int ff_compute_dist_compensator(const TransferFn *plant,
                                const TransferFn *dist,
                                TransferFn *comp);

int ff_compute_closed_loop_tf(const TransferFn *plant,
                               const TransferFn *fb_ctrl,
                               const TransferFn *ff_ref,
                               TransferFn *cl_tf);

int ff_is_realizable(const TransferFn *ff);

double ff_dc_gain(const TransferFn *ff, int is_dt);

void ff_init(const TransferFn *plant, FeedforwardCtrl *ff);

void ff_free(FeedforwardCtrl *ff);

/* ==========================================================================
 * L2: Polynomial Operations
 * ========================================================================== */

Poly poly_create(const double *coeff, int degree);
void poly_free(Poly *p);
double poly_eval(const Poly *p, double x);
Poly poly_mul(const Poly *p, const Poly *q);
Poly poly_add(const Poly *p, const Poly *q);
Poly poly_sub(const Poly *p, const Poly *q);
Poly poly_derivative(const Poly *p);

/* ==========================================================================
 * L3: Transfer Function Analysis
 * ========================================================================== */

void tf_freq_response(const TransferFn *tf, double omega,
                      double *mag, double *phase);

double tf_step_response(const TransferFn *tf, double t, int n_approx);

int tf_find_poles(const TransferFn *tf, double *poles, int n);

int tf_is_minimum_phase(const TransferFn *tf, int is_dt);

int tf_to_state_space(const TransferFn *tf,
                      double *A, double *B, double *C, double *D);

void tf_free(TransferFn *tf);

TransferFn tf_create(const double *num, int num_deg,
                     const double *den, int den_deg,
                     double gain);

/* ==========================================================================
 * L3: Transfer Function Algebra (implemented in feedforward_core.c)
 * ========================================================================== */

TransferFn tf_parallel(const TransferFn *g1, const TransferFn *g2);
TransferFn tf_series(const TransferFn *g1, const TransferFn *g2);
TransferFn tf_feedback(const TransferFn *g, const TransferFn *h);
TransferFn tf_negate(const TransferFn *tf);
TransferFn tf_scale(const TransferFn *tf, double k);
int tf_equals(const TransferFn *a, const TransferFn *b, double tol);

#endif /* FEEDFORWARD_CORE_H */
