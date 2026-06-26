/**
 * @file feedforward_design.h
 * @brief Feedforward controller design: model inversion, ZPETC,
 *        Diophantine design, properness enforcement.
 *
 * L3 --- Mathematical structures: transfer function inversion, polynomial
 * long division, Diophantine equations for feedforward synthesis.
 *
 * L4 --- Fundamental theorems: perfect feedforward condition,
 * internal model principle, causality constraints.
 *
 * L5 --- Computational methods: model inverse approximation,
 * zero-phase error tracking, proper inverse construction.
 *
 * Course alignment:
 *   MIT 6.245 Multivariable Control ? model inversion
 *   Stanford ENGR210B Robust Control ? Youla parameterization + FF
 *   Berkeley ME234 Nonlinear ? feedforward linearization
 *   ETH 151-0563 Robust ? 2-DOF design
 *   Caltech CDS 212 ? robust feedforward
 */

#ifndef FEEDFORWARD_DESIGN_H
#define FEEDFORWARD_DESIGN_H

#include "feedforward_core.h"

/* ==========================================================================
 * L4: Perfect Feedforward Theorem and Model Inversion
 * ========================================================================== */

int ff_model_inverse(const TransferFn *plant, TransferFn *inv,
                     int filter_order, double filter_bw);

int ff_causal_inverse(const TransferFn *plant, TransferFn *inv,
                      int max_order);

int ff_zpetc_design(const TransferFn *plant_dt, TransferFn *zp_ff,
                    double Ts);

int ff_diophantine_design(const TransferFn *plant,
                          const TransferFn *desired,
                          Poly *ff_num, Poly *ff_den);

int ff_enforce_properness(TransferFn *tf, double filter_bw);

/* ==========================================================================
 * L5: Reference Feedforward Design
 * ========================================================================== */

int ff_design_ref_feedforward(const TransferFn *plant,
                              const TransferFn *fb_ctrl,
                              const TransferFn *desired,
                              TransferFn *ff_ref);

void ff_design_prefilter_1st(double rise_time, Prefilter *prefilter);

void ff_design_prefilter_2nd(double rise_time, double damping,
                             Prefilter *prefilter);

/* ==========================================================================
 * L5: Disturbance Feedforward Design
 * ========================================================================== */

int ff_design_static_dist_ff(const TransferFn *plant,
                             const TransferFn *dist_model,
                             double *ff_gain);

int ff_design_dynamic_dist_ff(const TransferFn *plant,
                              const TransferFn *dist_model,
                              TransferFn *comp,
                              double filter_bw);

int ff_eval_dist_rejection(const TransferFn *plant,
                           const TransferFn *fb_ctrl,
                           const TransferFn *dist_ff,
                           const TransferFn *dist_model,
                           double omega, double *attenuation);

#endif /* FEEDFORWARD_DESIGN_H */
