/**
 * @file cascade_analysis.h
 * @brief Stability and frequency-domain analysis of cascade control systems
 *
 * L3 --- Mathematical Structures: Transfer function algebra for cascaded loops.
 * L4 --- Fundamental Laws: Internal stability conditions for cascade systems.
 * L5 --- Computational Methods: Frequency response computation, stability margins.
 *
 * Core mathematical structure:
 *   Inner closed loop:   Gi_cl = Ci*Gi / (1 + Ci*Gi)
 *   Equivalent plant:    Geq = Gi_cl * Go
 *   Outer open loop:     Lo = Co * Geq
 *   Overall closed loop: G_cl = Lo / (1 + Lo) = Co*Gi_cl*Go / (1 + Co*Gi_cl*Go)
 *
 * Internal stability requires all 4 transfer functions (from [r,d1,d2] to
 * [y,u_i,u_o]) to be stable. For cascade, this adds constraints beyond
 * simple single-loop Nyquist.
 */

#ifndef CASCADE_ANALYSIS_H
#define CASCADE_ANALYSIS_H

#include "cascade_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L3: Transfer Function Operations for Cascade Systems
 * ========================================================================== */

/** Compute frequency response G(jw) of a transfer function.
 *  Evaluates the complex rational function at s = j*omega.
 *  Complexity: O(deg_num + deg_den) using Horner's method.
 *
 *  @param tf    Transfer function G(s)
 *  @param omega Frequency (rad/s)
 *  @param mag   [out] |G(jw)| (linear magnitude)
 *  @param phase [out] arg(G(jw)) in radians, range [-pi, pi]
 *  @return 0 on success, -1 if den(jw) == 0 (pole on jw-axis) */
int cascade_tf_freq_response(const CascadeTF *tf, double omega,
                              double *mag, double *phase);

/** Compute closed-loop transfer function from open-loop.
 *  G_cl = L / (1 + L)  for unity negative feedback.
 *  Polynomial method: num_cl = num_L, den_cl = den_L + num_L
 *
 *  @param L   Open-loop transfer function L(s)
 *  @param cl  [out] Closed-loop transfer function
 *  @return 0 on success */
int cascade_tf_close_loop(const CascadeTF *L, CascadeTF *cl);

/** Form equivalent plant for outer loop: Geq(s) = Gi_cl(s) * Go(s).
 *  This is the transfer function seen by the outer controller.
 *  Key step in sequential loop closure design.
 *
 *  @param Gi        Inner process TF
 *  @param Ci        Inner controller PID
 *  @param Go        Outer process TF
 *  @param eq_plant  [out] Equivalent plant Geq(s)
 *  @return 0 on success */
int cascade_form_equivalent_plant(const CascadeTF *Gi,
                                   const CascadePID *Ci,
                                   const CascadeTF *Go,
                                   CascadeTF *eq_plant);

/** Compute inner closed-loop transfer function: Gi_cl = Ci*Gi/(1+Ci*Gi).
 *  Converts PID controller to transfer function form first.
 *
 *  @param Gi      Inner process TF
 *  @param Ci      Inner PID controller
 *  @param Gi_cl   [out] Inner closed-loop TF
 *  @return 0 on success */
int cascade_inner_closed_loop(const CascadeTF *Gi,
                               const CascadePID *Ci,
                               CascadeTF *Gi_cl);

/** Compute overall cascade closed-loop: y/r = Co*Gi_cl*Go/(1+Co*Gi_cl*Go).
 *
 *  @param Gi_cl    Inner closed-loop TF
 *  @param Co       Outer PID controller
 *  @param Go       Outer process TF
 *  @param overall  [out] Overall closed-loop TF from r to y
 *  @return 0 on success */
int cascade_overall_closed_loop(const CascadeTF *Gi_cl,
                                 const CascadePID *Co,
                                 const CascadeTF *Go,
                                 CascadeTF *overall);

/* ==========================================================================
 * L3: Polynomial Operations for Transfer Function Algebra
 * ========================================================================== */

/** Create a polynomial from coefficient array.
 *  @param coeff  Coefficient array [a0, a1, ..., an]
 *  @param degree Polynomial degree n
 *  @return Allocated polynomial (caller must cascade_poly_free) */
CascadePoly cascade_poly_create(const double *coeff, int degree);

/** Free polynomial memory. */
void cascade_poly_free(CascadePoly *p);

/** Evaluate polynomial p(s) at complex value s = sigma + j*omega.
 *  @param p      Polynomial
 *  @param sigma  Real part of evaluation point
 *  @param omega  Imaginary part of evaluation point
 *  @param re     [out] Real part of p(sigma + j*omega)
 *  @param im     [out] Imaginary part of p(sigma + j*omega) */
void cascade_poly_eval_complex(const CascadePoly *p,
                                double sigma, double omega,
                                double *re, double *im);

/** Multiply two polynomials: r = p * q.
 *  Complexity: O(deg(p) * deg(q)).
 *  @param p, q  Input polynomials
 *  @return Product polynomial (caller must free) */
CascadePoly cascade_poly_mul(const CascadePoly *p, const CascadePoly *q);

/** Add two polynomials: r = p + q.
 *  @param p, q  Input polynomials
 *  @return Sum polynomial (caller must free) */
CascadePoly cascade_poly_add(const CascadePoly *p, const CascadePoly *q);

/* ==========================================================================
 * L4: Internal Stability Analysis
 * ========================================================================== */

/** Verify internal stability of cascaded system.
 *
 *  Theorem (Internal Stability of Cascade): The cascade system is
 *  internally stable iff the following 6 transfer functions are stable:
 *    (1) y/r, (2) y/d1, (3) y/d2
 *    (4) u_i/r, (5) u_i/d1, (6) u_i/d2
 *  where d1 disturbs outer PV, d2 disturbs inner PV.
 *
 *  Practically, this reduces to checking that:
 *    - Inner closed loop is stable (poles of Gi_cl in LHP)
 *    - Outer closed loop is stable (poles of overall CL in LHP)
 *    - No unstable pole-zero cancellation between Ci and Gi
 *    - No unstable pole-zero cancellation between Co and Geq
 *
 *  @param Gi    Inner process TF
 *  @param Ci    Inner controller
 *  @param Go    Outer process TF
 *  @param Co    Outer controller
 *  @return 1 if internally stable, 0 otherwise */
int cascade_verify_internal_stability(const CascadeTF *Gi,
                                       const CascadePID *Ci,
                                       const CascadeTF *Go,
                                       const CascadePID *Co);

/** Check if all poles of G(s) are in the open left half-plane.
 *  For continuous-time stability: Re(p_i) < 0 for all i.
 *  Uses the Routh-Hurwitz criterion for polynomial root location.
 *
 *  @param tf  Transfer function
 *  @return 1 if all poles stable, 0 otherwise */
int cascade_is_stable(const CascadeTF *tf);

/** Compute the Routh array for a polynomial and count RHP roots.
 *  Theorem (Routh-Hurwitz): Number of RHP roots equals sign changes
 *  in the first column of the Routh array.
 *
 *  @param den    Denominator polynomial coefficients [a0..an]
 *  @param degree Polynomial degree
 *  @param n_rhp  [out] Number of roots in right half-plane
 *  @return 0 on success, -1 if degenerate case encountered */
int cascade_routh_hurwitz(const double *den, int degree, int *n_rhp);

/* ==========================================================================
 * L5: Frequency-Domain Analysis Methods
 * ========================================================================== */

/** Compute complete Bode frequency response over specified range.
 *  Evaluates G(jw) at logarithmically spaced frequencies.
 *
 *  @param tf        Transfer function
 *  @param freq_min  Minimum frequency (rad/s)
 *  @param freq_max  Maximum frequency (rad/s)
 *  @param n_points  Number of frequency points
 *  @param resp      [out] Frequency response data (caller must free points)
 *  @return 0 on success */
int cascade_bode_analysis(const CascadeTF *tf,
                           double freq_min, double freq_max,
                           int n_points,
                           CascadeFreqResponse *resp);

/** Compute gain and phase margins from frequency response.
 *  Gain margin: GM = -20*log10(|L(jw_pc)|) where arg(L(jw_pc)) = -180 deg
 *  Phase margin: PM = 180 + arg(L(jw_gc)) where |L(jw_gc)| = 1
 *
 *  @param resp  Frequency response of open-loop L(s)
 *  @param gm    [out] Gain margin in dB
 *  @param pm    [out] Phase margin in degrees
 *  @return 0 on success, -1 if margins undefined */
int cascade_stability_margins(const CascadeFreqResponse *resp,
                               double *gm, double *pm);

/** Compute maximum sensitivity Ms = max_w |1/(1+L(jw))|.
 *  Ms <= 2 ensures gain margin >= 6 dB and phase margin >= 30 deg.
 *  Ms <= 1.6 ensures gain margin >= 9 dB and phase margin >= 45 deg.
 *
 *  @param resp  Frequency response of open-loop L(s)
 *  @return Ms value (>= 1) */
double cascade_max_sensitivity(const CascadeFreqResponse *resp);

/** Free frequency response data.
 *  @param resp  Response data to free */
void cascade_freq_response_free(CascadeFreqResponse *resp);

/** Compute overall cascade open-loop frequency response.
 *  L_cascade(jw) = Co(jw) * Gi_cl(jw) * Go(jw)
 *
 *  @param Gi      Inner process
 *  @param Ci      Inner controller
 *  @param Go      Outer process
 *  @param Co      Outer controller
 *  @param omega   Frequency (rad/s)
 *  @param mag     [out] |L_cascade(jw)|
 *  @param phase   [out] arg(L_cascade(jw)) in radians
 *  @return 0 on success */
int cascade_loop_freq_response(const CascadeTF *Gi,
                                const CascadePID *Ci,
                                const CascadeTF *Go,
                                const CascadePID *Co,
                                double omega,
                                double *mag, double *phase);

/** Compute bandwidth (frequency where closed-loop magnitude = -3 dB).
 *  Solves |G_cl(jw)| = 1/sqrt(2) via bisection.
 *
 *  @param cl_tf  Closed-loop transfer function
 *  @param bw     [out] -3 dB bandwidth (rad/s)
 *  @return 0 on success */
int cascade_closed_loop_bandwidth(const CascadeTF *cl_tf, double *bw);

/** Free a CascadeTF structure (frees internal polynomial arrays).
 *  @param tf  Transfer function to free */
void cascade_tf_free(CascadeTF *tf);

/** Create a transfer function from coefficient arrays.
 *  @param num_coeff  Numerator coefficients [b0, b1, ..., bm]
 *  @param num_deg    Numerator degree m
 *  @param den_coeff  Denominator coefficients [a0, a1, ..., an]
 *  @param den_deg    Denominator degree n
 *  @param gain       Overall gain K
 *  @return Transfer function (caller must cascade_tf_free) */
CascadeTF cascade_tf_create(const double *num_coeff, int num_deg,
                              const double *den_coeff, int den_deg,
                              double gain);

/** Convert PID controller to transfer function form.
 *  C(s) = Kp + Ki/s + Kd*s/(N*s+1)
 *       = (Kp*N*s^2 + (Kp + Ki*N + Kd)*s + Ki) / (N*s^2 + s)
 *
 *  @param pid  PID controller
 *  @param tf   [out] Transfer function C(s)
 *  @return 0 on success */
int cascade_pid_to_tf(const CascadePID *pid, CascadeTF *tf);

/** Find dominant poles of a transfer function via companion matrix method.
 *  Uses the eigenvalues of the companion matrix (Francis, 1961).
 *  For continuous-time: Re(pole) < 0 means stable.
 *
 *  @param tf      Transfer function
 *  @param poles   [out] Array of complex poles (2*degree doubles: re,im pairs)
 *  @param max_n   Maximum number of poles to find
 *  @return Number of poles found (equal to deg(den)), -1 on error */
int cascade_find_poles(const CascadeTF *tf, double *poles, int max_n);

/** Compute step response of a transfer function via inverse Laplace.
 *  Uses partial fraction expansion for rational TFs.
 *  For TFs with delay, includes Pade approximation.
 *
 *  @param tf       Transfer function
 *  @param t        Time point (s)
 *  @param n_terms  Number of terms in partial fraction expansion
 *  @return Step response value y(t) */
double cascade_step_response(const CascadeTF *tf, double t, int n_terms);

/** Simulate step response over a time vector.
 *  @param tf       Transfer function
 *  @param t        Time array (s)
 *  @param n        Number of time points
 *  @param y        [out] Output response (must be pre-allocated, size n)
 *  @return 0 on success */
int cascade_step_response_vector(const CascadeTF *tf,
                                  const double *t, int n,
                                  double *y);

/* ==========================================================================
 * L1: System Lifecycle
 * ========================================================================== */

/** Initialize a cascade system with safe default parameters.
 *  @param sys  Cascade system to initialize
 *  @param name Human-readable system name */
void cascade_system_init(CascadeSystem *sys, const char *name);

/** Free all memory associated with a cascade system.
 *  @param sys  Cascade system to free */
void cascade_system_free(CascadeSystem *sys);

/** Set default design specifications.
 *  @param spec  Specification struct to populate */
void cascade_set_default_spec(CascadeDesignSpec *spec);

/* ==========================================================================
 * L6: Canonical System Model Creation
 * ========================================================================== */

int cascade_create_dc_motor_model(const CascadeDCMotor *motor,
                                   CascadeTF *vel_tf, CascadeTF *pos_tf);

int cascade_create_reactor_model(const CascadeReactor *rx,
                                  CascadeTF *jacket_tf, CascadeTF *reactor_tf);

int cascade_create_flow_pressure_model(const CascadeFlowPressure *fp,
                                        CascadeTF *pressure_tf,
                                        CascadeTF *flow_tf);

int cascade_create_level_tank_model(const CascadeLevelTank *tank,
                                     CascadeTF *flow_tf, CascadeTF *level_tf);

/* ==========================================================================
 * L5: Performance Metrics Computation
 * ========================================================================== */

int cascade_compute_performance(const CascadeTF *cl_tf,
                                 double t_final, int n_pts,
                                 CascadePerformance *perf);

int cascade_compare_single_vs_cascade(const CascadeTF *inner_model,
                                       const CascadeTF *outer_model,
                                       const CascadePID *ci,
                                       const CascadePID *co,
                                       double t_final, int n_pts,
                                       double *improvement_factor);

double cascade_compute_ise(const double *error, const double *time, int n);
double cascade_compute_iae(const double *error, const double *time, int n);
double cascade_compute_itae(const double *error, const double *time, int n);
double cascade_compute_tv(const double *control, int n);

int cascade_extract_rise_time(const double *t, const double *y, int n,
                               double y_final, double *rise_time);
int cascade_extract_settling_time(const double *t, const double *y, int n,
                                   double y_final, double band_pct,
                                   double *settle_time);
double cascade_extract_overshoot(const double *y, int n, double y_final);

/* ==========================================================================
 * L7: Application Assessment
 * ========================================================================== */

int cascade_assess_dc_motor(const CascadeDCMotor *motor,
                             CascadePerformance *perf);
int cascade_assess_reactor(const CascadeReactor *rx,
                            CascadePerformance *perf);
int cascade_assess_level_tank(const CascadeLevelTank *tank,
                               CascadePerformance *perf);

#ifdef __cplusplus
}
#endif

#endif /* CASCADE_ANALYSIS_H */
