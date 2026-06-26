/**
 * @file lag_optimization.c
 * @brief Parameter optimization for lag compensators
 *
 * Finds optimal lag compensator parameters (Kc, T, beta) to minimize
 * a cost function combining steady-state error, phase margin loss,
 * and bandwidth reduction.
 *
 * L5: Numerical optimization — golden-section search, gradient descent
 * L8: Multi-objective optimization for compensator tuning
 *
 * Textbook: Ogata Ch. 7; Nocedal & Wright, "Numerical Optimization" (2006)
 */

#include "lag_compensator.h"
#include "lag_types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ==========================================================================
 * L5: Objective function for lag compensator optimization
 * ========================================================================== */

/**
 * Computes a composite cost for a lag compensator design.
 *
 * The cost function balances three competing objectives:
 *   J1: Steady-state error (lower is better)
 *        J1 = K_required / (Kc * plant_dc_gain) - 1
 *        J1 >= 0, with J1 = 0 when ESS requirement is exactly met
 *
 *   J2: Phase margin reduction (lower is better)
 *        J2 = max(0, PM_req - PM_achieved) / PM_req
 *        Penalizes designs that sacrifice too much phase margin
 *
 *   J3: Bandwidth reduction penalty
 *        J3 = max(0, w_bw_min - w_bw_achieved) / w_bw_min
 *        Penalizes designs that reduce bandwidth too much
 *
 * Total cost: J = w1*J1 + w2*J2 + w3*J3
 *
 * Weights default: w1=1.0, w2=0.5, w3=0.3
 * (ESS is typically the primary design driver for lag compensation)
 *
 * Complexity: O(1)
 *
 * @param Kc       compensator gain to evaluate
 * @param T        time constant to evaluate
 * @param beta     separation factor to evaluate
 * @param plant_dc_gain  plant DC gain
 * @param w_gc     current gain crossover frequency
 * @param K_req    required error constant
 * @param PM_req   required phase margin (degrees)
 * @param w_bw_min minimum acceptable bandwidth
 * @return composite cost (lower is better)
 */
static double lag_cost_function(double Kc, double T, double beta,
                                 double plant_dc_gain, double w_gc,
                                 double K_req, double PM_req,
                                 double w_bw_min) {
    double total_cost = 0.0;

    /* J1: Steady-state error cost */
    double K_achieved = Kc * plant_dc_gain;
    double J1;
    if (K_achieved >= K_req) {
        J1 = 0.0;  /* ESS requirement met */
    } else {
        J1 = (K_req / K_achieved) - 1.0;
    }
    total_cost += 1.0 * J1;

    /* J2: Phase margin cost
     *
     * Estimate phase margin of compensated system.
     * phi_c at w_gc = arctan(w_gc*T) - arctan(w_gc*beta*T)
     * PM_achieved ~= PM_uncomp + phi_c (phi_c is negative for lag)
     *
     * For a typical system with PM_uncomp ~= 90 deg,
     * PM_achieved = 90 + phi_c_deg
     */
    double wgc_T = w_gc * T;
    double phase_contrib_rad = atan(wgc_T) - atan(wgc_T * beta);
    double phase_contrib_deg = phase_contrib_rad * 180.0 / M_PI;
    double PM_achieved = 90.0 + phase_contrib_deg;  /* assumes PM_uncomp=90 */
    if (PM_achieved < 0.0) PM_achieved = 0.0;

    double J2;
    if (PM_achieved >= PM_req) {
        J2 = 0.0;
    } else {
        J2 = (PM_req - PM_achieved) / PM_req;
    }
    total_cost += 0.5 * J2;

    /* J3: Bandwidth cost
     *
     * Estimate new bandwidth: w_bw_new ~= w_gc_new
     * For a lag compensator with Kc=beta, |G_c(jw_gc)| ~= 1,
     * so w_gc doesn't shift much. The bandwidth ~= w_gc.
     *
     * More precisely, when |G_c(jw)| has settled to its HF value:
     * w_bw_new ~= w_gc (approximately preserved if T is chosen well) */
    double w_bw_achieved = w_gc;  /* approximate */
    double J3;
    if (w_bw_achieved >= w_bw_min) {
        J3 = 0.0;
    } else {
        J3 = (w_bw_min - w_bw_achieved) / w_bw_min;
    }
    total_cost += 0.3 * J3;

    return total_cost;
}

/* ==========================================================================
 * L5: Golden-section search for optimal beta
 * ========================================================================== */

/**
 * Optimize beta parameter using golden-section search.
 *
 * Given Kc and T fixed, find beta in [beta_min, beta_max]
 * that minimizes the cost function.
 *
 * Golden-section search is an efficient method for unimodal
 * functions that does not require derivatives.
 * Convergence rate: each iteration reduces interval by factor ~0.618.
 *
 * Complexity: O(log(1/epsilon)) function evaluations
 *
 * @param Kc       fixed compensator gain
 * @param T        fixed time constant
 * @param beta_min minimum beta to consider (> 1)
 * @param beta_max maximum beta to consider
 * @param tol      convergence tolerance
 * @param params   cost function parameters
 * @return optimal beta value
 */
static double golden_section_beta(double Kc, double T,
                                   double beta_min, double beta_max,
                                   double tol,
                                   double plant_dc_gain, double w_gc,
                                   double K_req, double PM_req,
                                   double w_bw_min, int *iterations) {
    const double phi = (sqrt(5.0) - 1.0) / 2.0;  /* golden ratio conjugate */
    const int max_iter = 100;

    double a = beta_min;
    double b = beta_max;
    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    double fc = lag_cost_function(Kc, T, c, plant_dc_gain, w_gc,
                                   K_req, PM_req, w_bw_min);
    double fd = lag_cost_function(Kc, T, d, plant_dc_gain, w_gc,
                                   K_req, PM_req, w_bw_min);

    int iter;
    for (iter = 0; iter < max_iter && (b - a) > tol; iter++) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = lag_cost_function(Kc, T, c, plant_dc_gain, w_gc,
                                    K_req, PM_req, w_bw_min);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = lag_cost_function(Kc, T, d, plant_dc_gain, w_gc,
                                    K_req, PM_req, w_bw_min);
        }
    }

    *iterations = iter;
    return (a + b) / 2.0;
}

/* ==========================================================================
 * L5: Gradient descent for (Kc, T) optimization
 * ========================================================================== */

/**
 * Optimize Kc and T parameters simultaneously using gradient descent.
 *
 * While beta is optimized via golden-section search, Kc and T
 * are optimized via gradient descent with finite-difference gradients.
 *
 * Cost function J(Kc, T) is evaluated with beta fixed at its current optimum.
 *
 * Algorithm:
 *   1. Initialize Kc = K_req / plant_dc_gain, T = 10/w_gc
 *   2. For each iteration:
 *      a. Find optimal beta given (Kc, T) via golden-section search
 *      b. Compute gradient of J w.r.t. (Kc, T) using central differences
 *      c. Update: Kc <- Kc - alpha * dJ/dKc, T <- T - alpha * dJ/dT
 *   3. Stop when gradient norm < tolerance or max iterations reached
 *
 * Complexity: O(max_iter * log(1/tol_beta) * 6) function evaluations
 *
 * @param plant_dc_gain  plant DC gain
 * @param w_gc           current gain crossover frequency
 * @param spec           design specification
 * @param[out] result    optimized compensator and convergence info
 * @return 0 on success, negative on error
 */
int lag_optimize_parameters(double plant_dc_gain, double w_gc,
                             const LagDesignSpec *spec,
                             LagOptResult *result) {
    if (!spec || !result || plant_dc_gain <= 0 || w_gc <= 0) return -1;

    /* Initialize */
    double K_req = 1.0 / spec->ess_target;
    double Kc = K_req / plant_dc_gain;
    if (Kc < 1.0) Kc = 1.0;

    double T = 10.0 / w_gc;  /* standard decade placement */
    double beta_opt = Kc;    /* initial guess: beta = Kc */

    /* Gradient descent parameters */
    const double alpha = 0.01;   /* learning rate */
    const double grad_tol = 1e-6;/* gradient tolerance */
    const double delta = 1e-6;   /* finite difference step */
    const int max_iter = 500;

    double best_cost = INFINITY;
    double best_Kc = Kc, best_T = T, best_beta = beta_opt;
    int converged = 0;
    int iter;

    for (iter = 0; iter < max_iter; iter++) {
        /* Step a: optimize beta for current (Kc, T) */
        int beta_iters;
        double beta_min = 1.01;
        double beta_max = 100.0;
        beta_opt = golden_section_beta(Kc, T, beta_min, beta_max, 1e-4,
                                        plant_dc_gain, w_gc,
                                        K_req, spec->phase_margin_target,
                                        spec->bandwidth_min, &beta_iters);

        double cost = lag_cost_function(Kc, T, beta_opt,
                                         plant_dc_gain, w_gc,
                                         K_req, spec->phase_margin_target,
                                         spec->bandwidth_min);

        /* Track best */
        if (cost < best_cost) {
            best_cost = cost;
            best_Kc = Kc;
            best_T = T;
            best_beta = beta_opt;
        }

        /* Check convergence */
        if (cost < 1e-8) {
            converged = 1;
            break;
        }

        /* Step b: compute gradient via central differences */
        double cost_Kc_plus = lag_cost_function(Kc + delta, T, beta_opt,
                                                  plant_dc_gain, w_gc,
                                                  K_req, spec->phase_margin_target,
                                                  spec->bandwidth_min);
        double cost_Kc_minus = lag_cost_function(Kc - delta, T, beta_opt,
                                                   plant_dc_gain, w_gc,
                                                   K_req, spec->phase_margin_target,
                                                   spec->bandwidth_min);
        double dJ_dKc = (cost_Kc_plus - cost_Kc_minus) / (2.0 * delta);

        double cost_T_plus = lag_cost_function(Kc, T + delta, beta_opt,
                                                plant_dc_gain, w_gc,
                                                K_req, spec->phase_margin_target,
                                                spec->bandwidth_min);
        double cost_T_minus = lag_cost_function(Kc, T - delta, beta_opt,
                                                 plant_dc_gain, w_gc,
                                                 K_req, spec->phase_margin_target,
                                                 spec->bandwidth_min);
        double dJ_dT = (cost_T_plus - cost_T_minus) / (2.0 * delta);

        double grad_norm = sqrt(dJ_dKc * dJ_dKc + dJ_dT * dJ_dT);
        if (grad_norm < grad_tol) {
            converged = 1;
            break;
        }

        /* Step c: gradient descent update */
        Kc -= alpha * dJ_dKc;
        T  -= alpha * dJ_dT;

        /* Clamp to valid ranges */
        if (Kc < 0.01) Kc = 0.01;
        if (Kc > 1000.0) Kc = 1000.0;
        if (T < 1e-6) T = 1e-6;
        if (T > 1e6) T = 1e6;
    }

    /* Return results */
    result->optimal_Kc = best_Kc;
    result->optimal_T = best_T;
    result->optimal_beta = best_beta;
    result->objective_value = best_cost;
    result->iterations = iter;
    result->converged = converged;

    return 0;
}

/* ==========================================================================
 * L5: Grid search for robust initialization
 * ========================================================================== */

/**
 * Perform a coarse grid search to find a good initial guess
 * for gradient-based optimization.
 *
 * Evaluates cost on a 20x20 grid of (Kc, beta) values
 * and returns the best combination found.
 *
 * Complexity: O(400) cost function evaluations
 *
 * @param plant_dc_gain  plant DC gain
 * @param w_gc           gain crossover frequency
 * @param spec           design specification
 * @param[out] best_Kc   best Kc found
 * @param[out] best_beta best beta found
 * @param[out] best_cost cost at best point
 */
int lag_grid_search_initial(double plant_dc_gain, double w_gc,
                             const LagDesignSpec *spec,
                             double *best_Kc, double *best_beta,
                             double *best_cost) {
    if (!spec || !best_Kc || !best_beta || !best_cost) return -1;

    double K_req = 1.0 / spec->ess_target;
    double Kc0 = K_req / plant_dc_gain;
    if (Kc0 < 1.0) Kc0 = 1.0;

    *best_cost = INFINITY;
    *best_Kc = Kc0;
    *best_beta = Kc0;

    const int grid_size = 20;
    double T_fixed = 10.0 / w_gc;

    /* Search Kc from 0.5*Kc0 to 2*Kc0 */
    for (int i = 0; i < grid_size; i++) {
        double Kc = Kc0 * (0.5 + 1.5 * i / (grid_size - 1));

        /* Search beta from max(1.01, 0.5*Kc0) to 2*Kc0 */
        for (int j = 0; j < grid_size; j++) {
            double beta_min_g = Kc0 * 0.5;
            if (beta_min_g < 1.01) beta_min_g = 1.01;
            double beta = beta_min_g + (2.0*Kc0 - beta_min_g) * j / (grid_size - 1);

            double cost = lag_cost_function(Kc, T_fixed, beta,
                                             plant_dc_gain, w_gc,
                                             K_req, spec->phase_margin_target,
                                             spec->bandwidth_min);

            if (cost < *best_cost) {
                *best_cost = cost;
                *best_Kc = Kc;
                *best_beta = beta;
            }
        }
    }

    return 0;
}

/* ==========================================================================
 * L8: Multi-objective Pareto optimization
 * ========================================================================== */

/**
 * Compute the Pareto frontier for lag compensator design.
 *
 * A design (Kc, T, beta) is Pareto-optimal if no other design
 * improves one objective without worsening another.
 *
 * Objectives:
 *   f1 = steady-state error (minimize)
 *   f2 = phase margin reduction (minimize)
 *   f3 = settling time (minimize)
 *
 * This function samples the design space and identifies
 * non-dominated solutions.
 *
 * Complexity: O(n_samples^2) for dominance checking
 *
 * @param plant_dc_gain  plant DC gain
 * @param w_gc           gain crossover frequency
 * @param n_samples      number of design samples to evaluate
 * @param[out] n_pareto  number of Pareto-optimal designs found
 * @param[out] pareto_Kc array of Kc values for Pareto designs
 * @param[out] pareto_T  array of T values for Pareto designs
 * @param[out] pareto_beta array of beta values for Pareto designs
 * @return 0 on success, negative on error
 */
int lag_pareto_optimization(double plant_dc_gain, double w_gc,
                             int n_samples,
                             int *n_pareto,
                             double **pareto_Kc,
                             double **pareto_T,
                             double **pareto_beta) {
    if (n_samples < 10 || !n_pareto || !pareto_Kc || !pareto_T || !pareto_beta)
        return -1;

    /* Allocate sample arrays */
    double *Kc_samples = (double*)malloc((size_t)n_samples * sizeof(double));
    double *T_samples  = (double*)malloc((size_t)n_samples * sizeof(double));
    double *beta_samples = (double*)malloc((size_t)n_samples * sizeof(double));
    double *f1_samples = (double*)malloc((size_t)n_samples * sizeof(double));
    double *f2_samples = (double*)malloc((size_t)n_samples * sizeof(double));
    double *f3_samples = (double*)malloc((size_t)n_samples * sizeof(double));
    int *dominated = (int*)calloc((size_t)n_samples, sizeof(int));

    if (!Kc_samples || !T_samples || !beta_samples ||
        !f1_samples || !f2_samples || !f3_samples || !dominated) {
        free(Kc_samples); free(T_samples); free(beta_samples);
        free(f1_samples); free(f2_samples); free(f3_samples);
        free(dominated);
        return -2;
    }

    /* Generate samples using Latin Hypercube-inspired spacing */
    for (int i = 0; i < n_samples; i++) {
        double u = (i + 0.5) / n_samples;
        Kc_samples[i] = 1.0 + u * 99.0;     /* Kc in [1, 100] */
        T_samples[i] = 0.01 + u * 10.0;      /* T in [0.01, 10] */
        beta_samples[i] = 1.01 + u * 98.99;  /* beta in (1, 100] */
    }

    /* Shuffle for better coverage (Fisher-Yates) */
    for (int i = n_samples - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        double tmp;
        tmp = Kc_samples[i]; Kc_samples[i] = Kc_samples[j]; Kc_samples[j] = tmp;
        tmp = T_samples[i]; T_samples[i] = T_samples[j]; T_samples[j] = tmp;
        tmp = beta_samples[i]; beta_samples[i] = beta_samples[j]; beta_samples[j] = tmp;
    }

    /* Evaluate objectives */
    (void)plant_dc_gain;  /* used implicitly through K_req_temp */
    double K_req_temp = 50.0;  /* nominal */
    double PM_req_temp = 45.0;
    (void)PM_req_temp;  /* phase margin reference for Pareto ranking */
    for (int i = 0; i < n_samples; i++) {
        f1_samples[i] = 1.0 / (Kc_samples[i] * plant_dc_gain);  /* ~ess */
        double phase = atan(w_gc * T_samples[i]) -
                       atan(w_gc * beta_samples[i] * T_samples[i]);
        phase = phase * 180.0 / M_PI;
        f2_samples[i] = fabs(phase);  /* phase lag magnitude */
        f3_samples[i] = T_samples[i]; /* proxy for settling time */
    }

    /* Find Pareto-optimal designs (non-dominated sorting) */
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_samples; j++) {
            if (i == j || dominated[i]) continue;
            /* Check if j dominates i */
            if (f1_samples[j] <= f1_samples[i] &&
                f2_samples[j] <= f2_samples[i] &&
                f3_samples[j] <= f3_samples[i] &&
                (f1_samples[j] < f1_samples[i] ||
                 f2_samples[j] < f2_samples[i] ||
                 f3_samples[j] < f3_samples[i])) {
                dominated[i] = 1;
                break;
            }
        }
    }

    /* Count and collect Pareto designs */
    int count = 0;
    for (int i = 0; i < n_samples; i++) {
        if (!dominated[i]) count++;
    }

    *n_pareto = count;
    *pareto_Kc = (double*)malloc((size_t)count * sizeof(double));
    *pareto_T = (double*)malloc((size_t)count * sizeof(double));
    *pareto_beta = (double*)malloc((size_t)count * sizeof(double));

    if (!*pareto_Kc || !*pareto_T || !*pareto_beta) {
        free(*pareto_Kc); free(*pareto_T); free(*pareto_beta);
        free(Kc_samples); free(T_samples); free(beta_samples);
        free(f1_samples); free(f2_samples); free(f3_samples);
        free(dominated);
        return -2;
    }

    int idx = 0;
    for (int i = 0; i < n_samples; i++) {
        if (!dominated[i]) {
            (*pareto_Kc)[idx] = Kc_samples[i];
            (*pareto_T)[idx] = T_samples[i];
            (*pareto_beta)[idx] = beta_samples[i];
            idx++;
        }
    }

    /* Cleanup */
    free(Kc_samples); free(T_samples); free(beta_samples);
    free(f1_samples); free(f2_samples); free(f3_samples);
    free(dominated);

    return 0;
}