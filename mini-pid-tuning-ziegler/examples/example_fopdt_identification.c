/**
 * example_fopdt_identification.c — FOPDT Model Identification + PID Tuning
 *
 * Demonstrates the complete workflow:
 *   1. Generate synthetic process data (step response)
 *   2. Identify FOPDT model using multiple methods
 *   3. Compare identification accuracy
 *   4. Tune PID using the identified model
 *
 * Knowledge: L5 (identification methods), L6 (FOPDT benchmark),
 *             L7 (process automation workflow).
 */

#include "pid_tuning.h"
#include "fopdt_model.h"
#include "ziegler_nichols.h"
#include "cohen_coon.h"
#include "imc_tuning.h"
#include "advanced_tuning.h"
#include "application_tuning.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    printf("FOPDT Model Identification + PID Tuning Pipeline\n");
    printf("================================================\n\n");

    /* Step 1: Define the "true" unknown process */
    fopdt_model_t true_model = { 1.5, 8.0, 2.0 };
    printf("True process (unknown to tuning):\n");
    printf("  G(s) = %.2f * exp(-%.1f*s) / (%.1f*s + 1)\n\n",
           true_model.K, true_model.L, true_model.T);

    /* Step 2: Generate step response data (with added noise) */
    int N = 300;
    double Ts = 0.1;
    double *t = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        t[i] = i * Ts;
        u[i] = 1.0;
        if (t[i] < true_model.L) {
            y[i] = 0.0;
        } else {
            y[i] = true_model.K * (1.0 - exp(-(t[i] - true_model.L) / true_model.T));
        }
        /* Add 2% random noise */
        y[i] += 0.02 * true_model.K * ((double)rand() / RAND_MAX - 0.5);
    }

    step_response_data_t data;
    data.time     = t;
    data.input    = u;
    data.output   = y;
    data.N        = N;
    data.Ts       = Ts;
    data.step_mag = 1.0;

    /* Step 3: Identify using 4 different methods */
    printf("Identification Results:\n");
    printf("-----------------------\n");

    fopdt_id_result_t results[4];
    const char *method_names[] = {
        "Tangent (Z-N)    ",
        "Two-point (28/63)",
        "Area             ",
        "Sundaresan (35/85)"
    };

    /* Method 1: Graphical/Tangent */
    fopdt_identify_graphical(&data, &results[0]);

    /* Method 2: Two-point */
    fopdt_identify_two_point(&data, &results[1]);

    /* Method 3: Area */
    fopdt_identify_area(&data, &results[2]);

    /* Method 4: Sundaresan-Krishnaswamy */
    fopdt_identify_sundaresan(&data, &results[3]);

    printf("Method              K        T       L       Fit%%\n");
    printf("------------------  -------  ------  ------  -----\n");
    printf("True model          %.4f   %.2f   %.2f   100.0\n",
           true_model.K, true_model.T, true_model.L);

    int best_idx = 0;
    double best_fit = 0.0;

    for (int m = 0; m < 4; m++) {
        printf("%s   %.4f   %.2f    %.3f   %.1f%%\n",
               method_names[m],
               results[m].model.K,
               results[m].model.T,
               results[m].model.L,
               results[m].fit_pct);
        if (results[m].fit_pct > best_fit) {
            best_fit = results[m].fit_pct;
            best_idx = m;
        }
    }
    printf("\nBest identification: %s(R²=%.3f)\n\n",
           method_names[best_idx], results[best_idx].r_squared);

    /* Step 4: Also try Least-Squares */
    printf("Least-Squares Refinement:\n");
    fopdt_id_result_t ls_result;
    fopdt_identify_ls(&data, 0.5, 5.0, 20, 10, &ls_result);
    printf("  LS: K=%.4f, T=%.3f, L=%.3f, R²=%.4f\n",
           ls_result.model.K, ls_result.model.T,
           ls_result.model.L, ls_result.r_squared);

    /* Step 5: Compute model distance to true model */
    double dist_best = fopdt_distance(&true_model, &results[best_idx].model);
    double dist_ls   = fopdt_distance(&true_model, &ls_result.model);
    printf("  Distance to true: Best=%.4f, LS=%.4f\n\n", dist_best, dist_ls);

    /* Step 6: Tune PID using the best identified model */
    printf("PID Tuning using best identified model:\n");
    printf("---------------------------------------\n");

    fopdt_model_t *id_model = &results[best_idx].model;

    /* Z-N Step */
    pid_params_t zn_p;
    zn_step_tune(id_model, ZN_CONTROLLER_PID, &zn_p);
    printf("  Z-N Step:    Kp=%.3f, Ti=%.3f, Td=%.3f\n",
           zn_p.Kp, zn_p.Ti, zn_p.Td);

    /* Cohen-Coon */
    pid_params_t cc_p;
    cohen_coon_tune(id_model, CC_CONTROLLER_PID, &cc_p);
    printf("  Cohen-Coon:  Kp=%.3f, Ti=%.3f, Td=%.3f\n",
           cc_p.Kp, cc_p.Ti, cc_p.Td);

    /* SIMC */
    pid_params_t simc_p;
    imc_simc_tune(id_model, id_model->L, 1, &simc_p);
    printf("  SIMC:        Kp=%.3f, Ti=%.3f, Td=%.3f\n",
           simc_p.Kp, simc_p.Ti, simc_p.Td);

    /* Compare tuning quality using prediction */
    printf("\nPredicted closed-loop performance of each tuning:\n");
    pid_params_t *tunings[] = { &zn_p, &cc_p, &simc_p };
    const char *tune_names[] = { "Z-N Step", "Cohen-Coon", "SIMC" };

    for (int i = 0; i < 3; i++) {
        double decay, overshoot, omega;
        cohen_coon_predict_response(id_model, tunings[i],
                                     &decay, &overshoot, &omega);
        printf("  %-12s: overshoot=%.1f%%, decay ratio=%.3f, ω_osc=%.3f rad/s\n",
               tune_names[i], overshoot, decay, omega);
    }

    /* Step 7: Recommend tuning method */
    pid_tune_method_t rec_method = app_recommend_method(id_model);
    const char *rec_name;
    switch (rec_method) {
        case TUNE_METHOD_IMC: rec_name = "IMC/SIMC"; break;
        case TUNE_METHOD_ZN_STEP: rec_name = "Z-N Step Response"; break;
        case TUNE_METHOD_COHEN_COON: rec_name = "Cohen-Coon"; break;
        case TUNE_METHOD_GM_PM: rec_name = "Gain/Phase Margin"; break;
        default: rec_name = "Default"; break;
    }
    printf("\nRecommended tuning method: %s\n", rec_name);

    /* Step 8: Controllability assessment */
    double dtr = fopdt_dead_time_ratio(id_model);
    double ci  = fopdt_controllability_index(id_model);
    printf("Process characteristics:\n");
    printf("  Dead time ratio τ = %.3f ", dtr);
    if (dtr < 0.1) printf("(lag-dominant: easy)\n");
    else if (dtr < 0.3) printf("(balanced)\n");
    else if (dtr < 0.7) printf("(delay-dominant: harder)\n");
    else printf("(very delay-dominant: needs advanced control)\n");
    printf("  Controllability index κ = %.3f\n", ci);

    /* Cleanup */
    free(t); free(u); free(y);

    printf("\nPipeline complete.\n");
    return 0;
}
