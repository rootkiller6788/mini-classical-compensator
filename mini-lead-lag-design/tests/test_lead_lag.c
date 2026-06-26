#include "lead_design.h"
#include "lag_design.h"
#include "lead_lag_design.h"
#include "compensator_spec.h"
#include "compensator_freq_design.h"
#include "compensator_rlocus_design.h"
#include "compensator_digital.h"
#include "compensator_sensitivity.h"
#include "compensator_analytical.h"
#include "compensator_loopshaping.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-6
static int run = 0, pass = 0;

int main(void) {
    printf("=== mini-lead-lag-design Tests ===\n");

    /* Lead compensator tests */
    run++; if (fabs(lead_max_phase_deg(0.1) - 54.9032) < 0.01) pass++;
    else printf("FAIL: lead_max_phase_deg\n");
    run++; if (lead_alpha_from_phase(55.0) > 0.05 && lead_alpha_from_phase(55.0) < 1.0) pass++;
    else printf("FAIL: lead_alpha_from_phase\n");

    LeadCompensator lc;
    run++; if (lead_design_direct(30.0, 55.0, 10.0, 5.0, &lc) == 0 && lead_is_valid(&lc)) pass++;
    else printf("FAIL: lead_design_direct\n");

    LeadFreqPoint fp = lead_evaluate(&lc, 0.0);
    run++; if (fabs(fp.real - lc.Kc) < EPS) pass++;
    else printf("FAIL: lead_evaluate_dc\n");

    double z, p;
    lead_pole_zero(&lc, &z, &p);
    run++; if (z < 0.0 && p < z) pass++;
    else printf("FAIL: lead_pole_zero\n");

    double t[50], y[50];
    lead_step_response(&lc, 1.0, 50, t, y);
    run++; if (y[49] > 0.0) pass++;
    else printf("FAIL: lead_step_response\n");

    double b[2], a1;
    lead_discretize_tustin(&lc, 0.01, b, &a1);
    run++; if (fabs(a1) < 1.0) pass++;
    else printf("FAIL: lead_tustin_stable\n");

    /* Lag compensator tests */
    run++; if (fabs(lag_attenuation_db(10.0) + 20.0) < 0.01) pass++;
    else printf("FAIL: lag_attenuation\n");

    LagCompensator gc;
    run++; if (lag_design_direct(0.2, 0.02, 5.0, 10.0, &gc) == 0 && lag_is_valid(&gc)) pass++;
    else printf("FAIL: lag_design\n");

    run++; if (gc.beta >= 1.0) pass++;
    else printf("FAIL: lag_beta\n");

    /* Lead-lag tests */
    LeadLagCompensator llc;
    run++; if (lead_lag_design_direct(30.0, 55.0, 10.0, 0.2, 0.02, 5.0, &llc) == 0 && lead_lag_is_valid(&llc)) pass++;
    else printf("FAIL: lead_lag_design\n");

    LeadLagFreqPoint llfp = lead_lag_evaluate(&llc, 0.0);
    run++; if (fabs(llfp.real - llc.Kc) < EPS) pass++;
    else printf("FAIL: lead_lag_dc_gain\n");

    /* Spec conversion tests */
    run++; if (fabs(po_to_damping(16.3) - 0.5) < 0.02) pass++;
    else printf("FAIL: po_to_damping\n");

    double zeta = pm_to_damping(45.0);
    run++; if (zeta > 0.4 && zeta < 0.5) pass++;
    else printf("FAIL: pm_to_damping\n");

    /* Bode tests */
    BodeData *bd = bode_data_create(50);
    run++; if (bd != NULL && bd->num_points == 50) pass++;
    else printf("FAIL: bode_data_create\n");
    if (bd) {
        FirstOrderPlant fop = {10.0, 0.1};
        bode_first_order(&fop, 0.1, 100.0, bd);
        run++; if (fabs(bd->mag_db[0] - 20.0) < 1.0) pass++;
        else printf("FAIL: bode_first_order_dc\n");
        bode_data_free(bd);
    }

    /* Root locus tests */
    OpenLoopPZ ol;
    ol.num_poles = 1; ol.poles[0] = -1.0;
    ol.num_zeros = 0;
    double ang = rl_angle_condition(&ol, -2.0, 0.0);
    run++; if (fabs(fmod(fabs(ang), 360.0) - 180.0) < 5.0) pass++;
    else printf("FAIL: rl_angle_condition\n");

    /* Digital filter tests */
    DigitalFilter df;
    run++; if (digital_lead_tustin(&lc, 0.01, &df) == 0 && digital_is_stable(&df)) pass++;
    else printf("FAIL: digital_lead_tustin\n");

    double out = digital_apply(&df, 1.0);
    run++; if (!isnan(out)) pass++;
    else printf("FAIL: digital_apply\n");

    /* Sensitivity tests */
    SensitivityData *sd = sens_data_create(30);
    run++; if (sd != NULL) pass++;
    else printf("FAIL: sens_data_create\n");
    if (sd) sens_data_free(sd);

    /* Analytical design tests */
    FOPDTModel plant = {5.0, 2.0, 0.5};
    double Kp, Ki;
    run++; if (imc_design_pi(&plant, 1.0, &Kp, &Ki) == 0 && Kp > 0.0) pass++;
    else printf("FAIL: imc_design_pi\n");

    run++; if (simc_pi(&plant, &Kp, &Ki) == 0 && Kp > 0.0) pass++;
    else printf("FAIL: simc_pi\n");

    /* Loop shaping tests */
    LoopShape *ls = loopshape_create(30);
    run++; if (ls != NULL && ls->n == 30) pass++;
    else printf("FAIL: loopshape_create\n");
    if (ls) loopshape_free(ls);

    printf("\n=== %d/%d tests passed ===\n", pass, run);
    return (pass == run) ? 0 : 1;
}
