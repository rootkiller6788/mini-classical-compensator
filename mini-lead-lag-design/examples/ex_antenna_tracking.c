#include "compensator_freq_design.h"
#include "lead_lag_design.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Antenna Tracking - Frequency-Domain Lead-Lag ===\n");
    printf("Plant: G(s)=100/(s*(s+2)*(s+20))\n");
    double pm_cur = 15.0, pm_tgt = 50.0, wc = 12.0, ess_r = 10.0;
    LeadLagCompensator c;
    freq_design_lead_lag(pm_cur, pm_tgt, wc, ess_r, 5.0, &c);
    printf("Designed: Kc=%.3f, T_lead=%.4f, alpha=%.3f\n",
           c.Kc, c.T_lead, c.alpha);
    printf("T_lag=%.4f, beta=%.1f\n", c.T_lag, c.beta);
    printf("Lead phase boost=%.1f deg\n", lead_max_phase_deg(c.alpha));
    double wm = 1.0/(c.T_lead * sqrt(c.alpha));
    printf("Lead center=%.1f rad/s, Lag corner=%.3f rad/s\n", wm, 1.0/c.T_lag);
    LeadLagFreqPoint fp = lead_lag_evaluate(&c, wc);
    printf("|C(j*wc)|=%.1f dB, arg=%.1f deg\n", fp.magnitude_db, fp.phase_deg);
    printf("Design complete.\n");
    return 0;
}
