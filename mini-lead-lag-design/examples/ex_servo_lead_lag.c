#include "lead_lag_design.h"
#include <stdio.h>

int main(void) {
    printf("=== Servo Lead-Lag Compensator ===\n");
    LeadLagCompensator c;
    lead_lag_design_direct(20.0, 55.0, 8.0, 0.15, 0.05, 5.0, &c);
    printf("Kc=%.2f, Lead:T=%.4f,alpha=%.3f, Lag:T=%.4f,beta=%.1f\n",
           c.Kc, c.T_lead, c.alpha, c.T_lag, c.beta);
    printf("DC gain=%.2f, Est BW=%.2f rad/s\n",
           lead_lag_dc_gain(&c), lead_lag_bandwidth_estimate(&c));
    double zl, pl, zg, pg;
    lead_lag_pole_zero(&c, &zl, &pl, &zg, &pg);
    printf("Lead zero=%.1f pole=%.1f, Lag zero=%.3f pole=%.3f\n", zl, pl, zg, pg);
    printf("Design complete.\n");
    return 0;
}
