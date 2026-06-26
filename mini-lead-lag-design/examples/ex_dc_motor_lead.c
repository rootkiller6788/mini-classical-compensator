#include "lead_design.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== DC Motor Lead Compensator ===\n");
    LeadCompensator c;
    lead_design_direct(25.0, 55.0, 3.0, 7.0, &c);
    printf("C(s)=%.2f*(1+%.3fs)/(1+%.3fs)\n", c.Kc, c.T, c.alpha*c.T);
    printf("Max phase=%.1f deg, DC gain=%.2f\n",
           lead_max_phase_deg(c.alpha), lead_dc_gain(&c));
    double t[30], y[30];
    lead_step_response(&c, 0.3, 30, t, y);
    for (int i = 0; i < 5; i++)
        printf("step[%d]=%.3f\n", i, y[i]);
    printf("...\nstep[29]=%.3f (ss)\n", y[29]);
    printf("Design complete.\n");
    return 0;
}
