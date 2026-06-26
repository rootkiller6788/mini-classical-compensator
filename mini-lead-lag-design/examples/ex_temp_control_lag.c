#include "lag_design.h"
#include <stdio.h>

int main(void) {
    printf("=== Temperature Lag Compensator ===\n");
    LagCompensator c;
    lag_design_direct(0.20, 0.02, 0.5, 10.0, &c);
    printf("C(s)=%.2f*(1+%.2fs)/(1+%.2fs)\n", c.Kc, c.T, c.beta*c.T);
    printf("beta=%.1f, DC gain=%.2f, HF atten=%.1f dB\n",
           c.beta, lag_dc_gain(&c), lag_attenuation_db(c.beta));
    double z, p;
    lag_pole_zero(&c, &z, &p);
    printf("Zero at %.3f, Pole at %.3f\n", z, p);
    double t[30], y[30];
    lag_step_response(&c, 30.0, 30, t, y);
    printf("t=0: y=%.3f, t=30: y=%.3f\n", y[0], y[29]);
    printf("Design complete.\n");
    return 0;
}
