/**
 * @file ex_disturbance_rejection.c
 * @brief Example: Disturbance feedforward for continuous stirred-tank reactor (CSTR).
 *
 * Demonstrates L6-L7: chemical process control with measurable disturbance
 * feedforward. The CSTR has an exothermic reaction where feed temperature
 * fluctuations are a major disturbance.
 *
 * L7 Application: Industrial process control (chemical/Boeing/NASA/Dow/DuPont).
 */

#include <stdio.h>
#include <math.h>
#include "feedforward_core.h"
#include "feedforward_design.h"

int main(void)
{
    printf("=== CSTR Temperature Control with Disturbance Feedforward ===\n\n");

    /* CSTR model (simplified first-order + delay):
     *
     * Reactor temperature T [K] responds to:
     *   - Jacket temperature Tj [K] (manipulated variable)
     *   - Feed temperature Tf [K] (measured disturbance)
     *   - Feed flow rate F [m^3/s]
     *   - Reaction heat Q_rxn [W]
     *
     * Linearized model around operating point:
     *   dT/dt = -(F/V + UA/(rho*Cp*V)) * T + (UA/(rho*Cp*V)) * Tj + (F/V) * Tf
     *
     * Parameters (for a 1 m^3 reactor with water-like fluid):
     */
    double V   = 1.0;      /* Reactor volume, m^3 */
    double rho = 1000.0;   /* Density, kg/m^3 */
    double Cp  = 4200.0;   /* Heat capacity, J/(kg*K) */
    double UA  = 5000.0;   /* Heat transfer coefficient * area, W/K */
    double F   = 0.0001;   /* Feed flow rate, m^3/s (0.1 L/s) */
    double rhoCpV = rho * Cp * V;
    double tau = rhoCpV / (F * rho * Cp + UA);
    double K_Tj = UA / (F * rho * Cp + UA);
    double K_Tf = (F * rho * Cp) / (F * rho * Cp + UA);

    printf("Reactor parameters:\n");
    printf("  Volume V     = %.2f m^3\n", V);
    printf("  UA           = %.0f W/K\n", UA);
    printf("  Feed flow F  = %.6f m^3/s\n", F);
    printf("  Time constant tau = %.1f s\n", tau);
    printf("  Gain to jacket T: K_Tj = %.4f\n", K_Tj);
    printf("  Gain to feed T:   K_Tf = %.4f\n\n", K_Tf);

    /* Transfer functions:
     * P(s) = K_Tj / (tau*s + 1)  ? jacket temperature ? reactor temperature
     * G_d(s) = K_Tf / (tau*s + 1) ? feed temperature ? reactor temperature */
    double np[] = {K_Tj};
    double dp[] = {1.0, tau};
    TransferFn plant = tf_create(np, 0, dp, 1, 1.0);

    double nd[] = {K_Tf};
    double dd[] = {1.0, tau};
    TransferFn dist_model = tf_create(nd, 0, dd, 1, 1.0);

    /* Feedback controller: PI control
     * C(s) = Kc * (1 + 1/(Ti*s)) = Kc*(Ti*s + 1)/(Ti*s) */
    double Kc = 2.0;
    double Ti = tau / 2.0;
    double nc_pi[] = {Kc, Kc * Ti};
    double dc_pi[] = {0.0, Ti};
    TransferFn fb_ctrl = tf_create(nc_pi, 1, dc_pi, 1, 1.0);

    printf("Feedback controller: PI, Kc=%.2f, Ti=%.1f s\n", Kc, Ti);

    /* Disturbance feedforward design */
    double ff_gain_static;
    ff_design_static_dist_ff(&plant, &dist_model, &ff_gain_static);
    printf("\nStatic disturbance FF gain: %.4f\n", ff_gain_static);
    printf("  (u_ff = %.4f * (T_feed - T_feed_nom))\n\n", ff_gain_static);

    /* Dynamic disturbance feedforward */
    TransferFn dist_ff;
    ff_design_dynamic_dist_ff(&plant, &dist_model, &dist_ff, 1.0);

    double mag_ff, ph_ff;
    tf_freq_response(&dist_ff, 0.0, &mag_ff, &ph_ff);
    printf("Dynamic FF DC gain: %.4f (should match static: %.4f)\n",
           mag_ff, -ff_gain_static);

    /* Evaluate disturbance rejection at various frequencies */
    printf("\n--- Disturbance Rejection Performance ---\n");
    printf("  Freq (rad/s)  |S(jw)|    With FF    Attenuation\n");
    printf("  ------------  -------    -------    -----------\n");

    double freqs[] = {0.0, 0.01, 0.1, 0.5, 1.0, 2.0};
    for (int i = 0; i < 6; i++) {
        double w = freqs[i];
        double attenuation;
        ff_eval_dist_rejection(&plant, &fb_ctrl, &dist_ff,
                               &dist_model, w, &attenuation);

        /* Sensitivity |S(jw)| */
        TransferFn pc = tf_series(&plant, &fb_ctrl);
        double mag_pc, ph_pc;
        tf_freq_response(&pc, w, &mag_pc, &ph_pc);
        tf_free(&pc);
        double re_pc = mag_pc * cos(ph_pc);
        double im_pc = mag_pc * sin(ph_pc);
        double mag_S = 1.0 / sqrt((1.0+re_pc)*(1.0+re_pc) + im_pc*im_pc);

        printf("  %-12.4f  %-7.4f    %-7.4f    %-7.4f\n",
               w, mag_S, attenuation,
               (mag_S > 1e-10) ? attenuation / mag_S : -1.0);
    }

    /* Time-domain simulation */
    printf("\n--- Time-Domain Simulation (Feed temp disturbance) ---\n");

    double dt = 0.5;
    double t_final = 200.0;
    int n = (int)(t_final / dt) + 1;

    /* Operating point */
    double T_op = 350.0;      /* Reactor temperature setpoint, K */
    double Tf_nom = 300.0;    /* Nominal feed temperature, K */
    double Tj_nom = 330.0;    /* Jacket temperature at steady state */

    double T_fb = T_op, T_ff = T_op;
    double Tj_fb = Tj_nom, Tj_ff = Tj_nom;

    double ise_fb = 0.0, ise_ff = 0.0;

    printf("  Time(s)   T_feed   T_reactor(FB)  T_reactor(FF)  Tj_FF\n");
    printf("  --------  ------   -------------  -------------  -----\n");

    for (int k = 0; k < n; k++) {
        double t = k * dt;

        /* Feed temperature disturbance: sinusoidal + step change */
        double Tf = Tf_nom;
        if (t > 30.0) Tf += 10.0;                     /* Step +10K at t=30 */
        if (t > 100.0) Tf += 5.0 * sin(0.05 * t);     /* Sinusoidal after t=100 */

        /* --- Feedback only --- */
        double err_fb = T_op - T_fb;
        /* PI controller: u = Kc*(err + 1/Ti * integral(err)) */
        double Tj_fb_new = Tj_nom + Kc * err_fb;
        /* Clamp */
        if (Tj_fb_new > 400.0) Tj_fb_new = 400.0;
        if (Tj_fb_new < 280.0) Tj_fb_new = 280.0;
        Tj_fb = Tj_fb_new;

        /* Reactor dynamics */
        double dT_fb = (-T_fb + K_Tj * Tj_fb + K_Tf * Tf) / tau;
        T_fb += dT_fb * dt;
        ise_fb += (T_op - T_fb) * (T_op - T_fb) * dt;

        /* --- With feedforward --- */
        double err_ff = T_op - T_ff;
        /* FF: compensate feed temperature disturbance */
        double T_ff_correction = ff_gain_static * (Tf - Tf_nom);
        double Tj_ff_new = Tj_nom + Kc * err_ff + T_ff_correction;
        if (Tj_ff_new > 400.0) Tj_ff_new = 400.0;
        if (Tj_ff_new < 280.0) Tj_ff_new = 280.0;
        Tj_ff = Tj_ff_new;

        double dT_ff = (-T_ff + K_Tj * Tj_ff + K_Tf * Tf) / tau;
        T_ff += dT_ff * dt;
        ise_ff += (T_op - T_ff) * (T_op - T_ff) * dt;

        /* Print every 10 seconds */
        if (k % 20 == 0 || k == n - 1) {
            printf("  %-8.1f  %-6.1f  %-13.2f  %-13.2f  %-5.1f\n",
                   t, Tf, T_fb, T_ff, Tj_ff);
        }
    }

    printf("\n--- Summary ---\n");
    printf("  Feedback only ISE:  %.2f K^2*s\n", ise_fb);
    printf("  FB + Feedforward ISE: %.2f K^2*s\n", ise_ff);
    if (ise_fb > 1e-10) {
        printf("  Improvement:        %.1f%%\n",
               100.0 * (1.0 - ise_ff / ise_fb));
    }
    printf("  Static FF gain:      %.4f\n", ff_gain_static);

    tf_free(&plant);
    tf_free(&dist_model);
    tf_free(&fb_ctrl);
    tf_free(&dist_ff);

    printf("\n=== Done ===\n");
    return 0;
}
