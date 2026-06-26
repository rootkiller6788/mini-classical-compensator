/**
 * @file ex_2dof_temperature.c
 * @brief Example: 2-DOF temperature control with prefilter + disturbance FF.
 *
 * Demonstrates L6-L7: combined reference prefiltering and disturbance
 * feedforward for a building HVAC zone temperature control system.
 * L7 Application: smart grid / building automation / ISO 50001 energy management.
 */

#include <stdio.h>
#include <math.h>
#include "feedforward_core.h"
#include "feedforward_design.h"
#include "feedforward_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== 2-DOF HVAC Zone Temperature Control ===\n\n");

    /* Room thermal model (ISO 13790 simplified):
     * C * dT/dt = (T_out - T) / R_wall + Q_heat + Q_solar + Q_occupants
     *
     * Parameters for a medium office zone (20 m^2):
     */
    double C_thermal = 500000.0;   /* Thermal capacitance, J/K */
    double R_wall    = 0.005;      /* Wall thermal resistance, K/W */
    double tau_room = C_thermal * R_wall;  /* Time constant ~2500s ~ 42min */

    printf("Room thermal model:\n");
    printf("  Thermal capacitance C = %.0f J/K\n", C_thermal);
    printf("  Wall resistance R     = %.4f K/W (K_heat=%.4f)\n", R_wall, R_wall);
    printf("  Time constant tau     = %.0f s (%.1f min)\n\n", tau_room, tau_room/60.0);

    /* Transfer functions:
     * P(s) = K_heat / (tau_room*s + 1)  ? heating power ? room temperature
     * G_d(s) = K_dist / (tau_room*s + 1) ? outdoor temp ? room temperature
     * where K_dist = R_wall / R_wall = 1.0 (simplified) */
    double K_heat_val = 0.005;
    double np[] = {K_heat_val};
    double dp[] = {1.0, tau_room};
    TransferFn plant = tf_create(np, 0, dp, 1, 1.0);

    double nd[] = {0.8};  /* 80% of outdoor temperature affects indoor */
    double dd[] = {1.0, tau_room};
    TransferFn dist_model = tf_create(nd, 0, dd, 1, 1.0);

    /* Feedback: PI controller
     * C(s) = Kc * (1 + 1/(Ti*s)) */
    double Kc = 500.0;    /* Proportional gain */
    double Ti = tau_room / 3.0;  /* Integral time */
    double nc_pi[] = {Kc, Kc * Ti};
    double dc_pi[] = {0.0, Ti};
    TransferFn fb_ctrl = tf_create(nc_pi, 1, dc_pi, 1, 1.0);

    /* Reference prefilter: second-order, tr=600s (10 min) */
    Prefilter prefilter;
    ff_design_prefilter_2nd(600.0, 0.7, &prefilter);
    printf("Prefilter: 2nd order, tr=600s, zeta=0.7, bw=%.4f rad/s\n",
           prefilter.bandwidth);

    /* Disturbance feedforward */
    double ff_gain;
    ff_design_static_dist_ff(&plant, &dist_model, &ff_gain);
    printf("Disturbance FF gain: %.2f W/K (outdoor temperature compensation)\n\n",
           ff_gain);

    /* Simulation: 24-hour period with occupancy schedule */
    double dt = 60.0;        /* 1 minute time step */
    double t_final = 86400.0; /* 24 hours */
    int n = (int)(t_final / dt) + 1;

    double T_occupied = 21.0;  /* Occupied setpoint, C */
    double T_room_fb = 20.0;   /* Initial room temp, C */
    double T_room_ff = 20.0;

    /* Prefilter state for smooth setpoint transitions */
    double T_sp_filtered = T_room_fb;
    double alpha = dt / (dt + 600.0);  /* Smoothing ~10min */

    double heat_fb = 0.0, heat_ff = 0.0;
    double integral_fb = 0.0, integral_ff = 0.0;

    double energy_fb = 0.0, energy_ff = 0.0;
    double comfort_violation_fb = 0.0, comfort_violation_ff = 0.0;

    printf("--- 24-hour Simulation ---\n");
    printf("  Hour  T_out  T_set  T_room(FB)  T_room(FF)  Heat(FB)  Heat(FF)\n");
    printf("  ----  -----  -----  ----------  ----------  --------  --------\n");

    for (int k = 0; k < n; k++) {
        double t = k * dt;
        double hour = t / 3600.0;

        /* Outdoor temperature: sinusoidal diurnal cycle */
        double T_out = 10.0 + 8.0 * sin(2.0 * M_PI * (hour - 8.0) / 24.0);

        /* Occupancy schedule: setback at night */
        double T_set;
        if (hour >= 7.0 && hour < 19.0) {
            T_set = T_occupied;  /* Occupied: 7AM-7PM */
        } else {
            T_set = 17.0;  /* Night setback */
        }

        /* Prefilter the setpoint change */
        T_sp_filtered = alpha * T_set + (1.0 - alpha) * T_sp_filtered;

        /* Solar gain (simplified: peak at noon) */
        double Q_solar = 0.0;
        if (hour >= 6.0 && hour < 18.0) {
            Q_solar = 500.0 * sin(M_PI * (hour - 6.0) / 12.0);
        }

        /* Occupant heat gain */
        double Q_occ = (hour >= 7.0 && hour < 19.0) ? 300.0 : 50.0;

        /* Total disturbance: outdoor temp effect + solar + occupants */
        double T_out_effect = 0.8 * (T_out - 20.0);  /* Deviation from 20C */
        double Q_dist = (T_out_effect / R_wall) + Q_solar + Q_occ;

        /* --- Feedback only --- */
        double err_fb = T_sp_filtered - T_room_fb;
        integral_fb += err_fb * dt;
        heat_fb = Kc * (err_fb + integral_fb / Ti);
        if (heat_fb < 0.0) heat_fb = 0.0;  /* Heating only */
        if (heat_fb > 5000.0) heat_fb = 5000.0;

        double dT_fb = (Q_dist + heat_fb - T_room_fb / R_wall) / C_thermal;
        T_room_fb += dT_fb * dt;
        energy_fb += heat_fb * dt / 3600000.0;  /* kWh */
        if (T_room_fb < T_sp_filtered - 1.0) comfort_violation_fb += dt;

        /* --- With feedforward --- */
        double err_ff = T_sp_filtered - T_room_ff;
        integral_ff += err_ff * dt;
        /* Feedforward: compensate for measured outdoor temp */
        double ff_correction = ff_gain * (T_out - 20.0);
        /* Feedforward for solar: simplified estimate based on time */
        double solar_ff = (hour >= 6.0 && hour < 18.0) ?
                          -100.0 * sin(M_PI * (hour - 6.0) / 12.0) : 0.0;
        heat_ff = Kc * (err_ff + integral_ff / Ti) + ff_correction + solar_ff;
        if (heat_ff < 0.0) heat_ff = 0.0;
        if (heat_ff > 5000.0) heat_ff = 5000.0;

        double dT_ff = (Q_dist + heat_ff - T_room_ff / R_wall) / C_thermal;
        T_room_ff += dT_ff * dt;
        energy_ff += heat_ff * dt / 3600000.0;
        if (T_room_ff < T_sp_filtered - 1.0) comfort_violation_ff += dt;

        /* Print every hour */
        if (k % 60 == 0 || k == n - 1) {
            printf("  %4.0f  %5.1f  %5.1f  %10.2f  %10.2f  %8.0f  %8.0f\n",
                   hour, T_out, T_sp_filtered, T_room_fb, T_room_ff,
                   heat_fb, heat_ff);
        }
    }

    printf("\n--- 24-hour Performance Summary ---\n");
    printf("  Energy consumption:\n");
    printf("    Feedback only:       %.2f kWh\n", energy_fb);
    printf("    FB + Feedforward:    %.2f kWh\n", energy_ff);
    printf("    Energy savings:      %.1f%%\n",
           100.0 * (1.0 - energy_ff / (energy_fb + 1e-15)));
    printf("  Comfort violations (hours below setpoint-1):\n");
    printf("    Feedback only:       %.1f h\n", comfort_violation_fb / 3600.0);
    printf("    FB + Feedforward:    %.1f h\n", comfort_violation_ff / 3600.0);
    printf("\n  (ISO 50001 / Smart Grid energy management context)\n");

    tf_free(&plant);
    tf_free(&dist_model);
    tf_free(&fb_ctrl);
    tf_free(&prefilter.filter);

    printf("\n=== Done ===\n");
    return 0;
}
