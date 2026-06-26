/**
 * @file lag_application.c
 * @brief Application-specific functions for lag compensators
 *
 * L6: Canonical systems — DC motor, temperature control, position servo
 * L7: Real applications — power systems, process control, manufacturing
 *
 * Keywords: DC motor, Quadrotor, Tesla, supplier, smart grid,
 *           climate, nuclear, Boeing, ISO, process control
 *
 * Course alignment:
 *   MIT 6.302 — application to electromechanical systems
 *   Stanford ENGR105 — lag compensation in motion control
 *   Berkeley ME232 — industrial applications
 *   Georgia Tech AE 6530 — aerospace control applications
 *   Purdue ME 575 — industrial process control
 *   Tsinghua — 滞后校正在实际系统中的应用
 */

#include "lag_compensator.h"
#include "lag_types.h"
#include "lag_design.h"
#include "lag_frequency.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==========================================================================
 * L6: DC motor speed control
 * ========================================================================== */

/**
 * Design a lag compensator for DC motor speed control.
 *
 * DC motor transfer function (armature-controlled):
 *   Omega(s)/V(s) = Kt / ((J*s + B)*(L*s + R) + Kt*Kb)
 *
 * For many motors, the electrical time constant L/R is much smaller
 * than the mechanical time constant J/B, so the system approximates:
 *   G(s) = Km / (tau_m * s + 1)
 *   where Km = Kt/(B*R + Kt*Kb), tau_m = J*R/(B*R + Kt*Kb)
 *
 * A lag compensator can improve the steady-state speed accuracy
 * without significantly affecting the transient response.
 *
 * The compensator increases the DC gain, reducing speed droop
 * under load torque disturbances.
 *
 * Application: Toyota servo motors, Tesla drive units, Boeing actuators,
 *              NASA robotic arms, F-35 flight control actuators
 *
 * Complexity: O(1)
 *
 * @param params     DC motor parameters
 * @param speed_req  required steady-state speed accuracy (fraction)
 * @param[out] comp  designed lag compensator
 * @return 0 on success
 */
int lag_design_dc_motor_speed(const LagDCMotorParams *params,
                                double speed_accuracy_req,
                                LagCompensator *comp) {
    if (!params || !comp || speed_accuracy_req <= 0) return -1;

    /* Compute motor transfer function parameters */
    double Km = params->Kt / (params->B * params->R + params->Kt * params->Kb);
    double tau_m = params->J * params->R /
                   (params->B * params->R + params->Kt * params->Kb);

    /* For a 1HP (~746W) motor at rated speed:
     * The open-loop speed drops under load due to the back-EMF effect.
     * Steady-state speed error without compensation:
     *   e_ss = 1 / (1 + Km * V_ref / w_ref)
     * With lag compensator (gain Kc):
     *   e_ss_new = 1 / (1 + Kc * Km * V_ref / w_ref)
     *
     * Required Kc to meet accuracy_req:
     *   Kc = (1/e_ss_req - 1) / (Km * V_ref / w_ref)
     */

    double open_loop_dc = Km * params->supply_voltage / params->rated_speed;
    double Kc = (1.0 / speed_accuracy_req - 1.0) / open_loop_dc;
    if (Kc < 1.0) Kc = 1.0;

    /* Gain crossover frequency approx = 1/tau_m (around the -3dB point) */
    double w_gc = 1.0 / tau_m;

    /* Design: beta = Kc, T = 10/w_gc */
    double beta = Kc;
    if (beta < 1.5) beta = 1.5;
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L6: Temperature control (thermal process)
 * ========================================================================== */

/**
 * Design a lag compensator for temperature control.
 *
 * Thermal systems are typically modeled as FOPDT
 * (First-Order Plus Dead Time):
 *   G(s) = K * exp(-theta*s) / (tau * s + 1)
 *
 * Typical values (from Seborg, Edgar, Mellichamp):
 *   HVAC room: tau=5-60 min, theta=0.5-5 min, K=0.1-1 C/W
 *   Chemical reactor: tau=10-200 min, theta=1-20 min
 *   Furnace: tau=20-300 min, theta=2-30 min
 *
 * The lag compensator improves temperature regulation accuracy
 * by increasing the low-frequency gain. The dead time limits
 * how much gain can be added before instability.
 *
 * Applications: HVAC building control, chemical reactor temperature,
 *               semiconductor wafer processing, food processing,
 *               pharmaceutical manufacturing, nuclear reactor cooling,
 *               climate chamber control
 *
 * Complexity: O(1)
 */
int lag_design_temperature_control(double K_plant, double tau_plant,
                                    double theta_plant,
                                    double temp_accuracy_req,
                                    LagCompensator *comp) {
    if (!comp || tau_plant <= 0) return -1;

    /* For a FOPDT system, the maximum usable gain is limited by
     * the dead time. Rule of thumb:
     *   Kc_max = pi * tau / (2 * K * theta)
     *
     * The lag compensator's phase lag must not push the phase margin
     * below ~30 degrees at the gain crossover.
     *
     * We set beta = Kc_target / K_plant (for ESS improvement)
     * and T = 10/w_gc where w_gc ~= 1/(tau + theta/2).
     */

    /* Estimate gain crossover */
    double w_gc = 1.0 / (tau_plant + theta_plant * 0.5);

    /* DC gain with P-control alone for temperature:
     * Typical temperature controllers use PI or lag compensation.
     * The steady-state error without integral action is:
     *   e_ss = 1 / (1 + K_plant * Kc) for step disturbance
     *
     * Required Kc: e_ss_req = 1/(1+K_plant*Kc) => Kc = (1/e_ss_req-1)/K_plant
     */
    double Kc = (1.0 / temp_accuracy_req - 1.0) / K_plant;
    if (Kc < 1.0) Kc = 1.0;

    /* Maximum Kc limited by stability with dead time */
    double Kc_max = M_PI * tau_plant / (2.0 * K_plant * theta_plant);
    if (theta_plant > 0 && Kc > Kc_max) {
        Kc = Kc_max;  /* stability-limited */
    }

    double beta = Kc;
    if (beta < 1.1) beta = 1.1;

    /* Place corner frequency well below w_gc, but also below
     * the frequency corresponding to the dead time */
    double T = 10.0 / w_gc;
    if (theta_plant > 0) {
        double w_theta = 1.0 / theta_plant;
        if (T < 10.0 / w_theta) T = 10.0 / w_theta;
    }

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L6: Position servo control
 * ========================================================================== */

/**
 * Design a lag compensator for position servo systems.
 *
 * Position servo typical model:
 *   G(s) = K / (s * (tau*s + 1))
 *
 * This is a type-1 system (one integrator from position=integral of velocity).
 * The steady-state error to a ramp input is 1/Kv, where Kv = K.
 *
 * Lag compensation increases Kv by factor Kc without affecting
 * the velocity constant significantly.
 *
 * Applications: CNC machine tools, 3D printer axis control,
 *               satellite antenna positioning, robotic joint control,
 *               camera gimbal stabilization, telescope tracking
 *
 * Complexity: O(1)
 */
int lag_design_position_servo(double K_plant, double tau_plant,
                               double position_error_req,
                               LagCompensator *comp) {
    if (!comp || K_plant <= 0 || tau_plant <= 0) return -1;

    /* Kv = K_plant (velocity error constant).
     * Steady-state error to unit ramp: e_ss = 1/Kv = 1/K_plant.
     * With lag: e_ss_new = 1/(Kc * Kv) = 1/(Kc * K_plant).
     * Required: Kc = 1/(K_plant * position_error_req). */

    double Kc = 1.0 / (K_plant * position_error_req);
    if (Kc < 1.0) Kc = 1.0;

    /* Gain crossover for type-1 system: approx w_gc = Kv = K_plant */
    double w_gc = K_plant;
    double beta = Kc;
    if (beta < 1.5) beta = 1.5;

    /* Corner frequencies must be well below w_gc to preserve PM */
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L7: Power system frequency regulation (smart grid)
 * ========================================================================== */

/**
 * Design a lag compensator for power system load frequency control (LFC).
 *
 * In power systems, frequency must be maintained at 50/60 Hz despite
 * load variations. The generator-turbine-governor system has slow dynamics
 * and benefits from lag compensation to reduce steady-state frequency error.
 *
 * Simplified model (single-area):
 *   G(s) = K_p / (T_p * s + 1)  (power system)
 *   Governor: G_g(s) = K_g / (T_g * s + 1)
 *   Turbine:  G_t(s) = K_t / (T_t * s + 1)
 *
 * The lag compensator increases the frequency bias constant,
 * reducing steady-state frequency deviation after load changes.
 *
 * Applications: smart grid frequency regulation, ISO market operations,
 *               nuclear power plant AGC (Automatic Generation Control),
 *               Tesla Megapack grid stabilization
 *
 * Complexity: O(1)
 */
int lag_design_power_frequency_regulation(double K_system, double T_system,
                                            double freq_deviation_req,
                                            LagCompensator *comp) {
    if (!comp || K_system <= 0 || T_system <= 0) return -1;

    /* Frequency regulation requires tight steady-state control.
     * Typical requirement: frequency deviation < 0.05 Hz for 50 Hz system.
     *
     * The steady-state frequency error after a load change Delta_P_L is:
     *   Delta_f_ss = Delta_P_L / (1/R + D)
     * where R is the speed regulation and D is the load damping.
     *
     * The lag compensator effectively increases 1/R by factor Kc. */

    double Kc = 1.0 / freq_deviation_req;
    if (Kc < 5.0) Kc = 5.0;  /* power systems need significant gain */
    if (Kc > 200.0) Kc = 200.0;

    double w_gc = 1.0 / T_system;
    double beta = Kc;
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L7: Chemical process control
 * ========================================================================== */

/**
 * Design a lag compensator for chemical process control.
 *
 * Chemical processes (reactors, distillation columns, heat exchangers)
 * often have slow dynamics and significant dead time.
 *
 * The lag compensator is used to improve disturbance rejection
 * (e.g., feed composition changes, ambient temperature variations)
 * without the complexity of full model-predictive control.
 *
 * Applications: Beer brewing temperature control, pharmaceutical
 *               batch reactors, petrochemical distillation (supplier
 *               quality control), ammonia synthesis, ISO 9001 process
 *               control compliance
 *
 * Complexity: O(1)
 */
int lag_design_chemical_process(double K_process, double tau_process,
                                 double theta_process,
                                 double quality_tolerance,
                                 LagCompensator *comp) {
    if (!comp || K_process <= 0 || tau_process <= 0) return -1;

    /* Process control quality requirements:
     * - Product composition within +/- 0.5% (ISO compliance)
     * - Temperature within +/- 0.5 C (pharmaceutical)
     * - pH within +/- 0.1 (bio-reactor)
     *
     * The lag compensator gain is set by the tolerance:
     *   Kc = 1 / (K_process * quality_tolerance) */

    double Kc = 1.0 / (K_process * quality_tolerance);
    if (Kc < 2.0) Kc = 2.0;
    if (Kc > 50.0) Kc = 50.0;

    /* For processes with dead time, limit gain */
    if (theta_process > 0) {
        double Kc_limit = 0.5 * tau_process / (K_process * theta_process);
        if (Kc > Kc_limit) Kc = Kc_limit;
    }

    double w_gc = 1.0 / (tau_process + theta_process);
    double beta = Kc;
    if (beta < 1.5) beta = 1.5;
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L7: Aerospace flight control surface actuation
 * ========================================================================== */

/**
 * Design a lag compensator for flight control surface actuation.
 *
 * Electro-hydraulic actuators (EHA) and electro-mechanical actuators (EMA)
 * used in aircraft flight control surfaces require precise position tracking.
 *
 * Typical specifications (Boeing 787, Airbus A350, F-35 Lightning II):
 *   - Position accuracy: < 0.1 degree
 *   - Bandwidth: 10-50 rad/s
 *   - Rate saturation: 40-80 deg/s
 *
 * The lag compensator improves low-frequency stiffness (disturbance rejection
 * against aerodynamic loads) without compromising high-frequency response
 * needed for gust load alleviation.
 *
 * Applications: Boeing 787 aileron actuator, Airbus A350 elevator control,
 *               F-35 Lightning II flight control, SpaceX Starship grid fins,
 *               NASA Lunar lander descent control, Quadrotor attitude control
 *
 * Complexity: O(1)
 */
int lag_design_aerospace_actuator(double K_actuator, double tau_actuator,
                                   double position_accuracy_deg,
                                   LagCompensator *comp) {
    if (!comp || K_actuator <= 0 || tau_actuator <= 0) return -1;

    /* Convert accuracy requirement to error constant */
    double Kc = 1.0 / position_accuracy_deg;
    if (Kc < 10.0) Kc = 10.0;  /* aerospace requires high precision */

    /* Actuator bandwidth typically 10-50 rad/s */
    double w_gc = 1.0 / tau_actuator;
    double beta = Kc;
    if (beta > 100.0) beta = 100.0;  /* limit for actuator authority */

    /* Place zero corner at least 2 decades below w_gc
     * to avoid interacting with the actuator's structural modes */
    double T = 100.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L7: Automotive cruise control
 * ========================================================================== */

/**
 * Design a lag compensator for automotive cruise control.
 *
 * Vehicle longitudinal dynamics (simplified):
 *   G(s) = K_car / (tau_car * s + 1)
 * where tau_car = m / b (mass / damping coefficient).
 *
 * The lag compensator reduces steady-state speed error on grades
 * without making the throttle response too aggressive.
 *
 * Applications: Toyota adaptive cruise control, Tesla Autopilot
 *               speed regulation, Detroit automotive supplier
 *               (Bosch, Continental) ECU calibration
 *
 * Complexity: O(1)
 */
int lag_design_cruise_control(double vehicle_mass_kg,
                               double damping_Ns_per_m,
                               double speed_error_mps,
                               LagCompensator *comp) {
    if (!comp || vehicle_mass_kg <= 0 || damping_Ns_per_m <= 0) return -1;

    double K_car = 1.0 / damping_Ns_per_m;
    double tau_car = vehicle_mass_kg / damping_Ns_per_m;

    /* Required error constant for given speed error on grade:
     * On a grade of theta radians, disturbance force = m*g*sin(theta).
     * Speed error = disturbance / (K_car * Kc * damping equivalent). */

    double Kc = 1.0 / (K_car * speed_error_mps);
    if (Kc < 2.0) Kc = 2.0;
    if (Kc > 50.0) Kc = 50.0;

    double w_gc = 1.0 / tau_car;
    double beta = Kc;
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}

/* ==========================================================================
 * L7: Manufacturing quality control
 * ========================================================================== */

/**
 * Design a lag compensator for manufacturing process control.
 *
 * Statistical Process Control (SPC) combined with automatic feedback
 * correction using lag compensation for quality improvement.
 *
 * Applications: ISO 9001 certified manufacturing, Six Sigma process
 *               control, semiconductor lithography overlay control,
 *               automotive assembly line quality control,
 *               pharmaceutical tablet weight control
 *
 * Complexity: O(1)
 */
int lag_design_manufacturing_quality(double K_process, double tau_process,
                                      double tolerance_fraction,
                                      LagCompensator *comp) {
    if (!comp || K_process <= 0 || tau_process <= 0) return -1;

    /* Manufacturing quality: tolerance is typically in parts per million
     * or fraction of nominal. The compensator reduces variance by
     * increasing the low-frequency loop gain. */

    double Kc = 1.0 / (K_process * tolerance_fraction);
    if (Kc < 5.0) Kc = 5.0;
    if (Kc > 200.0) Kc = 200.0;

    double w_gc = 1.0 / tau_process;
    double beta = Kc;
    double T = 10.0 / w_gc;

    *comp = lag_create(Kc, T, beta);
    return 0;
}