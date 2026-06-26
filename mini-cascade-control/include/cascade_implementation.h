/**
 * @file cascade_implementation.h
 * @brief Digital implementation, anti-windup, and bumpless transfer for cascade control
 *
 * L2 --- Core Concepts: Anti-windup strategies for cascade PID, bumpless transfer.
 * L5 --- Computational Methods: Discrete-time PID, velocity form, digital cascade.
 * L7 --- Applications: Industrial DCS implementation patterns for cascade control.
 *
 * Industrial cascade controllers operate in DCS/PLC environments where:
 * - Inner loop executes at faster rate (e.g., 100 ms)
 * - Outer loop executes at slower rate (e.g., 1 s)
 * - Anti-windup is critical due to actuator saturation
 * - Bumpless transfer required for auto/manual/cascade mode switching
 */

#ifndef CASCADE_IMPLEMENTATION_H
#define CASCADE_IMPLEMENTATION_H

#include "cascade_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L5: Discrete-Time PID Implementation
 * ========================================================================== */

/** Initialize a PID controller with given parameters.
 *  Resets all internal states (integrator, previous error).
 *
 *  @param Kp  Proportional gain
 *  @param Ki  Integral gain (1/s)
 *  @param Kd  Derivative gain (s)
 *  @param N   Derivative filter coefficient
 *  @param Ts  Sampling period (s), 0 = continuous
 *  @param pid [out] Initialized PID controller
 *  @return 0 */
int cascade_pid_init(double Kp, double Ki, double Kd, double N,
                      double Ts, CascadePID *pid);

/** Set PID controller output limits for anti-windup.
 *  @param u_min  Lower saturation limit
 *  @param u_max  Upper saturation limit
 *  @param pid    PID controller to configure */
void cascade_pid_set_limits(double u_min, double u_max, CascadePID *pid);

/** Configure anti-windup method for PID controller.
 *  @param method  Anti-windup method
 *  @param Tt      Tracking time constant (for back-calculation)
 *  @param pid     PID controller */
void cascade_pid_set_antiwindup(CascadeAWMethod method, double Tt,
                                 CascadePID *pid);

/** Execute one step of discrete-time PID controller (velocity form).
 *
 *  Velocity (incremental) form:
 *    du(k) = Kp*(e(k)-e(k-1)) + Ki*Ts*e(k) + (Kd/Ts)*(e(k)-2e(k-1)+e(k-2))
 *    u(k) = u(k-1) + du(k)
 *
 *  Advantages: inherently anti-windup (no integrator accumulation),
 *  bumpless parameter changes, smooth mode transitions.
 *
 *  @param setpoint  Reference r(k)
 *  @param pv        Process variable y(k)
 *  @param pid       PID controller state (updated in-place)
 *  @return Control signal u(k), clamped to [u_min, u_max] */
double cascade_pid_velocity(double setpoint, double pv, CascadePID *pid);

/** Execute one step of discrete-time PID controller (position form).
 *
 *  Position form with back-calculation anti-windup:
 *    P = Kp * (b*setpoint - pv)
 *    I = I_prev + Ki*Ts*e + (Ts/Tt)*(u_sat - u)
 *    D = Kd * (c*setpoint - pv - prev_y) / (Ts + N)  [filtered derivative]
 *    u = P + I + D
 *
 *  @param setpoint  Reference r(k)
 *  @param pv        Process variable y(k)
 *  @param pid       PID controller state (updated in-place)
 *  @return Saturated control signal */
double cascade_pid_position(double setpoint, double pv, CascadePID *pid);

/** Reset PID controller internal states.
 *  Sets integrator to specified value (for bumpless transfer).
 *
 *  @param init_integrator  Initial value for integrator
 *  @param init_u           Initial control signal
 *  @param pid              PID controller to reset */
void cascade_pid_reset(double init_integrator, double init_u, CascadePID *pid);

/** Compute PID frequency response C(jw).
 *
 *  @param pid    PID controller
 *  @param omega  Frequency (rad/s)
 *  @param mag    [out] |C(jw)|
 *  @param phase  [out] arg(C(jw)) in radians
 *  @return 0 */
int cascade_pid_freq_response(const CascadePID *pid, double omega,
                               double *mag, double *phase);

/* ==========================================================================
 * L2: Anti-Windup Implementation
 * ========================================================================== */

/** Initialize anti-windup state for a cascade system.
 *  Configures both inner and outer loop anti-windup.
 *
 *  @param inner_method  Anti-windup method for inner loop
 *  @param outer_method  Anti-windup method for outer loop
 *  @param sys           Cascade system to configure */
void cascade_aw_init(CascadeAWMethod inner_method,
                      CascadeAWMethod outer_method,
                      CascadeSystem *sys);

/** Apply clamping anti-windup: freeze integrator when saturated.
 *  Simpler than back-calculation but can be slower to recover.
 *
 *  @param pid  PID controller with clamping AW state */
void cascade_aw_clamping(CascadePID *pid);

/** Apply back-calculation anti-windup.
 *  Feeds saturation error back through tracking time constant Tt.
 *  di/dt = Ki*e + (1/Tt)*(u_sat - u_unsat)
 *
 *  @param pid      PID controller
 *  @param u_unsat  Unsaturated control signal
 *  @param u_sat    Saturated (applied) control signal */
void cascade_aw_back_calculation(CascadePID *pid, double u_unsat,
                                  double u_sat);

/** Cascade-specific combined anti-windup.
 *  When inner loop saturates: freeze outer integrator too (cascade AW).
 *  This prevents outer loop windup when inner actuator is saturated.
 *
 *  Reference: Astrom & Hagglund, "Advanced PID Control" (2006), Sec. 7.5.
 *
 *  @param sys  Cascade system */
void cascade_aw_combined(CascadeSystem *sys);

/* ==========================================================================
 * L2: Bumpless Transfer Implementation
 * ========================================================================== */

/** Initialize bumpless transfer state machine.
 *
 *  @param initial_mode  0=manual, 1=auto(single-loop), 2=cascade
 *  @param transition_t  Time for bumpless transition (s)
 *  @param bt            [out] Bumpless state */
void cascade_bt_init(int initial_mode, double transition_t,
                      CascadeBumpless *bt);

/** Initiate bumpless transfer to target mode.
 *  Tracks current controller output so target controller starts
 *  from the current operating point.
 *
 *  @param target_mode   Target mode (0=manual, 1=auto, 2=cascade)
 *  @param bt            Bumpless state
 *  @param sys           Cascade system
 *  @param current_time  Current absolute time
 *  @return 0 on success, -1 if already in transition */
int cascade_bt_start(int target_mode, CascadeBumpless *bt,
                      CascadeSystem *sys, double current_time);

/** Update bumpless transfer state (call each sample).
 *  During transition: linearly blends from old output to new output.
 *  After transition_time: switches mode and copies integrator values.
 *
 *  @param bt           Bumpless state
 *  @param sys          Cascade system
 *  @param current_time Current absolute time
 *  @param auto_u       Current auto-mode control signal
 *  @param cascade_u    Current cascade-mode control signal
 *  @return Effective control signal during/after transition */
double cascade_bt_update(CascadeBumpless *bt, CascadeSystem *sys,
                          double current_time,
                          double auto_u, double cascade_u);

/** Check if bumpless transfer is complete.
 *  @param bt  Bumpless state
 *  @return 1 if transfer complete, 0 if in progress */
int cascade_bt_is_complete(const CascadeBumpless *bt);

/* ==========================================================================
 * L5: Cascade Loop Execution
 * ========================================================================== */

/** Execute one sample of the complete cascade control system.
 *
 *  Sequence:
 *    1. Read inner PV (fast measurement)
 *    2. Execute inner PID -> inner control signal u_i
 *    3. Apply inner anti-windup if saturated
 *    4. If new outer sample due:
 *       a. Read outer PV (slow measurement)
 *       b. Execute outer PID -> inner setpoint r_i
 *       c. Apply outer anti-windup
 *    5. Apply bumpless transfer logic
 *
 *  @param sys          Cascade system
 *  @param setpoint     Outer loop setpoint r_o
 *  @param inner_pv     Inner loop process variable y_i
 *  @param outer_pv     Outer loop process variable y_o
 *  @param bt           Bumpless transfer state (NULL if not used)
 *  @param inner_sp     [out] Inner loop setpoint (from outer controller)
 *  @param inner_u      [out] Inner control signal u_i
 *  @return 0 on success, -1 on error */
int cascade_execute(CascadeSystem *sys,
                     double setpoint,
                     double inner_pv,
                     double outer_pv,
                     CascadeBumpless *bt,
                     double *inner_sp,
                     double *inner_u);

/** Toggle cascade mode on/off.
 *  When switching off, outer controller tracks inner PV (bumpless).
 *  When switching on, outer integrator initialized for smooth transition.
 *
 *  @param sys           Cascade system
 *  @param cascade_on    1 to enable cascade, 0 for single-loop
 *  @param bt            Bumpless transfer state */
void cascade_set_mode(CascadeSystem *sys, int cascade_on,
                       CascadeBumpless *bt);

/* ==========================================================================
 * L5: Digital Discretization
 * ========================================================================== */

/** Discretize a continuous-time transfer function using Tustin (bilinear)
 *  approximation: s <- (2/Ts) * (z-1)/(z+1)
 *
 *  Preserves stability: LHP s-plane -> unit circle z-plane.
 *  Frequency warping: w_a = (2/Ts)*tan(w*Ts/2)
 *
 *  @param ct_tf   Continuous-time transfer function G(s)
 *  @param Ts      Sampling period (s)
 *  @param dt_tf   [out] Discrete-time transfer function G(z)
 *  @return 0 on success */
int cascade_tustin_discretize(const CascadeTF *ct_tf, double Ts,
                               CascadeDTF *dt_tf);

/** Discretize a PID controller for digital implementation.
 *  Uses: integrator -> forward Euler, derivative -> backward Euler + filter.
 *
 *  @param ct_pid  Continuous-time PID parameters
 *  @param Ts      Sampling period
 *  @param dt_pid  [out] Discrete-time PID (Kp, Ki, Kd, N same; Ts set) */
void cascade_pid_discretize(const CascadePID *ct_pid, double Ts,
                             CascadePID *dt_pid);

/** Free a discrete-time transfer function.
 *  @param dtf  Discrete TF to free */
void cascade_dtf_free(CascadeDTF *dtf);

/** Execute one step of discrete-time transfer function.
 *  y(k) = b0*u(k) + b1*u(k-1) + ... - a1*y(k-1) - a2*y(k-2) - ...
 *
 *  @param dtf       Discrete transfer function
 *  @param input     Current input u(k)
 *  @param u_history Input history (length = num_order+1, updated in-place)
 *  @param y_history Output history (length = den_order, updated in-place)
 *  @return Current output y(k) */
double cascade_dtf_execute(const CascadeDTF *dtf, double input,
                            double *u_history, double *y_history);

#ifdef __cplusplus
}
#endif

#endif /* CASCADE_IMPLEMENTATION_H */
