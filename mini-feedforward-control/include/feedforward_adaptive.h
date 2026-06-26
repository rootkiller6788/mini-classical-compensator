/**
 * @file feedforward_adaptive.h
 * @brief Adaptive and learning feedforward control.
 *
 * L5 --- Computational methods: LMS adaptive feedforward,
 * iterative learning control (ILC).
 *
 * L8 --- Advanced topics: adaptive feedforward, nonlinear
 * feedforward compensation, iterative learning control.
 *
 * Course alignment:
 *   MIT 6.243 Nonlinear Control ? adaptive feedforward
 *   Caltech CDS 243 Adaptive Control
 *   Berkeley ME234 Nonlinear ? feedback linearization
 *   Stanford ENGR209A ? nonlinear feedforward
 *   ETH 151-0567 ? learning control
 */

#ifndef FEEDFORWARD_ADAPTIVE_H
#define FEEDFORWARD_ADAPTIVE_H

#include "feedforward_core.h"

/* ==========================================================================
 * L1: Adaptive Feedforward Structures
 * ========================================================================== */

/** LMS (Least Mean Squares) adaptive filter for feedforward.
 *
 *  Structure: adaptive FIR filter that learns the feedforward
 *  compensator online to minimize a performance criterion.
 *
 *  Application: Active Noise Control (ANC), vibration cancellation,
 *  repetitive disturbance rejection.
 */
typedef struct {
    double  *weights;       /**< adaptive FIR filter coefficients */
    int      n_taps;        /**< number of filter taps */
    double   mu;            /**< step size / learning rate */
    double  *buffer;        /**< input buffer (circular) */
    int      buf_idx;       /**< current buffer index */
    double   leakage;       /**< leakage factor (0..1, prevents drift) */
    double   norm_factor;   /**< normalization factor */
} LMSFilter;

/** Iterative Learning Control (ILC) structure.
 *
 *  ILC improves feedforward from trial to trial:
 *    u_{k+1}(t) = Q(s) * [u_k(t) + L(s) * e_k(t)]
 *
 *  where e_k(t) = r(t) - y_k(t) is the tracking error on trial k,
 *  L(s) is the learning filter, Q(s) is the Q-filter (robustness).
 */
typedef struct {
    TransferFn  learning_filter;   /**< L(s) learning gain */
    TransferFn  q_filter;          /**< Q(s) robustness filter */
    double     *trial_data;        /**< stored trial signal data */
    double     *error_data;        /**< error from previous trial */
    int         n_samples;         /**< number of time samples */
    double      dt;                /**< sample time (s) */
    int         trial_num;         /**< current trial number */
    double      convergence_rate;  /**< estimated convergence rate */
} ILCController;

/** Nonlinear feedforward via computed torque / feedback linearization.
 *
 *  For a nonlinear system: M(q)*q_ddot + C(q,q_dot)*q_dot + G(q) = tau
 *  Feedforward: tau_ff = M(q_d)*q_ddot_d + C(q_d,q_dot_d)*q_dot_d + G(q_d)
 *
 *  This cancels the nonlinear dynamics along the desired trajectory.
 */
typedef struct {
    /** Inertia matrix function pointer: M(q) */
    void (*inertia)(const double *q, double *M_out, int n_dof, void *params);
    /** Coriolis/centripetal: C(q,qd)*qd */
    void (*coriolis)(const double *q, const double *qd, double *C_out,
                     int n_dof, void *params);
    /** Gravity vector: G(q) */
    void (*gravity)(const double *q, double *G_out, int n_dof, void *params);
    /** Friction model */
    void (*friction)(const double *qd, double *F_out, int n_dof, void *params);
    int    n_dof;          /**< number of degrees of freedom */
    void  *params;         /**< user-defined model parameters */
} NonlinearFFModel;

/* ==========================================================================
 * L5: LMS Adaptive Feedforward
 * ========================================================================== */

/**
 * Initialize an LMS adaptive filter.
 *
 * @param filter   Filter to initialize
 * @param n_taps   Number of FIR taps
 * @param mu       Learning rate (0 < mu < 2/lambda_max)
 * @param leakage  Leakage factor (0 = no leakage, ~0.999 typical)
 */
void lms_init(LMSFilter *filter, int n_taps, double mu, double leakage);

/**
 * Process one sample through the LMS adaptive filter.
 *
 * y[n] = sum(w_k[n] * x[n-k])
 *
 * @param filter  LMS filter
 * @param x       Input sample (reference signal)
 * @param desired Desired output (for weight update)
 * @return Filter output
 */
double lms_process(LMSFilter *filter, double x, double desired);

/**
 * Normalized LMS update (prevents gradient amplification with large inputs).
 *
 * @param filter  LMS filter
 * @param x       Input sample
 * @param desired Desired output
 * @param epsilon Regularization (small positive)
 * @return Filter output
 */
double lms_process_normalized(LMSFilter *filter, double x,
                               double desired, double epsilon);

/**
 * Get current filter weights.
 */
void lms_get_weights(const LMSFilter *filter, double *weights_out);

/**
 * Free LMS filter memory.
 */
void lms_free(LMSFilter *filter);

/* ==========================================================================
 * L8: Iterative Learning Control
 * ========================================================================== */

/**
 * Initialize an ILC controller.
 *
 * @param ilc      ILC controller to initialize
 * @param n_samples Number of time samples per trial
 * @param dt       Sample time
 * @param learn_q  Q-filter numerator (low-pass)
 * @param learn_l  L-filter (learning gain)
 */
void ilc_init(ILCController *ilc, int n_samples, double dt,
              const TransferFn *learn_l, const TransferFn *learn_q);

/**
 * Compute the feedforward signal for the next trial.
 *
 * u_next = Q(s) * [u_current + L(s) * e_current]
 *
 * @param ilc          ILC controller
 * @param u_current    Current trial feedforward (n_samples)
 * @param error        Current trial tracking error (n_samples)
 * @param u_next       Output next trial feedforward (n_samples)
 */
void ilc_update(ILCController *ilc,
                const double *u_current,
                const double *error,
                double *u_next);

/**
 * Check ILC convergence criterion.
 *
 * |Q(jw) * (1 - L(jw)*P(jw))| < 1 for all frequencies
 *
 * @param ilc     ILC controller
 * @param plant   Plant transfer function P(s)
 * @param omega   Frequency array (rad/s)
 * @param n_freq  Number of frequencies
 * @param margin  Output convergence margin (< 1 means convergence)
 */
void ilc_convergence_check(const ILCController *ilc,
                           const TransferFn *plant,
                           const double *omega, int n_freq,
                           double *margin);

/**
 * Free ILC controller memory.
 */
void ilc_free(ILCController *ilc);

/* ==========================================================================
 * L8: Nonlinear Feedforward
 * ========================================================================== */

/**
 * Compute feedforward torque via inverse dynamics.
 *
 * tau_ff = M(q)*qdd_d + C(q,qd)*qd + G(q) + F(qd)
 *
 * @param model     Nonlinear model
 * @param q         Desired joint positions (n_dof)
 * @param qd        Desired joint velocities (n_dof)
 * @param qdd       Desired joint accelerations (n_dof)
 * @param tau_ff    Output feedforward torques (n_dof)
 */
void nl_ff_computed_torque(const NonlinearFFModel *model,
                           const double *q, const double *qd,
                           const double *qdd, double *tau_ff);

/**
 * Compute feedforward for a simple pendulum (single DOF).
 *
 * Dynamics: J*qdd + b*qd + m*g*l*sin(q) = tau
 * Feedforward: tau_ff = J*qdd_d + b*qd_d + m*g*l*sin(q_d)
 *
 * @param J     Inertia (kg*m^2)
 * @param b     Damping coefficient (N*m*s/rad)
 * @param mgl   m*g*l product (N*m)
 * @param q_d   Desired position (rad)
 * @param qd_d  Desired velocity (rad/s)
 * @param qdd_d Desired acceleration (rad/s^2)
 * @return Feedforward torque (N*m)
 */
double nl_ff_pendulum(double J, double b, double mgl,
                       double q_d, double qd_d, double qdd_d);

/**
 * Compute feedforward for a 2-DOF planar robot arm.
 *
 * Standard two-link manipulator dynamics:
 *   tau1 = (m1+m2)*a1^2*qdd1 + m2*a2^2*(qdd1+qdd2)
 *          + m2*a1*a2*(2*qdd1+qdd2)*cos(q2)
 *          - m2*a1*a2*qd2*(2*qd1+qd2)*sin(q2)
 *          + (m1+m2)*g*a1*cos(q1) + m2*g*a2*cos(q1+q2)
 *
 * @param m1   Mass of link 1 (kg)
 * @param m2   Mass of link 2 (kg)
 * @param a1   Length of link 1 (m)
 * @param a2   Length of link 2 (m)
 * @param q    Joint positions [q1, q2] (rad)
 * @param qd   Joint velocities [qd1, qd2] (rad/s)
 * @param qdd  Joint accelerations [qdd1, qdd2] (rad/s^2)
 * @param tau  Output torques [tau1, tau2] (N*m)
 */
void nl_ff_2link_arm(double m1, double m2, double a1, double a2,
                     const double q[2], const double qd[2],
                     const double qdd[2], double tau[2]);

/**
 * Compute velocity and acceleration feedforward for a servo axis.
 *
 * Standard industrial motion control feedforward:
 *   u_ff = Kv_ff * v_d + Ka_ff * a_d + Ks_ff * sign(v_d)
 *
 * where Kv_ff = velocity FF gain (ideally = 1/Kv if Kv is velocity loop gain)
 *       Ka_ff = acceleration FF (ideally = J/Kt for motor inertia/torque const)
 *       Ks_ff = static friction compensation
 *
 * @param Kv_ff   Velocity feedforward gain
 * @param Ka_ff   Acceleration feedforward gain
 * @param Ks_ff   Static friction compensation gain
 * @param v_d     Desired velocity
 * @param a_d     Desired acceleration
 * @return Feedforward control signal
 */
double ff_servo_motion(double Kv_ff, double Ka_ff, double Ks_ff,
                        double v_d, double a_d);

#endif /* FEEDFORWARD_ADAPTIVE_H */
