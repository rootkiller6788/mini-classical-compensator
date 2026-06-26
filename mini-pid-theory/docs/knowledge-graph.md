# Knowledge Graph — mini-pid-theory

## L1: Definitions — COMPLETE
- pid_form_t: Standard (ISA), Parallel, Series, 2DOF PID forms
- pid_action_t: Direct/reverse acting controller
- pid_deriv_mode_t: Derivative on error / measurement / filtered
- pid_error_type_t: Absolute, squared, ITAE-weighted, saturated error
- pid_antiwindup_t: None, clamping, back-calculation, incremental
- pid_params_t: Complete PID parameter structure (Kc, Ti, Td, Ts, N, b, c, limits)
- pid_state_t: Internal state (integrator, derivative memory, performance accumulators)
- pid_2dof_params_t: Two-degree-of-freedom setpoint weights
- pid_performance_t: IAE, ISE, ITAE, ITSE, TV, overshoot, settling/rise time, margins

## L2: Core Concepts — COMPLETE
- PID initialization and validation (pid_init)
- PID control law computation (pid_compute) with all forms
- State reset and bumpless transfer (pid_reset, pid_bumpless_transfer)
- Frequency response (pid_freq_magnitude, pid_freq_phase)
- Form conversion (pid_convert_form)
- Relay feedback for auto-tuning (pid_relay_feedback)
- Closed-loop pole computation (pid_closed_loop_poles)
- Performance evaluation on FOPDT (pid_evaluate_fopdt)

## L3: Mathematical Structures — COMPLETE
- Laplace-domain PID: C(s) = Kc*(1 + 1/(Ti*s) + Td*s)
- Frequency response: |C(jw)|, arg(C(jw))
- Loop transfer function: L(s) = G(s)*C(s)
- Sensitivity function: S(s) = 1/(1+L(s))
- Complementary sensitivity: T(s) = L(s)/(1+L(s))
- Maximum sensitivity Ms via golden-section search
- Polynomial root-finding: quadratic (discriminant), cubic (Cardano), quartic (Ferrari)

## L4: Fundamental Laws — COMPLETE
- Routh-Hurwitz stability criterion for PID-controlled FOPDT
- Gain margin and phase margin via binary frequency search
- Nyquist stability concept via loop transfer analysis
- Bode integral constraints (implicit in Ms computation)

## L5: Computational Methods — COMPLETE
- Ziegler-Nichols step response tuning (1942)
- Ziegler-Nichols ultimate sensitivity tuning (1942)
- Cohen-Coon open-loop tuning (1953)
- AMIGO robust tuning (Astrom & Hagglund, 2004)
- IMC Lambda tuning (Rivera, Morari & Skogestad, 1986)
- SIMC tuning (Skogestad, 2003)
- Tyreus-Luyben conservative tuning (1992)
- Chien-Hrones-Reswick setpoint/disturbance tuning (1952)
- Tuning robustness analysis (gain/deadtime margins)
- Normalized deadtime and controller type recommendation
- Form conversion (standard/parallel/series/2DOF)

## L6: Canonical Systems — COMPLETE
- DC motor speed control (example_dc_motor.c)
- Temperature control of heated tank (example_temperature.c)
- Tank level control (example_level_control.c)
- Flow control loop (example_flow_control.c)

## L7: Applications — COMPLETE
- Industrial motor drives (DC motor example)
- Process temperature control (temperature example)
- Tank level control (level example)
- Flow control loop (flow example)

## L8: Advanced Topics — COMPLETE
- Anti-windup: back-calculation, clamping, conditional integration
- Cascade control (inner/outer loop architecture)
- Feedforward control (static + dynamic lead-lag)
- Gain scheduling (linear interpolation in lookup table)
- Setpoint filtering and rate limiting
- Setpoint weighting optimization (b, c weights)
- Bumpless parameter change for online gain updates
- Smith Predictor for deadtime-dominant processes
- MIT rule adaptive PID (Model Reference Adaptive Control)
- Extremum-seeking auto-tuner (sinusoidal perturbation)
- Iterative Feedback Tuning (IFT) foundation
- Self-Tuning Regulator (RLS + SIMC)

## L9: Research Frontiers — PARTIAL
- Auto-tuning: relay feedback method (Astrom-Hagglund) — IMPLEMENTED
- Adaptive PID: MIT rule, extremum-seeking — IMPLEMENTED
- Data-driven tuning: IFT framework — IMPLEMENTED
- Fractional-order PID: documented only
- Model-free tuning: documented only
