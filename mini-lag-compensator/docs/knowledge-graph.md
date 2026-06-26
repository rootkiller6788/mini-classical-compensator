# Knowledge Graph — mini-lag-compensator

## L1: Definitions (15 struct/enum types)

| Item | C Implementation | Description |
|------|-----------------|-------------|
| LagCompensator | struct { Kc, T, beta, zero, pole, dc_gain, hf_gain, max_phase_lag_rad, max_lag_freq } | Core lag compensator |
| LagComplex | struct { re, im } | Complex number |
| LagFreqPoint | struct { omega, magnitude, magnitude_db, phase_rad, phase_deg, real_part, imag_part } | Frequency response point |
| LagDesignSpec | struct { PM_target, GM_target, ess_target, ess_type, ... } | Design specification |
| LagStepResponse | struct { time, output, error, control_signal, ... } | Step response trajectory |
| LagRampResponse | struct { time, output, error, ... } | Ramp response trajectory |
| LagTransferFunction | struct { numerator, denominator, dc_gain } | Rational TF |
| LagPolynomial | struct { order, coeff } | Polynomial |
| LagBodeData | struct { points, num_points, phase_margin, gain_margin, ... } | Bode plot |
| LagNyquistData | struct { points, encirclements, is_stable } | Nyquist plot |
| LagOptResult | struct { optimal_Kc, optimal_T, optimal_beta, converged } | Optimization result |
| LagDCMotorParams | struct { R, L, Kb, Kt, J, B, ... } | DC motor model |
| LagESSType | enum { STEP, RAMP, PARABOLIC } | ESS type |
| LagClosedLoop | struct { plant, compensator, open_loop, closed_loop, sensitivity, comp_sensitivity } | Closed-loop system |
| LagDigital | struct { b0, b1, a1, e_prev, u_prev, Ts } | Digital compensator |

## L2: Core Concepts (10 items)

| Concept | Implementation | Source |
|---------|---------------|--------|
| DC gain | lag_get_dc_gain() | lag_compensator.c |
| HF gain | lag_get_hf_gain() | lag_compensator.c |
| Corner frequencies | lag_get_corner_frequencies() | lag_compensator.c |
| Maximum phase lag | lag_get_max_phase_lag() | lag_compensator.c |
| ESS improvement | lag_ess_improvement() | lag_compensator.c |
| Bandwidth | lag_compute_bandwidth() | lag_frequency.c |
| Low-freq asymptote | lag_low_freq_asymptote_db() | lag_compensator.c |
| High-freq asymptote | lag_high_freq_asymptote_db() | lag_compensator.c |
| Lag principle | Design docs + implementation | lag_design.c |
| Separation factor beta | lag_get_beta() | lag_compensator.c |

## L3: Mathematical Structures (7 items)

| Structure | Implementation | Source |
|-----------|---------------|--------|
| Complex s-domain eval | lag_eval_s() | lag_compensator.c |
| Rational TF | LagTransferFunction, lag_to_transfer_function() | lag_compensator.c |
| Polynomial convolution | build_open_loop_tf() (static) | lag_frequency.c |
| Frequency response | lag_eval_frequency() | lag_compensator.c |
| Bilinear transform | lag_to_digital_bilinear() | lag_digital.c |
| ZOH transform | lag_to_digital_zoh() | lag_digital.c |
| Matched pole-zero | lag_to_digital_matched() | lag_digital.c |

## L4: Fundamental Theorems (8 theorems)

| Theorem | Implementation | Test |
|---------|---------------|------|
| Nyquist Stability (Z=N+P) | lag_compute_nyquist() | test_nyquist_stable |
| Phase Margin definition | lag_compute_phase_margin() | test_phase_margin |
| Gain Margin definition | lag_compute_gain_margin() | Implemented |
| Final Value Theorem for ESS | lag_compute_steady_state_error() | lag_frequency.c |
| Maximum Phase Lag Theorem | lag_get_max_phase_lag() | test_max_phase_lag |
| Routh-Hurwitz (compensator) | lag_is_stable() | test_stability |
| Minimum-phase check | lag_is_minimum_phase() | test_stability |
| Error constant estimation | lag_estimate_error_constant() | lag_frequency.c |

## L5: Computational Methods (29 algorithms)

Bode/Nyquist computation (6), design algorithms (6), optimization (4),
digital discretization (3), time-domain simulation (4), system identification (6).

## L6: Canonical Systems (3 systems)

- DC motor speed control — lag_design_dc_motor_speed() / ex_dc_motor.c
- Temperature control (FOPDT) — lag_design_temperature_control() / ex_temp_control.c
- Position servo (type-1) — lag_design_position_servo() / ex_position_servo.c

## L7: Applications (6 domains)

- Power system frequency regulation (smart grid, nuclear, ISO)
- Chemical process control (Beer, supplier, ISO)
- Aerospace flight actuation (Boeing, F-35, Tesla, SpaceX, NASA, Quadrotor)
- Automotive cruise control (Toyota, Tesla, Detroit)
- Manufacturing quality (ISO, Six Sigma)
- DC motor control (Tesla, Toyota)

## L8: Advanced Topics (Partial+)

- Digital implementation (3 methods) — Complete
- Anti-windup — Complete
- Multi-objective Pareto optimization — Complete
- Pre-warping for bilinear transform — Complete
- Sample rate selection — Complete
- Gradient-based optimization — Complete

## L9: Research Frontiers (Partial)

- Safe RL with lag compensator baselines — Documented
- Cyber-physical system compensation — Documented
- Adaptive lag compensation — Documented
- Quantum control applications — Documented
- Human-machine shared control — Documented