# Knowledge Graph — mini-feedforward-control

## L1: Definitions (Complete)

| # | Definition | C Type | Location |
|---|-----------|--------|----------|
| 1 | Polynomial p(s) | `Poly` | feedforward_core.h |
| 2 | Transfer function G(s) | `TransferFn` | feedforward_core.h |
| 3 | Feedforward controller | `FeedforwardCtrl` | feedforward_core.h |
| 4 | 2-DOF controller | `TwoDOF` | feedforward_core.h |
| 5 | Prefilter | `Prefilter` | feedforward_core.h |
| 6 | Disturbance FF model | `DistFFModel` | feedforward_core.h |
| 7 | FF performance metrics | `FFPerformance` | feedforward_core.h |
| 8 | Input shaper | `InputShaper`, `Impulse` | feedforward_input_shaping.h |
| 9 | Shaper type enum | `ShaperType` | feedforward_input_shaping.h |
| 10 | LMS adaptive filter | `LMSFilter` | feedforward_adaptive.h |
| 11 | ILC controller | `ILCController` | feedforward_adaptive.h |
| 12 | Nonlinear FF model | `NonlinearFFModel` | feedforward_adaptive.h |
| 13 | FIR filter | `FIRFilter` | feedforward_filter.h |
| 14 | IIR filter | `IIRFilter` | feedforward_filter.h |
| 15 | Notch filter | `NotchFilter` | feedforward_filter.h |
| 16 | Butterworth LPF | `ButterworthLP` | feedforward_filter.h |

**Struct count: 16** (>=5 required)

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Feedforward compensation principle | ff_compute_ideal_filter, ff_compute_closed_loop_tf |
| 2 | 2-DOF separation | twodof_init, TwoDOF |
| 3 | Disturbance rejection via FF | ff_compute_dist_compensator |
| 4 | Model-inverse principle | ff_model_inverse |
| 5 | Causality/realizability | ff_is_realizable |
| 6 | Steady-state FF gain | ff_dc_gain |
| 7 | Reference prefiltering | ff_design_prefilter_1st/2nd |
| 8 | Internal Model Principle | twodof_check_imp |
| 9 | Robustness analysis | twodof_robustness_margin |
| 10 | FB-FF complementary roles | twodof_closed_loop_tfs |

**Count: 10 core concepts**

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Polynomial algebra | poly_add, poly_mul, poly_sub, poly_derivative |
| 2 | TF algebra | tf_parallel, tf_series, tf_feedback |
| 3 | Polynomial long division | ff_causal_inverse |
| 4 | Dominant pole expansion | tf_step_response |
| 5 | State-space realization | tf_to_state_space |
| 6 | Frequency response | tf_freq_response |
| 7 | Minimum-phase detection | tf_is_minimum_phase |
| 8 | Diophantine equation | ff_diophantine_design |
| 9 | Horner polynomial evaluation | poly_eval |
| 10 | Bilinear transform | notch_biquad_design |

**Count: 10 structures**

## L4: Fundamental Theorems (Complete)

| # | Theorem | Implementation |
|---|---------|---------------|
| 1 | Perfect Feedforward Theorem | ff_model_inverse |
| 2 | Internal Model Principle | twodof_check_imp |
| 3 | Disturbance Observer Principle | ff_design_dynamic_dist_ff |
| 4 | Small Gain Theorem | twodof_robustness_margin |
| 5 | Causality Constraint | ff_is_realizable |
| 6 | ZV Shaper Theorem | shaper_design_zv |
| 7 | ZVD Derivative Constraint | shaper_design_zvd |
| 8 | ILC Convergence Criterion | ilc_convergence_check |

**Count: 8 theorems**

## L5: Algorithms/Methods (Complete)

24 algorithms implemented: model inverse, ZPETC, causal inverse, Diophantine design,
properness enforcement, reference FF design, static/dynamic disturbance FF, ZV/ZVD/EI
shapers, two-mode shaper, negative shaper, sensitivity curve, time-optimal shaper,
LMS adaptive FF, normalized LMS, ILC update, FIR filter design (window method),
Butterworth LPF, notch filter, velocity-limiting prefilter, S-curve trajectory,
exponential smoothing, gain-scheduled FF.

## L6: Canonical Systems (Complete)

9 canonical systems: DC motor position/velocity servo, temperature control,
crane anti-sway, flexible transmission, pendulum, 2-link robot arm, HVAC 2-DOF,
industrial servo motion. 4 end-to-end examples (>30 lines + main).

## L7: Applications (Complete)

7 applications: CSTR process control, boiler steam pressure, pH neutralization,
satellite attitude control, quadrotor altitude control, HVAC building automation
(ISO 50001 / smart grid), DC motor industrial servo.

## L8: Advanced Topics (Partial+)

5 advanced topics: Adaptive FF (LMS), ILC, nonlinear FF (computed torque),
gain-scheduled FF, robust FF (small-gain analysis).

## L9: Research Frontiers (Partial)

Documented: AI/ML-based FF, safe learning FF, quantum control FF,
multi-agent FF consensus, cyber-physical FF (smart grid context).
Foundation implementations exist (LMS, ILC).
