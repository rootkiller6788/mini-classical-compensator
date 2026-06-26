# mini-pid-theory — PID Controller Theory

**Module Status: COMPLETE ✅**

Proportional-Integral-Derivative (PID) controller theory — the most widely used
control algorithm in industry (>90% of all control loops).

## Knowledge Coverage

| Level | Name | Status | Details |
|-------|------|--------|---------|
| L1 | Definitions | **Complete** | 9 struct/enum types covering all PID forms and parameters |
| L2 | Core Concepts | **Complete** | PID computation, frequency domain, form conversion |
| L3 | Math Structures | **Complete** | Laplace/Fourier domain, polynomial solvers (quadratic-cubic-quartic) |
| L4 | Fundamental Laws | **Complete** | Routh-Hurwitz, gain/phase margins, sensitivity analysis |
| L5 | Algorithms | **Complete** | 8 tuning methods (ZN, CC, AMIGO, IMC, SIMC, TL, CHR) |
| L6 | Canonical Systems | **Complete** | DC motor, temperature, level, flow control examples |
| L7 | Applications | **Complete** | Industrial motor drives, HVAC/process, tank level, flow loops |
| L8 | Advanced Topics | **Complete** | Anti-windup, cascade, FF, gain scheduling, Smith, adaptive PID |
| L9 | Research Frontiers | **Partial** | Auto-tuning, adaptive PID, IFT implemented; fractional PID documented |

**Score: 17/18 — COMPLETE**

## Core Definitions

- **Standard (ISA) form**: u = Kc·[e + (1/Ti)∫e·dt + Td·de/dt]
- **Parallel form**: u = Kp·e + Ki·∫e·dt + Kd·de/dt
- **Series (interacting) form**: C(s) = Kc'·(1 + 1/(Ti'·s))·(1 + Td'·s)
- **2DOF form**: u = Kc·(b·ysp − y) + Ki·∫(ysp−y)dt + Kd·(c·dysp/dt − dy/dt)

## Core Theorems

1. **Routh-Hurwitz Stability**: For PID-controlled 1st-order plant, closed-loop stable iff a₂a₁ > a₃a₀
2. **Bode Sensitivity Integral**: ∫₀^∞ log|S(jω)| dω = 0 for open-loop stable plants (waterbed effect)
3. **Maximum Sensitivity**: Ms = max_ω|S(jω)| — key robustness metric (target: 1.2–2.0)
4. **Internal Model Principle**: Perfect steady-state disturbance rejection requires integrator in loop

## Core Algorithms

1. **pid_compute()** — Discrete PID with derivative filtering & anti-windup
2. **pid_tune_zn_step()** — Ziegler-Nichols step response method (1942)
3. **pid_tune_zn_ultimate()** — Ziegler-Nichols ultimate sensitivity (1942)
4. **pid_tune_simc()** — Skogestad IMC tuning (2003)
5. **pid_tune_amigo()** — Robust Ms-constrained tuning (2004)
6. **pid_tune_imc_lambda()** — Internal Model Control tuning (1986)
7. **pid_closed_loop_poles()** — Analytic pole computation via Cardano/Ferrari
8. **pid_evaluate_fopdt()** — Time-domain performance simulation

## Classic Problems

1. **DC Motor Speed Control** — Second-order integrating plant
2. **Temperature Control** — Lag-dominant with deadtime (tau >> theta)
3. **Level Control** — Integrating process requiring conservative tuning
4. **Flow Control** — Fast first-order process with measurement noise

## Course Alignment

| School | Course | Key PID Topics |
|--------|--------|---------------|
| MIT | 6.302 Feedback Systems | PID design, loop shaping |
| Stanford | ENGR105 Feedback Control | PID tuning, frequency response |
| Berkeley | ME132 Dynamic Systems | PID implementation |
| ETH | 151-0591 Control I | PID theory and practice |
| Cambridge | 3F2 Systems & Control | Classical compensator design |

## Building

```bash
make          # Build library, tests, and examples
make test     # Run comprehensive test suite (27 tests)
make examples # Build all 4 example programs
make run-examples  # Run all examples
```

## File Structure

```
mini-pid-theory/
├── Makefile
├── README.md                    ← This file
├── include/
│   ├── mini-pid-theory.h        Core PID structures and API
│   ├── pid_tuning.h             Tuning method declarations
│   ├── pid_advanced.h           Advanced features (cascade, FF, Smith)
│   └── pid_adaptive.h           Adaptive PID (auto-tune, MIT, ES, IFT)
├── src/
│   ├── pid_core.c               PID computation, frequency domain, stability
│   ├── pid_tuning.c             All 8 tuning method implementations
│   ├── pid_advanced.c           Cascade, FF, gain scheduling, Smith predictor
│   └── pid_adaptive.c           Auto-tuning, MIT rule, extremum-seeking, STR
├── tests/
│   └── test_pid.c               Comprehensive test suite (27 test functions)
├── examples/
│   ├── example_dc_motor.c       DC motor speed control
│   ├── example_temperature.c    Temperature control with method comparison
│   ├── example_level_control.c  Tank level control (integrating process)
│   └── example_flow_control.c   Flow control with setpoint filtering
└── docs/
    ├── knowledge-graph.md       L1-L9 knowledge map
    ├── coverage-report.md       Completion assessment
    ├── gap-report.md            Remaining gaps and priorities
    ├── course-alignment.md      Nine-school curriculum mapping
    └── course-tree.md           Prerequisite and dependency tree
```

## References

1. Astrom, K.J. & Hagglund, T. (1995). *PID Controllers: Theory, Design, and Tuning*. ISA.
2. Ziegler, J.G. & Nichols, N.B. (1942). Optimum Settings for Automatic Controllers. *Trans. ASME*.
3. Astrom, K.J. & Murray, R.M. (2010). *Feedback Systems*. Princeton University Press.
4. Skogestad, S. (2003). Simple analytic rules for model reduction and PID controller tuning. *J. Process Control*.
5. Rivera, D.E., Morari, M. & Skogestad, S. (1986). Internal Model Control. 4. PID Controller Design. *IEC Proc. Des. Dev.*
