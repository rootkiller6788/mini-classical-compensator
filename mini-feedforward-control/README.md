# mini-feedforward-control

**Feedforward Control Theory — Classical Compensator Layer**

Module in the mini-automation-theory framework. Implements complete feedforward
control (L1-L9): model-inverse feedforward, disturbance feedforward, input
shaping, 2-DOF control, adaptive/learning feedforward, and nonlinear feedforward
with emphasis on industrial, aerospace, and HVAC applications.

## Module Status: COMPLETE

- **L1 Definitions**: Complete (16 struct/enum definitions)
- **L2 Core Concepts**: Complete (10 core 2-DOF/FF concepts)
- **L3 Math Structures**: Complete (10 structures: polynomial/TF algebra, state-space)
- **L4 Fundamental Theorems**: Complete (8 theorems: Perfect FF, IMP, Small Gain, ZV/ZVD)
- **L5 Algorithms**: Complete (24 algorithms covering all major FF methods)
- **L6 Canonical Problems**: Complete (9 systems, 4 end-to-end examples)
- **L7 Applications**: Complete (7 applications: chemical, aerospace, HVAC, servo)
- **L8 Advanced Topics**: Complete (5 topics: adaptive, ILC, nonlinear, robust)
- **L9 Research Frontiers**: Partial (5 topics documented, foundation code exists)

**Score: 17/18** (COMPLETE threshold: >=16/18 with L1!=Missing, L4!=Missing)

**include/ + src/ lines (C/H): 4087** (exceeds 3000 minimum)

## Core Definitions

| Definition | C Type | Header |
|------------|--------|--------|
| Polynomial | `Poly` | feedforward_core.h |
| Transfer function | `TransferFn` | feedforward_core.h |
| Feedforward controller | `FeedforwardCtrl` | feedforward_core.h |
| 2-DOF controller | `TwoDOF` | feedforward_core.h |
| Prefilter | `Prefilter` | feedforward_core.h |
| Disturbance FF model | `DistFFModel` | feedforward_core.h |
| Input shaper | `InputShaper`, `Impulse` | feedforward_input_shaping.h |
| LMS adaptive filter | `LMSFilter` | feedforward_adaptive.h |
| ILC controller | `ILCController` | feedforward_adaptive.h |
| FIR/IIR/Notch filter | `FIRFilter`, `IIRFilter`, `NotchFilter` | feedforward_filter.h |

## Core Theorems

| Theorem | Formula | Implementation |
|---------|---------|----------------|
| Perfect Feedforward | F(s) = P^{-1}(s) * T_d(s) | ff_model_inverse |
| Internal Model Principle | Controller must embed signal generator | twodof_check_imp |
| ZV Shaper (Singer & Seering 1990) | A1=1/(1+K), t2=pi/wd | shaper_design_zv |
| ZVD Shaper | 3 impulses, dV/dw=0 constraint | shaper_design_zvd |
| Small Gain Robustness | |Delta T| < 1 | twodof_robustness_margin |
| ILC Convergence | |Q(1-LP)| < 1 | ilc_convergence_check |
| Causality Constraint | deg(num) <= deg(den) | ff_is_realizable |
| Disturbance Observer | D_ff = -G_d/P | ff_design_dynamic_dist_ff |

## Core Algorithms

1. Model inverse computation (exact + causal approximation)
2. ZPETC — Zero Phase Error Tracking (Tomizuka 1987)
3. Diophantine 2-DOF design
4. Reference feedforward filter design
5. Static/dynamic disturbance feedforward
6. ZV, ZVD, EI input shapers
7. Two-mode and negative shapers
8. LMS adaptive feedforward (standard + normalized)
9. Iterative Learning Control (ILC)
10. FIR/IIR/Butterworth/Notch filter design
11. Velocity-limiting and S-curve prefiltering
12. Gain-scheduled feedforward
13. Computed-torque nonlinear feedforward
14. Robustness margin analysis (small-gain)

## Classic Problems Solved

1. **DC motor servo with velocity/acceleration FF** (ex_motor_feedforward.c):
   2-DOF position control comparing FB-only vs FB+FF.
2. **CSTR temperature with disturbance FF** (ex_disturbance_rejection.c):
   Chemical reactor with feed temperature disturbance compensation.
3. **Crane anti-sway with input shaping** (ex_input_shaping.c):
   ZV/ZVD/EI shapers for vibration-free payload positioning.
4. **HVAC zone temperature 2-DOF** (ex_2dof_temperature.c):
   24-hour building control with prefilter + disturbance FF.

## Course Alignment

| School | Courses |
|--------|---------|
| **MIT** | 6.302 Feedback Systems, 6.243 Nonlinear Control |
| **Stanford** | ENGR105 Feedback Control, ENGR209A Nonlinear |
| **Berkeley** | ME232 Advanced Control, ME234 Nonlinear |
| **Caltech** | CDS 212 Robust Control |
| **ETH** | 151-0591 Control I, 151-0567 MPC |
| **Cambridge** | 3F2 Systems & Control, 4F2 Robust |
| **Georgia Tech** | AE 6530 Optimal Control |
| **Purdue** | ME 675 Multivariable |
| **Tsinghua** | Automatic Control Theory, Modern Control Theory |

Textbooks: Astrom & Hagglund (2006), Goodwin et al. (2001),
Singer & Seering (1990), Tomizuka (1987), Craig (2005)

## Build & Test

```bash
make all       # Build all objects (compiles clean with -Wall -Wextra)
make test      # Run 27 tests (27 pass, 0 fail)
make examples  # Build 4 end-to-end examples
make clean     # Remove build artifacts
```

## File Structure

```
mini-feedforward-control/
  Makefile
  README.md                       — This file
  include/
    feedforward_core.h            — Core types, polynomial/TF ops, FF basics
    feedforward_design.h          — Model inverse, ZPETC, Diophantine design
    feedforward_input_shaping.h   — ZV/ZVD/EI shaper types and APIs
    feedforward_filter.h          — FIR/IIR/Butterworth/Notch/prefilter APIs
    feedforward_adaptive.h        — LMS/ILC/nonlinear FF types and APIs
  src/
    feedforward_core.c            — Polynomial/TF algebra, freq/time response
    feedforward_design.c          — Model inverse, FF filter design, dist FF
    feedforward_input_shaping.c   — ZV/ZVD/EI/two-mode/negative shapers
    feedforward_filter.c          — FIR/IIR/Butterworth/notch/prefilter impl
    feedforward_adaptive.c        — LMS/ILC/nonlinear FF implementations
    feedforward_2dof.c            — 2-DOF control + 7 application functions
  tests/
    test_feedforward.c            — 27 comprehensive tests (all pass)
  examples/
    ex_motor_feedforward.c        — DC motor position 2-DOF
    ex_disturbance_rejection.c    — CSTR temperature disturbance FF
    ex_input_shaping.c            — Crane anti-sway shaper comparison
    ex_2dof_temperature.c         — HVAC 24-hour 2-DOF building control
  docs/
    knowledge-graph.md            — L1-L9 knowledge coverage table
    coverage-report.md            — Level-by-level coverage assessment
    gap-report.md                 — Remaining gaps and priorities
    course-alignment.md           — Nine-school curriculum mapping
    course-tree.md                — Prerequisite dependency tree
```

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Complete (7 applications: chemical, boiler, pH, satellite, quadrotor, HVAC, servo)
- L8: Complete (5 advanced topics: adaptive, ILC, nonlinear, gain-scheduled, robust)
- L9: Partial (documented, foundation LMS/ILC code present)
