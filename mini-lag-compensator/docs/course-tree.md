# Course Prerequisite Tree — mini-lag-compensator

## Dependency Graph

```
mini-control-mathematics (0)
  +-- mini-complex-analysis
  +-- mini-laplace-transform
      |
      v
mini-system-modeling (1)
  +-- Transfer functions
      |
      v
mini-time-domain-analysis (2)
  +-- Step/ramp response, stability
      |
      v
mini-root-locus-method (3)
  +-- Root locus construction
      |
      v
mini-frequency-domain (4)
  +-- Bode plots, Nyquist criterion, stability margins
      |
      v
mini-classical-compensator (5)  <-- YOU ARE HERE
  +-- mini-lag-compensator      <-- THIS MODULE
  +-- mini-lead-compensator
  +-- mini-lead-lag-design
  +-- mini-pid-theory
  +-- mini-pid-tuning-ziegler
      |
      v
mini-state-space-theory (6)
mini-pole-placement-observer (7)
mini-kalman-estimation (8)
mini-stochastic-control (9)
mini-automation-application-systems (10)
```

## Internal Dependencies

```
lag_types.h
  +-- lag_compensator.h/.c (core)
        +-- lag_design.h/.c (design)
        |     +-- lag_frequency.h/.c (Bode/Nyquist)
        +-- lag_simulation.c (time-domain)
        +-- lag_optimization.c (optimization)
        +-- lag_digital.c (digital impl)
        +-- lag_application.c (applications)
  +-- lag_identification.h/.c (system ID)
```

## Prerequisites for This Module

1. Complex analysis (L3): s-domain evaluation
2. Transfer functions (L1): rational polynomial representation
3. Time-domain response (L2): step response, stability
4. Root locus (L5): pole-zero configuration effects
5. Frequency domain (L4): Bode, Nyquist, stability margins