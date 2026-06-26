# Course Dependency Tree — mini-pid-theory

## Prerequisites
This module depends on the following prior knowledge:

```
mini-control-mathematics (0)
  └── Laplace transform, complex analysis
      └── mini-system-modeling (1)
          └── Transfer functions, FOPDT models
              └── mini-time-domain-analysis (2)
                  └── Step response, performance metrics
                      └── mini-root-locus-method (3)
                          └── Pole placement intuition
                              └── mini-frequency-domain (4)
                                  ├── Bode plots, Nyquist criterion
                                  └── mini-classical-compensator (5)
                                      └── **mini-pid-theory ← YOU ARE HERE**
```

## What This Module Provides To

```
mini-pid-theory
  ├──→ mini-state-space-theory (6)
  │     PID as output feedback; comparison with state feedback
  ├──→ mini-pole-placement-observer (7)
  │     Observer-based control vs PID
  ├──→ mini-kalman-estimation (8)
  │     PID with Kalman-filtered measurements
  ├──→ mini-stochastic-control (9)
  │     PID performance under stochastic disturbances
  └──→ mini-automation-application-systems (10)
        Real-world PID applications
```

## Internal Dependency Tree

```
pid_core.c (foundation)
  ├── pid_tuning.c (uses frequency response from core)
  ├── pid_advanced.c (uses PID compute + form conversion)
  └── pid_adaptive.c (uses tuning + core compute)

All headers are self-contained with forward declarations.
```

## Key Concept Dependencies

| Concept | Depends On | Used By |
|---------|-----------|---------|
| PID Standard Form | Transfer functions (L3) | All tuning methods |
| Routh-Hurwitz | Characteristic polynomial (L3) | Stability analysis, robustness |
| Gain/Phase Margin | Frequency response (L4) | Robustness evaluation |
| ZN Tuning | Relay feedback (L5) | All auto-tuning methods |
| Anti-windup | Saturation nonlinearity (L2) | Cascade, adaptive PID |
| Smith Predictor | FOPDT model (L6) | Deadtime compensation |
| Adaptive PID | Online identification (L8) | Self-tuning regulators |
