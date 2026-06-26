# Course Tree — mini-feedforward-control

## Prerequisites (Dependency Tree)

```
mini-feedforward-control depends on:
|
+-- 0. mini-control-mathematics
|   +-- mini-complex-analysis (frequency response, transfer functions)
|   +-- mini-laplace-z-transform (Laplace domain, poles/zeros)
|   +-- mini-difference-equations (discrete-time FF, ZPETC)
|
+-- 1. mini-system-modeling
|   +-- Transfer function models (plant identification)
|   +-- State-space models (for nonlinear FF)
|
+-- 2. mini-time-domain-analysis
|   +-- Step/impulse response
|   +-- Performance specifications (ISE, settling time)
|
+-- 3. mini-root-locus-method
|   +-- Pole placement intuition
|   +-- Dominant pole concept
|
+-- 4. mini-frequency-domain
|   +-- Bode/Nyquist analysis
|   +-- Stability margins
|   +-- Bandwidth concepts
|
+-- 5. mini-classical-compensator (sibling modules)
    +-- mini-pid-theory (feedback controller design)
    +-- mini-lead-compensator (phase lead for feedback)
    +-- mini-lag-compensator (steady-state improvement)
    +-- mini-lead-lag-design
    +-- mini-cascade-control (nested loops)
    +-- mini-ratio-control (ratio FF is a form of FF)
    +-- mini-pid-tuning-ziegler
```

## Feedforward-Specific Dependency Tree

```
Feedforward Control Theory
|
+-- Feedforward Fundamentals
|   +-- Model inverse principle
|   +-- 2-DOF architecture
|   +-- Causality/realizability
|
+-- Reference Feedforward
|   +-- Prefilter design (1st order, 2nd order)
|   +-- Model-reference FF
|   +-- Input shaping (ZV, ZVD, EI)
|   +-- S-curve trajectory generation
|
+-- Disturbance Feedforward
|   +-- Static disturbance FF
|   +-- Dynamic disturbance FF
|   +-- Measurable disturbance models
|   +-- Process control applications
|
+-- Advanced Feedforward
    +-- Adaptive FF (LMS)
    +-- Iterative Learning Control (ILC)
    +-- Nonlinear FF (computed torque)
    +-- Gain-scheduled FF
    +-- Robust FF (small-gain)
```

## What This Module Provides to Later Modules

| Later Module | What feedforward provides |
|-------------|--------------------------|
| 6. mini-state-space-theory | State-space 2-DOF design, LQR + FF |
| 7. mini-pole-placement-observer | Observer-based FF, disturbance observer |
| 8. mini-kalman-estimation | Kalman filter for feedforward estimation |
| 9. mini-stochastic-control | Stochastic FF, robust FF under uncertainty |
| 10. mini-automation-application-systems | Full 2-DOF industrial applications |
