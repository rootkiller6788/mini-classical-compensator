# Course Tree — mini-ratio-control

## Prerequisites
- mini-system-modeling (transfer functions, Laplace domain)
- mini-time-domain-analysis (step response, settling time)
- mini-frequency-domain (Bode plots, lead-lag design)
- mini-root-locus-method (gain margin, stability)

## Dependency Graph

```
mini-system-modeling
  └── mini-time-domain-analysis
        └── mini-frequency-domain
              └── mini-root-locus-method
                    └── mini-ratio-control ← this module
                          ├── ratio_control_core (L1-L2)
                          ├── ratio_station (L1-L5)
                          ├── blending_control (L1-L6)
                          ├── cross_limiting (L1-L4)
                          ├── ratio_feedforward (L1-L8)
                          └── ratio_cascade (L1-L8)
```

## Postrequisites (modules that depend on this)
- mini-state-space-theory (MIMO ratio control)
- mini-pole-placement-observer (ratio with observer)
- mini-automation-application-systems (process control systems)

## L9 Research Frontiers
1. **AI-Based Ratio Optimization**: Reinforcement learning for ratio setpoint in blending
2. **Digital Twin for Ratio Systems**: Real-time ratio monitoring with predictive maintenance
3. **Quantum Control**: Ratio control for quantum state preparation (theoretical)
