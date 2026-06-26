# mini-lead-lag-design

**Classical compensator | mini-lead-lag-design | mini-classical-compensator (5/6)**

## Module Status: COMPLETE

- **L1-L6**: Complete
- **L7**: Complete (4 applications: DC motor, temperature, servo, antenna tracking)
- **L8**: Partial (H-infinity loop shaping, mixed sensitivity; QFT documented)
- **L9**: Partial (adaptive tuning via extremum seeking documented, not implemented)

---

## Overview

Lead-Lag compensator combines a phase-lead network (transient improvement)
and a phase-lag network (steady-state accuracy) in cascade:

    C(s) = Kc * (T_lead*s + 1)/(alpha*T_lead*s + 1) * (T_lag*s + 1)/(beta*T_lag*s + 1)

Full pipeline: specification -> analytical design -> frequency-domain verification
-> root locus analysis -> digital implementation -> sensitivity analysis.

## Knowledge Coverage (L1-L9)

### L1 - Core Definitions
| Type | Description |
|------|-------------|
| LeadLagCompensator | Combined lead-lag: Kc, T_lead, alpha, T_lag, beta |
| LeadLagDesignSpec | Design spec: PM_desired, PM_current, wc, ess |
| LeadLagFreqPoint | Complex frequency response: w, mag(dB), phase |
| CompensatorSpec | General compensator spec with freq/time domain |
| FrequencyDomainSpec | PM, GM, bandwidth, crossover requirements |
| TimeDomainSpec | Overshoot %, settling time, steady-state error |
| LoopShape | Target loop transfer function shape |
| SensitivityData | |S|, |T|, peak values, bandwidth metrics |
| DigitalFilter | Tustin/bilinear transformed filter coefficients |
| FOPDTModel / SOPDTModel | First-order / second-order plant models |

### L2 - Core Concepts
- **Phase-lead principle**: pole farther from origin than zero, giving phase boost near wm
- **Phase-lag principle**: zero farther from origin than pole, raising DC gain
- **Lead-Lag combination**: independent design with frequency separation (w_lead >> w_lag)
- **Bode gain-phase relationship**: minimum-phase systems, waterbed effect (S+T=1)
- **Loop shaping**: performance/robustness trade-off via open-loop shaping
- **Internal Model Control (IMC)**: model-based design for robustness
- **Direct Synthesis (DS)**: specify desired closed-loop response, solve for controller

### L3 - Mathematical Structures
- Complex frequency response: C(jw), plant(jw), L(jw) = C(jw)*plant(jw)
- Bode magnitude (dB) and phase (deg) on logarithmic frequency axis
- Pole-zero maps on s-plane: lead zero/pole pair, lag zero/pole pair
- FOPDT and SOPDT transfer function models with Pade delay approximation
- Sensitivity S(s) = 1/(1+L(s)), complementary sensitivity T(s) = L(s)/(1+L(s))
- Z-domain digital transfer functions via Tustin and ZOH transforms

### L4 - Fundamental Laws and Formulas
| Formula | Description |
|---------|-------------|
| phi_max = arcsin((1-alpha)/(1+alpha)) | Maximum phase lead from alpha |
| w_m = 1/(T*sqrt(alpha)) | Frequency of maximum phase lead |
| alpha = (1-sin(phi_max))/(1+sin(phi_max)) | Required alpha for desired phase boost |
| Lag attenuation = 20*log10(1/beta) dB | High-frequency attenuation from lag |
| zeta approx PM/100 | Engineering damping approximation |
| PO = 100*exp(-pi*zeta/sqrt(1-zeta^2)) | Percent overshoot from damping |
| ts approx 4/(zeta*wn) | 2% settling time |
| Integral(ln|S(jw)|) dw = 0 | Bode sensitivity integral (waterbed) |
| angle(L(jwgc)) + 180 deg = PM | Phase margin definition |
| Skogestad half-rule | Model reduction: tau1+tau2/2 to effective delay |

### L5 - Algorithms and Methods
| Algorithm | Method | Reference |
|-----------|--------|-----------|
| lead_design_ogata() | Ogata phase-margin targeting lead design | Ogata Ch.7 |
| lag_design_ess() | Lag design for steady-state error specification | Ogata Ch.7 |
| lead_lag_design_direct() | Independent lead+lag with freq separation | Franklin Ch.6 |
| lead_lag_bode_data() | Log-spaced Bode frequency response generation | - |
| compensator_freq_design() | Frequency-domain compensator synthesis | Nise Ch.9 |
| compensator_rlocus_design() | Root locus angle/magnitude condition design | Ogata Ch.7 |
| lead_lag_discretize_tustin() | Bilinear (Tustin) transform s->2/Ts*(z-1)/(z+1) | Franklin Ch.8 |
| imc_design_fopdt() | IMC-based design for FOPDT plants | Rivera-Morari-Skogestad 1986 |
| simc_tune() | SIMC tuning rules for PID | Skogestad 2003 |
| ds_lambda_tune() | Direct synthesis lambda-tuning | Chen and Seborg 2002 |
| loop_shape_fit() | Fit compensator to target loop shape | - |
| sensitivity_compute() | S(s), T(s) from loop frequency response | - |

### L6 - Canonical Problems (examples/)
| Example | Plant | Compensator | Domain |
|---------|-------|-------------|--------|
| ex_dc_motor_lead | DC motor (2nd order) | Lead compensator | Motion control |
| ex_temp_control_lag | Temperature process (FOPDT) | Lag compensator | Process control |
| ex_servo_lead_lag | Servo mechanism | Lead-lag | Precision positioning |
| ex_antenna_tracking | Antenna azimuth | Frequency-domain | Tracking systems |

### L7 - Applications
1. **Industrial motion control**: DC motor speed regulation with lead compensator
2. **Process control**: Temperature loop with lag + lead-lag disturbance rejection
3. **Aerospace servo systems**: Position servos with fast response + zero steady-state error
4. **Antenna/radar tracking**: Frequency-domain design for tracking bandwidth requirements

### L8 - Advanced Topics
- **H-infinity loop shaping**: McFarlane-Glover framework for robust stability margin
- **Mixed sensitivity design**: Simultaneous shaping of S(s) and T(s)
- **Adaptive loop shaping targets**: Performance-adaptive compensator re-design
- QFT (Quantitative Feedback Theory) - documented, not implemented
- Gain scheduling with lead-lag - documented, not implemented

### L9 - Industry Frontiers (Documented)
- Online loop shaping with real-time system identification
- Learning-based compensator tuning via RL/ILC
- Distributed compensator coordination in networked control
- Cyber-physical compensator verification

## Course Alignment

| School | Course | Topics Covered |
|--------|--------|---------------|
| MIT | 6.302 Feedback Systems | Lead/lag/lead-lag design, Bode methods, root locus |
| Stanford | ENGR105 Feedback Control | Frequency-domain compensator synthesis |
| Berkeley | ME132 Dynamic Systems | Bode-based lead design, loop shaping |
| Caltech | CDS 110 Intro to Control | Classical compensator analysis and synthesis |
| ETH | 151-0591 Control I | Phase-lead/lag networks (Korrekturglieder) |
| Cambridge | 3F2 Systems and Control | Compensator design methods |
| Georgia Tech | ECE 6550 | Classical lead-lag design procedures |
| Purdue | ECE 602 Lumped Systems | Compensator synthesis |
| Tsinghua | Automatic Control Theory | Lead/lag correction |

## Building

```bash
cd mini-lead-lag-design
make all      # build everything
make test     # run all 23 tests
make examples # build example programs
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-lead-lag-design/
+-- include/
|   +-- lead_lag_design.h          # Core lead-lag compensator struct and API
|   +-- lead_design.h              # Lead-only design algorithms
|   +-- lag_design.h               # Lag-only design algorithms
|   +-- compensator_spec.h         # Frequency/time domain specifications
|   +-- compensator_freq_design.h  # Frequency-domain design methods
|   +-- compensator_rlocus_design.h# Root locus design methods
|   +-- compensator_analytical.h   # Analytical compensator formulas
|   +-- compensator_loopshaping.h  # Loop shaping design
|   +-- compensator_sensitivity.h  # S/T sensitivity analysis
|   +-- compensator_digital.h      # Digital implementation (Tustin/ZOH)
+-- src/                           # Implementation (10 .c files)
+-- tests/                         # test_lead_lag.c (23 tests)
+-- examples/                      # 4 end-to-end examples
+-- docs/
|   +-- knowledge-graph.md         # L1-L9 knowledge map
|   +-- coverage-report.md         # Level-by-level assessment
|   +-- gap-report.md              # Missing items and priority
|   +-- course-alignment.md        # 9-school course mapping
|   +-- course-tree.md             # Prerequisite dependency tree
+-- Makefile
```

## License

MIT