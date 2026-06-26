# Knowledge Graph — mini-lead-compensator

## L1: Definitions
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Lead compensator transfer function C(s) = K_c*(s+z_c)/(s+p_c) | lead_compensator.h | Complete |
| 2 | K-T-alpha parameterization C(s) = K*(Ts+1)/(alpha*Ts+1) | lead_compensator.h | Complete |
| 3 | Lead ratio alpha = z_c/p_c, 0 < alpha < 1 | lead_compensator.h | Complete |
| 4 | lead_compensator_t struct | lead_compensator.h | Complete |
| 5 | lead_tf_t transfer function struct | lead_compensator.h | Complete |
| 6 | lead_system_t plant model struct | lead_compensator.h | Complete |
| 7 | lead_specs_t design specifications struct | lead_compensator.h | Complete |
| 8 | lead_design_result_t design result struct | lead_compensator.h | Complete |
| 9 | lead_performance_t performance metrics struct | lead_compensator.h | Complete |
| 10 | lead_bode_data_t Bode plot data struct | lead_compensator.h | Complete |
| 11 | DC gain K = K_c * alpha | lead_compensator.c | Complete |
| 12 | HF gain = K_c = K / alpha | lead_compensator.c | Complete |

## L2: Core Concepts
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Phase lead phi(w) = atan(wT) - atan(alpha*wT) | lead_compensator.c | Complete |
| 2 | Magnitude |C(jw)| = K*sqrt(1+(wT)^2)/sqrt(1+(alpha*wT)^2) | lead_compensator.c | Complete |
| 3 | Compensated phase: total = angle(C) + angle(G) | lead_compensator.c | Complete |
| 4 | Compensated magnitude: product of magnitudes | lead_compensator.c | Complete |
| 5 | Bandwidth extension via lead | lead_analysis.c | Complete |
| 6 | Noise amplification trade-off | lead_analysis.c | Complete |
| 7 | Control effort analysis | lead_analysis.c | Complete |
| 8 | Disturbance rejection via sensitivity | lead_analysis.c | Complete |
| 9 | Phase efficiency: phi/phi_max | lead_compensator.c | Complete |
| 10 | Phase bandwidth (effective lead range) | lead_compensator.c | Complete |

## L3: Mathematical Structures
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Complex number operations (lead_complex_t) | lead_compensator.c | Complete |
| 2 | Horner polynomial evaluation | lead_frequency.c | Complete |
| 3 | Transfer function evaluation at complex s | lead_frequency.c | Complete |
| 4 | PM-to-damping relationship | lead_compensator.c | Complete |
| 5 | Damping-to-PM relationship | lead_compensator.c | Complete |
| 6 | zeta/wn/PO/ts relationships | lead_compensator.c | Complete |
| 7 | Bandwidth from zeta, wn | lead_compensator.c | Complete |
| 8 | Dominant pole from zeta, wn | lead_design.c | Complete |

## L4: Fundamental Laws
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Maximum phase lead theorem: phi_max = asin((1-alpha)/(1+alpha)) | lead_compensator.c | Complete |
| 2 | Inverse: alpha = (1-sin(phi_m))/(1+sin(phi_m)) | lead_compensator.c | Complete |
| 3 | omega_m = 1/(T*sqrt(alpha)) | lead_compensator.c | Complete |
| 4 | |C(j*omega_m)| = K/sqrt(alpha) | lead_compensator.c | Complete |
| 5 | Nyquist stability criterion (compensated) | lead_frequency.c | Complete |
| 6 | Routh-Hurwitz stability criterion | lead_frequency.c | Complete |
| 7 | Characteristic polynomial: 1+C(s)G(s)=0 | lead_frequency.c | Complete |
| 8 | Sensitivity function S = 1/(1+L) | lead_analysis.c | Complete |
| 9 | Complementary sensitivity T = L/(1+L) | lead_analysis.c | Complete |
| 10 | S+T=1 identity | lead_analysis.c | Complete |

## L5: Computational Methods
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Bode design method (Ogata 7-4) | lead_design.c | Complete |
| 2 | Gain from steady-state error | lead_design.c | Complete |
| 3 | New crossover frequency bisection | lead_design.c | Complete |
| 4 | Root-locus design method | lead_design.c | Complete |
| 5 | Angle deficiency computation | lead_design.c | Complete |
| 6 | Zero/pole placement by bisector | lead_design.c | Complete |
| 7 | Analytic design from PO%, ts | lead_design.c | Complete |
| 8 | Bode plot generation | lead_frequency.c | Complete |
| 9 | Phase margin computation | lead_frequency.c | Complete |
| 10 | Gain margin computation | lead_frequency.c | Complete |
| 11 | Step response simulation (RK4) | lead_analysis.c | Complete |
| 12 | Step response metrics | lead_analysis.c | Complete |
| 13 | Ramp response simulation | lead_analysis.c | Complete |
| 14 | Performance prediction from PM/w_gc | lead_design.c | Complete |

## L6: Canonical Systems
| # | Topic | Example | Status |
|---|-------|---------|--------|
| 1 | DC motor position servo | example_lead_basic.c | Complete |
| 2 | DC motor velocity control | example_lead_dc_motor.c | Complete |
| 3 | Aircraft pitch control (Boeing 747) | example_lead_aero.c | Complete |

## L7: Applications
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Aerospace pitch control (Boeing 747, NASA) | example_lead_aero.c | Complete |
| 2 | Manufacturing automation (DC motor) | example_lead_dc_motor.c | Complete |
| 3 | Servo position control | example_lead_basic.c | Complete |

## L8: Advanced Topics
| # | Topic | Code Location | Status |
|---|-------|---------------|--------|
| 1 | Multi-stage lead compensation | lead_advanced.c | Complete |
| 2 | Discrete-time (Tustin) lead | lead_advanced.c | Complete |
| 3 | Robust lead design (GM+PM) | lead_advanced.c | Complete |
| 4 | Sensitivity-based robust stability | lead_advanced.c | Complete |
| 5 | GM degradation analysis | lead_advanced.c | Complete |

## L9: Research Frontiers
| # | Topic | Status |
|---|-------|--------|
| 1 | Adaptive lead compensation | Partial (documented) |
| 2 | Learning-based lead tuning | Partial (documented) |
| 3 | Quantum control applications | Partial (documented) |
