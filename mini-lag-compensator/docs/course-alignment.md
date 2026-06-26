# Course Alignment — mini-lag-compensator

## Nine-School Curriculum Mapping

### MIT — 6.302 Feedback Systems
- Frequency-response design methods -> lag_design_bode_method()
- Steady-state errors -> lag_design_for_steady_state_error()
- Root-locus design -> lag_design_root_locus()

### Stanford — ENGR105 Feedback Control
- Lag/Lead Design: Frequency-domain loop shaping -> lag_frequency.c
- Digital Implementation: Discrete-time equivalents -> lag_digital.c

### Berkeley — ME132/232 Dynamic Systems
- Classical compensator design -> all design functions
- Multi-objective trade-offs -> lag_pareto_optimization()

### Caltech — CDS 110 Introduction to Control
- Bode/Nyquist/Root-locus methods -> lag_frequency.c
- State-space connection -> lag_simulation.c

### ETH Zurich — 151-0591 Control I
- Frequenzkennlinienverfahren (Bode-based design)
- Korrekturglieder (correction elements)

### Cambridge — 3F2 Systems & Control
- Compensator Synthesis: frequency-domain methods
- Robustness: gain/phase margin analysis

### Georgia Tech — ECE 6550 / AE 6530 / ME 6401
- Aerospace Applications -> lag_design_aerospace_actuator()

### Purdue — ME 575 / ME 675
- Industrial Process Control -> lag_design_temperature_control()
- Multivariable Extensions -> multiple SISO lag loops

### Tsinghua — 自动控制原理
- 滞后校正 (Lag compensator design methodology)
- 频域设计 (Frequency-domain design)