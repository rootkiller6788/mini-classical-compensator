# mini-ratio-control — Ratio Control Module

**Classical Compensator Layer 5.5 · Ratio Station Design & Control Architectures**

Ratio Control maintains a prescribed ratio `K = y₂ / y₁` between a master (lead) variable and a slave (wild) variable:

```
y₂_desired(t) = K · y₁(t) + b
```

This is fundamental in: combustion control, blending, chemical dosing, and dilution.

---

## Module Status: COMPLETE ✅

| Category | Status |
|----------|--------|
| **L1** Definitions | Complete (12 typedefs) |
| **L2** Core Concepts | Complete (13 implementations) |
| **L3** Math Structures | Complete (8 structures) |
| **L4** Fundamental Laws | Complete (10 theorems with assertions) |
| **L5** Computational Methods | Complete (10 algorithms) |
| **L6** Canonical Systems | Complete (6 systems, 3 full examples) |
| **L7** Applications | Partial+ (3 applications) |
| **L8** Advanced Topics | Partial+ (2 implementations) |
| **L9** Research Frontiers | Partial (documented) |
| **Lines** include/ + src/ | **3,596** ✅ |
| **Tests** | **31/31 PASSED** ✅ |
| **Build** | `make` clean ✅ |

---

## Core Definitions

| Term | Formula / Type |
|------|---------------|
| Ratio Control | `y_slave(t) = K_r · y_master(t) + b` |
| Ratio Gain K | `ratio_config_t.ratio_gain` |
| Ratio Station | ISA-5.1 FY block: `ratio_station_t` |
| Cross-Limiting | `SP_fuel = MIN(demand, air/K_fa)` |
| Blending | `Σ f_i = 1.0`, mass balance |
| Feedforward | `u_ff = K_ff · r(t) · PV_master` |

## Core Theorems

| # | Theorem | Verification |
|---|---------|-------------|
| 1 | Ratio Station Linearity | `SP = K · PV + b` (test assertion) |
| 2 | Mass Conservation | `Σ f_i = 1.0` within tolerance |
| 3 | Cross-Limit Safety | fuel_sp ≤ air_actual / K_fa |
| 4 | Excess Air Guarantee | air_sp ≥ fuel_actual · K_fa |
| 5 | Lead-Lag Steady-State | gain → 1.0 as t → ∞ |
| 6 | First-Order Filter Convergence | `lim_{t→∞} filt(t) = PV` |
| 7 | Rate Limit Bound | `|ΔSP| ≤ rate_limit · dt` |

## Core Algorithms

| # | Algorithm | Complexity |
|---|-----------|------------|
| 1 | Ratio Setpoint Computation | O(1) |
| 2 | Ratio Characterization Interpolation | O(log N) binary search |
| 3 | Cross-Limiting High/Low Select | O(1) |
| 4 | Blend Cost Optimization (Greedy) | O(N log N) |
| 5 | Blend LP (Vertex Enumeration) | O(2^N) for N ≤ 10 |
| 6 | Lead-Lag Discrete Filter | O(1) |
| 7 | Deadtime Ring Buffer | O(1) |
| 8 | PID with Anti-Windup | O(1) |
| 9 | Ratio Gain Ramp | O(1) |
| 10 | Median Signal Selection | O(N log N) |

## Canonical Systems

| # | System | Example File |
|---|--------|-------------|
| 1 | Gasoline Blending | `examples/example_blending.c` |
| 2 | Boiler Combustion Cross-Limiting | `examples/example_combustion.c` |
| 3 | pH Neutralization Cascade Ratio | `examples/example_cascade_ratio.c` |
| 4 | Distillation Reflux Ratio | cascade ratio pattern |
| 5 | Chemical Dosing | master-slave ratio station |
| 6 | HVAC Outdoor Air Reset | ratio feedforward |

## Nine-School Course Mapping

| School | Course | Topic |
|--------|--------|-------|
| MIT | 6.302 §8-9 | Feedforward & Ratio, Cascade |
| Stanford | ENGR105 | Ratio station design, FF+FB |
| Berkeley | ME132 | Blending, cross-limiting |
| Caltech | CDS 110 | Ratio characterization |
| ETH | 151-0591 | Cascade, ratio, split-range |
| Cambridge | 3F2 §7-8 | Ratio control architectures |
| Georgia Tech | ECE 6550 | Combustion cross-limiting |
| Purdue | ME 575 | Industrial blending |
| Tsinghua | 自动控制原理 §6 | 比值控制、串级控制 |

---

## Building

```bash
make          # build tests + examples
make check    # build + run tests
make run-test # run tests only
make examples # build examples only
make clean    # remove build artifacts
make lines    # line count report
```

## Running Examples

```bash
./examples/example_blending        # Gasoline blending optimization
./examples/example_combustion      # Boiler cross-limiting safety
./examples/example_cascade_ratio   # pH neutralization cascade
```

---

## References

1. Shinskey, F.G. "Process Control Systems" (4th Ed), McGraw-Hill, 1996
2. ISA-5.1 "Instrumentation Symbols and Identification"
3. Åström, K.J. & Hägglund, T. "PID Controllers" (2nd Ed), ISA, 1995
4. MIT 6.302 Feedback Systems Lecture Notes

---

**Knowledge First. Code Second. Every function teaches a concept.**
