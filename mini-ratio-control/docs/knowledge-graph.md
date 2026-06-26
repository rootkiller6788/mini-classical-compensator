# Knowledge Graph — mini-ratio-control

## Module: Ratio Control
Classical Compensator Layer 5.5 — Ratio Station Design & Control Architectures

---

## L1: Definitions (Complete)

| # | Definition | Implementation |
|---|-----------|----------------|
| 1 | Ratio Control: y2(t) = K * y1(t) + b | `ratio_control.h`: ratio_config_t |
| 2 | Ratio Gain K: the proportional factor between master and slave | `ratio_config_t.ratio_gain` |
| 3 | Ratio Bias b: offset applied to ratio relationship | `ratio_config_t.ratio_bias` |
| 4 | Ratio Station: ISA-S5.1 FY function block | `ratio_station.h`: ratio_station_t |
| 5 | Ratio Mode: master-slave, cross-limiting, blending, etc. | `ratio_mode_t` enum |
| 6 | Ratio Formula: linear, square-root, polynomial, exp, log | `ratio_formula_t` enum |
| 7 | Ratio Characterization: nonlinear ratio as function of load | `ratio_char_table_t` |
| 8 | Cross-Limiting: high/low select for combustion safety | `cross_limit_t` struct |
| 9 | Blending System: multi-component mixing with mass balance | `blend_system_t` struct |
| 10 | Feedforward Compensator: ratio-based disturbance rejection | `feedforward_compensator_t` |
| 11 | Cascade Ratio: nested PID loops with ratio station | `cascade_ratio_t` struct |
| 12 | Ratio Integrity: loop health status enumeration | `ratio_integrity_t` enum |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Master-Slave Architecture | `ratio_set_master()`, `ratio_set_slave()` |
| 2 | Ratio Filtering | First-order low-pass on master PV |
| 3 | Ratio Rate Limiting | Rate-of-change limits on slave SP |
| 4 | Ratio Shedding | `ratio_shed()` safe fallback |
| 5 | Ratio Clamping | Min/max ratio enforcement |
| 6 | Dynamic Compensation (Lead-Lag) | `ratio_lead_lag()` |
| 7 | Deadtime Compensation | `deadtime_compensator()` |
| 8 | Square Root Extraction | `ratio_square_root_extract()` |
| 9 | High/Low/Median Signal Selection | `hi_select()`, `lo_select()`, `mid_select()` |
| 10 | Blend Recipe Management | `blend_load_recipe()`, `blend_save_recipe()` |
| 11 | PID with Anti-Windup | `pid_loop_t` with integral clamping |
| 12 | Bumpless Transfer | `pid_loop_bumpless_init()` |
| 13 | Split-Range Ratio Control | `CASCADE_MODE_SPLIT_RANGE` |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Linear Fractional Transform | Ratio = K * master + b |
| 2 | Piecewise-Linear Interpolation | `ratio_char_lookup()` via binary search |
| 3 | First-Order Discrete Filter | y[k] = y[k-1] + alpha*(x[k] - y[k-1]) |
| 4 | Discrete Transfer Function (Direct Form I) | `ratio_dynamic_compensator_step()` |
| 5 | Mass Conservation: sum(f_i) = 1.0 | `blend_mass_balance_error()` |
| 6 | Linear Programming (Vertex Enumeration) | `blend_linear_program_solve()` |
| 7 | Bilinear (Tustin) Discretization | `lead_lag_filter()` |
| 8 | Ring Buffer Deadtime Delay | `deadtime_compensator()` |

## L4: Fundamental Laws (Complete)

| # | Theorem/Law | Code Assertion |
|---|------------|----------------|
| 1 | Ratio Station Linearity: SP = K * PV + b | `test_ratio_set_master_multiply` |
| 2 | Mass Balance: sum(f_i) = 1.0 | `test_blend_mass_balance` |
| 3 | Cross-Limit Safety: fuel_sp <= air / K_fa | `test_cross_limit_fuel_limited` |
| 4 | Cross-Limit Excess Air: air_sp >= fuel * K_fa | `test_cross_limit_air_guaranteed` |
| 5 | Lead-Lag Steady-State Gain = 1.0 | `test_feedforward_lead_lag` |
| 6 | Filter Convergence: lim_{t->inf} filt(t) = PV | `test_ratio_set_master_filter` |
| 7 | Rate Limit Bound: |ΔSP| <= rate_limit * dt | `test_ratio_station_rate_limit` |
| 8 | Cascade Primary Anti-Windup | `pid_loop_execute()` integral clamping |
| 9 | Ratio Deviation: K_actual = slave_pv / master_pv | `test_ratio_set_slave` |
| 10 | Optimal Blend: LP optimal at vertex | `test_blend_optimization` |

## L5: Computational Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Ratio Setpoint Computation | `ratio_set_master()` |
| 2 | Ratio Characterization Interpolation | Binary search + linear interp |
| 3 | Blend Cost Optimization (Greedy) | `blend_optimize_cost()` |
| 4 | Blend LP (Vertex Enumeration) | `blend_linear_program_solve()` |
| 5 | Cross-Limit High/Low Select Algorithm | `cross_limit_execute()` |
| 6 | Signal Selection (Lo/Hi/Median) | `signal_selector_execute()` |
| 7 | Dynamic Lead-Lag Compensation | `lead_lag_filter()` |
| 8 | Deadtime Ring Buffer | `deadtime_compensator()` |
| 9 | Ratio Gain Ramp | `ratio_set_gain()` with rate limiting |
| 10 | Blend Recipe Normalization | `blend_normalize_fractions()` |

## L6: Canonical Systems (Complete)

| # | System | Example |
|---|--------|---------|
| 1 | Gasoline Blending (Octane/RVP) | `example_blending.c` |
| 2 | Boiler Combustion Cross-Limiting | `example_combustion.c` |
| 3 | pH Neutralization Cascade Ratio | `example_cascade_ratio.c` |
| 4 | Distillation Reflux Ratio | cascade ratio with temperature/flow |
| 5 | Chemical Dosing Ratio | master-slave ratio station |
| 6 | HVAC Outdoor Air Reset | feedforward ratio based on outdoor temp |

## L7: Applications (Partial+ — 3 implemented)

| # | Application | Reference |
|---|------------|-----------|
| 1 | Refinery Gasoline Blending | `example_blending.c` — ISO 8217 fuel standards |
| 2 | Power Plant Burner Management | `example_combustion.c` — Fukushima safety lesson |
| 3 | Water Treatment pH Control | `example_cascade_ratio.c` — NHS/Detroit Water |
| 4 | Pharmaceutical Blending | blend_system_t with quality constraints |
| 5 | HVAC Outdoor Reset | ratio_feedforward with gain scheduling |

## L8: Advanced Topics (Partial+ — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Adaptive Ratio Control | Config hooks via `ratio_set_gain()` ramp |
| 2 | Gain-Scheduled Feedforward | `feedforward_set_gain_schedule()` |
| 3 | Model-Predictive Ratio Control | Documented approach |
| 4 | Fuzzy Ratio Control | Documented approach |
| 5 | Multi-variable Ratio Decoupling | `cascade_ratio_t` for dual-output |

## L9: Research Frontiers (Partial — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | AI-Based Ratio Optimization | Documented in gap-report.md |
| 2 | Digital Twin for Ratio Systems | Documented in course-tree.md |
| 3 | Reinforcement Learning Ratio Tuning | Research direction noted |
