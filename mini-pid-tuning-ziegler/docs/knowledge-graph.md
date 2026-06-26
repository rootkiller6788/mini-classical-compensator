# Knowledge Graph — mini-pid-tuning-ziegler

## L1: Definitions (Complete ✅)

| ID | Definition | C Implementation | Lean Formalization |
|----|-----------|-----------------|-------------------|
| 1.1 | PID controller forms (ideal/parallel/series) | `pid_form_t` enum, `pid_params_t` struct | `PIDForm`, `PIDParams` |
| 1.2 | FOPDT model definition | `fopdt_model_t` struct | `FOPDTModel` |
| 1.3 | SOPDT model definition | `sopdt_model_t` struct | — |
| 1.4 | IPDT model definition | `ipdt_model_t` struct | — |
| 1.5 | Tuning method enumeration | `pid_tune_method_t` enum | — |
| 1.6 | Anti-windup mode enumeration | `antiwindup_mode_t` enum | `AntiWindupMode` |
| 1.7 | Performance metrics (IAE, ISE, ITAE) | `pid_perf_t` struct | `computeIAE`, `computeISE`, `computeITAE` |
| 1.8 | Z-N step response rule table | `zn_rule_entry_t` struct | `ZNStepRule` |
| 1.9 | Z-N frequency response rule table | `zn_freq_rule_entry_t` struct | `ZNFreqRule` |
| 1.10 | Gain/Phase margin specification | `margin_spec_t` struct | — |
| 1.11 | IMC tuning parameters | `imc_tuning_result_t` struct | — |
| 1.12 | Relay auto-tuner configuration | `relay_config_t`, `relay_result_t` | — |
| 1.13 | Gain schedule entry/table | `gs_entry_t`, `gs_table_t` | — |
| 1.14 | Cohen-Coon decay ratio | `cc_controller_type_t` enum | — |

## L2: Core Concepts (Complete ✅)

| ID | Concept | Implementation |
|----|---------|---------------|
| 2.1 | Feedback principle (error-based control) | `pid_update()` — error computation |
| 2.2 | Integral windup phenomenon | `aw_back_calculation()`, `aw_clamping()` |
| 2.3 | Derivative kick on setpoint change | `PID_FORM_DERIVATIVE_ON_PV`, setpoint weighting |
| 2.4 | Bumpless transfer | `pid_reset()`, `bumpless_config_t` |
| 2.5 | Setpoint tracking vs regulation | CHR rules: `CHR_RULES_SP` vs `CHR_RULES_DIST` |
| 2.6 | Stability margin concept | `margin_spec_t`, `gm_pm_compute_margins()` |
| 2.7 | Quarter-amplitude damping (Z-N target) | Z-N rule coefficient derivation |
| 2.8 | Describing function (relay analysis) | `relay_describing_function()` |
| 2.9 | Process reaction curve | `zn_identify_fopdt()` — tangent method |
| 2.10 | Robustness-performance tradeoff (λ tuning) | `imc_lambda_from_Ms()` |
| 2.11 | Controllability (dead time ratio) | `fopdt_dead_time_ratio()`, `fopdt_controllability_index()` |
| 2.12 | Cascade control structure | `cascade_tune_pid()` |
| 2.13 | Feedforward compensation | `feedforward_pid_design()` |

## L3: Mathematical Structures (Complete ✅)

| ID | Structure | Implementation |
|----|-----------|---------------|
| 3.1 | Laplace domain: PID transfer function | `pid_transfer_function_coeffs()`, `pid_freq_response()` |
| 3.2 | FOPDT in Laplace domain | `fopdt_freq_response()` |
| 3.3 | Complex frequency response analysis | PID + FOPDT magnitude/phase computation |
| 3.4 | Padé approximation of dead time | Used in IMC derivation (documented) |
| 3.5 | Optimal parameter estimation (least squares) | `fopdt_identify_ls()` |
| 3.6 | Prediction error minimization | `fopdt_identify_pem()` |
| 3.7 | Linear interpolation (gain scheduling) | `gs_linear_interpolate()` |
| 3.8 | Exponential weighting (smooth scheduling) | `gs_exponential_interpolate()` |

## L4: Fundamental Laws & Theorems (Complete ✅)

| ID | Theorem/Law | C Verification | Lean Statement |
|----|------------|---------------|----------------|
| 4.1 | Ziegler-Nichols step response rules (1942) | `zn_step_tune()` | `zn_step_pid_rule_spec` |
| 4.2 | Ziegler-Nichols frequency response rules (1942) | `zn_freq_tune()` | `zn_freq_pid_rule_spec` |
| 4.3 | Cohen-Coon tuning formulas (1953) | `cohen_coon_tune()` | `cohen_coon_pid_formula` |
| 4.4 | Chien-Hrones-Reswick rules (1952) | `zn_chien_hrones_reswick()` | — |
| 4.5 | IMC-PID equivalence (Rivera-Morari-Skogestad 1986) | `imc_tune_fopdt()`, `imc_tune_sopdt()` | — |
| 4.6 | SIMC rules (Skogestad 2003) | `imc_simc_tune()` | — |
| 4.7 | Nyquist stability criterion (margin computation) | `gm_pm_compute_margins()` | `gain_margin_definition`, `phase_margin_definition` |
| 4.8 | Relay oscillation condition: N(a)*G(jω)=-1 | `relay_describing_function()`, `relay_extract_ultimate()` | `relay_describing_function_magnitude`, `relay_ultimate_gain` |
| 4.9 | Skogestad half-rule model reduction | `fopdt_skogestad_half_rule()` | — |
| 4.10 | Sørensen ultimate gain method | `zn_sorensen_ultimate_gain()` | — |
| 4.11 | Åström-Hägglund auto-tuning (1984) | `relay_autotune_complete()` | — |
| 4.12 | IDEAL↔PARALLEL↔SERIES form equivalence | `pid_convert_form()` | `ideal_parallel_equivalence`, `series_to_parallel_conversion` |
| 4.13 | Ho-Hang-Cao margin-based PI/PID design (1995) | `gm_pm_tune_pi()`, `gm_pm_tune_pid()` | — |
| 4.14 | Anti-windup back-calculation (Fertik & Ross 1967) | `aw_back_calculation()` | `antiwindup_back_calculation` |
| 4.15 | Gain schedule interpolation convexity | `gs_linear_interpolate()` | `gain_schedule_interpolation_bounds` |

## L5: Computational Methods & Algorithms (Complete ✅)

| ID | Algorithm | Implementation | Complexity |
|----|-----------|---------------|------------|
| 5.1 | Tangent method FOPDT identification | `fopdt_identify_graphical()` | O(N) |
| 5.2 | Two-point method (28.3%/63.2%) identification | `fopdt_identify_two_point()` | O(N) |
| 5.3 | Area method identification | `fopdt_identify_area()` | O(N) |
| 5.4 | Least-squares FOPDT identification | `fopdt_identify_ls()` | O(N·L_steps) |
| 5.5 | PEM FOPDT identification | `fopdt_identify_pem()` | O(N·iter) |
| 5.6 | Sundaresan-Krishnaswamy identification | `fopdt_identify_sundaresan()` | O(N) |
| 5.7 | Z-N step response tuning | `zn_step_tune()`, `zn_step_tune_modified()` | O(1) |
| 5.8 | Z-N frequency response tuning | `zn_freq_tune()`, `zn_freq_tune_custom()` | O(1) |
| 5.9 | Ultimate gain numerical search | `zn_find_ultimate_gain()` | O(log(ω_max/ω_min)) |
| 5.10 | Cohen-Coon analytical tuning | `cohen_coon_tune()`, variants | O(1) |
| 5.11 | IMC-PID for FOPDT/SOPDT/IPDT | `imc_tune_fopdt()` etc. | O(1) |
| 5.12 | SIMC (Skogestad) PI/PID tuning | `imc_simc_tune()` | O(1) |
| 5.13 | Relay auto-tuning simulation | `relay_simulate_fopdt()` | O(N_steps) |
| 5.14 | Variable hysteresis relay method | `relay_variable_hysteresis()` | O(attempts·N_steps) |
| 5.15 | Gain/Phase margin computation | `gm_pm_compute_margins()` | O(log(ω_max/ω_min)²) |
| 5.16 | Maximum sensitivity (Ms) calculation | `gm_pm_max_sensitivity()` | O(N_ω) |
| 5.17 | Ms-constrained tuning | `gm_pm_tune_by_Ms()` | O(N_ω) |
| 5.18 | Margin-based PI/PID design | `gm_pm_tune_pi()`, `gm_pm_tune_pid()` | O(iter) |
| 5.19 | PID update with anti-windup | `pid_update()` | O(1) |
| 5.20 | Velocity-form PID | `aw_velocity_form()` | O(1) |
| 5.21 | Setpoint filtering | `sp_filter_update()` | O(1) |
| 5.22 | Derivative filtering | `pid_derivative_filter()` | O(1) |
| 5.23 | Linear gain schedule interpolation | `gs_linear_interpolate()` | O(n) |
| 5.24 | Exponential gain schedule interpolation | `gs_exponential_interpolate()` | O(n) |
| 5.25 | Cascade PID tuning procedure | `cascade_tune_pid()` | O(1) |
| 5.26 | Feedforward parameter design | `feedforward_pid_design()` | O(1) |
| 5.27 | PID method comparison/ranking | `pid_compare_tuning_methods()` | O(n_methods) |
| 5.28 | Optimal setpoint weight computation | `sp_compute_optimal_weights()` | O(1) |
| 5.29 | FOPDT response prediction (decay/overshoot) | `cohen_coon_predict_response()` | O(1) |
| 5.30 | FOPDT simulation (step & arbitrary input) | `fopdt_simulate_step()`, `fopdt_simulate_arbitrary()` | O(N) |

## L6: Canonical Systems (Complete ✅)

| ID | System | Implementation |
|----|--------|---------------|
| 6.1 | Temperature control (lag-dominant) | `app_tune_temperature()` |
| 6.2 | Flow control (fast, noisy) | `app_tune_flow()` |
| 6.3 | Pressure control | `app_tune_pressure()` |
| 6.4 | Level control (integrating, surge) | `app_tune_level()` |
| 6.5 | pH control (nonlinear) | `app_tune_ph()` |
| 6.6 | DC motor speed control | `app_tune_dc_motor()` |
| 6.7 | HVAC zone temperature | `app_tune_hvac()` |
| 6.8 | Chemical reactor (exothermic CSTR) | `app_tune_chemical_reactor()` |
| 6.9 | Paper machine headbox | `app_tune_paper_headbox()` |

## L7: Applications (Complete ✅)

| ID | Application | Implementation |
|----|------------|---------------|
| 7.1 | Industrial process control (temperature/flow/level/pressure) | `app_auto_tune()` dispatcher |
| 7.2 | PID auto-tuner (relay-based) | `relay_autotune_complete()` |
| 7.3 | Adaptive re-tune trigger | `relay_needs_retune()` |
| 7.4 | Model identification pipeline | `example_fopdt_identification.c` |
| 7.5 | Tuning method recommendation engine | `app_recommend_method()`, `app_recommend_adaptation()` |
| 7.6 | Quick FOPDT estimation from bump test | `app_estimate_fopdt()` |

## L8: Advanced Topics (Partial ⚠️)

| ID | Topic | Implementation | Status |
|----|-------|---------------|--------|
| 8.1 | Cascade control tuning | `cascade_tune_pid()` | ✅ Complete |
| 8.2 | Feedforward + PID design | `feedforward_pid_design()` | ✅ Complete |
| 8.3 | IMC for unstable processes | `imc_tune_unstable()` | ✅ Complete |
| 8.4 | Maximum sensitivity (Ms) design | `gm_pm_tune_by_Ms()`, `gm_pm_max_sensitivity()` | ✅ Complete |
| 8.5 | Iso-damping (flat-phase) PID | `gm_pm_tune_pid()` with `use_flat_phase` | ✅ Partial |
| 8.6 | Gain scheduling (linear + exponential) | `gs_linear_interpolate()`, `gs_exponential_interpolate()` | ✅ Complete |
| 8.7 | Fractional-order PID | — | ❌ Not implemented |
| 8.8 | Fuzzy-PID | — | ❌ Not implemented |

## L9: Research Frontiers (Partial ⚠️)

| ID | Topic | Documentation | Status |
|----|-------|--------------|--------|
| 9.1 | Self-tuning regulators (STR) | `SelfTuningRegulator` in Lean | Partial |
| 9.2 | Model-free adaptive PID | — | Documented only |
| 9.3 | Extremum-seeking PID tuning | — | Documented only |
| 9.4 | Multi-agent consensus control | — | Not covered |
