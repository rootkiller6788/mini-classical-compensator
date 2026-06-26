# Course Alignment — mini-pid-tuning-ziegler

Mapping to nine-university curriculum for PID tuning and Ziegler-Nichols methods.

| University | Course | Relevant Content | Implementation Mapping |
|-----------|--------|-----------------|----------------------|
| **MIT** | 6.302 Feedback Systems | PID control, Ziegler-Nichols, loop shaping | `zn_step_tune()`, `gm_pm_tune_pid()` |
| **Stanford** | ENGR105 Feedback Control | PID tuning, anti-windup, gain scheduling | `aw_back_calculation()`, `gs_linear_interpolate()` |
| **Berkeley** | ME132 Dynamic Systems | Process identification, FOPDT models | `fopdt_identify_*()` methods |
| **Caltech** | CDS 110 Intro Control | PID structures, stability margins | `gm_pm_compute_margins()`, `pid_freq_response()` |
| **ETH** | 151-0591 Control I | Classical PID tuning, Z-N methods | `zn_step_tune()`, `zn_freq_tune()` |
| **Cambridge** | 3F2 Systems & Control | Frequency response, describing functions | `relay_describing_function()`, `fopdt_freq_response()` |
| **Georgia Tech** | ECE 6550 Nonlinear | Relay auto-tuning, describing function | `relay_simulate_fopdt()`, `relay_autotune_complete()` |
| **Purdue** | ME 575 Industrial Control | Cohen-Coon, IMC, process applications | `cohen_coon_tune()`, `app_tune_temperature()` |
| **Tsinghua** | 自动控制原理 | PID参数整定, 齐格勒-尼科尔斯法 | All zn_*() and cohen_coon_*() functions |

## Topic-by-Course Matrix

| Topic | MIT | Stan | Berk | Calt | ETH | Camb | GT | Purd | THU |
|-------|-----|------|------|------|-----|------|----|------|-----|
| PID forms | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Z-N step response | ✓ | ✓ | ✓ | ✓ | ✓ | — | — | ✓ | ✓ |
| Z-N frequency response | ✓ | ✓ | ✓ | — | ✓ | ✓ | ✓ | — | ✓ |
| Cohen-Coon | — | — | — | — | — | — | — | ✓ | ✓ |
| IMC-based tuning | — | ✓ | — | — | — | — | — | ✓ | — |
| SIMC | — | — | — | — | ✓ | — | — | — | — |
| Relay auto-tuning | — | — | — | — | — | — | ✓ | ✓ | ✓ |
| Gain/Phase margin design | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Anti-windup | — | ✓ | — | — | — | ✓ | — | ✓ | ✓ |
| Gain scheduling | — | ✓ | — | ✓ | ✓ | — | — | ✓ | ✓ |
| FOPDT identification | — | — | ✓ | — | — | — | — | ✓ | ✓ |
