# Course Tree — mini-pid-tuning-ziegler

## Prerequisite Dependencies

```
mini-control-mathematics (0)
    ├── Laplace transform → PID transfer function analysis
    ├── Complex analysis → Frequency response, Nyquist
    └── Numerical methods → Optimization, root-finding

mini-system-modeling (1)
    ├── Transfer function models → FOPDT, SOPDT, IPDT
    └── Step response analysis → Process reaction curve

mini-time-domain-analysis (2)
    ├── Stability definitions → Margin verification
    ├── Time-domain specifications → Overshoot, settling time
    └── Step response analysis → Tuning evaluation

mini-root-locus-method (3)
    └── Pole placement intuition → Cohen-Coon derivation

mini-frequency-domain (4)
    ├── Bode/Nyquist plots → Gain/phase margin analysis
    ├── Describing function → Relay auto-tuning theory
    └── Stability margins → Gm/Pm specification design

mini-classical-compensator (5) ← WE ARE HERE
    ├── mini-pid-theory → PID structure foundations
    ├── mini-pid-tuning-ziegler ← Current module
    ├── mini-lead-compensator
    ├── mini-lag-compensator
    ├── mini-lead-lag-design
    ├── mini-feedforward-control
    ├── mini-cascade-control
    └── mini-ratio-control
```

## Internal Dependency Tree

```
pid_tuning.h (core types)
├── ziegler_nichols.h ──────────┐
├── cohen_coon.h ───────────────┤
├── imc_tuning.h ───────────────┼── all depend on pid_tuning.h
├── relay_autotune.h ───────────┤
├── gain_margin_tuning.h ───────┤
├── fopdt_model.h ──────────────┤
└── advanced_tuning.h ──────────┘

src/ dependencies:
├── pid_tuning_core.c ─── pid_tuning.h
├── ziegler_nichols.c ─── ziegler_nichols.h, pid_tuning.h
├── cohen_coon.c ──────── cohen_coon.h, ziegler_nichols.h (zn_verify_margins)
├── imc_tuning.c ──────── imc_tuning.h
├── relay_autotune.c ──── relay_autotune.h, ziegler_nichols.h (zn_freq_tune)
├── gain_margin_tuning.c ─ gain_margin_tuning.h, ziegler_nichols.h, imc_tuning.h
├── fopdt_model.c ─────── fopdt_model.h
├── advanced_tuning.c ─── advanced_tuning.h, ziegler_nichols.h, cohen_coon.h, imc_tuning.h
└── application_tuning.c ─ all headers (dispatcher uses all methods)

examples/ dependencies:
├── example_temperature_control.c ─ all headers
├── example_auto_tune_sim.c ──────── pid_tuning.h, ziegler_nichols.h, relay_autotune.h, fopdt_model.h, gain_margin_tuning.h
└── example_fopdt_identification.c ─ pid_tuning.h, fopdt_model.h, ziegler_nichols.h, cohen_coon.h, imc_tuning.h, advanced_tuning.h
```

## Knowledge Flow

```
Process Data (step response)
    │
    ▼
FOPDT Identification (fopdt_model.h/c)
    │  tangent, two-point, area, LS, PEM, Sundaresan
    │
    ▼
Process Model (K, T, L)
    │
    ├──► Z-N Step Method ──────────► PID Parameters
    ├──► Cohen-Coon ────────────────► PID Parameters
    ├──► IMC / SIMC ────────────────► PID Parameters
    │
    └──► Relay Experiment ──► Ku, Pu ──► Z-N Freq Method ──► PID Parameters
    │
    └──► Specified Gm/Pm ───────────► Gm/Pm Design ───────► PID Parameters
                │
                ▼
    Anti-Windup + Setpoint Filtering + Gain Scheduling
                │
                ▼
        Closed-Loop Performance Evaluation
```
