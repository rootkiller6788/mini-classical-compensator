# Mini Classical Compensator

A collection of **from-scratch, zero-dependency C implementations** of classical compensator design and PID control theory. Each module maps to MIT (and other top-tier university) courses, bridging textbook transfer functions and frequency-domain design methods into runnable C code.

## Module Status: COMPLETE ✅

| Sub-Module | include/+src/ Lines | L1-L6 | L7 | L8 | L9 | Status |
|------------|---------------------|-------|----|----|----|--------|
| mini-cascade-control | 3,194 | Complete | Complete | Partial | Partial | ✅ |
| mini-feedforward-control | 4,061 | Complete | Complete | Partial | Partial | ✅ |
| mini-lag-compensator | 5,362 | Complete | Complete | Partial | Partial | ✅ |
| mini-lead-compensator | 3,404 | Complete | Complete | Partial | Partial | ✅ |
| mini-lead-lag-design | 3,009 | Complete | Complete | Partial | Partial | ✅ |
| mini-pid-theory | 3,357 | Complete | Complete | Partial | Partial | ✅ |
| mini-pid-tuning-ziegler | 6,771 | Complete | Complete | Partial | Partial | ✅ |
| mini-ratio-control | 3,607 | Complete | Complete | Partial | Partial | ✅ |

All 8 sub-modules exceed the 3,000-line threshold and pass `make test` with zero failures.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-cascade-control](mini-cascade-control/) | Cascade control architecture, nested loops, bandwidth separation, sequential closure | MIT 6.302, Stanford ENGR105, Berkeley ME232 |
| [mini-feedforward-control](mini-feedforward-control/) | Feedforward compensation, 2-DOF architecture, disturbance rejection, input shaping | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-lag-compensator](mini-lag-compensator/) | Lag compensator, steady-state error reduction, DC gain boost, phase-lag analysis | MIT 6.302, Stanford ENGR105, Berkeley ME132 |
| [mini-lead-compensator](mini-lead-compensator/) | Lead compensator, phase margin improvement, transient response, Bode design | MIT 6.302, Stanford ENGR105, Caltech CDS 110 |
| [mini-lead-lag-design](mini-lead-lag-design/) | Lead-lag synthesis, root locus design, loop shaping, frequency-domain methods, sensitivity | MIT 6.302, Georgia Tech ECE 6550, ETH 151-0591 |
| [mini-pid-theory](mini-pid-theory/) | PID forms, anti-windup, 2-DOF PID, derivative filtering, bumpless transfer | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-pid-tuning-ziegler](mini-pid-tuning-ziegler/) | Ziegler-Nichols tuning, FOPDT identification, relay autotune, Cohen-Coon, IMC | MIT 6.302, Chalmers (Åström & Hägglund) |
| [mini-ratio-control](mini-ratio-control/) | Ratio control, cross-limiting, blending, master-slave, mass balance integrity | Purdue ECE 602, Industrial Process Control |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module translates classical control transfer functions and design procedures into runnable algorithms
- **Practical demos** — compensator design tools, PID autotuners, cascade tuning, ratio station simulators, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-cascade-control
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-classical-compensator/
├── mini-cascade-control/        # Cascade control with nested primary/secondary loops
├── mini-feedforward-control/    # Feedforward and 2-DOF control architecture
├── mini-lag-compensator/        # Phase-lag compensator for steady-state accuracy
├── mini-lead-compensator/       # Phase-lead compensator for transient improvement
├── mini-lead-lag-design/        # Combined lead-lag design and loop shaping
├── mini-pid-theory/             # PID control theory, forms, and advanced features
├── mini-pid-tuning-ziegler/     # Ziegler-Nichols and classic PID tuning methods
└── mini-ratio-control/          # Ratio control with cross-limiting and blending
```

## License

MIT
