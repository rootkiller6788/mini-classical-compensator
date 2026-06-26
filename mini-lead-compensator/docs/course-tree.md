# Course Tree — mini-lead-compensator

## Prerequisites
1. **Transfer functions** (L1) — mini-system-modeling
2. **Frequency response** (L3) — mini-frequency-domain
3. **Bode plots** (L3) — mini-bode-plot
4. **Nyquist criterion** (L4) — mini-nyquist-criterion
5. **Root locus** (L5) — mini-root-locus-method
6. **Phase/gain margin** (L4) — mini-gain-phase-margin

## Dependencies
`
mini-system-modeling
    └── mini-frequency-domain
            ├── mini-bode-plot
            ├── mini-nyquist-criterion
            ├── mini-gain-phase-margin
            └── mini-lead-compensator ← THIS MODULE
                    ├── mini-lag-compensator
                    ├── mini-lead-lag-design
                    └── mini-pid-theory
`

## Sequential Build Order
0 → 1 → 2 → 3 → 4 → **5 (here)** → 6 → 7 → 8 → 9 → 10
