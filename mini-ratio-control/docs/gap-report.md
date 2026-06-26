# Gap Report — mini-ratio-control

## Current Gaps

### L7: Applications (Partial → need 2 more)

Currently implemented:
- ✅ Refinery blending (example_blending.c)
- ✅ Combustion control (example_combustion.c)
- ✅ Water treatment pH (example_cascade_ratio.c)

Gaps:
- ⚠️ Pharmaceutical blending application code
- ⚠️ HVAC outdoor air reset application code

### L8: Advanced Topics (Partial → need 1 more)

Currently implemented:
- ✅ Gain-scheduled feedforward (`feedforward_set_gain_schedule`)
- ✅ Adaptive ratio via `ratio_set_gain()` ramp

Gaps:
- ⚠️ Model-predictive ratio control implementation
- ⚠️ Fuzzy ratio control implementation

### L9: Research Frontiers (Partial)

All L9 items are documented only. Acceptable per SKILL.md §6.1.

## Priority Actions
1. Add pharmaceutical blending example (L7)
2. Add HVAC outdoor reset example (L7)
3. Add MPC ratio control implementation (L8)
