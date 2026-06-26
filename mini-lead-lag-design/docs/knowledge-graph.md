# Knowledge Graph ！ mini-lead-lag-design

## L1 ！ Definitions (Complete)
- LeadCompensator: C(s)=Kc*(1+T*s)/(1+alpha*T*s), 0<alpha<1
- LagCompensator: C(s)=Kc*(1+T*s)/(1+beta*T*s), beta>1
- LeadLagCompensator: combined lead+lag cascade
- FrequencyDomainSpec, TimeDomainSpec, CompensatorSpec
- BodeData, PlantFreqData, FirstOrderPlant, SecondOrderPlant
- OpenLoopPZ, DominantPoleSpec
- DigitalFilter, DigitalPID
- SensitivityData, LoopShape, LoopSpec
- FOPDTModel, SOPDTModel

## L2 ！ Core Concepts (Complete)
- Phase-lead principle: pole farther left than zero on real axis
- Phase-lag principle: zero farther left than pole (dipole near origin)
- Lead-Lag combination: independent design with frequency separation
- Bode gain-phase relationship, waterbed effect (S+T=1)
- Loop shaping for performance/robustness trade-off
- Internal Model Control (IMC) principle
- Direct Synthesis method

## L3 ！ Mathematical Structures (Complete)
- Complex frequency response: C(j*omega)
- Bode magnitude/phase, dB/logarithmic scales
- Pole-zero maps on s-plane
- First-order and second-order transfer function models
- Rational function approximations (Pade for delays)
- Sensitivity S(s) and complementary sensitivity T(s)
- Digital z-domain transfer functions

## L4 ！ Fundamental Laws (Complete)
- Maximum phase lead: phi_max = arcsin((1-alpha)/(1+alpha))
- Frequency of max phase: omega_m = 1/(T*sqrt(alpha))
- Alpha from phase: alpha = (1-sin(phi))/(1+sin(phi))
- Lag attenuation: 20*log10(1/beta) dB
- Root locus angle/magnitude conditions
- PO-damping relation: zeta = -ln(PO/100)/sqrt(pi^2+ln^2(PO/100))
- PM-damping relation: PM ~ 100*zeta (engineering approx)
- Sensitivity peak to GM/PM bounds
- Bode integral (waterbed): integral of log|S| = 0
- Skogestad half-rule for model reduction

## L5 ！ Computational Methods (Complete)
- Lead design: Ogata method (phase margin targeting)
- Lag design: steady-state error reduction
- Lead-lag design: independent lead+lag with freq separation
- Bode data generation (log-spaced frequencies)
- Stability margin computation (PM, GM, crossovers)
- Root locus: angle condition, magnitude condition, breakaway
- Tustin/bilinear discretization
- Zero-order hold discretization
- Digital filter application (difference equation)
- Sensitivity function computation from loop data
- IMC design for FOPDT plants
- SIMC tuning rules (Skogestad 2003)
- DS/lambda tuning for PID
- Loop shape design and fitting

## L6 ！ Canonical Problems (Complete)
- DC motor speed control with lead compensator
- Temperature process with lag compensator
- Servo position control with lead-lag
- Antenna tracking with frequency-domain design

## L7 ！ Applications (Partial+)
- Industrial motion control (DC motor)
- Process control (temperature)
- Aerospace servo systems
- Antenna/radar tracking systems

## L8 ！ Advanced Topics (Partial+)
- H-infinity loop shaping principles
- Mixed sensitivity design
- McFarlane-Glover stability margin
- Adaptive loop shaping targets

## L9 ！ Research Frontiers (Partial)
- Documented: online loop shaping, learning-based compensator tuning
