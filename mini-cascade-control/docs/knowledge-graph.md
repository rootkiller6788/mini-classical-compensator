# Knowledge Graph -- mini-cascade-control

## L1: Definitions COMPLETE
- Cascade system structure (inner/secondary + outer/primary loops)
- CascadePoly: polynomial ring R[s] representation
- CascadeTF: rational transfer function G(s) = K*N(s)/D(s)
- CascadeDTF: discrete-time transfer function G(z)
- CascadePID: parallel-form PID controller with anti-windup
- CascadeInner/CascadeOuter: inner/outer loop specifications
- CascadeSystem: complete two-loop cascade system
- CascadeAWMethod: anti-windup method enumeration
- CascadeBumpless: bumpless transfer state machine
- CascadeDisturbance: disturbance input model
- CascadeSensor: sensor/measurement dynamics
- CascadeDesignSpec: design specifications container
- CascadeFreqPoint/CascadeFreqResponse: Bode/Nyquist data
- CascadePerformance: time-domain performance metrics
- CascadeDCMotor/CascadeReactor/CascadeFlowPressure/CascadeLevelTank: canonical models
- CascadeIMCParams/CascadeGainSchedule/CascadeSmithPredictor: advanced types (L8)

## L2: Core Concepts COMPLETE
- Bandwidth separation principle (wbi >= 3*wbo)
- Cascade disturbance rejection mechanism
- Equivalent plant concept (Geq = Gi_cl * Go)
- Sequential loop closure design philosophy
- Anti-windup: clamping, back-calculation, incremental, combined
- Bumpless transfer: manual/auto/cascade mode switching
- Cascade-specific anti-windup: freeze outer I when inner saturates

## L3: Mathematical Structures COMPLETE
- Polynomial ring R[s]: addition, multiplication (convolution)
- Horner method for complex polynomial evaluation O(n)
- Rational function field: transfer function algebra
- Complex frequency domain: s = sigma + j*omega
- Routh array: structured stability test matrix

## L4: Fundamental Laws COMPLETE
- Routh-Hurwitz stability criterion (1874/1895)
- Internal stability theorem for cascade systems
- Sequential loop closure theorem
- Bandwidth separation condition: wbi >= 3*wbo
- Ms-sensitivity bounds

## L5: Computational Methods COMPLETE
- Frequency-domain PI design (magnitude + phase conditions)
- Skogestad SIMC tuning rule for inner PI
- Direct Synthesis tuning for outer PI
- Nelder-Mead simplex optimization (4D)
- Velocity/position form PID
- Tustin (bilinear) discretization
- Bode analysis, stability margins, Ms computation
- Step response: exact 1st/2nd order, dominant-pole approximation
- Performance indices: ISE, IAE, ITAE, TV
- Rise time, settling time, overshoot extraction

## L6: Canonical Systems COMPLETE
- DC motor position/velocity cascade
- Jacketed CSTR temperature cascade
- Flow-pressure pipeline cascade
- Level-on-flow surge tank cascade

## L7: Applications COMPLETE
- Servo motor control (robotics, CNC, aerospace)
- Chemical reactor temperature control
- Surge tank level control (boiler, distillation)

## L8: Advanced Topics PARTIAL
- IMC cascade parameters (struct defined)
- Gain-scheduled cascade (struct defined)
- Smith predictor in cascade (struct defined)
- IMC cascade design algorithm (pending)
- Adaptive cascade (pending)

## L9: Research Frontiers PARTIAL
- Distributed cascade over networks
- Learning-based cascade optimization
- MPC-cascade architectures
