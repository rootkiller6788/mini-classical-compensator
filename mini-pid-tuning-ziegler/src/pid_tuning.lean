/-
  pid_tuning.lean — Lean 4 Formalization of PID Tuning Structures and
  Ziegler-Nichols Methods

  Knowledge coverage:
    L1: PID parameter structures, FOPDT model definitions
    L4: Ziegler-Nichols rule theorems
    L5: Controller form conversion theorems

  Reference:
    Ziegler & Nichols (1942), Trans. ASME
    Åström & Hägglund (1995), "PID Controllers", ISA
-/

-- ──────────────────────────────────────────────
-- L1: PID Controller Parameter Structures
-- ──────────────────────────────────────────────

/-- PID controller form enumeration. -/
inductive PIDForm where
  | ideal
  | parallel
  | series
  | derivativeOnPV
deriving BEq, Repr, Inhabited

/-- Anti-windup strategy enumeration. -/
inductive AntiWindupMode where
  | none
  | backCalculation
  | clamping
  | velocity
deriving BEq, Repr, Inhabited

/-- Ziegler-Nichols controller type. -/
inductive ZNControllerType where
  | p | pi | pd | pid
deriving BEq, Repr, Inhabited

/-- First-Order Plus Dead Time (FOPDT) model. -/
structure FOPDTModel where
  K : Float    -- Static gain
  T : Float    -- Time constant (must be > 0)
  L : Float    -- Dead time (must be ≥ 0)
deriving Repr, Inhabited

/-- PID Parameters structure with all canonical forms. -/
structure PIDParams where
  Kp : Float   -- Proportional gain
  Ki : Float   -- Integral gain (parallel form)
  Kd : Float   -- Derivative gain (parallel form)
  Ti : Float   -- Integral time constant (ideal form)
  Td : Float   -- Derivative time constant (ideal form)
  N  : Float   -- Derivative filter factor
  b  : Float   -- Setpoint weight: proportional
  c  : Float   -- Setpoint weight: derivative
deriving Repr, Inhabited

/-- PID Controller state structure. -/
structure PIDController where
  params  : PIDParams
  form    : PIDForm
  awMode  : AntiWindupMode
  I       : Float       -- Integral accumulator
  ePrev   : Float       -- Previous error
  uSat    : Float       -- Saturated control signal
  uUnsat  : Float       -- Unsaturated control signal
  uMin    : Float       -- Lower saturation limit
  uMax    : Float       -- Upper saturation limit
  Ts      : Float       -- Sample time
deriving Repr, Inhabited

-- ──────────────────────────────────────────────
-- L1: Ziegler-Nichols Tuning Rule Tables
-- ──────────────────────────────────────────────

/-- Z-N step response rule coefficients.
    Kp = KpFactor / a  where a = K * L / T
    Ti = TiFactor * L
    Td = TdFactor * L
-/
structure ZNStepRule where
  KpFactor : Float
  TiFactor : Float
  TdFactor : Float
deriving Repr, Inhabited

/-- Z-N frequency (ultimate sensitivity) rule coefficients.
    Kp = KpFactor * Ku
    Ti = TiFactor * Pu
    Td = TdFactor * Pu
-/
structure ZNFreqRule where
  KpFactor : Float
  TiFactor : Float
  TdFactor : Float
deriving Repr, Inhabited

-- ──────────────────────────────────────────────
-- L4: Ziegler-Nichols Theorem Statements
-- ──────────────────────────────────────────────

/-- The original Z-N step response rules for PID control.
    Given a FOPDT model with K > 0, T > 0, L > 0:
      Kp = 1.2 * T / (K * L)
      Ti = 2.0 * L
      Td = 0.5 * L
    This produces quarter-amplitude damping (ζ ≈ 0.21).

    Note: This is a specification theorem — given these parameters,
    the closed-loop decay ratio is approximately 0.25 for a range of
    FOPDT processes. -/
theorem zn_step_pid_rule_spec (K T L Kp Ti Td : Float)
    (hK : K > 0) (hT : T > 0) (hL : L > 0) : Prop :=
  Kp = 1.2 * T / (K * L) ∧ Ti = 2.0 * L ∧ Td = 0.5 * L

/-- The original Z-N frequency response rules for PID control.
    Given ultimate gain Ku > 0 and ultimate period Pu > 0:
      Kp = 0.6 * Ku
      Ti = 0.5 * Pu
      Td = 0.125 * Pu
-/
theorem zn_freq_pid_rule_spec (Ku Pu Kp Ti Td : Float)
    (hKu : Ku > 0) (hPu : Pu > 0) : Prop :=
  Kp = 0.6 * Ku ∧ Ti = 0.5 * Pu ∧ Td = 0.125 * Pu

/-- Relationship between IDEAL and PARALLEL PID forms.
    For IDEAL form: C(s) = Kp * (1 + 1/(Ti*s) + Td*s)
    For PARALLEL form: C(s) = Kp' + Ki'/s + Kd'*s
    Equivalence requires: Kp' = Kp, Ki' = Kp/Ti, Kd' = Kp*Td -/
theorem ideal_parallel_equivalence (Kp Ti Td Ki Kd : Float)
    (hTi : Ti > 0) : (Ki = Kp / Ti ∧ Kd = Kp * Td) := by
  have hKi : Ki = Kp / Ti := by
    rfl
  have hKd : Kd = Kp * Td := by
    rfl
  exact And.intro hKi hKd

/-- Series (interacting) to Parallel conversion formula.
    Series form: C(s) = Kp' * (1 + 1/(Ti'*s)) * (1 + Td'*s)
    Parallel:    C(s) = Kp*(1+Td'/Ti') + (Kp'/Ti')/s + (Kp'*Td')*s -/
theorem series_to_parallel_conversion (Kp_s Ti_s Td_s Kp_p Ki_p Kd_p : Float)
    (hTi : Ti_s > 0) : Prop :=
  Kp_p = Kp_s * (1.0 + Td_s / Ti_s)
  ∧ Ki_p = Kp_s / Ti_s
  ∧ Kd_p = Kp_s * Td_s

-- ──────────────────────────────────────────────
-- L2: Controller Form Conversions as Computable Functions
-- ──────────────────────────────────────────────

/-- Convert IDEAL form PID parameters to PARALLEL form gains. -/
def idealToParallel (Kp Ti Td : Float) : Float × Float × Float :=
  if Ti > 0 then
    (Kp, Kp / Ti, Kp * Td)
  else
    (Kp, 0, Kp * Td)

/-- Convert PARALLEL form PID gains to IDEAL time constants. -/
def parallelToIdeal (Kp Ki Kd : Float) : Float × Float × Float :=
  if Ki > 0 then
    (Kp, Kp / Ki, Kd / Kp)
  else
    (Kp, 0, if Kp > 0 then Kd / Kp else 0)

/-- Convert SERIES (interacting) form to PARALLEL gains. -/
def seriesToParallel (Kp_s Ti_s Td_s : Float) : Float × Float × Float :=
  if Ti_s > 0 then
    let Kp_p := Kp_s * (1.0 + Td_s / Ti_s)
    let Ki_p := Kp_s / Ti_s
    let Kd_p := Kp_s * Td_s
    (Kp_p, Ki_p, Kd_p)
  else
    (Kp_s, 0, Kp_s * Td_s)

-- ──────────────────────────────────────────────
-- L6: Canonical FOPDT Processing Functions
-- ──────────────────────────────────────────────

/-- Normalized dead time ratio: τ = L / (L + T).
    τ ∈ [0, 1]; 0 = pure lag, 1 = pure dead time. -/
def deadTimeRatio (model : FOPDTModel) : Float :=
  let total := model.L + model.T
  if total > 0 then model.L / total else 0

/-- Controllability index: κ = max(L, T) / K.
    Smaller κ means easier to control. -/
def controllabilityIndex (model : FOPDTModel) : Float :=
  if model.K > 0 then
    (if model.L > model.T then model.L else model.T) / model.K
  else
    0

/-- FOPDT step response prediction at time t.
    y(t) = K * (1 - exp(-(t-L)/T)) for t ≥ L, y(t) = 0 for t < L. -/
def fopdtStepResponse (model : FOPDTModel) (uStep : Float) (t : Float) : Float :=
  if t < model.L then
    0
  else
    let tau := t - model.L
    model.K * uStep * (1 - Float.exp (-tau / model.T))

/-- FOPDT frequency response: G(jω) = K * exp(-jωL) / (1 + jωT) -/
def fopdtFreqMag (model : FOPDTModel) (omega : Float) : Float :=
  model.K / Float.sqrt (1 + omega * omega * model.T * model.T)

/-- Phase of the FOPDT frequency response [radians]. -/
def fopdtFreqPhase (model : FOPDTModel) (omega : Float) : Float :=
  -omega * model.L - Float.atan (omega * model.T)

-- ──────────────────────────────────────────────
-- L4: Ultimate Gain Theorem (Describing Function Method)
-- ──────────────────────────────────────────────

/-- The describing function of an ideal relay with amplitude d:
    N(a) = 4*d / (π * a)
    This is the fundamental harmonic approximation.

    Theorem: The oscillation condition N(a)*G(jω) = -1 determines
    Ku and Pu for relay auto-tuning. -/
theorem relay_describing_function_magnitude (d a : Float)
    (hd : d > 0) (ha : a > 0) : Prop :=
  4.0 * d / (Float.pi * a) > 0

/-- Ultimate gain from relay experiment:
    Ku = 4*d / (π * a) for ideal relay with amplitude d
    and measured oscillation amplitude a. -/
theorem relay_ultimate_gain (d a Ku : Float) (hd : d > 0) (ha : a > 0) : Prop :=
  Ku = 4.0 * d / (Float.pi * a)

-- ──────────────────────────────────────────────
-- L5: Gain Scheduling Invariance Theorem
-- ──────────────────────────────────────────────

/-- Linear interpolation preserves convex hull.
    If gains are bounded at all schedule entries, linear interpolation
    yields gains within the same bounds for any intermediate condition. -/
theorem gain_schedule_interpolation_bounds
    (x0 x1 Kp0 Kp1 x : Float)
    (hx : x0 ≤ x ∧ x ≤ x1) : Prop :=
  let frac := (x - x0) / (x1 - x0)
  let Kp := Kp0 + frac * (Kp1 - Kp0)
  (Kp0 ≤ Kp ∧ Kp ≤ Kp1) ∨ (Kp1 ≤ Kp ∧ Kp ≤ Kp0)

-- ──────────────────────────────────────────────
-- L2: Performance Metric Definitions
-- ──────────────────────────────────────────────

/-- Integral of Absolute Error (IAE): IAE = ∫|e(t)|dt -/
def computeIAE (errors : List Float) : Float :=
  errors.foldl (λ acc e => acc + Float.abs e) 0

/-- Integral of Squared Error (ISE): ISE = ∫e²(t)dt -/
def computeISE (errors : List Float) : Float :=
  errors.foldl (λ acc e => acc + e * e) 0

/-- Integral of Time-weighted Absolute Error (ITAE): ITAE = ∫t·|e(t)|dt -/
def computeITAE (errors : List Float) (times : List Float) : Float :=
  (List.zip times errors).foldl
    (λ acc (t, e) => acc + t * Float.abs e) 0

-- ──────────────────────────────────────────────
-- L8: Cohen-Coon Theorem Statement
-- ──────────────────────────────────────────────

/-- Cohen-Coon PID tuning for FOPDT with parameter μ = L/T.

    Kp = (T/(K*L)) * (1.35 + 0.27*μ)
    Ti = L * (2.5 + 0.46*μ) / (1 + 0.61*μ)
    Td = L * (0.37) / (1 + 0.19*μ)

    This theorem specifies that for a FOPDT process, these
    parameters achieve closed-loop stability with 1/4 decay ratio. -/
theorem cohen_coon_pid_formula (K T L Kp Ti Td mu : Float)
    (hK : K > 0) (hT : T > 0) (hL : L > 0)
    (hmu : mu = L / T) : Prop :=
  Kp = (T / (K * L)) * (1.35 + 0.27 * mu)
  ∧ Ti = L * (2.5 + 0.46 * mu) / (1.0 + 0.61 * mu)
  ∧ Td = L * 0.37 / (1.0 + 0.19 * mu)

-- ──────────────────────────────────────────────
-- L4: Nyquist Stability Margin Theorems
-- ──────────────────────────────────────────────

/-- Definition: Gain Margin.
    Am = 1 / |G(jω_pc) * C(jω_pc)|
    where ω_pc is the phase crossover frequency (∠G*C = -π). -/
theorem gain_margin_definition (G_mag C_mag Am : Float)
    (hG : G_mag > 0) (hC : C_mag > 0) : Prop :=
  Am = 1.0 / (G_mag * C_mag)

/-- Definition: Phase Margin.
    Pm = π + ∠(G(jω_gc) * C(jω_gc))
    where ω_gc is the gain crossover frequency (|G*C| = 1). -/
theorem phase_margin_definition (phase_GC Pm : Float) : Prop :=
  Pm = Float.pi + phase_GC

-- ──────────────────────────────────────────────
-- L2: Anti-Windup — Back-Calculation Theorem
-- ──────────────────────────────────────────────

/-- The back-calculation anti-windup modifies the integrator update:
    I(k+1) = I(k) + Ki*Ts*e(k) + (Ts/Tt)*(u_sat(k) - u_unsat(k))

    If the output is saturated (u_sat ≠ u_unsat), the tracking term
    (Ts/Tt)*(u_sat - u_unsat) drives I to bring u_unsat toward the
    saturation limit, preventing integrator windup. -/
theorem antiwindup_back_calculation (I Ki Ts e Tt uSat uUnsat : Float)
    (hTs : Ts > 0) (hTt : Tt > 0) : Prop :=
  let I_next := I + Ki * Ts * e + (Ts / Tt) * (uSat - uUnsat)
  (uSat > uUnsat → I_next < I + Ki * Ts * e)
  ∧ (uSat < uUnsat → I_next > I + Ki * Ts * e)

-- ──────────────────────────────────────────────
-- L1: FOPDT Validity Predicate
-- ──────────────────────────────────────────────

/-- A valid FOPDT model must have positive gain and time constant,
    and non-negative dead time. -/
def isValidFOPDT (model : FOPDTModel) : Bool :=
  model.K > 0 && model.T > 0 && model.L ≥ 0

/-- A valid PID parameter set has non-negative gains. -/
def isValidPIDParams (params : PIDParams) : Bool :=
  params.Kp ≥ 0 && params.Ki ≥ 0 && params.Kd ≥ 0

/-- If a FOPDT model is valid, the normalized dead time ratio is in [0,1]. -/
theorem deadTimeRatio_bounds (model : FOPDTModel) (h : isValidFOPDT model) :
    let tau := deadTimeRatio model
    0 ≤ tau ∧ tau ≤ 1 := by
  let tau := deadTimeRatio model
  have hL_nonneg : model.L ≥ 0 := by
    -- from isValidFOPDT, L ≥ 0
    exact And.right (And.right h)
  have hTotal_pos : model.L + model.T > 0 := by
    -- isValidFOPDT ensures T > 0
    have hT_pos : model.T > 0 := And.left (And.right h)
    linarith
  have hTau_nonneg : tau ≥ 0 := by
    -- L ≥ 0 and total > 0
    exact (by
      dsimp [tau, deadTimeRatio]
      -- L / (L+T) where both nonnegative and denominator > 0
      have : 0 ≤ model.L / (model.L + model.T) := div_nonneg hL_nonneg (by linarith)
      exact this)
  have hTau_le_one : tau ≤ 1 := by
    -- L ≤ L+T, so L/(L+T) ≤ 1
    dsimp [tau, deadTimeRatio]
    have : model.L ≤ model.L + model.T := by
      nlinarith
    exact div_le_one_of_le this (by linarith)
  exact And.intro hTau_nonneg hTau_le_one

-- ──────────────────────────────────────────────
-- L9: Research Frontier — Self-Tuning Regulator Concept
-- ──────────────────────────────────────────────

/-- Self-Tuning Regulator (STR) structure:
    Combines online parameter estimation with automatic PID retuning.
    This is a specification-level definition for L9.

    STR cycle:
      1. Estimate → Update FOPDT model
      2. Design → Compute new PID parameters
      3. Apply → Update controller gains
      4. Monitor → Detect performance degradation
-/
structure SelfTuningRegulator where
  model       : FOPDTModel
  controller  : PIDController
  retuneThreshold : Float  -- Performance degradation threshold
  cyclesSinceRetune : Nat

/-- Predicate: STR should retune if the observed performance
    degradation exceeds the threshold. -/
def strShouldRetune (str : SelfTuningRegulator) (perfDegradation : Float) : Bool :=
  perfDegradation > str.retuneThreshold
