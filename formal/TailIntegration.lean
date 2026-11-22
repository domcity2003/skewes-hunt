import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Tactic.Ring

open Filter MeasureTheory Set
open scoped Interval Topology

namespace Skewes

private theorem zero_at_endpoints {f : ℝ → ℂ} {a b : ℝ}
    (hc : Continuous f) (hz : ∀ x ∉ Icc a b, f x = 0) :
    f a = 0 ∧ f b = 0 := by
  have hl : EqOn f (fun _ => 0) (Iio a) :=
    fun x hx => hz x (fun h => (not_le_of_gt hx) h.1)
  have hr : EqOn f (fun _ => 0) (Ioi b) :=
    fun x hx => hz x (fun h => (not_le_of_gt hx) h.2)
  exact ⟨hl.closure hc continuous_const (by simp),
    hr.closure hc continuous_const (by simp)⟩

private theorem derivative_zero_outside {f g : ℝ → ℂ} {a b : ℝ}
    (hd : ∀ x, HasDerivAt f (g x) x) (hz : ∀ x ∉ Icc a b, f x = 0) :
    ∀ x ∉ Icc a b, g x = 0 := by
  intro x hx
  have hneighborhood : (Icc a b)ᶜ ∈ 𝓝 x := isClosed_Icc.isOpen_compl.mem_nhds hx
  have hlocal : f =ᶠ[𝓝 x] (fun _ => (0 : ℂ)) :=
    mem_of_superset hneighborhood fun y hy => hz y hy
  exact (hd x).unique ((hasDerivAt_const x (0 : ℂ)).congr_of_eventuallyEq hlocal)

private theorem integrate_step {a b : ℝ} {z : ℂ} (hz : z ≠ 0)
    {f g : ℝ → ℂ} (hd : ∀ x, HasDerivAt f (g x) x) (hg : Continuous g)
    (ha : f a = 0) (hb : f b = 0) :
    (∫ x in a..b, f x * Complex.exp (z * (x : ℂ))) =
      (-z⁻¹) * (∫ x in a..b, g x * Complex.exp (z * (x : ℂ))) := by
  let e : ℝ → ℂ := fun x => z⁻¹ * Complex.exp (z * (x : ℂ))
  have he : ∀ x : ℝ, HasDerivAt e (Complex.exp (z * (x : ℂ))) x := by
    intro x
    have hcast : HasDerivAt (fun t : ℝ => (t : ℂ)) 1 x := by
      have hlinear : HasFDerivAt (fun t : ℝ => (t : ℂ)) Complex.ofRealCLM x :=
        Complex.ofRealCLM.hasFDerivAt
      simpa using hlinear.hasDerivAt
    have h := ((hcast.const_mul z).cexp).const_mul z⁻¹
    have hfactor : z⁻¹ * (Complex.exp (z * (x : ℂ)) * (z * 1)) =
        Complex.exp (z * (x : ℂ)) := by
      rw [mul_one, mul_comm (Complex.exp (z * (x : ℂ))) z, ← mul_assoc,
        inv_mul_cancel₀ hz, one_mul]
    rwa [hfactor] at h
  have hexp : Continuous (fun x : ℝ => Complex.exp (z * (x : ℂ))) :=
    Complex.continuous_exp.comp (continuous_const.mul Complex.continuous_ofReal)
  have h := intervalIntegral.integral_mul_deriv_eq_deriv_mul
    (fun x _ => hd x) (fun x _ => he x)
    (hg.intervalIntegrable a b) (hexp.intervalIntegrable a b)
  rw [ha, hb, zero_mul, zero_mul, sub_self, zero_sub] at h
  rw [h]
  have hp : (fun x => g x * e x) =
      (fun x => z⁻¹ * (g x * Complex.exp (z * (x : ℂ)))) := by
    funext x
    dsimp [e]
    rw [mul_left_comm]
  rw [hp, intervalIntegral.integral_const_mul, neg_mul]

/-- A continuous derivative tower supported in the integration interval permits
repeated integration by parts against a complex exponential. -/
theorem integral_exp_parts {a b : ℝ} {z : ℂ} (hz : z ≠ 0)
    (n : ℕ) (f : ℕ → ℝ → ℂ)
    (hcontinuous : ∀ j ≤ n, Continuous (f j))
    (hderiv : ∀ j < n, ∀ x, HasDerivAt (f j) (f (j + 1) x) x)
    (hsupport : ∀ x ∉ Icc a b, f 0 x = 0) :
    (∫ x in a..b, f 0 x * Complex.exp (z * (x : ℂ))) =
      (-z⁻¹) ^ n * (∫ x in a..b, f n x * Complex.exp (z * (x : ℂ))) := by
  induction n generalizing f with
  | zero => simp
  | succ n ih =>
    have hpos : 0 < n + 1 := Nat.zero_lt_succ n
    obtain ⟨ha, hb⟩ := zero_at_endpoints (hcontinuous 0 hpos.le) hsupport
    rw [integrate_step hz (hderiv 0 hpos) (hcontinuous 1 (Nat.succ_le_succ (Nat.zero_le n))) ha hb]
    have h := ih (fun j => f (j + 1))
      (fun j hj => hcontinuous (j + 1) (Nat.succ_le_succ hj))
      (fun j hj => hderiv (j + 1) (Nat.succ_lt_succ hj))
      (derivative_zero_outside (hderiv 0 hpos) hsupport)
    change (∫ x in a..b, f 1 x * Complex.exp (z * (x : ℂ))) = _ at h
    rw [h, pow_succ', mul_assoc]

/-- The last derivative and the exponential envelope give the inverse-frequency bound
for one integral. The same estimate applies when the exponent has a nonzero real part. -/
theorem norm_integral_exp_le {a b : ℝ} {z : ℂ} (hz : z ≠ 0)
    (n : ℕ) (f : ℕ → ℝ → ℂ)
    (hcontinuous : ∀ j ≤ n, Continuous (f j))
    (hderiv : ∀ j < n, ∀ x, HasDerivAt (f j) (f (j + 1) x) x)
    (hsupport : ∀ x ∉ Icc a b, f 0 x = 0)
    {M : ℝ} (hM : ∀ x ∈ uIcc a b, ‖f n x * Complex.exp (z * (x : ℂ))‖ ≤ M) :
    ‖∫ x in a..b, f 0 x * Complex.exp (z * (x : ℂ))‖ ≤ ‖z‖⁻¹ ^ n * (M * |b - a|) := by
  rw [integral_exp_parts hz n f hcontinuous hderiv hsupport, norm_mul, norm_pow,
    norm_neg, norm_inv]
  exact mul_le_mul_of_nonneg_left
    (intervalIntegral.norm_integral_le_of_norm_le_const fun x hx => hM x (uIoc_subset_uIcc hx))
    (pow_nonneg (inv_nonneg.mpr (norm_nonneg z)) n)

/-- Including the prefactor in a single zero term gives decay in the imaginary
ordinate to power `n + 1`. -/
theorem norm_zero_term_le {a b : ℝ} {ρ c : ℂ} (hρ : ρ.im ≠ 0)
    (n : ℕ) (f : ℕ → ℝ → ℂ)
    (hcontinuous : ∀ j ≤ n, Continuous (f j))
    (hderiv : ∀ j < n, ∀ x, HasDerivAt (f j) (f (j + 1) x) x)
    (hsupport : ∀ x ∉ Icc a b, f 0 x = 0)
    {M : ℝ} (hM0 : 0 ≤ M)
    (hM : ∀ x ∈ uIcc a b, ‖f n x * Complex.exp ((ρ - ((1 / 2 : ℝ) : ℂ)) * (x : ℂ))‖ ≤ M) :
    ‖(c / ρ) * (∫ x in a..b, f 0 x * Complex.exp ((ρ - ((1 / 2 : ℝ) : ℂ)) * (x : ℂ)))‖ ≤
      ‖c‖ * (M * |b - a|) / |ρ.im| ^ (n + 1) := by
  let z : ℂ := ρ - ((1 / 2 : ℝ) : ℂ)
  have him : z.im = ρ.im := by simp [z]
  have hz : z ≠ 0 := by
    intro h
    exact hρ (by simpa [h] using him.symm)
  have hi : 0 < |ρ.im| := abs_pos.mpr hρ
  have hnorm : |ρ.im| ≤ ‖z‖ := by simpa only [him] using Complex.abs_im_le_norm z
  have hpow : ‖z‖⁻¹ ^ n ≤ |ρ.im|⁻¹ ^ n :=
    pow_le_pow_left₀ (inv_nonneg.mpr (norm_nonneg z)) (inv_anti₀ hi hnorm) n
  have hfactor : ‖c‖ / ‖ρ‖ ≤ ‖c‖ / |ρ.im| := by
    exact div_le_div_of_nonneg_left (norm_nonneg c) hi (Complex.abs_im_le_norm ρ)
  have hwidth : 0 ≤ M * |b - a| := mul_nonneg hM0 (abs_nonneg _)
  change ‖(c / ρ) * (∫ x in a..b, f 0 x * Complex.exp (z * (x : ℂ)))‖ ≤ _
  rw [norm_mul, norm_div]
  calc
    _ ≤ (‖c‖ / |ρ.im|) * (|ρ.im|⁻¹ ^ n * (M * |b - a|)) := by
      exact mul_le_mul hfactor
        ((norm_integral_exp_le hz n f hcontinuous hderiv hsupport hM).trans
          (mul_le_mul_of_nonneg_right hpow hwidth))
        (norm_nonneg _) (div_nonneg (norm_nonneg c) (abs_nonneg _))
    _ = ‖c‖ * (M * |b - a|) / |ρ.im| ^ (n + 1) := by
      simp only [pow_succ, div_eq_mul_inv, mul_inv_rev, inv_pow]
      ring

end Skewes

#print axioms Skewes.integral_exp_parts
#print axioms Skewes.norm_integral_exp_le
#print axioms Skewes.norm_zero_term_le
