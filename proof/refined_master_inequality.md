# Refined master inequality (Epoch 1 — original work)

Standing hypotheses (H): H finite-dim; Ω smooth bounded convex, W(M)⁻ ⊂ Ω; γ = holomorphic
functional calculus; Φ = conjugate Cauchy transform (unital, antilinear, contractive);
P(σ) ⪰ 0 the double-layer potential, ∮P ds = 2I; K = ‖γ‖ > 1; (f₀,x₀) extremal pair,
μ the extremal (probability) measure of (f₀,x₀); g₀ := Φ(f₀).
[All ingredients PROVED in literature: S–dV Props 1,3, CP17, CGL18.]

## L3 (exact energy identity) — status: PROVED (elementary, from S–dV proof internals)
K² = ⟨γ_Φ(f₀)x₀, γ(f₀)x₀⟩ − conj(∫ g₀ f₀ dμ).
Proof: ⟨γ(f₀)x₀, γ(f₀)x₀⟩ = ⟨γ(f₀)x₀, γ_Φ(f₀)x₀⟩ − ⟨γ(f₀)x₀, γ(g₀)*x₀⟩ and
⟨γ(f₀)x₀, γ(g₀)*x₀⟩ = ⟨γ(g₀)γ(f₀)x₀, x₀⟩ = ⟨γ(g₀f₀)x₀,x₀⟩ = ∫g₀f₀ dμ. Take adjoint/conj as needed. ∎

## L4 (free scalar subtraction) — status: PROVED
For every ω₀ ∈ C: ⟨γ(f₀)x₀, γ_Φ(f₀)x₀⟩ = ⟨γ(f₀)x₀, (γ_Φ(f₀) − ω₀I)x₀⟩, since ⟨γ(f₀)x₀,x₀⟩ = 0 (S–dV Prop 1, needs K>1). ∎

## Definitions (new quantities)
- ν := ⟨P(·)x₀, x₀⟩ ds — positive measure on ∂Ω, total mass 2 ("x₀-localized double-layer measure").
- Exact relation: ∮ f dν = ∫f dμ + conj(∫Φ(f) dμ) for all f ∈ A(Ω). (Both sides = ⟨γ_Φ(f)x₀,x₀⟩.)
- W := ‖γ_Φ(f₀)x₀‖ ∈ [0, 2].
- β := ⟨γ_Φ(f₀)x₀, x₀⟩ = conj(∫ g₀ dμ)   [since ⟨γ(f₀)x₀,x₀⟩ = 0].
- c := |∫ g₀ f₀ dμ| (the S–dV interaction term).

## L5 (unimodularity pressure) — status: PROVED
W² ≤ 2 ∮ |f₀|² dν. In particular W = 2 forces |f₀| = 1 ν-a.e.
Proof: ‖∮f₀ P x₀ ds‖ = sup_{‖y‖=1} |∮ f₀ ⟨P x₀, y⟩ ds| ≤ sup_y ∮|f₀| ⟨Px₀,x₀⟩^{1/2}⟨Py,y⟩^{1/2} ds
≤ (∮|f₀|²⟨Px₀,x₀⟩ds)^{1/2} (∮⟨Py,y⟩ds)^{1/2} = (∮|f₀|²dν)^{1/2}·√2. ∎
(Cauchy–Schwarz for the PSD-operator-valued measure, twice.)

## L6 (refined master inequality) — status: PROVED
K² ≤ K·√(W² − |β|²) + c.
Proof: L3 + L4 with optimal ω₀ = β: |⟨γ(f₀)x₀, (γ_Φ(f₀)−β)x₀⟩| ≤ K‖(γ_Φ(f₀)−β)x₀‖
= K√(W²−|β|²) (orthogonal decomposition of γ_Φ(f₀)x₀ along x₀), plus |∫g₀f₀dμ| = c. ∎

## L7 (target sufficient condition) — status: UNRESOLVED (the new frontier)
Claim to investigate:  √(W² − |β|²) ≤ 2 − c/2  at every extremal configuration.
- If true ⟹ K ≤ 2 always ⟹ CROUZEIX'S CONJECTURE. [K² ≤ K(2−c/2)+c ⟹ (K−2)(K+c/2) ≤ 0.]
- TIGHTNESS CHECK (Crabb, K=2): Ω=disk ⟹ g₀ = Φ(z^k) = 0 ⟹ β=0, c=0, γ_Φ(f₀)x₀ = γ(f₀)x₀ ⟹ W=2.
  Condition reads 2 ≤ 2: EQUALITY. The candidate inequality is sharp exactly at known extremals. ✓
- Note c ≤ |free κ-shift|: c = |∫(g₀−κ)f₀dμ| ∀κ (∫f₀dμ=0) ⟹ c ≤ dist_{L∞(supp μ)}(g₀, C)·∫|f₀|dμ.
- Note L6 unused slack: Cauchy–Schwarz |⟨γ(f₀)x₀, ·⟩| ≤ K‖·‖ discards angle information; the
  even sharper form keeps ⟨u, (γ_Φ(f₀)−β)x₀⟩ with u = γ(f₀)x₀/K.

## L8 (direct-sum obstruction to operator-norm routes) — status: HEURISTIC (to formalize)
For M = M₁ ⊕ M₂ with M₂ normal, spectrum spread in Ω: L(M) := ∮λ_min(P_M) ds → can be ~0
while extremal data (K, c) of M₁ unchanged. Hence ANY sufficient condition phrased through
global operator quantities like ‖γ_Φ − ω‖ or ∮λ_min(P)ds (CGL/(⋆)) is NOT ⊕-stable and cannot
prove the conjecture. The refined chain (L6/L7) uses only x₀-localized quantities (W, β, c, ν) —
these are ⊕-stable: P_M x₀ = P_{M₁}x₀ for x₀ in the first block. **This is why L7 is the right
shape of target and (⋆) was not.**

## Immediate actions
1. Numerically evaluate (W, β, c) at genuinely near-extremal pairs (extremality certificate:
   |⟨γ(f₀)x₀,x₀⟩| ≈ 0 and K stationary). Test L7. [experiments/refined_test.py]
2. Try to prove L7 from L5 + potential theory: W near 2 forces |f₀|=1 ν-a.e.; understand how
   unimodular f₀ + analytic constraint couples ∫g₀f₀dμ with ∫g₀dμ.
3. Check literature for anything equivalent to L6/L7 (MMOR IMRN 2025, "DLP revisited") — avoid rediscovery.
