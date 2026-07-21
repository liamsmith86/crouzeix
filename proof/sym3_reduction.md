# Closed-form reduction: doubly-symmetric 3×3 family (VERIFIED numerically, 2026-07-20)

Family: A = [[0,a,0],[b,0,c],[0,d,0]] real; σ(A) = {0, ±e}, e = √(ab+cd); A³ = e²A.
Right/left null vectors: χ = (c, 0, −b), ϱ = (d, 0, −a); obliquity κ = ‖χ‖‖ϱ‖/|ϱ*χ| ≥ 1.
Critical domain Ω = int W(A): doubly symmetric; Riemann map ψ: D → Ω odd; τ := φ(e) ∈ (0,1).

Extremal (numerics): f₀ = B∘φ, B(w) = (w²−α²)/(1−α²w²), α real (symmetric even Blaschke).

## Verified formulas (match pipeline to all printed digits; test case tri(2, .2, 1.4, .1))
- v := B(τ) = (τ²−α²)/(1−α²τ²); m := α²+v.
- f₀(A) = v·I − m·P₀ (P₀ = χϱ*/(ϱ*χ) rank-one oblique).
- **K = ‖[[−α², −m√(κ²−1)],[0, v]]‖** (2×2 triangular; maximize over α).
  [verified: 1.867263 = pipeline 1.867263]
- Conformal data: z₁ = ψ(α), ψ'(α); B'(α) = 2α/(1−α⁴); r₁ = ψ'(α)/B'(α).
- g₀ = Φ(f₀) even with: g₀(0) = −1/α² + 2r₁/z₁; g₀(e) = 1/v − 2r₁z₁/(e²−z₁²).
- g₀(A) = g₀(e)I + (g₀(0)−g₀(e))P₀ ⟹ G_m = [[g₀(0), (g₀(0)−g₀(e))√(κ²−1)],[0, g₀(e)]].
- **ρ = Re⟨G_m M x₀, x₀⟩**, x₀ = top right singular vector of M = [[−α², −m√(κ²−1)],[0, v]].
  [verified: +0.000262 = pipeline +0.000262]

## Structural observations
- At the extremal, g₀(0) ≈ −7e-4, g₀(e) ≈ +2e-3 — Φ(f₀) NEARLY ANNIHILATES σ(A).
  So ρ ≈ 0 to first order; positivity is second-order in the annihilation defects.
- H-r for this family ⟺ inequality among (α*, τ, z₁, ψ'(α), e, κ) where α* = argmax K and the
  conformal data belong to the CRITICAL domain W(A(a,b,c,d)) (n=3 confocality constraint:
  Kippenhahn cubic-class curve with foci {0,±e} and shape tied to κ).
- 2×2 analog: ρ = 1 − √k·r₁-type with positivity ⟸ K(m) ≥ π/2. Here expect an analogous
  "capacity" inequality; the two equality mechanisms (disk: g const; critical: g(σ)≈0) must be
  unified limits.

## Next
1. α-stationarity dK/dα = 0 explicitly (one real equation) → eliminate α.
2. Express the criticality constraint: for the family, parameterize W(A) (Kippenhahn) and get
   (τ, z₁, ψ'(α)) as functions of (a,b,c,d) — numerically map ρ(family) landscape; locate ρ = 0 set.
3. Attempt: prove g₀(0), g₀(e) → 0 identities/estimates at criticality (annihilation lemma), then
   show the second-order form is PSD. This is the concrete n=3 H-r program.
