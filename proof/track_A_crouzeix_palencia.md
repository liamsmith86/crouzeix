# Track A — The Crouzeix–Palencia mechanism, dissected

Sources verified: Ransford–Schwenninger (SIMAX 2018, arXiv:1708.08633v2, read in full);
Schwenninger–de Vries "On Abstract Spectral Constants" (arXiv:2302.05389, read in full).

## Setup (all verified from the papers)
- Ω ⊂ C smooth bounded open, A(Ω) uniform algebra of functions holomorphic on Ω, continuous on Ω̄.
- M bounded on Hilbert H, σ(M) ⊂ Ω. γ: A(Ω) → L(H) the Riesz–Dunford calculus, γ(f) = f(M). K := ‖γ‖ ≥ 1.
- Φ: A(Ω) → A(Ω) conjugate Cauchy transform Φ(f)(z) = (1/2πi)∮ f(σ)* (σ−z)^{-1} dσ. Antilinear, UNITAL (Φ(1)=1). If Ω convex: Φ contractive [CP17 Lemma 2.1].
- Symmetrized calculus γ_Φ(f) := γ(f) + γ(Φ(f))*. NOT multiplicative.
- Double-layer potential P(σ) := (v(σ)/2π)(σ−M)^{-1} + [same]* (v = i·unit tangent = outward normal rotation). Hermitian-valued.
- Identity: γ_Φ(f) = ∮_{∂Ω} f(σ(s)) P(σ(s)) ds, and ∮ P(σ(s)) ds = 2I.
- If Ω convex and W(M) ⊆ Ω: P(σ) ⪰ 0 ∀σ ∈ ∂Ω ⟹ ‖γ_Φ‖ ≤ 2 (Cauchy–Schwarz for PSD operator measure).

## Extremal machinery (Schwenninger–de Vries)
Assume H finite-dim (suffices: Prop 9 — Krylov compression argument reduces general H to finite-dim, NO dimension constant lost). Then an extremal pair (f₀,x₀) exists: ‖f₀‖_Ω = ‖x₀‖ = 1, ‖γ(f₀)x₀‖ = K.
- **Prop 1**: K > 1 ⟹ ⟨γ(f₀)x₀, x₀⟩ = 0.
- **Prop 3**: ∃ probability measure μ on ∂Ω (extremal measure): ⟨γ(f)x₀,x₀⟩ = ∫ f dμ ∀f ∈ A(Ω). Corollary: ∫ f₀ dμ = 0 when K > 1.
- **Lemma 2**: γ(f₀)*γ(f₀)x₀ = K²x₀ (x₀ is the top right-singular vector of f₀(M)).

## Theorem 5 (the master bound)
For any bounded antilinear Φ, K > 1, extremal (f₀,x₀):
K ≤ d + sqrt(d² + c),  where
d := ½ inf_{ω ∈ A(Ω)'} ‖γ_Φ − ω‖   (ω scalar functional, embedded f ↦ ω(f)·I — rank-one-type perturbation),
c := |⟨γ(Φ(f₀)f₀)x₀, x₀⟩| = |∫ Φ(f₀)·f₀ dμ|.

Proof core: K² = ⟨γ(f₀)x₀, γ(f₀)x₀⟩ = ⟨γ(f₀)x₀, (γ_Φ(f₀) − ω(f₀))x₀⟩ − ⟨γ(Φ(f₀)f₀)x₀,x₀⟩* + ω(f₀)⟨γ(f₀)x₀,x₀⟩ and the last term dies by Prop 1; hence K² ≤ 2dK + c. [Uses multiplicativity of γ: γ(Φ(f₀))γ(f₀) = γ(Φ(f₀)f₀).]

### Where each constant enters
- C–P 2017: d ≤ 1 (take ω = 0; ‖γ_Φ‖ ≤ 2), c ≤ ‖Φ(f₀)‖‖f₀‖ ≤ 1 ⟹ K² ≤ 2K + 1 ⟹ K ≤ 1+√2. **The entire loss from 2 to 1+√2 is the crude bound c ≤ 1.**
- Okubo–Ando (Ω disk): Φ(f₀) is CONSTANT (range of Φ on disk = constants) ⟹ c = |Φ(f₀)| · |∫f₀ dμ| = 0 by Prop 3 ⟹ K ≤ 2. Sharp.
- R–S sharpness example: T=[[1,1],[0,0]], Ω = two tiny disks (NOT ⊇ W(T), not connected). Their α is antilinear but NOT unital (1 ↦ −1). c = 1 attained, K = 1+√2. **The sharpness example violates unitality and the geometric containment — both are real structure still on the table.**

## The quantified frontier
**Remark 13 (S–dV)**: Crouzeix's conjecture FOLLOWS if for every finite-dim M, convex smooth Ω ⊃ W(M)⁻, extremal (f₀,x₀,μ):
  inf_ω ‖γ_Φ − ω‖ ≤ 2 − c/2,  i.e.  2d + c/2 ≤ 2.
[Check: d ≤ 1 − c/4 gives d + √(d²+c) ≤ (1−c/4) + √((1−c/4)²+c) = (1−c/4)+(1+c/4) = 2.]

**Prop 8 (Caldwell–Greenbaum–Li)**: inf_ω ‖γ_Φ − ω‖ ≤ 2 − L, where L := ∮ λ_min(P(σ(s))) ds ≥ 0, via ω(f) := ∮ f(σ(s)) λ_min(P(σ(s))) ds.

**⟹ SUFFICIENT CONDITION FOR CROUZEIX (new working target, "Inequality (⋆)"):**
    c ≤ 2L + (c²/8 rounding? — exact algebra):  K ≤ 2 iff (2−L)/2 + sqrt((2−L)²/4 + c) ≤ 2 iff c ≤ 2L. 
    **(⋆):  |∫ Φ(f₀) f₀ dμ| ≤ 2 ∮ λ_min(P(σ)) ds**  at every extremal configuration with K > 1.

## Sanity checks of (⋆)
- Disk: LHS = 0 ✓ (any L ≥ 0 works).
- If K = 2 exactly (Crabb/Jordan equality, W = disk): LHS = 0, consistent even with L = 0.
- Danger zone: Ω shrinking to non-circular W(M): λ_min(P) → 0 where ∂Ω touches "flat"/extremal parts of W(M)? And c > 0 possible. **This is exactly where the conjecture lives.**

## Notes on known follow-up work to check before attacking (⋆)
- MMOR 2024 (C_N < 1+√2): likely uses exactly this slack with compactness for fixed N. Read to see why not uniform.
- "The Double-Layer Potential for Spectral Constants Revisited" (Schwenninger et al., ~2024).
- MMOR IMRN 2025 "Double-Layer Potentials, Configuration Constants, and Applications to Numerical Ranges".
- Clouâtre–Ostermann–Ransford, "An Abstract Approach to the Crouzeix Conjecture" (arXiv:2011.10422, JOT) — abstract framework [COR20]; also studies limitations.
- R–S Question 4.1 (unitality forces 2?): status to determine — if still open, it's a clean attack target: EITHER prove it (⇒ Crouzeix) OR find abstract counterexample (⇒ this whole route capped; pivot).

## Immediate Track A actions
1. Numerically compute c and L for: Crabb matrices with Ω = enlarged disk (control), non-disk W examples, and the local maximizers from Track C. Test (⋆) empirically. [experiments/star_inequality.py]
2. Determine status of Question 4.1 (search: Bickel–Gorkin–et al. CMFT 2020 survey [BGG+20] states it as open problem; later work?).
3. Understand the range of Φ for non-disk convex Ω (for disk = constants; the failure of c = 0 measures non-circularity — quantify via Neumann–Poincaré spectrum: Φ relates to the NP operator whose spectrum on convex smooth domains lies in [0? ...]; configuration constant literature).
