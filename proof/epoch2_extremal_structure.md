# Epoch 2 — Structure of true extremal pairs (theory notes)

Sources: BGG+20 survey (arXiv:2006.04901, §1.2, §3, §4 read); Crouzeix 2004 Thm 2.1 (cited);
own derivations marked (*).

## Established facts (literature)
- E1 (Crouzeix 04, Thm 2.1): extremal f₀ for (A,Ω) (Ω simply connected, smooth, ⊇ σ(A)) is
  B_A∘φ, φ:Ω→D conformal, B_A Blaschke of degree ≤ n−1. In particular |f₀| ≡ 1 on ∂Ω.
- E2 (orthogonality, [4] Thm 5.1 via BGG+20): K>1 ⟹ ⟨f₀(A)x₀,x₀⟩ = 0.
- E3 (BGG+20 Thm 4.1, cancellation): f₀ = f₁f₂ (each in ball) ⟹ K²⟨f₁(A)x₀,x₀⟩ = ⟨f₁(A)x₀, f₂(A)*f₂(A)x₀⟩.
- E4 (BGG+20 Thm 4.5): extremal measure μ ≥ 0 probability on ∂Ω: ⟨h(A)x₀,x₀⟩ = ∫h dμ ∀h ∈ A(Ω).
- E5 (compressed shifts): for (S_Θ, D), extremal functions are exactly Blaschke of degree < deg Θ (K=1 case).
- Open (their Q1.2): degree of B_A for Ω = W(A)°.

## Derived (*) — to re-verify carefully before use
- D1: Left-side measure. u₀ := f₀(A)x₀/K is extremal for the adjoint problem; first-order
  stationarity of f ↦ ‖f(A)x₀‖ under multiplicative perturbations f₀(1+tq) yields a probability
  measure μ̃ on ∂Ω with ⟨q(A)u₀,u₀⟩ = ∫q dμ̃ ∀q ∈ A(Ω). [Also follows by symmetry: x₀ extremal
  vector of f̃(Ã) for adjoint configuration.] f₀(A)*u₀ = Kx₀.
- D2: c-symmetry: c = |∫f₀ Φ(f₀) dμ| = |∫f₀ Φ(f₀) dμ̃|. [Via c = |⟨(Φ(f₀)f₀)(A)x₀,x₀⟩| and
  x₀ = f₀(A)*u₀/K, multiplicativity.]
- D3: Energy decomposition: W² = K² + 2Re C + G², where C := ∫g₀f₀dμ (c=|C|), G := ‖g₀(A)*x₀‖,
  g₀ := Φ(f₀). And β = conj(∫g₀ dμ).
  ⟹ L7@ext ⟺ K² + 2ReC + G² − |∫g₀dμ|² ≤ (2 − |C|/2)².
- D4: For unimodular f₀ (E1!): on ∂Ω, f̄₀ = 1/f₀, so g₀ = Φ(f₀) = [1/f₀]_analytic = 1/f₀ − PP,
  PP := sum of principal parts of the meromorphic 1/f₀ at the zeros of f₀ in Ω.
  ⟹ **c = |1 − ∫ f₀·PP dμ|.** (Disk: PP = everything, c = 0 — Okubo–Ando mechanism.)
  ⟹ G = ‖(f₀g₀)(A)*u₀‖/K with f₀g₀ = 1 − f₀·PP.
- D5: Crabb control case: g₀ = 0, G = 0, C = 0, W = K, β = 0.

## Key open determinations (numerics first)
- N1: size/sign of C at true extremals (is ReC ≤ 0 forced? C ≈ 0?).
- N2: is G ≈ 0 at true extremals (indicated: W ≈ K at near-extremal configs)? Exact rate?
- N3: test L7@ext on exact extremal pairs with Ω = ellipse (explicit conformal map,
  Blaschke parameterization — exact extremal machinery for ANY matrix A with W(A) ⊆ Ω).
- N4: where in (A,Ω)-space is K largest for fixed Ω-eccentricity? (K → 2 only as Ω → disk?)
  Quantify: K_max(eccentricity) curve — if K ≤ 2 − ρ(ecc) with ρ > 0 for non-disk, and c ≤ C(ecc),
  a two-regime argument (near-disk: perturbative; far-disk: K bounded away) could close everything.

## THE TWO-REGIME STRATEGY (Epoch 2 master plan)
1. Near-disk regime: Ω close to disk (after conformal normalization). Perturbative: c = O(ε),
   K ≤ 2 − ?·ε or K ≤ 2 + O(ε²)-danger — need sharp first-order analysis of K and c in domain
   perturbation. Target: show d/dε[bound] ≤ 0 at disk (the disk is a strict local max of the
   Crouzeix functional in domain space — is it?).
2. Far-from-disk regime: use a(Ω)-type or elliptic estimates (MMOR Cor 5: explicit K(a,b) for
   ellipses!) to keep K ≤ 2 outright. MMOR: K(a,b) = 1+√(1+(2/π)arctan(...)) ≤ 2 iff
   (2/π)arctan(e²/(2√(1−e²))) ≤ ... K(a,1)→2 as a→1: for NEAR-disk ellipses MMOR ALREADY gives ≤ 2−o(1)?
   Check: K(a,b) ≤ 2 ⟺ a(Ω) ≤ ... a(Ω_ell) = (2/π)arctan(e²/(2√(1−e²))): K ≤ 2 ⟺ a(Ω) ≤ 0?? No:
   1+√(1+a) ≤ 2 ⟺ a ≤ 0. So MMOR gives ≤ 2 ONLY for disks. The gap: their route pays c ≤ a(Ω)
   with no K-feedback. Ours has K-feedback via L6. Combining L6 with c ≤ a(Ω)·(localized) may beat
   both: K² ≤ K√(W²−|β|²) + min(a(Ω), ...) and W-control...
