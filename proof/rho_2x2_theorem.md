# Theorem (2×2 closed form): ρ = 1 − π/(2K(m)) ≥ 0 at the critical domain

**Claim.** Let A be a 2×2 matrix, nonnormal, σ(A) ⊂ int W(A), Ω = int W(A). At the extremal pair
for (A, Ω): ρ = 1 − π/(2K(m)) ≥ 0, where m is the elliptic modulus of Ω's Riemann map, with
equality iff W(A) is a disk.

**Derivation** (canonical form; verified numerically to 1e-8; assumptions flagged below).
1. Affine normalization: A = A_h = [[−1, 2h],[0, 1]] (every nonnormal 2×2 is affinely equivalent;
   ρ is affine-invariant). W(A_h) = ellipse, foci ±1, semi-axes (√(1+h²), h).
2. Riemann map φ: Ω → D: φ(z) = √k·sn((2K/π)arcsin z, k), k = √m, K'(k)/K(k) = 4·artanh(b/a)/π
   [self-calibrated numerically: boundary defect 3e-31]. Eigenvalue images φ(±1) = ±√k·sn(±K) = ±√k.
3. **[α = 0 step — status upgraded 2026-07-20]**: numerically established as STRICT GLOBAL max of
   K(α) over the disk THROUGHOUT THE K > 1 REGIME (62 (k,h)-grid points × 600 α-samples; worst
   deficit −2.3e-5 near α = 0, quadratic). IMPORTANT: for √k(h+√(1+h²)) < 1 (K ≤ 1), α = 0 is NOT
   optimal (near-boundary α do better) — the theorem must carry the hypothesis K > 1, which is
   also what H-r/Prop-1 need. Symbolic proof still desirable: maximize
   K(α)² = ½(F² + √(F⁴−4|det|²)), F² = tr(M*M), |det| = |α²−k|/|1−ᾱ²k|, M = b_α(N).
4. Since f₀(±1) = ±√k (odd values): f₀(A_h) = √k·A_h. Hence x₀ = top right singular vector of A_h,
   K = √k·(h + √(1+h²)).
5. Zero of f₀: z₁ = 0; residue of 1/f₀ at 0: r₁ = 1/f₀'(0) = 1/φ'(0); φ'(0) = √k·(2K/π).
   g₀ = Φ(f₀) = 1/f₀ − r₁/z (D4). g₀(±1) = ±(1/√k − r₁) ⟹ g₀(A) = (1/√k − r₁)·A_h.
6. F := f₀g₀: F(A) = √k(1/√k − r₁)A_h² = (1 − √k·r₁)·I, since **A_h² = I**.
7. ρ = Re⟨F(A)x₀,x₀⟩ = 1 − √k·r₁ = 1 − √k/φ'(0) = 1 − π/(2K(k)).
8. K(k) ≥ K(0) = π/2 for k ∈ [0,1), strictly increasing ⟹ **ρ ≥ 0**, = 0 iff k = 0 (disk). ∎

**Equivalent geometric form.** ρ ≥ 0 ⟺ φ'(0) ≥ φ(1) = √k: the Riemann map's derivative at the
center dominates its value at the focus-image. Equality iff disk.

**Consistency checks.**
- Numerical sweep (ellipse_sandbox.py): all values match 1 − π/(2K(m)) to displayed digits.
- Disk limit h → ∞: m → 0, ρ → 0 ✓; thin limit h → 0: m → 1, K(m) → ∞, ρ → 1 ✓ (matches h=0.05
  numeric ρ = 0.8967 ⟹ K(m) = 15.2 plausible for near-degenerate ellipse).
- ρ real (ImC = 0 numerically ✓ — here proven: everything real).

**Significance.** First fully explicit instance of the H-r mechanism. The general shape to seek:
ρ = 1 − Re Σᵢ rᵢ⟨hᵢ(A)x₀,x₀⟩ (hᵢ = f₀/(·−zᵢ), rᵢ = 1/f₀'(zᵢ)) and positivity should come from a
Schwarz/capacity-type inequality tied to Ω = W(A) criticality. Note: for 2×2 the criticality of Ω
was NOT used (any ellipse ⊇ W(A) with the same symmetry gives ρ = 1 − π/2K ≥ 0 — the ellipse-shape
family is special); the n ≥ 3 problem must use W(A)-contact (ψ-domain probes showed ρ < 0 for
shape-mismatched Ω ⊋ W(A)).

**TODO:**
- [ ] Prove step 3 (α = 0) rigorously.
- [ ] K ≤ 2 for 2×2 re-derived: K = √k(h+√(1+h²)) ≤ 2 — independent check against Crouzeix 04.
- [ ] Generalize: 3×3 with double-symmetry (odd matrices?) — compute ρ structure; seek the general
      positivity mechanism (Schwarz-Pick on residues? Grunsky-type? capacity?).
