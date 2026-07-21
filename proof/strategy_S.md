# Strategy S′ — reduction of Crouzeix to a domain-monotonicity lemma (Epoch 3)

## Chain (each link's status marked)
1. [PROVED, SV24 Thm 6.1 — verified] For any smooth convex Ω ⊇ W(A)⁻, extremal pair (f₀,x₀):
   K(A,Ω) ≤ 1 + √(1 − ρ), ρ := Re∫Φ(f₀)f₀ dμ₀.
2. [PROVED, classical] ρ ≥ 0 at Ω would give K(A,Ω) ≤ 2; applying along Ω ↓ int W(A) gives the
   conjecture (shrinking argument; SV24 Prop A.2 handles approximation).
   Boundary-eigenvalue matrices reduce (normal ⊕ rest) — handle separately.
3. [ANCHOR — provable pieces] Disk case: Φ(f₀) constant ⟹ (K−1)ρ = 0 [SV24]; so at an exact disk
   with K > 1: ρ = 0. If K = 1 at the disk end: unimodular-constant extremal pairs give ρ = 1 ≥ 0
   (branch-selection issue to handle carefully).
4. **[NEW LEMMA NEEDED — "M"] Monotonicity: along the Minkowski family
   Ω_s := (1−s)·W(A) + s·D_R (D_R = a disk ⊇ W(A); support functions interpolate linearly:
   h_s = (1−s)h_{W(A)} + s h_{D_R}), the extremal value ρ(s) satisfies ρ(0) ≥ ρ(1).**
   (Weaker than pointwise dρ/ds ≤ 0 — an endpoint comparison suffices!)
5. 1+3+4 ⟹ ρ(0) ≥ 0 at Ω = W(A)-limit ⟹ K ≤ 2 ⟹ CROUZEIX.

## Evidence for M (2026-07-20)
- Offset families (= Minkowski with disk, reparameterized): ρ strictly ↑ as domain shrinks, 3/3
  matrices (n=3), all t ∈ [0.002, 0.05]. (experiments/Hr_test.py output.)
- Exact 2×2 complete sweep at Ω = W(A): ρ ∈ (0, 0.9] > 0 always. ρ → (2−K) in round limit
  (saturating SV24 equality asymptotically).
- ψ-domain probes with SHAPE-MISMATCHED Ω ⊋ W(A) gave ρ < 0 — consistent: M is claimed only along
  the Minkowski-to-disk family, not all domain paths. (H-r needs only ONE good path per A.)

## Attack plan for M
- (a) Numerics: support-function-driven Theodorsen (Minkowski families are linear in support
  functions); test ρ(s) monotonicity/endpoint comparison over battery of A, n = 3..5. FIRST.
- (b) Theory: envelope/Dini derivative of s ↦ ρ(s) through the extremal pair; domain-variation of
  the conjugate Cauchy transform Φ_Ω (kernel varies with ∂Ω_s); DLP degeneracy structure at s = 0
  (kernel field P(σ_θ)(σ_θ−A)x_θ = 0 — derived, to verify). The variation of ρ has an explicit
  boundary-integral form; positivity of the DLP and the contact structure at W(A) should enter.
- (c) Watch for failure mode: extremal-pair branch jumps (K or f₀ discontinuous in s) — use upper
  Dini derivatives / interval argument; ρ need only be compared at endpoints through a chain of
  branch segments.

## Fallbacks if M false
- P2 (q ≤ 2 + ρ/2) remains: strictly weaker demand than ρ ≥ 0 when q < 2. The β-subtraction gives
  q < 2 strictly whenever β ≠ 0 (generic). So even ρ slightly negative is tolerable:
  need ρ ≥ −2(2−q)... precisely q ≤ 2 + ρ/2.
- Angle form of P1 (unused phase of the inner product).
