# P2 — the campaign's central target (Epoch 3)

Adopted from external audit (proof/initial_strategy_audit.md) after independent verification
(2026-07-20). Supersedes L7/L7@ext as the target; L6 superseded by P1.

## Setting (as in refined_master_inequality.md; conventions: ⟨a,b⟩ linear in first arg)
Extremal pair (f₀,x₀) for (A,Ω), K = ‖γ‖ > 1, μ extremal measure; g₀ = Φ(f₀) = K_Ω-transform;
C := ⟨(g₀f₀)(A)x₀,x₀⟩ = ∫g₀f₀dμ; r := Re C; W := ‖γ_Φ(f₀)x₀‖; β := ⟨γ_Φ(f₀)x₀,x₀⟩;
q := √(W²−|β|²).

## P1 (PROVED — verified by hand 2026-07-20)
K² + r ≤ K·q.
Proof: K² + conj(C) = ⟨(γ_Φ(f₀)−β)x₀, f₀(A)x₀⟩ (exact identity; β-subtraction free since
⟨f₀(A)x₀,x₀⟩=0). Take Re, Cauchy–Schwarz; ‖(γ_Φ(f₀)−β)x₀‖ = q by orthogonal decomposition. ∎
Relation to literature: Schwenninger–de Vries arXiv:2409.15954 §6 (VERIFIED against PDF) proves
K² + ρ ≤ 2K with ρ = r (their Theorem 6.1: K ≤ 1+√(1−ρ)); their proof already uses W ≤ 2.
**P1's new content over SV24 = replacing 2 by q ≤ W (β-subtraction).** They also prove |ρ| ≤ a(Ω) < 1
(Remark 6.2) and ρ ≥ −1 for convex Ω, and disk ⟹ (K−1)ρ = 0.

## P2 (CENTRAL TARGET — UNRESOLVED)
Claim: q ≤ 2 + r/2 at every extremal pair (K > 1, Ω convex smooth ⊇ W(A)).
- P1 + P2 ⟹ (K−2)(K − r/2) ≤ 0 and K − r/2 ≥ K − 1/2 > 0 ⟹ **K ≤ 2 ⟹ CROUZEIX** (verified).
- If r ≥ 0: P2 trivial (q ≤ W ≤ 2). The difficult case is r < 0 — which IS generic at extremals
  per Epoch-2 data (ReC < 0 predominantly off-disk).
- Sharp at Crabb (q = 2, r = 0).
- Weaker than dead L7@ext (which used |C| ≥ |r|); ImC costs nothing here.

## Stinespring defect form (verified identities)
T(h) := ½∮hP ds is positive unital; Stinespring T = V*π(·)V, ξ := Vx₀, U := π(f₀) unitary
(|f₀|≡1 on ∂Ω by Crouzeix's Blaschke theorem), Q := VV*. Then:
  W² = 4‖QUξ‖² = 4(1 − ‖(I−Q)Uξ‖²);  β = 2⟨Uξ,ξ⟩;
  4 − q² = 4[‖(I−Q)Uξ‖² + |⟨Uξ,ξ⟩|²].
**P2 (r<0 case) ⟺ ‖(I−Q)Uξ‖² + |⟨Uξ,ξ⟩|² ≥ |r|/2·(1/2) − r²/16 ... precisely: ≥ (−2r − r²/4)/4.**
Interpretation: negative interaction r forces quantified leakage of Uξ out of ran V or self-overlap.
This is the quantitative L9 and explains the L10 failure (scalar data can't see leakage).

## Epoch-3 plan (audit order, adopted)
1. Re-score all stored extremal records with s_phase := 2 + r/2 − q (P2 slack). [Existing Epoch-2
   data: P2 holds a fortiori since L7@ext held; margins larger.]
2. Falsification-first: exact 2×2 sandbox, Ω = W(A) (ellipse; elliptic-function Riemann map),
   deg-1 Blaschke extremals; adversarial ascent on −s_phase over eccentricity + matrix params.
   Also ψ-family domains at n=2,3 (pullback machinery) targeting r < 0, q → 2.
3. If P2 survives: derive full Euler–Lagrange at extremals (zero-variations of the Blaschke product;
   left/right singular vector stationarity; D1 adjoint measure done rigorously via conjugate algebra)
   → prove the defect bound (Stinespring stability route).
4. If P2 fails at a genuine extremal: fall back to the pre-Cauchy–Schwarz angle form
   K² + conj(C) = ⟨(γ_Φ(f₀)−β)x₀, Ku₀⟩ — the phase of the inner product is still unused slack.

## Notes
- D2 exactness: ⟨Hx₀,x₀⟩ = ⟨Hu₀,u₀⟩ for H = (f₀g₀)(A) (proof: Hf₀(A)x₀ = f₀(A)Hx₀, f₀(A)*u₀ = Kx₀). Verified.
- Caution (audit item): arXiv:2603.15536 (q-numerical ranges) Theorem 3.3's K<2-type sentence may
  have a gap; do not build on it without independent verification.
