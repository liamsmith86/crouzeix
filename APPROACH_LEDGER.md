# APPROACH_LEDGER.md

Every attempted approach; why it succeeded/failed. Includes known-failed approaches from literature (do not repeat blindly).

## Known dead ends / pitfalls (from literature — do not repeat)
1. **Improving the abstract C–P lemma alone.** Ransford–Schwenninger: sharp at 1+√2 abstractly. Any improvement MUST inject concrete structure of (A, W(A), Cauchy transform).
2. **Invalid degree reduction via Cayley–Hamilton.** Replacing p by its remainder mod char poly changes max over W(A) (can increase it), so "wlog deg p < n" is FALSE for the sup side. Valid only as: ‖p(A)‖ computable from remainder.
3. **Pointwise-to-operator-norm leaps.** Scalar estimates on W(A) do not transfer to ‖·‖ for nonnormal A.
4. **Proving the completely-bounded version as a shortcut.** Possibly strictly harder; a cb-proof failure says nothing about scalar version.
5. **Von Neumann on a disk containing W(A).** Gives constant depending on disk vs W(A) geometry; radius arguments alone max out well above 2 for eccentric shapes.

## Campaign attempts
| ID | Track | Idea | Status | Outcome / failure condition |
|---|---|---|---|---|
| A1 | A | Reconstruct C–P line-by-line; locate constant entry | DONE | Constant enters as c ≤ 1 in K² ≤ 2dK + c (S–dV Thm 5). Full dissection in proof/track_A_crouzeix_palencia.md |
| C1 | C | Randomized + gradient search max R, dims 3–5, poly deg ≤ 6 | DONE (round 1) | Max R_inner ≈ 1.9747 (n=3,d=3). No violation. Optimizer stalls below 2 at nonsmooth points; landscape matches G–O |
| A2 | A | (⋆): c ≤ 2L with L = ∮λ_min(P)ds (CGL slack) | ABANDONED — WRONG SHAPE | Not ⊕-stable: M⊕(normal junk) sends L→0, keeps c. Global operator quantities cannot work. Lesson → localized quantities |
| A3 | A | L6 refined master inequality (vector-localized) + L7′ target | ACTIVE — CENTRAL | L6 proved. L7′ = the campaign's main open lemma; numerics support; sharp at Crabb |
| C2 | C | Adversarial optimization of L7′ violation over (M,p,x₀) | ACTIVE | Requires self-certifying quadrature (see pitfalls below) |

## Methodological pitfalls discovered (Track G)
- P1. Offset-curve construction silently breaks on matrices with nonsmooth W(M) (near-normal, corners/flats): FFT ringing → garbage integrals that still satisfy holomorphic-identity checks. FIX: certify positivity of P(σ) pointwise + convex tangent winding + mass 2I + Cauchy reconstruction; reject otherwise. COVERAGE GAP: such matrices currently unsearchable — needs smoothed support-function domains.
- P2. Background process hygiene: stale processes with old code contaminated logs → phantom "violations". FIX: harness-tracked jobs only; verify claims by re-evaluating saved records with current code.
- P3. Catastrophic cancellation: optimizer drives poly coefficients to ~1e50 with boundary cancellation; consistency certificates fooled (both sides share garbage). FIX: theorem guards — computed ‖p(M)‖ > 1+√2 or ‖γ_Φ(p)x₀‖ > 2 under certified hypotheses PROVES numerics lie (they contradict established theorems) → reject. Plus coefficient-magnitude cap.
- P4. Near-singular layer-potential quadrature (outer→inner contour, gap ε/2): N=1024 insufficient; passed all other certificates while Φ values were wrong (fake viol +0.358 at n=3, killed at N=4096). FIX: unitality certificate Φ(1)=1 through the same kernel + N-vs-2N stability check.
- **Meta-pattern (important):** every relaxation of extremality so far admits sharp counterconfigurations: abstract lemma (R–S, non-unital α), domain-only constants (MMOR, thin quadrilaterals), scalar localization (L10, odd-symmetric mechanism — ours). Conjecture-strength inequalities must engage TRUE extremal structure (global maximality of f₀ + positive extremal measure), not just its first-order shadows.

## Epoch 2 outcome
2×2/ψ-domain exact extremal machinery built (extremal_pullback.py). Empirical laws at true extremals: c small (0 on disk, ~ε²-ish off), G ≈ |β|, ReC ≤ 0 generic, L7@ext holds broadly. L12-strong false.

## Epoch 3 (current): P2 target, audit-adopted (chatgpt/FABLE_RESEARCH_AUDIT.md — verified before adoption)
P1 (proved): K² + ReC ≤ Kq. P2 (target): q ≤ 2 + ReC/2 at extremals ⟹ Crouzeix. See proof/P2_target.md.
Attacks: (i) falsification-first: exact 2×2, Ω = W(A) ellipse (elliptic Riemann map) + ψ-family adversarial ascent on −s_phase = q − 2 − r/2; (ii) if P2 survives: Euler–Lagrange at extremals (Blaschke-zero variations + singular-vector stationarity) → Stinespring leakage bound; (iii) fallback: pre-Cauchy–Schwarz angle form of P1.
Audit adoptions verified by hand: P1/P2 algebra, Stinespring identities (D), D2-exactness, SV24 §6 baseline (ρ = ReC, K²+ρ≤2K — our β-subtraction is the increment). Audit hygiene items fixed (ledgers synced, 2410.10678 corrected, run tables filled).
