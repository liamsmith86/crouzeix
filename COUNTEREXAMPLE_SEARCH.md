# COUNTEREXAMPLE_SEARCH.md

Search spaces, algorithms, seeds, near misses, certification status.

## Certification standard (from goal.txt)
A candidate R(A,p) > 2 counts ONLY after: increased precision, conditioning analysis, interval/exact certification, independent implementation. Floating point slightly above 2 = NOT a counterexample.

## Ratio computation method
- ‖p(A)‖₂: SVD of p(A) (Horner). Cross-check with mpmath at ≥50 digits.
- max_{W(A)}|p|: max on ∂W(A) suffices (max-modulus: p analytic everywhere, W(A) compact convex — max of |p| over W(A) is attained on ∂W(A) since |p| subharmonic; NOTE: if p constant, trivial).
- ∂W(A) via θ-sweep: λ_max(Re(e^{-iθ}A)) with eigvec x(θ) → boundary point x*Ax. Handle crossings by dense θ grid + refinement; flat portions appear as arcs collapsing to segments — max of |p| on segment endpoints + interior grid.
- Overestimate check: W(A) is intersection of half-planes → polygon enclosure P_out ⊇ W(A) from finitely many θ gives max_{P_out}|p| ≥ max_{W(A)}|p| → certified LOWER bound on R needs max over W(A) from OUTER approx (max_{P_out}|p| ≥ true max ⇒ R_reported = ‖p(A)‖/max_{P_out}|p| ≤ true R). For claiming R > 2 use outer polygon (safe direction).

## Runs (R = Crouzeix ratio; violation = R > 2 certified — NONE FOUND)
| Run | Space | Algorithm | Best found | Notes |
|---|---|---|---|---|
| results_n3_d3_s1 | n=3, deg 3, 40 restarts | L-BFGS on R_inner (θ=64), verify θ=1024 | R = 1.9747 | Landscape matches Greenbaum–Overton |
| results_n3_d5_s2 | n=3, deg 5, 40 | same | R = 1.8398 | |
| results_n4_d4_s3 | n=4, deg 4, 30 | same | R = 1.9463 | |
| results_n4_d6_s4 | n=4, deg 6, 30 | same | R = 1.3320 | optimizer stalls (nonsmooth) |
| results_n5_d5_s5 | n=5, deg 5, 20 | same | R = 1.2550 | |
| advL7 runs (n=2,3,4) | L7′ violation search | L-BFGS + certified quadrature | max verified viol < 0 | ALL positive prints were artifacts (pitfalls P1–P4 in APPROACH_LEDGER); strict re-eval killed each |
| scalar_adv (b=.6,.3,.15) | scalar L10 | Nelder-Mead | none certified | superseded by targeted construction |
| odd_symmetric_test | scalar L10, odd-symmetric | targeted (theory-driven) | **L10 violated +0.33** (not a Crouzeix c.e. — scalar lemma only) | mechanism documented in LEMMA_LEDGER |
| epoch2/L12 extremal probes | true extremal pairs, ψ-domains | Blaschke-zero opt | L7@ext/P2 hold everywhere; c ≤ 0.006 | data in RESEARCH_STATE |
| direct n=6 (2026-07-21) | n=6, deg ≤ 6, multistart | L-BFGS on R_inner, outer-verify | R = 1.148 | nothing near 2 |
| direct n=7 (2026-07-21) | n=7, deg ≤ 6, multistart | same | R = 1.465 | nothing near 2 |
| Hr_adversarial n=3 | min ρ over certified extremal pairs | multistart + certificate battery | min ρ = +2e-4 | H-r floor POSITIVE |
| Hr_adversarial n=4 | same | same | min ρ = +1.7e-2 (run2: +7.2e-2) | POSITIVE |
| Hr_adversarial n=5 | same | same | min ρ = +1.6e-2 (run2: +1.8e-2) | POSITIVE (weaker search) |
| Hr_adversarial n=6 s42 | same | same | KILLED 2026-07-21 after 2h04m, ZERO output | FUTILE BY DESIGN (pitfall P7): init acceptance loop is uncapped and the diag < 1e-5 extremality certificate is unreachable at n=6 dense — probe reproduced: old solver 14s/eval diag 1.7e-3 REJECT (~500 rejects burned); even best_extremal 98s diag 1.9e-3 REJECT (its deg-3 K 1.58723 beat old deg-5 1.58701). Redesign: structured n=6 families (zero-diag tridiagonal, graded collapse → block certification) or stronger solver |

## Near misses
None above 2. Best R = 1.9747 (converging to Crabb-type equality configuration, R → 2⁻).

## Structural counterexample (not a Crouzeix violation)
The claimed level-4 parity classification is numerically false: exact elliptic-slice weights
`(.8,2.4,1.1,.3)` have a shifted degree-one Möbius candidate with K=1.569762, above the computed
odd 1.511514 and proved-global even 1.489136 sector values. Its rho is +0.15504, so it supports
H-r; it only refutes the reduction to parity sectors. `experiments/slice_phase_audit.py`.

## Coverage gaps (honest accounting, refreshed 2026-07-21)
- Matrices with nonsmooth W(A) (flat portions / corners): rejected by certificates in operator
  pipeline; not yet searched. Plan: smoothed support-function domains or near-polygon ψ-families.
- Dimensions ≥ 8, polynomial degrees > 6: unsearched (n=6,7 now covered at low degree).
- Complex-eigenvalue sectors of the sym-tower (imaginary pairs, probe case C): treated separately,
  lightly probed only.
- Ω strictly = W(A) (not a neighborhood): exact machinery exists for 2×2 (ellipse_sandbox) and
  the elliptic sym4 slice (slice_exact); general domains still use inflate > 0 (theodorsen) —
  K(inflate→0) convergence verified on the slice (monotone from below, ~1%/0.005 inflate).
