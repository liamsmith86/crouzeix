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

## Near misses
None above 2. Best R = 1.9747 (converging to Crabb-type equality configuration, R → 2⁻).

## Coverage gaps (honest accounting)
- Matrices with nonsmooth W(A) (flat portions / corners): rejected by certificates in operator
  pipeline; not yet searched. Plan: smoothed support-function domains or near-polygon ψ-families.
- Dimensions ≥ 6, high polynomial degrees: unsearched.
- Ω strictly = W(A) (not a neighborhood): extremal machinery currently uses Ω ⊋ σ(A) smooth; the
  2×2 sandbox (Epoch 3) uses Ω = W(A) ellipse directly.
