# COUNTEREXAMPLE_SEARCH.md

Search spaces, algorithms, seeds, near misses, certification status.

## Certification standard (from goal.txt)
A candidate R(A,p) > 2 counts ONLY after: increased precision, conditioning analysis, interval/exact certification, independent implementation. Floating point slightly above 2 = NOT a counterexample.

## Ratio computation method
- ‖p(A)‖₂: SVD of p(A) (Horner). Cross-check with mpmath at ≥50 digits.
- max_{W(A)}|p|: max on ∂W(A) suffices (max-modulus: p analytic everywhere, W(A) compact convex — max of |p| over W(A) is attained on ∂W(A) since |p| subharmonic; NOTE: if p constant, trivial).
- ∂W(A) via θ-sweep: λ_max(Re(e^{-iθ}A)) with eigvec x(θ) → boundary point x*Ax. Handle crossings by dense θ grid + refinement; flat portions appear as arcs collapsing to segments — max of |p| on segment endpoints + interior grid.
- Overestimate check: W(A) is intersection of half-planes → polygon enclosure P_out ⊇ W(A) from finitely many θ gives max_{P_out}|p| ≥ max_{W(A)}|p| → certified LOWER bound on R needs max over W(A) from OUTER approx (max_{P_out}|p| ≥ true max ⇒ R_reported = ‖p(A)‖/max_{P_out}|p| ≤ true R). For claiming R > 2 use outer polygon (safe direction).

## Runs
| Run | Space | Algorithm | Best R found | Notes |
|---|---|---|---|---|
| (pending) | | | | |

## Near misses
(none yet)
