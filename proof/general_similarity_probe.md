# General contraction-similarity probe (2026-07-21)

## Question

L21 identifies the least similarity square

\[
t_*(T)=\min\{t:I\preceq P\preceq tI,\ T^*PT\preceq P\}.
\]

For `T=phi(A)`, a uniform `t_*(T)<=4` would give a condition-two contraction
similarity and hence a completely bounded Crouzeix theorem.  This may be
strictly stronger than the scalar conjecture, so it needs a general-matrix
falsification gate before more slice-specific certification.

## Protocol

`experiments/general_similarity_sdp.py` generated 120 centered and
numerical-range-normalized complex matrices (`n=3,...,8`, five samples from
each of four families: dense, perturbed Crabb, upper triangular, and nearly
normal; seed `20260721`).  It used the outer domain
`Omega=W(A)+0.01 D`, constructed its Riemann map by Theodorsen iteration, and
evaluated `phi(A)` by Cauchy functional calculus.

A record was accepted only after:

1. adaptive boundary resolutions from 512 through at most 4096 agreed after
   phase alignment;
2. scalar and matrix Cauchy unitality, double-layer positivity/mass, spectral
   radius, and eigenvalue-interpolation cross-checks passed;
3. the primal SDP and L21 dual trace ratio agreed within `2e-3`, and the
   contraction LMI slack was at least `-2e-7`.

These gates test numerical reliability; they do not make the output a proof.

## Results

- 102 of 120 records passed every gate.  None had `t_*>4`.
- The largest accepted value was `3.839929574938` for the `n=3`, near-Crabb,
  sample-zero case.  The accepted median was `2.902898996977`.
- Accepted/max by family: dense `30/30`, `3.400009239257`; near-Crabb
  `29/30`, `3.839929574938`; triangular `30/30`, `3.717197170825`; nearly
  normal `13/30`, `3.218205404227`.
- Seventeen nearly normal cases were rejected because the conformal boundary
  calculation did not converge through the certificate gates at resolution
  4096.  One near-Crabb case was rejected because its SDP primal/dual gap was
  `0.0665`.  Rejected cases are not counted as passes or failures.

The largest case was then recomputed while shrinking the outer offset:

| offset | accepted `t_*` |
|---:|---:|
| 0.02000 | 3.689563892643 |
| 0.01000 | 3.839929574938 |
| 0.00500 | 3.917984555735 |
| 0.00250 | 3.957755651620 |
| 0.00125 | 3.977830404106 |

Every recalculation passed at boundary resolution 1024 with primal/dual gaps
below `4e-8`.  The sequence approaches four from below rather than crossing
it, exactly the behavior expected near the sharp Crabb configuration.

## Verdict and limitations

The similarity route survives its first reliable test away from the elliptic
slice and is now a credible general attack, not merely a slice artifact.  No
counterexample to the stronger `t_*<=4` target was found.  This is still only
finite numerical evidence: the sample is not exhaustive, the calculations use
strict outer offsets, and the nearly normal/polygonal cases expose a current
map-resolution limitation.  A future apparent value above four must survive
the same map, primal/dual, resolution, and offset checks before it is logged as
evidence against the route.

## Exact margin for the next adversarial gate

For a normalized witness `tr(Z)=1`, put

\[
 D=Z-TZT^*.
\]

Writing \(p=\operatorname{tr}D_+\) and
\(n=\operatorname{tr}D_-\), L21's target \(n\leq4p\) is exactly

\[
 \boxed{3\lVert D\rVert_1+5\operatorname{tr}D\geq0}, \tag{1}
\]

because the left side is (8p-2n).  This is one compact nonsmooth
objective for a joint search over general \(A\) and \(Z\succeq0\).  A negative
floating-point value is not a counterexample: the conformal map, the strict
outer domain, positivity of (Z), and the trace-norm eigenvalue split would
all need independent interval enclosures.

Hartz--McCarthy (arXiv:2606.02922, June 2026) supplies a second exact lens.
For the disk-algebra functional calculus \(\theta_T\) and every scalar
functional \(\beta\), their theorem and Paulsen similarity give

\[
 \sqrt{t_*(T)}=\lVert\theta_T\rVert_{cb}
 \leq\max\{1,\lVert\theta_T+\beta I\rVert_{cb}\}. \tag{2}
\]

Thus a geometrically natural scalar shift with cb norm at most two would
prove the general L21 target.  The theorem does not replace the
operator-valued conjugate-Cauchy correction in the Crouzeix--Palencia map by
a scalar functional, so (2) is a route, not a solution.  The first gate
remains adversarial optimization of (1), especially on the previously
rejected nearly normal/flat-range cases and reducible matrices with small
coupling.  A certified violation would kill only the stronger completely
bounded route; the scalar H-r program would remain live.

Reproduction:

```bash
.venv/bin/python -u experiments/general_similarity_sdp.py \
  --per-family 5 --resolution 512 --max-resolution 4096 --inflate 0.01 \
  --output experiments/general_similarity_sweep_s20260721.jsonl
.venv/bin/python -u experiments/general_similarity_sdp.py \
  --per-family 5 --case 3:near_crabb:0 \
  --inflates 0.02 0.01 0.005 0.0025 0.00125 \
  --resolution 512 --max-resolution 8192 \
  --output experiments/general_similarity_inflate_sensitivity_s20260721.jsonl
```
