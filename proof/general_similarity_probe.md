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

The three strongest accepted triangular `n=3` samples were also recomputed
on the same offset ladder.  Fourteen of fifteen records passed every gate;
the single rejection was an inaccurate SDP at offset `.005`, not a value
above four.  The largest accepted triangular value was `3.807797994233` at
offset `.00125`.  Thus the best non-Crabb family in the initial sweep also
remained below four as the outer domain shrank.

## Repeated-equality-block gate (2026-07-22)

The first general sweep did not directly attack multiplicity at the exact
equality locus.  `experiments/general_similarity_equality_probe.py` now tests
Crabb blocks of sizes three and four with multiplicities two and three.  It
uses three independently seeded transverse perturbation types: full dense,
cross-copy coupling, and noncommuting operator weights/diagonal blocks.  The
perturbation ladder is `1e-4,3e-4,1e-3,3e-3,1e-2`, at offsets zero and
`.000625`.  In addition to every earlier gate, a 4096-angle minimum gap between
the two top support eigenvalues checks that a perturbation has left the
support-multiplicity locus.

The main sweep has 128 records: eight exact repeated-block calibrations and
120 perturbed cases.  All calibrations pass the gates.  At zero offset their
primal values are `4.0000000008` through `4.0000000924`, while their dual
ratios remain below four; this known exact-equality overshoot measures the SDP
noise floor and is not a violation.  Of the perturbed records, 118 pass every
gate and none exceeds four:

| perturbation | accepted | largest `t_*` |
|---|---:|---:|
| full dense | 40/40 | 3.999751968486 |
| cross-copy | 40/40 | 3.999844312562 |
| operator weight | 38/40 | 3.999799090562 |

The two rejected records are one operator-weight matrix evaluated at both
offsets; its primal/dual gaps are too large, and neither record is counted.
Every perturbed matrix has a positive sampled support gap.  The strongest case
is the size-three, multiplicity-two cross coupling at `delta=1e-4`, zero
offset.

A second 48-record run recomputes every `delta=1e-4` direction with boundary
resolution 512→1024, support resolution 16384, and four offsets.  Forty-six
records pass; no accepted value exceeds four.  The strongest case agrees with
the main run within `5e-9`.  Finally, a 65536-angle, boundary-resolution
1024→2048 run follows the two strongest cross directions down to `delta=1e-5`:

| Crabb block size | `delta` | `t_*` | `(4-t_*)/delta` |
|---:|---:|---:|---:|
| 3 | 1e-5 | 3.999984429891 | 1.557011 |
| 3 | 3e-5 | 3.999953290384 | 1.556987 |
| 3 | 1e-4 | 3.999844312558 | 1.556874 |
| 4 | 1e-5 | 3.999983630407 | 1.636959 |
| 4 | 3e-5 | 3.999950891493 | 1.636950 |
| 4 | 1e-4 | 3.999836310481 | 1.636895 |

Thus these genuinely transverse directions leave the equality value with a
stable first-order *decrease*.  This is strong numerical evidence that the
repeated Crabb locus is a local maximum in the tested directions, not a proof
of local maximality and not an exhaustive general-matrix gate.

### Analytic first-order follow-up

L61 subsequently explains and proves the sign of these slopes.  Degenerate
support perturbation and the Schwarz integral reduce the conformal operator
tangent to finitely many Fourier coefficients.  The active-kernel L21 tangent
SDP then has a closed dual: its value is `-4` times the Jensen gap between the
mean of the maximum eigenvalue of the repeated support compression and the
maximum eigenvalue of its mean.  It is therefore nonpositive for every
perturbation, not just the sampled families.  Under a uniform conformal
expansion, a strict feasible-metric lift gives the corresponding upper Dini
derivative bound for `t_*`.  See `proof/general_similarity_tangent.md`.

This settles only first order.  A full neighbourhood theorem still requires
classifying zero-Jensen-gap directions and treating second order there, as
well as closing conformal regularity at compression-eigenvalue crossings.

L62 now performs the next reduction on the foundational single-block face,
where the Jensen gap vanishes in every direction.  First-order complementarity
and second-order Schur complements give three finite affine block LMIs.  Twelve
`p=3,4` directions all have negative quadratic coefficients, including two
nearly flat structured cases.  Second support perturbation and the Schwarz
integral now give the conformal coefficient in finite analytic form; the
arbitrary-`p` SDP sign remains open.  See
`proof/general_similarity_second_order.md`.

For single `3×3` and `4×4` Crabb blocks, L63–L64 subsequently eliminate that SDP
exactly.  The `p=3` value is `-2(Re(E01-E12))^2-21|E20|^2/4`; the `p=4`
value is the five-term negative sum of squares (28).  Hence the stronger
similarity square has nonpositive upper second-order change in every direction
at both low-order blocks.  Equality still leaves large spaces requiring higher
order or an exact-orbit argument.

L65 subsequently closes the same sign in every single-block dimension.  The
support defect is a path Laplacian, and circle grading turns every nonzero mode
into an explicit negative kernel; the grade-zero mode is an exact weighted
shift.  The resulting quadratic form has rank `p(p-2)`.  This remains a local
second-order theorem, not a neighbourhood or general-matrix theorem.

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

If \(\lVert\theta_T\rVert_{cb}>1\), (2) and the choice \(\beta=0\) in fact
give the exact identity

\[
 \inf_\beta\lVert\theta_T+\beta I\rVert_{cb}
 =\lVert\theta_T\rVert_{cb}. \tag{3}
\]

Consequently, finding a scalar shift of cb norm at most two is equivalent to
the L21 target, not a relaxation of it.  Crouzeix--Palencia does not supply
such a shift in general: after pullback it controls

\[
 \theta(f)+R(f),\qquad R(f)=\theta(\alpha(f))^*, \tag{4}
\]

where \(R\) is linear but operator-valued.  At matrix level its entries are
full operator blocks, while Hartz--McCarthy allows only
\([\beta(f_{ij})]\otimes I\).  The two coincide when
\(\theta\mathbin\circ\alpha\) has scalar range, including the disk case, but
not for a general numerical-range domain.  Nor can a unital complete
contraction postprocess (4) while fixing \(\theta\): such a map is completely
positive and star-preserving, hence it fixes \(R\) as well.

A reproducible numerical probe also rules out the most natural attempted
collapse of (4).  For every state \(Q\succeq0\), \(\operatorname{tr}Q=1\), set
\(\beta_Q(f)=\operatorname{tr}(QR(f))\).  If
\((\theta+\beta_QI)/2\) were completely contractive, all its
operator-valued Toeplitz moment matrices would be positive semidefinite.  On
the normalized dense \(3\times3\) seed-`20260721` sample with outer offset
`.01`, maximizing the least eigenvalue of the order-three matrix over *all*
states gave:

| boundary resolution | best possible least eigenvalue |
|---:|---:|
| 512 | -0.01552653 |
| 1024 | -0.01553114 |
| 2048 | -0.01553144 |
| 4096 | -0.01553149 |
| 8192 | -0.01553153 |

SCS and Clarabel agree at resolution 4096 to \(2\times10^{-9}\); the conformal
map diagnostics and an independent numerical rebuild are recorded with every
row.  This is stable numerical evidence, not an interval certificate.  It
closes positive-state scalarization (including normalized trace) as a
promising shortcut.  Hartz--McCarthy currently explains exactly why the
scalar-range disk case works, rather than reducing the general problem.

The live general gate remains adversarial optimization of (1), while the
local equality-block gate has moved to zero-Jensen-gap classification and
second order.  A
certified violation would kill only the stronger completely bounded route;
the scalar H-r program would remain live.  A possible CP-based alternative is
to retain the full correction moments in (4) and ask whether L21's trace
inequality follows from their block-Toeplitz positivity.

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
.venv/bin/python -u experiments/general_similarity_sdp.py \
  --per-family 5 \
  --case 3:triangular:0 --case 3:triangular:1 --case 3:triangular:2 \
  --inflates 0.02 0.01 0.005 0.0025 0.00125 \
  --resolution 512 --max-resolution 8192 \
  --output experiments/general_similarity_triangular_sensitivity_s20260721.jsonl
.venv/bin/python -u experiments/general_similarity_scalarization_probe.py \
  --resolutions 512 1024 2048 4096 8192 --solver CLARABEL \
  --output experiments/general_similarity_scalarization_s20260721.jsonl
.venv/bin/python -u experiments/general_similarity_scalarization_probe.py \
  --resolutions 4096 --solver SCS \
  --output experiments/general_similarity_scalarization_cross_solver_s20260721.jsonl
.venv/bin/python -u experiments/general_similarity_equality_probe.py \
  --output experiments/general_similarity_equality_s9173401.jsonl
.venv/bin/python -u experiments/general_similarity_equality_probe.py \
  --deltas .0001 --inflates 0 .00015625 .0003125 .000625 \
  --no-include-base --resolution 512 --max-resolution 4096 \
  --support-resolution 16384 \
  --output experiments/general_similarity_equality_sensitivity_s9173401.jsonl
.venv/bin/python -u experiments/general_similarity_equality_probe.py \
  --block-sizes 3 4 --multiplicities 2 --families cross \
  --deltas .00001 .00003 .0001 --inflates 0 --no-include-base \
  --resolution 1024 --max-resolution 8192 --support-resolution 65536 \
  --output experiments/general_similarity_equality_ultralocal_s9173401.jsonl
```
