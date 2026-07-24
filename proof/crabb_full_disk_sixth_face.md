# The recentered full-disk sixth face (A126, 2026-07-24)

## 1. Status

This note records the discovery and fixed-size audit trail after
L176--L177.  L183--L184 have since proved the sixth-order inequality
in every size; see `proof/crabb_full_disk_sixth_hardy_factor.md`.

The initial finite-scale data suggested that the true circular-normal
gradient was always fourth order after full-disk recentering.  The new
exact series audit disproves that extrapolation: the cubic response
first survives at `p=7`.  The correct next target is therefore a
sixth-order Schur inequality.

The remaining problem is no longer sixth-order positivity.  It is the
classification of the equality variety \({\cal R}_L(z)=0\), the
eighth-order fallback there, and the nonlinear singular blow-up.

## 2. Exact recentered disk path

Write `p=L+1`, let

\[
 z=(z_0,\ldots,z_{L-2}),\qquad
 W=z\wedge J\overline z,
\]

and use L176's interval-pulse correction `B(W)`.  The exact disk path is

\[
 H(s)=\frac12I+sZ(z)+s^2B(W).                       \tag{1}
\]

For L122's canonical rank-one Stein metric, let
\(\lambda_-(s),\lambda_+(s)\) be its two simple generalized endpoint
levels and put

\[
 \delta(s)=\lambda_+(s)-4\lambda_-(s).              \tag{2}
\]

The sign of \(\delta\) is the sign of the condition-square excess over
four.  On the generic exact rays tested through `L=9`, the recurrence
gives

\[
 [s^j]\delta=0\quad(0\le j\le5),\qquad
 [s^6]\delta<0.                                     \tag{3}
\]

There are important higher-order strata.  In `L=3`, the first term is
degree eight.  In the real `L=4` slice
\(z=(a_1,a_2,a_3)\),

\[
\boxed{
[s^6]\delta
=-\frac{128}{9}a_2^2(a_1-a_3)^2(2a_1^2+a_2^2).
}                                                    \tag{4}
\]

Thus the sixth face also vanishes when `a_2=0`, even away from the
phase-palindromic equality cone.  On that terminal-only stratum the
first term is degree eight.  A proof cannot divide globally by a
strict sixth-order form.

## 3. First cubic true-normal response

The exact characteristic-Blaschke/Riemann engine was generalized from
quadratic to arbitrary finite series order and evaluated on the fixed
true circular-normal support modes `3,...,p`.

After (1), the quadratic response vanishes in every tested size, as
required by L177.  In the exact dimensions tested, the cubic response
has the following onset:

\[
\begin{array}{c|c}
L&\text{active cubic true-normal modes}\\ \hline
3,4,5&\varnothing\\
6&\{3\}\\
7&\{3,4\}\\
8&\{3,4,5\}\\
9&\{3,4,5,6\}.
\end{array}                                         \tag{5}
\]

For the first case `L=6`, index \(z\) from zero and put

\[
 C(z,W)=6z_3W_{03}-5z_4W_{02}+2z_3W_{12}.           \tag{6}
\]

The complex mode-three response is

\[
\boxed{G_{6,3}^{(3)}=-\frac{22}{45}C(z,W).}          \tag{7}
\]

Equation (7) was first derived on the full symbolic real slice and
checked on five unrelated Gaussian-rational complex directions.
L180's slower upstream checker now propagates five generic complex
variables through all ten true-normal polarizations: mode three is the
only cubic mode and (7) holds coefficientwise.

L173's exact positive curvature on either mode-three support Riesz row
is

\[
 b_{6,3}=\frac{145}{2592}.
\]

Consequently its completed-square gain is

\[
\boxed{
 {\cal G}_6(z)
 =\frac{|G_{6,3}^{(3)}|^2}{4b_{6,3}}
 =\frac{3872}{3625}|C(z,W)|^2.
}                                                    \tag{8}
\]

## 4. The polynomial inequality exposed in `p=7`

On the real `L=6` slice, exact endpoint propagation gives

\[
[s^6]\delta=-\frac{32}{225}P_6(z),                  \tag{9}
\]

where \(P_6\) is a homogeneous real sextic with 68 terms.  The base
condition deficit is therefore \(64P_6/225\).  Combining (8)--(9), the
sixth-order Schur face is nonpositive exactly when

\[
\boxed{
P_6(z)\ge\frac{1089}{290}|C(z,W)|^2.
}                                                    \tag{10}
\]

A degree-three Gram SDP, used only as a discovery tool, returns the
much stronger sharp-looking bound

\[
\boxed{P_6(z)\ge9|C(z,W)|^2,}                        \tag{11}
\]

with equality when the first three Toeplitz coefficients vanish.
Two thousand random complex samples followed by twenty BFGS searches
also converge to the value `9` from above.

At this stage equation (11) was only numerical evidence.  L184 now
proves it at `p=7` as the sharp anti-diagonal Cauchy--Schwarz case.

L179 proves the weaker inequality (10), which is exactly the one needed
for the Schur face, on the complete real `p=7` slice by an exact
rational rank-seven Gram certificate.  L180 now proves (10) on the
complete **complex** `p=7` slice by an independent rational ten-square
certificate on \(z\otimes(z\wedge J\overline z)\).  L184 subsequently
proves the sharper inequality (11) and the arbitrary-size block theorem.

L181 closes the next complex size `p=8`.  It first proves by exact
sparse polarization that modes three and four are the complete cubic
response, then represents their full Schur residual by fifteen positive
rational squares on 37 Pluecker-tensor coordinates.  Thus the first two
active complex sizes supplied exact precursors to L184.

L182 proves the cubic response in arbitrary size.  With
\(S_t=\sum_{i<t-i}(t-2i)W_{i,t-i}\), every active \(G_{L,k}^{(3)}\)
is the explicit triangular interval-flux sum in
`proof/crabb_full_disk_cubic_response.md`; modes above \(L-3\)
vanish.  L183--L184 then close the base side: the canonical disk Hardy
residual has a cubic skew coefficient \({\cal R}_L\), with
\(D_{6,L}=8\|{\cal R}_L\|^2\).  L182's response rows are disjoint
weighted anti-diagonal projections of this matrix, so
Cauchy--Schwarz proves the flux-only Schur inequality and L173's
positive null lift proves the actual one.

L178 subsequently resolves the equality edge of (11) without needing
the full SOS.  When only the last two coefficients remain, the actual
Schur gain/base ratio is exactly

\[
{6(4k-1)^2\over
 6(4k-1)^2+169k(k-1)(k-2)}<1,\qquad k=L-3.
\]

Thus the stronger constant-nine inequality is sharp there, but the
actual Schur face is strict because L173's null-lift curvature retains
the positive `169` term.

## 5. Exact and nonlinear guards

For the canonical deterministic complex rays, the exact ratios
`cubic Schur gain / sixth-order base deficit` are

\[
\begin{array}{c|c}
p&\text{ratio}\\ \hline
7&1612809/35259650\approx0.04574\\
8&377063523/3272726627\approx0.1152\\
9&11112691589457/61490382896306\approx0.1807\\
10&18866054150719637/79550535749464164\approx0.2372.
\end{array}                                         \tag{12}
\]

Five unrelated exact complex rays at `p=7` have maximum ratio about
`0.277`.  A separate optimized-rank-one nonlinear probe uses ten
directions and scales `0.12,0.09,0.0675`.  All 30 residuals are
positive; the largest finite-scale ratio is `0.3375`.

On the first deterministic ray, the optimized defect differs from the
canonical defect approximately at order six and changes the condition
only around order twelve.  This supports, but does not prove, that the
canonical and optimized sixth jets agree.

## 6. Correct next target

The useful all-size statement was not “the cubic response vanishes.”
It was the block inequality

\[
\boxed{
D_{6,L}(z)
-\frac14\sum_{k\ge3}
 (G_{L,k}^{(3)})^*B_{L,k}^{-1}G_{L,k}^{(3)}
\ge0.
}                                                    \tag{13}
\]

L184 proves (13) and identifies its exact kernel as
\({\cal R}_L(z)=0\).  The remaining algebraic program is:

1. classify the cubic equations \({\cal R}_L(z)=0\);
2. compute and prove the eighth-order face on every component;
3. apply the singular-cone analytic blow-up used in L174.

L178 completes this program on the last-two-coefficient edge and
identifies the flux/null-lift split which the general blocks should
retain.  L180 completes it in the first active size and shows that the
right complex Gram organization has separate global-phase sectors; a
formal Hermitian lift of the real certificate is false.  L181 completes
the second active size with the same grade-block organization and
literal extensions of several L180 factor vectors.  L182 closes their
common response transform and L184 closes the sixth face, so no
further fixed-size sixth-order SOS reconstruction is needed.

This is the current load-bearing route to the nonlinear full-disk
tube.  The elliptic and compact merger remains downstream.

## 7. Regeneration

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_base_jet_exact.py \
  --minimum-length 3 --maximum-length 6 \
  --output experiments/crabb_full_disk_base_jet_exact_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_cubic_normal_exact.py \
  --minimum-length 3 --maximum-length 9 \
  --output experiments/crabb_full_disk_cubic_normal_exact_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_cubic_normal_exact.py \
  --minimum-length 6 --maximum-length 6 --direction-count 5 \
  --output \
  experiments/crabb_full_disk_cubic_normal_L6_five_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_sixth_ratio_adversary.py \
  --seed 70226 --random-count 2000 --start-count 20 \
  --output \
  experiments/crabb_full_disk_sixth_ratio_adversary_s70226.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_real_sos_probe.py \
  --solver CLARABEL \
  --output experiments/crabb_full_disk_real_sos_probe_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_normal_schur_probe.py \
  --minimum-length 6 --maximum-length 6 --resolution 4096 \
  --scales 0.12 0.09 0.0675 --direction-count 10 --seed 70225 \
  --output experiments/crabb_full_disk_normal_L6_s70225.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_cubic_response.py \
  --output \
  experiments/crabb_full_disk_complex_cubic_response_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_sixth_certificate.py \
  --output \
  experiments/crabb_full_disk_complex_sixth_certificate_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_L7_cubic_response.py \
  --workers 8 \
  --output \
  experiments/crabb_full_disk_complex_L7_cubic_response_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_L7_sixth_certificate.py \
  --output \
  experiments/crabb_full_disk_complex_L7_sixth_certificate_s70224.jsonl
```

The generalized series routines also regenerate the pre-existing L122
and L173 datasets byte-for-byte after refactoring.
