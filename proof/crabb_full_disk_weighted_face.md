# The full-disk weighted face and recentered normal audit (L175, 2026-07-24)

## 1. Status

This note records a reproducible **numerical structural conjecture**, not
an all-size proof.

L174 proves a tube over the Toeplitz part of L122's disk chart.  The
remaining circular-range coordinates cannot simply be included among
L174's coercive normals: some are tangent to the full general-H disk
manifold and cancel the entire Toeplitz quartic at leading order.

The new calculations identify that cancellation much more precisely.
They also test the first genuine circular-normal face after the
corresponding full-disk recentering.  No adverse sign was found.

## 2. Full weighted-face conjecture

Put \(p=L+1\), let \(C_p\) be the Crabb matrix, and write

\[
 H_0=\frac12I_L,\qquad
 H(s)=H_0+sZ(h).
\]

Let \({\cal C}_p\) be the positive curvature \(-e_p\) from L65 on
physical real matrix coordinates.  Let

\[
 {\cal T}_p=D_{H_0}X(\operatorname{Herm}_L)
\]

be the tangent image of L122's physical general-H disk chart

\[
 X(H)=2K^{-1/2}\widehat H R K^{-1/2},\qquad
 K=\widehat H+R^*\widehat H R.
\]

Finally, let \(g_2(h)\) be the even quadratic coefficient of the
ambient gradient of L157's characteristic-Blaschke branch, after the
first Riemann correction.  The leading full-strong face is

\[
 {\cal F}_4(h,E)
 =-32{\cal Q}(h)+g_2(h)^TE-E^T{\cal C}_pE.             \tag{1}
\]

The data support the exact identities

\[
\boxed{
\begin{aligned}
g_2(h)&\perp\ker{\cal C}_p,\\
g_2(h)&\in\operatorname{range}({\cal C}_p{\cal T}_p),\\
\frac14g_2(h)^T{\cal C}_p^\dagger g_2(h)
  &=32{\cal Q}(h).
\end{aligned}}                                        \tag{2}
\]

Consequently the maximum of (1) is exactly zero, and a maximizing
class modulo \(\ker{\cal C}_p\) is tangent to the exact disk manifold.
Thus the full-H face should be semidefinite, not strictly negative:

\[
 \boxed{\max_E{\cal F}_4(h,E)=0.}                      \tag{3}
\]

This is the precise geometric form of L160 equation (7).  It does not
contradict L173.  L173 maximizes only over the true circular-normal
quotient and gets a strict gain smaller than \(32{\cal Q}\); the extra
gain in (2) comes from disk-tangent curvature directions.

## 3. Rank identities and numerical evidence

The checker reconstructs the L65 Hessian independently, differentiates
the physical general-H chart using divided differences for
\(K^{-1/2}\), and extracts \(g_2\) by symmetric Richardson
extrapolation.  In every tested size it finds

\[
\begin{aligned}
\operatorname{rank}{\cal C}_p&=p(p-2),\\
\operatorname{rank}{\cal T}_p&=L^2-1=p(p-2),\\
\operatorname{rank}({\cal C}_p{\cal T}_p)&=(p-2)^2.
\end{aligned}                                        \tag{4}
\]

Tests cover \(p=4,5,6,7\), one structured and three seeded generic
complex directions per size.  Across all 16 records:

* the gain ratio in (2) differs from one by at most \(7.2\,10^{-10}\);
* the relative component of \(g_2\) in \(\ker{\cal C}_p\) is at most
  \(1.3\,10^{-9}\);
* the relative residual in
  \(2{\cal C}_p{\cal T}_pb=g_2\) is at most \(1.4\,10^{-9}\).

Two clean runs are byte-identical.  These residuals are much smaller
than the effect and stable across unrelated directions, but they are
still floating-point evidence.

## 4. Recentered exact disk paths

L176 now gives an explicit, all-size Hermitian correction
\(B_2(h)\), linear in \(h\wedge J\overline h\), and proves exactly
that its L65 curvature energy is \(32{\cal Q}(h)\).  The numerical
ambient-response equation in (2) selects this same correction.  This
gives the exact-disk path

\[
 H_h(s)=\frac12I+sZ(h)+s^2B_2(h).                     \tag{5}
\]

Every point of (5) is a disk matrix by L122.  On this path the
optimized rank-one envelope loses its quartic term.  The measured base
deficit is generically

\[
 4-\Gamma(H_h(s))=
 \begin{cases}
   O(s^8),&p=4,\\
   O(s^6),&p=5,6.
 \end{cases}                                         \tag{6}
\]

The extra size-four cancellation is a low-dimensional phenomenon and
must not be extrapolated.

For each point, the second checker computes the envelope's exact
defect-vector gradient, the ambient Riemann derivative, the transported
support-mode normal frame, and the L65 normal Schur gain.  Three
directions per size and four scales \(0.15,0.10,0.075,0.05\) give:

* every residual `base deficit - normal gain` is positive;
* the largest normal/base ratio is \(0.223\), at the coarsest scale;
* at scale \(0.05\), the largest ratios are \(1.2\,10^{-4}\),
  \(1.30\,10^{-2}\), and \(2.42\,10^{-2}\) for \(p=4,5,6\);
* for \(p=5,6\), the base is approximately sixth order, the true-normal
  gradient fourth order, and its Schur gain eighth order.  Hence the
  ratio decays approximately as \(s^2\).

The formal completion in the whole L65 range is again close to the
base deficit, but the gradient has a nonzero kernel component of
relative order \(s\).  That full-range number is only a recentering
diagnostic; it is not a valid Schur bound.  The true-normal projection
above is the relevant favorable test.

## 5. What remains to prove

L176 proves the curvature-energy half of (2).  The load-bearing next
target is the response-matching half.  A useful proof should:

1. derive the complete all-mode formula for \(g_2(h)\), extending
   L173's true-normal projection;
2. verify
   \(2{\cal C}_pD X[B_2(h)]=g_2(h)\) for L176's explicit correction.
   The exact energy then follows from L176.

After that identity, the correct nonlinear strategy is a splitting
over the **full** disk manifold:

* use Berger--Okubo--Ando/L122 for the exact disk base;
* use L65/L173 for the genuine coercive circular normals;
* prove the recentered normal response has the extra ideal order
  indicated by (6);
* then merge L117/L149/L150 and L160/L163 for the elliptic and compact
  variables.

One must not claim a strict quartic after all strong variables are
allowed, infer a theorem from the finite-scale exponents, or assume
that L123 classifies every equality disk matrix near \(C_p\).

## 6. Regeneration

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_full_disk_weighted_face.py \
  --minimum-length 3 --maximum-length 6 --resolution 4096 \
  --scale 0.002 --direction-count 4 --seed 70224 \
  --output experiments/crabb_full_disk_weighted_face_s70224.jsonl

OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_full_disk_normal_schur_probe.py \
  --minimum-length 3 --maximum-length 5 --resolution 4096 \
  --scales 0.15 0.1 0.075 0.05 \
  --direction-count 3 --seed 70224 \
  --output experiments/crabb_full_disk_normal_schur_probe_s70224.jsonl
```

The second script uses the analytic real gradient of the rank-one
condition objective.  Its seeded central finite-difference audits have
relative residual at most \(7.9\,10^{-10}\) in the recorded sizes.
