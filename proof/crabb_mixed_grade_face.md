# Distinct-grade first-face polarization target (2026-07-23)

## 1. Exact finite finding

Let `Q_L(u;c)=[a^2]\Gamma_L(a,u,c)` be L118's optimized rank-one
amplitude Hessian on the phase-palindromic equality/ellipse stratum.
For two distinct independent grades `1<=k<l<=floor(L/2)`, define the
real polarization

\[
B_{k,l}(c)
=\frac12\{Q_L(e_k+e_l;c)-Q_L(e_k;c)-Q_L(e_l;c)\}.    \tag{1}
\]

Exact formal arithmetic gives

\[
\boxed{
[c^r]B_{k,l}=0\qquad(0\le r\le k+l)
}                                                     \tag{2}
\]

in every audited pair through size seven.  Separate higher-order
runs found the whole recorded jet zero through order 11 in size five
and through order nine in sizes six and seven.  In particular, the
same-parity pair `(k,l)=(1,3)` vanishes; parity alone does not explain
the observation.

This note does **not** label (2) an all-size lemma.  It states the
precise polarization theorem still to prove.

### A stronger falsification target

The finite data point to a substantially stronger identity:

\[
\boxed{B_{k,l}(c)\equiv0\quad(k\ne l).}              \tag{2a}
\]

For size five, the exact `(1,2)` cross jet vanishes through `c^17`.
In size seven all three grade pairs vanish through `c^11`, including
the same-parity pair `(1,3)`.
Independent finite-amplitude calculations at
`c=0.12,0.25,0.4` find the mixed Hessian at optimizer tolerance; the
apparent residual scales quadratically with the amplitude
finite-difference step and hence is quartic contamination, not an
amplitude-Hessian cross term.

Equation (2a) is a conjecture, not a result used below.  It suggests
that the correct proof may diagonalize the complete Hessian rather
than commute a single Newton face.  The likely coordinates are the
DCT-I modes from L117's elliptic Szegő kernel: the Crabb-axis Stein
kernel is Toeplitz-plus-Hankel and exactly DCT diagonal, while the
linearized characteristic/Faber data are discrete cosine grades.
The missing derivation is to carry the *optimized defect* quadratic
through that transform.  Merely diagonalizing the scalar Faber rows
does not establish (2a).

Here “diagonal” should mean **block diagonal by grade**, not
phase-isotropic inside each complex grade at every fixed `c`.
Complex finite-amplitude probes also put real/imaginary cross-grade
terms at quartic-contamination scale, but the real and imaginary
diagonal Hessians themselves differ away from `c=0` (for example,
in length six, grade one, at `c=.4`).  L134 proves their equality only
on the first Newton face, which is exactly what (4) needs.

## 2. Why `k+l` is the decisive coefficient

L125's inverse-map filtration says that a grade-`k` equality
coefficient cannot reach its reflected elliptic normal below `c^k`.
Therefore a bilinear interaction between grades `k` and `l` cannot
appear below `c^(k+l)`.  The candidate Hardy face is diagonal exactly
when its coefficient at that first allowed order also vanishes:

\[
\boxed{B_{k,l}(c)=O(c^{k+l+1}).}                     \tag{3}
\]

If (3) holds and the one-grade diagonal coefficient is `-64`, then
polarization gives

\[
[a^2]\Gamma_L(a,u,c)
=-64\sum_j|u_j|^2c^{2j}
+\text{terms above the first Newton face}.            \tag{4}
\]

L133--L134 prove the complex diagonal coefficient for offset one;
L132 transfers it to divisor grades.  The diagonal coefficient for a
general nondivisor grade remains a separate gate.

## 3. The associated-graded proof target

Three proved results already identify the ingredients of (3).

1. **Scalar filtration (L125).**  In the associated graded ring, the
   direct ellipse map is
   `z/(1+cz^2)`.  A coefficient of circle grade `j` first enters with
   elliptic weight `c^j`.
2. **Faber reflection (L131).**  At that weight, the grade-`j`
   amplitude row is the coordinate vector
   `4u_jc^j e_(L-j)^*`.  Distinct grades are exactly orthogonal by
   ordinary Fourier/coordinate Parseval.
3. **Metric mode splitting (L65/L118).**  At the Crabb point, the
   optimized rank-one Hessian is L65's circle-mode-diagonal quadratic
   form, and its defect-variable Hessian is invertible.

The missing statement is that taking L118's metric Schur complement
commutes with this associated-graded projection:

\[
\operatorname{gr}_{k+l}
\bigl(\operatorname{Schur}_{\rm defect}{\cal Q}\bigr)
=\operatorname{Schur}_{\rm defect}
\bigl(\operatorname{gr}_{k+l}{\cal Q}\bigr).          \tag{5}
\]

Once (5) is proved, L131 makes the right side diagonal and (3)
follows immediately.

The danger is terminal aliasing: lower-order defect propagation runs
along both ends of the finite Crabb chain before the reflected row is
reached.  A valid proof of (5) must show that every mixed path either

* remains in distinct L65 circle modes and is orthogonal; or
* wraps/folds at a terminal boundary and therefore gains strictly
  more than weight `k+l`.

It is not enough to cite Parseval before checking the defect Schur
correction; L131's falsified naive-defect shortcut shows why.

## 4. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_mixed_grade_face.py \
  --output experiments/crabb_mixed_grade_face_s70223.jsonl
```

The checker uses the guarded rational ellipse/Stein engine.  For each
grade pair it reconstructs the optimized Hessians of both basis
directions and their sum, polarizes them exactly, and audits every
coefficient through `c^(k+l)`.  The calculation is a regression for
(3), not a substitute for (5).

Optional `--first-grade`, `--second-grade`, and `--audit-order`
arguments run beyond the first face.  The persisted deep audits are

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_mixed_grade_face.py \
  --minimum-size 5 --maximum-size 5 \
  --first-grade 1 --second-grade 2 --audit-order 18 \
  --output experiments/crabb_mixed_grade_deep_p5_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_mixed_grade_face.py \
  --minimum-size 7 --maximum-size 7 --audit-order 12 \
  --output experiments/crabb_mixed_grade_deep_p7_s70223.jsonl
```
