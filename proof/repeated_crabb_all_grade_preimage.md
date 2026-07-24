# A universal polynomial preimage for every transfer grade

## 1. Result (L212, 2026-07-24)

Let \(S\) be a spectrally stable finite pure partial isometry with
equal orthogonal defect frames

\[
 I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0,
\]

and put

\[
\begin{aligned}
 Q&=I-VV^*,\\
 B_k&=W^*(S^*)^kV,\\
 {\mathfrak C}(K)&=\sum_{n\ge1}B_nKB_n^*.
\end{aligned}                                      \tag{1}
\]

Use the physical equality metric and operator

\[
 P=2I-VV^*+2WW^*,\qquad T=P^{-1/2}SP^{1/2}.        \tag{2}
\]

For every \(k\ge1\), define the balanced column

\[
\boxed{
\widehat C_k=-\frac72Q\left\{
 S^kWB_k+
 \sum_{j=1}^{k-1}(S^*)^{k-j}VB_k^*B_j
\right\},\qquad C_k=P^{1/2}\widehat C_k.}           \tag{3}
\]

Then \(V^*C_k=0\), and L204's physical endpoint map obeys

\[
\boxed{
 {\cal M}_T(C_k)
 =28\left\{
 {\mathfrak C}(B_k^*B_k)-B_kB_k^*
 \right\}.}                                        \tag{4}
\]

Thus no copy projection, pseudoinverse, or constant-rank hypothesis
is needed at any transfer grade.  The extra terms in (3) are exactly
the lower-grade contamination found in L208; retaining rather than
projecting them gives a universal polynomial preimage.

By linearity, for any real summable weights \(r_k\),

\[
\boxed{
\begin{aligned}
\widehat C(r)&=\sum_{k\ge1}r_k\widehat C_k,\\
{\cal M}_T(P^{1/2}\widehat C(r))
&=28\left\{
{\mathfrak C}\!\left(\sum_{k\ge1}r_kB_k^*B_k\right)
-\sum_{k\ge1}r_kB_kB_k^*
\right\}.
\end{aligned}}                                     \tag{5}
\]

In particular \(r_k=|c|^{2k}\) is the exact Faber/Hardy reflection
weight from L201 and L210.  If the prepared physical base endpoint is

\[
12\sum_{k\ge1}|c|^{2k}B_kB_k^*
-28{\mathfrak C}\!\left(
\sum_{k\ge1}|c|^{2k}B_k^*B_k\right),               \tag{6}
\]

then (5) changes it to the coercive endpoint

\[
-16\sum_{k\ge1}|c|^{2k}B_kB_k^*.                  \tag{7}
\]

L212 proves the correction and range statement (3)--(5).  It does
**not** assume or prove the physical base identity (6).  Deriving
(6), or a one-sided endpoint no larger than it, from the complete
Riemann/metric preparation is now the sole elliptic endpoint debt.

## 2. The all-grade intertwining defect

For an arbitrary Hermitian copy matrix \(Y\), let

\[
 H-SHS^*=WYW^*,\qquad
 A=V^*HV,\qquad R=QHV.                              \tag{8}
\]

L208 proves the exact ordered identity

\[
\boxed{
YB_k-B_kA
=W^*(S^*)^kR+
\sum_{j=1}^{k-1}
B_jR^*(S^*)^{k-j}V.}                              \tag{9}
\]

No copy factors commute in (9).

## 3. Dual proof of the polynomial preimage

Stein adjointness for the transfer channel gives

\[
\operatorname{tr}\{Y{\mathfrak C}(K)\}
=\operatorname{tr}(AK).                            \tag{10}
\]

Apply (10) with \(K=B_k^*B_k\).  Since both resulting traces are real,

\[
\begin{aligned}
\operatorname{tr}Y
\left\{{\mathfrak C}(B_k^*B_k)-B_kB_k^*\right\}
&=-\operatorname{Re}\operatorname{tr}
B_k^*(YB_k-B_kA).                                  \tag{11}
\end{aligned}
\]

Substitute (9).  Its first term is

\[
\operatorname{Re}\operatorname{tr}
B_k^*W^*(S^*)^kR
=\operatorname{Re}\operatorname{tr}
(QS^kWB_k)^*R.                                    \tag{12}
\]

For \(1\le j<k\), cyclicity of the scalar trace and taking the real
part give

\[
\begin{aligned}
&\operatorname{Re}\operatorname{tr}
B_k^*B_jR^*(S^*)^{k-j}V\\
&\qquad=
\operatorname{Re}\operatorname{tr}
\left\{Q(S^*)^{k-j}VB_k^*B_j\right\}^*R.           \tag{13}
\end{aligned}
\]

L207--L208's physical endpoint normalization says that every balanced
perpendicular column \(\widehat C\) satisfies

\[
\operatorname{tr}\{Y{\cal M}_T(P^{1/2}\widehat C)\}
=8\operatorname{Re}\operatorname{tr}(\widehat C^*R). \tag{14}
\]

Equations (11)--(14), with the factor \(-7/2\) in (3), prove equality
of the pairings of both sides of (4) against every Hermitian \(Y\).
This proves (4).  Summation proves (5).

## 4. Analyticity and the Crabb apex

For each fixed finite state dimension, \(S\) is algebraic and its
transfer is rational, so only a convergent finite-dimensional series
is involved in (5).  Formula (3) uses multiplication and adjoint
only.  It is therefore jointly real analytic in the equality
colligation data, with no rank-changing inverse.

At the repeated Crabb apex,

\[
B_k=0\quad(k<L),\qquad B_L\ \hbox{unitary}.
\]

The sole possible terminal term also vanishes:

\[
QS^LWB_L=0,
\]

because \(S^LW\mathbb C^m=V\mathbb C^m\), while every lower
contamination factor is zero.  Hence the complete correction (5)
vanishes at the apex, exactly as the scalar elliptic-axis metric
requires.  L201's invertible \(B_L\) makes the left Gram in (7)
positive definite for every nonzero small \(c\).

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_all_grade_preimage.py \
  --output \
  experiments/repeated_crabb_all_grade_preimage_s70224.jsonl
```

The checker uses unstructured stable partial isometries and
independently gauged heterogeneous delay sums.  It tests (4) through
grade six, verifies that the lower-grade contamination is genuinely
nonzero in the unstructured cases, checks the weighted sum (5), and
checks vanishing of (3) on repeated monomial shifts.  The equations
above, not the floating audit, prove L212.  The tracked dataset
SHA-256 is
`bf1d90206cada4fb7983235bd2399b961b1b7ee0e5f6bb3f6f886112b5a126b1`.
