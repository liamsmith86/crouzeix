# The exact delayed face of the boundary-layer metric

## 1. Result (L223, 2026-07-24)

Retain L219's balanced pure partial isometry

\[
\begin{gathered}
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\\
P=2I-E+2F,\qquad
B_j=W^*(S^*)^jV,
\end{gathered}                                      \tag{1}
\]

and, with \(q=c^2\),

\[
P_{\rm bl}(q)
=I-\sum_{j\ge1}\frac{q^j}{1+q^j}(S^*)^jES^j
 \sum_{j\ge1}q^jS^jF(S^*)^j.                       \tag{2}
\]

Suppose

\[
B_1=\cdots=B_{k-1}=0.                               \tag{3}
\]

Let

\[
P_{\rm phys}(q)=P^{1/2}P_{\rm bl}(q)P^{1/2},
\]

and take the Schur complement of
\(P_{\rm phys}(q)-4I\) to the upper endpoint
\(W\mathbb C^m\).  Then

\[
\boxed{
\operatorname {Schur}_W\!
\left(P_{\rm phys}(q)-4I\right)
=-4q^kB_kB_k^*+O(q^{k+1}).}                        \tag{4}
\]

Equivalently, in the ellipse parameter,

\[
\boxed{
[c^{2k}]\operatorname {Schur}_W
\left(P_{\rm phys}(c)-4I\right)
=-4B_kB_k^*.}                                      \tag{5}
\]

The orientation is the left Gram, including when \(B_k\) is singular.
If \(B_k=0\), the face vanishes and the same statement advances to
the next transfer grade.

L223 is an exact all-grade metric theorem.  Combined with L219 it
shows that every delayed flag already has a canonical negative upper
budget of the correct reflected weight.  It does **not** make the
balanced elliptic pullback contractive: L219's Stein slack can be
indefinite.  A repair must preserve enough of (4), or replace it by a
stronger negative endpoint.

## 2. Positive upper gap

Put

\[
E_j=(S^*)^jES^j,\qquad F_j=S^jF(S^*)^j.
\]

L219's left-orbit resolution is

\[
\sum_{j\ge0}F_j=I.
\]

Since

\[
P^{-1}=\frac12I+\frac12E-\frac14F,
\]

direct subtraction from (2) gives the positive upper gap

\[
\boxed{
U(q):=4P^{-1}-P_{\rm bl}(q)
=2E+\sum_{j\ge1}\frac{q^j}{1+q^j}E_j
 \sum_{j\ge1}(1-q^j)F_j.}                          \tag{6}
\]

At \(q=0\),

\[
U(0)=2E+I-F.                                       \tag{7}
\]

Its kernel is exactly \(W\mathbb C^m\), and its restriction to
\(W^\perp\) is positive definite.

## 3. The first delayed compression

For every \(j\ge1\),

\[
W^*F_jW=0,                                         \tag{8}
\]

because the left wandering orbits
\(W,SW,S^2W,\ldots\) are mutually orthogonal.
Also

\[
\begin{aligned}
W^*E_jW
&=W^*(S^*)^jVV^*S^jW\\
&=(V^*S^jW)^*(V^*S^jW)
=B_jB_j^*.                                        \tag{9}
\end{aligned}
\]

Compressing (6) to \(W\) and using (3) therefore gives

\[
\boxed{
W^*U(q)W
=q^kB_kB_k^*+O(q^{k+1}).}                         \tag{10}
\]

The cross row starts at the same order.  Indeed, the constant gap
\(U(0)\) annihilates \(W\), every \(F_j\) annihilates \(W\), and

\[
E_jW=(S^*)^jVB_j^*.
\]

Hence

\[
(I-WW^*)U(q)W=O(q^k).                              \tag{11}
\]

Let \(J\) be an isometry onto \(W^\perp\).  By (7),

\[
J^*U(q)J=J^*U(0)J+O(q)
\]

has an analytic inverse.  Equations (10)--(11) show that the Schur
cross-square is

\[
W^*U(q)J\{J^*U(q)J\}^{-1}J^*U(q)W
=O(q^{2k})=O(q^{k+1}).                              \tag{12}
\]

Thus

\[
\operatorname {Schur}_W U(q)
=q^kB_kB_k^*+O(q^{k+1}).                           \tag{13}
\]

## 4. Return to physical coordinates

The physical upper gap is

\[
4I-P_{\rm phys}(q)=P^{1/2}U(q)P^{1/2}.             \tag{14}
\]

The endpoint is an eigenspace of \(P\):

\[
P^{1/2}W=2W.                                       \tag{15}
\]

Schur complements are covariant under a block-diagonal congruence.
The two endpoint factors in (15) multiply (13) by four:

\[
\operatorname {Schur}_W
\left(4I-P_{\rm phys}(q)\right)
=4q^kB_kB_k^*+O(q^{k+1}).                          \tag{16}
\]

Negating (16) proves (4)--(5).

## 5. Consequence for the Stein repair

At the repeated monomial apex, \(k=L\) and \(B_L\) is unitary, so

\[
\operatorname {Schur}_W
\left(P_{\rm phys}(c)-4I\right)
=-4c^{2L}I+O(c^{2L+2}).                            \tag{17}
\]

The optimal axis metric has the larger leading margin
\(-16c^{2L}I\), but (17) is already strict.  More generally, (5)
gives a negative copy-space budget on the range of every first active
\(B_k\).  L222 says that only reducing delayed separators can obstruct
redistribution of a repaired endpoint.  The remaining quantitative
problem is to construct a Stein repair whose cost in (16) is
strictly smaller than this budget, uniformly through rank changes.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_boundary_metric_face.py \
  --output \
  experiments/repeated_crabb_boundary_metric_face_s70224.jsonl
```

The deterministic audit covers grades one through seven,
unstructured fully delayed colligations, rank-changing heterogeneous
shifts, and repeated monomial apices.  It checks every earlier upper
Schur coefficient and the matrix identity (5).  Equations
(6)--(16), not the floating audit, prove L223.  The tracked dataset
SHA-256 is
`9ad734c74ed6cea07d07c32fb2a5cd62249345e4a8ab041ccb28ceb682d6acfe`.
