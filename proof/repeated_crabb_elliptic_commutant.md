# The commutant theorem for the repeated elliptic cokernel

## 1. Result (L206, 2026-07-24)

At every fixed L193 inverse-block-Toeplitz equality anchor, L203's
oriented second face is solvable:

\[
\boxed{
 X-T^*XT=F_2+VC_2^*+C_2V^*,\qquad
 V^*XV=0,\qquad W^*XW=-ZZ^*.}                      \tag{1}
\]

Here

\[
 Z=W^*P^{1/2}E_1P^{-1/2}V=4B_1.                   \tag{2}
\]

Thus the grade-one upper loss is the correctly oriented left Gram

\[
 \boxed{-ZZ^*=-16B_1B_1^*.}                       \tag{3}
\]

The proof closes every cokernel condition in L204.  It is
finite-dimensional and order-safe.  It does **not** yet produce a
bounded real-analytic choice of \(X,C_2\) through the rank-jumping
Crabb apex, and it does not cover higher reflected grades.

The new structural statement is that L204's cokernel is exactly the
self-adjoint commutant of L201's partial-isometry colligation.

## 2. Balanced observability Gramian

Retain L202's notation

\[
 S=P^{1/2}TP^{-1/2},\qquad
 I-S^*S=VV^*,\qquad I-SS^*=WW^*.                  \tag{4}
\]

For \(Y=Y^*\in M_m\), let the physical dual Gramian be

\[
 {\cal Z}_Y-T{\cal Z}_YT^*=WYW^*.                  \tag{5}
\]

Normalize it in balanced coordinates:

\[
 H_Y=\frac14P^{1/2}{\cal Z}_YP^{1/2}.              \tag{6}
\]

Since \(P^{1/2}W=2W\), equations (5)--(6) give

\[
\boxed{
 H_Y-SH_YS^*=WYW^*,\qquad
 H_Y=\sum_{n\ge0}S^nWYW^*(S^*)^n.}                \tag{7}
\]

Let

\[
 Q=I-VV^*=S^*S.
\]

Because \(P\) commutes with \(VV^*\), L204's adjoint-kernel condition

\[
 (I-VV^*){\cal Z}_YV=0
\]

is equivalent to

\[
 \boxed{QH_YV=0.}                                  \tag{8}
\]

Hence there is a Hermitian \(A\in M_m\) with

\[
 H_YV=VA,\qquad A=V^*H_YV.                         \tag{9}
\]

## 3. Cokernel equals commutant

Multiply (7) on the right by \(S\).  Since \(W^*S=0\),

\[
 H_YS=SH_YS^*S=SH_YQ.
\]

Using \(Q=I-VV^*\) and (9),

\[
 H_YS-SH_Y=-SH_YVV^*=-SVAV^*=0.                   \tag{10}
\]

Thus

\[
 \boxed{H_YS=SH_Y.}                                \tag{11}
\]

Conversely, if the Hermitian Gramian \(H_Y\) commutes with \(S\), then

\[
 SH_YV=H_YSV=0,
\]

so \(H_YV\in\ker S=V\mathbb C^m\), which is (8).  Therefore (8) and
(11) are equivalent.

Since \(H_Y\) is Hermitian, (11) also gives commutation with \(S^*\).
The left defect is then invariant:

\[
 S^*H_YW=H_YS^*W=0.
\]

The \(n=0\) term in (7), together with \(W^*S^n=0\) for \(n\ge1\),
gives \(W^*H_YW=Y\).  Hence

\[
\boxed{
 H_YW=WY,\qquad H_YV=VA.}                          \tag{12}
\]

This proves the state-space commutant description.

## 4. Equivalent transfer intertwining

Let

\[
 B_n=W^*(S^*)^nV,\qquad
 B_H(z)=\sum_{n\ge1}B_nz^n                        \tag{13}
\]

be L201's genuine matrix-inner transfer.  Equations (11)--(12) give

\[
\boxed{
 YB_n=B_nA\quad(n\ge1),\qquad
 YB_H(z)=B_H(z)A.}                                 \tag{14}
\]

There is also a converse.  Purity of \(S\) and
\(I-SS^*=WW^*\) give the telescoping observability identity

\[
 \sum_{n\ge0}S^nWW^*(S^*)^n=I.                    \tag{15}
\]

If (14) holds, then

\[
\begin{aligned}
 H_YV
 &=\sum_{n\ge0}S^nWY\,W^*(S^*)^nV\\
 &=\sum_{n\ge0}S^nWYB_n\\
 &=\sum_{n\ge0}S^nWB_nA
 =VA,
\end{aligned}                                      \tag{16}
\]

where the last equality is (15) applied to \(V\).  Thus

\[
\boxed{
\ker{\cal M}_T^*
\cong
\{(Y,A):Y=Y^*,A=A^*,\ YB_H=B_HA\}.}                \tag{17}
\]

At the Crabb apex \(B_H(z)=Uz^L\), every Hermitian \(Y\) occurs, with
\(A=U^*YU\).  At an irreducible nearby colligation, only scalar
intertwiners remain.  Reducible copy strata retain one scalar on each
reducing summand.  This explains exactly the rank jump observed in
L203--L204.

## 5. Weighted L203 identity

It remains to prove L204's Fredholm condition

\[
 \operatorname{tr}(YD_T)=0
\quad(Y\in\ker{\cal M}_T^*),                        \tag{18}
\]

where

\[
 D_T=-ZZ^*-W^*X_0W,\qquad
 X_0={\cal G}_T(F_2-VKV^*),\quad K=V^*F_2V.        \tag{19}
\]

The commutant theorem makes this a direct sum of scalar L203
identities.

Indeed, take the spectral decomposition

\[
 H_Y=\sum_\lambda\lambda E_\lambda.                \tag{20}
\]

Every \(E_\lambda\) commutes with \(S,S^*\), hence with

\[
 S^*S,\quad SS^*,\quad VV^*,\quad WW^*,\quad P.
\]

It also commutes with \(T=P^{-1/2}SP^{1/2}\), and therefore with the
physical ellipse coefficients \(E_1,E_2\).  The canonical first defect
motion respects the same decomposition: if
\(p_\lambda^V=V^*E_\lambda V\), then its formula in L203 gives
\(E_\lambda C_1=C_1p_\lambda^V\).  Consequently
\(C_1C_1^*\), \(F_2\), the balanced first tangent, and its normal
corner \(Z\) all reduce blockwise.

The restriction of \(S\) to \(E_\lambda\mathcal H\) is again a pure
partial isometry.  Its left and right defect dimensions are equal
(both equal the rank loss of the square restricted operator).
L203's trace calculation used only these partial-isometry relations
and the three-eigenvalue formula for \(P\), so its exact scalar
identity applies on every such reducing block.

Multiply the block identity by \(\lambda\) and sum.  In physical
coordinates this gives

\[
\boxed{
\operatorname{tr}\{
 {\cal Z}_Y(F_2-VKV^*)\}
=-\operatorname{tr}(YZZ^*).}                       \tag{21}
\]

The subtraction \(VKV^*\) is essential: on each block it changes
\(4P^{-1}\) into L203's dual
\(4(P^{-1}-VV^*)\).

By Stein adjointness and the definition of \(X_0\),

\[
\operatorname{tr}\{YW^*X_0W\}
=\operatorname{tr}\{
 {\cal Z}_Y(F_2-VKV^*)\}.                          \tag{22}
\]

Combining (19), (21), and (22) proves

\[
\operatorname{tr}(YD_T)
=-\operatorname{tr}(YZZ^*)
 +\operatorname{tr}(YZZ^*)=0.
\]

This is (18).  The finite-dimensional Fredholm alternative applied to
L204's real-linear map now proves existence of \(C\), hence (1), at
every fixed equality anchor.

## 6. Audit

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_elliptic_commutant.py \
  --output \
  experiments/repeated_crabb_elliptic_commutant_s70224.jsonl
```

The checker constructs the complete Hermitian cokernel at
noncommuting equality anchors of lengths \(2,\ldots,5\), copy
multiplicities \(2,3\), and two equality amplitudes.  It verifies:

1. the adjoint-kernel equation (8);
2. the state commutator (11);
3. both defect intertwining equations (12);
4. transfer intertwining (14) through degree \(2L+8\); and
5. the weighted compatibility (18).

The deterministic multiplicity-two colligations have a
one-dimensional scalar cokernel.  The multiplicity-three examples
have a two-dimensional reducing cokernel, so the audit genuinely
tests more than L203's trace direction.  The tracked dataset SHA-256
is
`d3f81984b52e7c6261928e95a73b4ba8992c722f1027af9db26e5165ec4feb49`.

## 7. Subsequent development

L206 upgrades the oriented grade-one matrix face from numerical to
pointwise proved.  L207 subsequently removes the selection debt as
well: it gives the explicit analytic preimage

\[
 \widehat C=-\frac72(I-VV^*)SWB_1.
\]

Thus the remaining problem is the higher associated-graded iteration
through \(B_2,\ldots,B_L\), followed by the disk/circular merger.
