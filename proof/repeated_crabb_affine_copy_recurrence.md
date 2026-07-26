# Recursive retightening has an exact partial-scale Gram reserve

## 1. Result (L307, 2026-07-26)

Retain L306's balanced stable partial isometry, moving pair, and copy
closure

\[
I-S^*S=VV^*,\qquad A=S+\Delta,\qquad D=V+H,
\qquad
\Psi_S(Z)=V^*\sum_{n\ge0}S^nZ(S^*)^nV.           \tag{1}
\]

Let \(X=X^*\), \(C\), and a Hermitian forcing \(R\) obey the affine
fixed-base equation

\[
\boxed{X-S^*XS=R+VC^*+CV^*.}                    \tag{2}
\]

This is the equation encountered after a new correction is chosen to
cancel the preceding residual \( -R\).  At the moving pair,

\[
\begin{aligned}
&X-A^*XA-DC^*-CD^*-CC^*\\
&\qquad=R+{\cal N}(X,C),                          \tag{3}\\
{\cal N}(X,C)
&=-\Delta^*XS-S^*X\Delta-\Delta^*X\Delta\\
&\quad-HC^*-CH^*-CC^*.                           \tag{4}
\end{aligned}
\]

Thus adding (3) after the preceding residual \(-R\) cancels it
exactly and leaves only the same six-term successor (4).  The
recurrence is affine; L306 is not restricted to its first
homogeneous use.

The initial-copy boundary term is exact:

\[
\boxed{
\Psi_S(VC^*+CV^*)=C^*V+V^*C.}                    \tag{5}
\]

Consequently

\[
\boxed{
\begin{aligned}
&\Psi_S(X-A^*XA)
 -(C^*D+D^*C+C^*C)\\
&\qquad=\Psi_S(R)+{\mathfrak q}(X,C),             \tag{6}\\
{\mathfrak q}(X,C)
&=-\Psi_S(\Delta^*XS+S^*X\Delta+\Delta^*X\Delta)\\
&\quad-(C^*H+H^*C+C^*C).
\end{aligned}}                                    \tag{7}
\]

In particular L306's apparently separate fixed-base copy term is
identically zero when \(R=0\), not merely response-equivalent to
zero.

There is also an exact reason to use a fixed partial scale.  For real
\(0\le\theta\le1\), scale the entire affine correction
\((R,X,C)\) by \(\theta\).  Then

\[
\boxed{
\begin{aligned}
{\cal N}(\theta X,\theta C)
&=\theta{\cal N}(X,C)
  +\theta(1-\theta)CC^*,\\
{\mathfrak q}(\theta X,\theta C)
&=\theta{\mathfrak q}(X,C)
  +\theta(1-\theta)C^*C .
\end{aligned}}                                    \tag{8}
\]

Both reserve terms are positive semidefinite.  At
\(\theta=1/2\), the state and copy recurrences retain exactly one
quarter of the correction-frame Gram relative to one half of the
full-scale successor.

There is an important scope restriction.  If \(R\ne0\), scaling the
entire affine correction cancels only \(\theta R\) of the preceding
\(-R\); their combination also retains \(-(1-\theta)R\).  Therefore
(8) is an exact positive term in the **successor scaling identity**,
not by itself a positive-margin induction through arbitrary affine
steps.  For the first homogeneous L298 direction, \(R=0\), so no such
leftover occurs.

L307 supplies the arbitrary-step algebra and an exact scaling reserve.
It does **not** prove that the full-scale successor
\({\mathfrak q}(X,C)\) is favorable.  L309 later disproves global
prepared odd-response parity.  The remaining gate is to show that
each lower-neutral odd retained class is a bounded hereditary factor
through the cumulative active transfer row, and that each even
retained part is quadratic in prior transfer rows, so the reserve and
the earlier strict margins can absorb them.

## 2. Proof

Expand \(A=S+\Delta\) and \(D=V+H\) in the left side of (3).  The
constant terms are

\[
X-S^*XS-VC^*-CV^*=R
\]

by (2), and the six remaining terms are exactly (4).

For (5), the \(n=0\) term of (1) is \(C^*V+V^*C\).
Every \(n\ge1\) term vanishes: \(S^nV=0\), and its adjoint gives
\(V^*(S^*)^n=0\).  Applying \(\Psi_S\) to (2) therefore gives

\[
\Psi_S(X-S^*XS)=\Psi_S(R)+C^*V+V^*C.
\]

Subtract the moving operator and frame increments to obtain (6)--(7).

Finally, every term of (4) and (7) is linear in \((X,C)\) except the
negative squares \(-CC^*\) and \(-C^*C\).  Hence

\[
\theta L-\theta^2G
=\theta(L-G)+\theta(1-\theta)G
\]

with \(G=CC^*\) or \(C^*C\), proving (8).

## 3. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_affine_copy_recurrence.py \
  --output \
  experiments/repeated_crabb_affine_copy_recurrence_s70226.jsonl
```

The checker starts with L298's simultaneous weighted correction,
uses its L299 residual as the forcing canceled by a second
bridge-valued correction, and verifies (3), (6), and (8) independently
at three partial scales.  Its unstructured, rank-collapse, and
nonvacuous complete-delay cases exercise correction grades through
three.  The tracked dataset SHA-256 is

```text
0b5f54dfb20b9fba53971ab8d03f7daad7833f4b020b7ee0c8fcf80b9dc4e385
```
