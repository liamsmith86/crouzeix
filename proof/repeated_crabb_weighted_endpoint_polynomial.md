# Physical reverse-edge polynomials have an exact weighted transfer factor

## 1. Result (L313, 2026-07-26)

Let

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

and put

\[
B_j=W^*(S^*)^jV,\qquad
J=(I+F)S^*(I+E).                                  \tag{1}
\]

The letter \(J\) is the physical reverse edge in the pulled ellipse
pencil.  For a word \(w\) in the alphabet \(S,J\), let \(r(w)\) be
the number of \(J\)'s and set

\[
K_w=W^*w(S,J)V.
\]

There are explicit polynomial copy matrices
\(C_1(w),\ldots,C_{r(w)}(w)\) such that

\[
\boxed{K_w=\sum_{j=1}^{r(w)}B_jC_j(w)}            \tag{2}
\]

and, if \(d=|w|\),

\[
\boxed{
\left\|\begin{bmatrix}C_1(w)\\ \vdots\\
C_{r(w)}(w)\end{bmatrix}\right\|
\le(d+1)4^{r(w)}.}                                \tag{3}
\]

In particular, every physical word with at most \(r\) reverse edges
factors boundedly through the first \(r\) transfer blocks.  This is
stronger weight information than applying L312 after expanding the
defect projections in \(J\): that expansion can create several raw
\(S^*\)'s from one physical reverse edge and hides the cancellation
proved below.

The factorization has an exact weighted form.  Let

\[
{\cal P}(c)=\sum_{\nu=1}^N
\alpha_\nu c^{d_\nu}w_\nu(S,J),\qquad
d_\nu\ge r_\nu:=r(w_\nu),
\quad 0<c\le1,                                    \tag{4}
\]

and let \(R=\max_\nu r_\nu\).  Then

\[
\boxed{
W^*{\cal P}(c)V
=[cB_1\ \ c^2B_2\ \cdots\ c^RB_R]\,{\cal C}(c),} \tag{5}
\]

where

\[
\boxed{
\|{\cal C}(c)\|
\le\sum_{\nu=1}^N
|\alpha_\nu|(|w_\nu|+1)4^{r_\nu}
c^{d_\nu-r_\nu}.}                                \tag{6}
\]

Consequently

\[
\boxed{
\{W^*{\cal P}(c)V\}\{W^*{\cal P}(c)V\}^*
\preceq \|{\cal C}(c)\|^2
\sum_{j=1}^Rc^{2j}B_jB_j^*.}                     \tag{7}
\]

This is polynomial and rank-stable.  It uses no transfer inverse,
pseudoinverse, Schur projection, or minimum singular value.

## 2. Wandering-chain recursion

It is clearest to factor the adjoint endpoint.  Apply \(w^*\) to
\(Wy\), in the same chronological letter order as \(w\).  The two
adjoint letters are

\[
S^*,\qquad J^*=(I+E)S(I+F).
\]

For the forward chain \(S^hW\), partial-isometry algebra gives

\[
\boxed{
\begin{aligned}
S^*W&=0,\\
S^*S^hW
 &=S^{h-1}W-VB_{h-1}^*, &&h\ge1,\\
J^*W
 &=2SW+2VB_1^*,\\
J^*S^hW
 &=S^{h+1}W+VB_{h+1}^*, &&h\ge1.
\end{aligned}}                                    \tag{8}
\]

For the second line, use

\[
S^*S^hW=(I-E)S^{h-1}W
\]

and \(V^*S^{h-1}W=B_{h-1}^*\).  For the last two
lines, use \(FW=W\),

\[
FS^hW=0\quad(h\ge1),
\]

and \(ES^hW=VB_h^*\).

Starting with height zero, follow only the first, chain-valued term
on the right of (8).  Each \(J^*\) raises the height by one and each
\(S^*\) lowers it by one; an \(S^*\) at height zero kills the main
path.  Every discarded term is exactly

\[
\lambda\,VB_j^*y,\qquad 1\le j\le r(w).           \tag{9}
\]

Apply the remaining adjoint-letter suffix to (9), then compress by
\(V^*\).  It contributes

\[
\lambda V^*R_{\rm suffix}V\,B_j^*y.              \tag{10}
\]

If the main path survives to height \(h>0\), its final contribution
is \(\lambda B_h^*y\), again with \(h\le r(w)\).
Collecting (10) by \(j\), then taking adjoints, constructs (2)
explicitly.

There are at most \(d+1\) emitted terms.  The main coefficient can
only double when \(J^*\) leaves height zero, and

\[
\|J^*\|\le\|I+E\|\,\|S\|\,\|I+F\|\le4.
\]

Thus every emitted multiplier has norm at most \(4^{r(w)}\).
The triangle inequality for the stacked column proves (3).

## 3. Weighted assembly

Apply (2) to one term in (4):

\[
\alpha c^dK_w
=\sum_{j=1}^{r}c^jB_j
\{\alpha c^{d-j}C_j(w)\}.                         \tag{11}
\]

Because \(j\le r\le d\) and \(0<c\le1\),

\[
c^{d-j}\le c^{d-r}.
\]

Stack the coefficients in (11), use (3), and then sum over
\(\nu\).  This gives (5)--(6).  Equation (7) follows from
\(XX^*\preceq\|{\cal C}\|^2{\mathbb B}{\mathbb B}^*\).

## 4. Direct elliptic-map corollary

L125 writes the direct map as

\[
\phi_c(z)=\sum_{h\ge0}a_h(c)z^{2h+1},
\qquad \operatorname {ord}_c a_h\ge h.            \tag{12}
\]

The complete physical pulled operator is

\[
A(c)=\phi_c(S+cJ).                                \tag{13}
\]

Expand one power in (12).  A word containing \(r\) copies of \(J\)
carries \(c^r\); multiplying by a coefficient of \(c\)-degree \(t\)
in \(a_h\) gives total degree

\[
d=t+r\ge r.                                       \tag{14}
\]

For a fixed degree \(d\), (12) has only finitely many contributing
powers because \(t\ge h\).  Therefore every finite jet of
\(A(c)-S\), and every endpoint polynomial built from its physically
grouped \(S+cJ\) expansion without a Stein closure, satisfies
(5)--(7).

This closes the finite conformal-map part of the weighted
L306/L307 multiplier problem.  It also explains why expanding
\(E,F\) inside \(J\) and counting raw \(S^*\)'s is the wrong
filtration.

## 5. Exact scope guard

L313 does **not** complete the prepared recurrence.  The raw defect
frame is obtained after a stable Stein solve and factorization.
L217's early future-transfer images show that its coefficients do
not separately obey (14).  Hence one may now localize the remaining
weighted debt:

> prove that the frame terms
> \(C^*H+H^*C\) cancel their noncausal part against the
> \(\Psi_S(\Delta^*XS+S^*X\Delta)\) terms before L313 is applied.

The favorable squares \(-\Psi_S(\Delta^*X\Delta)\) when their sign
is available and \(-C^*C\) must remain grouped as in L306/L307.
One must not apply L313 termwise to the raw Hermitian-gauge frame,
nor infer the all-series recurrence from the direct-map corollary.

Thus the live obstruction is no longer all of L306: it is the
paired **Stein/frame covariance** that removes the premature images.
L220/L236/L244 are the natural coordinates for that pairing.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_weighted_endpoint_polynomial.py \
  --output \
  experiments/repeated_crabb_weighted_endpoint_polynomial_s70226.jsonl
```

The checker constructs (2) by the recursion (8)--(10), not by a
pseudoinverse.  It audits random causal physical polynomials and the
complete direct-map word jet through degree six at three ellipse
parameters, on unstructured, rank-chain, and completely delayed
colligations of multiplicities one to three.  Each record also checks
all 510 physical words through length eight individually.  A separate
comparison reconstructs the physical direct-map coefficients from the
repository's balanced implementation and agrees through degree six.
All 42 records pass.  The exact identities above, not the floating
audit, prove L313.  The tracked dataset SHA-256 is

```text
5e20fe4db61c1b80582436faa3c2c98b2d6699fd329d6b006f895196523e0252
```
