# The flagged transfer channel is covariant under deflation

## 1. Result (L211, 2026-07-24)

Retain L209--L210's pure partial isometry \(S\), defect frames \(V,W\),
transfer coefficients and transfer channel

\[
\begin{aligned}
 B_n&=W^*(S^*)^nV,\\
 {\mathfrak C}_S(K)&=\sum_{n\ge1}B_nKB_n^*,
\end{aligned}
\]

and an isometry \(U:\mathbb C^d\to\mathbb C^m\) satisfying

\[
 B_j^*U=0\qquad(1\le j<k).                         \tag{1}
\]

Let \(S_{k,U}\) be the L209 restriction obtained by removing

\[
 WU,\ SWU,\ldots,S^{k-2}WU,
\]

and use its left defect frame

\[
 \widetilde W=[WU_\perp,S^{k-1}WU].
\]

Write its transfer coefficients and channel as
\(\widetilde B_n\) and \({\mathfrak C}_{S_{k,U}}\).  L210 proves

\[
\boxed{
 \widetilde B_n=
 \begin{bmatrix}
  U_\perp^*B_n\\
  U^*B_{n+k-1}
 \end{bmatrix}.}                                   \tag{2}
\]

Let \(J_U:\mathbb C^d\to
\mathbb C^{m-d}\oplus\mathbb C^d\) be the bottom-coordinate
inclusion.  Then for every copy matrix \(K\),

\[
\boxed{
 J_U^*{\mathfrak C}_{S_{k,U}}(K)J_U
 =U^*{\mathfrak C}_S(K)U.}                         \tag{3}
\]

Also,

\[
\boxed{
 J_U^*\widetilde B_1\widetilde B_1^*J_U
 =U^*B_kB_k^*U.}                                   \tag{4}
\]

In particular L208's proved flagged correction is exactly the bottom
compression of the corresponding deflated grade-one channel
coboundary:

\[
\boxed{
\begin{aligned}
&28J_U^*\left\{
 {\mathfrak C}_{S_{k,U}}
  (\widetilde B_1^*J_UJ_U^*\widetilde B_1)
 -\widetilde B_1\widetilde B_1^*
 \right\}J_U\\
&\qquad
=28\,U^*\{
   {\mathfrak C}_S(B_k^*UU^*B_k)-B_kB_k^*
   \}U.
\end{aligned}}                                     \tag{5}
\]

Thus L208's dual preimage and L209--L210's state-space deflation are
the same correction mechanism, including the completely positive
channel term and every multiplication order.

L211 still does **not** identify the complete prepared physical
Riemann/metric **base endpoint** after least-squares removal of prior
active rows.  The full L207 base contains the right Gram
\(\widetilde B_1^*\widetilde B_1\), including its top-row
contribution.  Proving that the later Schur quotient replaces this by
\(\widetilde B_1^*J_UJ_U^*\widetilde B_1\) is the remaining physical
endpoint-functoriality gate.  L211 must not be read as assuming that
replacement.

## 2. Channel compression

By (2), the bottom channel block is

\[
\begin{aligned}
 J_U^*{\mathfrak C}_{S_{k,U}}(K)J_U
 &=\sum_{n\ge1}
 U^*B_{n+k-1}KB_{n+k-1}^*U\\
 &=U^*\left\{
 {\mathfrak C}_S(K)
 -\sum_{j=1}^{k-1}B_jKB_j^*
\right\}U.
\end{aligned}                                      \tag{6}
\]

Every subtracted term vanishes because \(U^*B_j=0\), proving (3).
The bottom row of

\[
 \widetilde B_1=
 \begin{bmatrix}
 U_\perp^*B_1\\
 U^*B_k
 \end{bmatrix}
\]

immediately proves (4).  Finally,

\[
 \widetilde B_1^*J_UJ_U^*\widetilde B_1
 =B_k^*UU^*B_k,
\]

so substituting this \(K\) into (3) and using (4) proves (5).

## 3. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_transfer_channel_covariance.py \
  --output \
  experiments/repeated_crabb_transfer_channel_covariance_s70224.jsonl
```

The checker independently re-audits L210's coefficient and resolvent
shift on gauged heterogeneous delay sums and on unstructured partial
isometries enlarged by non-coordinate wandering chains.  It verifies
(3) by independent Stein solves and (4)--(5) directly.  The tracked
dataset SHA-256 is
`4b5f7d4d9bf7c926eba47ffd3361b26c741a4422c2e082678275d787e207ba59`.
The exact equations above, not the floating tests, prove L211.
