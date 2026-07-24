# Wandering-chain deflation and exact reflected transfer pullback

## 1. Results (L209--L210, 2026-07-24)

Let \(S\) be a finite pure partial isometry with equal orthogonal
defect frames

\[
 I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0,
\]

and genuine transfer coefficients

\[
 B_j=W^*(S^*)^jV,\qquad B_0=0.                    \tag{1}
\]

Fix an isometry \(U:\mathbb C^d\to\mathbb C^m\) and an integer
\(k\ge2\).  Suppose

\[
\boxed{B_j^*U=0\quad(1\le j<k).}                  \tag{2}
\]

Then the state columns

\[
 WU,\ SWU,\ \ldots,\ S^{k-1}WU                  \tag{3}
\]

are isometries with pairwise orthogonal ranges.  Remove the first
\(k-1\) of these ranges:

\[
\begin{aligned}
 {\cal L}_{k,U}
 &=\bigoplus_{j=0}^{k-2}S^jWU\mathbb C^d,\\
 {\cal H}_{k,U}&={\cal L}_{k,U}^{\perp}.           \tag{4}
\end{aligned}
\]

The space \({\cal H}_{k,U}\) is invariant for \(S\), and

\[
 S_{k,U}=S|_{{\cal H}_{k,U}}                       \tag{5}
\]

is again a partial isometry with right defect \(V\mathbb C^m\).
If \(U_\perp\) is any orthonormal complement of \(U\), its left defect
frame is

\[
\boxed{
 \widetilde W
 =\bigl[WU_\perp,\ S^{k-1}WU\bigr].}               \tag{6}
\]

In the output coordinates
\(\mathbb C^{m-d}\oplus\mathbb C^d\), the first coefficient of the
deflated transfer is

\[
\boxed{
 \widetilde B_1
 =
 \begin{bmatrix}
 U_\perp^*B_1\\
 U^*B_k
 \end{bmatrix}.}                                   \tag{7}
\]

In fact every coefficient of the smaller transfer satisfies

\[
\boxed{
 \widetilde B_n
 =
 \begin{bmatrix}
 U_\perp^*B_n\\
 U^*B_{n+k-1}
 \end{bmatrix}
 \qquad(n\ge1).}                                    \tag{7a}
\]

Thus a grade-\(k\) coefficient on the common left kernel of the
earlier grades becomes a genuine grade-one coefficient of a smaller
colligation, and the entire surviving reflected series is its exact
monomial shift.

L209 is a colligation theorem.  It does **not** yet identify the
physical grade-\(k\) ellipse jet with L207 applied to (5); that
weighted-jet lifting is the remaining interface.  L210 closes the
transfer-series part of that interface, but not the physical endpoint
covariance.

## 2. The delay line is orthonormal

For every state vector \(x\),

\[
 \|Sx\|^2=\|x\|^2-\|V^*x\|^2.                     \tag{8}
\]

Also

\[
 V^*S^jWU=B_j^*U.                                  \tag{9}
\]

Equations (2), (8), and \(B_0=0\) show inductively that every column
in (3) is an isometry through level \(k-1\).

If \(0\le i<j\le k-1\), the same no-loss identities permit \(i\)
successive cancellations of \(S^*S\):

\[
\begin{aligned}
 (S^iWU)^*S^jWU
 &=U^*W^*S^{j-i}WU=0.                              \tag{10}
\end{aligned}
\]

The last equality follows from \(W^*S=0\), the adjoint of
\(S^*W=0\).  This proves the orthonormality assertion without
commuting any copy matrices.

## 3. Compression and its defects

For \(j\ge1\), equations (2) and (8) also give

\[
 S^*S^jWU=S^{j-1}WU,\qquad S^*WU=0.                \tag{11}
\]

Hence

\[
 S^*{\cal L}_{k,U}\subseteq{\cal L}_{k,U}.
\]

Taking orthogonal complements proves that
\({\cal H}_{k,U}\) is invariant for \(S\).

The old right defect lies in the retained space.  Indeed,

\[
 V^*S^jWU=B_j^*U=0\quad(0\le j\le k-2).            \tag{12}
\]

For \(x\in{\cal H}_{k,U}\), invariance and (8) therefore give

\[
 S_{k,U}^*S_{k,U}
 =I_{{\cal H}_{k,U}}-VV^*.                         \tag{13}
\]

Thus the compressed right defect is exactly the old \(m\)-dimensional
one.

Both parts of (6) lie in the retained space, are mutually orthogonal,
and are orthogonal to \(V\): for the promoted block the last fact is
exactly \(B_{k-1}^*U=0\).  Moreover,

\[
\begin{aligned}
 S_{k,U}^*WU_\perp&=0,\\
 S_{k,U}^*S^{k-1}WU
 &=P_{{\cal H}_{k,U}}S^{k-2}WU=0.                 \tag{14}
\end{aligned}
\]

So (6) supplies \(m\) orthonormal left-defect columns.  A square
finite matrix has equal left and right nullities; (13) shows the
right nullity is \(m\).  Therefore (6) is the complete left defect,
and \(S_{k,U}\) is a partial isometry.

## 4. Promotion of \(B_k\)

Use the left frame (6) and the unchanged right frame \(V\).  The top
row of the first deflated coefficient is

\[
 U_\perp^*W^*S^*V=U_\perp^*B_1.
\]

For the promoted row,

\[
\begin{aligned}
 (S^{k-1}WU)^*S^*V
 &=U^*W^*(S^*)^kV\\
 &=U^*B_k.                                         \tag{15}
\end{aligned}
\]

The compression projections may be omitted in (15) because both
defect frames lie in \({\cal H}_{k,U}\).  Stacking the two rows proves
(7).

At the repeated Crabb apex \(B_j=0\) for \(j<L\) and \(B_L\) is
unitary.  More generally, L201 shows that \(B_L\) remains invertible
on a fixed equality neighbourhood.  Consequently the descending
left-kernel flag

\[
 U_k=\bigcap_{j=1}^{k}\ker B_j^*                   \tag{16}
\]

must terminate by grade \(L\).  There is no infinite reflected
copy-space recursion.

## 5. Exact transfer and reflected-Gram pullback (L210)

Let

\[
 \widetilde B(z)=\sum_{n\ge1}\widetilde B_nz^n
\]

be the genuine matrix-inner transfer of the deflated colligation.
Because \({\cal H}_{k,U}\) is invariant for \(S\), for either part
of the left frame and every \(n\ge1\),

\[
 \widetilde B_n^*
 =V^*S_{k,U}^n\widetilde W
 =V^*S^n\widetilde W.                              \tag{17}
\]

The first block in (6) gives \(B_n^*U_\perp\), while its second block
gives \(B_{n+k-1}^*U\).  Taking adjoints proves (7a) without inserting
or commuting any copy matrices.

Let \(J_U:\mathbb C^d\to
\mathbb C^{m-d}\oplus\mathbb C^d\) denote the last-coordinate
inclusion.  Assumption (2) and (7a) give the exact analytic identities

\[
\boxed{
\begin{aligned}
 U^*B_H(z)
 &=z^{k-1}J_U^*\widetilde B(z),\\
 U^*B_H(c/z)
 &=(c/z)^{k-1}J_U^*\widetilde B(c/z).
\end{aligned}}                                     \tag{18}
\]

Consequently the left-oriented reflected Hardy Gram pulls back with
the exact weight shift

\[
\boxed{
\begin{aligned}
U^*\!\left(
\sum_{j\ge1}|c|^{2j}B_jB_j^*
\right)U
={}&|c|^{2(k-1)}
J_U^*\!\left(
\sum_{n\ge1}|c|^{2n}
\widetilde B_n\widetilde B_n^*
\right)J_U .
\end{aligned}}                                     \tag{19}
\]

This is the precise transfer-level weight-\(k\) pullback suggested by
L209.  In particular its first nonzero term is
\(|c|^{2k}U^*B_kB_k^*U\): the bottom grade-one Gram of the deflated
model multiplied by \(|c|^{2(k-1)}\).  Matrix innerness of
\(\widetilde B\) follows from the partial-isometry theorem above, so
(18)--(19) are genuine Hardy identities rather than formal coefficient
matches.

What remains is narrower: prove that the **prepared physical
similarity endpoint**, not merely its reflected transfer Gram, obeys
the same monomial-shift covariance.  That is the rank-\(m\) analogue
of L149--L150's model-complement/residue step.

## 6. Analytic use along a flag

Along a real-analytic arc, apply the finite Schur-orthogonal procedure
of L197 to the left Grams \(B_jB_j^*\).  On every constant-rank
segment, the resulting kernel bundle admits a local analytic frame
\(U(s)\).  Equations (3)--(7) then use only multiplication, adjoint,
and orthogonal projection, so the deflated colligation is analytic
on that segment.

The remaining higher-grade statement is now precise:

> show that the associated physical ellipse endpoint on
> \(U_{k-1}\) has the same monomial-shift covariance already proved
> for its transfer Gram in (18)--(19).

Once that weighted-jet identity is established, L207 supplies the
negative left Gram \(-16B_kB_k^*\) on each successive quotient, and
invertibility of \(B_L\) closes the elliptic flag.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_transfer_deflation.py \
  --output \
  experiments/repeated_crabb_transfer_deflation_s70224.jsonl
```

The checker uses heterogeneous direct sums of finite shifts, random
state conjugacies, and independent random left/right defect gauges.
It verifies (2)--(7a) and (19) for non-coordinate flag subspaces at several
intermediate and terminal grades.  The tests audit the construction;
equations (8)--(19) are the all-size proof.  The tracked dataset
SHA-256 is
`89ba19a7f78ea471b8605c5aade4a69308eddabce59b845918bb0dd655daae5a`.
