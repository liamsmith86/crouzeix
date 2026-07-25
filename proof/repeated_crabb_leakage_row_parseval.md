# The complete future leakage row has exactly its diagonal energy

## 1. Result (L261, 2026-07-25)

Retain L259's square matrix-inner transfer

\[
B(z)=\sum_{j\geq1}B_jz^j
\]

and causal Toeplitz multiplier \({\cal T}_B\).  Put

\[
{\cal L}_B={\cal T}_B{\cal T}_B^*.
\]

Matrix innerness makes \({\cal T}_B\) an isometry, so
\({\cal L}_B\) is an orthogonal projection.  If

\[
B_1=\cdots=B_{k-1}=0,
\]

then L259 gives the first nonzero row

\[
P_k{\cal L}_BP_j=B_kB_j^*,\qquad j\geq k.
\]

The norm of this **entire future row**, not merely its diagonal
block, is exactly the desired active energy:

\[
\boxed{
\operatorname {tr}
\{P_k{\cal L}_B^2P_k\}
=\operatorname {tr}
\{P_k{\cal L}_BP_k\}
=\|B_k\|_F^2.}                                    \tag{1}
\]

Equivalently, in copy coordinates,

\[
\boxed{
\sum_{j\geq k}
B_kB_j^*B_jB_k^*
=B_kB_k^*.}                                      \tag{2}
\]

Thus L243's future coefficients need not cancel term by term.
They may close as one Hardy-row norm, after which matrix-inner
Parseval collapses the whole row to its diagonal energy.

L261 does **not** prove that L258's physical return is the row norm
in (1).  It sharpens A194's remaining statement to

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathfrak D_{\rm ret}
=4\operatorname {tr}
\{P_k{\cal L}_B^2P_k\},}                           \tag{3}
\]

together with the lower vanishings.  Compared with L259's
diagonal-selection formulation, (3) is the form naturally suggested
by a closed entry/return path: keep every future block, prove that
the return supplies its adjoint row with the unweighted Hardy
metric, and use projection/Parseval only at the end.

## 2. Proof

Since \(B\) is matrix inner, multiplication by \(B\) is an isometry:

\[
{\cal T}_B^*{\cal T}_B=I.
\]

Therefore

\[
{\cal L}_B^2
={\cal T}_B({\cal T}_B^*{\cal T}_B){\cal T}_B^*
={\cal L}_B.
\]

Compressing this projection identity to the individual Hardy row
\(P_k\) proves the first equality in (1).  L259 gives

\[
P_k{\cal L}_BP_k=B_kB_k^*,
\]

whose copy-space trace proves the second equality.

Expanding the middle row/column product in
\(P_k{\cal L}_B^2P_k\) gives

\[
\sum_j
(P_k{\cal L}_BP_j)(P_j{\cal L}_BP_k).
\]

Complete delay kills \(j<k\); L259 identifies the two remaining
blocks as \(B_kB_j^*\) and \(B_jB_k^*\).  This proves (2), including
the noncommutative multiplication order.

## 3. Interface with the live proof

L251 proves that physical analytic ports occur in opposite pairs,
and L258 resums every final-row visit into a closed return.  L261
identifies the only acceptable final simplification of that paired
return: it must become the ordinary Hardy norm of L259's full
leakage row.  No individual future block is zero, and L260 shows
that the two ellipse orientations cannot be separated before this
closure.

The remaining proof obligation is therefore a metric statement, not
an autocorrelation guess:

> after L244 removes the half-line coisometric graph, show that the
> first physical return uses the standard Hardy inner product on
> the leakage row, with the doubled remote amplitude on both sides.

If this is established, (1) supplies the future-coefficient
cancellation and the factor \(2^2=4\) gives (3).
