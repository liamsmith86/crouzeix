# Every Schur transport is a graph gauge plus its endpoint residual

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

> **Route update.**  L292 subsequently replaces the fixed-rank
> pseudoinverse existence factor below by an exact analytic
> Smith-valuation criterion along every L197 failure arc.  Pointwise
> flag zero is not enough for bounded selection; the physical
> endpoint residual must satisfy the relative valuation bounds.

## 1. Result (L291, 2026-07-26)

Let \(H,\widetilde H\) be Hermitian state matrices split into a fixed
endpoint of dimension \(m\) and its interior:

\[
H=\begin{bmatrix}A&B\\B^*&D\end{bmatrix},\qquad
\widetilde H=
\begin{bmatrix}\widetilde A&\widetilde B\\
\widetilde B^*&\widetilde D\end{bmatrix}.          \tag{1}
\]

Assume \(D,\widetilde D\succ0\), and write their Schur endpoints as

\[
{\cal U}(H)=A-BD^{-1}B^*,\qquad
{\cal U}(\widetilde H)
=\widetilde A-\widetilde B\widetilde D^{-1}
 \widetilde B^*.                                  \tag{2}
\]

Then there is an invertible block-lower-triangular matrix

\[
R=\begin{bmatrix}I&0\\K&C\end{bmatrix}            \tag{3}
\]

such that

\[
\boxed{
\widetilde H
=R^*HR+
\begin{bmatrix}
{\cal U}(\widetilde H)-{\cal U}(H)&0\\
0&0
\end{bmatrix}.}                                   \tag{4}
\]

The first term in (4) is a pure graph gauge: the identity block and
zero upper-right block in (3) make its Schur endpoint exactly
\({\cal U}(H)\).  Thus every state-space complication away from the
literal endpoint residual can be absorbed into one bounded
endpoint-quotient-fixing congruence.

Now collect the first \(r\) transfer rows into

\[
{\mathbb B}=[B_1\ \cdots\ B_r]:
\mathbb C^{rm}\longrightarrow\mathbb C^m,\qquad
{\cal K}_r=\ker{\mathbb B}^*
=\bigcap_{j=1}^r\ker B_j^*.                       \tag{5}
\]

For a Hermitian endpoint matrix \(Q\), the following are equivalent:

\[
\boxed{
\begin{aligned}
P_{{\cal K}_r}QP_{{\cal K}_r}=0
\quad\Longleftrightarrow\quad
Q=X{\mathbb B}^*+{\mathbb B}X^*
\quad\hbox{for some }X.
\end{aligned}}                                    \tag{6}
\]

Combining (4) and (6), a complete physical transport preserves the
surviving transfer flag if and only if it has the exact normal form

\[
\boxed{
\widetilde H
=R^*HR+
\begin{bmatrix}
X{\mathbb B}^*+{\mathbb B}X^*&0\\0&0
\end{bmatrix}.}                                   \tag{7}
\]

This sharpens L290's architecture.  A hereditary **state** change is
a convenient sufficient certificate, but it is stronger than
necessary.  The exact invariant lives at the mixed-graph Schur
endpoint.  The large interior metric witnesses need not themselves
belong to a delay ideal.

L291 does not prove the physical arbitrary-grade endpoint
factorization or boundedness of \(X\) through rank changes.  It
removes the nonendpoint state algebra from that debt:

> A178 should construct or estimate the finite endpoint factors
> \(X\), while using the graph gauge \(R\) only as an automatically
> bounded coordinate change.

The block normal form and (6) are elementary linear algebra and are
not claimed as new general matrix theorems.  Their use to isolate the
only non-gauge content of the repeated-Crabb transport is the new
campaign reduction.

## 2. Exact graph-gauge construction

Define the unit lower elimination factors

\[
L_H=
\begin{bmatrix}
I&0\\D^{-1}B^*&I
\end{bmatrix},\qquad
L_{\widetilde H}=
\begin{bmatrix}
I&0\\
\widetilde D^{-1}\widetilde B^*&I
\end{bmatrix}.                                    \tag{8}
\]

Direct block multiplication gives

\[
\boxed{
\begin{aligned}
H&=L_H^*
\begin{bmatrix}{\cal U}(H)&0\\0&D\end{bmatrix}
L_H,\\
\widetilde H&=L_{\widetilde H}^*
\begin{bmatrix}{\cal U}(\widetilde H)&0\\
0&\widetilde D\end{bmatrix}
L_{\widetilde H}.
\end{aligned}}                                    \tag{9}
\]

Choose any invertible \(C\) satisfying

\[
C^*DC=\widetilde D.                               \tag{10}
\]

For positive interiors, one canonical choice is

\[
\boxed{
C=D^{-1/2}
\left(D^{1/2}\widetilde D D^{1/2}\right)^{1/2}
D^{-1/2}.}                                        \tag{11}
\]

Put

\[
\boxed{
R=L_H^{-1}
\begin{bmatrix}I&0\\0&C\end{bmatrix}
L_{\widetilde H}.}                                \tag{12}
\]

Every factor on the right of (12) is lower triangular with identity
endpoint block, so \(R\) has form (3).  Equations (9)--(12) give

\[
R^*HR
=L_{\widetilde H}^*
\begin{bmatrix}{\cal U}(H)&0\\0&\widetilde D\end{bmatrix}
L_{\widetilde H}.                                 \tag{13}
\]

Subtract (13) from the second line of (9).  Since the elimination
factor has an identity endpoint block, the remaining diagonal
endpoint term is unchanged by its congruence.  This proves (4).

If \(H,\widetilde H\) are analytic families with uniformly positive
interiors, the principal square roots in (11) are analytic.  Hence
\(C\), \(R\), and their inverses are locally analytic and bounded.
The graph gauge creates no rank-change singularity.

The same proof works for formal Hermitian series whose interior
constant coefficients are positive definite: inverse and principal
square-root series are unique coefficientwise.

## 3. Endpoint hereditary factorization is equivalent to flag zero

Let \({\cal R}=\operatorname{ran}{\mathbb B}\) and let
\(\Pi\) be the orthogonal projection onto \({\cal R}\).  Then
\({\cal K}_r={\cal R}^{\perp}\).

If

\[
Q=X{\mathbb B}^*+{\mathbb B}X^*,
\]

both terms vanish after compression to \({\cal K}_r\), proving the
reverse implication in (6).

Conversely, assume

\[
(I-\Pi)Q(I-\Pi)=0.                                \tag{14}
\]

Let \({\mathbb B}^{\dagger}\) be the Moore--Penrose inverse.  It is
used only for this fixed-matrix existence proof, and

\[
{\mathbb B}{\mathbb B}^{\dagger}=\Pi.
\]

Set

\[
\boxed{
X=
\left\{
\frac12\Pi Q\Pi+(I-\Pi)Q\Pi
\right\}
({\mathbb B}^{\dagger})^*.}                       \tag{15}
\]

Then

\[
\begin{aligned}
X{\mathbb B}^*
&=\frac12\Pi Q\Pi+(I-\Pi)Q\Pi,\\
{\mathbb B}X^*
&=\frac12\Pi Q\Pi+\Pi Q(I-\Pi).
\end{aligned}
\]

Their sum is \(Q\) by (14), proving (6).

Equation (15) is **not** an admissible campaign selection across a
rank-changing family: the pseudoinverse can be discontinuous and
unbounded.  It proves only pointwise equivalence.  A178 still needs
one of:

1. an explicit polynomial/analytic endpoint factor \(X\), such as
   L288's factor at sextic order; or
2. L282's uniform energy estimate, which supplies a bounded physical
   correction without selecting a rank projection.

## 4. Consequence for the transport campaign

L285 shows that raw state forcing does not remain in a delay word
ideal.  L290 shows that exposed hereditary state channels would be
sufficient.  L291 now identifies the exact middle ground:

1. use L285 and L194 to construct an admissible transported state
   \(\widetilde H\);
2. use L289 to compute its complete endpoint difference, including
   every graph cross;
3. factor only that endpoint difference through the earlier
   transfer row as in (6); and
4. treat all remaining state motion as the bounded gauge \(R\).

In particular, failure of a metric witness to lie in a state delay
ideal is no longer an obstruction.  Conversely, membership of a raw
state polynomial in a two-sided ideal is not a substitute for the
endpoint factorization.

The normal form does not allow an arbitrary \(R\) to be inserted
into the physical construction.  The two states must first arise
from valid lower-tight Stein/factor data.  Equation (4) then
classifies their difference; it does not manufacture admissibility.

The remaining mathematical target is therefore precise:

\[
\boxed{
[c^n]\{
{\cal U}(\widetilde H)-{\cal U}(H)
\}
=X_n{\mathbb B}_{r(n)}^*
+{\mathbb B}_{r(n)}X_n^*,}                        \tag{16}
\]

with \(X_n\) bounded by an arbitrary-grade rule and with the retained
even direct-Gram/Schur margin positive through the terminal grade.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_gauge_normal_form.py \
  --output \
  experiments/repeated_crabb_endpoint_gauge_normal_form_s70226.jsonl
```

The checker uses complex positive interiors, independently moving
cross blocks, and transfer rows of ranks one through four.  It
verifies:

1. both exact elimination factorizations (9);
2. the canonical interior congruence (11);
3. the quotient-fixing triangular form of \(R\);
4. the normal form (4);
5. the literal endpoint difference;
6. flag compression; and
7. the converse factorization (15).

The tracked dataset has SHA-256

```text
c7eea10f5d4e58b25a13fe67040d671816245dc0be8b5fafccc276a61b997bb0
```
