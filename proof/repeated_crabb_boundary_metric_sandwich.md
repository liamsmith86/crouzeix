# The delayed boundary metric is strictly inside the \(1\)--\(4\) sandwich

## 1. Result (L219, 2026-07-24)

Let \(S\) be a finite pure partial isometry with equal orthogonal
defect projections

\[
 E=VV^*=I-S^*S,\qquad F=WW^*=I-SS^*,\qquad EF=0.
\]

Use the repeated-Crabb equality metric

\[
 P=2I-E+2F
\]

and, for \(0<c<1\), put \(q=c^2\) and

\[
\boxed{
 P_{\rm bl}(c)
 =I-\sum_{j\ge1}\frac{q^j}{1+q^j}(S^*)^jES^j
   +\sum_{j\ge1}q^jS^jF(S^*)^j.}                 \tag{1}
\]

Then the physical boundary-layer metric

\[
 P_{\rm phys}(c)=P^{1/2}P_{\rm bl}(c)P^{1/2}
\]

obeys the strict operator sandwich

\[
\boxed{I<P_{\rm phys}(c)<4I.}                     \tag{2}
\]

In particular,

\[
\kappa(P_{\rm phys}(c))<4.                         \tag{3}
\]

This proves that the explicit metric in the all-grade one-image
candidate has exactly the required strict condition-number geometry,
for every pure partial isometry and without any transfer-delay
hypothesis.  The remaining obstruction is solely the contraction
condition: in general

\[
P_{\rm bl}-\widehat T_c^*P_{\rm bl}\widehat T_c
\]

is indefinite for the balanced elliptic pullback \(\widehat T_c\).
Thus (1) is not by itself a similarity certificate.  Its Stein slack
must still be factored or repaired by the ordered defect-frame
construction in L214--L218.

## 2. The two orbit resolutions

Purity gives norm-convergent telescoping identities

\[
\boxed{
\begin{aligned}
 I&=\sum_{j\ge0}(S^*)^jES^j,\\
 I&=\sum_{j\ge0}S^jF(S^*)^j.
\end{aligned}}                                    \tag{4}
\]

Indeed, the first partial sum is

\[
\sum_{j=0}^{N}(S^*)^j(I-S^*S)S^j
=I-(S^*)^{N+1}S^{N+1},
\]

and the second is its left-defect analogue.  Since a finite pure
contraction has spectral radius below one, both remainders converge
to zero in norm.

Write

\[
 E_j=(S^*)^jES^j,\qquad F_j=S^jF(S^*)^j.
\]

Every \(E_j,F_j\) is positive semidefinite, and (4) says that each
family resolves the identity.

## 3. Exact positive decompositions of both gaps

Because \(E,F\) are orthogonal projections,

\[
\boxed{P^{-1}=\frac12I+\frac12E-\frac14F.}          \tag{5}
\]

Set

\[
a_j=\frac{q^j}{1+q^j},\qquad b_j=q^j.
\]

Using the first resolution in (4), direct subtraction gives

\[
\boxed{
P_{\rm bl}-P^{-1}
=\sum_{j\ge1}\left(\frac12-a_j\right)E_j
 +\frac14F+\sum_{j\ge1}b_jF_j.}                   \tag{6}
\]

Using the second resolution gives

\[
\boxed{
4P^{-1}-P_{\rm bl}
=2E+\sum_{j\ge1}a_jE_j
 +\sum_{j\ge1}(1-b_j)F_j.}                        \tag{7}
\]

For \(0<q<1\),

\[
0<a_j<\frac12,\qquad 0<b_j<1.
\]

All coefficients in (6)--(7) are therefore positive.  Formula (6)
contains the complete left-orbit resolution \(F_0,F_1,\ldots\) with
positive coefficients, while (7) contains the complete right-orbit
resolution \(E_0,E_1,\ldots\) with positive coefficients.  Hence
both operators in (6)--(7) are positive definite.

Finally,

\[
\begin{aligned}
P_{\rm phys}-I
 &=P^{1/2}(P_{\rm bl}-P^{-1})P^{1/2},\\
4I-P_{\rm phys}
 &=P^{1/2}(4P^{-1}-P_{\rm bl})P^{1/2}.
\end{aligned}
\]

Equations (6)--(7) prove (2), and (3) follows immediately.

## 4. Consequence for the all-grade gate

L214--L215 obtained the lower and upper endpoint inequalities only
after solving a coefficientwise rank-\(m\) Stein factorization.  L219
separates that work into two logically independent parts:

1. the full boundary-layer metric already lies strictly between the
   desired endpoints, by the elementary orbit sums (6)--(7);
2. only its elliptic Stein slack can fail positivity.

Therefore a future all-grade proof need not re-prove endpoint metric
positivity coefficient by coefficient.  It is enough to construct an
ordered correction \(X(c)\) such that

\[
\begin{aligned}
&(P_{\rm bl}+X)
 -\widehat T_c^*(P_{\rm bl}+X)\widehat T_c\succeq0,\\
&-(P_{\rm bl}-P^{-1})<X
  <4P^{-1}-P_{\rm bl}.                             \tag{8}
\end{aligned}
\]

The one-image frame is one candidate way to produce such an \(X\).
Equations (6)--(8) give an alternative target: control a Stein repair
directly inside two explicit positive orbit budgets.  This is weaker
than preserving a rank-\(m\) Stein defect and may admit a simpler
construction.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_boundary_metric_sandwich.py \
  --output \
  experiments/repeated_crabb_boundary_metric_sandwich_s70224.jsonl
```

The deterministic audit covers unstructured pure partial isometries
and fully delayed grades two through five.  It verifies both orbit
resolutions, both positive decompositions, and both physical metric
gaps.  It also records negative elliptic Stein eigenvalues, confirming
that the metric-only shortcut is false.  Equations (4)--(7), not the
floating audit, prove L219.  The tracked dataset SHA-256 is
`05088a6efb82245dd4a76ba9290ee0cd946b9e508696dbbe51721017d7e2b8b9`.
