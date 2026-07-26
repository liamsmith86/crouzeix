# Every first moving successor is a gap-free Markov response

## 1. Result (L301, 2026-07-26)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\]

and let \(X=X^*\) and \(C\) be an arbitrary fixed-base Stein
direction:

\[
\boxed{X-S^*XS=VC^*+CV^*.}                       \tag{1}
\]

No endpoint condition on \(X\) or perpendicularity condition on
\(C\) is assumed.

Let the first moving elliptic operator/frame tangents be

\[
\begin{aligned}
A_1&=(I+F)S^*(I+E)-S^3,\\
D_1&=-2(S^*)^2V.
\end{aligned}                                     \tag{2}
\]

The first successor forcing produced by inserting \((X,C)\) into
the moving pair is

\[
\boxed{
{\cal T}(X,C)
=-\{A_1^*XS+S^*XA_1+D_1C^*+CD_1^*\}.}            \tag{3}
\]

Its balanced upper endpoint

\[
\tau(X,C)=W^*{\cal G}_S({\cal T}(X,C))W           \tag{4}
\]

is always in the complete L204/L280 free-column response range.
More precisely, there is a perpendicular column \(C^+\) such that

\[
\boxed{
{\cal M}_S(C^+)=\tau(X,C),\qquad
\|C^+\|_F
\le10\|X\|_F+2\|C\|_F.}                           \tag{5}
\]

Here

\[
{\cal M}_S(K)
=W^*{\cal G}_S(VK^*+KV^*)W.
\]

Thus transport through one moving order creates no scalar,
commutant, Markov-gap, or rank-selection obstruction.  The result is
uniform and uses no pseudoinverse.

For L300's cubic repair this has an immediate consequence.  If two
bounded cubic columns cancel the same \(Q_3\), their difference is an
endpoint-null fixed-base direction.  L301 says their quartic
successors differ only by another bounded Markov response.
Therefore:

\[
\boxed{\text{the quartic commutant/separator pairing is independent
of the chosen cubic preimage.}}                   \tag{6}
\]

This removes the need to optimize the cubic and quartic selections
simultaneously.  Choose any bounded L300 cubic preimage; at quartic
order retain the selection-independent commutant component and use
L301/L280 to redistribute the response component.

L301 does not prove that the retained quartic component has the
favorable sign required by the metric sandwich.  Computing that
physical invariant and bounding the new quartic selection on ordered
flags are the next gates.

## 2. Exact successor decomposition

L300 introduced the universal frame-eliminated expression

\[
\begin{aligned}
{\cal L}(X)={}&
2\{(S^*)^2X+XS^2\}\\
&-\{(S^*)^3XS+S^*XS^3\}\\
&-(I+E)S(I+F)XS\\
&-S^*X(I+F)S^*(I+E).                             \tag{7}
\end{aligned}
\]

For L298's special frame, \((S^*)^2C=0\), and L300 had
\({\cal T}={\cal L}\).  For an arbitrary direction (1), use

\[
(S^*)^2VC^*
=(S^*)^2(X-S^*XS)-(S^*)^2CV^*.
\]

Adding the adjoint gives the exact general formula

\[
\boxed{
{\cal T}(X,C)
={\cal L}(X)
-2\{(S^*)^2CV^*+VC^*S^2\}.}                      \tag{8}
\]

This identifies precisely why L300 did not automatically iterate:
the displayed remainder is absent for L298 but is present for a
generic selected response frame.

## 3. Dual pairing and gap-free bound

For a Hermitian endpoint test \(Y\), put

\[
\begin{aligned}
H_Y-SH_YS^*&=WYW^*,\\
Z_Y&=(I-E)H_YV.
\end{aligned}                                     \tag{9}
\]

L300 proves

\[
\langle Y,W^*{\cal G}_S({\cal L}(X))W\rangle
=\operatorname {tr}(XK_Y),\qquad
\|K_Y\|_F\le20\|Z_Y\|_F.                          \tag{10}
\]

It remains to pair the two-term remainder in (8).  L208's iterated
commutator identity gives

\[
H_Y(S^*)^2
=(S^*)^2H_Y+VZ_Y^*(S^*)^2+S^*VZ_Y^*S^*.
\]

Since \(V^*(S^*)^2=0\) and \(V^*S^*=0\),

\[
\boxed{V^*H_Y(S^*)^2=Z_Y^*(S^*)^2.}              \tag{11}
\]

Stein adjointness, (8), and (11) now give the exact pairing

\[
\boxed{
\langle Y,\tau(X,C)\rangle
=\operatorname {tr}(XK_Y)
-4\operatorname {Re}\operatorname {tr}
 Z_Y^*(S^*)^2C.}                                 \tag{12}
\]

Consequently

\[
\boxed{
|\langle Y,\tau(X,C)\rangle|
\le
(20\|X\|_F+4\|C\|_F)\|Z_Y\|_F.}                 \tag{13}
\]

If \(Z_Y=0\), the pairing vanishes.  L206 identifies these tests with
the entire endpoint cokernel, so (13) proves
\(\tau(X,C)\in\operatorname {ran}{\cal M}_S\), including on
reducible strata.

The adjoint of \({\cal M}_S\), on perpendicular columns, is

\[
{\cal M}_S^*Y=2Z_Y.
\]

Equation (13) is therefore

\[
|\langle Y,\tau(X,C)\rangle|
\le
(10\|X\|_F+2\|C\|_F)\|{\cal M}_S^*Y\|_F.
\]

The finite-dimensional Hilbert quotient/Riesz theorem produces a
preimage with the bound (5), without choosing a singular subspace.

L281 also turns (13) into the Markov-energy form

\[
|\langle Y,\tau(X,C)\rangle|
\le
(10\sqrt2\|X\|_F+2\sqrt2\|C\|_F)
\sqrt{\langle Y,(I-\Phi\Phi^*)Y\rangle}.          \tag{14}
\]

## 4. Trace and selection invariance

Taking \(Y=I\) in (12) gives \(Z_I=0\), hence

\[
\boxed{\operatorname {tr}\tau(X,C)=0.}            \tag{15}
\]

There is also a direct Green identity behind (15).  The left defect
frame tangent is \(L_1=-2S^2W\).  Pairing the left slack tangent with
\(X\) reduces the operator part of (3) to

\[
2\operatorname {Re}\operatorname {tr}(L_1^*XW).
\]

Multiplying (1) by \((S^*)^2\) and taking the trace gives

\[
\operatorname {tr}\{W^*(S^*)^2XW\}
=\operatorname {tr}\{C^*(S^*)^2V\},
\]

which cancels the frame part exactly.

For (6), let \((X,C)\) and \((\widetilde X,\widetilde C)\) be two
cubic solutions with the same endpoint.  Their difference obeys
(1), and its fixed-base endpoint response is zero.  The difference
of their first successors is exactly

\[
\tau(\widetilde X-X,\widetilde C-C).
\]

By (5) this lies in the response range with a bound depending only on
the difference direction.  Every L206 commutant pairing is therefore
the same for the two quartic successors.

## 5. Scope in the all-series recurrence

L285 says each prepared coefficient creates one new endpoint
homology equation at the next order.  L301 closes the part of that
equation produced by the **first** moving tangent:

\[
(X_n,C_n)
\longmapsto
{\cal T}(X_n,C_n)
\]

always has a bounded-energy response.

Higher successors also contain:

1. higher fixed operator/frame coefficients \(A_j,D_j\);
2. products with earlier metric directions; and
3. quadratic frame terms.

L301 does not claim those pieces are responses or favorable.
However, it removes the universal nearest-neighbour transport from
the arbitrary-grade debt.  The retained even invariant must come
from the remaining physical terms, not from a choice of cubic
response representative.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_successor_markov_transport.py \
  --output \
  experiments/repeated_crabb_successor_markov_transport_s70226.jsonl
```

The checker verifies (8), (11), and (12) exactly on a rational
noncommuting direction.  It then audits unstructured and reducible
colligations and repeated monomial apices, where the response map and
successor endpoint both collapse.  It independently synthesizes the
minimum-norm response column and checks (5), (13), and L281's energy
identity.  All 10 records pass.  The tracked dataset has SHA-256

```text
24ca5eade9348ee61821c24604abd8dcbb8ad1156e62cb021f89e9e193343d7c
```
