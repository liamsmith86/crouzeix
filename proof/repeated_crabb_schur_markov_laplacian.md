# The first Markov gap is the Schur-parameter commutator Laplacian

> **Route update.**  L294 subsequently proves the exact nonlinear
> numerator/Dirichlet Cauchy--Schwarz estimate for L212's canonical
> channel-minus-Gram target.  L295 then proves that the complete
> raw L227/L285 repair upper face cannot be compared with that
> lower-tight target: their apex values are respectively `-12I` and
> `-16I`.  L297 identifies the gap as lower retightening, so the
> direct comparison is not obstructed.  L298 subsequently bypasses
> that comparison and closes the desired complete-delay lower-tight
> face by an exact orbit-Gram/L212 cancellation.  Its moving-series
> defect remains open.  If the raw favorable lower budget is retained,
> pair only its mixed transport remainder; L296 handles arbitrary
> polarized responses once their column energy is controlled.

## 1. Result (L293, 2026-07-26)

Retain L218's matrix-Schur chart for a square inner transfer near the
repeated monomial

\[
B_0(z)=z^L U.
\]

Right-normalize the transfer by its terminal unitary and hence take
\(U=I\).  This loses no generality even when the terminal unitary
varies: if \(\widetilde B=B U_B^*\), then
\(\widetilde\Phi\widetilde\Phi^*=\Phi\Phi^*\) pointwise.  Let all
nonterminal Schur parameters vary along one real parameter:

\[
\Gamma_j(\varepsilon)=\varepsilon\Delta_j+O(\varepsilon^2),
\qquad 1\le j<L.                                  \tag{1}
\]

Let

\[
\Phi_\varepsilon(H)
=\sum_{n\ge1}B_n(\varepsilon)H B_n(\varepsilon)^*,
\qquad
K_\varepsilon=I-\Phi_\varepsilon\Phi_\varepsilon^*              \tag{2}
\]

be L280's bistochastic transfer channel and Markov Laplacian.  Then

\[
\boxed{
K_\varepsilon(H)
=\varepsilon^2\sum_{j=1}^{L-1}
\left(
 [\Delta_j^*,[\Delta_j,H]]
 +[\Delta_j,[\Delta_j^*,H]]
\right)
+O(\varepsilon^3).}                               \tag{3}
\]

For every Hermitian \(H\), its Dirichlet form is therefore

\[
\boxed{
\langle H,K_\varepsilon H\rangle_{\rm HS}
=2\varepsilon^2\sum_{j=1}^{L-1}
\|[H,\Delta_j]\|_F^2+O(\varepsilon^3).}            \tag{4}
\]

In particular, the quadratic Markov gap is positive semidefinite and
its kernel is exactly the common self-adjoint commutant of the first
Schur jets:

\[
\boxed{
K^{(2)}_\Delta(H)=0
\quad\Longleftrightarrow\quad
[H,\Delta_j]=0\quad(1\le j<L).}                   \tag{5}
\]

Along an arbitrary real-analytic arc, first make this analytic
terminal normalization and let

\[
q=\min_j\operatorname{ord}_s\Gamma_j(s),\qquad
\Delta_j=[s^q]\Gamma_j(s),
\]

with \(\Delta_j=0\) for parameters of higher order.  Then (3)--(5)
hold with \(\varepsilon=s^q\).  L281's physical observability column
has the matching first size

\[
\boxed{
\|2R_H(s)\|_F^2
=4s^{2q}\sum_j\|[H,\Delta_j]\|_F^2
+O(s^{2q+1}).}                                    \tag{6}
\]

This is an arbitrary-length, all-Schur-grade result.  It removes one
possible failure mode from L282: a closing Markov gap cannot have an
uncontrolled first exponent.  Its leading square is exactly the
simultaneous commutator defect of the first nonzero Schur jets.

L293 alone does **not** prove L282's physical flux inequality.  L294
subsequently proves it for every transfer flux \(G_k\).  L295 shows
that one proposed raw-to-lower-tight comparison was false, while
L297 identifies its four-unit gap as lower retightening.  The
remaining task is now sharper:

> either retain L283's favorable raw lower budget, separate its
> \(-12B_kB_k^*\) upper face, and control the remaining mixed
> L285/L289 transport; or prove the lower-tight L212 physical base in
> arbitrary grade.

L294 then gives L282's desired gap-free bound for such a flux.  L293
also does not replace L292's endpoint-valuation route; it gives the
exact leading geometry of its Markov alternative.

## 2. First transfer jet in Schur coordinates

L218 proves, after the terminal normalization above,

\[
[\,\varepsilon\,]B(z)
=\sum_{j=1}^{L-1}
\left(z^j\Delta_j-z^{2L-j}\Delta_j^*\right).       \tag{7}
\]

The \(2(L-1)\) displayed degrees are pairwise distinct and none is
\(L\).  Write their coefficient matrices as

\[
C_j=\Delta_j,\qquad C_{2L-j}=-\Delta_j^*.          \tag{8}
\]

Thus

\[
B_L(\varepsilon)=I+\varepsilon^2A+O(\varepsilon^3),
\qquad
B_n(\varepsilon)=\varepsilon C_n+O(\varepsilon^2)
\quad(n\ne L).                                    \tag{9}
\]

Matrix innerness gives both Parseval identities

\[
\sum_nB_nB_n^*=I,\qquad
\sum_nB_n^*B_n=I.                                 \tag{10}
\]

At order two, either identity gives

\[
\boxed{
A+A^*=-R,\qquad
R=\sum_j(\Delta_j\Delta_j^*
          +\Delta_j^*\Delta_j).}                  \tag{11}
\]

The skew-Hermitian part of \(A\) is not determined by Parseval.  It
will cancel from \(K_\varepsilon\), so no second-order Schur-recursion
formula is needed.

## 3. Channel expansion

Equations (8)--(11) give

\[
\Phi_\varepsilon
=I+\varepsilon^2\Psi+O(\varepsilon^3),             \tag{12}
\]

where

\[
\Psi(H)=AH+HA^*
+\sum_j\{\Delta_jH\Delta_j^*
          +\Delta_j^*H\Delta_j\}.                 \tag{13}
\]

Consequently

\[
K_\varepsilon
=-\varepsilon^2(\Psi+\Psi^*)+O(\varepsilon^3).     \tag{14}
\]

Put

\[
T(H)=\sum_j(\Delta_jH\Delta_j^*
             +\Delta_j^*H\Delta_j).
\]

Using (11) in (13)--(14) cancels the skew part of \(A\) and gives

\[
\begin{aligned}
K^{(2)}_\Delta(H)
&=RH+HR-2T(H)\\
&=\sum_j\left(
 [\Delta_j^*,[\Delta_j,H]]
 +[\Delta_j,[\Delta_j^*,H]]
\right),
\end{aligned}                                      \tag{15}
\]

which proves (3).

## 4. Exact Dirichlet square

For Hermitian \(H\), cyclicity gives

\[
\begin{aligned}
\left\langle H,K^{(2)}_\Delta(H)\right\rangle
&=\sum_j\left(
\|[H,\Delta_j]\|_F^2
+\|[H,\Delta_j^*]\|_F^2
\right)\\
&=2\sum_j\|[H,\Delta_j]\|_F^2.                    \tag{16}
\end{aligned}
\]

The last equality uses
\([H,\Delta_j^*]=-[H,\Delta_j]^*\).
This proves (4).  Equality in (16) holds exactly when all the
commutators vanish, proving (5).

For an analytic arc, discard parameters whose order is above the
minimum \(q\) and substitute \(\varepsilon=s^q\).  Terms involving a
higher jet have strictly larger order, proving the arc statement.
Finally L281 gives

\[
\|2R_H\|_F^2=2\langle H,K_\varepsilon H\rangle,
\]

which combined with (4) proves (6).

## 5. Consequence for the live flux gate

The denominator in L282's exact dual ratio is

\[
\langle\widehat Y,K_\varepsilon\widehat Y\rangle.
\]

At the first nonzero Schur order, L293 identifies it with

\[
2\varepsilon^2
\sum_j\|[\widehat Y,\Delta_j]\|_F^2.               \tag{17}
\]

There are therefore only two cases on a hypothetical failure arc.

1. Some first Schur jet fails to commute with the positive test.
   Then (17) supplies the exact quadratic denominator needed for a
   Cauchy--Schwarz flux estimate.
2. Every first Schur jet commutes with the test.  Then the test
   reduces the first perturbed copy decomposition, the quadratic
   denominator vanishes, and the analysis advances to the next
   Schur order.  At a fully reducing first active block, L279 gives
   the strict negative endpoint trace.

L294 subsequently supplies the exact numerator identity for every
\(G_k\)-flux, at full amplitude rather than only in this first jet.
L295 proves that the raw L285/L289 upper endpoint does **not** lie
below L212's lower-tight target: at the apex the former is \(-12I\)
and the latter is \(-16I\).  L297 shows that lower retightening
changes the former to \(-16I\), so it leaves the effective
partial-flag comparison open.  In the raw normalization, what is
still missing is the physical remainder identity after the favorable
L283 Gram is separated.  Merely invoking compactness would not prove
that identity, and A179 forbids raw all-series superposition.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_schur_markov_laplacian.py \
  --output \
  experiments/repeated_crabb_schur_markov_laplacian_s70226.jsonl
```

The exact SymPy checker uses several lengths, multiplicities, and
simultaneously active noncommuting Schur jets.  It verifies:

1. the two order-two Parseval identities;
2. cancellation of an arbitrary skew-Hermitian terminal jet;
3. the superoperator identity (15) on a full Hermitian basis;
4. the Dirichlet identity (16); and
5. the pairwise-disjoint tangent degrees in (7).

Every residual is exactly zero; no floating tolerance is used.
The tracked dataset has SHA-256

```text
2e5c20fabc9fed72a752a802cee12374e751327585144afb04883b93f6a17059
```
