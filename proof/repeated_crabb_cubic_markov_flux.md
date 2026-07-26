# The moving cubic defect has a gap-free Markov preimage

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

## 1. Result (L300, 2026-07-26)

Retain L299's pure partial-isometry colligation

\[
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\]

its first transfer coefficient

\[
B_1=W^*S^*V,
\]

and L298's grade-one metric/frame direction

\[
{\mathscr R}-S^*{\mathscr R}S
=V{\mathscr F}^*+{\mathscr F}V^*.                 \tag{1}
\]

L299 found the cubic defect \({\cal N}_3\) left when this
fixed-base direction is inserted into the moving elliptic pair.  Let

\[
Q_3=W^*{\cal G}_S({\cal N}_3)W.                  \tag{2}
\]

Then \(Q_3\) belongs to the full free-column/Markov endpoint-response
range at every colligation, including reducible and rank-changing
ones.  More precisely, there is a perpendicular state column
\(C_3\), \(V^*C_3=0\), such that

\[
\boxed{
W^*{\cal G}_S(VC_3^*+C_3V^*)W=Q_3,\qquad
\|C_3\|_F\le10\|{\mathscr R}\|_F.}                \tag{3}
\]

This is a pointwise bounded-energy selection theorem, not a
pseudoinverse formula.

Let \(\Gamma_S\) be the Frobenius operator norm of the stable Stein
inverse \({\cal G}_S\).  L298's explicit formula gives

\[
\|{\mathscr R}\|_F
\le2\Gamma_S\|B_1\|_F,                            \tag{4}
\]

and hence

\[
\boxed{\|C_3\|_F\le20\Gamma_S\|B_1\|_F.}          \tag{5}
\]

On every fixed repeated-Crabb neighbourhood, \(\Gamma_S\) is
uniformly bounded.  Thus the cubic correction does not acquire a
hidden Markov spectral-gap denominator as transfer rank changes.

Equivalently, for every Hermitian endpoint test \(Y\), let

\[
\begin{aligned}
H_Y-SH_YS^*&=WYW^*,\\
Z_Y&=(I-E)H_YV.
\end{aligned}                                     \tag{6}
\]

Then

\[
\boxed{
|\langle Y,Q_3\rangle|
\le20\|{\mathscr R}\|_F\|Z_Y\|_F
\le20\sqrt2\,\Gamma_S\|B_1\|_F
\sqrt{\langle Y,(I-\Phi\Phi^*)Y\rangle}.}         \tag{7}
\]

The second inequality uses L281's exact identity

\[
\langle Y,(I-\Phi\Phi^*)Y\rangle=2\|Z_Y\|_F^2.
\]

In particular, \(Q_3\) annihilates the entire L206 commutant
cokernel, not only the scalar trace direction.  L299's nonzero
compression to \(\ker B_1^*\) is therefore a genuine Markov
coboundary, not an obstruction.

L300 closes the **full cubic range and bounded-energy debt** created
by L299.  It does not yet prove the sharper ordered-flag estimate
needed after the cubic is canceled.  Inserting \(C_3\) into the
moving frame creates a new **quartic** cross with the raw first frame
coefficient; its quadratic self-cost appears later at degree six.
The quartic successor must be computed and charged to the next
surviving transfer row before any sextic conclusion is drawn.

## 2. Eliminate the moving frame

L299 writes

\[
{\cal N}_3
=-\{A_1^*{\mathscr R}S+S^*{\mathscr R}A_1
     +D_1{\mathscr F}^*+{\mathscr F}D_1^*\},      \tag{8}
\]

where

\[
A_1=(I+F)S^*(I+E)-S^3.
\]

Direct partial-isometry reduction gives the much shorter right and
left defect-frame tangents

\[
\boxed{D_1=-2(S^*)^2V,\qquad L_1=-2S^2W.}         \tag{9}
\]

The following special identity of L298's grade-one frame is
load-bearing:

\[
\boxed{(S^*)^2{\mathscr F}=0.}                    \tag{9a}
\]

Indeed, L298 gives

\[
{\mathscr F}
=-\frac12VB_1^*B_1-\frac12(I-E)SWB_1.
\]

Since \(VB_1^*=ESW\), \(E+(I-E)=I\), and
\((S^*)^2SW=S^*(S^*S)W=S^*W=0\), equation (9a)
follows.  This is not a generic consequence of \(S^*W=0\) for an
arbitrary correction frame.

Now use (1) in the frame cross in (8).  Equation (9a) gives

\[
(S^*)^2(V{\mathscr F}^*+{\mathscr F}V^*)
=(S^*)^2({\mathscr R}-S^*{\mathscr R}S)
=(S^*)^2{\mathscr R}-(S^*)^3{\mathscr R}S.
\]

Adding the adjoint and expanding \(A_1\) gives

\[
\boxed{
\begin{aligned}
{\cal N}_3={}&
2\{(S^*)^2{\mathscr R}+{\mathscr R}S^2\}\\
&-\{(S^*)^3{\mathscr R}S+S^*{\mathscr R}S^3\}\\
&-(I+E)S(I+F){\mathscr R}S\\
&-S^*{\mathscr R}(I+F)S^*(I+E).
\end{aligned}}                                   \tag{10}
\]

Thus the moving frame disappears completely.  The remaining defect
is one universal linear expression in the already bounded L298
metric direction.

## 3. The scalar coboundary is identically zero

Stein adjointness gives

\[
\langle Y,Q_3\rangle=\operatorname {tr}(H_Y{\cal N}_3).       \tag{11}
\]

Cyclically move \({\mathscr R}\) to the left in (11).  The
coefficient of a scalar test \(H_Y=I\) is

\[
\begin{aligned}
K_I={}&2(S^*)^2+2S^2-S(S^*)^3-S^3S^*\\
&-S(I+E)S(I+F)-(I+F)S^*(I+E)S^*.                 \tag{12}
\end{aligned}
\]

This vanishes as an operator.  Indeed

\[
\begin{aligned}
S(I+E)S(I+F)&=2S^2-S^3S^*,\\
(I+F)S^*(I+E)S^*&=2(S^*)^2-S(S^*)^3,
\end{aligned}
\]

using only \(E=I-S^*S\), \(F=I-SS^*\), and
\(SS^*S=S\).

Equation (12) proves \(\operatorname {tr}Q_3=0\), but the
commutator calculation below proves the stronger full-cokernel
compatibility.

## 4. Every dual pairing is controlled by the Markov defect

L208 gives

\[
[S,H_Y]=SZ_YV^*.                                  \tag{13}
\]

Put \(C_Y=SZ_YV^*\).  Since \(H_Y\) commutes with \(F\),

\[
\begin{aligned}
[E,H_Y]&=VZ_Y^*-Z_YV^*,\\
[S^*,H_Y]&=-VZ_Y^*S^*.
\end{aligned}                                     \tag{14}
\]

Subtract \(H_YK_I=0\) from the cyclic coefficient in (11), and use
(13)--(14).  One obtains

\[
\boxed{\langle Y,Q_3\rangle
=\operatorname {tr}({\mathscr R}K_Y),}            \tag{15}
\]

where

\[
\begin{aligned}
K_Y={}&2(SC_Y+C_YS)-C_Y(S^*)^3\\
&-(S^2C_Y+SC_YS+C_YS^2)S^*\\
&-C_Y(I+E)S(I+F)\\
&-(I+F)\bigl[
 S^*(VZ_Y^*-Z_YV^*)
 -VZ_Y^*S^*(I+E)
\bigr]S^* .
\end{aligned}                                     \tag{16}
\]

Every term now contains \(Z_Y\).  Since all projections and \(S\)
are contractions,

\[
\begin{aligned}
\|K_Y\|_F
&\le(4+1+3+4+8)\|Z_Y\|_F\\
&=20\|Z_Y\|_F.                                    \tag{17}
\end{aligned}
\]

Cauchy--Schwarz in (15) proves the first inequality in (7).
If \(Z_Y=0\), then the pairing is zero.  By L206 this is exactly
orthogonality to the full endpoint cokernel, so finite-dimensional
Fredholm duality proves existence in (3).

For the norm, the adjoint of the balanced response map

\[
{\cal M}_S(C)=W^*{\cal G}_S(VC^*+CV^*)W
\]

is \(Y\mapsto2Z_Y\).  Inequality (7) in its first form says

\[
|\langle Y,Q_3\rangle|
\le10\|{\mathscr R}\|_F\|{\cal M}_S^*Y\|_F.
\]

The Hilbert-space quotient/Riesz theorem therefore supplies a
preimage with norm at most \(10\|{\mathscr R}\|_F\), proving (3)
without selecting a singular subspace.

## 5. Uniform local transfer bound

L298 has

\[
\begin{aligned}
{\mathscr R}
&=-{\cal G}_S(VB_1^*B_1V^*)
  +{\cal G}_S(VD^*+DV^*),\\
D&=-\frac12(I-E)SWB_1.
\end{aligned}                                     \tag{18}
\]

Because \(\|B_1\|\le1\),

\[
\|B_1^*B_1\|_F\le\|B_1\|_F,\qquad
2\|D\|_F\le\|B_1\|_F.
\]

Applying \(\|{\cal G}_S\|_{F\to F}=\Gamma_S\) to both terms in
(18) proves (4).  Stability makes \(\Gamma_S\) finite, and
continuity makes it uniformly bounded after shrinking any fixed
repeated-block neighbourhood.

This is the correct scope: no global dimension-free Stein bound is
claimed for arbitrary pure partial isometries approaching the unit
circle.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_cubic_markov_flux.py \
  --output \
  experiments/repeated_crabb_cubic_markov_flux_s70226.jsonl
```

The checker verifies (9)--(16) exactly on L299's rational rank-one
flag example.  It then audits unstructured, rank-chain, and reducible
two-block colligations, including a nonscalar commutant cokernel.  It
independently synthesizes the minimum-norm response column and checks
both bounds in (3) and (7).  All 13 records pass.  The tracked
dataset has SHA-256

```text
1c8bf694c4f0783ce9d3b38da44e046ad3a5b9a406888c1a60ebe95db8eb064d
```
