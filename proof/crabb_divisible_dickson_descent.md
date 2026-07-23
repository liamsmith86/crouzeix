# Exact rank-one metric lift for divisible Crabb grades (2026-07-23)

## 1. Result

Let `L=qk` with `q>=2`.  In size `L+1`, take the phase-one equality
direction supported at coefficient grades `k` and `L-k`, and write

\[
T_{L,k}(a,c)
=\phi_c\{S_{L,k}(a,c)\}.
\]

Let `T_{q,1}(a,r)` denote the same equality/ellipse construction in
size `q+1`, now at first offset, with `r=c^k`, and let
`\Gamma_(q+1)` be L118's locally optimized rank-one Stein envelope
(optimized within that analytic rank-one branch).
Then, for every fixed pair `(q,k)`,

\[
\boxed{
t_*(T_{q,1}(a,c^k))
\leq t_*(T_{L,k}(a,c))
\leq\Gamma_{q+1}(T_{q,1}(a,c^k))
}
\qquad (|a|+c\ \hbox{sufficiently small}).            \tag{1}
\]

The neighbourhood may depend on `(q,k)`.  Both comparisons in (1) are
exact, not only amplitude-Hessian inequalities.  L130 is the special
case `q=2`, where the outer size-three rank-one branch is known to be
optimal and the sandwich collapses to equality.

Consequently, once the offset-one coefficient is computed in every
size, every grade dividing `L` follows automatically.  In particular
an upper face `-64a^2c^(2k)` for a divisible grade is the offset-one
rank-one face after the substitution `r=c^k`.

It would be incorrect to claim equality of the two `t_*` values in
general from the present argument.  L118 proves that its rank-one
branch is feasible and has the correct second variation, but it does
not prove that branch is the globally optimal SDP metric at every
larger-size off-axis point.

## 2. Exact polynomial and coordinate reduction

Let

\[
P_0=2,\qquad P_1=z,\qquad P_m=zP_{m-1}-cP_{m-2}
\]

and put

\[
V=\operatorname{span}\{e_0,e_k,e_{2k},\ldots,e_{qk}\}.
\]

Dickson path cancellation gives

\[
\boxed{
P_k(S_{L,k}(a,c))
\simeq S_{q,1}(a,c^k)\oplus S_{\rm in}(a,c),}         \tag{2}
\]

where `V` and `V^\perp` reduce the displayed matrix and the first
summand acts on `V`.  The symbol `simeq` in (2) is only the coordinate
ordering; the block equality itself is exact.  Moreover,

\[
\boxed{
\Pi_VP_m(S_{L,k}(a,c))\Pi_V=0,\qquad1\leq m<k.}       \tag{3}
\]

Both identities include every amplitude power.  In fact `P_k` is
exactly affine in `a` on this pencil.

Here is a direct recurrence proof.  The finite Dickson model used in
L126 works with terminal length `L`, not only `L=2k`: take the basis
formed by the normalized `P_j`, impose
`P_(L+1)=cP_(L-1)`, and use

\[
 P_mP_n=P_{m+n}+c^nP_{m-n}\quad(m>n),\qquad
 P_n^2=P_{2n}+2c^n.
\]

Multiplication by `z` is `C+cJCJ`.  Multiplication by `P_k` changes a
basis grade only by `+k` or `-k`, with the same terminal folding.
It therefore preserves every residue class modulo `k`; on residue
zero the product law gives exactly
`C_(q+1)+c^kJC_(q+1)J`.  The equality tangent is

\[
2(e_0-e_L)
\{e_{k+1}^*+e_{L-k+1}^*
-c e_{k-1}^*-c e_{L-k-1}^*\}.                       \tag{4}
\]

In the differentiated Dickson recurrence, the two terms at each end
of (4) pair the first reversed path exactly as in L126 and L131.
They leave the residue-zero columns of the size-`q+1` first-offset
tangent.  Every cross-residue path has an unmatched distance strictly
below `k` and cancels.  More explicitly, between two amplitude
insertions one encounters a factor
`v^*(C+cJCJ)^j(e_0-e_L)`, where `v` is the row vector in (4);
this is zero for `0<=j<=k-2`.  A degree-`k` Dickson word containing
two insertions has at most `k-2` intervening matrix factors, so every
such word vanishes.  This proves exact affinity and (2).  The same
residue/distance argument before the first possible return proves
(3).

L123's Toeplitz coordinate Gramian reduces `V` as well, and

\[
\boxed{
K_{L,k}(a)|_V=K_{q,1}(a).}                            \tag{5}
\]

This follows entrywise: the only nonzero Toeplitz offsets are `k` and
`(q-1)k`, which become offsets one and `q-1` after division by `k`.
Thus (2) is an orthogonal physical reduction, not merely a coordinate
invariant subspace.

## 3. Descended inner blocks are inactive

At `(a,c)=(0,0)`, the residue-zero block in (2) is the size-`q+1`
Crabb operator and has similarity square four.  Every other residue
class has `q` coordinates and carries the unweighted nilpotent shift

\[
\begin{pmatrix}
0&1&&\\ &0&\ddots&\\ &&\ddots&1\\ &&&0
\end{pmatrix},
\]

which is already a contraction.  The coordinate Gramian is the
identity on those classes.  Strict contraction feasibility persists
for fixed `(q,k)` under a sufficiently small perturbation, while the
outer value remains near four.  Hence

\[
t_*\!\left(B_{k,c}(T_{L,k})\right)
=t_*(T_{q,1}(a,c^k))                                 \tag{6}
\]

locally, where

\[
B_{k,c}\circ\phi_c=\phi_{c^k}\circ P_k.
\]

## 4. The metric lift has no condition cost

L130's metric-lift proof does not use that the active outer block has
dimension three.  For the degree-`k` finite Blaschke product, factor
the outside critical points and put `F=Q/D` as there.  Its residue
theorem gives, for every `f` in the model space,

\[
{1\over k}\sum_{B(t)=u}{f(t)\over F(t)}
={\langle f,F\rangle\over\|F\|^2}.                   \tag{7}
\]

Equations (3) and polynomial division transfer (7) to the
`(q+1)`-dimensional physical outer space:

\[
O^*(f/F)(T)O
={\langle f,F\rangle\over\|F\|^2}I_{q+1}.            \tag{8}
\]

Now take the defect `v` selected by L118's locally optimized
rank-one-envelope branch, define `J=F(T)^{-*}O`, and set `d=Jv`.
The identical reconstruction calculation from L130 gives
`sum gamma_jw_j=0`, so the lifted full Stein metric reduces the outer
space and its compression is a scalar copy of the outer rank-one
envelope metric.

At the apex, that full metric is

\[
\operatorname{diag}(1,2,\ldots,2,4).
\]

The outer residue-zero levels include the endpoints one and four;
every omitted level equals two.  Continuity therefore keeps all inner
levels strictly between the active outer endpoints.  This proves the
right inequality in (1), with the outer rank-one envelope.  The
finite-Blaschke dual/lower transfer and (6) give the left inequality.

## 5. Scope and next step

This is a genuine enlargement of L130, but it does not yet cover
`k` when `k` does not divide `L`.  For a nondivisible grade the two
endpoint residue chains have unequal lengths and `P_k` has cross
blocks.  L131 shows that their only first-face invariant is still the
reflected coefficient `u_kc^k`; the remaining task is an
unequal-residue Schur localization, not another divisible descent.

The other immediate debt is a written all-size calculation of the
offset-one coefficient.  Its target order is only `c^2`, and A85's
exact records already show the dimension-stable endpoint pair
`(48,128)`, hence `128-4(48)=-64`.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_divisible_dickson_descent.py \
  --output experiments/crabb_divisible_dickson_descent_s70223.jsonl
```

The checker verifies (2)--(5), exact amplitude linearity, and the
inactive apex shifts on the complete grid `1<=k<=6`, `2<=q<=6`.
The independent floating checker constructs the critical-factor lift
itself and verifies scalar compression, outer/inner metric reduction,
preservation of the outer rank-one condition, and strict inner
endpoint margins:

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_divisible_metric_lift.py \
  --output experiments/crabb_divisible_metric_lift_s70223.jsonl
```
