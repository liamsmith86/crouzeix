# Exact analytic metric chart near repeated `C3` (2026-07-22)

## 1. Statement

Put

\[
 A_0=I_m\otimes C_3,\qquad
 M=I_m\otimes\operatorname{diag}(1,2,4),
\]

and order the space by the three physical levels, each with copy
multiplicity `m`.  There are neighbourhoods of `A_0`, zero, and `M` with
the following property.

For every operator `T` near `A_0`, every free complex block

\[
 B=(P_{01}\ P_{02})\in M_{m,2m}
\]

near zero, and every Hermitian Stein-slack parameter `H` near zero, there
is a unique Hermitian range block

\[
 C=\begin{bmatrix}P_{11}&P_{12}\\P_{21}&P_{22}\end{bmatrix}
\]

near `diag(2I,4I)` such that the lower metric constraint has zero Schur
complement and the Stein constraint has Schur complement `H` on its base
kernel.  The resulting metric `P=P(T,B,H)` is real analytic in `(T,B,H)`.

Moreover, `P>=I` holds automatically, and `P-T^*PT>=0` holds whenever
`H>=0`.  The remaining condition `P<=4I` is equivalent to one explicit
analytic `m x m` endpoint inequality.

## 2. Tighten the lower constraint exactly

Given `B` and `C`, define

\[
 P_{00}=I+B(C-I)^{-1}B^*,\qquad
 P=\begin{bmatrix}P_{00}&B\\B^*&C\end{bmatrix}.       \tag{1}
\]

Because `C-I` is positive near `diag(I,3I)`, the Schur complement gives

\[
 P-I\succeq0,\qquad
 (P-I)/(C-I)=0.                                       \tag{2}
\]

Thus (1) is the exact nonlinear version of the lower-metric recursions used
in L76--L104.  The whole level-zero/range cross block `B` remains free.

## 3. Tighten the Stein constraint by the implicit-function theorem

Let

\[
 S(T,B,C)=P-T^*PT
\]

and split its physical levels as range `R={0}` and kernel `K={1,2}`.  Near
the base, `S_RR` is positive definite.  For a Hermitian parameter `H`,
consider the `2m x 2m` equation

\[
 {\cal F}(T,B,C)
 =S_{KK}-S_{KR}S_{RR}^{-1}S_{RK}=H.                  \tag{3}
\]

At `(T,B,C,H)=(A_0,0,diag(2I,4I),0)`, equation (3) holds.  If

\[
 \dot C=\begin{bmatrix}X_{11}&X_{12}\\
 X_{12}^*&X_{22}\end{bmatrix},
\]

direct multiplication by the two weight-`sqrt(2)` links of `C_3` gives

\[
 D_C{\cal F}(\dot C)
 =\begin{bmatrix}
 X_{11}&X_{12}\\
 X_{12}^*&X_{22}-2X_{11}
 \end{bmatrix}.                                      \tag{4}
\]

This real-linear map on the Hermitian matrices is invertible, with inverse

\[
 (Y_{11},Y_{12},Y_{22})
 \longmapsto(Y_{11},Y_{12},Y_{22}+2Y_{11}).           \tag{5}
\]

The finite-dimensional real-analytic implicit-function theorem therefore
gives a unique analytic solution `C=C(T,B,H)` to (3).

Since `S_RR` remains positive, the Schur complement in the other direction
now gives

\[
 S(T,B,C(T,B,H))\succeq0
 \quad\Longleftrightarrow\quad H\succeq0.             \tag{6}
\]

Equations (1)--(6) prove the asserted analytic parameterization of every
lower-tight certificate with prescribed nearby Stein Schur complement.

## 4. One exact upper endpoint

For `Q=4I-P`, split off physical level two as its base kernel and let
`L={0,1}`.  The range block `Q_LL` is positive near `diag(3I,2I)`.  Hence

\[
 P\preceq4I
\quad\Longleftrightarrow\quad
 {\cal E}(T,B,H)\preceq0,                             \tag{7}
\]

where

\[
 \boxed{
 {\cal E}(T,B,H)
 =P_{22}-4I+
 P_{2L}(4I-P_{LL})^{-1}P_{L2}.}                       \tag{8}
\]

Thus the repeated-`C3` complete-similarity problem is locally reduced to
choosing `B` and `H>=0` for which the single analytic endpoint (8) is
nonpositive.  No infinite forced-metric recursion remains: for any analytic
choice of `T,B,H`, all the other metric coefficients are the convergent
Taylor series of `C(T,B,H)`.

L104's homogeneous block is exactly a Taylor coefficient of `B`; its
propagated `P_{12}=2P_{01}` is the linearization (5).

## 5. Consequence and remaining gate

This closes one part of the all-orders debt: divergence of the forced lower
or Stein metric coefficients cannot occur.  It does **not** yet prove that
there is a convergent choice of the free block `B` making (8) nonpositive
on the terminal tube.

The unresolved scalar issue is now isolated from metric feasibility:

1. use the off-diagonal freedom in `B` and L104's surjective mixed
   derivative to center the traceless part of (8) for `d!=0`;
2. prove that the resulting scalar endpoint retains L100's negative
   one-sided Jensen slope uniformly as `d->0`;
3. merge with L99 when the transverse and normal scales are comparable.

The only nonsmooth object in this remaining step is the map
`A -> T=phi_A(A)` at repeated support branches.  The exact metric chart
itself is analytic in the matrix `T`.

## 6. Slack-compatible frozen-normal bridge

The parameter `H` resolves the compatibility issue between this chart and
L106.  Let `P_N=diag(P_j)` be L106's exact metric for
`T_N=f(N)`.  Scale each single-block metric so that its smallest eigenvalue
is one.  Near the Crabb metric that eigenvalue is simple, while the
physical-level range block of `P_j-I` stays positive.  Since `P_j-I` is
singular positive semidefinite, its level-zero Schur complement is exactly
zero.  The same is true for their direct sum.

Let `B_N` be the level-zero/range block of `P_N`, and let

\[
 H_N=S_{KK}-S_{KR}S_{RR}^{-1}S_{RK}\succeq0           \tag{9}
\]

be the Stein Schur complement of `S=P_N-T_N^*P_NT_N`.
Then `P_N` solves (1), (3) with data `(T_N,B_N,H_N)`.  Local uniqueness gives

\[
 \boxed{P_N=P(T_N,B_N,H_N).}                          \tag{10}
\]

Thus no preliminary tightening is needed to place the inherited metric in
the chart.  Holding `H_N` fixed would preserve lower and Stein feasibility,
but need not preserve the upper bound after a transverse perturbation.  The
next section proves that on the normal anchor one can instead move safely to
the sharper zero-slack branch.  L106 then supplies a uniform transverse
expansion for `T-T_N`.

## 7. Zero-slack tightening preserves the normal upper bound

On the normal direct-sum stratum, `T_N,B_N,H_N` are block diagonal in copy
space.  The chart and endpoint therefore split into independent
single-copy problems.  For one copy, write

\[
 H_j=\begin{bmatrix}h_{11}&h_{12}\\\bar h_{12}&h_{22}
 \end{bmatrix}\succeq0 .
\]

At the Crabb base, (5) and (8) give

\[
 D_H{\cal E}(H_j)=h_{22}+2h_{11}.                    \tag{11}
\]

The upper Schur penalty has zero derivative there because its cross block
vanishes at `M`.  For a positive `2 x 2` matrix,

\[
 h_{22}+2h_{11}\ge\operatorname{tr}H_j\ge\|H_j\|.     \tag{12}
\]

Real analyticity now gives, uniformly over the single-block normal
neighbourhood,

\[
 {\cal E}(T_j,B_j,H_j)-{\cal E}(T_j,B_j,0)
 \ge(1-Cs)\|H_j\|,                                    \tag{13}
\]

where
`s=||T_j-C_3||+||B_j||+||H_j||` controls the distance to the base along the
segment from zero to `H_j`.  After shrinking so that `Cs<1`, the right side
is nonnegative.  Since the inherited normal metric has

\[
 {\cal E}(T_j,B_j,H_j)\le0,
\]

equation (13) proves

\[
 \boxed{{\cal E}(T_j,B_j,0)\le
 {\cal E}(T_j,B_j,H_j)\le0.}                          \tag{14}
\]

Taking direct sums proves the same statement in arbitrary copy
multiplicity.  Thus the frozen normal certificate may be tightened all the
way to **zero Stein Schur slack** without losing the upper bound.  This is a
local statement on the block-diagonal normal stratum; no monotonicity is
claimed for a general coupled `H`.

Consequently the transverse construction should use the zero-slack analytic
branch `P(T,B,0)`.  This is exactly the branch on which L98--L104 computed
the sharp cubic and fourth-order endpoints.

## 8. Weighted terminal slack transfer

The zero-slack tightening does more than preserve the normal upper bound: it
supplies exactly the quadratic margin which the transverse matrix consumes.
Work in the real canonical two-copy chart and write

\[
 D=\begin{bmatrix}d&a\\0&-d\end{bmatrix},\qquad
 w_{\rm cen}={3d^2\over8\sqrt2}.
\]

Along the weighted path of L100, the second inverse-map coefficient is

\[
 F_2(\zeta)=m_2\zeta+c\zeta^3,\qquad
 m_2={5(2d^2+a^2)\over128},\quad c={9d^2\over64}.     \tag{15}
\]

Put `Delta m=5a^2/128`, and let `m_3` be the mean top-support
coefficient from L100.  If the actual domain map is frozen and evaluated on
the block-diagonal normal part, the Stein Schur slack of the inherited
normal metric has coefficients

\[
\begin{aligned}
 H_2&=\operatorname{diag}(4\Delta m,8\Delta m,
                          4\Delta m,8\Delta m),\\
 H_3&=\begin{bmatrix}
 4m_3&\gamma&0&0\\
 \gamma&8m_3&0&0\\
 0&0&4m_3&-\gamma\\
 0&0&-\gamma&8m_3
 \end{bmatrix},
 \qquad
 \gamma={15\sqrt2\,a^2d\over256}.                    \tag{16}
\end{aligned}
\]

The ordering in (16) is the two active physical levels of copy one followed
by those of copy two.  Notice that the coefficient `H_3` need not itself be
positive; it is a coefficient of the positive analytic slack
`H(epsilon)`.  Exact propagation shows that (16) makes the second and third
metric coefficients identical to those of the unperturbed normal
certificate.

Retain a fraction `theta` of both coefficients in (16).  For the frozen
normal operator `T_N=f(N)` and the actual operator `T_A=f(A)`, respectively,
the upper endpoints are

\[
\begin{array}{c|cc}
 &{\cal E}_2&{\cal E}_3\\ \hline
 T_N&-\frac58(1-\theta)a^2I&
       -16(1-\theta)m_3I\\[1mm]
 T_A& \frac58\theta a^2I&
       -16(1-\theta)m_3I .
\end{array}                                           \tag{17}
\]

Thus freezing the inherited slack (`theta=1`) creates the positive
quadratic transverse endpoint seen in the numerical adversary.  Tightening
to `theta=0` gives

\[
 {\cal E}_2(T_N)=-{5a^2\over8}I,\qquad
 {\cal E}_2(T_A)=0,\qquad
 {\cal E}_3(T_N)={\cal E}_3(T_A)=-16m_3I.             \tag{18}
\]

The positive transverse quadratic change is therefore canceled by the
normal tightening margin, while the sharp negative cubic survives
unchanged.  L100 gives

\[
 m_3\ge {a(a^2+2d^2)\over64},
\]

so at physical scales `r=epsilon|d|`,
`delta=epsilon a`, the surviving coefficient is at most

\[
 -{\delta(\delta^2+2r^2)\over4}I.                    \tag{19}
\]

This is the missing **joint finite-jet cancellation**.  It is not yet the
full terminal-tube theorem: one must still prove that the exact chart
remainder, after subtracting the frozen normal endpoint, is uniformly
`o(delta(r^2+delta^2))`.  L106 supplies the transverse factor, while the
analytic chart supplies the Taylor majorant; the remaining work is to
organize those two facts without assuming differentiability of the varying
Riemann map.

## 9. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_exact_metric_chart.py
.venv/bin/python -u experiments/repeated_p3_slack_transfer.py
```

The checker uses a general complex `2 x 2` block `X12` and general Hermitian
`X11,X22`, constructs the repeated `C3` Stein defect, and proves (4) entry by
entry.  The additive parameter `H` does not change this Jacobian.  The same
block multiplication is independent of copy multiplicity.

The second checker independently constructs the weighted actual and
frozen-normal jets, derives (16), verifies both rows of (17), and confirms
that full inherited slack reproduces the unperturbed normal metric through
third order.
