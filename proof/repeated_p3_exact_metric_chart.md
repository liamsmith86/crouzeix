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

For every operator `T` near `A_0` and every free complex block

\[
 B=(P_{01}\ P_{02})\in M_{m,2m}
\]

near zero, there is a unique Hermitian range block

\[
 C=\begin{bmatrix}P_{11}&P_{12}\\P_{21}&P_{22}\end{bmatrix}
\]

near `diag(2I,4I)` such that the lower metric constraint and the Stein
constraint have zero Schur complements on their base kernels.  The resulting
metric `P=P(T,B)` is real analytic in `(T,B)`.

Moreover, `P>=I` and `P-T^*PT>=0` hold automatically after shrinking the
neighbourhood.  The remaining condition `P<=4I` is equivalent to one
explicit analytic `m x m` endpoint inequality.

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
the base, `S_RR` is positive definite.  Define the Hermitian `2m x 2m`
equation

\[
 {\cal F}(T,B,C)
 =S_{KK}-S_{KR}S_{RR}^{-1}S_{RK}=0.                  \tag{3}
\]

At `(T,B,C)=(A_0,0,diag(2I,4I))`, equation (3) holds.  If

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
gives a unique analytic solution `C=C(T,B)` to (3).

Since `S_RR` remains positive, the Schur complement in the other direction
now gives

\[
 S(T,B,C(T,B))\succeq0.                               \tag{6}
\]

Equations (1)--(6) prove the asserted analytic parameterization of every
tight lower/Stein certificate in this chart.

## 4. One exact upper endpoint

For `Q=4I-P`, split off physical level two as its base kernel and let
`L={0,1}`.  The range block `Q_LL` is positive near `diag(3I,2I)`.  Hence

\[
 P\preceq4I
\quad\Longleftrightarrow\quad
 {\cal E}(T,B)\preceq0,                               \tag{7}
\]

where

\[
 \boxed{
 {\cal E}(T,B)
 =P_{22}-4I+
 P_{2L}(4I-P_{LL})^{-1}P_{L2}.}                       \tag{8}
\]

Thus the repeated-`C3` complete-similarity problem is locally reduced to
finding the free block `B` for which the single analytic endpoint (8) is
nonpositive.  No infinite forced-metric recursion remains: for any analytic
choice of `T` and `B`, all the other metric coefficients are the convergent
Taylor series of `C(T,B)`.

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

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_exact_metric_chart.py
```

The checker uses a general complex `2 x 2` block `X12` and general Hermitian
`X11,X22`, constructs the repeated `C3` Stein defect, and proves (4) entry by
entry.  The same block multiplication is independent of copy multiplicity.
