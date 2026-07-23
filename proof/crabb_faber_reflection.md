# Faber reflection of a palindromic equality tangent (2026-07-23)

## 1. Exact endpoint law

Put `p=L+1`, let `C=C_p`, and let `J` reverse coordinates.  For

\[
 S_0=C+cJCJ
\]

write `P_m=P_(m,c)` for the Dickson polynomials

\[
 P_0=2,\qquad P_1=z,\qquad P_m=zP_{m-1}-cP_{m-2}.
\]

Their two endpoint rows are

\[
\boxed{
 e_0^*P_m(S_0)=2e_m^*,\qquad
 e_L^*P_m(S_0)=2c^m e_{L-m}^*
 \quad(0\leq m\leq L).}                              \tag{1}
\]

For `1<=j<L`, put

\[
 E_j=2(e_0-e_L)e_{j+1}^*,\qquad
 H_j=E_j+cJE_{L-j}J
    =2(e_0-e_L)(e_{j+1}^*-ce_{j-1}^*).              \tag{2}
\]

Then the compatible equality/adjoint perturbation is invisible in the
two endpoint rows of the top Dickson polynomial:

\[
\boxed{
 e_0^*DP_L(S_0)[H_j]=
 e_L^*DP_L(S_0)[H_j]=0.}                             \tag{3}
\]

Equations (1)--(3) hold over the polynomial ring `C[c]`; they are not
asymptotic statements.

To prove them, differentiate the Dickson recurrence:

\[
 D_0=0,\quad D_1=H_j,\quad
 D_m=H_jP_{m-1}(S_0)+S_0D_{m-1}-cD_{m-2}.            \tag{4}
\]

A simultaneous induction on `m` and distance from the two endpoints
gives (1).  The same induction makes the cancellation in (3)
coefficientwise explicit.  If `D_m=DP_m(S_0)[H_j]`, then

\[
 {1\over2}e_0^*D_m=A_{m,j}-B_{m,j},\qquad
 e_L^*D_m=-e_0^*D_m,                                 \tag{4a}
\]

where

\[
\begin{aligned}
A_{m,j}
&=\begin{cases}
e_{j+m}^*,&j+m\leq L,\\
c^{j+m-L}e_{2L-j-m}^*,&j+m>L,
\end{cases}\\
B_{m,j}
&=\begin{cases}
c^m e_{j-m}^*,&m\leq j,\\
c^j e_{m-j}^*,&m>j.
\end{cases}                                         \tag{4b}
\end{aligned}
\]

These formulas start at `D_1=H_j`; substituting them in (4), with
(1) at the two boundary crossings, proves the next value of `m`.
At `m=L`, both terms in (4b) are
`c^j e_(L-j)^*`, so they cancel.  This proves (3) without an
asymptotic or genericity argument.  Equivalently, expand

\[
 P_L(z)=\sum_{r=0}^{\lfloor L/2\rfloor}
 (-1)^r{L\over L-r}\binom{L-r}{r}c^rz^{L-2r};
\]

inserting `H_j` once pairs the summand with `r` reverse steps from
`E_j` with the summand having `r-1` reverse steps from
`cJE_(L-j)J`.  The binomial recurrence is exactly (4), including the
two Crabb endpoint weights.

## 2. The reflected equality polynomial

Let

\[
 g_a(\xi)=\xi^L+2a\sum_{j=1}^{L-1}u_j\xi^j,\qquad
 u_j=\overline {u_{L-j}}.                            \tag{5}
\]

In coefficient coordinates the physical-adjoint tangent is the sum of
`u_j H_j`; the conjugation in the adjoint is essential for complex
directions.  Define its Faber transform

\[
 {\cal F}_{c}g_a(S)
 =P_L(S)+2a\sum_{j=1}^{L-1}u_jP_j(S).                \tag{6}
\]

Combining (1) and (3), the amplitude derivative at the Crabb axis has
the exact endpoint rows

\[
\boxed{
\begin{aligned}
e_0^*D_a{\cal F}_{c}g_a(S)|_{a=0}
 &=4\sum_j u_j e_j^*,\\
e_L^*D_a{\cal F}_{c}g_a(S)|_{a=0}
 &=4\sum_j u_jc^j e_{L-j}^*.
\end{aligned}}                                      \tag{7}
\]

The second row is the coefficient vector reflected across the chain
and weighted by its circle grade.  Its squared Euclidean norm is

\[
16\sum_j|u_j|^2c^{2j}.                               \tag{8}
\]

The same identity is visible on the scalar Joukowski boundary:

\[
P_j(\zeta+c/\zeta)=\zeta^j+c^j\zeta^{-j},
\]

so

\[
({\cal F}_{c}g_a)(\zeta+c/\zeta)
=g_a(\zeta)+c^L\zeta^{-L}
+2a\sum_j u_jc^j\zeta^{-j}.                          \tag{9}
\]

Parseval turns the negative-frequency part of (9) into

\[
c^{2L}+4a^2\sum_j|u_j|^2c^{2j},                     \tag{10}
\]

with no mixed circle grades.  Multiplication by the size-three
coefficient `-16` is exactly A84's candidate face.

## 3. What this proves and what remains

This closes two algebraic gaps in the proposed localization:

1. the grade-weighted vector is forced by an exact endpoint identity,
   not inferred from fitted Hessians; and
2. orthogonality of distinct grades is the ordinary Fourier/coordinate
   orthogonality in (8)--(10).

It does **not yet prove** that the optimized Stein envelope depends to
first Newton order only on this reflected endpoint vector.  That is the
remaining metric-normal-form statement.  A proof must show that the
other rows of (6), together with the lifted model-space metric, are
inactive or removable to strictly higher Newton weight.  L130 proves
that statement when only one central grade is present.  The next attack
is therefore a block perturbation/localization theorem for (6), rather
than another raw coefficient computation.

A tempting shortcut has already been falsified.  If `d_c` is the axis
rank-one defect and `U(u)` is lower triangular Toeplitz, the naive
transport

\[
 d_c+2aU(u)d_c                                      \tag{11}
\]

does reproduce every optimizer coefficient below the first reflected
grade.  It is not sharp at that grade.  For `k=3`, in sizes
`L=6,7,8`, (11) gives condition coefficient `+192c^6`; the true
stationarity equation adds `8c^3` in the relevant near-endpoint defect
coordinate and changes `+192` to `-64`.  Thus Parseval identifies the
normal vector, but the metric Schur correction is essential.  Any proof
that simply transports the disk defect without solving that correction
is wrong.

## 4. Exact regeneration

Run

```bash
.venv/bin/python -u experiments/crabb_faber_reflection.py \
  --output experiments/crabb_faber_reflection_s70223.jsonl
```

The checker uses exact SymPy polynomial arithmetic.  It verifies
(1), (3), and (7) through `L=14`, for every independent real and
imaginary phase-one coefficient direction.
