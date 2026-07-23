# Second-order Stein sign on the repeated-`C3` cross quotient (2026-07-22)

## 1. Statement and scope

Let `A_0=C_3 direct-sum C_3`, use the L74 cross perturbation

\[
 E=\begin{bmatrix}0&Y\\X&0\end{bmatrix},\qquad
 X=\sum_{j=0}^2\bar\alpha_jX_j,\quad
 Y=\sum_{j=0}^2\alpha_jY_j,                         \tag{1}
\]

and let `T_epsilon` be the numerical-range Riemann pullback.  L75 gives, for
positive `epsilon`,

\[
 T_\epsilon=A_0+\epsilon E-\epsilon^2k_0A_0+o(\epsilon^2), \tag{2}
\]

where

\[
\begin{aligned}
k_0={}&\frac5{128}|\alpha_0|^2+\frac14|\alpha_1|^2
 +\frac5{72}|\alpha_2|^2+\frac{\sqrt2}{12\pi}|d|,\\
d={}&3\alpha_0\bar\alpha_1+4\alpha_1\bar\alpha_2 .   \tag{3}
\end{aligned}
\]

There is an explicit feasible second-order metric certificate whose endpoint
coefficient is

\[
 \boxed{e=-8|\alpha_1|^2-\frac{4\sqrt2}{3\pi}|d|\le0.} \tag{4}
\]

Consequently the L21 similarity square satisfies the one-sided bound

\[
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(T_\epsilon)-4}{\epsilon^2}\le e\le0.        \tag{5}
\]

Equality in this certificate occurs exactly when `alpha_1=0`; then `d=0`
automatically.  Formula (5) handles one multiplicity-two cross pair with zero
diagonal perturbation blocks.  L77 subsequently proves strict cubic descent
on the flat `(alpha_0,alpha_2)` plane, and L78 proves the arbitrary-multiplicity
star sign at second order.  Mixtures with diagonal directions and the
higher-multiplicity equality set remain open.

## 2. First metric tangent

Work in copy-major coordinates and put

\[
 P_0=\operatorname{diag}(1,2,4,1,2,4).
\]

The lower, upper, and Stein active kernels are respectively the two copy
vectors at levels zero, two, and levels one-and-two.  Define a Hermitian first
metric tangent

\[
 {\cal X}=\begin{bmatrix}0&Z\\Z^*&0\end{bmatrix},\qquad
 Z=\begin{bmatrix}
 0&-3\sqrt2\alpha_0/8&0\\
 \alpha_2/\sqrt2&-2\sqrt2\alpha_1&3\sqrt2\alpha_0/4\\
 0&-\sqrt2\alpha_2&0
 \end{bmatrix}.                                       \tag{6}
\]

Direct substitution shows that all three active first compressions vanish.
Thus (6) lies on the first-order tangent face required by the L62
second-order PSD lemma.

## 3. Construct the complete second metric

Let `D_1` be the first Stein tangent and let `K` and `R` denote its active
kernel and level-zero range.  Since the positive range block of the base Stein
defect is the identity, its Schur penalty is

\[
 C=(D_1)_{KR}(D_1)_{RK}.                                \tag{7}
\]

The lower endpoint penalty obtained from `P_0-I` is

\[
 L=\operatorname{diag}\left(
 \frac9{32}|\alpha_0|^2,\frac12|\alpha_2|^2\right),    \tag{8}
\]

and the upper penalty obtained from `4I-P_0` is

\[
 U=\operatorname{diag}\left(
 |\alpha_2|^2,\frac9{16}|\alpha_0|^2\right).           \tag{9}
\]

Write the positive forcing subtracted from `Y-A_0^*YA_0` in the second Stein
tangent as

\[
 F=H^*P_0A_0+A_0^*P_0H+E^*P_0E
   +E^*{\cal X}A_0+A_0^*{\cal X}E,\qquad H=-k_0A_0.    \tag{10}
\]

Regard the second metric `cal Y` as `2 x 2` copy blocks
`cal Y_ij` indexed by Crabb levels `i,j=0,1,2`.  Set

\[
 {\cal Y}_{00}=L,\qquad {\cal Y}_{0j}={\cal Y}_{j0}=0
 \quad(j>0),                                           \tag{11}
\]

and, with `w_0=w_1=sqrt(2)`, recursively define

\[
 {\cal Y}_{ij}=w_{i-1}w_{j-1}{\cal Y}_{i-1,j-1}
                  +(F+C)_{ij},\qquad i,j\ge1.           \tag{12}
\]

Here `C` is extended by zero off the Stein kernel.  The forcing is Hermitian,
so (11)--(12) construct a Hermitian full matrix.  They make the lower and
Stein second-order Schur complements equal to zero.  In particular,

\[
 {\cal Y}_{22}=4L+2(F+C)_{11}+(F+C)_{22}.               \tag{13}
\]

## 4. Endpoint collapse and sign

Exact multiplication in (6)--(13) gives the striking scalar identity

\[
 {\cal Y}_{22}+U=\left(
 \frac58|\alpha_0|^2-4|\alpha_1|^2
 +\frac{10}{9}|\alpha_2|^2-16k_0\right)I_2.            \tag{14}
\]

Choosing the upper metric coefficient `e` equal to the scalar in (14) makes
the upper Schur complement zero as well.  Substitution of (3) cancels the
`alpha_0` and `alpha_2` squares exactly and yields (4).  The strict-lift and
limiting argument of L62 then proves (5).

The cancellation explains the computations after L75: generator 1 decreases
the certificate quadratically, whereas generators 0 and 2 are genuinely
second-order flat.  L77 computes their first nonzero term at cubic order.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_stein_sign.py
```

The checker constructs (6), verifies the active first Stein compression,
builds the complete second metric by (11)--(12), checks all three second-order
Schur complements exactly, derives the scalar endpoint (14), and proves the
symbolic reduction to (4) for arbitrary complex coefficients.
