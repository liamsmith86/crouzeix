# Radial-gap coupling theorem for repeated `C3` (2026-07-22)

## 1. Model and theorem

Let `d>0` and consider two repeated Crabb copies with one selected radial
winner:

\[
 A_\epsilon=
 \begin{bmatrix}(1+d\epsilon)C_3&\epsilon Y\\
                 \epsilon X&(1-d\epsilon)C_3\end{bmatrix},             \tag{1}
\]

where `(X,Y)` is the L74 cross quotient with coefficients
`alpha_0,alpha_1,alpha_2`.  The selected copy is the unique first-order top
support branch for every boundary angle.  The numerical-range Riemann
pullback satisfies

\[
 \boxed{
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(T_\epsilon)-4}{\epsilon^2}
 \le-\frac{25}{8}|\alpha_0|^2-8|\alpha_1|^2
       -\frac{50}{9}|\alpha_2|^2<0}                   \tag{2}
\]

for every nonzero cross direction.  The coefficient is independent of the
fixed positive gap `d`; the size of the neighbourhood in which the expansion
is useful is not uniform as `d` tends to zero.

For one selected copy and any number of equally shrinking orthogonal copies,
(2) holds with each square replaced by the squared norm of the corresponding
coefficient vector.  Thus this is the first exact diagonal--cross coupling
theorem in the repeated-block campaign.

## 2. Pulled operator through second order

The first support derivative is the constant `d`.  The cross compression is
zero by L74.  At second order the selected effective support branch is the
non-maximized diagonal branch from L75.  Its mean is

\[
 m_0=\frac5{128}|\alpha_0|^2+\frac14|\alpha_1|^2
       +\frac5{72}|\alpha_2|^2,                       \tag{3}
\]

and its first Fourier contribution evaluated at `C3` is

\[
 \frac{\sqrt2}{24}
 (3\alpha_0\bar\alpha_1+4\alpha_1\bar\alpha_2)C_3^2. \tag{4}
\]

If `E` denotes the cross matrix in (1), inversion of the domain map gives

\[
 T_\epsilon=A_0+\epsilon G+\epsilon^2H+o(\epsilon^2), \tag{5}
\]

where

\[
 G=\begin{bmatrix}0&Y\\X&-2dC_3\end{bmatrix}.        \tag{6}
\]

The diagonal blocks of `H` are minus (3)--(4) on the selected copy and
`2d^2C_3` minus (3)--(4) on the other copy; its cross blocks are `-d(X,Y)`.
These terms are all retained in the exact checker.

## 3. A first metric with strict orthogonal slack

On the nonselected copy, add

\[
 \operatorname{diag}(d,0,-d)                           \tag{7}
\]

to the first metric tangent.  This is positive on the lower metric endpoint
and negative on the upper one, as required.  The first Stein compression on
its two positive Crabb levels becomes

\[
 \operatorname{diag}(6d,15d)\succ0.                   \tag{8}
\]

The selected compression and every selected/cross compression remain zero.
Hence only the selected copy remains active at second order; the other copy
absorbs all sufficiently small second-order remainders through (7)--(8).

The complete cross first metric tangent has two free real scalars `x,y`:

\[
 Z(x,y)=\begin{bmatrix}
0&x\alpha_0&0\\
y\alpha_2&-2\sqrt2\alpha_1&(2x+3\sqrt2/2)\alpha_0\\
0&(2y-2\sqrt2)\alpha_2&0
\end{bmatrix}.                                        \tag{9}
\]

The non-free entries in (9) are exactly the first active Stein equations.

## 4. Selected endpoint

Insert the selected lower and contraction Schur penalties, recursively build
the selected second metric through its two Crabb levels, and add the upper
penalty.  Exact multiplication gives

\[
\begin{aligned}
e(x,y)={}&\frac{96x^2-48\sqrt2x-63}{24}|\alpha_0|^2
-8|\alpha_1|^2\\
&+\frac{96y^2-256\sqrt2y+208}{24}|\alpha_2|^2.        \tag{10}
\end{aligned}
\]

Every occurrence of the radial gap and the first Fourier term (4) cancels
from this endpoint.  The two scalar quadratics have independent minimizers

\[
 x=\frac{\sqrt2}{4},\qquad y=\frac{4\sqrt2}{3}.       \tag{11}
\]

Substitution gives exactly the right side of (2).  The orthogonal constraints
are already strict at first order, and the selected endpoint is strict at
second order, so the usual strict-lift argument proves (2).

For several orthogonal copies, use (7) and (9) in each selected/cross block.
All selected scalar penalties add, producing the coefficient-vector version
of (2); the strict first-order blocks remain a direct sum.

## 5. Checks, scope, and regeneration

For `d=1`, nonlinear Riemann-map/SDP quotients converge numerically to
`-25/8`, `-8`, and `-50/9` on the three pure generators.  The convergence is
slower for larger `d` because the asymptotic regime requires
`epsilon` small compared with the fixed first-order gap.

Run

```bash
.venv/bin/python -u experiments/repeated_p3_radial_gap.py
```

The checker reconstructs (3)--(6), verifies every first active compression
and the strict slack (8), eliminates the complete selected second metric, and
proves (10)--(11) symbolically for arbitrary complex coefficients.

The theorem treats radial diagonal splitting plus pure star blocks.  General
diagonal perturbations, blocks internal to the losing copy space, ties between
several first-order winners, and the weighted transition `d -> 0` remain.
