# The Gau--Wu cyclic corner is the top weighted support mode

> **Status and scope.**  The identity below is exact at every finite
> nondegenerate Gau--Wu model.  It identifies L346's southwest
> cyclic-closing generator with one coefficient of L343's polynomial
> shape.  It does not prove that this coefficient vanishes on L345's
> dual lift, so phase covariance and the final Hessian sign remain
> open.

## 1. Exact corner formula

Let \(S=S_\phi\) be the upper-triangular compressed shift for
\(\phi=zf\), with

\[
 \operatorname {diag}S=(0,b_1,\ldots,b_{n-2},0),
 \qquad p=e_0,\quad q=e_L,\quad L=n-1.             \tag{1}
\]

Assume that the interior zeros are nonzero, and write

\[
 f(z)=z\prod_{j=1}^{n-2}
       {z-b_j\over1-\overline b_jz},
 \qquad
 f'(0)=\prod_{j=1}^{n-2}(-b_j)\ne0.               \tag{2}
\]

Let \(A=XSX^{-1}\) be the Gau--Wu disk matrix.  The diagonal
similarity fixes the two endpoint vectors, so the same notation
\(p,q\) may be used in either coordinate system.

For an ambient first direction \(E\), let \(G_E\) be its
inverse-Riemann normalized first operator direction, and let

\[
 Y_1(E,\dot b)
 =Df(A)[G_E]+\dot f(A)                            \tag{3}
\]

be the first moving-Blaschke image.  Then

\[
\boxed{
 \langle q,Y_1(E,\dot b)p\rangle
 =4\,{f'(0)\over\gamma}\,
   \widehat{\sigma_E}(-n),}                       \tag{4}
\]

where

\[
 \sigma_E=\omega_\phi s_E,\qquad
 \omega_\phi=|\delta(\overline\zeta)|^2D_f,
 \tag{5}
\]

is L343's weighted polynomial support motion, and
\(\gamma\ne0\) is determined by the endpoint columns of L343's
polynomial frame:

\[
 {\cal C}e_0=\gamma p,\qquad {\cal C}e_L=q.        \tag{6}
\]

In the standard normalized Takenaka frame used by the experiments,
\(\gamma=1\).  Thus

\[
\boxed{
 (Y_1)_{L0}=4f'(0)\widehat{\sigma_E}(-n).}         \tag{7}
\]

In particular, motion of the zeros of \(f\) cannot alter the cyclic
corner.  The corner vanishes exactly when the top negative Fourier
coefficient of the weighted shape vanishes.

## 2. Highest coefficient of the support polynomial

Put

\[
 {\bf f}_\zeta=(1,\zeta,\ldots,\zeta^L)^T,\qquad
 M_E={\cal C}^*E{\cal C}.
\]

L343 gives

\[
 \sigma_E(\zeta)
 ={1\over4}{\bf f}_\zeta^*
   \{\overline\zeta M_E+\zeta M_E^*\}
   {\bf f}_\zeta.                                 \tag{8}
\]

The exponent \(-n=-L-1\) in the first term of (8) can only come
from the entry \((L,0)\).  The adjoint term supplies the conjugate
positive mode.  Therefore

\[
 \widehat{\sigma_E}(-n)={1\over4}(M_E)_{L0}.       \tag{9}
\]

The final column of L343's frame is
\({\cal C}e_L=Lr_0=q\).  The adjugate recurrence gives
\(Sr_L=0\), so its first column is
\({\cal C}e_0=\gamma p\); invertibility of \({\cal C}\) makes
\(\gamma\ne0\).  Hence

\[
 (M_E)_{L0}
 =\langle{\cal C}e_L,E{\cal C}e_0\rangle
 =\gamma\langle q,Ep\rangle.                     \tag{10}
\]

Equations (9)--(10) identify the support coefficient without a
pseudoinverse, Schur complement, or numerical phase choice.

## 3. The functional-calculus corner

The first inverse-Riemann correction has the analytic form
\(AH_{s_E}(A)\).  It is upper triangular, so

\[
 \langle q,G_Ep\rangle=\langle q,Ep\rangle.       \tag{11}
\]

For the same reason, every zero-motion term \(\dot f(A)\) is upper
triangular and has zero \((L,0)\) entry.  Finally, the divided
difference of \(f\) between the two endpoint eigenvalues \(0,0\) is
\(f'(0)\).  No other path can reach the southwest corner of an
upper-triangular base matrix.  Thus

\[
 \langle q,Df(A)[G_E]p\rangle
 =f'(0)\langle q,G_Ep\rangle.                    \tag{12}
\]

Combining (9)--(12) proves (4).

## 4. Consequence for the live flag

L346 reduced the observed upper-Schur response to the diagonal
conditions and the single cyclic condition
\((Y_1)_{L0}=0\).  Equation (4) now rewrites the latter as

\[
\boxed{\widehat{\sigma_{W_\phi}}(-n)=0.}          \tag{13}
\]

This is a useful separation:

1. the interior diagonal conditions are moving-root/interpolation
   equations involving the zero velocities;
2. the cyclic condition is purely physical and is exactly the last
   negative Hardy coefficient of L343's shape polynomial.

The next proof must derive (13) from L345's Euler equation in
L344's fixed Toeplitz port.  It may not use the zero velocities to
cancel the corner, because (4) proves that they do not enter it.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_cyclic_corner_shape.py
```

The checker independently constructs L343's polynomial frame,
forms \(\sigma_E\) from sampled support data, compares its highest
negative coefficient with \(({\cal C}^*E{\cal C})_{L0}/4\), and
compares (7) with a separate block functional-calculus jet.  It
also verifies directly that every pure zero-velocity direction has
zero cyclic corner.

Twelve models, two in each dimension \(3,\ldots,8\), give maximum
relative highest-mode and cyclic-corner residuals
\(6.36\cdot10^{-15}\) and \(6.35\cdot10^{-15}\), respectively.
The dataset SHA-256 is
`a4bad7eb40ff2aa41d65c8b7cee515acb4f8ae21b1ea4b0caa3eb353f8646829`.
