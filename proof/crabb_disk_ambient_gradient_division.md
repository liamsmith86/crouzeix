# The proposed `O(Q)` ambient-gradient division is false (2026-07-23)

## 1. Falsified target

Let `A_z` be L122's normalized Toeplitz disk chart and

\[
{\cal Q}(z)=\|z\|^4-|z^TJz|^2.
\]

The proposed shortcut after L162 was

\[
\left|D\Gamma(A_z)[E-h_{z,E}(A_z)]\right|
\le C_L{\cal Q}(z)\|E\|,                            \tag{1}
\]

uniformly in every ambient direction `E`.  A weaker version asked for
(1) only after projecting onto L115's `2p-4` coercive
circular-normal covectors.

The full ambient target is exactly false already in size four.  The
coercive-projected target is numerically false with stable
high-resolution controls.  Along deterministic transverse disk paths

\[
z(\delta)=u+\delta v,\qquad {\cal Q}(z(\delta))\asymp\delta^2,
\]

the characteristic-Blaschke ambient gradient and its complete
coercive-normal projection satisfy

\[
\boxed{
\|\nabla_{\rm amb}R\|\asymp\sqrt{\cal Q},\qquad
\|\Pi_{\rm coercive}\nabla_{\rm amb}R\|
\asymp\sqrt{\cal Q},
}                                                     \tag{2}
\]

where

\[
R(z,E)=D\|B_z(T)\|_{K_z}^2
       [E-h_{z,E}(A_z)],\qquad
B_z={g_z\over g_z^\sharp}.
\]

Thus cone division cannot replace the sharp joint normal
Schur/covariant calculation.

## 2. Direct finite calculation

For every disk coordinate, construct

\[
K=H+R^*HR,\qquad A=2K^{-1}HR,
\qquad \det(\xi I-A)=\xi g(\xi).
\]

The fixed-inner Frechet derivative is evaluated without root
tracking:

\[
DB(A)[F]
=Dg(A)[F]\,g^\sharp(A)^{-1}
-B(A)Dg^\sharp(A)[F]\,g^\sharp(A)^{-1}.             \tag{3}
\]

Horner recurrences evaluate both polynomial derivatives exactly up to
floating arithmetic.  The normalized ambient direction is

\[
F=E-h_E(A),
\]

where the support Fourier coefficients of

\[
s_E(w)
={\operatorname{Re}\{\bar w f(w)^*KEf(w)\}
  \over f(w)^*Kf(w)}
\]

give `h_E`, exactly as in L162.  The top singular value of
`K^(1/2)B(A)K^(-1/2)` is simple, so ordinary singular differentiation
gives every real ambient coordinate of `R`.

### 2.1 Exact size-four counterexample

At the Crabb apex take

\[
\begin{aligned}
z_1&=x(37/1000+2i/125),\\
z_2&=x(11/250+31i/1000),
\end{aligned}
\]

and the fixed real ambient direction `E_(1,2)`.  Exact
Gaussian-rational Laurent expansion of the support correction and
Horner expansion of (3) give

\[
\boxed{
{\cal Q}(z)
={25281\over15625000000}x^4,\qquad
R(z,E_{12})=-{424\over15625}x^2.
}                                                     \tag{4}
\]

The nonzero quadratic coefficient against a positive quartic
coefficient proves that the full ambient estimate (1) is false; no
floating inference is involved.

The transported circular-normal covector space in the broader
adversarial audit is computed
independently from the same support map.  Its rows are the real and
imaginary Fourier coefficients in modes `3,...,p`.  Their rank is
exactly

\[
2(p-2)=2p-4,
\]

and orthogonal projection of the full gradient onto their row space
gives the second quantity in (2).  This is the correct transported
normal space, not a fixed Crabb-point basis.

## 3. Adversarial coercive-projection results

The persisted audit covers dimensions `p=4,...,7` and transverse
scales

\[
\delta=.03,\ .01,\ .003,\ .001.
\]

For each fixed dimension:

* `Q/delta^2` is stable;
* the full gradient divided by `sqrt(Q)` is stable;
* the coercive projection divided by `sqrt(Q)` is stable and nonzero;
* the equality control at `delta=0` is between `2.1e-14` and
  `3.5e-14`; and
* the projected normal rank is exactly `2p-4`.

For example, at `p=5` the projected ratios are

\[
1.2761,\ 1.2895,\ 1.2942,\ 1.2956,
\]

while at `p=7` they are

\[
6.6185,\ 6.6190,\ 6.6196,\ 6.6198.
\]

Changing the boundary resolution through
`1024,2048,4096,8192` leaves the displayed `p=5,delta=.01` ratio
unchanged in the first twelve reported digits.  This rules out FFT
aliasing as the explanation.

## 4. Why the same obstruction reaches the upper envelope

Let `widehat Gamma=4+Gamma` be the feasible optimized rank-one upper
certificate and put

\[
\Delta(z,Y)
=\widehat\Gamma(A_z+Y)-\|B_z(T(A_z+Y))\|^2\ge0.
\]

L157 and the explicit model upper certificate give

\[
0\le\Delta(z,0)\le C_L{\cal Q}(z)^2.
\]

The nonnegative-gradient estimate from L158 therefore still proves

\[
\|D_Y\Delta(z,0)\|\le C_L{\cal Q}(z).                \tag{5}
\]

Equation (5) says the upper and dual ambient gradients differ only by
`O(Q)`.  It cannot cancel the nonzero `sqrt(Q)` term in (2).
Consequently (1) fails for the optimized upper envelope as well.

This also explains why L156's scalar endpoint residual could not prove
(1): the missing full endpoint vectors carry a genuine first normal
jet rather than a removable proof artifact.

## 5. Reproduction and consequence

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/crabb_off_equality_dual_gradient.py \
  --output experiments/crabb_off_equality_dual_gradient_s70223.jsonl

OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/crabb_off_equality_dual_gradient_exact.py
```

The result does not contradict the conjecture.  It falsifies only the
attempt to absorb circular normals using the scalar disk deficit and
an `O(Q)` gradient.  The campaign must retain the sharp joint
quadratic/covariant cancellation isolated in L160--L163.
