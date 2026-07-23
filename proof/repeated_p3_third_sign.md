# Third-order descent on the repeated-`C3` flat plane (2026-07-22)

## 1. The L76 equality plane

Keep the multiplicity-two cross perturbation of L74--L76 and set
`alpha_1=0`.  Write

\[
 \alpha=\alpha_0,\qquad \beta=\alpha_2,
 \qquad E=E(\alpha,\beta).                              \tag{1}
\]

L76 gives zero at second order for every such pair.  The next coefficient is
strictly negative except at the origin.  More precisely, define on `|q|=1`

\[
\begin{aligned}
c(q)={}&-\frac{(3\alpha q^2+4\beta)}{1728q^3}\Bigl(
-9|\alpha|^2q^2+18\alpha\bar\beta q^4
-18\bar\alpha\beta+16|\beta|^2q^2\Bigr),              \tag{2}\\
m_3={}&\frac1{2\pi}\int_0^{2\pi}|c(e^{i\theta})|\,d\theta. \tag{3}
\end{aligned}
\]

Then the L21 similarity square of the numerical-range Riemann pullback obeys

\[
 \boxed{
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(T_\epsilon)-4}{\epsilon^3}\le-16m_3<0
 }
 \quad\text{when }(\alpha,\beta)\ne(0,0).              \tag{4}
\]

This closes the pure-pair equality plane through its first nonzero order.  It
does not yet treat diagonal/cross mixtures or simultaneous directions in more
than one orthogonal copy.

## 2. Third effective support split

Let `H(q)=Re(q^{-1}C_3)`, let `r(q)` be its top unit vector, and let

\[
 R(q)=(I-H(q))^\dagger=I-\frac14H(q)-\frac34H(q)^2.
\]

The first compression of the support perturbation `V(q)` is zero.  On the
L76 flat plane, exact multiplication gives a scalar second effective matrix

\[
 PVRVP=M_2(q)I_2,                                      \tag{5}
\]

where

\[
 M_2(q)=m_2+\frac1{32}
 (\alpha\bar\beta q^2+\bar\alpha\beta q^{-2}),
 \qquad
 m_2=\frac5{128}|\alpha|^2+\frac5{72}|\beta|^2.       \tag{6}
\]

Because (5) is scalar, the third degenerate perturbation matrix is simply

\[
 PVRVRVP=
 \begin{bmatrix}0&c(q)\\\overline{c(q)}&0\end{bmatrix}. \tag{7}
\]

Here the lower entry means the circle adjoint, so (7) is Hermitian on
`|q|=1`.  Its top eigenvalue is `|c(q)|`.  Uniformly in the boundary angle,

\[
 h_\epsilon(q)=1+\epsilon^2M_2(q)
                    +\epsilon^3|c(q)|+O(\epsilon^4).   \tag{8}
\]

The Laurent polynomial `c` contains only odd powers of `q`; hence `|c|` is
`pi`-periodic and has zero first Fourier coefficient.

## 3. The third Riemann pullback

No normal-angle correction enters before order four because the first domain
motion is order two.  The second inverse-map term obtained from (6) is

\[
 F_2(w)=m_2w+\frac{\alpha\bar\beta}{16}w^3.            \tag{9}
\]

Let `F_3` be the Schwarz transform of `|c|`.  Periodicity and `C_3^3=0` give

\[
 F_3(C_3)=m_3C_3.                                      \tag{10}
\]

Expanding the inverse map at the perturbed operator therefore gives

\[
\begin{aligned}
T_\epsilon={}&A_0+\epsilon E-\epsilon^2m_2A_0
 +\epsilon^3H_3+o(\epsilon^3),\\
H_3={}&-m_2E-\frac{\alpha\bar\beta}{16}
 (A_0^2E+A_0EA_0+EA_0^2)-m_3A_0.                     \tag{11}
\end{aligned}
\]

The absolute value in (8) is the only nonsmooth ingredient.  It is continuous
and the finite support cluster remains uniformly separated from the other
four support eigenvalues, so the expansion and disk shape derivative are
uniform in `q`.

## 4. Third metric construction

Use the L76 first metric tangent `cal X` and its recursively constructed
second metric `cal Y`, now with `alpha_1=0`.  For any tight PSD block

\[
 \begin{bmatrix}0&0\\0&R_0\end{bmatrix}
 +\epsilon M_1+\epsilon^2M_2+\epsilon^3M_3,
\]

the order-three Schur endpoint forced after the order-two endpoint is tight is

\[
 B_2R_0^{-1}B_1^*+B_1R_0^{-1}B_2^*
 -B_1R_0^{-1}(M_1)_{RR}R_0^{-1}B_1^*,                \tag{12}
\]

where `B_j=(M_j)_{KR}`.  Apply (12) to the lower metric, Stein, and upper
metric constraints.

For the lower constraint the forced third endpoint vanishes.  Set the other
level-zero blocks of the third metric `cal Z` to zero.  On positive Crabb
levels recursively define

\[
 {\cal Z}_{ij}=2{\cal Z}_{i-1,j-1}
                  +(F_3^{\rm Stein}+C_3^{\rm Schur})_{ij},
 \qquad i,j\ge1.                                      \tag{13}
\]

This makes the complete third Stein Schur complement zero.  Exact propagation
to the upper endpoint then collapses to

\[
 {\cal Z}_{22}+U_3=-16m_3I_2.                         \tag{14}
\]

All terms involving `m_2`, the cubic coefficient in (9), and the lower-order
metric cancel.  A strict perturbation of the finite tangent LMIs followed by
the L62 limiting argument proves (4).

## 5. Strictness and checks

The product in (2) cannot vanish identically unless both `alpha` and `beta`
vanish.  Indeed, the first factor is identically zero only at the origin.  If
the second factor is identically zero, its `q^4`, `q^2`, and constant
coefficients force

\[
 \alpha\bar\beta=0,\qquad 9|\alpha|^2=16|\beta|^2,
\]

and again both coefficients vanish.  Thus `m_3>0` away from the origin.  It
is continuous and homogeneous of degree three, so the negativity is uniform
on the unit sphere of this four-real-dimensional plane.

The pure axes provide simple checks:

\[
 \beta=0:\quad -16m_3=-\frac14|\alpha|^3,
 \qquad
 \alpha=0:\quad -16m_3=-\frac{16}{27}|\beta|^3.       \tag{15}
\]

For the fresh mixed check `alpha=.7+.2i`, `beta=-.3+.8i`, the predicted
coefficient is `-.5822766772`; the numerical quotients at
`epsilon=.1,.05,.025` are `-.5696246,-.5757497,-.5783211`.  No numerical
estimate is used in the proof.

Run

```bash
.venv/bin/python -u experiments/repeated_p3_third_sign.py
```

The checker reconstructs (5)--(7), (9)--(11), every lower/Stein/upper
third-order Schur term, the full recursive third metric, and the scalar
endpoint (14) exactly for arbitrary complex `alpha,beta`.
