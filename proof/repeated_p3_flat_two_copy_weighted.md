# Weighted two-copy terminal theorem (2026-07-22)

## 1. Weighted chart

Let

\[
 D=\begin{bmatrix}d&a\\0&-d\end{bmatrix},\qquad a>0,
\]

and use L88's physical perturbation `E_D`.  Consider the weighted path

\[
 A_\epsilon=I_2\otimes C_3+\epsilon E_D
 +\epsilon^2\{E_{zI+D_1}+I_2\otimes W(w)\},           \tag{1}
\]

where

\[
 W(w)=w(E_{10}+E_{21}).                               \tag{2}
\]

This is the transition left after L94--L97: the first direction lies on the
trace-zero nonnormal terminal face, while an arbitrary trace-zero Schur
tangent `D_1`, the scalar trace coordinate, and the common flat mode enter one
order later.  Differentiating the canonical first metric at
`D+epsilon D_1` puts the free copy block

\[
 U_1=-\frac{3\sqrt2}{8}D_1                              \tag{3}
\]

in the second metric.

## 2. Scalar second support

Let `V_1(q),V_2(q)` be the Hermitian support perturbations from the two terms
in (1), and let `R(q)` be the Crabb reduced resolvent.  Exact multiplication
gives

\[
 P\{V_2+V_1RV_1\}P=M_2(q)I_2,                         \tag{4}
\]

where

\[
\begin{aligned}
M_2(q)={}&m_2-\frac3{128}
 \{d^2q^2+\bar d^2q^{-2}\}\\
&+\frac{\sqrt2}{4}\{wq^{-2}+\bar wq^2\},\\
m_2={}&\frac5{128}(2|d|^2+a^2).                       \tag{5}
\end{aligned}
\]

The full second direction has scalar first support, so the
copy degeneracy survives through order two.  The second inverse-map term is

\[
 F_2(\zeta)=m_2\zeta+
 \left(-\frac3{64}d^2+\frac{\bar w}{\sqrt2}\right)\zeta^3. \tag{6}
\]

## 3. Third effective support

The complete third degenerate support matrix is

\[
 Q_3(q)=P\{V_1RV_1RV_1+V_1RV_2+V_2RV_1\}P.           \tag{7}
\]

It is Hermitian on the circle, but need not be traceless because `D_1`
changes the scalar second support at order three.  Scalar terms cannot affect
the matrix Jensen gap.  Put

\[
 Q_3^0(q)=Q_3(q)-\frac12\operatorname{tr}Q_3(q)I.
\]

The tangent `D_1` does not alter the following three traceless coefficients:

\[
\begin{aligned}
\overline {Q_3^0}
 &=\frac5{64}(\bar zD+zD^*),\\
\widehat {Q_3^0}(2)&=-\frac3{64}zD,\\
[q^1](Q_3^0)_{01}
 &=\frac{a}{128}(2a^2+4|d|^2).                        \tag{8}
\end{aligned}
\]

The last coefficient is independent of both weighted coordinates `z,w`; it
is the same coercive mode as L94.

## 4. Complete weighted third metric

Let

\[
 m_3=\mathop{\rm mean}_{|q|=1}\lambda_{\max}Q_3(q).
\]

Keep the possibly nonzero first conformal coefficient of that top branch.
Insert (6), its full Fréchet derivative on `E_D`, the second physical
direction from (1), and the free block (3) into L92's tight
lower/Stein/upper third metric.  Every deterministic term and the first
conformal coefficient cancel, leaving

\[
\boxed{
 {\cal E}_3=16\{\overline Q_3-m_3I_2\}\preceq0.}       \tag{9}
\]

This is a matrix Jensen identity one order above L83/L85, not merely the
scalar L92 endpoint.

## 5. Strictness

If `z=0`, the traceless mean vanishes.  Scalar Fourier terms cancel from the
matrix Jensen endpoint, while the last coefficient in (8) is nonzero.
Therefore (9) is negative definite.

Suppose `z!=0` and (9) has a null vector.  Equality in the averaged Hermitian
order would make that vector a common top eigenvector of every `Q_3(q)`.
It must therefore reduce both `overline {Q_3^0}` and every traceless Fourier
coefficient, including `hat Q_3^0(2)`.  But

\[
 [\overline {Q_3^0},\widehat {Q_3^0}(2)]
 =-\frac{15z^2}{4096}[D^*,D]\ne0,                    \tag{10}
\]

because `a>0` makes `D` nonnormal.  Two `2 x 2` matrices with a common
reducing line are simultaneously block diagonal and commute, contradicting
(10).  Therefore (9) is again negative definite.

We have proved:

\[
\boxed{\text{Every bounded weighted chart (1) at a nonnormal
two-copy terminal block descends strictly at order three.}}          \tag{11}
\]

On normalized blocks with `a` bounded below relative to `||D||`, and bounded
`D_1,z,w`, compactness makes the strict coefficient uniform.  The remaining
terminal issue is the approach `a->0` to L88's exact normal direct-sum
stratum, not another bounded weighted cancellation.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_two_copy_weighted.py
```

The checker uses arbitrary complex `d,z,w,D_1` and real Schur edge/tangent
coordinates.  It derives (4), every identity in (8), retains an arbitrary
complex first conformal coefficient, inserts the free block (3), rebuilds the
full third metric through the shared L77/L92 propagation module, and proves
(9) entry by entry.
