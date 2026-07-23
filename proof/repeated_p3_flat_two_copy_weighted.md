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

## 6. The sharp normal center and its transverse gap

The degeneration in the last paragraph can now be located exactly.  On the
normal face

\[
 a=0,\qquad z=0,\qquad D_1=0,
\]

write `Q_3^0=diag(h,-h)`.  Direct Laurent extraction gives

\[
\begin{aligned}
 [q^3]h&=\frac{d(3d^2-8\sqrt2\bar w)}{128},\\
 [q^1]h&=-\frac{\bar d(3d^2-8\sqrt2\bar w)}{128}.       \tag{12}
\end{aligned}
\]

The negative modes are their Hermitian conjugates.  Consequently, for
`d!=0`,

\[
 Q_3^0(q)\equiv0
 \quad\Longleftrightarrow\quad
 w=w_{\rm cen}(d):=\frac{3\bar d^2}{8\sqrt2}.           \tag{13}
\]

This is the unique weighted third-order center of the normal chart.  It lies
on L88's exact normal direct-sum manifold, so zero cubic endpoint there is
expected rather than a new equality obstruction.

More importantly, (8) is independent of `w` and `D_1`.  If `z=0` and `a>0`,
then `mean Q_3^0=0`, and the mean top eigenvalue dominates the modulus of the
displayed cross Fourier coefficient.  Thus even at (13),

\[
 {\cal E}_3\preceq
 -\frac{a(a^2+2|d|^2)}4 I_2.                           \tag{14}
\]

Hence the center is transversely strict, with leading margin
`-a|d|^2/2` as `a->0` for fixed nonzero `d`.  Equations (13)--(14) identify
the correct normal/tangential split for a tubular compactness argument.
They do **not** by themselves absorb analytic remainders when `a` tends to
zero on a later scale; that remains the uniformity debt.

## 7. Fourth-order persistence of the center

There is no hidden quartic copy splitting at (13).  Circle covariance makes
`d` real, and a diagonal copy unitary then makes the upper Schur edge real
without changing `d`; these two symmetries recover the general complex
case.  Let `N(q)` be the first support generated by `diag(d,-d)`, let `E(q)`
be the first support generated by the upper Schur edge, let `V_2(q)` be the
common support at (13), and let `R(q)` be the reduced Crabb resolvent.  Put
`M_2(q)I` for the scalar second effective support.

The fourth Feshbach coefficient is

\[
\begin{aligned}
K_4=P\{&
 NRNRNRN+V_2RNRN+NRV_2RN+NRNRV_2+V_2RV_2\\
&{}-M_2NR^2N\}P.                                      \tag{15}
\end{aligned}
\]

The last term is the energy-dependence correction in the reduced
resolvent; omitting it would not be a valid fourth effective Hamiltonian.
Exact multiplication gives

\[
\boxed{
 K_4-\frac12\operatorname{tr}(K_4)I=0,\qquad
 \left.\partial_a\left(
 K_4(a)-\frac12\operatorname{tr}(K_4(a))I
 \right)\right|_{a=0}=0.}                              \tag{16}
\]

Thus the weighted normal center stays copy-scalar through order four, and
the quartic traceless coefficient begins at least quadratically in the
nonnormal edge.  In physical normal/transverse scales `(r,delta)`, its size
is therefore

\[
 O\{\delta^2(r+\delta)^2\},                            \tag{17}
\]

whereas (14) has size comparable to
`delta(r^2+delta^2)`.  Their ratio tends to zero with `r+delta`.
Higher coefficients may again be linear in `delta`.  At the raw Feshbach
matrix level, subtracting the exact normal value forces at least one
transverse factor.  The remaining terminal-tube step is to preserve that
factor uniformly through L101's conformal response and through one
simultaneous feasible-metric selection.

## 8. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_two_copy_weighted.py
.venv/bin/python -u experiments/repeated_p3_flat_two_copy_fourth.py
```

The checker uses arbitrary complex `d,z,w,D_1` and real Schur edge/tangent
coordinates.  It derives (4), every identity in (8), retains an arbitrary
complex first conformal coefficient, inserts the free block (3), rebuilds the
full third metric through the shared L77/L92 propagation module, and proves
(9) entry by entry.  It also extracts (12), verifies the conjugated center in
(13), and substitutes it into the full traceless support.

The second checker independently constructs the fourth Feshbach coefficient
(15), including its energy correction, and proves both identities in (16)
entry by entry as rational Laurent identities.

The non-load-bearing nonlinear regression is

```bash
.venv/bin/python -u experiments/repeated_p3_normal_center_probe.py
```

It evaluates seven two-scale charts at map resolutions 2048 and 4096 and
compares the observed similarity deficit with
`delta(r^2+delta^2)`.  It is evidence only and is not used in (16).
