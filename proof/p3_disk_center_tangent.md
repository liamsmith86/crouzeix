# Exact stationarity of the `p=3` disk-center certificate (2026-07-22)

## 1. Statement and scope

Let `D(e)` be the analytic disk-matrix curve from L71, let `c(e)` and `r(e)>0`
be its disk center and radius, and normalize

\[
 B_e=\frac{D(e)-c(e)I}{r(e)}.
\]

For every sufficiently small `e`, the numerical range of `B_e` is the unit disk.
This note proves that the explicit rank-one Stein certificate of condition four on
this curve is stationary under **every complex matrix perturbation after the
first-order Riemann-map correction**.  Equivalently, the feasible-certificate
function has zero first derivative in every ambient direction along the disk
curve.

This is the critical-manifold statement suggested by L70--L71.  By itself it is
not a local-neighbourhood theorem; L73 subsequently combines it with the
weighted leading form and symmetry to obtain the required uniform sign.

## 2. Canonical real Schur form

L71 gives the disk factorization and a double eigenvalue at the center.  A further
exact invariant calculation gives

\[
 \sum_{|I|=2}\det\big((D-cI)^T(D-cI)\big)_{I,I}=4r^4.       \tag{1}
\]

The regeneration script proves (1) modulo the same exact quadratic relation that
defines the disk curve.

Put `l=tr(B_e)=-3c/r`.  A real Schur form, ordered with the two zero eigenvalues
at the endpoints, is

\[
 T=\begin{bmatrix}0&a&b\\0&l&d\\0&0&0\end{bmatrix}.
\]

Signs may be chosen so that `a,d>0` near `e=0`.  Comparing its Kippenhahn
polynomial with `(z+lx)(z^2-x^2-y^2)` gives

\[
 a^2+b^2+d^2=4,\qquad l(a^2+d^2)+abd=0.                    \tag{2}
\]

The normalized version of (1) is

\[
 (ad-bl)^2=4.                                               \tag{3}
\]

The sign is `ad-bl=2` by continuity from `e=0`.  Set `X=a^2+d^2` and
`Y=ad`.  Equations (2)--(3) give

\[
 b=-\frac{lX}{Y},\qquad Y+\frac{l^2X}{Y}=2,
 \qquad X+\frac{l^2X^2}{Y^2}=4.
\]

Writing `t=X/Y`, the last two equations reduce to `2t=4`.  Hence `X=2Y`,
so `a=d`; (2) then gives the exact canonical form

\[
 \boxed{T_l=
 \begin{bmatrix}
 0&a&-2l\\0&l&a\\0&0&0
 \end{bmatrix},\qquad a^2=2(1-l^2).}                       \tag{4}
\]

Thus (4) is not an assumed classification of all `3 x 3` disk matrices: the
additional singular-product identity (1) selects this one-parameter subfamily.

## 3. Rank-one Stein certificate and its endpoint derivative

For `|l|<1`, put

\[
 P_0=\operatorname{diag}(1,2,4),\qquad c_0=(1,0,0)^T.
\]

Direct multiplication using (4) gives

\[
 P_0-T_l^TP_0T_l=c_0c_0^T.                                 \tag{5}
\]

Consequently `P_0` is a contraction metric of condition four.  Perturb the
operator by `G` and the rank-one defect by `d`, fixing the harmless defect scale
with `d_0=0`.  If `X` is the metric derivative, differentiation of (5) gives

\[
 X-T_l^*XT_l
 =G^*P_0T_l+T_l^*P_0G+c_0d^*+dc_0^*.                     \tag{6}
\]

The endpoint eigenvalues of `P_0` are simple, so the derivative of
`lambda_max(P)/lambda_min(P)` is `X_22-4X_00`.  Solving the finite Stein
system (6) exactly gives

\[
 X_{22}-4X_{00}=4\,\ell(\operatorname{Re}G),               \tag{7}
\]

Here `Re G` means the entrywise real part in this real Schur basis (not the
Hermitian part of `G`).

where

\[
 \ell(E)=a(E_{01}+E_{12})+2lE_{11}-lE_{02}                 \tag{8}
\]

for a real matrix `E`.  All real and imaginary components of `d` cancel from
(7).  A pure imaginary operator direction also contributes zero.  The exact
symbolic solve in `experiments/p3_disk_center_tangent.py` audits these claims for
all matrix entries simultaneously.

## 4. First conformal correction

Write `q=e^{i theta}`.  Parameterize

\[
 l=\frac{2\rho}{1+\rho^2},\qquad
 a=\frac{\sqrt2(1-\rho^2)}{1+\rho^2},\qquad |\rho|<1.
\]

A normalized top support eigenvector for
`Re(q^{-1}T_l)` has numerator

\[
 n(q)=(q^{-1}-l,\ a,\ q-l)^T                             \tag{9}
\]

and squared norm

\[
 N(q)=4\left(1-\frac l2(q+q^{-1})\right).                  \tag{10}
\]

Both the eigenvector equation and (10) are exact identities in the audit script.
For a perturbation `E`, its support derivative is therefore

\[
 s_E(q)=\frac{\bar n(q)^T(E/q+qE^*)n(q)}{2N(q)}.           \tag{11}
\]

Let `s_hat(k)` be its Fourier coefficients and define the positive analytic
value

\[
 \mathcal S_E(l)=\widehat s(0)+\sum_{k\ge1}\widehat s(k)l^k
 =\frac1{2\pi i}\int_{|q|=1}\frac{s_E(q)}{q-l}\,dq.        \tag{12}
\]

The rational integrand in (12) has only the poles `0,rho,l` inside the circle.
Taking their residues gives the two real-linear identities

\[
 \mathcal S_E(l)=\frac14\ell(E)\quad(E\text{ real}),
 \qquad \operatorname{Re}\mathcal S_{iE}(l)=0
 \quad(E\text{ real}).                                    \tag{13}
\]

The script checks the first identity with nine symbolic real entries and the
second on the nine imaginary matrix units, which is equivalent by real linearity.
The calculation is made for nonzero `rho`; both identities extend to `rho=0` by
continuity (and agree there with the Crabb calculation in L61).

The first shape derivative of the disk-to-numerical-range map is the Schwarz
integral

\[
 h(z)=\widehat s(0)z+2\sum_{k\ge1}\widehat s(k)z^{k+1}.    \tag{14}
\]

Hence the pulled-back operator perturbation is

\[
 G=E-h(T_l).                                                \tag{15}
\]

Finally, direct multiplication gives

\[
 \ell(T_l)=4,\qquad \ell(T_l^{k+1})=2l^k\quad(k\ge1).      \tag{16}
\]

The second formula follows for every `k` from the checked base case and
`T_l^3=lT_l^2`.  Combining (13)--(16) yields

\[
 \ell(\operatorname{Re}h(T_l))
 =4\operatorname{Re}\mathcal S_E(l)
 =\ell(\operatorname{Re}E).                                \tag{17}
\]

Substitution into (7) proves

\[
 \boxed{\frac d{d\epsilon}\Big|_{0}
 \frac{\lambda_{\max}P_\epsilon}{\lambda_{\min}P_\epsilon}=0}
\]

for every complex ambient perturbation after conformal pullback.

## 5. Consequence and remaining debt

The L71 disk curve is therefore an exact critical manifold of this analytic
rank-one feasible certificate, not merely a curve on which several Taylor
coefficients happened to vanish.  Together with L70's negative-definite weighted
leading form, this supplies the two ingredients expected in a weighted
Morse--Bott argument.

L73 justifies analytic dependence of the optimized certificate and uses the
negative strong Hessian plus the soft quartic normal form to prove a uniform
sign.  Thus L70--L72 are promoted there to a full neighbourhood theorem around
the single `C_3` block.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/p3_crabb_disk_center.py
.venv/bin/python -u experiments/p3_disk_center_tangent.py
```

Both scripts use exact SymPy identities; the tangent audit completes in about one
second and makes no floating-point sign decision.
