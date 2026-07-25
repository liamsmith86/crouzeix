# The one-sided elliptic model has an exact rank-\(m\) Stein factor

## 1. Result (L240, 2026-07-24)

Let \(\mathcal H_+=\ell^2(\mathbb N_0;\mathbb C^m)\), let \(L^*\)
be the backward shift, and let
\(J_n:\mathbb C^m\to\mathcal H_+\) be the \(n\)-th coordinate
isometry.  Put

\[
S_\infty=L^*\otimes I_m,\qquad
E_\infty=J_0J_0^*.
\]

Thus

\[
S_\infty S_\infty^*=I,\qquad
I-S_\infty^*S_\infty=E_\infty.
\]

Work in formal power series at \(c=0\), put \(q=c^2\), and define the
right-balanced ellipse pencil

\[
\Xi_\infty(c)
=S_\infty+cS_\infty^*(I+E_\infty),                \tag{1}
\]

and let

\[
A_\infty(c)=\phi_c(\Xi_\infty(c)),
\]

where \(\phi_c\) is the centered normalized Riemann map of the ellipse
parameterized by \(z=\zeta+c/\zeta\).  Here the expression is defined
coefficientwise from the finite sections below.  Every coefficient is
a finite-band operator.  This formal scope is exactly what the
associated-graded L228 calculation needs; no boundary \(H^\infty\)
calculus is being assumed.

Put

\[
D_R(q)=\operatorname {diag}
\left(1,\frac1{1+q},\frac1{1+q^2},\ldots\right)
\otimes I_m.                                      \tag{2}
\]

Then the complete one-sided Stein slack has formal rank \(m\):

\[
\boxed{
D_R-A_\infty^*D_RA_\infty=d_\infty d_\infty^*.}   \tag{3}
\]

Here \(d_\infty:\mathbb C^m\to\mathcal H_+\) is the scalar Hardy
column

\[
\boxed{
d_\infty(c)
=\frac1{\vartheta_3(c^2)}
\left\{
J_0
+2\sum_{j\geq1}
\frac{(-c)^j}{1+c^{4j}}J_{2j}
\right\},}                                        \tag{4}
\]

The explicit column series in (4) also converges in operator norm for
\(|c|<1\), although only its formal identity is used below.

Since \(J_0^*d_\infty=\vartheta_3(c^2)^{-1}I_m\) is invertible, the
Schur complement of (3) away from the right defect is exactly zero:

\[
\boxed{
Q_\infty H_\infty Q_\infty
-Q_\infty H_\infty J_0
 (J_0^*H_\infty J_0)^{-1}
 J_0^*H_\infty Q_\infty=0,}                       \tag{5}
\]

where \(H_\infty=D_R-A_\infty^*D_RA_\infty\) and
\(Q_\infty=I-E_\infty\).

Equation (5) is the exact half-line cancellation required after
L237--L239.  It explains why retained-block calculations alone can
look wrong: the Stein product also includes amplitude that exits into
the eliminated clean chain.  On the full one-sided model that flux is
part of the rank-\(m\) Gram (3), and its Schur residual vanishes.

L240 still does **not** prove L228.  A finite left endpoint changes
the half-line model by L239's reflection
\(\lambda_r(t)=t^r/(1+t^r)\).  One must compute the first active
coefficient of the **exact reflected transform** of (5), including
the state-lift columns, and show that its first active Hardy cell is

\[
E_1F_r+F_rE_1.
\]

The bulk term is now closed exactly; only that reflected boundary
coefficient remains.  One must not replace the exact transform by
its term linear in \(\lambda_r\): depending on \(r\), higher powers
can still reach the requested coefficient, and L235 already warned
that terminal-crossing truncations require proof.

## 2. Derivation from the finite elliptic Crabb axis

Let \(S_N\) be the backward shift on
\(\mathbb C^{N+1}\otimes\mathbb C^m\).  Its two defect projections
are the first and last coordinate blocks.  The balanced finite
ellipse pencil is

\[
\Xi_N
=S_N+c(I+F_N)S_N^*(I+E_N).                        \tag{6}
\]

L117's arbitrary-size elliptic Crabb-axis theorem constructs an exact
physical metric \(P_N(c)\) with rank-\(m\) Stein slack for
\(\phi_c(\Xi_N)\).  In the axis notation, with
\(\ell=-\log c\), its diagonal entries are

\[
p_{N,n}(c)
=\frac{
\sum_{k\in\mathbb Z}
\operatorname {sech}((n+2Nk)\ell)
}{
\sum_{k\in\mathbb Z}
\operatorname {sech}(2Nk\ell)
}\,c^{-n}.                                        \tag{7}
\]

For every fixed \(n\),

\[
\begin{aligned}
p_{N,0}(c)&\longrightarrow1,\\
p_{N,n}(c)&\longrightarrow
\operatorname {sech}(n\ell)c^{-n}
=\frac2{1+c^{2n}}\qquad(n\geq1).                  \tag{8}
\end{aligned}
\]

The repeated equality metric is \(1\) at the right defect and \(2\)
on every fixed interior coordinate.  Passing from physical to
balanced coordinates therefore turns (8) precisely into (2).

For every fixed Taylor degree in \(c\) and every fixed finite
coordinate window, (6) stabilizes to (1).  Indeed, L125's direct-map
coefficient of that degree is a finite polynomial in the banded
pencil, so once the remote left endpoint is farther away than its
maximum path length it cannot enter the coefficient.  The same
finite-path stabilization applies to the Stein product.  We therefore
define \(A_\infty(c)\) coefficientwise by these stable Taylor
coefficients and pass coefficientwise in L117's finite rank-\(m\)
Stein identity.

Since this stabilization works at every degree, the result is an
all-order formal identity, not an extrapolation from any fixed jet.
This proves (3) once the limiting defect column is identified.

## 3. The limiting defect column

The scalar Szegő factor in L117 is the normalized Jacobi
\(\operatorname {nd}\) Fourier column.  Its standard nome expansion,
in the present \(c\)-normalization, has coordinate coefficients

\[
\begin{aligned}
[d_\infty]_0&=\vartheta_3(c^2)^{-1},\\
[d_\infty]_{2j}
&=\frac{2(-c)^j}
 {\vartheta_3(c^2)(1+c^{4j})}\qquad(j\geq1),\\
[d_\infty]_{2j+1}&=0.                             \tag{9}
\end{aligned}
\]

This is (4).  Absolute summability follows geometrically from
\(c^j/(1+c^{4j})\leq c^j\).  The finite DCT-I defect columns in
L117 converge coefficientwise to (9), so their rank-\(m\) Grams
converge to \(d_\infty d_\infty^*\).  This completes the limiting
proof of (3).

The first terms are

\[
\begin{aligned}
d_\infty(c)
={}&(1-2c^2+4c^4-\cdots)J_0\\
&+(-2c+4c^3-6c^5+\cdots)J_2\\
&+(2c^2-6c^4+\cdots)J_4+\cdots ,
\end{aligned}                                     \tag{10}
\]

agreeing with the independently derived finite-axis columns in
L217.

## 4. Schur residual

Equation (3) is a Gram with invertible first block

\[
J_0^*d_\infty
=\vartheta_3(c^2)^{-1}I_m.
\]

For any block column \(d=\binom{d_0}{d_1}\) with invertible \(d_0\),

\[
d_1d_1^*
-d_1d_0^*(d_0d_0^*)^{-1}d_0d_1^*=0.
\]

Apply this identity to the splitting
\(E_\infty\mathcal H_+\oplus Q_\infty\mathcal H_+\).  This proves
(5) without a coefficient expansion.

## 5. Exact formal audit

The proof above is analytic and all-order.  An independent checker
works in the coisometric word algebra

\[
SS^*=I,\qquad E=I-S^*S.
\]

It generates the full L125 theta/ODE direct map, the metric (2), and
the column (4) as rational truncated series.  Through degree eight,
every coefficient of

\[
D_R-A_\infty^*D_RA_\infty-d_\infty d_\infty^*
\]

reduces exactly to zero before any floating evaluation.

Regenerate with

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_half_line_stein_factor.py \
  --output \
  experiments/repeated_crabb_half_line_stein_factor_s70224.jsonl
```

The tracked SHA-256 is
`4ff399c6bf44361f56be24b6c9cd29c932529fe637eb8b9e7b3ff0a28a93f4eb`.
