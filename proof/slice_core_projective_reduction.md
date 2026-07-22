# Projective reduction of the remaining elliptic transfer core

**Status (2026-07-21):** exact reduction proved; global positivity is still open.
The fixed-parameter Bernstein results below are falsification evidence, not a
certificate over the continuous nome parameter.

## 1. The two minors that remain

Use the polynomial transfer core from `slice_coupled_defects.md` (27):

\[
 Q=\begin{bmatrix}E_o&-3J\\-3J^T&E_e\end{bmatrix}.
\]

In the open parameter square, the already-proved positivity of (E_o) is
strict.  Its Schur complement is the real symmetric (2\times2) matrix

\[
 S=E_e-9J^TE_o^{-1}J.
\]

Consequently (Q\succeq0) is equivalent to only

\[
 \det(E_o)S_{11}\ge0,\qquad \det Q\ge0. \tag{1}
\]

The first expression is one leading (3\times3) principal minor.  This is an
exact Sylvester/Schur reduction; all other leading minors are already supplied
by (E_o\succ0).  Parameter edges follow by continuity from L41/L48/L50.

## 2. Cancellation-free modal coordinates

Put (T=\tan^2u\) and (N=(1+T)(1+r^2T)).  The common modal matrix is

\[
 M={\sqrt{k}\over\sqrt N}
 \begin{bmatrix}
  1+prT&-(p-r)\sqrt T\\
  (1-pr)\sqrt T&p+rT
 \end{bmatrix}. \tag{2}
\]

With (D=\operatorname{diag}(1,c)),

\[
 B=c^{-1/2}DMD^{-1},\qquad C=c^{1/2}DM^TD^{-1}. \tag{3}
\]

Substitution of (2)--(3) into (1), followed only by multiplication by positive
powers of (N,c,1-a^2,1-b^2), gives two integer polynomials.  Simultaneously
changing the signs of (a,b) is a block congruence, so it is enough to take
(a\ge0) and treat (b\ge0) and (b\le0) separately.

## 3. The projective cancellation

Compactify orientation by (o=T/(1+T)), and set

\[
 P=p^2,\qquad t=1-o.
\]

Write the full rigorous cubic envelope as

\[
 r=p(s+\gamma P),\qquad
 \gamma=(1-y){s(1+k^2-s^2)\over6}+y(1-s),\quad0\le y\le1. \tag{4}
\]

After (4), every power of (p) in both normalized minors is even.  Exact
coefficient collection gives

\[
 \operatorname{ord}_{(P,t)} Q_3=3,qquad
 \operatorname{ord}_{(P,t)}\det Q=4. \tag{5}
\]

Thus two polynomial charts cover the complete ((P,t))-square:

\[
 (P,t)=(w,wv),\qquad (P,t)=(wv,w),\qquad0\le w,v\le1, \tag{6}
\]

after division by (w^3) for (Q_3) and (w^4) for (det Q).  No square
roots or limiting numerical cancellations remain.  This is the correct
blow-up of the intersection (p=0,o=1); raw boxes converge indefinitely onto
that intersection.

## 4. Evidence and what it does not prove

For each fixed

\[
 c\in\{.001,.003,.01,.03,.05,.1,.3,.6,.629\},
\]

the two charts in (6), with (y,a,|b|\) retained as free Bernstein variables,
certified both signs of (b) for both expressions in (1).  Typical positive-
sign determinant charts required at most 277 leaves; the negative-sign charts
usually required at most three.  A separate 100,000-sample sign-variation
audit also found the (3\times3) minor increasing through the cubic envelope
and no determinant minimum below its envelope endpoints.  Neither finite
observation is used as a theorem.

Directed interval/Taylor enclosure in (c) was then attempted.  It correctly
refused to certify boxes meeting further exact zero intersections:

* in the second chart, (v=0, |b|=1) is a zero of the (3\times3) minor;
* at (P=1,o=0), the determinant is (c^7F(a,b,c,k)^2), so its zero curve
  moves with (c);
* the coordinate axes and (p=0) face are already closed by L48/L50, but
  their intersections still appear as zero Bernstein coefficients.

These are genuine algebraic zeros, not evidence of a negative determinant.
Any rigorous continuation must factor or blow up these intersections, or use
an exact moving coordinate such as (F).  It must not accept a small negative
interval bound as rounding error.

## 5. Shortcuts falsified during this reduction

Three tempting simplifications fail numerically and should not be retried:

1. Replacing every block operator norm in (BE) by its Frobenius norm reaches
   (2), not (1), as (c\to0,p\to1,a,b\to0).
2. Replacing each block norm by the rational trace/determinant majorant used in
   L47 still reaches (3/2).
3. The block energy is not maximized at the two orientation endpoints; sampled
   interior improvements over both endpoints are large.

The exact polynomial core, not another uncoupled norm majorant, is therefore
the preferred route.

## 6. Next exact target

Center the second chart simultaneously at (v=0) and (1-|b|=0), divide the
resulting common first-order factor, and rerun the directed (c)-certificate
away from (P=1,o=0).  Treat that last corner in coordinates

\[
 1-P,quad o,quad F(a,b,c,k),
\]

using the exact square on its face.  Success for these charts proves (RT) and
therefore completes the elliptic (4\times4) slice, but not the general
Crouzeix conjecture.
