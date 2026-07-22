# Projective reduction of the remaining elliptic transfer core

**Status (2026-07-21):** exact reductions and the two orientation faces are
proved; global interior positivity is still open. The fixed-parameter
Bernstein results below are falsification evidence, not a certificate over the
continuous nome parameter.

## 1. The two minors that remain

Use the polynomial transfer core from `slice_coupled_defects.md` (27):

\[
 Q=\begin{bmatrix}E_o&-3J\\-3J^T&E_e\end{bmatrix}.
\]

In the open parameter square, the already-proved positivity of $E_o$ is
strict. Its Schur complement is the real symmetric $2\times2$ matrix

\[
 S=E_e-9J^TE_o^{-1}J.
\]

Consequently, $Q\succeq0$ is equivalent to only

\[
 \det(E_o)S_{11}\ge0,\qquad \det Q\ge0. \tag{1}
\]

The first expression is one leading $3\times3$ principal minor. This is an
exact Sylvester/Schur reduction; all other leading minors are already supplied
by $E_o\succ0$. Parameter edges follow by continuity from L41/L48/L50.

## 2. Cancellation-free modal coordinates

Put $T=\tan^2u$ and $N=(1+T)(1+r^2T)$. The common modal matrix is

\[
 M={\sqrt{k}\over\sqrt N}
 \begin{bmatrix}
  1+prT&-(p-r)\sqrt T\\
  (1-pr)\sqrt T&p+rT
 \end{bmatrix}. \tag{2}
\]

With $D=\operatorname{diag}(1,c)$,

\[
 B=c^{-1/2}DMD^{-1},\qquad C=c^{1/2}DM^TD^{-1}. \tag{3}
\]

Substitution of (2)--(3) into (1), followed only by clearing positive
factors involving $N,c,1-a^2,1-b^2$, gives two integer polynomials.
Simultaneously changing the signs of $a,b$ is a block congruence, so it is
enough to take $a\ge0$ and treat $b\ge0$ and $b\le0$ separately.

## 3. The projective cancellation

Compactify orientation by $o=T/(1+T)$, and set

\[
 P=p^2,\qquad t=1-o.
\]

Write the full rigorous cubic envelope as

\[
 r=p(s+\gamma P),\qquad
 \gamma=(1-y){s(1+k^2-s^2)\over6}+y(1-s),\quad0\le y\le1. \tag{4}
\]

After (4), every power of $p$ in both normalized minors is even. Exact
coefficient collection gives

\[
 \operatorname{ord}_{(P,t)} Q_3=3,\qquad
 \operatorname{ord}_{(P,t)}\det Q=4. \tag{5}
\]

Thus two polynomial charts cover the complete $(P,t)$-square:

\[
 (P,t)=(w,wv),\qquad (P,t)=(wv,w),\qquad0\le w,v\le1. \tag{6}
\]

Divide by $w^3$ for $Q_3$ and by $w^4$ for $\det Q$. No square roots or
limiting numerical cancellations remain. This is the correct blow-up of the
intersection $p=0,o=1$; raw boxes converge indefinitely onto that
intersection.

## 4. The two orientation faces are closed

The apparent square obstruction at $P=1,o=0$ is not an interior zero. It can
be resolved without a certificate. For $0\le x\le k$, define

\[
\begin{aligned}
 G_x(a,b)={}&4a^2b^2c^2x-16a^2b^2cx^2-a^2b^2c+4a^2b^2x\\
 &-4(a^2+b^2)c^2x+4(a^2+b^2)cx^2
   +4(a^2+b^2)c-4(a^2+b^2)x\\
 &+18abcx+4c^2x-cx^2-16c+4x,\\
 H_x(a,b)={}&4a^2b^2cx-a^2b^2-a^2cx+a^2
              -4b^2cx+4b^2+cx-4.
\end{aligned} \tag{7}
\]

The elementary nome bound gives $x\le k<4c$, and $x<1$. Then

\[
 \boxed{G_x(a,b)<0,\qquad H_x(a,b)\le0.} \tag{8}
\]

For $H_x$, regard it as bilinear in $a^2,b^2$. Its four corner values are
$cx-4,-3,-3cx,0$, so it is nonpositive.

For $G_x$, changing the sign of exactly one of $a,b$ only decreases it. Take
$a,b\ge0$, put $q=ab$, and $z=a^2+b^2$. Then

\[
 G_x=C_0+C_1z+C_2q^2+18cxq,
\]

where

\[
 C_0=-(4c-x)(4-cx)<0,\quad
 C_1=4(c-x)(1-cx),\quad
 C_2=4c^2x-16cx^2-c+4x. \tag{9}
\]

The elementary bounds $2q\le z\le1+q^2$ finish the proof. If $x\le c$,
then $C_1\ge0$, and replacing $z$ by $1+q^2$ gives

\[
 G_x(a,b)\le G_x(1,q)
 =-3c(2qx-q+x-2)(2qx+q-x-2)\le0. \tag{10}
\]

If $x\ge c$, replace $z$ by $2q$ and put $r=1-q$. The resulting upper
bound is

\[
 -9c(1-x)^2+6c(4x+1)(x-1)r+C_2r^2. \tag{11}
\]

The first two terms are nonpositive. If $C_2\le0$, this is immediate; if
$C_2\ge0$, the expression is convex in $r$, so its maximum is at an endpoint,
where it equals either $-9c(1-x)^2$ or $C_0$. Both are strictly negative in
the nondegenerate range. This proves (8).

For a cancellation-free statement of the endpoint factorization, apply the
positive diagonal congruence

\[
 \widetilde Q=\operatorname{diag}(I,\sqrt c I)\,
 Q\,\operatorname{diag}(I,\sqrt c I).
\]

At $o=0$, exact substitution factors its two expressions in (1) as

\[
 \widetilde Q_3=(1-a^2)(1-b^2)H_{kp^2}G_k,\qquad
 \det\widetilde Q=(1-a^2)^2(1-b^2)^2G_kG_{kp^2}. \tag{12}
\]

At $o=1$, the modal labels interchange, so the leading minor contains
$H_kG_{kp^2}$ and the determinant is unchanged. Hence both expressions are
nonnegative on both orientation faces. Together with the positive diagonal
block, this proves (RT) at $o=0$ and $o=1$. Boundary values of $a,b$ follow by
continuity (and were already closed before L52).

`experiments/slice_orientation_factor_audit.py` regenerates the scalar core
determinant and every polynomial identity in (7)--(12) with exact symbolic
arithmetic.

In particular, after removal of the displayed positive defect factors, the
determinant at $p=1,o=0$ is the square $G_k^2$. It has no moving interior zero
curve. Its small margin near degenerate nome limits explains why a
low-precision interval box looked singular.

## 5. Evidence and what it does not prove

For each fixed

\[
 c\in\{.001,.003,.01,.03,.05,.1,.3,.6,.629\},
\]

the two charts in (6), with $y,a,|b|$ retained as free Bernstein variables,
returned nonnegative subdivisions for both signs of $b$ and both expressions
in (1). Typical positive-sign determinant charts required at most 277 leaves;
the negative-sign charts usually required at most three. A separate
100,000-sample sign-variation audit also found the $3\times3$ minor increasing
through the cubic envelope and no determinant minimum below its envelope
endpoints. Neither finite observation is used as a theorem.

Directed interval/Taylor enclosure in $c$ was then attempted. It correctly
refused to certify boxes meeting further degenerate intersections:

* in the second chart, $v=0,|b|=1$ is a zero of the $3\times3$ minor;
* at $P=1,o=0$, the determinant is, up to known positive factors, the strictly
  positive square $G_k(a,b)^2$, but its margin degenerates as $c\downarrow0$ and at the
  disk/coalescence limits;
* the coordinate axes and $p=0$ face are already closed by L48/L50, but their
  intersections still appear as zero Bernstein coefficients.

The first intersection begins with an exact order-one secondary blow-up. Put
$|b|=1-q$ in the second chart and use

\[
 (v,q)=(w,wh),\qquad(v,q)=(wh,w),
\]

after division by $w$. Both charts pass every fixed-$c$ test above, for both
signs of $b$ and the complete cubic envelope. This is still finite evidence
until the continuous-$c$ enclosure is regenerated in these coordinates.

The first secondary chart retains one smaller exact intersection. On its
$w=0$ face, setting the secondary ratio $h=0$ leaves a factor $1-a^2$; in
particular, the common zero is $(h,1-a)=(0,0)$. Put $z=1-a$ and apply

\[
 (h,z)=(w,wj),\qquad(h,z)=(wj,w),
\]

again dividing by $w$. Exact coefficient collection gives order one. The two
tertiary charts pass the same complete fixed-$c$ tests. For the positive sign,
the worst low-nome test at $c=.001$ uses 1,408 leaves; the negative sign uses
at most four. Thus the previously observed interval stall is explained by a
specific nested boundary intersection, not a negative minor.

`experiments/slice_projective_core.py` regenerates the complete sparse
polynomials from (2)--(4), verifies the main orders $3,4$, and verifies both
subsequent order-one blow-ups exactly. The cache is deliberately generated and
git-ignored rather than checked in.

A corrected Taylor prototype also identified a separate numerical issue in
the first continuous-$c$ attempt: it enclosed the zeroth Taylor term between
zero and its value instead of retaining the value exactly. Preserving the
Taylor coordinate as a common Bernstein axis removes that dependency loss.
With this correction, all determinant charts and the secondary/tertiary minor
charts certify on the continuous test intervals $[.01,.0101]$ and
$[.1,.101]$. This remains evidence until the directed-rounding implementation
and a finite cover of the entire nome range are committed and audited.

These degeneracies are not evidence of a negative determinant. Any rigorous
continuation must factor or blow them up and must not accept a small negative
interval bound as rounding error.

## 6. Shortcuts falsified during this reduction

Three tempting simplifications fail numerically and should not be retried:

1. Replacing every block operator norm in (BE) by its Frobenius norm reaches
   $2$, not $1$, as $c\to0,p\to1,a,b\to0$.
2. Replacing each block norm by the rational trace/determinant majorant used in
   L47 still reaches $3/2$.
3. The block energy is not maximized at the two orientation endpoints; sampled
   interior improvements over both endpoints are large.

The exact polynomial core, not another uncoupled norm majorant, is therefore
the preferred route.

## 7. Next exact target

Build the directed $c$-certificate from the reproducible record generator,
keeping $c$ as a shared Bernstein coordinate. The compact range away from
$c=0$ is now locally certified. At $c=0$, normalize
$k/c$, $\gamma_-/c^2$, and $(\gamma_+-\gamma_-)/c^4$ before applying a final
asymptotic chart cover; raw determinant records have two additional cancelling
powers of $c$. Success proves (RT) and therefore completes the elliptic
$4\times4$ slice, but not the general Crouzeix conjecture.
