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

A corrected Taylor implementation also identified a separate numerical issue in
the first continuous-$c$ attempt: it enclosed the zeroth Taylor term between
zero and its value instead of retaining the value exactly. Preserving the
Taylor coordinate as a common Bernstein axis removes that dependency loss.
`experiments/slice_projective_interval_certificate.py` now implements the
result with Arb balls for the cancellation-sensitive scalar coefficients and
one-ulp outward binary64 arithmetic for the large Bernstein tensors and every
de Casteljau subdivision. Theta tails and Taylor remainders are enclosed
explicitly. More precisely, if the first omitted theta-series exponent is
$E_M$, its order-$j$ Taylor coefficient is bounded by the first omitted term
times the geometric factor

\[
 \left(1-4^j c_+^{4M+2}\right)^{-1}. \tag{13}
\]

This dominates both exponent sequences used by $\theta_3$ and the auxiliary
series. Arb automatic differentiation propagates those coefficient balls
through $k,s,\gamma_-,\gamma_+-\gamma_-$. On a rational nome box with midpoint
$c_0$ and radius $h$, the implementation retains Taylor orders $0$ through
$9$ as one Bernstein coordinate and encloses the order-$10$ remainder by
$F^{(10)}([c_0-h,c_0+h])h^{10}/10!$. Exact rational change-of-basis constants
are rounded outwards to binary64; the same is true of every subsequent sum,
product, and de Casteljau half-sum. Thus a nonnegative terminal Bernstein box
is an enclosure proof, not a floating-point sign test.

All determinant and final minor charts, both signs, certify on the four-box
cover

\[
 [.01,.012]\cup[.012,.0144]\cup[.0144,.01728]
 \cup[.01728,.020736]=[.01,.020736]. \tag{14}
\]

On the first box, the positive-sign leaf counts are respectively
$4,127,7,129,3,436,144,575$; for the negative sign they are
$1,3,5,2,4$. On the second box the corresponding counts are
$3,049,4,985,2,215,130,384$ and again $1,3,5,2,4$. On the third and fourth
boxes the positive-sign counts are respectively
$2,668,4,215,1,718,88,334$ and $1,932,2,902,1,286,84,272$; the negative-sign
counts remain $1,3,5,2,4$ on both. Thus (RT) is proved on the
continuous nome interval (14), not
merely sampled there. This local certificate is L54. It is not yet a finite
cover of the entire nome range.

The certificate is reproduced by running the script twice, with
`--sign 1` and `--sign -1`, on each pair of exact decimal endpoints.

### The low-nome exceptional hierarchy

The raw scalar powers obscure two exact cancellations at $c=0$. Introduce
analytic deviations by

\[
\begin{aligned}
 k&=c(4+c^2K_2),&s&=1+c^2S_2,\\
 \gamma_-&=c^2(4+c^2L_2),&
 \gamma_+-\gamma_-&=c^4(32+c^2D_2).
\end{aligned} \tag{15}
\]

Cancellation-free theta formulas make all four deviations regular at zero,
with respective values $-16,-4,-44,-272$. Exact coefficient collection then
gives order $c^9$ for both determinant charts and order $c^2$ for all three
final minor charts. This is stronger than the naive monomial orders $7$ and
$0$ and explains the catastrophic cancellation in direct low-nome arithmetic.

The only sharp determinant corner has main-chart coordinates

\[
 u=0,\qquad v=1,\qquad a=1,\qquad b=0. \tag{16}
\]

Let $H$ denote the determinant after division by $c^9$, and let
$\sigma\in\{1,-1\}$ be the retained sign of $b$. In either main chart, make
the first corner scaling

\[
 u=cU,\quad1-v=cV,\quad1-a=cA,\quad |b|=cB.
\]

The exact first exceptional coefficient is

\[
 H=1536c(A+3U)+O(c^2). \tag{17}
\]

On its equality face put $U=cU_2$ and $A=cA_2$. Exact collection through the
next order gives the sum of squares

\[
 H=192c^2\{8A_2+24U_2+3V^2+18(B-2\sigma)^2\}+O(c^3). \tag{18}
\]

For the negative sign the bracket is strict. For the positive sign its sole
remaining ridge is $A_2=U_2=V=0,B=2$. Direct substitution into the exact
orientation-face factorization gives

\[
 H=46080c^4-1622016c^6+O(c^8). \tag{19}
\]

Thus the apparent negative low-nome sliver is another nested equality ridge,
not a negative determinant. Equations (15)--(19), including equality of the
two main charts and cancellation of the envelope variable, are regenerated by
`experiments/slice_low_nome_certificate.py --audit-ridge`. They are L55. A
finite interval enclosure of the remainders in these final coordinates is
still required; L55 alone does not prove a low-nome interval.

Two further exact hierarchies were exposed while enclosing that remainder.
First choose the main variable in the largest-coordinate chart and center its
upper orientation corner. With weighted variables $(C,S,D,B)$ of weights
$(1,2,1,1)$, and the remaining ratio $Q\in[0,1]$, the two determinant charts
start respectively with

\[
16(X+4QS)(X+4S(1-Q)),\qquad
64(1-Q)(X+4S(1-Q)),\quad
X=16C^2+3D^2+3B^2. \tag{20}
\]

For the second chart, centering $P=1-Q$ gives the next forms

\[
16X(X+4S+4P),\qquad
16S(4P+7S)(4P+4S)\quad(X=0). \tag{21}
\]

All displayed factors are nonnegative. For det0, the weighted $S$-chart has
two endpoint zeros; exact affine endpoint neighborhoods certify directly, and
the remaining chart closes in 22 leaves. For det1, (21) identifies the endpoint
hierarchy exactly. In its $P$-dominant chart the next transverse form is

\[
64(16C^2+3D^2+3B^2+4P^2S). \tag{22}
\]

Only the $S$-dominant chart at $P=0$ survives. If $s$ is its radial coordinate,
the final quadratic is

\[
64(1+s^2)\{16C^2+3B^2+3(D+sP)^2+4(1+s^2)P^2\}. \tag{23}
\]

Thus the last det1 face is positive definite in its four transverse variables.
The analytic-tail enclosure on $0\leq c\leq.005$ certifies this complete
hierarchy for both signs: four first endpoint charts; ten centered $P$-models
for each second chart; and four half-box models for each final chart. All 59
charts per sign pass. The positive sign takes about 84 minutes and uses at most
nine leaves; the negative sign takes about 76 minutes and every chart passes
at the root. Peak memory is 11.25 GB without swap. Thus the complete det1 main
chart is proved on that interval for both signs.

The other cancellation-sensitive face occurs in the top `one-minus-ratio`
chart at the already-proved parameter face $a=0$. After centering its three
upper coordinates, give $(Z,U,R,A,B)$ weights $(1,2,2,2,1)$. Its exact leading
form is

\[
16(3B^2+16Z^2)(3B^2+16Z^2+4R+4U). \tag{24}
\]

In the $U$-dominant chart, after replacing the even radial variable by its
square, the transverse leading form is

\[
64(16Z^2+3B^2+3sA^2+4sR). \tag{25}
\]

With the analytic tail enclosure on $0\leq c\leq.01$, the nome-, $R$-, $A$-,
and $B$-dominant charts certify even after all ratios are enlarged to $[0,2]$.
The $U$-axis equality line has one further finite hierarchy. In the
$R$-dominant transverse chart its equality face vanishes identically, and the
next form is

\[
64(1+r^2)\{16Z^2+3B^2+s(4+7r^2+6rA+3A^2)\}. \tag{26}
\]

In the $A$-dominant chart the corresponding form is

\[
64\{(1+ra^2)(16Z^2+3B^2)
 +s(3+4r+6ra+3ra^2+11r^2a^2+6r^2a^3+7r^3a^4)\}. \tag{27}
\]

Both are coefficient-positive on the unit cube and strict in their three
transverse variables. On $0\leq c\leq.005$, 26 centered-model charts certify
the complete positive-sign $U$-dominant hierarchy: four nome, four $B$, six
nested $R$, and twelve nested $A$ boxes, all at the root. Regeneration takes
about 47 minutes and peaks at 12.72 GB without swap. The negative-sign run and
integration of this local chart into the global top chart remain pending.
Equations (20)--(27) are checked as exact rational sparse-map identities by
`--audit-ridge`.

For the analytic tails, exact rational coefficients are retained through
order 15. Since each normalized scalar factor is $c^p$ times an even analytic
function, even and odd records use directed Arb remainder orders 16 and 17
respectively. This preserves the parity needed by the weighted charts. Only
the genuine remainder is interval-valued; exact Taylor cancellations are no
longer performed between Arb balls.

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

Extend L54's four-box adaptive cover across the compact range, keeping $c$ as
a shared Bernstein coordinate. At $c=0$, regenerate the det0 $U$-hierarchy
for the negative sign and integrate its five local corner charts into the
global top chart. Then certify the three low minor charts and bridge
$[.005,.01]$. Completing all three tasks proves (RT) and therefore the elliptic
$4\times4$ slice, but not the general Crouzeix conjecture.
