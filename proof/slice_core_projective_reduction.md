# Projective reduction of the remaining elliptic transfer core

**Status (2026-07-22):** the exact reductions, both orientation faces, the
complete continuous interval $0\leq c\leq.020736$, and the high-nome range
$c\geq2^{-2/3}$ are proved. In the compact interval, the negative sign is
certified through $.63$ and the positive sign through
$c_+=.5566585294072849\ldots$. Thus only the positive-sign window
$c_+<c<2^{-2/3}$ remains open.

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

For the positive sign, the analytic tail enclosure on $0\leq c\leq.01$
certifies the nome-, $R$-, $A$-, and $B$-dominant charts even after all ratios
are enlarged to $[0,2]$. The four negative-sign counterparts have now also
been regenerated on $0\leq c\leq.005$ and pass at the root.
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
the complete $U$-dominant hierarchy for both signs: four nome, four $B$, six
nested $R$, and twelve nested $A$ boxes, all at the root. The negative run took
$3032.81$ seconds (50m33s), peaked at 11.02 GB, and used no swap. Equations
(20)--(27) are checked as exact rational sparse-map identities by
`--audit-ridge`.

The Cartesian cutoff is not intrinsic. The exact transverse quadratic (25)
has no zero at $s=1/2$; that apparent transition was only a centered-model
seam. Two scale-free arm charts remove the cutoff completely. First, a
nome-dominant chart allows each remaining normalized ratio up to $1/2$ and
certifies the entire nome radial interval for both signs. Second, the
$U$-dominant chart retains its radial variable in exact degree-73 Bernstein
form, factors the four transverse variables, and models only their selected
radial and ratio coordinates. Twenty boxes cover the complete positive-sign
largest-coordinate chart; for the negative sign the same boxes cover the
ratio cap $1/2$, which is enough for the equality tube while the strict
complement remains global. All local boxes pass (all negative boxes at the
root).

After admitting the proved Cartesian corner and the two arm regions, the top
det0 global complement certifies in 235 leaves at depth 8 for the positive
sign and 399 leaves at depth 9 for the negative sign. Hence the formerly open
det0 `one-minus-ratio` chart, and therefore every low determinant chart, is
proved for both signs on $0\leq c\leq.005$.

### Complete low-nome minor certificates

The three normalized final minors have two sharp mechanisms: the original
corner and the midpoint face where the main coordinate is $1/2$. At the first
corner, use $C$ for the scaled nome. For the secondary minor put
$U=1-u$, $V=v$, $Y=1-y$, and $A=a$. Its exact transverse quadratic is

\[
 256C^2+64U^2+64UV+64UY+48V^2+48A^2. \tag{28}
\]

For the tertiary minors put $R=1-r$ and retain $A$ for the lower ratio. The
common form is

\[
 256C^2+64U^2+64UV+48V^2+96VR+96R^2+96RA+48A^2, \tag{29}
\]

and tertiary chart one has the additional nonnegative term $96VA$. The sole
secondary equality line after selecting $Y$ has the positive-definite
transverse form

\[
 256C^2+64U+48V^2+48A^2. \tag{30}
\]

Two more secondary boxes are needed to cover the ends of that line. Their
leading forms are respectively

\[
 256C^2+48V^2+48A^2+64U(1-q),\qquad
 256C^2+64U+48V^2+64q+48A^2. \tag{31}
\]

At the midpoint write $M=u-1/2$ and $S=1-v$. With the remaining chart
coordinates denoted as in the checker, the three exact forms are

\[
\begin{aligned}
 H_{\rm sec}&=144(A-2\sigma C)^2+144qC^2+384M^2+(64+96q)S,\\
 H_{\rm ter0}&=144(T+R-2\sigma C)^2+144C^2+384M^2+160S,\\
 H_{\rm ter1}&=144C^2+(1-R)
 \{144(T-2\sigma C)^2+384M^2+32(5-2R)S\}.
\end{aligned} \tag{32}
\]

The positive secondary chart has one further ridge, $q=0,A=2C$. In the actual
$c\leq1/200$ scaling its leading form is

\[
 {3\over12500000}Z^2+24M^2+4S+{9\over2500}q+9D^2, \tag{33}
\]

where $Z,D$ are the radial and $A-2C$ blow-up coordinates. Tertiary chart one
also degenerates when $R=1$; its separate ratio-zero chart has leading form
$144(T-2\sigma Cr)^2+144C^2+384M^2+(96+64r^2)S$.

The exact audit regenerates (28)--(33) for both signs. The interval checker
then certifies all local charts and the global complement on
$0\leq c\leq.005$. The global leaf counts for the secondary, tertiary-zero,
and tertiary-one minors are respectively $81,32,184$ for the positive sign and
$71,32,188$ for the negative sign, at maximum depths four, four, and five.
Complete self-contained regeneration took 10m19s and 8.13 GB for the positive
sign, and 3m52s and 5.77 GB for the negative sign, with no swap:

```text
slice_low_nome_certificate.py 0 .005 --certify-asymptotic-minors --sign 1
slice_low_nome_certificate.py 0 .005 --certify-asymptotic-minors --sign -1
```

For memory control, an irrelevant envelope axis is converted to Bernstein
form and collapsed before the remaining axes. This is rigorous: the envelope
control balls are replaced by their convex hull, and every subsequent
power-to-Bernstein conversion is a positive linear map, so it preserves that
inclusion. It is not a sampled or floating-point shortcut.

For the analytic tails, exact rational coefficients are retained through
order 15. Since each normalized scalar factor is $c^p$ times an even analytic
function, even and odd records use directed Arb remainder orders 16 and 17
respectively. This preserves the parity needed by the weighted charts. Only
the genuine remainder is interval-valued; exact Taylor cancellations are no
longer performed between Arb balls.

These degeneracies are not evidence of a negative determinant. Factoring the
arms, rather than accepting tiny negative interval bounds or moving a cutoff,
is exactly what closes the low chart. Combining the determinant result above
with the three complete minor certificates proves the full polynomial core on
$0\leq c\leq.005$ for both signs.

### The bridge and compact directed covers

The generic shared-coordinate checker closes the former bridge on the exact
four-box cover

\[
 [.005,.006]\cup[.006,.0072]\cup[.0072,.00864]\cup[.00864,.01]. \tag{34}
\]

Every one of the two determinant and three final-minor charts passes for both
signs.  For the positive sign, the first-box leaf counts are respectively
$12377,20230,9192,214,1757$; the other three boxes require
$(9071,14337,6388,196,1308)$,
$(6850,11753,4966,160,1104)$, and
$(5351,7554,3770,148,695)$.  The negative sign needs only the root or one
subdivision.  Hence L55, (34), and L54 prove (RT) continuously through
$c=.020736$; (34) is L56.

On the compact range, geometric rational boxes give the current rigorous
frontier.  The negative sign passes from $324/15625$ to $.63$ (coarse ratio
$6/5$, followed by ratio $11/10$).  The positive sign passes successive
ratio-$6/5,11/10,21/20,41/40,81/80$ covers through

\[
 c_+=
 {7471344308886696360166308338925321050191281742307131164807\over
  13421772800000000000000000000000000000000000000000000000000}.
 \tag{35}
\]

These are continuous interval certificates, not fixed-$c$ samples.  Since
$(63/100)^3=250047/10^6>1/4$, the negative cover overlaps L42.  Only the
positive-sign tail after (35) remains.

### Exact factorization of the sharp positive face

The remaining positive determinant stalls are confined numerically to
$a,b\to1$.  This edge is not a zero or a sign obstruction.  Let $X,R,Y$ denote
the radial, ratio, and envelope coordinates, and put

\[
 q=\gamma_-+Y(\gamma_+-\gamma_-),\quad C=s+Xq. \tag{36}
\]

Exact collection at $a=b=1$ gives, in determinant chart zero,

\[
 D_0=81c^9(1-k)^2(1-kX)^2
       [R+(1-RX)C^2]^4. \tag{37}
\]

In chart one, with $\widetilde C=s+RXq$,

\[
 D_1=81c^9(1-k)^2(1-kRX)^2
       [1+R(1-X)\widetilde C^2]^4. \tag{38}
\]

Both are manifestly nonnegative.  More is true.  Write $A=1-a,B=1-b$.
The coefficients of $A$ and $B$ are identical.  With $Z=X$ in chart zero
and $Z=RX$ in chart one, their common value is

\[
 54c^9(1-k)(1-kZ)
 [2+3k+3Zk-8Zk^2]\,{\cal B}^4, \tag{39}
\]

where ${\cal B}$ is the bracket in (37) or (38).  The remaining factor is
nonnegative: if $k\leq3/8$, minimize in $Z$ at $Z=0$; if $k\geq3/8$, minimize
at $Z=1$ and obtain $2(1-k)(4k+1)$.  Thus moving inward from the face first
*increases* the determinant.  This explains why raw binary64 subdivision was
stalling on a positive edge.

`experiments/slice_positive_face_audit.py` regenerates (37)--(39) from all
$197563$ integer records in each chart.  It also collects every term of total
$(A,B)$-degree at least two and obtains the same exact coefficient l1 norm
$41235531913$ in both charts.  The remaining certificate task is to retain
the correlated nome/spatial cancellation in that higher-order remainder;
using the global l1 norm directly is rigorous but much too coarse.

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

Close the positive-sign tail (35) to $2^{-2/3}$ by using the exact face value
and positive first variation (37)--(39), while keeping the higher-order
deficit coefficients correlated in the nome and the three surviving
projective variables.  Do not accept a tolerance or repeat unbounded raw
subdivision.  Completing this tail proves (RT), hence the complete elliptic
$4\times4$ slice; it still does not resolve the general Crouzeix conjecture.
