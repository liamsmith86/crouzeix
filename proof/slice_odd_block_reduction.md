# Sign-separated coupled faces and the odd-Blaschke block

This note identifies and proves the exact scalar theorem hidden in the
generic rank-one/full KKT face of `slice_coupled_defects.md`.  It turns the
coupled dual matrix problem into one explicit odd finite-Blaschke block,
removes the orientation and Blaschke parameters, and certifies the remaining
two-variable endpoint inequalities.  By itself it does **not** close the
separate rank-one/rank-one face; the downstream L40--L59 chain now closes that
face and proves L20.

## 1. A rank-one boundary is an orthogonal colligation

Write

\[
 T=\begin{pmatrix}0&B\\ C&0\end{pmatrix},\qquad F=BC.
\]

Consider the sign-separated boundary of a coupled dual face

\[
 Z_o=xx^T,\qquad
 D_o=Z_o-BZ_eB^T=-uu^T,\qquad
 D_e=Z_e-CZ_oC^T=yy^T. \tag{1}
\]

Eliminating \(Z_e=Cxx^TC^T+yy^T\) gives

\[
 xx^T+uu^T=(Fx)(Fx)^T+(By)(By)^T. \tag{2}
\]

Hence the two matrices with columns \([x,u]\) and \([Fx,By]\) have the
same row Gramian.  In the nondegenerate case they differ by a real orthogonal
matrix.  Writing one of its entries as \(a\in[-1,1]\), and eliminating the
other entry, gives

\[
 u={F-aI\over\sqrt{1-a^2}}x,\qquad
 y={aC-B^{-1}\over\sqrt{1-a^2}}x. \tag{3}
\]

The reflection component of \(O(2)\) gives the same quotient.  Therefore

\[
 {\operatorname{tr}(D_o)_-\over\operatorname{tr}(D_e)_+}
 ={\|u\|^2\over\|y\|^2}
 \le \left\|B\,R_a(CB)\right\|^2,\qquad
 R_a(X)=(X-aI)(I-aX)^{-1}. \tag{4}
\]

Conversely, every generalized Rayleigh quotient on the right of (4) is
obtained from (2), so the supremum over this boundary is exactly

\[
 \boxed{\sup_{-1\le a\le1}\|B R_a(CB)\|^2.} \tag{5}
\]

At a generic rank-one/full KKT point the smallest and largest eigenvalues of
the primal metric are simple.  Since the dual stationarity matrix is block
diagonal and has one positive and one negative eigenvalue, they lie in
opposite parity blocks, giving precisely (1), up to swapping parity.  Multiple
sandwich eigenvalues are limiting cases.  Thus the theorem proved below
closes the generic rank-one/full face and its closure.

## 2. The quotient is an odd finite-Blaschke block

Use the modal form

\[
 B=c^{-1/2}S_o\Sigma S_e^{-1},\qquad
 C=c^{1/2}S_e\Sigma S_o^{-1},
 \quad S_o=DU,\ S_e=DV.
\]

The push-through identity for rational functions gives

\[
 B R_a(CB)=c^{-1/2}S_o
 \operatorname{diag}(q_1,q_2)S_e^{-1}, \tag{6}
\]

where

\[
 q_i=\tau_i{\tau_i^2-a\over1-a\tau_i^2}. \tag{7}
\]

Equivalently, (6) is the upper parity block of

\[
 f_a(T),\qquad f_a(z)=z{z^2-a\over1-az^2}. \tag{8}
\]

For \(a\in[-1,1]\), \(f_a\) is an odd finite Blaschke product (with the
endpoint cancellations understood), and \(|q_i|\le\tau_i\).  Thus the live
rank-one/full target is the concrete scalar-looking inequality

\[
 \boxed{\|B R_a(CB)\|\le2\quad(-1\le a\le1).} \tag{OB}
\]

The swapped face gives the analogous lower block \(C R_a(BC)\).  The lower
block is closed directly later in this section, and the upper inequality
(OB) is proved by (13at)--(13bg).  In conceptual terms, (OB) is the complete
norm-two theorem for the odd two-node Pick sector.

## 3. The orientation variable disappears exactly

Put \(p=\tau_2/\tau_1\), \(r=H(p)\), and \(x=\tan u\), as in the modal
notes.  For arbitrary real \(q_1,q_2\), set

\[
 M(x)=c^{-1/2}DU\operatorname{diag}(q_1,q_2)V^TD^{-1}.
\]

The same direct expansion used in L24/L26 gives

\[
 d_x^2\det(4I-M(x)^TM(x))
 =E_0(1+r^2x^4)+E_1x^2, \tag{9}
\]

with

\[
 E_0=(4-q_1^2/c)(4-q_2^2/c)>0. \tag{10}
\]

Here (10) follows from \(|q_i|\le\tau_i\) and L23.  Consequently it is
enough to check the reciprocal orientation \(x^2=1/r\).  At this point

\[
 U\operatorname{diag}(1,r)V^T
 =\begin{pmatrix}\sqrt r&0\\1-r&\sqrt r\end{pmatrix},
\]

so the reciprocal orientation is exactly the palindromic weighted-shift
slice.  This explains why it is the sharp orientation rather than merely an
algebraic midpoint.  The odd block is

\[
 M_*={1\over\sqrt c(1+r)}
 \begin{pmatrix}
  \sqrt r(q_1+q_2)&(rq_1-q_2)/c\\
  c(q_1-rq_2)&\sqrt r(q_1+q_2)
 \end{pmatrix}. \tag{11}
\]

Define

\[
 \begin{aligned}
 \mathcal N={}&(1+r)^2(cq_1^2q_2^2+16c^3)\\
 &-4\{(c^2+r)^2q_1^2-2r(1-c^2)^2q_1q_2
 +(1+c^2r)^2q_2^2\}. \tag{12}
 \end{aligned}
\]

Then exactly

\[
 \det(4I-M_*^TM_*)={\mathcal N\over c^3(1+r)^2}. \tag{13}
\]

Also \(|\det M_*|=|q_1q_2|/c\le kp/c<4\), so nonnegativity of
\(\mathcal N\) puts both singular values below two.  Thus (OB) has reduced
to one explicit three-parameter scalar inequality: substitute (7),
\(\tau_1^2=k\), \(\tau_2^2=kp^2\), and \(r=H(p)\) into (12).

There is a substantially simpler equivalent test.  Write (11) as

\[
 M_*=\begin{pmatrix}\alpha&\beta\\\gamma&\alpha\end{pmatrix},
 \qquad \eta=\beta+\gamma,
 \qquad \delta=\det M_*={q_1q_2\over c}. \tag{13g}
\]

The two singular values of any real equal-diagonal matrix obey

\[
 \sigma_{\max}={|\eta|+\sqrt{\eta^2+4\delta}\over2}. \tag{13h}
\]

Here \(|\delta|<4\), by L23.  Squaring (13h), with the sign justified by
that determinant bound, shows the exact equivalence

\[
 \boxed{\|M_*\|\le2\quad\Longleftrightarrow\quad
 \delta+2|\eta|\le4,} \tag{13i}
\]

where

\[
 \eta={(r+c^2)q_1-(1+c^2r)q_2\over c\sqrt c(1+r)}. \tag{13j}
\]

Thus the quartic determinant numerator (12) can be replaced by two signed
bilinear inequalities.  Moreover, (13j) is monotone in \(r\): its derivative
has the sign of \((1-c^2)(q_1+q_2)\).  Consequently the envelope (14) again
reduces the full conformal displacement to its two cubic endpoints.

For certification it is better to replace the Blaschke zero \(a\) by its
value at the inner squared node,

\[
 t={kp^2-a\over1-akp^2}\in[-1,1],\qquad
 d={k(1-p^2)\over1-k^2p^2}. \tag{13e}
\]

The automorphism addition law then removes both rational denominators:

\[
 q_2=\sqrt k\,p t,\qquad
 q_1=\sqrt k\,{t+d\over1+dt}. \tag{13f}
\]

Thus the final compact variables can be taken as \((c,p,t)\), with only the
complete theta quantities \(k(c),s_0(c)\) left transcendental.  The sharp
ridge becomes \(p\to1/2\), \(t\sim-2c\).

In these variables put \(g=\sqrt{k/c}=2/\ell\),
\(X=q_1/\sqrt c=g(t+d)/(1+dt)\), and \(Y=q_2/\sqrt c=gpt\).  For
\(\varepsilon=\pm1\), define

\[
 J_\varepsilon=4-XY-2\varepsilon
 { (r+c^2)X-(1+c^2r)Y\over c(1+r)}. \tag{13k}
\]

Equation (13i) is exactly \(J_+,J_-\ge0\).  Multiplication by the positive
factor \(1+dt\) gives a quadratic

\[
 (1+dt)J_+=A t^2+B t+C, \tag{13l}
\]

with

\[
 \begin{aligned}
 A={}&{gp[2d(1+c^2r)-cg(1+r)]\over c(1+r)},\\
 B={}&{2c^2gpr-2c^2g-cdg^2p(1+r)+4cd(1+r)
          +2gp-2gr\over c(1+r)},\\
 C={}&4-{2dg(c^2+r)\over c(1+r)}. \tag{13m}
 \end{aligned}
\]

For \(J_-\), the quadratic coefficient is

\[
 -{gp[2d(1+c^2r)+cg(1+r)]\over c(1+r)}<0. \tag{13n}
\]

It is therefore concave in \(t\), and its minimum is at \(t=\pm1\), where
the block is \(\pm B\) and L26 applies.  Hence **the entire \(J_-\) half is
proved**.  The \(J_+\) quadratic also has nonnegative endpoints.  It can dip
inside \([-1,1]\) only when

\[
 A>0,\qquad |B|<2A. \tag{13o}
\]

On precisely that branch, the last upper-odd target is

\[
 \boxed{4AC-B^2\ge0} \tag{13p}
\]

at \(r=L(p)\) and \(r=U(p)\).  This removes both the orientation and the
Blaschke parameter from the remaining certificate: (13p) depends only on
\((c,p)\) and complete theta quantities.

The discriminant in (13p) has a much sharper square factorization.  Define

\[
 \begin{aligned}
 Q&=2g(1+c^2p)-d(kp+4c),\\
 S&=-2g(c^2+p)+d(kp+4c),
 \end{aligned} \tag{13q}
\]

where we used \(cg^2=k\).  Exact expansion of (13m) gives

\[
 \boxed{
 c^2(1+r)^2(4AC-B^2)
 =16rg^2p(1-c^2)^2(1-d^2)-(Qr-S)^2.} \tag{13r}
\]

This is not a numerical fit: comparison of the coefficients in \(r\) gives
the leading and constant terms \(-Q^2\) and \(-S^2\), while the difference
between half the middle coefficient and \(QS\) is exactly
\(8g^2p(1-c^2)^2(1-d^2)\).  The remaining positive factor also has the
elementary form

\[
 1-d^2=
 { (1-k^2)(1-k^2p^4)\over(1-k^2p^2)^2}. \tag{13s}
\]

Consequently the final target can equivalently be written as the sharp
two-variable distortion estimate

\[
 \boxed{|Qr-S|\le
 4g(1-c^2)\sqrt{rp(1-d^2)}}. \tag{13t}
\]

The right side of (13r) is a concave quadratic in \(r\).  Hence its minimum
on the conformal envelope interval is attained at \(L(p)\) or \(U(p)\),
which recovers the endpoint reduction without differentiating (13j).  The
only remaining content of the rank-one/full face is therefore (13t) at
those two cubic endpoints, restricted to the interior-vertex branch (13o).

There is no remaining inequality in the high-nome regime.  Recall
\(\ell=2/g\).  If

\[
 c\ell\ge1, \tag{13u}
\]

then the entire rectangle \(|q_1|\le\tau_1\),
\(|q_2|\le\tau_2\) maps into the norm-two ball.  By convexity it is enough
to check its four vertices.  The equal-sign vertices are \(\pm B\), already
closed by L26.  At an opposite-sign vertex the reciprocal block is, up to
global sign,

\[
 {\sqrt h\over1+r}
 \begin{pmatrix}
  \sqrt r(1-p)&(r+p)/c\\
  c(1+rp)&\sqrt r(1-p)
 \end{pmatrix}. \tag{13v}
\]

Its determinant numerator (12) factors as

\[
 -{16c\over\ell^4}F_-F_+, \tag{13w}
\]

where

\[
 F_\pm=\ell[c^2(1+pr)+p+r]
 \mathbin{\pm}c(1+r)(\ell^2+p). \tag{13x}
\]

Clearly \(F_+>0\).  Regard \(F_-\) as a concave quadratic in \(\ell\).
At \(\ell=1/c\),

\[
 F_-=-{(c^2-1)(p-1)\over c}\le0,
\]

and its derivative there is

\[
 c^2(1+pr)+p-r-2
 \le (p-1)(1+r)\le0.
\]

The derivative decreases thereafter, so (13u) implies \(F_-\le0\), and
(13w) is nonnegative.  As before \(|\det M|=hp<4\), and the reciprocal-
quadratic orientation identity propagates the vertex bound to all
orientations.  Therefore L29 is already proved whenever \(c\ell\ge1\).
The five-factor truncation of the product (3) in the L26 note gives this for
\(12599/20000\le c<1\).  Indeed, after clearing the positive denominator,
\(cP_5-1=(1-c)Q(c)/(\text{positive})\), where \(Q'(c)\) has nonnegative
coefficients and \(Q(12599/20000)>0\).  The exact polynomial audit is included
in `experiments/slice_odd_block_check.py`.  Thus the square target (13t)
need only be certified on

\[
 0<c<12599/20000. \tag{13y}
\]

For a finite certificate, (13r) has a second cancellation-free form.  Put

\[
 \lambda={d(kp+4c)\over2g},\qquad
 a=1+c^2p-\lambda,\qquad b=c^2+p-\lambda. \tag{13z}
\]

Thus \(Q=2ga\) and \(S=-2gb\).  Define

\[
 \begin{aligned}
 F_1&=-4c^2d+2cg(p+1)-dg^2p,\\
 F_2&=c^2dg^2p-2cg(p+1)+4d.
 \end{aligned} \tag{13aa}
\]

Direct factorization gives

\[
 p(1-c^2)^2(1-d^2)-ab={F_1F_2\over4g^2}. \tag{13ab}
\]

Consequently (13r) is equivalently

\[
 \boxed{
 c^2(1+r)^2(4AC-B^2)
 =4\{rF_1F_2-g^2D^2\},} \tag{13ac}
\]

where

\[
 D=ar-b=(r-p)-c^2(1-pr)+\lambda(1-r). \tag{13ad}
\]

At an envelope endpoint, evaluate \(p-r\) directly rather than subtracting:
it is \((1-s_0)p(1-p^2)\) for \(U\), and
\(p[1-s_0-a_3p^2]\) for \(L\).  Formula (13ac) exposes the two small
quantities separately and avoids the subtraction of two order-one squares
in (13r); it is the preferred interval-certificate form.

There is a further exact normalization which locates the only sharp
small-nome ridge.  Substitute

\[
 d={cg^2(1-p^2)\over1-c^2g^4p^2}
\]

in \(F_1\), and put

\[
 \begin{aligned}
 G_1(p)={}&g^3(1-2c^2g)p^2+(4c^2g-g^3)p+2-4c^2g,\\
 H={}&32c^4g^3-16c^4-24c^2g^2-g^4+8g.
 \end{aligned} \tag{13ae}
\]

Then direct cancellation gives

\[
 F_1={cg(1+p)\over1-c^2g^4p^2}G_1(p),
 \qquad \operatorname{disc}_pG_1=-g^2H. \tag{13af}
\]

Away from \(1-2c^2g=0\), this is the completed square

\[
 G_1(p)=g^3(1-2c^2g)(p-p_*)^2
       +{H\over4g(1-2c^2g)},
 \quad
 p_*={g^2-4c^2\over2g^2(1-2c^2g)}. \tag{13ag}
\]

The polynomial identity (13af), not the divided form (13ag), is used at
the exceptional parameter.  From the defining theta series,
\(p_*=1/2+3c^2/2+O(c^6)\).  Thus the previously observed sharp corner is
not an arbitrary guessed center: it is the exact vertex of the \(F_1\)
quadratic.  The remaining certificate should use \(p=p_*+c^4x\), rather
than axis-aligned boxes in \((c,p)\).  The identities (13af)--(13ag) are
exactly audited in `experiments/slice_odd_block_check.py`; no sign claim
about \(H\) is being assumed here.

The exact centered coordinate also makes the singular endpoint corner
finite.  Let \(\operatorname{Num}_U\) and \(\operatorname{Num}_L\) denote
the positive-denominator numerators of

\[
 rF_1F_2-g^2D^2
\]

at the upper and lower cubic envelopes.  Put

\[
 D_0=2g^2(1-2c^2g),\qquad p=p_*+c^4x. \tag{13ah}
\]

**Tiny-edge theorem.**  If \(0<c\le1/20\) and \(|x|\le8\), then

\[
 \operatorname{Num}_U>0,\qquad \operatorname{Num}_L>0. \tag{13ai}
\]

Here is the finite exact certificate.  The theta series and geometric
tails give

\[
 \begin{aligned}
 g={}&2-4c^2+10c^4-20c^6+36c^8-64c^{10}+110c^{12}+c^{14}P,
       &|P|\le200,\\
 s_0={}&1-4c^2+12c^4-32c^6+76c^8-168c^{10}+352c^{12}+c^{14}Q,
       &|Q|\le750.
 \end{aligned} \tag{13aj}
\]

Exact polynomial expansion, after clearing \(D_0^{12}\), gives

\[
 D_0^{12}\operatorname{Num}_U=c^{10}\{T(c^2,x)+R_U\},
 \quad K_U=9{,}895{,}604{,}649{,}984, \tag{13ak}
\]

where

\[
 \begin{aligned}
 T(t,x)={}&K_U(2x^2+4x+3)\\
 &-t(2{,}275{,}989{,}069{,}496{,}320x^2
      +4{,}690{,}516{,}604{,}092{,}416x
      +3{,}225{,}967{,}115{,}894{,}784)\\
 &+t^2(-13{,}194{,}139{,}533{,}312x^3
      +132{,}911{,}164{,}588{,}818{,}432x^2\\
 &\hspace{35mm}+280{,}995{,}589{,}640{,}945{,}664x
      +177{,}393{,}556{,}757{,}938{,}176).
 \end{aligned} \tag{13al}
\]

Exact tensor-product Bernstein coefficients on \(0\le t\le1/400\), with
the \(x\)-interval split into eighths, give

\[
 {T\over K_U}>0.8124\quad(-2\le x\le2),\qquad
 {T\over K_U}>2.3479\quad(-8\le x\le-2),\qquad
 {T\over K_U}>14.3177\quad(2\le x\le8). \tag{13am}
\]

Absolute coefficient domination over the terms of degree at least six in
\(c\) proves

\[
 |R_U|<0.1K_U\quad(|x|\le2),\qquad
 |R_U|<0.75K_U\quad(|x|\le8). \tag{13an}
\]

Thus (13am)--(13an) prove the upper assertion.  For the lower endpoint,
the shorter expansions

\[
 g=2-4c^2+10c^4-20c^6+c^8E,\quad |E|\le37,
 \qquad
 s_0=1-4c^2+12c^4+c^6J,\quad |J|\le33
\]

give

\[
 D_0^{12}\operatorname{Num}_L
 =c^8(K_L+R_L),\quad
 K_L=316{,}659{,}348{,}799{,}488,\quad |R_L|<0.9K_L. \tag{13ao}
\]

All coefficients, tail bounds, vanishing orders, and remainder sums are
reconstructed with rational arithmetic by
`experiments/slice_odd_tiny_edge_certificate.py`; the measured upper tail
bounds are \(0.085844K_U\) and \(0.741649K_U\), while the lower remainder
is \(0.843779K_L\).  The checker also
verifies that the discarded denominators are positive squares.  Thus this
is a proof, not a floating-point sweep.  This first certificate closes the
sharp core; the deeper centered completion below closes its complement.

The compact part away from the singular nome is also finite.

**Compact low-nome theorem.**  At either cubic-envelope endpoint, on the
branch (13o),

\[
 rF_1F_2-g^2D^2\ge0
 \qquad\left({1\over12}\le c\le{12599\over20000},\ 0\le p\le1\right).
 \tag{13ap}
\]

This is proved by the independently reproducible directed-interval checker
`experiments/slice_odd_compact_certificate.py`.  Its finite certificate uses
second-order Taylor models in two variables.  For \(1/12\le c\le1/2\), it
splits the \(p\)-interval at the exact \(p_*\) in (13ag) and uses

\[
 p=p_*v,\qquad p=p_*+(1-p_*)v,\qquad 0\le v\le1. \tag{13aq}
\]

The checker first interval-proves (0\le p_*\le1) on this range, so the two
charts cover every \(p\).  For \(1/2\le c\le12599/20000\), ordinary
\((c,p)\) boxes are already well conditioned.

No elliptic-function library or sampled modulus is used.  With (n=7), the
omitted theta tail and its first two derivatives are bounded, for example,
by

\[
 \begin{aligned}
 0\le T-T_6&\le {2c^{2n^2}\over1-c^{4n+2}},\\
 0\le T'-T_6'&\le {4n^2c^{2n^2-1}\over1-4c^{4n+2}},\\
 0\le T''-T_6''&\le
 {4n^2(2n^2-1)c^{2n^2-2}\over1-16c^{4n+2}},
 \end{aligned} \tag{13ar}
\]

with analogous bounds for \(R=\sum c^{2j(j+1)}\).  Every binary64
operation is expanded by one ulp outward.  On each box the resulting Taylor
intervals enclose \(A,2A+B,2A-B\), and the residual in (13ap).  A box is
accepted only if one branch form has upper endpoint at most zero, or the
residual has lower endpoint at least zero; otherwise it is bisected.
The chart value interval is intersected with its separately proved codomain
\([0,1]\); if a coarse dependency interval straddles a positive denominator,
the box is bisected rather than accepted.

The completed run has no unresolved boxes.  The ridge charts use 15,047
bisections (maximum depth 12; 25,600 branch-excluded and 56,647 proved box
evaluations).  The direct charts use 17,594 bisections (maximum depth 11;
36,160 branch-excluded and 7,434 proved evaluations).  Thus (13ap) is a
finite enclosure proof.  Combining it with L33 leaves only

\[
 0<c<{1\over12}. \tag{13as}
\]

Within (13as), L35 already removes the sharp tube
\(c\le1/20, |p-p_*|\le8c^4\).

The rest of (13as) can now be closed by two complementary exact
certificates.  First, deeper centered expansions remove the whole singular
band rather than only its core.

**Centered completion theorem.**  At both cubic-envelope endpoints,

\[
 \operatorname{Num}_U>0,\quad \operatorname{Num}_L>0
 \quad\left(0<c\le {1\over20},\quad
 \left|{p-p_*\over c^2}\right|\le50\right). \tag{13at}
\]

The same certificate also proves the bridge tube

\[
 \operatorname{Num}_U>0,\quad \operatorname{Num}_L>0
 \quad\left(0<c\le {1\over12},\ |p-p_*|\le4c^4\right). \tag{13au}
\]

For the upper endpoint, use the deeper exact enclosure

\[
 \begin{aligned}
 g={}&2-4c^2+10c^4-20c^6+36c^8-64c^{10}
       +110c^{12}-180c^{14}+c^{16}W,\\
 s_0={}&1-4c^2+12c^4-32c^6+76c^8-168c^{10}
       +352c^{12}-704c^{14}+c^{16}Z,
 \end{aligned}
 \qquad |W|\le300,\quad |Z|\le1500. \tag{13av}
\]

For the lower endpoint it is enough to use the shorter bounds

\[
 g=2-4c^2+10c^4-20c^6+36c^8-64c^{10}+c^{12}P,
 \quad |P|\le120,
\]
\[
 s_0=1-4c^2+12c^4-32c^6+76c^8+c^{10}N,
 \quad |N|\le200. \tag{13aw}
\]

All four remainder bounds hold on the larger interval \(c\le1/12\), by
the same rational theta tails used above.  With \(p=p_*+c^4x\), exact
expansion gives

\[
 D_0^{12}\operatorname{Num}_U=c^{10}\widehat U(c,x,W,Z),
 \qquad
 D_0^{12}\operatorname{Num}_L=c^8\widehat L(c,x,N,P). \tag{13ax}
\]

Tensor-product Bernstein bounds for the parts of \(\widehat U\) through
\(c^6\), and of \(\widehat L\) through \(c^4\), followed by absolute
coefficient bounds for the tails, give the following normalized data:

\[
\begin{array}{c|ccc|c}
 &\min U_{\le6}/K_U\ (|x|\le2)&
 \min U_{\le6}/K_U\ (x\le-2)&
 \min U_{\le6}/K_U\ (x\ge2)&|U_{\ge8}|/K_U\\ \hline
 c\le1/20,\ |x|\le16&0.81162&2.34015&14.23828&0.19648\\
 c\le1/12,\ |x|\le4&0.55453&1.48088&8.25522&0.42678/1.16046
\end{array} \tag{13ay}
\]

The two tail entries in the second row are for \(|x|\le2\) and
\(|x|\le4\), respectively.  For the lower endpoint,

\[
\begin{array}{c|cc}
 &\min L_{\le4}/K_L&|L_{\ge6}|/K_L\\ \hline
 c\le1/20,\ |x|\le16&0.76251&0.07344\\
 c\le1/12,\ |x|\le4&0.53031&0.22142
\end{array}. \tag{13az}
\]

Thus the two tubes in (13at)--(13au) are strict.  To finish (13at), put
\(t=c^2\) and \(y=(p-p_*)/c^2=tx\).  After the same exact substitutions,

\[
 D_0^{12}\operatorname{Num}_U=c^6Q_U(c,y),\qquad
 D_0^{12}\operatorname{Num}_L=c^6Q_L(c,y). \tag{13ba}
\]

On \(16t\le|y|\le50\), coefficient domination gives

\[
 Q_U=K_U(2y^2+4ty+3t^2)+E_U,\qquad
 |E_U|<1.474641K_Uy^2. \tag{13bb}
\]

But \(2y^2+4ty+3t^2=2(y+t)^2+t^2\ge(225/128)y^2\).
For the lower endpoint the leading coefficient is
\(A_L=712{,}483{,}534{,}798{,}848\); after discarding one positive
\(t\)-term,

\[
 Q_L=A_Ly^2+E_L+(\hbox{positive term}),\qquad
 |E_L|<0.968283A_Ly^2. \tag{13bc}
\]

This annulus meets the \(|x|\le16\) tube exactly and proves (13at).
`experiments/slice_odd_centered_completion_certificate.py` regenerates the
309,479-term and 565,425-term exact polynomials and checks every number in
(13av)--(13bc); saved expansion files are not inputs.

It remains to prove that (13at) contains every branch point not handled by
regular estimates.  The uncentered endpoint numerators have the common
leading polynomial

\[
 \operatorname{Num}_U=c^2L_0(p)+O(c^4),\qquad
 \operatorname{Num}_L=36c^2L_0(p)+O(c^4),
\]
\[
 L_0(p)=64p(p+1)^2(2p-1)^2(3-4p). \tag{13bd}
\]

The exact regular-region certificate uses four elementary pieces.

1. If \(p<c/4\), substitute \(p=cz\) in the positive-denominator
   numerator of the branch form \(2A-B\).  After removing \(c^2\), its
   leading terms are \(24z-12\) and \(144z-72\); the absolute remainders
   are only \(0.079427\) of their margins at \(z=1/4\).  Hence this strip
   does not meet the branch (13o).
2. On \(c/4\le p\le1/4\), \(L_0(p)\ge50p\).  Every higher term, bounded as
   a multiple of \(p\), is at most \(0.907878\) and \(0.907895\) of that
   leading lower bound at the upper and lower endpoints.
3. Put \(h=|2p-1|\).  On \(1/4\le p\le2/3\),
   \(L_0(p)\ge(25/3)h^2\).  If \(h\ge96c^2\), exact Bernstein bounds on
   the \(c^4/h\) term and coefficient domination on the \(c^6\) tail give
   total ratios \(0.675539\) and \(0.674419\).
4. On \(2/3\le p\le3/4\), put \(q=3/4-p\).  Here
   \(L_0(p)\ge(12800/243)q\), and the \(c^4\) constant is the helpful
   \(68943/256\) (36 times this below).  The \(q\)-correction ratios are
   \(0.709433,0.708693\), while the remaining \(c^2\)-ratios are
   \(0.095632,0.099831\).  If \(p\ge3/4\), the branch numerator instead
   starts with
   \(-42q-105c/16\) and \(-252q-315c/8\), where now \(q=p-3/4\); all
   potentially positive corrections use less than \(0.0890\) and \(0.695\)
   of these two margins.  Thus this second boundary strip is also outside
   (13o).

Finally, the exact identity

\[
 {p_*-1/2\over c^2}={g^3-2\over g^2(1-2c^2g)} \tag{13be}
\]

and the rational bounds \(19/10\le g\le2\) imply
\(1/2\le p_*\le1/2+2c^2\) for \(c\le1/20\).  Therefore the only middle
points not covered by item 3 have

\[
 |p-p_*|/c^2<48+2=50,
\]

and are in (13at).  This proves both endpoint residuals on the entire
branch for \(0<c\le1/20\).  The exact regenerating checker is
`experiments/slice_odd_small_edge_certificate.py`.

The last bridge \(1/20\le c\le1/12\) is finite.  Outside (13au), use the
two exact ridge-complement charts

\[
 p=(p_*-4c^4)v,\qquad
 p=p_*+4c^4+(1-p_*-4c^4)v,\qquad 0\le v\le1. \tag{13bf}
\]

The directed-interval checker first proves that their endpoints lie in
\([0,1]\), then applies the same second-order Taylor enclosure as L36.
It closes this bridge with 15,267 bisections (maximum depth 17), 952
branch-excluded and 17,715 residual-positive box evaluations, and no
unresolved boxes.  Together with (13au), this proves the endpoint residuals
through \(c=1/12\).

Combining the small-edge partition, bridge certificate, L36, and L33 covers
every \(0\le c<1\).  L30--L32 then propagate the endpoint square to every
Blaschke parameter and orientation.  Consequently

\[
 \boxed{\|B R_a(CB)\|\le2\quad(-1\le a\le1),} \tag{13bg}
\]

so **L29 is proved** and the generic rank-one/full KKT face in L28 is
closed.  By itself this is only a theorem for that coupled face of the
elliptic \(4\times4\) slice; the downstream L40--L59 chain now proves L20.
Neither result proves the general Crouzeix conjecture.

The **lower odd block is already proved**.  For fixed modal parameters its
matrix is linear in \((q_1,q_2)\), and the operator norm is convex.  Since
\(|q_i|\le\tau_i\), it is enough to check the four vertices of the rectangle
\([ -\tau_1,\tau_1]\times[-\tau_2,\tau_2]\).  Equal-sign vertices are L24.
For the opposite-sign vertex, the reciprocal orientation is, up to an
irrelevant global sign,

\[
 {\sqrt h\over1+r}
 \begin{pmatrix}
 c\sqrt r(1-p)&1+rp\\ c^2(r+p)&c\sqrt r(1-p)
 \end{pmatrix}. \tag{13a}
\]

Here \(c^2(r+p)\le1+rp\), and the triangle inequality reduces its norm to

\[
 {\sqrt h\over1+r}\,[c\sqrt r(1-p)+1+rp]. \tag{13b}
\]

By L23 it is enough to show the bracket is at most
\((1+c^2)(1+r)\).  With \(y=\sqrt r\), the difference is

\[
 c^2(1+y^2)+(1-p)y(y-c). \tag{13c}
\]

It is nonnegative when \(y\ge c\).  When \(y<c\), use
\(1-p\le1-r=1-y^2\) and \(y(c-y)\le c^2/4\); (13c) is then at least
\(c^2(1+y^2)-c^2/4>0\).  The reciprocal-quadratic reduction (9), now for
the lower block, propagates the midpoint bound to every orientation.  Hence

\[
 \boxed{\|C R_a(BC)\|\le2\quad(-1\le a\le1).} \tag{13d}
\]

## 4. A sharp algebraic envelope for the conformal displacement

Let

\[
 s_0=H'(0)={\pi\over2K(k^2)},\qquad
 a_3={s_0(1+k^2-s_0^2)\over6}.
\]

The Kanas--Sugawa coefficient theorem used in L26 says

\[
 H(p)=s_0p+a_3p^3+\sum_{j\ge2}a_{2j+1}p^{2j+1},
 \qquad a_{2j+1}\ge0.
\]

Since \(H(1)=1\), this yields the two-sided cubic envelope

\[
 \boxed{s_0p+a_3p^3\le H(p)
 \le s_0p+(1-s_0)p^3.} \tag{14}
\]

Before L30, numerical global optimization indicated the stronger raw
statement that (12) is nonnegative for **every** \(r\) in the interval
(14), not only for \(r=H(p)\).  That statement survived random and
differential-evolution searches, but is not separately claimed here.  The
proved route above first removes the Blaschke parameter and certifies the
smaller discriminant; it establishes (OB) without evaluating an incomplete
elliptic integral and makes this stronger raw formulation unnecessary.

For comparison with the original determinant route, as a polynomial in \(r\), write
\(\mathcal N=A_r r^2+B_r r+C_r\), where

\[
 \begin{aligned}
 A_r={}&16c^3+cq_1^2q_2^2-4q_1^2-4c^4q_2^2,\\
 B_r={}&32c^3+2cq_1^2q_2^2-8c^2(q_1^2+q_2^2)
       +8(1-c^2)^2q_1q_2. \tag{14a}
 \end{aligned}
\]

If \(A_r\le0\), the minimum on the envelope interval is at an endpoint.
If \(A_r>0\), all searches give \(B_r\ge0\), which makes the quadratic
increasing for \(r\ge0\) and again puts its minimum at an endpoint.  Thus a
particularly concrete proof package is now visible:

\[
 A_r>0\Longrightarrow B_r\ge0,
 \qquad \mathcal N(L(p))\ge0,
 \qquad \mathcal N(U(p))\ge0, \tag{14b}
\]

where \(L,U\) are the two cubics in (14).  This older sufficient package was
not needed: L30--L32 remove the Blaschke parameter first, and
(13at)--(13bg) certify the resulting smaller discriminant directly.  The
standalone implication in (14b) remains unaudited and should not be used as
an additional theorem.

The small-nome edge explains why the estimate is delicate.  With
\(a=Ac\),

\[
 k=4c+O(c^3),\qquad r=p-4c^2p(1-p^2)+O(c^4),
\]

and (11) tends to the one-entry matrix with norm

\[
 8p(1-p)\le2. \tag{15}
\]

Equality is at \(p=1/2\).  The next determinant term there is proportional
to \((A-3)^2c^2\), so the sharp corner follows the ridge
\(p\to1/2\), \(a\sim3c\).  Any interval proof must isolate this ridge rather
than expect a uniform positive margin.

## 5. Remaining proof debt

1. Treat L27's rank-one/rank-one coupled face.  The entire rank-one/full
   sector, including the singular \(c\to0\) ridge, is now closed by L29.
2. Audit whether the rank-one/rank-one KKT equations have an analogous
   colligation/Blaschke interpretation, rather than introducing four free
   defect-direction coordinates prematurely.

`experiments/slice_odd_block_check.py` audits (5)--(13), the conformal
envelope, and the symbolic reductions.  The three certificate scripts cited
above audit the complete parameter range.
