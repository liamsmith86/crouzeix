# Sign-separated coupled faces and the odd-Blaschke block

This note identifies the exact scalar theorem hidden in the generic
rank-one/full KKT face of `slice_coupled_defects.md`.  It does **not** yet
prove that theorem.  The reduction is useful because it turns a coupled dual
matrix problem into one explicit odd finite-Blaschke block and then removes
the orientation parameter.

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
 u={F-aI\over\sqrt{1-a^2}}x,qquad
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
sandwich eigenvalues are limiting cases.  Thus the theorem below would close
the generic rank-one/full face and its closure.

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
block can in fact be closed immediately (the proof is in §3); only the upper
inequality (OB) remains.  In conceptual terms, (OB) is the complete norm-two
theorem for the odd two-node Pick sector.

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
 q_2=\sqrt k\,p t,qquad
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
is a proof, not a floating-point sweep.  It closes the sharp tube only; the
rest of the low-nome rectangle remains below.

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
second-order Taylor models in two variables.  For (1/12\le c\le1/2), it
splits the (p)-interval at the exact (p_*) in (13ag) and uses

\[
 p=p_*v,\qquad p=p_*+(1-p_*)v,\qquad 0\le v\le1. \tag{13aq}
\]

The checker first interval-proves (0\le p_*\le1) on this range, so the two
charts cover every (p).  For (1/2\le c\le12599/20000), ordinary
((c,p)) boxes are already well conditioned.

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

with analogous bounds for (R=\sum c^{2j(j+1)}).  Every binary64
operation is expanded by one ulp outward.  On each box the resulting Taylor
intervals enclose (A,2A+B,2A-B), and the residual in (13ap).  A box is
accepted only if one branch form has upper endpoint at most zero, or the
residual has lower endpoint at least zero; otherwise it is bisected.

The completed run has no unresolved boxes.  The ridge charts use 15,203
bisections (maximum depth 12; 25,756 branch-excluded and 56,647 proved box
evaluations).  The direct charts use 17,594 bisections (maximum depth 11;
36,160 branch-excluded and 7,434 proved evaluations).  Thus (13ap) is a
finite enclosure proof.  Combining it with L33 leaves only

\[
 0<c<{1\over12}. \tag{13as}
\]

Within (13as), L35 already removes the sharp tube
\(c\le1/20, |p-p_*|\le8c^4\).

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

Numerical global optimization indicates the stronger algebraic statement
that (12) is nonnegative for **every** \(r\) in the interval (14), not only
for \(r=H(p)\).  This survived random and differential-evolution searches;
it is not yet a certificate.  Proving this envelope version would establish
(OB) without interval evaluation of an incomplete elliptic integral.

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

where \(L,U\) are the two cubics in (14).  The implication and the two
endpoint inequalities are still unproved; they are lower-dimensional
algebraic theta inequalities suitable for a finite certificate.

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

1. Prove the square distortion inequality (13t) on its branch (13o), at the
   two envelope endpoints, only for \(0<c<12599/20000\).  This is exactly the
   discriminant (13p); the high-nome complement is closed by (13u)--(13y).
   L36 and L33 reduce this to \(0<c<1/12\).  The sharp tube
   \(0<c\le1/20\), \(|(p-p_*)/c^4|\le8\), is closed by (13ah)--(13ao).
   Certify its complement, using multiscale coordinates near \(c=0\) rather
   than an axis-aligned box.
   The older equivalent route is \(\mathcal N\ge0\) under (14), but it
   retains an unnecessary Blaschke parameter.
2. Treat the rank-one/rank-one coupled face after this rank-one/full sector.

`experiments/slice_odd_block_check.py` audits (5)--(13), the conformal
envelope, and the current numerical status.
