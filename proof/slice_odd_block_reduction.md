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

1. Prove the two-variable discriminant (13p) on its branch (13o), at the two
   envelope endpoints.  The older equivalent route is \(\mathcal N\ge0\)
   under (14), but it retains an unnecessary Blaschke parameter.
2. Treat the rank-one/rank-one coupled face after this rank-one/full sector.

`experiments/slice_odd_block_check.py` audits (5)--(13), the conformal
envelope, and the current numerical status.
