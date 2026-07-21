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
enough to check the reciprocal orientation \(x^2=1/r\).  There

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

1. Prove \(\mathcal N\ge0\) under the envelope (14), preferably by minimizing
   the quadratic in \(r\) at the two envelope endpoints and finding the
   small-nome square structure.
2. Treat the rank-one/rank-one coupled face after this rank-one/full sector.

`experiments/slice_odd_block_check.py` audits (5)--(13), the conformal
envelope, and the current numerical status.
