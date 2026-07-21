# Globality of the level-4 even midpoint (Epoch 6, 2026-07-21)

This note closes the "global maximum versus stationary point" debt in L16.  The result is
domain-independent: it is a theorem about a two-node Schur/Pick problem and a complementary
two-idempotent frame.

## Theorem

Let (u_2<u_1) be real points of the unit disk, and let (Q_1,Q_2) be complementary
rank-one idempotents on a two-dimensional Hilbert space:

\[
 Q_1+Q_2=I,\qquad Q_iQ_j=\delta_{ij}Q_i.
\]

For a scalar Schur function (F), set

\[
 X_F=F(u_1)Q_1+F(u_2)Q_2.
\]

Let (v) be the real pseudo-hyperbolic midpoint of (u_1,u_2), and put

\[
 d=b_v(u_1)=-b_v(u_2)>0,\qquad R=Q_1-Q_2.
\]

Then

\[
 \sup_{\|F\|_\infty\leq1}\|X_F\|=\max\{1,\,d\|R\|\}. \tag{1}
\]

If (d\|R\|>1), every maximizer is, up to a unimodular scalar, the midpoint automorphism
(F=b_v).  In particular, the midpoint point is the **global** nonconstant maximizer, not
merely a stationary point.

## 1. Reduction to disk automorphisms

The two-node Pick body of pairs ((F(u_1),F(u_2))) is compact and contains the origin in its
interior.  A nonzero interior pair cannot maximize the norm: radial scaling reaches the boundary
and scales (X_F) by the same factor.  On the Pick boundary either:

1. one nodal value is unimodular, in which case both values are the same unimodular constant and
   (|X_F|=1); or
2. equality holds in Schwarz--Pick, in which case the pair is interpolated by a disk
   automorphism.

Thus a maximum larger than one is attained by an automorphism.  After composing with (b_v^{-1}),
every such automorphism has the form

\[
 F=e^{i\gamma}\varphi_a\circ b_v,\qquad
 \varphi_a(z)=\frac{z+a}{1+\overline a z},\qquad |a|<1.
\]

The irrelevant scalar (e^{i\gamma}) does not change the norm.  Functional calculus for the
involution (R^2=I) gives

\[
 X_a=\varphi_a(dR)=(dR+aI)(I+\overline a,dR)^{-1}. \tag{2}
\]

The midpoint is (a=0).

## 2. The automorphism orbit is maximized at the midpoint

Choose an orthonormal basis in which

\[
 R=\begin{bmatrix}1&2h\\0&-1\end{bmatrix},\qquad h\geq0,
\]

and write (T=dR), (lambda=\|T\|^2).  The singular-value polynomial at (a=0) is

\[
 \lambda^2-2d^2(1+2h^2)\lambda+d^4=0,
\]

so

\[
 h^2=\frac{(\lambda-d^2)^2}{4d^2\lambda}. \tag{3}
\]

Let (a=x+iy), (r=|a|^2).  Direct (2\times2) algebra using (3) yields

\[
 \det(\lambda I-X_a^*X_a)
 =\frac{(\lambda-1)\,J(a)}
 {|1+\overline a d|^2|1-\overline a d|^2}, \tag{4}
\]

where

\[
 J(a)=r\left[(\lambda+1)(d^4-1)r+2(1+d^2)(\lambda-d^2)\right]
      -4d^2(\lambda-1)x^2. \tag{5}
\]

Assume (lambda>1), the only regime relevant to a nonconstant extremal with norm larger than
one.  Since (x^2\leq r<1) and ((\lambda+1)(d^4-1)<0),

\[
 \frac{J(a)}r\geq
 (\lambda+1)(d^4-1)+2(1+d^2)(\lambda-d^2)-4d^2(\lambda-1)
 =(\lambda-1)(1-d^2)^2>0. \tag{6}
\]

Hence the determinant in (4) is strictly positive for every (a\ne0) in the disk.  As
(|a|\to1), (2) tends to the constant (aI), whose norm is one.  The punctured disk is
connected, and an eigenvalue of (X_a^*X_a) cannot cross (lambda) without making (4) zero.
It follows that

\[
 \|X_a\|^2<\lambda=\|dR\|^2\qquad(0<|a|<1). \tag{7}
\]

This proves (1) and the uniqueness statement when (d\|R\|>1).  If (d\|R\|\leq1), the
constant boundary functions give the value one, proving the other branch of (1).  (square)

## Consequence for L16

In the even sector of the level-4 reduction, every even disk Schur function has the form
(B(w)=F(w^2)).  The two collapsed nodes are (u_j=\tau_j^2), and the active block is exactly
(F(u_1)Q_1+F(u_2)Q_2).  Therefore, whenever the even-sector extremal norm is larger than one,
the pseudo-hyperbolic midpoint automorphism is the global even-sector maximizer.  The existing
(q_1=q_2=1/2) and EL4 formulas now apply without the former stationarity/globality caveat.

This theorem does **not** show that a global maximizer over all (not necessarily even) Schur
functions has definite parity.  That symmetry-breaking question, plus the odd and degree-one
phases, remains before the full 4-by-4 family is closed.

