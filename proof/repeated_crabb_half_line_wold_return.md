# The normalized half-line return has the ordinary Hardy metric

## 1. Result (L262, 2026-07-25)

Retain L244's normalized half-line operator and defect column

\[
U_c=D_R^{1/2}A_\infty D_R^{-1/2},\qquad
e_c=D_R^{-1/2}d_\infty .
\]

L244 proves, coefficientwise formally,

\[
U_cU_c^*=I,\qquad I-U_c^*U_c=e_ce_c^*,
\qquad e_c^*e_c=I_m.                              \tag{1}
\]

Define the Wold synthesis map on vector-valued Hardy space by

\[
\boxed{
{\cal W}_c^*(z^jv)=(U_c^*)^je_cv,\qquad j\geq0.}  \tag{2}
\]

Then \({\cal W}_c\) is a formal unitary:

\[
\boxed{
{\cal W}_c{\cal W}_c^*=I,\qquad
{\cal W}_c^*{\cal W}_c=I,\qquad
U_c{\cal W}_c^*={\cal W}_c^*L^*.}                 \tag{3}
\]

In particular, for every polynomial \(p\), and hence
coefficientwise for the formal analytic calculi used here,

\[
p(U_c){\cal W}_c^*
={\cal W}_c^*p(L^*),\qquad
{\cal W}_cp(U_c)^*=p(L^*)^*{\cal W}_c.            \tag{3a}
\]

Consequently every paired channel that enters and returns through
the normalized half-line has the **ordinary unweighted Hardy
metric**:

\[
\boxed{
({\cal W}_c^*Y)^*({\cal W}_c^*Y)
=Y^*Y
=\sum_{j\geq0}Y_j^*Y_j.}                          \tag{4}
\]

The order of copy matrices is unchanged because all half-line
coefficients are scalar multiples of \(I_m\).

L262 settles the metric of a genuine half-line Wold sandwich.  It
does not by itself prove that L258's first nonconstant physical
return is exactly such a sandwich.  The remaining A194 interface is
now:

> use L243's zero/one inverse-kernel filtration and L245/L251's
> paired state columns to place the first leakage row between
> \({\cal W}_c\) and \({\cal W}_c^*\), before coefficient
> extraction.

If that placement is proved, (4) gives the ordinary row norm and
L261 immediately collapses every future coefficient to
\(\|B_k\|_F^2\); L245's two remote factors supply the multiplier
four.

## 2. Orthogonality

Equation (1) implies

\[
U_ce_c=0.                                         \tag{5}
\]

For \(j\geq i\),

\[
\begin{aligned}
e_c^*U_c^i(U_c^*)^je_c
&=e_c^*(U_c^*)^{j-i}e_c\\
&=
\begin{cases}
I_m,&j=i,\\
0,&j>i,
\end{cases}
\end{aligned}
\]

where \(U_c^i(U_c^*)^i=I\) follows from coisometry and the last zero
uses (5).  Taking adjoints covers \(i>j\).  Thus the columns in (2)
are orthonormal and

\[
{\cal W}_c{\cal W}_c^*=I
\]

on \(H^2(\mathbb C^m)\).

## 3. Formal completeness

At \(c=0\),

\[
U_0=L^*,\qquad e_0=J_0.
\]

Therefore

\[
(U_0^*)^je_0=J_j
\]

and the constant coefficient of \({\cal W}_c^*\) is the identity on
Hardy space.  L240 defines every coefficient of \(U_c\) and \(e_c\)
as a finite-band operator, so

\[
{\cal W}_c^*=I+O(c)
\]

has a coefficientwise formal inverse.  Its orthonormal-column
identity says that this inverse is \({\cal W}_c\).  Hence

\[
{\cal W}_c^*{\cal W}_c=I
\]

as well.  This proves formal completeness without assuming a
boundary \(H^\infty\) calculus or excluding a unitary summand by
continuity.

Finally, (5) and the coisometry relation give

\[
U_c(U_c^*)^je_c=
\begin{cases}
0,&j=0,\\
(U_c^*)^{j-1}e_c,&j\geq1,
\end{cases}
\]

which is exactly the intertwining in (3).  Equation (4) is now
unitarity applied to an arbitrary finite Hardy column; formal
coefficient stabilization extends it to every coefficient required
by the campaign.  Iterating the intertwining and taking adjoints
proves (3a).
