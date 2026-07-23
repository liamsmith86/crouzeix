# Unique-winner reduction at repeated `C3` (2026-07-22)

## 1. First support decomposition

Let `A_0=I_m tensor C_3` and let `E` lie on L61's zero-Jensen face.  Choose
the common maximizing copy vector as `y=e_0`.  After subtracting the selected
copy's scalar boundary motion through the first Riemann correction, the pulled
first tangent has copy decomposition

\[
 G=\begin{bmatrix}G_0&Y\\X&G_-\end{bmatrix}.          \tag{1}
\]

Against the Crabb top support vector `r(q)`,

\[
 r(q)^*\operatorname{Re}(q^{-1}G_0)r(q)=0,            \tag{2}
\]

the selected/losing cross compression is zero by L74, and the losing
copy-space compression obeys

\[
 B_-(q):=V_-(q)^*\operatorname{Re}(q^{-1}G_-)V_-(q)
 \preceq0.                                             \tag{3}
\]

Call the face **unique winner** when

\[
 \bigcap_{|q|=1}\ker B_-(q)=\{0\}.                    \tag{4}
\]

Condition (4) allows the pointwise gap to vanish and allows eigenvalue ties at
isolated angles.  It excludes only a copy vector that ties the selected copy
for every angle; such a vector belongs in an enlarged winner sector instead.

## 2. Strict first-order slack on the losing sector

Since (3) is continuous and NSD, a vector in the kernel of its mean would lie
in every pointwise kernel.  Thus (4) implies

\[
 \overline B_-:=\frac1{2\pi}\int_0^{2\pi}B_-(e^{i\theta})\,d\theta
 \prec0.                                               \tag{5}
\]

The L61 mean identity identifies `overline B_-` with the adjacent-block matrix
in the tangent SDP.  Equivalently, restrict L61 equation (11) to the losing
copy space: its optimal first-order metric-bound coefficient is strictly
negative by (5).

Therefore one may choose a losing-sector first metric tangent for which all
three lower, upper, and Stein active compressions are strict while the global
bound coefficient remains zero.  To see strictness directly, start with a
losing tangent whose bound coefficient is negative, mix in an arbitrarily
small tangent-Slater point from L61 section 7, and use the remaining negative
margin.  This works even when (3) touches zero pointwise.

## 3. Additivity of the selected second endpoint

Let `a_0,a_1,a_2` collect the L74 star-quotient coefficients from the selected
copy to the losing copy space.  Use any L65 second-order metric certificate on
the selected diagonal block and the optimized L81 cross tangent

\[
 x=\frac{\sqrt2}{4},\qquad y=\frac{4\sqrt2}{3}.        \tag{6}
\]

The selected second endpoint splits **exactly** into diagonal and cross
pieces.  This is a block-parity fact:

- lower and upper Schur penalties are sums of squares over orthogonal copy
  sectors;
- the selected block of `G^*P_0G` is a sum over the intermediate copy index;
- the selected blocks of `G^*cal X A_0+A_0^*cal X G` have the same sector
  separation;
- every diagonal--cross term in the second pulled operator is off-diagonal in
  copy space and does not enter the selected endpoint.

Consequently no mixed diagonal--star term survives.  L81's exact algebra is
unchanged by the shape of (3), and the selected endpoint is

\[
 \boxed{
 e_{\rm total}\le e_{\rm selected}
 -\frac{25}{8}\|a_0\|^2-8\|a_1\|^2
 -\frac{50}{9}\|a_2\|^2,}                            \tag{7}
\]

where `e_selected<=0` is the single-`C3` coefficient from L63/L65.

The losing-sector constraints are already strict at first order by section 2,
so their second-order diagonal or internal blocks are absorbed.  The selected
constraints are controlled by (7).  Under the same uniform one-sided
conformal expansion used in L61--L62, the strict-lift argument proves the
corresponding upper second-order Dini bound.

## 4. Consequences

If the star quotient is nonzero, (7) is strictly negative.  Therefore every
unique-winner direction with nonzero selected/losing coupling descends at
second order, regardless of:

- the detailed nonconstant support-gap shape;
- isolated pointwise contacts in that gap;
- arbitrary blocks internal to the losing copy space.

If the star quotient is zero, the selected diagonal block is governed by the
single-block local theorem L73, while the losing sector has strict first-order
metric slack.  The remaining global difficulty at this stage was consequently
the **multiple-winner sector**, where (4) fails and all common kernel vectors
must be promoted into the winner multiplicity.  L83--L86 subsequently carry
out that promotion and close the complete fixed-direction second-order sign.

The scalar contact test

\[
 B_-(e^{i\theta})=-(1-\cos\theta)                     \tag{8}
\]

illustrates why uniform pointwise separation is unnecessary.  It is realized
by the copy difference `C_3-I`.  Nonlinear quotients for its three pure cross
generators converge to `-25/8,-8,-50/9`, exactly as (7) predicts.

## 5. Audit trail and scope

The universal cross calculation and its full second metric regenerate via

```bash
.venv/bin/python -u experiments/repeated_p3_radial_gap.py
```

L61 supplies the dimension-free tangent dual/Slater argument used in section
2; L63/L65 supply `e_selected<=0`.  The proof of additivity in section 3 is an
exact copy-block decomposition of those same three finite Schur complements,
not a numerical extrapolation from the radial model.

This theorem is a fixed-direction, unique-winner result.  L83--L86
subsequently close the multiple-winner fixed-direction sign.  A uniform
neighbourhood estimate as the mean gap in (5) tends to zero remains open.
