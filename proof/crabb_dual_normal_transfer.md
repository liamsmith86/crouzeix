# Positivity transfer of the reflected/normal row (L166, 2026-07-23)

## 1. Result

Fix a reflected grade `k>=1`, put `d=k+1`, and use the physical
weighted path from L160--L163.  Let

\[
U_\epsilon(a,s)
\]

be L118's locally optimized rank-one Stein envelope, where the
equality amplitude is `epsilon*a` and the eligible circular normal is
inserted as `epsilon^d*sN`.  Let

\[
R_\epsilon(a,s)
=\|B_{\epsilon,a}(T_{\epsilon,a,s})\|^2              \tag{1}
\]

be the sharp prepared Blaschke lower certificate used in
L140--L149, frozen as a polynomial with respect to `s` and evaluated
at the strongly perturbed normalized operator.

Then the leading reflected/normal coefficients of `U` and `R` agree:

\[
\boxed{
[\,\epsilon^{2d}\,]\,
\partial_a\partial_s U_\epsilon(0,0)
=
[\,\epsilon^{2d}\,]\,
\partial_a\partial_s R_\epsilon(0,0).
}                                                     \tag{2}
\]

The same statement holds after real/imaginary polarization and for a
central reflected fold.  Thus the **one-reflection part** of L163 may
be proved on the sharp prepared Blaschke side.  No differentiated
Stein optimizer is needed for this transfer.  Terms nonlinear in the
unreflected equality amplitude are not represented by
`partial_a|_(a=0)`; they remain governed by L162's equality-ridge
transport.

Equation (2) does **not** say that the common coefficient is zero.
That is the remaining dual companion/inner calculation.

## 2. The elementary positivity lemma

Let `U,R` be real `C^2` functions near the origin with

\[
U\ge R,\qquad U(0)=R(0).
\]

Then

\[
H=D^2(U-R)(0)\succeq0.                               \tag{3}
\]

If `H(x,x)=0`, positivity of every `2 x 2` principal restriction
gives

\[
\boxed{H(x,y)=0\quad\hbox{for every }y.}             \tag{4}
\]

Equivalently,

\[
|H(x,y)|^2\le H(x,x)H(y,y).                          \tag{5}
\]

This is the same zero-diagonal PSD polarization used in L146, now
with the second variable allowed to be an arbitrary ambient strong
normal rather than another reflected grade.

## 3. Application at each nonzero scale

For every sufficiently small fixed nonzero `epsilon`,

\[
U_\epsilon(a,s)\ge t_*(T_{\epsilon,a,s})
\ge R_\epsilon(a,s).                                \tag{6}
\]

At `a=s=0`, L117 and the axis Blaschke certificate give

\[
U_\epsilon(0,0)=R_\epsilon(0,0)=4.                  \tag{7}
\]

Hence the `2 x 2` Hessian of

\[
\Delta_\epsilon=U_\epsilon-R_\epsilon
\]

in `(a,s)` is positive semidefinite.

L142--L145 prove equality of the primal and dual pure reflected
principal faces.  In the present weighting this is

\[
\partial_a^2\Delta_\epsilon(0,0)
=o(\epsilon^{2d}).                                  \tag{8}
\]

The strong insertion already contains `epsilon^d`; analyticity of the
two finite-dimensional constructions gives

\[
\partial_s^2\Delta_\epsilon(0,0)
=O(\epsilon^{2d}).                                  \tag{9}
\]

Applying (5) to (8)--(9),

\[
|\partial_a\partial_s\Delta_\epsilon(0,0)|
\le
\sqrt{
 \partial_a^2\Delta_\epsilon(0,0)
 \partial_s^2\Delta_\epsilon(0,0)}
=o(\epsilon^{2d}).                                  \tag{10}
\]

Analyticity in `epsilon` makes (10) equivalent to the coefficient
identity (2).

The argument is unchanged for two real coordinates representing a
complex amplitude.  A central fold is already a real coordinate in
the fixed phase gauge.  No constant hidden in (9) needs to be uniform
in dimension; the campaign is local in each fixed `p`.

## 4. Audit of hypotheses

1. **Lower certificate.**  `B_(epsilon,a)` is a finite Blaschke
   product, so (1) is bounded above by the functional-calculus square
   `t_*` for every nearby normalized operator.  Freezing the
   polynomial in `s` is legitimate and preserves analyticity.
2. **Upper certificate.**  L118's analytic optimized rank-one Stein
   metric is feasible, giving the first inequality in (6).
3. **Touching axis.**  L117 and L140 give (7), including simple active
   endpoints.
4. **Zero reflected gap diagonal.**  L142--L145 construct matching
   primal and dual one-grade coefficients `-64`; L146 supplies the
   associated-graded PSD formulation and covers mixed phases and
   folds.
5. **Strong scaling.**  The physical insertion is exactly
   `epsilon^d sN`.  The inverse Riemann, matrix evaluation, Stein
   solve, and simple endpoint maps are analytic, proving (9).

Thus the positivity transfer is an all-size theorem conditional only
on the already proved L117--L118 and L140--L146 results.

## 5. Remaining scalar gate

By (2), the compact one-reflection part of L163 is reduced to

\[
[\,\epsilon^{2d}\,]\,
\partial_a\partial_s
\|B_{\epsilon,a}(T_{\epsilon,a,s})\|^2=0,
\qquad k\ge2.                                       \tag{11}
\]

The grade-one coefficient is nonzero and must remain an exception.
The companion decomposition in L165 shows that (11) is precisely a
balance between endpoint basis motion and bottom-row characteristic
motion.  L149's logarithmic-inner identity handles the one-reflection
leg, while L162's Schwarz transfer handles arbitrary normalized
ambient motion on the zero-reflection ridge.  The unresolved step is
to combine those two identities in one differentiated companion
formula and verify that its fixed-degree hypothesis fails exactly in
grade one.  A complete L163 proof must then combine this reflected
identity with L162 for the zero-reflection/equality-amplitude
components of the same weighted coefficient.

L168--L169 now identify the differentiated object that should perform
this combination.  Its first coefficient is
`k(z^(-m)-z^m)`, and its exact nonlinear continuation is the inner
colligation transfer

\[
\theta(\zeta)
=d+\zeta\sqrt\alpha\,q^T(I-\zeta T)^{-1}Jq.
\]

The feedthrough `d` retains the moving characteristic constant and
avoids A115's fixed-zero failure.  The remaining obligation at this
stage was precise: show that (11) is the corresponding marked
real-mean coefficient of `dot(theta)/theta`.

L171 subsequently separates reflection count and proves that the
one-reflection row has no ordinary quadratic singular-Hessian term.
Its direct endpoint/cofactor term is reduced to L156's fully
differentiated Hardy functional applied to L168's
inner/optimized-defect split.  The required kernel value of that
functional is still open.
