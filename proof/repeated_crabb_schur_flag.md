# Finite Schur-orthogonal flag for repeated Crabb disk residuals

## 1. Result (L197, 2026-07-24)

Let \({\cal E}(s)\) be the final copy-space upper endpoint of L193's
canonical Hardy certificate along any real-analytic block-disk path.
Then

\[
S(s)=-{\cal E}(s)\succeq0                           \tag{1}
\]

is an analytic Hermitian \(m\times m\) matrix.  There is a finite
descending copy-space flag

\[
\mathbb C^m=W_0\supsetneq W_1\supsetneq\cdots
\supsetneq W_r                                      \tag{2}
\]

and analytic invertible congruences which put \(S\) into hierarchical
block form

\[
\boxed{
S(s)\sim
\operatorname{diag}\bigl(
s^{2\nu_0}A_0(s),\ldots,
s^{2\nu_{r-1}}A_{r-1}(s),0_{W_r}\bigr),}            \tag{3}
\]

where

\[
1\le\nu_0<\cdots<\nu_{r-1},\qquad A_j(0)\succ0.     \tag{4}
\]

The flag has at most \(m\) strict steps.  Its first member is exactly
L195's common residual kernel, and L196 identifies its first
inverse-Gram promotion.  Later members are obtained by Schur
orthogonalizing residual columns against all earlier active Hardy
ranges.

Thus the repeated block-disk endpoint has no infinite-order metric
recursion or hidden positive remainder.  L197 closes the **disk-only**
Schur flag.  It does not merge true circular normals or elliptic
support branches.

## 2. One analytic Schur step

The following elementary lemma supplies the induction.  Let
\(S(s)\succeq0\) be analytic and not identically zero.  If its first
nonzero Taylor coefficient is

\[
S(s)=s^{2\nu}S_{2\nu}+O(s^{2\nu+1}),\qquad
S_{2\nu}\succeq0,                                   \tag{5}
\]

split

\[
V=\operatorname{ran}S_{2\nu},\qquad
W=\ker S_{2\nu}.                                    \tag{6}
\]

The exponent is even and the coefficient is positive semidefinite:
apply scalar nonnegativity to \(v^*S(s)v\) for every \(v\).

In the fixed decomposition \(V\oplus W\), write

\[
S=\begin{bmatrix}A&C\\C^*&D\end{bmatrix}.           \tag{7}
\]

Then

\[
A=s^{2\nu}\widehat A,\qquad \widehat A(0)\succ0.
\]

The leading coefficient has no cross row into its kernel, so
\(C=s^{2\nu+1}\widehat C\).  Positivity forces the same order gain in
the kernel scalar forms.  Consequently

\[
X=A^{-1}C=s\,\widehat A^{-1}\widehat C              \tag{8}
\]

is analytic, and the exact congruence

\[
\begin{bmatrix}I&-X\\0&I\end{bmatrix}^*
S
\begin{bmatrix}I&-X\\0&I\end{bmatrix}
=
\begin{bmatrix}
A&0\\0&D-C^*A^{-1}C
\end{bmatrix}                                      \tag{9}
\]

splits off the active block.  The reduced Schur endpoint

\[
S_W=D-C^*A^{-1}C\succeq0                            \tag{10}
\]

is analytic and has strictly higher valuation unless it vanishes
identically.

## 3. Finite induction

Apply the step to \(S_W\).  Each nonzero step removes the positive
rank of its leading coefficient, so the domain dimension strictly
drops.  After at most \(m\) steps the process reaches either zero
dimension or an identically zero analytic block.  Multiplying the
finitely many analytic triangular congruences proves (3).

At the first step, L195 gives

\[
S_{2\nu_0}=4\sum_{r,c}F_{r,c}^*F_{r,c},             \tag{11}
\]

so \(W_1=\bigcap\ker F_{r,c}\).  Formula (10) has the standard least
squares meaning: it removes from later residual columns their
orthogonal projection onto the earlier active Hardy range.  This is
the Schur-orthogonal residual required after L196.

An identically zero terminal block in (3) is an exact analytic kernel
bundle for the canonical upper endpoint.  No assertion that this
bundle reduces every noncommuting equality coefficient is needed for
the disk certificate; L193 already gives \(P\preceq4I\) exactly.

## 4. Why this is enough for curve selection

Every possible local failure after adding other faces would yield a
real-analytic arc by finite-dimensional curve selection.  Along that
arc, (3) provides a finite hierarchy of strictly positive disk losses
on the successive active quotients.  Other response terms must either:

1. be absorbed at the first active quotient of this hierarchy; or
2. vanish and descend to its terminal analytic kernel.

Thus later mergers can work flag by flag without assuming a uniform
inverse gap across rank changes.  The rank-changing strata are the
next members of the same finite analytic hierarchy.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_schur_flag.py \
  --output experiments/repeated_crabb_schur_flag_s70224.jsonl
```

The checker first uses an exact polynomial Gram family with three
copy columns.  Its first Schur endpoint has valuation four and rank
one; the second has valuation six and positive leading coefficient.
It then constructs actual canonical Hardy endpoints with prescribed
first kernel dimensions and verifies positivity of the active and
Schur-reduced blocks.

The standard dataset covers lengths two through five, multiplicities
three and four, and every nontrivial first-kernel dimension.  Equations
(5)--(10) prove the general finite induction.
