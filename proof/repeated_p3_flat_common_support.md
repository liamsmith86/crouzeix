# Common flat mode in copy-matrix coordinates (2026-07-22)

## 1. Setting

Use L88's complete flat copy matrix

\[
 E_Z=Z^*\otimes X_0+Z\otimes Y_0
\]

and add the common single-block flat mode

\[
 I\otimes W(w),\qquad W(w)=w(E_{10}+E_{21}).          \tag{1}
\]

The `Z` part has zero first top-support compression.  The common part has the
scalar compression

\[
 s_w(q)=\frac{\sqrt2}{4}(wq^{-2}+\bar wq^2),          \tag{2}
\]

so it preserves the tied copy eigenspace at first order.

## 2. Complete second effective support

Define

\[
 \begin{aligned}
 r_w(q)&=\frac1{16}
 \{2|w|^2-q^4\bar w^2-q^{-4}w^2\},\\
 \ell_w(q)&=\frac{\sqrt2}{16}
 \{q^{-1}w-q^3\bar w\}.                              \tag{3}
 \end{aligned}
\]

Exact reduced-resolvent multiplication gives

\[
 \boxed{
 Q_{Z,w}(q)=Q_Z(q)+r_w(q)I
 +\ell_w(q)Z+\overline{\ell_w(q)}Z^*,}                \tag{4}
\]

where

\[
 Q_Z(q)=\frac5{128}(ZZ^*+Z^*Z)
 -\frac3{128}\{q^2Z^2+q^{-2}(Z^*)^2\}.               \tag{5}
\]

Thus the entire L87 higher-order center—including the common flat mode—is an
explicit degree-two noncommutative Laurent polynomial in one copy matrix `Z`.
The scalar `r_w I` shifts every eigenvalue equally and is irrelevant to the
matrix-Jensen gap.  Only the Hermitian linear pencil

\[
 \ell_w(q)Z+\overline{\ell_w(q)}Z^*                  \tag{6}
\]

can change its common-top eigenspaces.

## 3. Consequence for the remaining classification

L88 already closes normal `Z` exactly.  In Schur coordinates, equations
(4)--(6) make the remaining task finite and algebraic:

- the diagonal of `Z` and `w` enter only through the linear pencil (6) and
  the diagonal contributions in (5);
- the strictly upper part enters through the anticommutator and square terms
  in (5);
- a common top vector is characterized by simultaneous eigenvector equations
  for the finitely many Laurent coefficients of (4).

This is the correct starting point for extending L89's zero-diagonal
three-copy classification to arbitrary Schur diagonal and common `w`; no
additional conformal-map fitting is needed at second order.

## 4. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_copy_matrix.py
```

The checker keeps a fully symbolic copy matrix, `w`, and all formal adjoints.
It verifies (2) and proves every coefficient of (4) exactly.
