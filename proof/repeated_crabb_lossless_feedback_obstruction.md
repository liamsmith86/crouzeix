# A single lossless feedback cannot transport the delayed quotient

## 1. Result (L268, 2026-07-25)

L267's abstract Redheffer defect identity does not directly implement
the one-delay recursion in A194.  Already on the scalar grade-two
monomial channel, the first full closed-return quotient has rank two,
whereas the deflated grade-one quotient has rank one.

More precisely, let

\[
S_2=
\begin{bmatrix}
0&0&0\\
1&0&0\\
0&1&0
\end{bmatrix},
\qquad
E_2=I-S_2^*S_2,\qquad F_2=I-S_2S_2^* .
\]

This is the monomial channel \(B(z)=z^2\).  Form the physical ellipse
map, delete the active metric coefficient as in L250, and eliminate
the final row as in L258.  On
\((I-F_2)\mathbb C^3=\operatorname {span}\{e_1,e_2\}\), exact series
arithmetic gives

\[
\boxed{
\mathfrak D_{\rm ret}^{S_2}(c)
=c^4
\begin{bmatrix}2&0\\0&2\end{bmatrix}
+O(c^5).}                                        \tag{1}
\]

Every coefficient below degree four is zero.

Removing the first clean layer leaves

\[
S_1=\begin{bmatrix}0&0\\1&0\end{bmatrix}.
\]

The independently balanced, edge-deleted, final-row-renewed tail has

\[
\boxed{
\mathfrak D_{\rm ret}^{S_1}(c)
=c^2
\begin{bmatrix}0&0\\0&4\end{bmatrix}
+O(c^3).}                                        \tag{2}
\]

Every coefficient below degree two is zero.  Thus the scalar
associated covariance is exact in this example:

\[
\operatorname {tr}[c^4]\mathfrak D_{\rm ret}^{S_2}
=4
=\operatorname {tr}[c^2]\mathfrak D_{\rm ret}^{S_1},
\]

but the operator faces have ranks two and one.

Consequently there is no analytic **single-port** lossless Redheffer
realization whose sole non-background defect is the deflated tail
quotient and which shifts that quotient by one clean layer.  In
particular, the proposed direct use of L267 cannot replace A194's
trace-only cancellation.

This does not rule out a larger multipoint network carrying an
additional independent defect channel.  Such a network would no
longer reduce the full face to the tail defect alone, however; its
extra channel would have to be evaluated and cancelled explicitly.

## 2. Rank obstruction

Suppose, contrary to the result, that an analytic unitary colligation
turns the deflated tail into the full physical map through degree
four, and that the tail Schur quotient is the only non-background
defect entering L267.  Factor out the moving rank-one background
defect and then take the fixed final-row quotient.  L263's
first-face calculation gives an effective analytic entrance
\(Q(c)\) such that

\[
K_{\rm full}(c)
=Q(c)^*K_{\rm tail}(c)Q(c)+O(c^5).               \tag{3}
\]

Equation (2) says that, on the one-dimensional tail quotient,

\[
K_{\rm tail}(c)=4c^2+O(c^3).                     \tag{4}
\]

The full quotient has no degree-two term.  Since the coefficient in
(4) is strictly positive, (3) forces \(Q(0)=0\).  Write

\[
Q(c)=cQ_1+O(c^2).
\]

The degree-four term in (3) is then

\[
4Q_1^*Q_1,                                       \tag{5}
\]

whose rank is at most one.  Equation (1) has degree-four coefficient
\(2I_2\), of rank two.  This is impossible.

The argument allows arbitrary analytic motion of the unitary
colligation and of the defect graph.  The obstruction is not that the
natural port was guessed incorrectly: any single effective entrance
from the one-dimensional tail quotient has the same rank bound.

## 3. Exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_lossless_feedback_obstruction.py \
  --output \
  experiments/repeated_crabb_lossless_feedback_obstruction_s70225.jsonl
```

The checker uses exact SymPy rational matrices.  It independently
constructs:

1. the full theta/ODE ellipse-map coefficients;
2. L219's boundary metric and the appropriate active deletion;
3. L258's output Gram and complete final-row renewal; and
4. the full and deflated closed-return quotients.

It verifies every earlier vanishing and the two exact matrices
(1)--(2), including their equal traces and unequal ranks.  No
floating-point tolerance or stored target is used.

The tracked dataset SHA-256 is
`fb53ef7433874e6be9887f68519bd467bf9c403e96864288b36d6aa2b4da4f65`.
The checker passes Ruff and `py_compile`.

## 4. Consequence for the live route

L267 remains a correct classical identity and may still organize
individual lossless subchannels.  What fails is A214's hoped-for
one-shot identification of the complete physical delay removal with
one transported tail defect.

The live proof must therefore remain trace-sensitive.  The two
explicit A213 stop conditions are again decisive:

1. place the complete L243/L251/L258 numerator in the two-sided
   leakage trace ideal; and
2. prove the post-contour boundary-depth/row filtration before using
   L266.

The rank created in (1) but absent from (2) is precisely the kind of
operator-level redistribution whose trace cancels correctly.  A
direct operator congruence is too rigid to capture it.
