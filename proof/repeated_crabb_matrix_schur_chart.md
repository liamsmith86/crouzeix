# Ordered matrix Schur coordinates for the repeated Crabb transfer

## 1. Result (L218, 2026-07-24)

Let \(U\in U(m)\), and consider a square rational inner function \(B\)
near

\[
B_0(z)=z^LU                                             \tag{1}
\]

whose determinant has winding number \(Lm\).  Assume \(B(0)=0\).
There are unique strict matrix Schur parameters

\[
\Gamma _0,\ldots,\Gamma _{L-1}\in M_m(\mathbb C)
\]

and a terminal unitary \(U_B\), all near \(0,\ldots,0,U\), such that

\[
\boxed{
F_j(z)=\Gamma_j+
D_{\Gamma_j^*}\,zF_{j+1}(z)
\bigl(I+\Gamma_j^*zF_{j+1}(z)\bigr)^{-1}D_{\Gamma_j},
}                                                       \tag{2}
\]

\[
F_0=B,\qquad F_L=U_B,\qquad
D_{\Gamma^*}=(I-\Gamma\Gamma^*)^{1/2},\quad
D_\Gamma=(I-\Gamma^*\Gamma)^{1/2}.                    \tag{3}
\]

The forced zero is

\[
\Gamma _0=B(0)=0.                                     \tag{4}
\]

Thus the local transfer chart has coordinates

\[
(\Gamma _1,\ldots,\Gamma _{L-1},U_B)
\]

and real dimension

\[
2(L-1)m^2+m^2=(2L-1)m^2.                              \tag{5}
\]

This exactly matches the dimension of L193's
inverse-block-Toeplitz equality manifold.

For a real variation \(\Gamma_j=\varepsilon\Delta\) at (1), with
all other nonterminal parameters zero and \(U_B=U\), the
linearization is

\[
\boxed{
\delta_j B(z)
=z^j\Delta-z^{\,2L-j}U\Delta^*U,
\qquad 1\le j\le L-1.
}                                                       \tag{6}
\]

Equation (6) is the noncommutative phase-palindromic reflection law.
It explains why a lower transfer coefficient and a later reflected
coefficient cannot be varied independently.

L201's transfer \(B_H\) lies in this chart at every sufficiently
nearby inverse-block-Toeplitz equality anchor.  Hence L218 supplies
an exact matrix Schur/Levinson generator for all of its transfer
coefficients.  It does **not** yet identify the nonlinear physical
ellipse metric or prove the all-grade prepared endpoint asserted in
the candidate one-image formula.

## 2. One ordered Schur step

For a strict contraction \(\Gamma\), put

\[
\mathcal S_\Gamma(Z)
=\Gamma+D_{\Gamma^*}Z(I+\Gamma^*Z)^{-1}D_\Gamma .
                                                               \tag{7}
\]

Direct multiplication gives the kernel identity

\[
\boxed{
\begin{aligned}
I-\mathcal S_\Gamma(Z(z))^*
       \mathcal S_\Gamma(Z(w))
={}&D_\Gamma(I+Z(z)^*\Gamma)^{-1}\\
&{}\cdot[I-Z(z)^*Z(w)]\\
&{}\cdot(I+\Gamma^*Z(w))^{-1}D_\Gamma .
\end{aligned}}                                             \tag{8}
\]

Therefore \(\mathcal S_\Gamma\) carries a unitary boundary value
\(Z\) to a unitary boundary value.

The inverse is also ordered:

\[
\boxed{
Z=D_{\Gamma^*}^{-1}(F-\Gamma)
  (I-\Gamma^*F)^{-1}D_\Gamma .
}                                                           \tag{9}
\]

Indeed, the defect intertwining

\[
\Gamma^*D_{\Gamma^*}=D_\Gamma\Gamma^*
\]

and (7) imply

\[
I-\Gamma^*F
=D_\Gamma(I+\Gamma^*Z)^{-1}D_\Gamma .
                                                               \tag{10}
\]

Substitution of (10) in (9) returns \(Z\).  In particular, if
\(\Gamma=F(0)\), the numerator in (9) vanishes at zero, so
\(F_{\rm next}=Z/z\) is analytic.

## 3. Termination and uniqueness

Apply (9) recursively with

\[
\Gamma_j=F_j(0),\qquad F_{j+1}(z)=Z_j(z)/z.          \tag{11}
\]

Equation (8) proves inductively that every \(F_j\) is rational inner.
The linear-fractional map in (7) is homotopic, through strict
parameters, to \(Z\mapsto Z\).  It therefore preserves determinant
winding.  Replacing \(Z\) by \(zF_{j+1}\) adds \(m\), and hence

\[
\operatorname{wind}\det F_{j+1}
=\operatorname{wind}\det F_j-m.                    \tag{12}
\]

Starting from \(Lm\), after \(L\) steps the winding is zero.  A
finite square rational inner function of zero determinant winding
has constant unimodular determinant.  Its values are contractions,
so the product of their singular values can have modulus one only
when every value is unitary; analyticity then makes the function
constant.  Thus \(F_L=U_B\).  Formula (11) makes every
parameter and the terminal unitary unique.  Conversely, (2)
reconstructs the unique rational inner function from these data.
Principal defect square roots are real analytic near zero, so this
is a real-analytic local chart.

For L201, the apex is (1).  Matrix innerness keeps the determinant
nonzero on the circle, so its winding number is locally constant and
equals \(Lm\).  Equation (4) follows independently from L201's
orthogonal defect frames.  This places the complete nearby transfer
family in the chart.

## 4. Tangent reflection

At the apex, the defect roots in (2) have zero first derivative.
At stage \(j\),

\[
zF_{j+1}(z)=z^{L-j}U.
\]

Differentiating (2) in the real direction
\(\Gamma_j=\varepsilon\Delta\) gives

\[
\Delta-z^{2(L-j)}U\Delta^*U.                       \tag{13}
\]

The \(j\) preceding zero-parameter steps multiply (13) by \(z^j\),
which proves (6).  The adjoint in (6) is forced by innerness; it is
not a commutative scalar analogy.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_matrix_schur_chart.py \
  --output \
  experiments/repeated_crabb_matrix_schur_chart_s70224.jsonl
```

The 16 deterministic records cover lengths two through five,
multiplicities two and three, both repeated apices and noncommuting
inverse-block-Toeplitz anchors.  They verify:

1. the forced zero and strictness of every recovered parameter;
2. termination at a unitary with zero tail;
3. coefficientwise inverse/forward round-trip;
4. boundary innerness and the ordered kernel identity; and
5. the tangent formula (6).

The exact identities (8)--(13), not the floating audit, prove L218.
The tracked dataset SHA-256 is
`c43f74ab4fdb14816232263a7b656688264d80ddab6317e746d6300b64262d51`.
