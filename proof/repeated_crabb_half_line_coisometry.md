# The normalized elliptic half-line is a coisometry

## 1. Result (L244, 2026-07-24)

Retain L240's formal half-line operator, metric, and defect column

\[
A_\infty=\phi_c(\Xi_\infty),\qquad
D_R=\operatorname {diag}
\left(1,\frac1{1+q},\frac1{1+q^2},\ldots\right),
\]

\[
D_R-A_\infty^*D_RA_\infty=d_\infty d_\infty^*,
\qquad q=c^2.                                     \tag{1}
\]

Then the normalized defect is exactly isometric:

\[
\boxed{d_\infty^*D_R^{-1}d_\infty=I_m.}           \tag{2}
\]

Consequently

\[
\widetilde A_\infty
=D_R^{1/2}A_\infty D_R^{-1/2}
\]

is a partial isometry with initial defect

\[
e_\infty=D_R^{-1/2}d_\infty.
\]

In fact it is a coisometry:

\[
\boxed{
\widetilde A_\infty\widetilde A_\infty^*=I,
\qquad
A_\infty D_R^{-1}A_\infty^*=D_R^{-1}.}            \tag{3}
\]

Thus L240's zero Schur baseline is stronger than an accidental
rank-\(m\) Gram cancellation.  After the metric similarity it is the
standard defect geometry of a coisometry.  In particular, a
reflected state column and its quadratic Gram are naturally absorbed
into the updated defect column and disappear under right-defect
Schur elimination.  This does not by itself evaluate L228's
remaining linear cross term, but it explains why the first
nonconstant inverse-kernel sector in L243 must be assembled with the
Schur square rather than treated as an independent slack term.

## 2. Theta normalization

L240 gives

\[
\begin{aligned}
[d_\infty]_0&=\vartheta_3(q)^{-1}I_m,\\
[d_\infty]_{2j}
&=\frac{2(-c)^j}
 {\vartheta_3(q)(1+c^{4j})}I_m,\qquad j\ge1,
\end{aligned}
\]

and all odd coordinates vanish.  Since

\[
[D_R^{-1}]_{2j}=1+c^{4j},
\]

we obtain

\[
d_\infty^*D_R^{-1}d_\infty
=\frac{
1+4\sum_{j\ge1}q^j/(1+q^{2j})
}{\vartheta_3(q)^2}I_m.                            \tag{4}
\]

Jacobi's Lambert-series identity

\[
\boxed{
\vartheta_3(q)^2
=1+4\sum_{j\ge1}\frac{q^j}{1+q^{2j}}}             \tag{5}
\]

proves (2).  The series converge absolutely for \(|q|<1\), and (5)
also proves the identity coefficientwise formally.

Conjugating (1) by \(D_R^{-1/2}\) gives

\[
I-\widetilde A_\infty^*\widetilde A_\infty
=e_\infty e_\infty^*.                             \tag{6}
\]

Equation (2) says that the right side is an orthogonal projection.
Hence \(\widetilde A_\infty\) is a partial isometry with
\(\ker\widetilde A_\infty=\operatorname {ran}e_\infty\).

## 3. Surjectivity

It remains to rule out a left defect.  Write \(S_\infty=L^*\) and
\(E_\infty=J_0J_0^*\).  The ellipse pencil is

\[
\Xi_\infty=S_\infty+cS_\infty^*(I+E_\infty).      \tag{7}
\]

For \(y=(y_0,y_1,\ldots)\), set \(x_0=0\) and solve

\[
x_1=y_0,\qquad x_2=y_1,\qquad
x_{n+1}=y_n-cx_{n-1}\quad(n\ge2).                 \tag{8}
\]

The resulting geometric convolution is bounded on \(\ell^2\) for
\(|c|<1\), and (8) gives \(\Xi_\infty x=y\).  Thus
\(\Xi_\infty\) is onto.  Its kernel has dimension \(m\), from the
same recurrence with \(y=0\).

The centered odd Riemann map has the form

\[
\phi_c(z)=z\,h_c(z^2),
\]

where \(h_c\) has no zero on the closed spectral ellipse.  Therefore
\(h_c(\Xi_\infty^2)\) is invertible and

\[
A_\infty
=\Xi_\infty h_c(\Xi_\infty^2)
\]

is also onto with an \(m\)-dimensional kernel.  Similarity by
\(D_R^{1/2}\) preserves surjectivity.  A surjective partial isometry
is a coisometry, proving (3).

The formal coefficient version follows by stabilization from finite
sections exactly as in L240.  The bounded-operator argument above
additionally justifies (3) for every real \(0<c<1\).

## 4. Independent exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_half_line_coisometry.py \
  --output \
  experiments/repeated_crabb_half_line_coisometry_s70224.jsonl
```

The checker independently:

1. inverts \(D_R\) as an exact noncommutative operator series;
2. verifies
   \(A_\infty D_R^{-1}A_\infty^*-D_R^{-1}=0\)
   in the coisometric word algebra; and
3. verifies (4)--(5) coefficientwise.

All coefficients through degree eight vanish exactly.  The tracked
SHA-256 is
`1393ba2a734305af91775317759d3b7b3e8812d658e40f5fe82473bca66a5388`.
