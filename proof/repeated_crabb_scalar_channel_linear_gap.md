# A sharp linear scalar-channel gap on repeated disk anchors

> **Campaign scope.**  This strengthens L320 on L193's repeated disk
> equality manifold.  It is a scalar Crouzeix estimate, not a
> complete local neighbourhood theorem.

## 1. Result (L322, 2026-07-26)

In L320's notation, let

\[
P=2I-VV^*+2WW^*,\qquad
C=P^{1/2}TP^{-1/2},
\]

and let \(\sigma(B)\) be the scalar-channel score of L201's
matrix-inner transfer.  Then every scalar Schur function \(f\)
satisfies the stronger estimate

\[
\boxed{
\|f(T)\|^2
\le
\frac{4+\sigma(B)+
\sqrt{\sigma(B)^2+8\sigma(B)}}{2}.}                \tag{1}
\]

In particular,

\[
\boxed{
4-\|f(T)\|^2
\ge\frac43\{1-\sigma(B)\}.}                        \tag{2}
\]

This replaces L320's quartic lower reserve
\((1-\sqrt{\sigma})^2\) by a reserve linear in \(1-\sigma\).
Combined with L321, an analytic matrix-inner curve through the
repeated monomial apex obeys

\[
\boxed{
4-\|f(T_s)\|^2
\ge\frac43\Lambda(D)s^{2q}+O(s^{2q+1}).}           \tag{3}
\]

Thus the available scalar reserve is quadratic, not quartic, in the
first channel-breaking amplitude.  L321 remains correct as an
expansion of L320's older bound; that bound was simply non-sharp.

The proof is an abstract sharp endpoint-angle inequality for a
contraction.  It uses no optimizer and no special property of \(f\)
beyond \(\|f(C)\|\le1\).

## 2. Sharp abstract endpoint-angle inequality

Let \(F\) be any contraction and put

\[
s=\|V^*FW\|\le1.
\]

Define

\[
\lambda(s)
=\frac{4+s^2+s\sqrt{s^2+8}}2.                     \tag{4}
\]

Then

\[
\boxed{
\|P^{-1/2}FP^{1/2}\|^2\le\lambda(s).}              \tag{5}
\]

The three eigenspaces of \(P\) are

\[
\operatorname{ran}V,\quad
{\cal R}=(\operatorname{ran}V\oplus
\operatorname{ran}W)^\perp,\quad
\operatorname{ran}W
\]

with eigenvalues \(1,2,4\), respectively.

Take a unit vector \(x\), and set

\[
y=P^{1/2}x,\qquad z=Fy,\qquad
a=\|WW^*y\|,\quad b=\|(I-WW^*)y\|.
\]

The input normalization gives the relaxed sharp constraint

\[
\frac{a^2}{4}+\frac{b^2}{2}\le1,                  \tag{6}
\]

because the \(V\)-part of \(y\) costs one rather than one half.
At the output,

\[
\begin{aligned}
\|P^{-1/2}z\|^2
&\le\frac12\|z\|^2+\frac12\|VV^*z\|^2\\
&\le\frac12(a^2+b^2)+\frac12\|V^*Fy\|^2.           \tag{7}
\end{aligned}
\]

To estimate the last term, fix a unit vector \(u\) in the \(V\)
defect coordinate and put

\[
t=\|WW^*F^*Vu\|\le s.
\]

Since \(F^*\) is a contraction,

\[
\|(I-WW^*)F^*Vu\|\le\sqrt{1-t^2}.
\]

Therefore

\[
|\langle Vu,Fy\rangle|
\le ta+\sqrt{1-t^2}\,b.
\]

Taking the output supremum gives

\[
\|V^*Fy\|
\le\max_{0\le t\le s}
\{ta+\sqrt{1-t^2}\,b\}.                            \tag{8}
\]

For a fixed \(t\), substitute

\[
a=2\xi,\qquad b=\sqrt2\,\eta.
\]

Equations (6)--(8) reduce the squared transformed norm to the largest
eigenvalue of

\[
M(t)=
\begin{bmatrix}
2+2t^2&\sqrt2\,t\sqrt{1-t^2}\\
\sqrt2\,t\sqrt{1-t^2}&2-t^2
\end{bmatrix}.                                    \tag{9}
\]

Its trace and determinant are

\[
\operatorname{tr}M(t)=4+t^2,\qquad
\det M(t)=4.
\]

Hence

\[
\lambda_{\max}M(t)
=\frac{4+t^2+t\sqrt{t^2+8}}2=\lambda(t).           \tag{10}
\]

This is increasing on \([0,1]\).  Taking the maximum over
\(0\le t\le s\) proves (5).

The bound is sharp for every \(s\).  In a three-dimensional
\(V\oplus{\cal R}\oplus W\) subspace, put

\[
F_s=
\begin{bmatrix}
0&-\sqrt{1-s^2}&s\\
0&s&\sqrt{1-s^2}\\
1&0&0
\end{bmatrix}.                                    \tag{11}
\]

This is unitary, \(V^*F_sW=s\), and the upper singular value of
\(P^{-1/2}F_sP^{1/2}\) is exactly \(\sqrt{\lambda(s)}\).

## 3. Insert the matrix-inner transfer

For \(F=f(C)\), von Neumann's inequality gives \(\|F\|\le1\).
L320's scalar Hardy argument proves

\[
\|V^*FW\|^2\le\sigma(B).                           \tag{12}
\]

Since \(\lambda\) is increasing, (5) and (12) give

\[
\|f(T)\|^2
\le\lambda(\sqrt{\sigma(B)}),
\]

which is exactly (1).

It remains to obtain the convenient linear form.  For
\(0\le\sigma\le1\),

\[
\frac{4+\sigma+\sqrt{\sigma^2+8\sigma}}2
\le4-\frac43(1-\sigma).                            \tag{13}
\]

After moving terms and squaring the two nonnegative sides, (13)
reduces to

\[
16(1-\sigma)^2\ge0.
\]

This proves (2).  Insert L321's
\(1-\sigma(B_s)=\Lambda(D)s^{2q}+O(s^{2q+1})\) to
obtain (3).

## 4. Consequence for the live merger

The strategic obstruction recorded after L321 is removed:
channel-breaking amplitude \(s^q\) supplies order \(s^{2q}\), not
order \(s^{4q}\), of scalar norm reserve.

The remaining repeated-block task is still nontrivial, but it is now
an ordinary associated-graded completion:

1. express the first circular-normal/disk response on a
   channel-adapted scalar compression;
2. use its vanishing on an exact L205 channel to factor it by the
   square root of \(1-\sigma\);
3. complete that response against L199's normal curvature while
   spending only a fixed fraction of (2); and
4. descend L321's common-eigenline flag when the leading factor
   vanishes.

No quartic endpoint calculation should be started unless this
quadratic factorization fails.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_scalar_channel_linear_gap.py \
  --output \
  experiments/repeated_crabb_scalar_channel_linear_gap_s70224.jsonl
```

The checker:

1. verifies equality in (5) on the unitary family (11);
2. tests random complex contractions with defect multiplicities one
   through three;
3. tests (1) on noncommuting L193 matrix-inner anchors and sampled
   scalar Schur polynomials; and
4. checks the global linear simplification (13).

The numerical work is an adversarial audit only.  Equations
(6)--(13) prove L322.  The tracked dataset regenerates byte for byte
with SHA-256

```text
52c16681fcddfcb21c08f4483093ea2f9e8a175e4dc8f22d6dfa2e90a46835b2
```
