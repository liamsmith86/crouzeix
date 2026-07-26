# Leading scalar-channel leakage at the repeated monomial apex

> **Campaign scope.**  This is a local statement about analytic
> matrix-inner transfer curves through the repeated monomial
> \(z^L I_m\).  It is not a neighbourhood theorem and does not claim
> that L320's quartic reserve already absorbs every transverse term.

## 1. Result (L321, 2026-07-26)

Let

\[
B_s(z)=\sum_{n\ge1}B_n(s)z^n
\]

be a real-analytic curve of \(m\times m\) matrix-inner functions near
\(s=0\), analytic as an \(H^2\)-valued curve, and choose constant
input/output frames so that

\[
B_0(z)=z^L I_m.
\]

Suppose its first nonzero transfer jet has order \(q\ge1\), in the
\(\ell^2\) space of coefficient sequences:

\[
\begin{aligned}
B_L(s)&=I_m+s^qD_L+O(s^{q+1}),\\
B_n(s)&=s^qD_n+O(s^{q+1})\quad(n\ne L).
\end{aligned}                                      \tag{1}
\]

For L320's scalar-channel score

\[
\sigma(B_s)=
\sup_{\|u\|=\|v\|=1}\sum_n|u^*B_n(s)v|^2,
\]

one has the exact leading valuation

\[
\boxed{
1-\sigma(B_s)
=s^{2q}\Lambda(D)+O(s^{2q+1}),}                   \tag{2}
\]

where

\[
\boxed{
\Lambda(D)=
\min_{\|v\|=1}
\sum_{n\ne L}
\left\|(I-vv^*)D_nv\right\|^2.}                   \tag{3}
\]

In particular,

\[
\Lambda(D)=0
\quad\Longleftrightarrow\quad
\text{the off-monomial jets \(D_n,\ n\ne L\), have a common
eigenvector}.                                      \tag{4}
\]

The coefficient \(D_L\) is absent from (3): it only rotates the
common output line to first order.  Matrix-inner Parseval also forces
\(D_L+D_L^*=0\), consistently with that interpretation.

Combining (2) with L320 gives

\[
\bigl(1-\sqrt{\sigma(B_s)}\bigr)^2
=\frac14\Lambda(D)^2s^{4q}+O(s^{4q+1})            \tag{5}
\]

when \(\Lambda(D)>0\).  Thus the scalar-channel reserve is generally
quartic in the first channel-breaking amplitude.  Any merger which
spends it as a quadratic reserve is invalid.

If \(\Lambda(D)=0\), (2) says only that the leakage begins later.
The common eigenline must be carried to the next nonzero jet, or
split exactly by L205 if it persists.  This is the finite flag
required by the channel-stratum induction.

## 2. Fixed-input covariance

Fix a unit input vector \(v\), and put

\[
x_v(s)=B_L(s)v,\qquad
A_v(s)=\sum_{n\ne L}B_n(s)vv^*B_n(s)^*.
\]

The best constant output line for this fixed \(v\) has score

\[
\sigma_v(s)=\lambda_{\max}R_v(s),\qquad
R_v(s)=x_v(s)x_v(s)^*+A_v(s).                     \tag{6}
\]

Matrix-inner Parseval gives

\[
\operatorname{tr}R_v(s)
=\sum_n\|B_n(s)v\|^2=1.                           \tag{7}
\]

Let

\[
a_v(s)=\operatorname{tr}A_v(s),\qquad
w_v(s)=\frac{x_v(s)}{\|x_v(s)\|}.
\]

Uniformly on the unit sphere, \(a_v(s)=O(s^{2q})\), \(w_v(s)=
v+O(s^q)\), and (7) rewrites (6) as

\[
R_v(s)=\{1-a_v(s)\}w_v(s)w_v(s)^*+A_v(s).         \tag{8}
\]

The leading eigenvalue of (8) is separated uniformly from the other
eigenvalues.  In the decomposition
\(\mathbb Cw_v\oplus w_v^\perp\), its upper-left entry is

\[
1-a_v+w_v^*A_vw_v.
\]

The off-diagonal block has norm at most \(a_v\), while the lower
block has norm at most \(a_v\).  The Schur eigenvalue equation
therefore gives, uniformly in \(v\),

\[
\lambda_{\max}R_v
=1-a_v+w_v^*A_vw_v+O(a_v^2).                      \tag{9}
\]

Consequently

\[
1-\sigma_v(s)
=\operatorname{tr}\{(I-w_vw_v^*)A_v(s)\}
O(s^{4q}).                                        \tag{10}
\]

## 3. Insert the first transfer jet

Equation (1), in the Hilbert space of square-summable coefficient
sequences, implies uniformly in \(v\)

\[
A_v(s)=
s^{2q}\sum_{n\ne L}D_nvv^*D_n^*
+O_{\mathfrak S_1}(s^{2q+1}).                     \tag{11}
\]

Replacing \(w_v\) by \(v\) in (10) costs
\(O(s^{3q})=O(s^{2q+1})\).  Since
\(O(s^{4q})=O(s^{2q+1})\) as well,

\[
1-\sigma_v(s)
=s^{2q}\sum_{n\ne L}
\|(I-vv^*)D_nv\|^2
+O(s^{2q+1}),                                     \tag{12}
\]

with a remainder uniform on the unit sphere.

Finally,

\[
1-\sigma(B_s)=\min_{\|v\|=1}\{1-\sigma_v(s)\}.
\]

Uniformity in (12) permits the minimum to pass through the expansion,
proving (2)--(3).  The summands in (3) are nonnegative, so their
minimum is zero exactly when one unit vector is an eigenvector of
every off-monomial \(D_n\), proving (4).

The first-order coefficient of
\(\sum_nB_n(s)^*B_n(s)=I\) is
\(D_L+D_L^*=0\).  This confirms directly that no missing norm loss
from \(D_L\) belongs in (3).

## 4. Matrix-Schur tangent check

In L218's ordered matrix-Schur chart at \(z^LI_m\), a tangent
\(\Delta_j\) in parameter \(j\), \(1\le j<L\), contributes

\[
D_j=\Delta_j,\qquad
D_{2L-j}=-\Delta_j^*,                             \tag{13}
\]

and no other first-order coefficient.  For \(m=2,L=3\), take

\[
\Delta_1=a\sigma_x,\qquad
\Delta_2=a\sigma_z.
\]

The two Pauli matrices have no common eigenvector, and (3) becomes

\[
\Lambda
=2a^2\min_{\|v\|=1}
\{\operatorname{Var}_v(\sigma_x)
 \operatorname{Var}_v(\sigma_z)\}
=2a^2.                                            \tag{14}
\]

At \(a=1/4\), the predicted coefficient is \(1/8\).  Exact nested
Schur reconstruction gives

\[
\frac{1-\sigma(B_s)}{s^2}\longrightarrow0.125.
\]

Commuting diagonal tangents and a pure terminal-unitary rotation give
\(\Lambda=0\) and retain an exact scalar channel, as required.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_scalar_channel_valuation.py \
  --output \
  experiments/repeated_crabb_scalar_channel_valuation_s70224.jsonl
```

The checker reconstructs exact matrix-inner curves from their Schur
parameters, independently finite-differences the tangent
coefficients, globally optimizes the two-dimensional input sphere as
a diagnostic, and verifies:

1. the Pauli coefficient \(1/8\);
2. a generic nonnormal positive coefficient;
3. an exact commuting-channel path; and
4. an exact terminal-rotation channel path.

The numerical optimizer is not used in the proof of (2).  Equations
(6)--(12) prove the result.  The tracked dataset regenerates byte for
byte with SHA-256

```text
72b24fe92a8c2ec97148015cc540fba0f3a7bd7df3455bd8c886c4576c2fc035
```
