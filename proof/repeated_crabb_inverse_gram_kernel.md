# Promoting the repeated first-residual kernel in inverse-Gram coordinates

## 1. Result (L196, 2026-07-24)

Let \(B(s)=\widehat H(s)^{-1}\) be an analytic path in L193's block
disk chart through \(B(0)=2I\).  Let \({\cal T}\) denote block
Toeplitz projection by averaging each block diagonal, and put

\[
N(s)=B(s)-{\cal T}B(s).                             \tag{1}
\]

Thus \(N\) is the exact transverse inverse-Gram coordinate and
\(N=0\) is L193's equality manifold.

Suppose

\[
N(s)=s^qN_q+O(s^{q+1}),\qquad N_q\ne0.              \tag{2}
\]

Let \({\cal G}_q\) be L195's first copy-space residual Gram and let

\[
W=\ker{\cal G}_q\subset\mathbb C^m.                 \tag{3}
\]

Then

\[
\boxed{
N_q(I_L\otimes P_W)=0,\qquad
(I_L\otimes P_W)N_q=0.}                             \tag{4}
\]

Consequently, after recentering at the explicit equality anchor

\[
B_{\rm eq}(s)={\cal T}B(s),\qquad
H_{\rm eq}(s)=B_{\rm eq}(s)^{-1},                   \tag{5}
\]

the whole level-by-copy subspace
\(\mathbb C^L\otimes W\) and its cross blocks have no transverse
coefficient through order \(q\).  The residual valuation on this copy
kernel is strictly higher.

L196 is the exact first stable-kernel promotion needed by the repeated
metric flag.  It does not by itself control the next Schur-orthogonal
residual, circular normals, or elliptic support branches.

## 2. Adjacent difference is an exact transverse isomorphism

For an \(Lm\times Lm\) Hermitian block matrix define

\[
\Delta N=N_{1:L,1:L}-N_{0:L-1,0:L-1}.              \tag{6}
\]

Its kernel is exactly the Hermitian block-Toeplitz space.  On the
complement selected in (1), it is an isomorphism onto the Hermitian
\((L-1)m\)-square matrices.

The inverse is elementary.  Fix one block diagonal and write its
blocks as \(n_0,\ldots,n_{d-1}\), with zero mean.  If

\[
\delta_j=n_{j+1}-n_j,
\]

then

\[
n_0=-{1\over d}\sum_{j=0}^{d-2}(d-1-j)\delta_j,
\qquad
n_i=n_0+\sum_{j=0}^{i-1}\delta_j.                  \tag{7}
\]

Apply (7) independently on every block diagonal.  It preserves every
fixed copy-space compression and cross block.

## 3. Residual kernel equals transverse support kernel

At the repeated Crabb base, L193's exact reflection factorization
linearizes to

\[
C_q=-{1\over4}\Delta N_q\,{\cal Q}_0,               \tag{8}
\]

where \(C_q\) is the actual interior-output/terminal-Krylov residual
matrix and

\[
{\cal Q}_0=J_{L-1}\otimes I_m                      \tag{9}
\]

is the invertible block reversal in coefficient gauge.  In particular
\({\cal Q}_0\) preserves
\(\mathbb C^{L-1}\otimes W\).

Write \(C_{r,c}\) for the copy blocks of \(C_q\).  L195 gives

\[
{\cal G}_q=\sum_{r,c}C_{r,c}^*C_{r,c}.              \tag{10}
\]

For \(w\in\mathbb C^m\),

\[
w\in\ker{\cal G}_q
\Longleftrightarrow
C_{r,c}w=0\quad\hbox{for every }r,c
\Longleftrightarrow
C_q(\mathbb C^{L-1}\otimes w)=0.                   \tag{11}
\]

Equations (8)--(9) turn the last condition into

\[
\Delta N_q(\mathbb C^{L-1}\otimes W)=0.             \tag{12}
\]

Because \(\Delta N_q\) is Hermitian, the corresponding rows vanish
as well.  Applying the copy-compression-preserving inverse (7) proves
(4).

This argument explains why the flag is a genuine copy-space flag,
despite arbitrary noncommuting block-Toeplitz equality coordinates:
those tangential coordinates are absorbed into \(B_{\rm eq}\), while
the first transverse coefficient is supported only on \(W^\perp\).

## 4. What remains for a repeated neighbourhood theorem

At the first active order, L195 closes \(W^\perp\) and L196 raises the
transverse order on \(W\).  At the next order, residual columns on
\(W\) can have components in the previously active Hardy range.
Taking the upper Schur complement orthogonalizes those columns against
the earlier range.  The next lemma must identify that
Schur-orthogonal residual with a recentered inverse-Gram coefficient,
or provide a direct analytic Gram-flag induction.

After that disk-residual induction, two mergers remain:

1. tensor L188's weighted anti-diagonal circular-normal response and
   strict null lift with the copy flag; and
2. prove that elliptic copy data either reduce the stable flag to
   scalar single-block branches or have a strict matrix-Jensen gap.

No claim about those two steps is included in L196.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_inverse_gram_kernel.py \
  --output \
  experiments/repeated_crabb_inverse_gram_kernel_s70224.jsonl
```

The checker:

1. constructs dense zero-Toeplitz transverse coefficients supported
   on arbitrary non-coordinate copy complements;
2. reconstructs them from (7);
3. verifies the exact base residual factorization (8); and
4. checks that the residual Gram kernel and both row/column support
   kernels agree.

The standard dataset covers lengths two through six, copy
multiplicities two through four, and every nontrivial kernel
dimension.  Equations (6)--(12) are the all-size proof.
