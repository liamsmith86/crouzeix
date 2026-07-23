# Metric flag and the arbitrary-copy flat core (2026-07-22)

## 1. Canonical endpoint and its common-top kernel

Set the common flat mode to zero and use L88's perturbation

\[
 E_Z=Z^*\otimes X_0+Z\otimes Y_0.
\]

Its second effective support is

\[
 Q_Z(q)=\frac5{128}(ZZ^*+Z^*Z)
 -\frac3{128}\{q^2Z^2+q^{-2}(Z^*)^2\}.              \tag{1}
\]

Let

\[
 k_0=\mathop{\rm mean}_{|q|=1}\lambda_{\max}Q_Z(q).
\]

The canonical L83--L85 metric has upper endpoint

\[
 {\cal E}_0=16\{\mathop{\rm mean}Q_Z-k_0I\}\preceq0. \tag{2}
\]

Its kernel is exactly the common top space

\[
 K_0=\{x:Q_Z(q)x=\lambda_{\max}Q_Z(q)x
              \text{ for every }q\}.                \tag{3}
\]

This follows by averaging the pointwise positive matrices
`lambda_max(Q_Z(q)) I-Q_Z(q)`, as in L87.

## 2. A free metric variation

In physical-level coordinates, the free level-zero/level-one block of the
canonical first metric is

\[
 U_0=-\frac{3\sqrt2}{8}Z.                             \tag{4}
\]

For an orthogonal copy-space projection `P`, put `R=I-P` and vary it by

\[
 \dot U_P=PZR-\frac34RZP.                             \tag{5}
\]

The other first-metric blocks are forced exactly as in L80:

\[
 X_{11}=(E_Z^*MA_0+A_0^*ME_Z)_{11},\qquad
 X_{12}=2U+(E_Z^*MA_0+A_0^*ME_Z)_{12}.               \tag{6}
\]

Propagate the lower, Stein, and upper Schur complements with
`U=U_0+tau dot U_P`.  Exact multiplication gives

\[
\boxed{
 P{\cal E}'(0)P=
 -5\sqrt2\,PZRZ^*P
 -\frac{15\sqrt2}{4}\,PZ^*RZP.}                      \tag{7}
\]

In particular (7) is negative semidefinite, and its kernel in `ran P` is

\[
 \{x\in\mathop{\rm ran}P:Zx\in\mathop{\rm ran}P,\
                         Z^*x\in\mathop{\rm ran}P\}. \tag{8}
\]

Formula (7) is the dimension-free form of L80's missed tangent.  Internal
blocks on `ran R` cancel from its compression, so it is not restricted to a
pure star.

## 3. The descending metric flag

Starting from (3), define

\[
 K_{j+1}=\{x\in K_j:Zx,Z^*x\in K_j\}.                \tag{9}
\]

If `P_j` projects onto `K_j`, equations (7)--(8) say that the derivative
from `dot U_{P_j}` is strictly negative on `K_j` modulo `K_{j+1}`.
The finite descending chain stabilizes.  Its stable space is reducing for
`Z`, because `K_{j+1}=K_j` means both `ZK_j` and `Z^*K_j` lie in `K_j`.

Consequently, if `Z` is unitarily irreducible and

\[
 0\ne K_0\ne\mathbb C^m,                              \tag{10}
\]

the chain ends at zero.  A standard hierarchical perturbation now makes the
whole endpoint strictly negative.  Explicitly, use (5) first for `K_0`.
On its derivative kernel use a sufficiently small positive combination of
the variations for `K_1,K_2,...`; their nested negative compressions give a
direction whose derivative is negative definite on `K_1`.  Insert that
direction at order `tau^2` with a large fixed coefficient.  It dominates the
quadratic remainder on `K_1`, while the order-`tau` gap on
`K_0 minus K_1` and the order-zero gap on `K_0^perp` remain negative.

Thus an irreducible block has strict second-order descent unless its common
top space is the whole block.  The cases `K_0=0` and (10) are both closed at
second order.

## 4. Scalar-support blocks have size at most two

It remains to classify an irreducible block with `K_0=C^m`.  Every vector is
a top eigenvector of (1), so each Laurent coefficient is scalar:

\[
 Z^2=\alpha I,\qquad ZZ^*+Z^*Z=\beta I.              \tag{11}
\]

### Nilpotent case

If `alpha=0` and `beta>0`, then

\[
 (Z^*Z)^2
 =Z^*(\beta I-Z^*Z)Z
 =\beta Z^*Z.                                        \tag{12}
\]

Therefore `Z^*Z/beta` is a projection, `ZZ^*/beta` is its complementary
projection, and `Z/sqrt(beta)` is a partial isometry between the two
orthogonal spaces.  A unitary basis decomposes `Z` into identical `2 x 2`
square-zero blocks.  If `beta=0`, then `Z=0`.

### Invertible case

If `alpha!=0`, choose `s^2=alpha` and set `X=Z/s`.  Then

\[
 X^2=I,\qquad XX^*+X^*X=\gamma I.                    \tag{13}
\]

Write `X=H+iK` with `H,K` Hermitian.  Equations (13) give

\[
 H^2-K^2=I,\qquad HK+KH=0,\qquad
 2(H^2+K^2)=\gamma I,                                \tag{14}
\]

hence

\[
 H^2=\frac{\gamma+2}{4}I,\qquad
 K^2=\frac{\gamma-2}{4}I.                            \tag{15}
\]

If `gamma=2`, then `K=0` and `Z` is normal.  If `gamma>2`, normalized
versions of `H` and `K` are two anticommuting Hermitian unitaries.
The second swaps the `+1` and `-1` eigenspaces of the first, so paired
orthonormal bases decompose the representation into `2 x 2` blocks.

We have proved

\[
\boxed{\text{every irreducible solution of (11) has dimension at most two.}}
                                                                    \tag{16}
\]

The one-dimensional blocks are normal and belong to L88/L73.  Every
nonnormal two-dimensional block has trace zero and belongs to L92, which
gives strict cubic descent.

## 5. Arbitrary multiplicity, fixed directions

Decompose an arbitrary `Z` into its minimal unitary reducing blocks.  The
full repeated matrix decomposes by the same copy unitary.

On each irreducible block:

- no common top vector gives a strict canonical second endpoint;
- a proper common top space is strict by the metric flag;
- a full common top space is one- or two-dimensional by (16), hence is
  covered by L88 or L92.

Domain monotonicity transfers the complete-`2` estimates from the summands
to their direct sum.  Therefore every fixed `w=0` direction in the complete
L87 flat copy core is closed at exact, second, or cubic order, in arbitrary
copy multiplicity.  Together with L91, the same fixed-direction conclusion
holds when `w!=0`.

This is not yet a uniform neighbourhood theorem.  The remaining issue is a
weighted compactness/blow-up argument as the Jensen gap, the flag gaps,
`w`, and the cubic two-copy margin tend to zero together.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_metric_flag.py
```

The checker reconstructs the full lower/Stein/upper second-order metric with
a symbolic `3 x 3` copy matrix.  It proves (7) for both a rank-one and a
corank-one coordinate projection with independent formal adjoints.  The
arbitrary-dimensional identity is the same block multiplication and
polarization.
