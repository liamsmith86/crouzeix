# Rigidity of repeated scalar-support blocks (2026-07-22)

## 1. Statement

L93 reduces a full common-top block to

\[
 Z^2=\alpha I,\qquad ZZ^*+Z^*Z=\beta I.              \tag{1}
\]

Every nonnormal irreducible solution is two-dimensional.  A possible
uniformity obstruction would be a direct sum of many equivalent
two-dimensional solutions: the multiplicity space might appear to create
new equality-preserving coupling directions not covered by the two-copy
tube L110.

It does not.  Let `Z_0=Z_ir tensor I_k`, where `Z_ir` is a nonnormal
two-dimensional solution of (1).  Suppose

\[
 Z(t)=Z_0+tX+O(t^2)                                  \tag{2}
\]

satisfies (1) to first order, allowing `alpha,beta` to vary.  Then

\[
 \boxed{X=[K,Z_0]+X_{\rm par}\otimes I_k,}            \tag{3}
\]

where `K^*=-K` and `X_par` is a tangent to the ordinary two-dimensional
solution family.  The parameter tangent has exactly three real dimensions.
Thus multiplicity contributes only unitary change of copy basis; it creates
no new local modulus.

## 2. Linearized relations

Differentiating (1) gives

\[
\begin{aligned}
 Z_0X+XZ_0&=\dot\alpha I,\\
 XZ_0^*+Z_0X^*+X^*Z_0+Z_0^*X&=\dot\beta I.           \tag{4}
\end{aligned}
\]

The first right side is complex and the second is real, accounting for the
three scalar parameter directions in (3).  The homogeneous system is (4)
with both right sides zero.

## 3. Square-zero stratum

After scaling and a unitary, the nilpotent irreducible and its multiplicity
are

\[
 Z_0=\begin{bmatrix}0&I_k\\0&0\end{bmatrix}.
\]

Write

\[
 X=\begin{bmatrix}A&B\\C&D\end{bmatrix}.
\]

The first equation in (4) says

\[
 C=cI_k,\qquad D=-A,                                  \tag{5}
\]

for one complex scalar `c`.  The second says

\[
 B+B^*=bI_k                                           \tag{6}
\]

for one real scalar `b`; its off-diagonal equation repeats `D=-A`.
Consequently

\[
 X=
 \begin{bmatrix}A&B-\frac b2I\\0&-A\end{bmatrix}
 +
 \begin{bmatrix}0&\frac b2I\\cI&0\end{bmatrix}.       \tag{7}
\]

The first matrix in (7) is `[K,Z_0]` for a skew-Hermitian block matrix `K`;
arbitrary `A` and arbitrary skew-Hermitian `B-bI/2` occur.  The second is
the same `2 x 2` parameter tangent on every multiplicity coordinate.  It
integrates explicitly:

\[
 \begin{bmatrix}0&1+tb/2\\tc&0\end{bmatrix}\otimes I_k
\]

still has scalar square and scalar `ZZ^*+Z^*Z`.  This proves (3) at the
nilpotent stratum, including its transition toward the invertible stratum.

## 4. Invertible nonnormal stratum

Normalize the phase and scale so that `Z_ir^2=I`.  Section 4 of L93 writes

\[
 Z_{\rm ir}=H+iK_0,\qquad
 H^2=a^2I,\quad K_0^2=b^2I,\quad HK_0+K_0H=0,         \tag{8}
\]

with `a,b>0`.  The normalized Hermitian matrices `H/a,K_0/b` are two
anticommuting symmetries and generate `M_2(C)`.  Hence `Z_0` is the
multiplicity-`k` representation of this matrix algebra.

Choose a single-block tangent `X_par` with the prescribed scalar
`dot alpha,dot beta` in (4), and subtract
`X_par tensor I_k`.  The remainder obeys the homogeneous equations.  In
write that remainder as `U+iV`, with `U,V` Hermitian, and choose the Pauli
normalization

\[
 H=a\sigma_z\otimes I_k,\qquad
 K_0=b\sigma_x\otimes I_k.
\]

The homogeneous square and anticommutator equations are equivalent to

\[
 \{H,U\}=0,\qquad \{K_0,V\}=0,\qquad
 \{H,V\}+\{K_0,U\}=0.                                \tag{9}
\]

Expand `U,V` in the Pauli basis with Hermitian `k x k` coefficients.  The
first equation leaves only the `sigma_x,sigma_y` coefficients of `U`; the
second leaves only the `sigma_y,sigma_z` coefficients of `V`; and the last
imposes the single Hermitian matrix relation
`a V_z+b U_x=0`.  Hence

\[
 \dim_{\mathbb R}\ker(4)_{\rm hom}=3k^2.              \tag{10}
\]

On the other hand, skew-Hermitian commutators `[K,Z_0]` lie in that kernel.
Their kernel consists exactly of the skew-Hermitian commutant
`I_2 tensor M_k`, of real dimension `k^2`.  Since the skew-Hermitian
matrices on `C^{2k}` have dimension `4k^2`, their commutator range has
dimension

\[
 4k^2-k^2=3k^2.                                      \tag{11}
\]

Equations (10)--(11) prove that every homogeneous tangent is an infinitesimal
unitary orbit direction.  Adding back the three scalar parameter directions
proves (3).  Equivalently, this is the elementary matrix-unit proof that
first-order deformations of a finite-dimensional `M_2` representation are
inner.

## 5. Consequence for the L93 flag

For either nonnormal stratum, the full relation-tangent dimension is

\[
 3k^2+3,                                              \tag{12}
\]

while the homogeneous/unitary-orbit dimension is `3k^2`.  Therefore a
weighted sequence approaching a multiplicity-`k` scalar-support block has
only two possibilities at its next scale:

1. it moves, after a small copy unitary, along the same repeated
   two-dimensional solution family, where L110 applies blockwise and domain
   monotonicity handles the direct sum;
2. it leaves the linearized scalar-support relations, in which case L93's
   Jensen/metric-flag endpoint has a genuine transverse negative
   coefficient.

There is no third, multiplicity-only equality direction.  This removes the
main algebraic obstruction to a uniform arbitrary-copy flag.  It does not
yet treat collisions with the normal stratum or merge the flag with L86's
losing-space gaps; those are the remaining steps.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_scalar_support_rigidity.py
```

The checker constructs the full real linearized systems for both a
square-zero irreducible and an invertible nonnormal irreducible.  For
multiplicities `k=1,...,4`, exact SymPy ranks prove

\[
 \dim T_{\rm full}=3k^2+3,\qquad
 \dim T_{\rm hom}=\dim T_{\rm orbit}=3k^2.
\]

The block calculations in Sections 3--4 are independent of `k`; the finite
rank checks are regression tests, not extrapolation.
