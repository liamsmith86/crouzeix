# Matrix Faber reflection is not a matrix Blaschke lift

## 1. Result (L200, 2026-07-24)

The scalar elliptic merger L190 cannot be transferred to a repeated Crabb
block merely by replacing its scalar characteristic coefficients by copy
matrices.

The obstruction already occurs at repeated \(C_3\).  A positive
inverse-block-Toeplitz equality anchor has an exact matrix companion
polynomial \(G\), and the Dickson/Faber reflection formula remains valid
coefficientwise.  Nevertheless neither ordered analogue

\[
 G(G^\sharp)^{-1},\qquad (G^\sharp)^{-1}G             \tag{1}
\]

is unitary on the circle for a generic nonnormal copy coefficient.  Thus
neither is a matrix inner transfer, so L190's logarithmic inner-variation
and model-complement argument cannot be applied to it.

This disproves a **proof route**, not the desired elliptic similarity
bound.  In fact L114 already proves a full neighbourhood theorem for every
fixed repeated \(C_3\) multiplicity by a different metric-stratification
argument.  The issue here is whether L190's shorter scalar reflected proof
can become an all-length repeated-block mechanism.  A valid version must
construct the genuine matrix inner transfer from the rank-\(m\) unitary
colligation, equivalently from a matrix Schur/Toeplitz factorization; raw
coefficientwise Faber reflection is insufficient.

## 2. Exact repeated-\(C_3\) equality companion

Let \(Z\in M_m(\mathbb C)\) and suppose

\[
 B=\begin{bmatrix}2I&Z\\Z^*&2I\end{bmatrix}\succ0,
 \qquad H=B^{-1}.                                    \tag{2}
\]

Extend \(H\) by one zero block, let \(S=S_2\otimes I_m\), and use L193's
coefficient-gauge disk construction

\[
 K=H+S^*HS,\qquad A=2K^{-1}HS.                       \tag{3}
\]

Direct multiplication in \(KA=2HS\) gives

\[
 \boxed{
 A=\begin{bmatrix}
 0&2I&-Z/2\\
 0&0&I\\
 0&0&Z^*/2
 \end{bmatrix}.}                                     \tag{4}
\]

For example, the last-column identities follow from the off-diagonal
blocks of \(HB=I\):

\[
 H_{00}Z+2H_{01}=0,\qquad
 2H_{10}+H_{11}Z^*=0.                                \tag{5}
\]

Put \(C=Z/2\).  Since the first block column of \(A\) is zero, block
triangular elimination proves

\[
 \det(zI-A)=z^m\det G(z),\qquad
 G(z)=z^2I-C^*z.                                     \tag{6}
\]

The reciprocal-adjoint polynomial is

\[
 G^\sharp(z)=z^2G(1/\overline z)^*=I-Cz.             \tag{7}
\]

Equations (4)--(7) are exact and use no commutativity.

## 3. What remains true: coefficientwise Faber reflection

Let \(D_0=2,D_1=x\), and

\[
 D_2(x;c)=x^2-2c.
\]

On the Joukowski boundary \(x=\zeta+c/\zeta\),

\[
 D_j(\zeta+c/\zeta;c)=\zeta^j+c^j\zeta^{-j}.
\]

Consequently the coefficientwise Faber transform of (6) satisfies

\[
 \boxed{
 ({\cal F}_cG)(\zeta+c/\zeta)
 =G(\zeta)+c^2\zeta^{-2}I-C^*c\zeta^{-1}.}           \tag{8}
\]

Thus the reflected legs can still be read off as matrices.  The failure is
not in (8); it is in treating the positive matrix polynomial \(G\) as if
the scalar quotient \(g/g^\sharp\) remained inner.

## 4. Exact square-zero counterexample

Take

\[
 C=tN,\qquad
 N=\begin{bmatrix}0&1\\0&0\end{bmatrix},\qquad
 0<t<1.                                              \tag{9}
\]

Then (2) is positive because \(Z=2tN\) has norm below two.  At \(z=1\),

\[
 G(1)=I-tN^*,\qquad G^\sharp(1)=I-tN.
\]

The right ordered quotient is

\[
 G(1)G^\sharp(1)^{-1}
 =\begin{bmatrix}1&t\\-t&1-t^2\end{bmatrix}.         \tag{10}
\]

Its first column has squared norm \(1+t^2\), so it is not unitary.  The
left ordered quotient is

\[
 G^\sharp(1)^{-1}G(1)
 =\begin{bmatrix}1-t^2&t\\-t&1\end{bmatrix},         \tag{11}
\]

whose first column has squared norm

\[
 (1-t^2)^2+t^2=1-t^2+t^4\ne1.                       \tag{12}
\]

No choice of multiplication order repairs the scalar quotient.

This is precisely where scalar commutativity entered L190.  Although on
the circle

\[
 G^\sharp(\zeta)=\zeta^2G(\zeta)^*,
\]

the factors \(G,G^*,G^{-1},G^{-*}\) cannot be reordered to cancel.

## 5. Correct next object

L193 already supplies a rank-\(m\) defect and hence a genuine
matrix-valued unitary colligation.  Its transfer function is matrix inner,
but equations (10)--(12) show that it is not obtained by the naïve raw
ratio in (1).  The positive block Toeplitz matrix (2) also naturally
supplies a finite matrix Schur/Levinson recursion.

The next elliptic attack should therefore:

1. construct that colligation/Schur transfer with its multiplication order
   fixed;
2. express the Joukowski reflected legs in its orthogonal matrix-polynomial
   coordinates; and
3. derive the copy-space endpoint as an oriented Gram before attempting a
   statewise or flag argument.

Using (8) alone and asserting matrix innerness would recreate the same
ordering error exposed by L198.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_matrix_faber_obstruction.py \
  --output \
  experiments/repeated_crabb_matrix_faber_obstruction_s70224.jsonl
```

The checker verifies:

1. the exact companion formula (4);
2. the determinant identity (6);
3. the valid coefficientwise Faber identity (8);
4. nonnormality of the tested copy coefficient;
5. a definite unitary defect for both ordered quotients in (1); and
6. the exact repeated elliptic support formula for
   \(X+cX^*\), where \(X=K^{1/2}AK^{-1/2}\).

The standard data use multiplicities two, three, and four.  The tracked
dataset regenerates byte for byte with SHA-256

```text
8ae795a804886f0e6d3dee30b2ce5bd9758ec0b689bd4f3e0640781e8a02f9e4
```

The exact counterexample (9)--(12), not the floating audit, disproves the
naïve matrix-Blaschke lift.
