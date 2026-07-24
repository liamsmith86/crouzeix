# Complete Jordan absorption of repeated Crabb circular normals

## 1. Result (L198, 2026-07-24)

Let \(F=(F_{ab})_{1\le a,b\le L-1}\), \(F_{ab}\in M_m\), be the
first block Hardy residual at a repeated Crabb disk point.  Its
block-reflection symmetry is

\[
JF=(JF)^*,\qquad
F_{ab}=F_{L-b,L-a}^*.                               \tag{1}
\]

For \(3\le k\le L-3\), put \(r=L-k\),

\[
D_{ab}={F_{ab}-F_{ba}\over2},\qquad
P_k=\sum_{\substack{a<b\\a+b=r}}(b-a)D_{ab},\qquad
C_k=\binom r3.                                      \tag{2}
\]

The correct complete amplification of L188 is **not** the one-sided
square \(P_k^*P_k\).  The two reflected circular characters give the
Jordan-symmetric square

\[
\boxed{P_k^*P_k+P_kP_k^*.}                          \tag{3}
\]

More precisely, let \(Y_k\in M_m\) be the complex copy coordinate in
L173's real/imaginary mode pair, and put

\[
\alpha_k={4(4k-1)\over L^2},\qquad G_k=\alpha_kP_k.
\]

Before scalarizing copy space, the residual/normal face is

\[
\begin{aligned}
{\cal N}_k(F,Y_k)
={}&-{b_{L,k}\over2}(Y_k^*Y_k+Y_kY_k^*)\\
&+{1\over4}\{
G_k^*Y_k+Y_k^*G_k+G_kY_k^*+Y_kG_k^*\}.              \tag{4}
\end{aligned}
\]

Here \(b_{L,k}\) is L173's exact curvature.  Completing the two
operator squares gives

\[
\boxed{
{\cal N}_k(F,Y_k)
\preceq {G_k^*G_k+G_kG_k^*\over8b_{L,k}}.}           \tag{5}
\]

L173's flux part is

\[
b_{L,k}^{\rm flux}={\alpha_k^2C_k\over64},
\qquad b_{L,k}>b_{L,k}^{\rm flux}.                  \tag{6}
\]

Consequently the operator-valued L195 Hardy loss absorbs every
circular-normal completion:

\[
\boxed{
\sum_k{\cal N}_k(F,Y_k)
\preceq4\sum_{a,b}F_{ab}^*F_{ab}.}                  \tag{7}
\]

The inequality is strict on every nonzero active normal component
after the L173 null lift is retained.  Thus noncommutativity creates
no new first-residual obstruction.  It does, however, make the
reflection pairing in (3) load-bearing.

L198 closes the **complete algebraic first-residual/circular-normal
face**.  It does not yet prove the full repeated circular-normal tube:
the later L197 Schur quotients must still be identified with reflected
Hardy residuals after their analytic recenterings.

## 2. Why the block jet is Jordan-symmetric

Repeat the L173/L188 coefficient paths before commuting their scalar
amplitudes.  The grade \(L-k\) path contains \(G_k^*Y_k+Y_k^*G_k\).
Reflection (1) sends it to the grade \(L+k\) path and reverses the
copy-product order, giving \(G_kY_k^*+Y_kG_k^*\).  The two paths have
the same scalar coefficient.  Dividing their sum by four gives the
second line of (4), which reduces to
\(\operatorname{Re}(\overline G_ky_k)\) when \(m=1\).

The same precommutative lift of L173's two reflected path chains gives

\[
-{b_{L,k}\over2}Y_k^*Y_k
-{b_{L,k}\over2}Y_kY_k^*.
\]

It reduces to \(-b_{L,k}|y_k|^2\) in one copy.  Thus (4) follows from
the actual ordered coefficient paths, not from guessing a tensor
product of the scalar answer.

At a repeated support eigenvalue, second-order numerical-range
normalization can additionally contribute the Jensen endpoint

\[
8L\left\{\overline Q-
\left(\mathop{\rm mean}_\theta\lambda_{\max}Q(\theta)\right)I\right\}
\preceq0,                                           \tag{8}
\]

where \(Q(\theta)\) is the effective support matrix.  This is L61's
tangent recurrence reapplied after the lower-order recentering.
Dropping (8) only weakens the upper estimate, so it is not needed for
(7).  It becomes useful when classifying equality and support-rank
changes.

## 3. Operator completion and weighted Cauchy

Split (4) into two halves.  The elementary identities

\[
\begin{aligned}
&-{b\over2}Y^*Y+{G^*Y+Y^*G\over4}
\preceq {G^*G\over8b},\\
&-{b\over2}YY^*+{GY^*+YG^*\over4}
\preceq {GG^*\over8b}
\end{aligned}                                      \tag{9}
\]

prove (5).  Using only the smaller flux curvature gives

\[
{G_k^*G_k+G_kG_k^*\over8b_{L,k}}
\preceq {8\over C_k}(P_k^*P_k+P_kP_k^*).            \tag{10}
\]

Matrix Cauchy--Schwarz, with
\(\sum_{a<b,\ a+b=r}(b-a)^2=C_k\), gives

\[
\begin{aligned}
P_k^*P_k&\preceq
C_k\sum_{\substack{a<b\\a+b=r}}D_{ab}^*D_{ab},\\
P_kP_k^*&\preceq
C_k\sum_{\substack{a<b\\a+b=r}}D_{ab}D_{ab}^*.
\end{aligned}                                      \tag{11}
\]

For one lower reflected pair write \(X=F_{ab}\), \(Y=F_{ba}\).
Its upper reflected pair is \(X^*,Y^*\), by (1).  Hence its
contribution to the full right Gram is

\[
X^*X+Y^*Y+XX^*+YY^*.                               \tag{12}
\]

Since

\[
\begin{aligned}
D^*D&\preceq\tfrac12(X^*X+Y^*Y),\\
DD^*&\preceq\tfrac12(XX^*+YY^*),
\end{aligned}                                      \tag{13}
\]

equations (10)--(13) show that the completed flux gain is at most four
times (12).  The active anti-diagonals are disjoint, and all unused
residual blocks add positive Gram terms.  This proves (7).

## 4. The one-sided amplification is false

If one incorrectly replaces (3) by \(2P_k^*P_k\), the proposed gain is

\[
{16\over C_k}P_k^*P_k.                              \tag{14}
\]

Take a nonnormal square-zero copy matrix \(X\), and on one lower
anti-diagonal put

\[
F_{ab}=(b-a)X,\qquad F_{ba}=-(b-a)X,
\]

then fill the upper anti-diagonal by (1).  Now

\[
P_k=C_kX,\qquad
4\sum F_{ab}^*F_{ab}=8C_k(X^*X+XX^*).               \tag{15}
\]

The symmetric flux gain in (10) equals (15) exactly, while the
one-sided residual is

\[
8C_k(XX^*-X^*X),                                    \tag{16}
\]

which is indefinite.  Thus a scalar-only tensorization would have
created a false lemma.  L173's null lift makes the correct symmetric
face strict.

## 5. Kernel compatibility with the first flag

Suppose a copy vector lies in the kernel of the complete leading
face.  The positive null-lift difference
\(b_{L,k}-b_{L,k}^{\rm flux}\) forces

\[
Y_kv=Y_k^*v=0
\]

for every active mode.  With the cross response gone, L195's Gram
forces \(F_{ab}v=0\) for every block.  Reflection (1) then also gives
\(F_{ab}^*v=0\).  Therefore the equality kernel reduces every leading
residual and every leading circular-normal copy matrix.

This is exactly the algebra needed to compress to the next copy-space
flag.  L196 supplies the first inverse-Gram promotion.  The remaining
proof debt is to show that L197's later Schur-orthogonal residuals,
after each analytic recentering, retain the reflected form (1) with
no additional mixed normal term.  Once that identification is made,
the argument above iterates in at most \(m\) steps.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_circular_jordan_absorption.py
```

The checker:

1. builds random noncommuting \(F\) with \(JF=(JF)^*\);
2. verifies the flux-only and actual L173 Loewner inequalities;
3. constructs the sharp square-zero example (15);
4. confirms that the symmetric flux face is exact there; and
5. confirms that the one-sided residual has a negative eigenvalue.

The standard dataset covers lengths \(6,\ldots,12\) and multiplicities
\(2,3,4\).  The proof is equations (1)--(13); the computation is an
indexing and noncommutativity audit.

Its SHA-256 hash is

```text
0843574dcbf656855f8a6aae0a2feb60bf8ddb4fc917970bbe1059d317113b49
```

An untracked broad audit covered every length \(6,\ldots,30\) and
multiplicity \(2,\ldots,5\) (104 records), with every check passing.
