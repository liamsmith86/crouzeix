# Jordan absorption guard for repeated Crabb circular normals

## 1. Conditional result and obstruction (A141, 2026-07-24)

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

The one-sided guess \(P_k^*P_k\) is not a valid automatic
amplification of L188.  A reflection-symmetric candidate instead uses
the Jordan square

\[
\boxed{P_k^*P_k+P_kP_k^*.}                          \tag{3}
\]

Suppose, conditionally, that a block endpoint calculation produces a
complex copy coordinate \(Y_k\in M_m\) in L173's real/imaginary mode
pair and the paired face below.  Put

\[
\alpha_k={4(4k-1)\over L^2},\qquad G_k=\alpha_kP_k.
\]

The candidate residual/normal face is

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

Equations (5)--(7) rigorously prove absorption **if** (4) is the
actual ordered block jet.  They do not prove that identification.
Scalar L173/L188 data do not by themselves determine every
noncommutative product order, and no independent full endpoint
derivation of (4) has been completed.

L199 avoids this unnecessary stronger claim.  A nonscalar normal copy
coefficient has an earlier strict support-Jensen gap.  On a
zero-Jensen winner it is scalar, so L188 amplifies state by state with
no product-order ambiguity.  Therefore this note is retained as an
adversarial guard and a useful conditional inequality, not as a
dependency of the repeated proof.

## 2. Why the candidate is natural but not yet a theorem

Formally lifting the two reflected L173/L188 character paths suggests
that the grade \(L-k\) path contains
\(G_k^*Y_k+Y_k^*G_k\), while grade \(L+k\) reverses the copy-product
order and gives \(G_kY_k^*+Y_kG_k^*\).  Their equal scalar
coefficients would give the second line of (4), reducing to
\(\operatorname{Re}(\overline G_ky_k)\) when \(m=1\).

The analogous formal lift of L173's two reflected path chains gives

\[
-{b_{L,k}\over2}Y_k^*Y_k
-{b_{L,k}\over2}Y_kY_k^*.
\]

It reduces to \(-b_{L,k}|y_k|^2\) in one copy.  However, establishing
these ordered terms from the full repeated Riemann/metric endpoint is
precisely the missing calculation.  The scalar restriction cannot
exclude additional noncommutative terms which vanish when \(m=1\).

At a repeated support eigenvalue, support normalization also has the
Jensen endpoint

\[
8L\left\{\overline Q-
\left(\mathop{\rm mean}_\theta\lambda_{\max}Q(\theta)\right)I\right\}
\preceq0,                                           \tag{8}
\]

where \(Q(\theta)\) is the effective support matrix.  L199 uses this
term first: it is strict for nonscalar normal coefficients, and its
zero set scalarizes the normal coefficient.  That route is both
weaker and fully justified.

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

## 5. Conditional kernel compatibility

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

This would be the algebra needed to compress to the next copy-space
flag if (4) were established.  L199 proves the required reducing
kernel directly on the scalarized zero-Jensen winner for the first raw
residual face.  A later Jordan lift may still be avoidable, but the
interaction of analytic normal response with L197's Schur-orthogonal
later residuals has not yet been derived.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_circular_jordan_absorption.py
```

The checker:

1. builds random noncommuting \(F\) with \(JF=(JF)^*\);
2. verifies the flux-only and L173-curvature conditional inequalities;
3. constructs the sharp square-zero example (15);
4. confirms that the symmetric flux face is exact there; and
5. confirms that the one-sided residual has a negative eigenvalue.

The standard dataset covers lengths \(6,\ldots,12\) and multiplicities
\(2,3,4\).  It audits equations (1)--(13), not the identification of
(4) with the actual repeated endpoint.

Its SHA-256 hash is

```text
0843574dcbf656855f8a6aae0a2feb60bf8ddb4fc917970bbe1059d317113b49
```

An untracked broad audit covered every length \(6,\ldots,30\) and
multiplicity \(2,\ldots,5\) (104 records), with every check passing.
