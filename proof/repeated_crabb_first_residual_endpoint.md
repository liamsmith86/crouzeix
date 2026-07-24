# The repeated first-residual endpoint and its metric flag

## 1. Result (L195, 2026-07-24)

Use L193's block disk chart and normalized rank-\(m\) Stein Gramian.
In physical coordinates put

\[
T=K^{1/2}AK^{-1/2},\qquad
P=K^{-1/2}MK^{-1/2}.                               \tag{1}
\]

Then:

1. **Exact zero-slack bridge.**  The metric \(P\) is exactly the
   zero-Stein-slack point of L194's analytic chart:

   \[
   \boxed{P=P(T,B,0),}                              \tag{2}
   \]

   where \(B\) is its level-zero/range row.

2. **Operator-valued first-residual square.**  Let
   \(\widehat H(s)\) be an analytic disk-Gram path through
   \(I/2\otimes I_m\).  Suppose the first nonzero terminal Hardy
   residual has order \(q\):

   \[
   P_I(A-S)A^rE_L
   =s^qF_r+O(s^{q+1}),\qquad 0\le r\le L-2.         \tag{3}
   \]

   Write \(F_{r,c}\in M_m\) for the level-\(c\) block of \(F_r\), and
   define

   \[
   {\cal G}(F)
   =\sum_{r=0}^{L-2}\sum_{c=1}^{L-1}F_{r,c}^*F_{r,c}
   \succeq0.                                       \tag{4}
   \]

   If \({\cal E}\) is L194's final-level upper endpoint, then

   \[
   \boxed{
   {\cal E}(s)
   =-4s^{2q}{\cal G}(F)+O(s^{2q+1}).}               \tag{5}
   \]

3. **First metric flag.**  The null space of the first face is

   \[
   \boxed{
   \ker{\cal G}(F)
   =\bigcap_{r,c}\ker F_{r,c}.}                     \tag{6}
   \]

Thus the scalar Frobenius deficit in L188 becomes a negative
semidefinite copy-space Gram matrix.  Strictly active copy directions
are closed immediately; only their common residual kernel advances to
the next flag.

L195 is not yet the repeated neighbourhood theorem.  One must prove
that a stable kernel in (6) can be promoted into a smaller-copy L193
equality block while controlling circular-normal and elliptic
couplings.

## 2. Why the canonical metric has zero Stein slack

The physical Stein defect is

\[
P-T^*PT=\widetilde Q\widetilde Q^*,\qquad
\widetilde Q=K^{-1/2}Q.                             \tag{7}
\]

Its level-zero block has rank \(m\) near the repeated Crabb point.
For any rank-\(m\) Gram column
\(\widetilde Q=(Q_0,Q_K)^T\) with invertible \(Q_0\),

\[
Q_KQ_K^*
-Q_KQ_0^*(Q_0Q_0^*)^{-1}Q_0Q_K^*=0.               \tag{8}
\]

Hence the Stein Schur complement in L194 is exactly zero.

Also \(P-I\succeq0\) locally by L193, its kernel has dimension \(m\),
and its range principal block is positive near the base.  Therefore
the level-zero lower Schur complement is zero.  The metric has both
defining properties of L194's zero-slack chart.  Local uniqueness
proves (2).

## 3. The upper endpoint as a Hardy quotient norm

Put

\[
G=4I-P\succeq0.                                    \tag{9}
\]

L194's endpoint is the negative Schur complement

\[
{\cal E}=-G/G_{JJ},\qquad J=\{0,\ldots,L-1\}.       \tag{10}
\]

Equivalently, for \(v\in\mathbb C^m\),

\[
-v^*{\cal E}v
=\min_x
\left\langle
G\binom{x}{v},\binom{x}{v}
\right\rangle.                                     \tag{11}
\]

L193's operator-valued Hardy defect realizes the quadratic form in
(11) as four times an orbit-complement \(H^2\) norm.  At the repeated
Crabb base every Hardy/model map is the corresponding scalar map
tensored with \(I_m\).  The scalar quotient calculation in L188 is an
orthogonal projection identity: the leading-coordinate orbit columns
span the model part, and the quotient coordinates are exactly the
finite residual blocks.  Tensoring that projection with \(I_m\)
preserves orthogonality and amplifies the scalar squared norm to the
matrix Gram.  Equivalently, apply the scalar quotient isometry to
each column \(F_{r,c}v\) and polarize in \(v\).

Thus the minimizing leading vector removes the model-space part and
leaves exactly the stacked coefficient column

\[
{\mathfrak F}v=(F_{r,c}v)_{r,c}.                   \tag{12}
\]

If the residual first occurs at order \(q\), orthogonal projection and
analyticity give

\[
\min_x
\left\langle
G(s)\binom{x}{v},\binom{x}{v}
\right\rangle
=4s^{2q}\|{\mathfrak F}v\|^2+O(s^{2q+1})\|v\|^2.   \tag{13}
\]

Polarization in \(v\) proves (5).  This is the operator-valued version
of L188's universal first-residual square; no trace or commutativity
is used.

Finally,

\[
v^*{\cal G}(F)v=\sum_{r,c}\|F_{r,c}v\|^2,
\]

which proves (6).

## 4. Consequence for the repeated campaign

Equations (2) and (5) align the two new charts:

- L193 describes the exact inverse-block-Toeplitz equality stratum;
- L194 supplies every convergent lower/Stein-tight metric nearby;
- L195 identifies the first upper loss as a copy-space Gram square.

There is no free metric choice to optimize at this first disk-residual
face: the canonical Hardy metric already occupies the zero-slack
branch and produces the correct negative endpoint.  The next
nontrivial issue is geometric, not analytic convergence:

1. split off the active range of \({\cal G}(F)\);
2. analyze the common kernel (6) as a smaller-copy equality problem;
3. retain L173/L188's strict circular-normal null lift on the active
   range; and
4. prove a matrix-Jensen alternative for elliptic support data on the
   kernel.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_first_residual_endpoint.py \
  --output \
  experiments/repeated_crabb_first_residual_endpoint_s70224.jsonl
```

The checker:

1. constructs the normalized physical Hardy metric at deterministic
   non-Toeplitz disk Gramians;
2. verifies its lower and Stein Schur complements vanish and its
   upper endpoint is nonpositive;
3. forms the exact linear residual blocks for a dense Hermitian
   direction; and
4. recovers (5) by a symmetric second difference, including
   noncommuting copy-space off-diagonal entries.

The standard dataset covers lengths two through five at
multiplicities two and three.  Equations (7)--(13) are the structural
proof.
