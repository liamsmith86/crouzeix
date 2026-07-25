# Every ordered boundary-metric flag has an exact transfer-Gram face

## 1. Result (L224, 2026-07-24)

Retain L223's balanced pure partial isometry

\[
\begin{gathered}
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\\
P=2I-E+2F,\qquad B_j=W^*(S^*)^jV,
\end{gathered}                                      \tag{1}
\]

and the positive upper gap

\[
U(q)=4P^{-1}-P_{\rm bl}(q),\qquad q=c^2.
\]

For \(k\ge1\), define the surviving left-copy flag

\[
{\cal K}_{k-1}
=\bigcap_{1\le j<k}\ker B_j^*,                     \tag{2}
\]

and let \(U_{k-1}:\mathbb C^{d_k}\to\mathbb C^m\)
be an isometry onto \({\cal K}_{k-1}\).

First Schur-eliminate the state complement of \(W\mathbb C^m\).
Then Schur-eliminate the earlier active copy directions
\({\cal K}_{k-1}^{\perp}\).  The remaining physical upper endpoint
has the exact first face

\[
\boxed{
\operatorname {Schur}_{WU_{k-1}}
\left(P_{\rm phys}(q)-4I\right)
=-4q^k
 (U_{k-1}^*B_k)(U_{k-1}^*B_k)^*
 +O(q^{k+1}).}                                     \tag{3}
\]

Here the Schur complement is ordered as just described; for every
fixed \(q>0\), associativity makes it the ordinary direct Schur
complement onto \(WU_{k-1}\).

Thus L223's fully delayed formula is not restricted to complete
copy-space delays.  It survives every rank-changing partial flag,
with the raw compressed transfer coefficient and no pseudoinverse or
extra whitening.  If \(U_{k-1}^*B_k=0\), the flag does not change and
the same assertion advances to the next nonzero grade.

L224 closes the endpoint-budget side of the ordered repeated flag.
It does **not** repair the indefinite elliptic Stein slack of
\(P_{\rm bl}\).  A contraction certificate must still be constructed
without spending all of the negative budget in (3).

## 2. Exact transfer-factor expansion

Put

\[
a_j(q)=\frac{q^j}{1+q^j},\qquad
E_j=(S^*)^jES^j,\qquad F_j=S^jF(S^*)^j.
\]

L219 gives

\[
U(q)=2E+\sum_{j\ge1}a_j(q)E_j
          +\sum_{j\ge1}(1-q^j)F_j.                \tag{4}
\]

Choose an isometry \(J\) onto \(W^\perp\), and set

\[
A(q)=J^*U(q)J,\qquad
C_j=J^*(S^*)^jV.
\]

The constant matrix \(A(0)\) is positive definite.  Moreover,
\(EW=0\), every \(F_jW=0\), and

\[
E_jW=(S^*)^jVB_j^*.
\]

Consequently

\[
\begin{aligned}
W^*U(q)W
 &=\sum_{j\ge1}a_j(q)B_jB_j^*,\\
J^*U(q)W
 &=\sum_{j\ge1}a_j(q)C_jB_j^*.
\end{aligned}                                      \tag{5}
\]

Let \(H(q)=\operatorname {Schur}_W U(q)\).  Equations
(5) give the exact expansion

\[
\boxed{
\begin{aligned}
H(q)
={}&\sum_{j\ge1}a_j(q)B_jB_j^*\\
&-\sum_{j,\ell\ge1}
 a_j(q)a_\ell(q)
 B_jC_j^*A(q)^{-1}C_\ell B_\ell^* .
\end{aligned}}                                     \tag{6}
\]

This formula is the mechanism behind the flag theorem.  The direct
term of grade \(j\) starts at \(q^jB_jB_j^*\), while every interaction
contains two transfer factors and starts at \(q^{j+\ell}\).

## 3. Compression to the surviving flag

Write \(K=U_{k-1}\) and let \(K_\perp\) be an isometry onto
\({\cal K}_{k-1}^{\perp}\).  By (2),

\[
K^*B_j=0\qquad(1\le j<k).                          \tag{7}
\]

Since \(a_j(q)=q^j+O(q^{2j})\), the \(K\)-compression of the direct
sum in (6) is

\[
q^k(K^*B_k)(K^*B_k)^*+O(q^{k+1}).                 \tag{8}
\]

Both transfer factors in the interaction term must have grade at
least \(k\) after compression on the two sides by \(K\).  Its
contribution is therefore \(O(q^{2k})=O(q^{k+1})\).
Thus

\[
K^*H(q)K
=q^k(K^*B_k)(K^*B_k)^*+O(q^{k+1}).                \tag{9}
\]

The cross block from \(K\) to \(K_\perp\) is \(O(q^k)\): every
direct term of smaller grade is killed by (7), and every interaction
with a surviving factor begins later.

The earlier flag steps split
\({\cal K}_{k-1}^{\perp}\) into positive blocks whose valuations are
at most \(k-1\).  Equivalently, the elementary analytic Schur
induction gives

\[
\left(K_\perp^*H(q)K_\perp\right)^{-1}
=O(q^{-(k-1)}).                                    \tag{10}
\]

For completeness, (10) follows inductively from (6): on the common
kernel through grade \(r-1\), the first direct term is
\(q^rB_rB_r^*\), while all two-factor interactions and the previous
Schur cross-square are \(O(q^{r+1})\).  Splitting its positive range
removes one nonzero copy block.  The finite copy dimension makes the
induction terminate.

The last copy-space Schur cross-square now has order

\[
O(q^k)\,O(q^{-(k-1)})\,O(q^k)=O(q^{k+1}).          \tag{11}
\]

It cannot change (9).  This proves the balanced version of (3).

## 4. Physical sign and factor

The physical upper gap is

\[
4I-P_{\rm phys}(q)=P^{1/2}U(q)P^{1/2}.
\]

Because \(P^{1/2}W=2W\), its complete left-copy Schur endpoint is
four times \(H(q)\).  The later copy-space Schur complements preserve
that common factor.  Finally, changing from the positive upper gap
to \(P_{\rm phys}-4I\) changes the sign.  This gives exactly (3).

## 5. Consequence for the Stein-repair problem

At every ordered flag, the available physical upper budget is the
negative Gram

\[
-4c^{2k}(U_{k-1}^*B_k)(U_{k-1}^*B_k)^*.
\]

This is stronger endpoint information than the pointwise
semidefinite alternative L222 requires, and it remains valid through
rank jumps.  What it does not control is the size of a correction
\(X\) satisfying

\[
(P_{\rm bl}+X)-\widehat T_c^*(P_{\rm bl}+X)\widehat T_c\succeq0.
\]

The next step is therefore quantitative: express a structured Stein
repair in the same transfer flag and show that its upper-endpoint cost
is strictly smaller than the factor \(4\) in (3).

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_boundary_metric_flag.py \
  --output \
  experiments/repeated_crabb_boundary_metric_flag_s70224.jsonl
```

The deterministic audit uses noncommuting matrix-Schur rank chains of
multiplicities two through five and rank-jumping heterogeneous shifts
through grade seven.  It constructs lossless partial-isometry
realizations, performs the ordered formal Schur reductions, and checks
the physical coefficient in (3).  Equations (4)--(11), not the
floating audit, prove L224.  The tracked dataset SHA-256 is
`641b1290380e3ed8efc46137c713334d1d8f5a2b7dd74d523c5cde8dde7f71ec`.
