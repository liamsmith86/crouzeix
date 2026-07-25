# Every ordered lower metric flag has the right transfer-Gram face

## 1. Result (L226, 2026-07-24)

Retain L219's balanced pure partial isometry

\[
\begin{gathered}
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\\
P=2I-E+2F,\qquad B_j=W^*(S^*)^jV,
\end{gathered}                                      \tag{1}
\]

and put

\[
L(q)=P_{\rm bl}(q)-P^{-1},\qquad q=c^2.
\]

For \(k\ge1\), define the surviving right-copy flag

\[
\mathcal R_{k-1}
=\bigcap_{1\le j<k}\ker B_j,                       \tag{2}
\]

and let \(K_{k-1}:\mathbb C^{d_k}\to\mathbb C^m\)
be an isometry onto \(\mathcal R_{k-1}\).

First Schur-eliminate the state complement of \(V\mathbb C^m\).
Then Schur-eliminate every earlier active right-copy direction
\(\mathcal R_{k-1}^{\perp}\).  The remaining physical lower endpoint
has first face

\[
\boxed{
\operatorname {Schur}_{VK_{k-1}}
\left(P_{\rm phys}(q)-I\right)
=q^k(B_kK_{k-1})^*(B_kK_{k-1})
O(q^{k+1}).}                                      \tag{3}
\]

For every fixed \(q>0\), associativity makes the ordered operation in
(3) the ordinary direct Schur complement onto
\(VK_{k-1}\mathbb C^{d_k}\).

Thus the lower endpoint has the exact right-Gram counterpart of
L224's left-Gram upper budget.  The statement remains valid through
rank changes without a pseudoinverse.  It is an endpoint theorem,
not a Stein theorem: it does not repair the indefinite Stein slack
of \(P_{\rm bl}\).

## 2. Exact lower transfer expansion

Put

\[
E_j=(S^*)^jES^j,\qquad F_j=S^jF(S^*)^j.
\]

L219 gives

\[
\boxed{
L(q)
=\sum_{j\ge1}\left(\frac12-\frac{q^j}{1+q^j}\right)E_j
 \frac14F+\sum_{j\ge1}q^jF_j.}                    \tag{4}
\]

Choose an isometry \(J\) onto \(V^\perp\), and set

\[
A(q)=J^*L(q)J,\qquad C_j=J^*S^jW.
\]

The constant matrix \(A(0)\) is positive definite.  Since

\[
E_jV=0,\qquad FV=0,\qquad
F_jV=S^jWB_j,
\]

the endpoint and cross blocks are

\[
\begin{aligned}
V^*L(q)V
 &=\sum_{j\ge1}q^jB_j^*B_j,\\
J^*L(q)V
 &=\sum_{j\ge1}q^jC_jB_j.                          \tag{5}
\end{aligned}
\]

Let \(H(q)=\operatorname {Schur}_V L(q)\).  Equation (5) gives the
exact expansion

\[
\boxed{
\begin{aligned}
H(q)
={}&\sum_{j\ge1}q^jB_j^*B_j\\
&-\sum_{j,\ell\ge1}q^{j+\ell}
 B_j^*C_j^*A(q)^{-1}C_\ell B_\ell .
\end{aligned}}                                     \tag{6}
\]

The direct grade-\(j\) term begins at \(q^jB_j^*B_j\); every
interaction contains two transfer legs.

## 3. Compression to the surviving right flag

Write \(K=K_{k-1}\), and let \(K_\perp\) be an isometry onto
\(\mathcal R_{k-1}^{\perp}\).  By (2),

\[
B_jK=0\qquad(1\le j<k).                            \tag{7}
\]

Compressing the direct sum in (6) to \(K\) gives

\[
q^k(B_kK)^*(B_kK)+O(q^{k+1}).                     \tag{8}
\]

Both transfer factors in the interaction term must have grade at
least \(k\) after compression on both sides by \(K\).  That term is
therefore \(O(q^{2k})=O(q^{k+1})\).  Hence

\[
K^*H(q)K
=q^k(B_kK)^*(B_kK)+O(q^{k+1}).                    \tag{9}
\]

The cross block from \(K\) to \(K_\perp\) is \(O(q^k)\).  The earlier
flag steps split \(\mathcal R_{k-1}^{\perp}\) into positive blocks
whose valuations are at most \(k-1\).  The same elementary analytic
Schur induction as in L224 gives

\[
\left(K_\perp^*H(q)K_\perp\right)^{-1}
=O(q^{-(k-1)}).                                    \tag{10}
\]

Consequently the final copy-space Schur cross-square has order

\[
O(q^k)\,O(q^{-(k-1)})\,O(q^k)=O(q^{k+1}).          \tag{11}
\]

It cannot change (9), which proves the balanced statement.

## 4. Physical coordinates and relation to the repair

The physical lower gap is

\[
P_{\rm phys}(q)-I=P^{1/2}L(q)P^{1/2}.
\]

At the right defect,

\[
P^{1/2}V=V.
\]

Schur complements are covariant under the resulting block-diagonal
congruence, so no scalar factor is introduced.  This proves (3).

Together, L224 and L226 give the uncorrected endpoint faces

\[
\begin{aligned}
\text{lower: }&
\quad+c^{2k}(B_kK)^*(B_kK),\\
\text{upper: }&
\quad-4c^{2k}(U^*B_k)(U^*B_k)^*.
\end{aligned}
\]

L225 shows that on a complete delay, re-tightening the lower face
and repairing the first Stein-rank failure would produce the desired
upper trace \(-16\|B_k\|_F^2\), conditional only on the delayed slack
trace covariance.  L226 makes the lower cost explicit on every
partial flag and removes another possible rank-selection ambiguity.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_lower_metric_flag.py \
  --output \
  experiments/repeated_crabb_lower_metric_flag_s70224.jsonl
```

The deterministic audit uses noncommuting matrix-Schur rank chains
of multiplicities two through five and rank-jumping heterogeneous
shifts through grade seven.  It performs every ordered formal Schur
reduction and verifies (3).  Equations (4)--(11), not the floating
audit, prove L226.  The tracked dataset SHA-256 is
`db3e68e6967b65a84580b55785555b75fc98886c4ecba4e9a588d382302fcfca`.
