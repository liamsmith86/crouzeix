# The delayed volume is a canonical one-frame model-space functional

## 1. Result (L257, 2026-07-25)

Retain a finite pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

with matrix-inner transfer

\[
B(z)=\sum_{j\geq1}B_jz^j.
\]

L236's left Hardy map is a unitary identification of the state space
with the model space

\[
\mathcal K=\mathcal K_{zB}
=H^2(\mathbb C^m)\ominus zB H^2(\mathbb C^m).
\]

In this single canonical frame, every object in L250's scalar volume
has an exact formula involving only \(B\), fixed Hardy shifts, and
fixed diagonal weights.  No state-space realization remains.

Let \(L\) be the forward shift, let \(P_0\) project onto the constant
row, and let

\[
\mathscr H=\mathcal O_R\mathcal O_L^*,\qquad
\mathscr H_{n,j}=B_{n+j}^*                         \tag{1}
\]

be the transfer Hankel operator, restricted to \(\mathcal K\).
Then

\[
\boxed{
\begin{aligned}
S_B^*&=L^*|_{\mathcal K},\\
S_B&=L(I-E_B)|_{\mathcal K}=P_{\mathcal K}L|_{\mathcal K},\\
F_B&=P_0|_{\mathcal K},\\
E_B&=\mathscr H^*P_0\mathscr H|_{\mathcal K}.
\end{aligned}}                                    \tag{2}
\]

For \(q=c^2\), put

\[
D_R=\operatorname{diag}
\left(1,\frac1{1+q},\frac1{1+q^2},\ldots\right),
\qquad
D_L=\operatorname{diag}(0,q,q^2,\ldots).
\]

L219's metric and the physical ellipse pencil become

\[
\boxed{
\begin{aligned}
R_B&=\mathscr H^*D_R\mathscr H
+P_{\mathcal K}D_LP_{\mathcal K},\\
\Xi_B&=S_B+c(I+P_0)L^*(I+E_B),\\
A_B&=\phi_c(\Xi_B).
\end{aligned}}                                    \tag{3}
\]

Consequently L250's edge-deleted mass is exactly

\[
\boxed{
\mu^\circ_{B,k}
=-m+\operatorname{tr}_{\mathcal K}
\left\{
(R^\circ_{B,k}-A_B^*P_0A_B)^{-1}
(R^\circ_{B,k}-A_B^*R^\circ_{B,k}A_B)
\right\},}                                        \tag{4}
\]

where

\[
R^\circ_{B,k}=R_B-q^k[q^k]R_B.                    \tag{5}
\]

Equations (1)--(5) are exact.  They reduce the remaining theorem to
one model-space shift covariance and expose why the stronger
whole-series statement is unnecessary.

## 2. Exact delay split

Assume

\[
B_1=\cdots=B_{k-1}=0
\]

and write

\[
B(z)=z^{k-1}\widetilde B(z),\qquad
\widetilde B_j=B_{j+k-1}.
\]

The standard model-kernel split from L221 gives the literal
orthogonal decomposition

\[
\boxed{
\mathcal K_{zB}
=\bigoplus_{j=0}^{k-2}z^j\mathbb C^m
\ \oplus\
z^{k-1}\mathcal K_{z\widetilde B}.}               \tag{6}
\]

Thus the first summand is exactly the clean wandering chain, while
the second is the arbitrary deflated grade-one model.  Under (6),
(2)--(3) reproduce L235's arrowhead pencil, L221's retained metric,
and L245's paired Green columns.  In particular, the tail is not a
new or approximate realization: it is the canonical model of
\(\widetilde B\).

Define the one-delay version by \(B=z\widetilde B\).  The desired
all-grade theorem is now precisely the associated coefficient
recursion

\[
\boxed{
[c^{2k}]\mu^\circ_{B,k}
=[c^{2k-2}]\mu^\circ_{\widetilde B,k-1}.}          \tag{7}
\]

Together with the corresponding earlier vanishings, iterating (7)
and using L256 gives

\[
[c^{2k}]\mu^\circ_{B,k}
=[c^2]\mu^\circ_{\widetilde B^{(k-1)},1}
=4\|B_k\|_F^2.                                    \tag{8}
\]

Conversely, the active part of (7) is exactly the scalar
one-delay covariance still missing from A194.  Therefore (7) is a
reformulation, not a claimed proof.

Separate floating diagnostics show the stronger relation
\(\mu^\circ_{B,k}-c^2\mu^\circ_{\widetilde B,k-1}
=O(c^{2k+2})\) in grades two through four, with a nonzero next-even
coefficient.  This is consistent with L235, but it is not used here.
The correct induction transports one associated face and nothing
more.

## 3. Proof of the canonical formulas

L236 proves that \(\mathcal O_L\) is an isometry and L254 identifies
its range with \(\mathcal K_{zB}\).  It is therefore unitary from the
finite state space onto \(\mathcal K\).

The forward intertwining

\[
\mathcal O_LS^*=L^*\mathcal O_L
\]

gives the first line of (2).  L255's opposite intertwining and defect
actions give

\[
\begin{aligned}
\mathcal O_LS
&=L\mathcal O_L-L\mathscr H^*P_0\mathcal O_R,\\
\mathcal O_LF&=P_0\mathcal O_L,\\
\mathcal O_LE&=\mathscr H^*P_0\mathcal O_R.
\end{aligned}
\]

Since

\[
\mathcal O_R=\mathscr H\mathcal O_L,
\]

these are exactly (2).  The equality
\(L(I-E_B)=P_{\mathcal K}L\) is the usual compressed-shift formula
and also follows by taking the adjoint of
\(L^*|_{\mathcal K}\).

Conjugating L236's metric pullback by \(\mathcal O_L\) gives

\[
\mathcal O_LR\mathcal O_L^*
=\mathscr H^*D_R\mathscr H
+P_{\mathcal K}D_LP_{\mathcal K},
\]

which is the first line of (3).  Substitution of (2) into
\(\Xi=S+c(I+F)S^*(I+E)\) gives its second line.  Holomorphic
functional calculus commutes with unitary conjugacy, proving the
third line.  L250 then gives (4)--(5).

Finally, \(zB=z^k\widetilde B\), so the elementary identity

\[
\mathcal K_{z^r\Theta}
=\bigoplus_{j=0}^{r-1}z^j\mathbb C^m
\oplus z^r\mathcal K_\Theta
\]

with \(r=k-1\) and \(\Theta=z\widetilde B\) proves (6).

## 4. Exact finite covariance audit

The exact rational word checker constructs both sides of (7) in the
same partial-isometry quotient.  The embedded tail uses

\[
I_{\rm tail}=I-F,\qquad
S_1=S(I-F),\qquad
E_1^{\rm tail}=E,\qquad
F_1^{\rm tail}=SFS^*.
\]

For grades two through five, the active cyclic difference reduces
identically to zero.  The unreduced quotient differences contain
\(36,82,159,268\) words, so this is not a comparison of two stored
target values.

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_volume_delay_covariance.py \
  --maximum-grade 5 \
  --output \
  experiments/repeated_crabb_volume_delay_covariance_s70225.jsonl
```

The tracked dataset SHA-256 is

```text
62981378ee78b5aa498def5882672976cd9c5b8e71d1046761b1690a3e259d9e
```

This audit proves neither (7) for arbitrary \(k\) nor Crouzeix's
conjecture.  Its value is adversarial: the same exact engine computes
the original and deflated faces rather than comparing two stored
copies of the target \(4\|B_k\|_F^2\).

## 5. Remaining proof gate

Prove (7) directly on (6), using:

1. L244 to remove the word-free half-line coisometric sector;
2. L245/L251 to keep the two orientations of every chain port paired;
3. L243 to exclude inverse-kernel powers above one in the active
   window; and
4. L254 to identify the surviving shifted-left complement.

This is now a one-frame block trace theorem.  No new state words,
metric roots, or choice of realization are permitted or required.
