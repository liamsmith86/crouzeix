# The delayed slack lives between two Hardy frames

## 1. Result (L236, 2026-07-24)

Retain a finite pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

with transfer coefficients

\[
B_j=W^*(S^*)^jV,\qquad B_0=0.
\]

Define the right and left Hardy analysis maps

\[
\begin{aligned}
({\cal O}_Rx)_n&=V^*S^nx,\\
({\cal O}_Lx)_j&=W^*(S^*)^jx,
\end{aligned}
\qquad n,j\ge0.                                    \tag{1}
\]

They are isometries from the state space into
\(\ell^2(\mathbb N_0;\mathbb C^m)\):

\[
\boxed{{\cal O}_R^*{\cal O}_R
={\cal O}_L^*{\cal O}_L=I.}                        \tag{2}
\]

Their cross Gram is exactly the block Hankel matrix of the inner
transfer:

\[
\boxed{
\left[{\cal O}_R{\cal O}_L^*\right]_{n,j}
=B_{n+j}^*.}                                      \tag{3}
\]

Let \(L^*\) be the backward shift on the coefficient space.  The two
state shifts intertwine exactly:

\[
\boxed{
{\cal O}_RS=L^*{\cal O}_R,\qquad
{\cal O}_LS^*=L^*{\cal O}_L.}                      \tag{4}
\]

Finally, L219's boundary metric is the pullback of two fixed scalar
diagonal weights.  For \(q=c^2\), put

\[
\begin{aligned}
D_R(q)&=\operatorname{diag}
\left(1,\frac1{1+q},\frac1{1+q^2},\ldots\right),\\
D_L(q)&=\operatorname{diag}(0,q,q^2,\ldots).
\end{aligned}
\]

Then

\[
\boxed{
P_{\rm bl}^S(c)
={\cal O}_R^*D_R(q){\cal O}_R
+{\cal O}_L^*D_L(q){\cal O}_L.}                    \tag{5}
\]

Thus all dependence on the colligation is concentrated in two
isometric Hardy frames and their Hankel cross Gram; the metric
weights themselves are scalar and universal.

## 2. Delay is an anti-diagonal Hankel filtration

Assume the complete delay

\[
B_1=\cdots=B_{k-1}=0.                              \tag{6}
\]

Equation (3) says precisely that every positive Hankel
anti-diagonal below \(k\) vanishes.  The first active anti-diagonal
has constant block \(B_k^*\):

\[
\left[{\cal O}_R{\cal O}_L^*\right]_{n,j}
=B_k^*
\quad\text{when }n+j=k.                            \tag{7}
\]

This identifies L228's candidate face without a word expansion.
Let \(P_j\) denote coefficient projection onto Hardy row \(j\).  The
two defect orbits in L228 are

\[
E_1={\cal O}_R^*P_1{\cal O}_R,\qquad
F_{k-1}={\cal O}_L^*P_{k-1}{\cal O}_L.             \tag{8}
\]

Using the single Hankel cell \((1,k-1)\),

\[
\begin{aligned}
E_1F_{k-1}
&={\cal O}_R^*P_1
  ({\cal O}_R{\cal O}_L^*)
  P_{k-1}{\cal O}_L\\
&={\cal O}_R^*P_1 B_k^*P_{k-1}{\cal O}_L.
\end{aligned}                                      \tag{9}
\]

Therefore

\[
\boxed{
E_1F_{k-1}+F_{k-1}E_1
=\operatorname{Herm}
\left\{
{\cal O}_R^*P_1B_k^*P_{k-1}{\cal O}_L
\right\},}                                         \tag{10}
\]

where \(\operatorname{Herm}(X)=X+X^*\).

So L228's arbitrary-grade assertion has the following exact
interpretation:

> after the scalar ellipse/metric operation and right-defect Schur
> elimination, the first surviving associated coefficient is the
> Hermitian lift of the first active cell of the transfer Hankel
> matrix, with scalar multiplier one.

This removes the apparent free-word complexity and preserves the
matrix multiplication order automatically.

## 3. Proof of the two isometries and shifts

Purity and defect telescoping give

\[
\sum_{n\ge0}(S^*)^nES^n=I,\qquad
\sum_{j\ge0}S^jF(S^*)^j=I.
\]

The two left sides are respectively
\({\cal O}_R^*{\cal O}_R\) and
\({\cal O}_L^*{\cal O}_L\), proving (2).

For every \(n\ge0\),

\[
({\cal O}_RSx)_n=V^*S^{n+1}x
=({\cal O}_Rx)_{n+1}.
\]

This is the first identity in (4).  The second is its left-defect
analogue:

\[
({\cal O}_LS^*x)_j=W^*(S^*)^{j+1}x
=({\cal O}_Lx)_{j+1}.
\]

For the cross Gram, the adjoint of Hardy row \(j\) of
\({\cal O}_L\) maps a copy vector to \(S^jW\).  Hence

\[
\begin{aligned}
\left[{\cal O}_R{\cal O}_L^*\right]_{n,j}
&=V^*S^nS^jW\\
&=V^*S^{n+j}W
=B_{n+j}^*,
\end{aligned}
\]

which proves (3), including its orientation.

## 4. Proof of the metric pullback

Write

\[
E_n=(S^*)^nES^n
={\cal O}_R^*P_n{\cal O}_R,
\]

\[
F_j=S^jF(S^*)^j
={\cal O}_L^*P_j{\cal O}_L.
\]

Since \(\sum_{n\ge0}E_n=I\),

\[
\begin{aligned}
{\cal O}_R^*D_R(q){\cal O}_R
&=E_0+\sum_{n\ge1}\frac1{1+q^n}E_n\\
&=I-\sum_{n\ge1}\frac{q^n}{1+q^n}E_n.
\end{aligned}
\]

Likewise,

\[
{\cal O}_L^*D_L(q){\cal O}_L
=\sum_{j\ge1}q^jF_j.
\]

Adding the two identities gives exactly L219's definition of
\(P_{\rm bl}^S(c)\), proving (5).

## 5. Relation to L220, L221, L235, and L237

The range of \({\cal O}_R\) is a finite backward-shift-invariant
realization of the characteristic model space.  L220's ordered Schur
features give an orthogonal coordinate system for the corresponding
transfer model without requiring a convention-dependent identification
of \(B\) versus its coefficientwise adjoint.  L221's delay splitting
is the statement that the Hankel matrix in (3) has the zero
anti-diagonal prefix (6) and that the model space splits off the
corresponding monomial layers.

L235 and L237 supply the complementary physical information: after the
complete delay is removed, the balanced ellipse resolvent resums all
terminal round trips into one scalar insertion.  L236 identifies
what those grouped terms must act on—the universal diagonal weights
in (5) and the anti-diagonal filtration of the Hankel cross Gram.

What remains unproved is the load-bearing scalar/physical interface.
One must show that L235/L237's grouped direct-map and Schur-complement
operations are causal for this anti-diagonal filtration and that
their first derivative on the active cell is exactly the identity
in (10).  The full theta/ODE coefficients are required; L235
disproves the Newton-edge-only shortcut.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_hardy_two_frame.py \
  --output \
  experiments/repeated_crabb_hardy_two_frame_s70224.jsonl
```

The eleven deterministic records cover unstructured grade-one
colligations and complete delays through grade five, with defect
multiplicities two and three.  They audit both Stein isometries,
the shift intertwinings, Hankel cells, zero anti-diagonals, the
metric pullback, and the target-cell lift (10).  The exact equations
above prove the result; the audit checks state/copy orientations and
implementation.  The tracked SHA-256 is
`c771b6c491e1d80f08c18e3f0df958bb2ef7d8899576ed1f35f0b12950c0b5f6`.

## 7. Next gate

Work entirely in the two-frame realization:

1. conjugate L237's resummed terminal reflection by
   \({\cal O}_R,{\cal O}_L\);
2. retain the full scalar theta/ODE coefficients;
3. prove that coefficients below the first active Hankel
   anti-diagonal vanish after right-defect Schur elimination; and
4. prove that the active cell is returned with multiplier one.

This would prove L228's one-delay associated-graded recursion and
its all-grade anticommutator in one step.
