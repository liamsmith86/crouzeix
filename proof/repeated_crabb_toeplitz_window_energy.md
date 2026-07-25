# The last radial moment is one causal Toeplitz-window energy

## 1. Result (L253, 2026-07-25)

Retain L252's finite pure partial isometry and matrix-inner transfer

\[
B(z)=\sum_{h\geq1}B_hz^h,\qquad
Q_j=(S^*)^jS^j .
\]

Let \({\cal T}_B\) be multiplication by \(B\) on
\(H^2(\mathbb C^m)\), written as the causal block Toeplitz matrix

\[
[{\cal T}_B]_{i,j}=
\begin{cases}
B_{i-j},&i\geq j,\\
0,&i<j,
\end{cases}
\qquad B_0=0.
\]

Let \(P_k\) project onto Hardy coordinates \(0,\ldots,k\).  Then

\[
\boxed{
\begin{aligned}
\operatorname {tr}
\{Q_{k+2}-(k+2)Q_1+(k+1)I\}
&=\|P_k{\cal T}_BP_k\|_{\rm HS}^2\\
&=\sum_{h=1}^{k}(k+1-h)\|B_h\|_F^2\\
&=(k+1)m-\operatorname {tr}(P_kP_{{\cal K}_B}P_k),
\end{aligned}}                                    \tag{1}
\]

where

\[
{\cal K}_B=H^2(\mathbb C^m)\ominus
B H^2(\mathbb C^m)
\]

and \(P_{{\cal K}_B}\) is its model-space projection.

Under the complete delay

\[
B_1=\cdots=B_{k-1}=0,                             \tag{2}
\]

the finite Toeplitz window has exactly one nonzero block:

\[
\boxed{
P_k{\cal T}_BP_k
=P_{\{k\}}B_kP_{\{0\}},\qquad
\|P_k{\cal T}_BP_k\|_{\rm HS}^2=\|B_k\|_F^2.}     \tag{3}
\]

Thus A198's remaining cyclic congruence has an intrinsic Hardy-space
target:

\[
\boxed{
[c^{2k}]
\left\{-m+\operatorname {tr}(N^{-1}H)\right\}
=4\|P_k{\cal T}_BP_k\|_{\rm HS}^2.}               \tag{4}
\]

Equation (1)--(3) is proved exactly in every grade.  Equation (4)
remains open: L253 identifies the scalar energy that L251's paired
channels must produce, but does not identify the physical volume with
that energy.

This formulation avoids two previous traps.  It is a positive
**scalar trace energy**, not the false positive-matrix strengthening
excluded by L250; and the factor two is still part of the complete
physical analytic port of L251, not a detached multiplier as in A195.

## 2. Coefficient proof

The finite causal window contains \(B_h\) in positions

\[
(h,0),(h+1,1),\ldots,(k,k-h).
\]

There are \(k+1-h\) such positions.  Orthogonality of matrix units in
the Hilbert--Schmidt inner product therefore gives

\[
\|P_k{\cal T}_BP_k\|_{\rm HS}^2
=\sum_{h=1}^{k}(k+1-h)\|B_h\|_F^2.                \tag{5}
\]

L252 proves

\[
\operatorname {tr}Q_{k+2}
=n-(k+2)m+
\sum_{h=1}^{k}(k+1-h)\|B_h\|_F^2,
\]

while \(\operatorname {tr}Q_1=n-m\).  Substitution cancels the
dimension and defect-rank terms and proves the first two expressions
in (1).

If (2) holds, every block in the window vanishes except the unique
copy of \(B_k\) at row \(k\), column \(0\).  This proves (3).

## 3. Model-space proof

Matrix innerness makes \({\cal T}_B\) an isometry, with range
\(BH^2(\mathbb C^m)\).  Consequently

\[
{\cal T}_B{\cal T}_B^*
=I-P_{{\cal K}_B}.                                \tag{6}
\]

Causality gives

\[
P_k{\cal T}_B=P_k{\cal T}_BP_k.                  \tag{7}
\]

Use (6)--(7) and cyclicity of the finite trace:

\[
\begin{aligned}
\|P_k{\cal T}_BP_k\|_{\rm HS}^2
&=\operatorname {tr}
\{P_k{\cal T}_BP_k{\cal T}_B^*P_k\}\\
&=\operatorname {tr}
\{P_k{\cal T}_B{\cal T}_B^*P_k\}\\
&=(k+1)m-\operatorname {tr}(P_kP_{{\cal K}_B}P_k).
\end{aligned}
\]

This proves the last expression in (1).  It also places the target in
the exact model space split by L221 and the two Hardy frames of L236.

## 4. Consequence for the live proof

The arbitrary-grade task is no longer to discover the output of the
paired trace.  Its required output is the finite-section leakage in
(4).  A suitable continuation of L251 should:

1. conjugate each closed paired channel into L236's Hardy frames;
2. use L244 to remove the full-space coisometric background;
3. show that the remaining chain truncation is exactly
   \(P_k{\cal T}_BP_k\), rather than a separately weighted
   \(\mathcal B^\sharp\) row; and
4. retain the physical factor \(4\) only after this closed-window
   identification.

Proving those four operations would establish A198 and A194.  L253
proves only the target equivalences used at the last step.

## 5. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_toeplitz_window_energy.py \
  --output \
  experiments/repeated_crabb_toeplitz_window_energy_s70224.jsonl
```

The checker uses general and completely delayed matrix-inner
colligations.  It compares the state radial trace, the explicit
Toeplitz finite section, the coefficient sum, and the model-kernel
compression for window sizes one through six.  The equations above,
not the floating audit, prove L253.
The tracked dataset regenerates byte for byte with SHA-256
`ab4fe9e48646843ac26772003a3c89a5847f68c9f1b0694a35c30eae656d843d`.
