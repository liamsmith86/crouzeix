# The physical ellipse pencil is a doubled Hardy-channel operator

## 1. Result (L255, 2026-07-25)

Retain L236's pure partial isometry, Hardy analysis maps, and transfer
Hankel cross Gram

\[
{\cal H}={\cal O}_R{\cal O}_L^*,\qquad
{\cal P}=P_{\{0\}}.
\]

Let \({\cal S}=L^*\) be the backward shift on
\(H^2(\mathbb C^m)\).  Besides L236's two forward intertwinings,
the opposite actions have exact Hankel corrections:

\[
\boxed{
\begin{aligned}
{\cal O}_RS&={\cal S}{\cal O}_R,&
{\cal O}_LS^*&={\cal S}{\cal O}_L,\\
{\cal O}_RS^*
&={\cal S}^*{\cal O}_R-{\cal S}^*{\cal H}{\cal P}{\cal O}_L,&
{\cal O}_LS
&={\cal S}^*{\cal O}_L-{\cal S}^*{\cal H}^*{\cal P}{\cal O}_R.
\end{aligned}}                                    \tag{1}
\]

The four defect actions are likewise

\[
\boxed{
\begin{aligned}
{\cal O}_RE&={\cal P}{\cal O}_R,&
{\cal O}_LF&={\cal P}{\cal O}_L,\\
{\cal O}_RF&={\cal H}{\cal P}{\cal O}_L,&
{\cal O}_LE&={\cal H}^*{\cal P}{\cal O}_R.
\end{aligned}}                                    \tag{2}
\]

Now put

\[
\Xi=S+c(I+F)S^*(I+E).
\]

Then the complete physical pencil has the exact doubled-Hardy
representation

\[
\boxed{
\begin{bmatrix}{\cal O}_R\Xi\\{\cal O}_L\Xi\end{bmatrix}
=
\begin{bmatrix}
X_{RR}&X_{RL}\\
X_{LR}&X_{LL}
\end{bmatrix}
\begin{bmatrix}{\cal O}_R\\{\cal O}_L\end{bmatrix},} \tag{3}
\]

where

\[
\boxed{
\begin{aligned}
X_{RR}
&={\cal S}
 c\{{\cal S}^*(I+{\cal P})
     +{\cal H}{\cal P}{\cal S}{\cal H}^*{\cal P}\},\\
X_{RL}
&=c\{-{\cal S}^*{\cal H}{\cal P}
       +{\cal H}{\cal P}{\cal S}\},\\
X_{LR}
&=-{\cal S}^*{\cal H}^*{\cal P}
 c(I+{\cal P}){\cal S}{\cal H}^*{\cal P},\\
X_{LL}
&={\cal S}^*+c(I+{\cal P}){\cal S}.
\end{aligned}}                                    \tag{4}
\]

At \({\cal H}=0\), the diagonal blocks are precisely L240's
half-line pencil and its adjoint:

\[
X_{RR}=\Xi_\infty
={\cal S}+c{\cal S}^*(I+{\cal P}),\qquad
X_{LL}=\Xi_\infty^*.                              \tag{5}
\]

Thus every dependence on the finite colligation is concentrated in
one Hankel channel:

1. the two off-diagonal blocks are linear in \({\cal H}\) or
   \({\cal H}^*\);
2. the sole diagonal correction is the ordered quadratic term
   \(c{\cal H}{\cal P}{\cal S}{\cal H}^*{\cal P}\); and
3. L236's metric is already the pullback of fixed diagonal weights
   \(D_R,D_L\) in these same two frames.

This is the missing exact interface among L236, L244, and L254.
It does **not** yet evaluate L250's volume.  The remaining A194
calculation is now a universal second reflected-channel response of
the block operator (3): apply the analytic functional calculus and
final-row whitening, use L244 on the diagonal background, and prove
that the surviving quadratic trace is four times
\({\cal T}_B{\cal T}_B^*
=I-\widehat{\cal O}_L\widehat{\cal O}_L^*\).
L243's reflection-power filtration is still required to exclude
higher channel powers at the active delayed grade.

## 2. Opposite intertwining proof

For \(n\ge1\),

\[
\begin{aligned}
V^*S^nS^*
&=V^*S^{n-1}(I-F)\\
&=V^*S^{n-1}-B_{n-1}^*W^*.
\end{aligned}                                     \tag{6}
\]

The row \(n=0\) is zero because \(SV=0\).  Since
\({\cal H}_{n,0}=B_n^*\), equation (6) is the third identity in
(1).  The fourth follows symmetrically from

\[
W^*(S^*)^nS
=W^*(S^*)^{n-1}-B_{n-1}V^*.                      \tag{7}
\]

The first two identities in (1) are L236.

For the defects, \(SV=0\) and \(S^*W=0\) give

\[
{\cal O}_RE={\cal P}{\cal O}_R,\qquad
{\cal O}_LF={\cal P}{\cal O}_L.
\]

Also

\[
[{\cal O}_RF]_n=B_n^*W^*,\qquad
[{\cal O}_LE]_j=B_jV^*,
\]

which, with the endpoint rows of the opposite frame, gives the last
two formulas in (2).

## 3. Pencil calculation

Write \(J=(I+F)S^*(I+E)\).  Equations (1)--(2), together with

\[
{\cal P}{\cal H}{\cal P}=B_0^*=0,                 \tag{8}
\]

give

\[
\begin{aligned}
{\cal O}_RJ
={}&\{
{\cal S}^*(I+{\cal P})
+{\cal H}{\cal P}{\cal S}{\cal H}^*{\cal P}
\}{\cal O}_R\\
&+\{
-{\cal S}^*{\cal H}{\cal P}
+{\cal H}{\cal P}{\cal S}
\}{\cal O}_L,                                     \tag{9}\\
{\cal O}_LJ
={}&(I+{\cal P}){\cal S}{\cal H}^*{\cal P}
  {\cal O}_R
+(I+{\cal P}){\cal S}{\cal O}_L.                 \tag{10}
\end{aligned}
\]

Add \({\cal O}_RS={\cal S}{\cal O}_R\) and the last identity in
(1), then multiply (9)--(10) by \(c\).  This gives exactly
(3)--(4).  No delay, truncation, commutativity, or functional
calculus has been used.

## 4. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_hardy_two_frame.py \
  --output experiments/repeated_crabb_hardy_two_frame_s70224.jsonl
```

The checker now additionally verifies (1)--(4) on general
colligations and complete delays through grade five, at defect
multiplicities two and three.  It uses finite Hardy prefixes only
where every displayed row is exact and discards the one artificial
terminal row introduced by truncating the backward shift.

All doubled-pencil errors are below \(3.8\times10^{-15}\).  The
tracked dataset regenerates with SHA-256
`5de72b1a2a1de56643f750dbe15d93645e04c0b81fafdca2832edcd0429186fa`.
