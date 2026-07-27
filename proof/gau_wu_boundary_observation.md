# The Gau--Wu loss is observed at the boundary

> **Status.**  The paired frame equations in Sections 1--2 are
> exact consequences of L127.  The rank and kernel statements in
> Section 3 are numerical evidence, not a theorem.  In particular,
> this note does not yet factor or dominate the L340 loss.

## 1. Paired transferred frames (L341/A288, 2026-07-26)

Use L336's shift coordinates

\[
 \phi=zf,\qquad S=S(\phi),\qquad p=f,\quad q=1,
\]

and let

\[
 T_\varepsilon=S+\varepsilon C+O(\varepsilon^2),
 \qquad
 Y_\varepsilon=f_\varepsilon(T_\varepsilon)
 =pq^*+\varepsilon Y_1+O(\varepsilon^2).
\]

For an orthonormal basis \(e_1,\ldots,e_{n-1}\) of
\({\cal K}_{f_\varepsilon}\), define the basis-independent positive
frame operators

\[
\begin{aligned}
 R_\varepsilon
 &=\sum_j e_j(T_\varepsilon)qq^*
                 e_j(T_\varepsilon)^*,\\
 L_\varepsilon
 &=\sum_j e_j(T_\varepsilon)^*pp^*
                 e_j(T_\varepsilon).
\end{aligned}                                      \tag{1}
\]

L127's two transfer identities give exactly

\[
\boxed{
\begin{aligned}
 R_\varepsilon-T_\varepsilon R_\varepsilon
 T_\varepsilon^*
 &=qq^*-Y_\varepsilon qq^*Y_\varepsilon^*,\\
 L_\varepsilon-T_\varepsilon^*L_\varepsilon
 T_\varepsilon
 &=pp^*-Y_\varepsilon^*pp^*Y_\varepsilon.
\end{aligned}}                                     \tag{2}
\]

At the Gau--Wu model,

\[
 R_0=I-pp^*,\qquad L_0=I-qq^*.                    \tag{3}
\]

Taking traces in (2) recovers the two L336 endpoint defects:

\[
\begin{aligned}
1-\|Y_\varepsilon q\|^2
 &=\operatorname{tr}\{(I-T_\varepsilon^*
 T_\varepsilon)R_\varepsilon\},\\
1-\|Y_\varepsilon^*p\|^2
 &=\operatorname{tr}\{(I-T_\varepsilon
 T_\varepsilon^*)L_\varepsilon\}.
\end{aligned}                                      \tag{4}
\]

This is the exact trace organization behind the observed
boundary-size loss.

## 2. Exact first frame response

Put

\[
 u=Y_1q,\qquad v=Y_1^*p,
\]

and write
\(R_\varepsilon=R_0+\varepsilon R_1+\cdots\) and
\(L_\varepsilon=L_0+\varepsilon L_1+\cdots\).
Differentiating (2) gives the two stable Stein equations

\[
\boxed{
\begin{aligned}
R_1-SR_1S^*
={}&C(I-pp^*)S^*+S(I-pp^*)C^*
      -up^*-pu^*,\\
L_1-S^*L_1S
={}&C^*(I-qq^*)S+S^*(I-qq^*)C
      -vq^*-qv^* .
\end{aligned}}                                     \tag{5}
\]

Thus the first frame forcing sees precisely the two endpoint
vectors \(u,v\), together with the normalized first operator
direction.  Since each frame in (1) has rank \(n-1\) near the base,
its moving kernel also gives the exact second-order identities

\[
 p^*R_2p=\|R_1p\|^2,\qquad
 q^*L_2q=\|L_1q\|^2.                              \tag{6}
\]

Equations (4)--(6) are the preferred starting point for deriving an
actual Gram for L340; expanding individual basis functions is
unnecessary.

## 3. Numerical boundary-observation law

Let \(E\) be a physical normal direction in the original Gau--Wu
coordinates, let \(s_E(\zeta)\) be its first support variation, and
perform L338's exact zero completion \(k=\kappa_E\).  Denote the
resulting first Blaschke image by \(Y_{1,E}^{\rm opt}\), and define

\[
\boxed{
 {\cal O}_\phi(E)=
 \left(
 s_E,\;
 Y_{1,E}^{\rm opt}q,\;
 (Y_{1,E}^{\rm opt})^*p
 \right).}                                        \tag{7}
\]

In twelve complete generic models, two per dimension \(4,\ldots,9\),

\[
\boxed{
\begin{aligned}
\operatorname{rank}(u,v)&=4n-8,\\
\operatorname{rank}{\cal O}_\phi&=6n-14,\\
\ker{\cal O}_\phi&=\ker{\cal G}_\phi,\\
\dim\ker{\cal O}_\phi&=(n-4)^2.
\end{aligned}}                                     \tag{8}
\]

The support coordinate contributes exactly \(2n-6\) new real
directions beyond the two endpoint vectors.  Moreover, on the
common kernel the **entire** optimized matrix \(Y_1^{\rm opt}\)
vanishes numerically, not merely its two endpoint vectors.

The loss therefore appears to factor through the canonical
boundary observation:

\[
 {\cal G}_\phi={\cal O}_\phi^*M_\phi{\cal O}_\phi
 \quad\text{for some }M_\phi\succeq0.             \tag{9}
\]

Equation (9) and positivity of \(M_\phi\) remain open.

## 4. Revised proof target

Use (5)--(6) and L339's support response
\(\gamma_E\) to derive \(M_\phi\) without taking a numerical square
root of \({\cal G}_\phi\).  The desired port estimate is

\[
 \|M_\phi^{1/2}{\cal O}_\phi(E)\|^2
 \leq2\|E\|_F^2+2\|\gamma_E\|_{L^2}^2.            \tag{10}
\]

This is exactly L338's remaining physical sign.  The two endpoint
frames must stay coupled; neither one-sided form is positive.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_boundary_observation.py \
  --output experiments/gau_wu_boundary_observation_s70224.jsonl
```

The checker rebuilds L340, applies its exact zero optimizer, forms
(7), and compares its null projector with the L340 loss kernel.  The
12-record dataset has SHA-256
`73b7d48403d8321ec89068ae6f492b5517c81153cf8b5383d37c99c32dda8e56`.
