# The boundary slack has a canonical repair and a compact first face

> **Closure note (2026-07-25).**  L283 proves the delayed residual
> lift in every completely delayed grade.  The canonical repair
> therefore has faces \(+B_k^*B_k\) and \(-12B_kB_k^*\) there and
> remains analytic through transfer-rank changes.  A172 still forbids
> iterating the unmodified repair through partial flags: mixed odd
> coefficients require the L230--L234/A178 preparation mechanism.

## 1. Result (L227, 2026-07-24)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

put \(Q=I-E=S^*S\), and let

\[
B_1=W^*S^*V.
\]

For L219's boundary metric, let

\[
H(c)=P_{\rm bl}(c)-\widehat T(c)^*
                    P_{\rm bl}(c)\widehat T(c)
\]

be its balanced Stein slack.  Schur-compress \(H(c)\) away from the
right defect \(V\), as in L225, and denote its second coefficient by
\(K_2\).  Then the entire state-space matrix collapses to

\[
\boxed{
K_2=2F-(S^*)^2S^2F-F(S^*)^2S^2.}                 \tag{1}
\]

Its upper Stein endpoint is the exact matrix Gram

\[
\boxed{
W^*{\cal G}_S(K_2)W=2B_1B_1^*.}                  \tag{2}
\]

Thus L225's scalar identity
\(\operatorname {tr}K_2=2\|B_1\|_F^2\) is the trace of a stronger
matrix endpoint identity.  No copy matrices are commuted.

There is also an exact nonlinear repair behind (1).  For any small
real \(c\), put

\[
\begin{aligned}
G(c)&=V^*H(c)V,\\
{\cal K}(c)
&=H(c)-H(c)V\,G(c)^{-1}V^*H(c),\\
X(c)&=-{\cal G}_{\widehat T(c)}({\cal K}(c)).
\end{aligned}                                      \tag{3}
\]

Then \(G(c)>0\), \(X(c)\) is analytic, and

\[
\boxed{
\begin{aligned}
&(P_{\rm bl}+X)
-\widehat T^*(P_{\rm bl}+X)\widehat T\\
&\qquad
=HVG^{-1}V^*H\succeq0.
\end{aligned}}                                     \tag{4}
\]

Consequently \(P_{\rm bl}+X\) is an exact local contraction metric.
This statement alone does **not** bound its condition number by four.

At grade one, however, (1)--(2) determine its two metric faces.
The leading lower physical face remains

\[
+c^2B_1^*B_1,                                     \tag{5}
\]

while the leading upper physical face becomes

\[
\boxed{-12c^2B_1B_1^*.}                           \tag{6}
\]

The boundary metric supplies \(-4B_1B_1^*\), and the canonical repair
supplies the additional \(-8B_1B_1^*\).

L227 makes the desired delay theorem sharper: it is enough to prove
that the first delayed Schur residual is the lift of (1).  L216 then
transports (2), so the same canonical repair gives the all-grade
faces \(+B_k^*B_k\) and \(-12B_kB_k^*\).  That delay covariance is
still the unproved part of L225.

## 2. The second slack Schur coefficient

Use L225's notation

\[
\begin{aligned}
\dot S&=J-S^3,\qquad J=(I+F)S^*(I+E),\\
\ddot S
&=2S-(S^2J+SJS+JS^2)+S^5,\\
R_2&=SFS^*-S^*ES.
\end{aligned}                                      \tag{7}
\]

The first two slack coefficients are

\[
\begin{aligned}
H_1&=-(\dot S^*S+S^*\dot S),\\
H_2&=R_2-\ddot S^*S-S^*\ddot S-\dot S^*\dot S
     -S^*R_2S.
\end{aligned}                                      \tag{8}
\]

Since \(V^*H_0V=I\), Schur complementation gives

\[
K_2=QH_2Q-AA^*,\qquad
A=QH_1V=-S^*\dot SV.                              \tag{9}
\]

Introduce the fixed grade-one forcing

\[
\widehat F_2
=AA^*+\ddot S^*S+S^*\ddot S+\dot S^*\dot S.       \tag{10}
\]

Because \(A=QA\), equations (8)--(10) give

\[
\boxed{
K_2=Q\{R_2-S^*R_2S-\widehat F_2\}Q.}              \tag{11}
\]

Substitute (7), retain every factor order, and use only

\[
S^*S=Q,\quad SS^*=I-F,\quad
SE=0,\quad S^*F=0,\quad EF=FE=0.
\]

The direct reduced word expansion is

\[
\begin{aligned}
K_2={}&2S^*S-2(S^*)^2S^2
-2S^*S^2(S^*)^2S\\
&+(S^*)^2S^3(S^*)^2S
 +S^*S^2(S^*)^3S^2.                               \tag{12}
\end{aligned}
\]

This five-word form has a short collapse.  Put
\(Q_2=(S^*)^2S^2\).  The two defect telescoping identities

\[
\begin{aligned}
S^*S^2(S^*)^2S&=Q-F,\\
(S^*)^2S^3(S^*)^2S&=Q_2(I-F)
\end{aligned}                                      \tag{13}
\]

follow from

\[
I-S^2(S^*)^2=F+SFS^*
\]

and \(QF=FQ=F\).  The last word in (12) is the
adjoint of the second line of (13).  Substitution in (12) gives
exactly (1).

## 3. Exact endpoint response

For a state forcing \(Y\), write

\[
\Phi(Y)=W^*{\cal G}_S(Y)W
=\sum_{n\ge0}W^*(S^*)^nYS^nW.                    \tag{14}
\]

Because \(W^*S=0\),

\[
FS^nW=0\qquad(n\ge1).
\]

It follows termwise that

\[
\begin{aligned}
\Phi(F)&=I,\\
\Phi(Q_2F)&=W^*Q_2W,\\
\Phi(FQ_2)&=W^*Q_2W.
\end{aligned}                                      \tag{15}
\]

The first forward step from \(W\) has no right-defect loss, while the
second loses exactly \(V^*SW=B_1^*\).  Therefore

\[
\begin{aligned}
W^*Q_2W
&=(S^2W)^*S^2W\\
&=I-(V^*SW)^*(V^*SW)
=I-B_1B_1^*.                                      \tag{16}
\end{aligned}
\]

Apply \(\Phi\) to (1) and insert (15)--(16):

\[
\Phi(K_2)
=2I-2(I-B_1B_1^*)
=2B_1B_1^*.
\]

This proves (2), including the left-Gram orientation.

## 4. Canonical exact Stein repair

The construction (3)--(4) is an elementary Schur-factor identity.
Since \(H(0)=VV^*\), the pivot \(G(c)\) remains positive and
invertible.  The residual \({\cal K}\) is Hermitian and annihilates
\(V\) on both sides.  Stability of \(\widehat T(c)\) makes its Stein
inverse analytic, and

\[
X-\widehat T^*X\widehat T=-{\cal K}.
\]

Thus the repaired slack is

\[
H-{\cal K}=HVG^{-1}V^*H
=(HVG^{-1/2})(HVG^{-1/2})^*,
\]

which proves (4).  Also \(X(0)=0\), so the repaired metric is positive
for all sufficiently small \(c\).

At order two,

\[
[c^2]X=-{\cal G}_S(K_2).                           \tag{17}
\]

Since \(K_2V=0\) and \(SV=0\),

\[
V^*[c^2]XV=0.                                     \tag{18}
\]

At the upper defect, (2) and the physical factor
\(P^{1/2}W=2W\) give

\[
4W^*[c^2]XW=-8B_1B_1^*.                           \tag{19}
\]

L225's boundary faces are \(+B_1^*B_1\) below and
\(-4B_1B_1^*\) above.  Equations (18)--(19) prove (5)--(6).

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_boundary_slack_repair.py \
  --output \
  experiments/repeated_crabb_boundary_slack_repair_s70224.jsonl
```

The deterministic audit covers defect dimensions one through three
and twelve unstructured pure partial isometries.  It independently
constructs the slack Schur coefficient, checks (1)--(2), solves the
canonical repair, and verifies the repaired faces (5)--(6).
Equations (7)--(19), not the floating audit, prove L227.  The tracked
dataset SHA-256 is
`193dca291f79120c1eea4605081b8db9bc542962dddf097088de109f84761801`.
