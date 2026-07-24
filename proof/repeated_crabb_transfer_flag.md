# The higher transfer flag has an explicit coboundary preimage

## 1. Result (L208, 2026-07-24)

L207 removes the grade-one rank jump.  The same state-space mechanism
extends to every genuine transfer coefficient once the earlier left
rows are removed.

Let

\[
\begin{aligned}
I-S^*S&=VV^*,\qquad I-SS^*=WW^*,\\
Q&=I-VV^*,\\
B_k&=W^*(S^*)^kV\qquad(k\ge0)
\end{aligned}                                      \tag{1}
\]

for a spectrally stable pure partial isometry with equal orthogonal
defects.  Thus \(B_0=0\), and define the transfer channel

\[
 {\mathfrak C}(K)
 =\sum_{n\ge1}B_nKB_n^*
 =W^*{\cal G}_S(VKV^*)W.                           \tag{2}
\]

For a grade \(k\ge1\), let \(P=P^*=P^2\) be any copy projection such
that

\[
 PB_j=0\qquad(1\le j<k).                            \tag{3}
\]

Then for every real scalar \(\eta\), the physical L204 endpoint map
has the explicit compressed preimage

\[
\boxed{
\begin{aligned}
 \widehat C_{k,P,\eta}
 &=-\frac{\eta}{8}QS^kWPB_k,\\
 C_{k,P,\eta}
 &=P_{\rm met}^{1/2}\widehat C_{k,P,\eta},\\
 P{\cal M}_T(C_{k,P,\eta})P
 &=\eta P\{
 {\mathfrak C}(B_k^*PB_k)-B_kB_k^*
 \}P .
\end{aligned}}                                    \tag{4}
\]

Here \(P_{\rm met}=2I-VV^*+2WW^*\) is the physical equality metric;
the letter \(P\) without a subscript always denotes the copy-space
flag projection.

For \(k=1\), \(P=I\), and \(\eta=28\), formula (4) is exactly L207:

\[
\widehat C_{1,I,28}=-\frac72QSWB_1.
\]

Thus a higher reflected face of the channel-minus-left-Gram form in
(4) has no new pointwise range or rank-jump obstruction on its active
flag.  What remains open is to derive the **actual physical
higher-grade target** from the complete Faber/Riemann/metric jet and
to organize the projections analytically through changing ranks.
L208 does not assume that target identification.

## 2. Arbitrary higher intertwining defect

Let \(Y=Y^*\) be arbitrary, and let

\[
 H-SHS^*=WYW^*,\qquad
 A=V^*HV,\qquad R=QHV.                              \tag{5}
\]

The Stein series immediately gives

\[
 HW=WY.                                             \tag{6}
\]

Multiplying (5) on the right by \(S\), using \(W^*S=0\), and
decomposing \(HV=VA+R\) gives

\[
 HS-SH=-SRV^*.                                     \tag{7}
\]

Taking adjoints,

\[
 HS^*=S^*H+VR^*S^*.                                \tag{8}
\]

Iteration of (8), without commuting any factors, yields

\[
H(S^*)^k
=(S^*)^kH+
\sum_{j=0}^{k-1}(S^*)^jVR^*(S^*)^{k-j}.            \tag{9}
\]

Multiply (9) by \(V\), compress on the left by \(W^*\), and use
\(HV=VA+R\) and \(B_0=W^*V=0\).  This proves the exact all-grade
identity

\[
\boxed{
\begin{aligned}
YB_k-B_kA
={}&W^*(S^*)^kR\\
&+\sum_{j=1}^{k-1}
B_jR^*(S^*)^{k-j}V .
\end{aligned}}                                    \tag{10}
\]

At grade one the sum is empty, recovering L207's intertwining
defect.  At higher grades it records every lower transfer channel and
explains why a blind replacement \(B_1\mapsto B_k\) would be false
before flag compression.

## 3. Earlier rows disappear on the flag

Under (3), left multiplication of (10) by \(P\) kills every summand
in its second line:

\[
\boxed{
P(YB_k-B_kA)=PW^*(S^*)^kR.}                        \tag{11}
\]

No commutation between \(P,Y,B_k\), or \(A\) is used.  This is the
precise recursive simplification supplied by the common left kernel

\[
\operatorname{ran}P
\subseteq\bigcap_{j<k}\ker B_j^*.
\]

## 4. Dual proof of the explicit preimage

It suffices to test the \(P\)-compression against arbitrary Hermitian
copy matrices \(Y=PYP\).

First, Stein adjointness in (2) gives

\[
\begin{aligned}
\operatorname{tr}\{Y{\mathfrak C}(B_k^*PB_k)\}
 &=\operatorname{tr}(AB_k^*PB_k).
\end{aligned}                                      \tag{12}
\]

Therefore

\[
\begin{aligned}
&\operatorname{tr}Y\,
 P\{{\mathfrak C}(B_k^*PB_k)-B_kB_k^*\}P\\
&\qquad
=-\operatorname{Re}\operatorname{tr}
\{B_k^*P(YB_k-B_kA)\}\\
&\qquad
=-\operatorname{Re}\operatorname{tr}
\{B_k^*PW^*(S^*)^kR\},                             \tag{13}
\end{aligned}
\]

where the last equality is (11).  The real part is harmless and makes
the Hermitian pairing explicit.

For a balanced perpendicular column \(\widehat C\), L207's physical
endpoint normalization gives

\[
\operatorname{tr}\{Y{\cal M}_T(P_{\rm met}^{1/2}\widehat C)\}
=8\operatorname{Re}\operatorname{tr}(\widehat C^*R). \tag{14}
\]

Substituting the first line of (4) into (14) gives \(\eta\) times
(13).  Since every Hermitian test supported on \(P\) agrees, the
compressed matrix identity (4) follows.

The column is admissible because

\[
V^*C_{k,P,\eta}
=V^*P_{\rm met}^{1/2}Q(\cdots)=0.                  \tag{15}
\]

This completes the proof.

## 5. Interpretation and remaining gate

The nested exact copy spaces

\[
{\cal K}_{k-1}
=\bigcap_{j=1}^{k-1}\ker B_j^*                    \tag{16}
\]

form the natural left reflected flag.  L201's invertible terminal
\(B_L\) makes \({\cal K}_L=\{0\}\) near the Crabb apex.  Equation
(10) shows that all apparent higher-grade order contamination comes
from coefficients already removed by this flag, while (4) supplies a
polynomial preimage for the surviving channel coboundary.

Two debts remain before this becomes a higher-grade elliptic theorem.

1. The complete prepared Faber/Riemann jet must be shown to have, on
   \({\cal K}_{k-1}\), the channel-minus-left-Gram target appearing
   in (4), with the correct scalar coefficient.
2. Exact kernel projections can jump rank.  A local proof must use
   analytic Schur-orthogonal flags along arcs, rather than insert a
   discontinuous Moore--Penrose projection into the metric formula.

L208 removes the state-space range problem after those two
identifications; it does not silently assume them.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_transfer_flag.py \
  --output experiments/repeated_crabb_transfer_flag_s70224.jsonl
```

The checker verifies (10) for grades one through five on
unstructured partial isometries.  It then uses independently gauged
direct sums of delay channels to produce nontrivial flags through
grade six and verifies (3)--(4).  These delayed models test changing
flag dimensions and nonnormal transfer matrices; the proof itself is
(5)--(15).  The tracked dataset SHA-256 is
`913902774e66fa0fee8238242e9e51736009d081c18e9294d8b84e996f8b7c5d`.
