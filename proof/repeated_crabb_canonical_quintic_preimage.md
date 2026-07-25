# A polynomial quintic preparation clears the second transfer flag

## 1. Result (L233, 2026-07-24)

Retain the balanced equality colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad Q=I-E,
\]

and the transfer coefficients

\[
B_j=W^*(S^*)^jV.
\]

Insert L230's cubic and L232's quartic columns into the exact
canonical positive Stein factor before expanding the fifth order.
There is a universal polynomial perpendicular column \(C_5\) such
that the resulting physical fifth upper endpoint satisfies

\[
\boxed{
U^*[c^5]{\cal U}_{\rm final}U=0,\qquad
\operatorname{ran}U=\ker B_1^*\cap\ker B_2^* .}       \tag{1}
\]

This is a full partial-flag statement, not merely the identity on
colligations with \(B_1=B_2=0\).  The column has no inverse,
pseudoinverse, flag projection, or selected singular vectors, and

\[
\boxed{\|C_5\|\le48.}                                \tag{2}
\]

Thus it remains bounded and real analytic through every rank change.
The fifth endpoint need not vanish away from the flag.  Its
off-diagonal coupling to the positive quartic range contributes to
the effective sextic Schur face, so L233 does not by itself prove the
next even face positive.

## 2. Explicit twenty-term column

Write \(a=S^*\) and \(s=S\), with products read from left to right.
Set

\[
\begin{aligned}
{\cal C}_5={}&
\tfrac12s^4a^2+\tfrac32a^3s^3a^2
+\tfrac32a^2s^3a^3-\tfrac12s^5a^3\\
&+3a^4s^4a^2+2a^3s^4a^3
+\tfrac52a^7s^3a^2-3a^6s^3a^3\\
&-a^4s^3a^5-a^3s^3a^6
+\tfrac32a^3s^7a^2-\tfrac12a^2s^3a^7\\
&+\tfrac12a^2s^7a^3+\tfrac12sa^3s^6a^2
-\tfrac32s^2a^5s^3a^2\\
&+\tfrac32s^2a^4s^3a^3+\tfrac12s^4a^3s^3a^2\\
&+17a^3s^3a^3s^3a^2
-\tfrac12a^2s^3a^4s^3a^2
+\tfrac{15}{2}a^2s^3a^3s^3a^3 .
\end{aligned}
\]

The required column is

\[
\boxed{C_5=Q{\cal C}_5V.}                            \tag{3}
\]

The coefficient \(\ell^1\)-norm of the braces is \(48\), proving
(2).  The leading \(Q\) proves \(V^*C_5=0\).  Consequently the
same-order lower endpoint motion of this free column is zero.

## 3. Exact Stein certificate

Let \({\cal F}_5\) be the exact fifth Stein forcing obtained by:

1. factoring L227's canonical slack in the Hermitian defect gauge;
2. adding \(c^3C_3+c^4C_4\) to that factor;
3. reconstructing the slack through order five; and
4. completing the variable metric through L230's \(Z_3\) and L232's
   \(Z_4\).

The exact checker obtains 158 reduced words in \({\cal F}_5\).
Put

\[
{\cal P}_5=C_5V^*.
\]

There is an explicit Hermitian 55-term polynomial \(Z_5\), with 110
reduced words, satisfying

\[
\boxed{
{\cal F}_5+{\cal P}_5+{\cal P}_5^*
=Z_5-S^*Z_5S.}                                      \tag{4}
\]

All coefficients of \(Z_5\) are integers or half-integers.  The
checker stores the complete certificate verbatim in
`quintic_coboundary_witness()`.  Exact rational word reduction leaves
zero words in (4), zero words in \(Z_5-Z_5^*\), and zero words in
\(E{\cal P}_5\).

## 4. Exact partial-flag factorization

Let

\[
R=I-F-\frac23E,\qquad M_j=[c^j]\widetilde P,
\]

and define the fifth upper Schur cross

\[
{\cal S}_5=M_2RM_3F+M_3RM_2F.
\]

The balanced endpoint lift after (3) is

\[
{\cal H}_5=FZ_5F+F{\cal S}_5.                       \tag{5}
\]

Define the two delayed channel columns

\[
G_1=ESF,\qquad G_2=ES^2F
\]

and the left multipliers

\[
\begin{aligned}
R_1&=F\left\{\frac72(S^*)^4S^3(S^*)^2
              +\frac52(S^*)^2S^3(S^*)^4\right\},\\
R_2&=F(S^*)^2S^4(S^*)^2.
\end{aligned}
\]

A second exact word identity is

\[
\boxed{
{\cal H}_5=R_1G_1+R_2G_2+(R_1G_1)^*+(R_2G_2)^*.}    \tag{6}
\]

This compact six-term Hermitian factorization is stronger than the
complete-delay reduction used to discover the column.

If \(x\in\ker B_1^*\cap\ker B_2^*\), then

\[
G_jWx=ES^jWx=VB_j^*x=0,\qquad j=1,2.
\]

Therefore \(x^*W^*{\cal H}_5Wx=0\), including all polarized
compressions inside the joint flag.  Physical congruence gives

\[
[c^5]{\cal U}_{\rm final}=-4W^*{\cal H}_5W,
\]

which proves (1).

## 5. Scope and the next gate

L232 and L233 now give the ordered upper expansion

\[
[c^4]{\cal U}_{\rm final}\succeq0,\qquad
\ker[c^4]{\cal U}_{\rm final}
=\ker B_1^*\cap\ker B_2^*,
\]

followed by a fifth coefficient whose compression to that kernel is
zero.  This removes a necessary odd-order obstruction through the
second transfer flag.

The next required calculation has two parts:

1. on \(B_1=B_2=0\), compute the sixth face in the new
   cubic--quartic--quintic gauge and compare it with L215/L228's
   grade-three prediction; and
2. at partial rank changes, include the sextic Schur cost generated
   by the fifth off-diagonal block against L232's positive quartic
   range.

The first is polynomial word algebra.  The second is the genuinely
rank-changing quantitative gate and must not be replaced by the
complete-delay special case.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_canonical_quintic_preimage.py \
  --output \
  experiments/repeated_crabb_canonical_quintic_preimage_s70224.jsonl
```

The checker proves (4) and (6) with exact rational reduction, then
reconstructs the full prepared factor and both physical endpoints on
eight unstructured, sixteen rank-changing, and three complete
double-delay colligations.  All 27 records pass.  The largest direct
fifth-endpoint reconstruction error is below \(1.4\times10^{-12}\);
the largest joint-flag compression is below
\(1.2\times10^{-12}\).

The tracked SHA-256 is
`b5558b357fecfdfc4a92697635e89b82873dabe61262302b1696841c878a40e3`.
