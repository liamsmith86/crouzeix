# A polynomial sextic preparation closes the third upper flag

## 1. Result (L234, 2026-07-24)

Retain the balanced equality colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad Q=I-E,
\]

and transfer coefficients

\[
B_j=W^*(S^*)^jV.
\]

Insert L230's cubic, L232's quartic, and L233's quintic
defect-factor columns before expanding order six.  There is a
universal perpendicular polynomial column

\[
\boxed{C_6=Q{\cal C}_6(S,S^*)V,\qquad \|C_6\|\le181} \tag{1}
\]

such that, on

\[
{\cal K}_2=\ker B_1^*\cap\ker B_2^*,
\]

the direct physical sixth upper face is exactly

\[
\boxed{
U^*[c^6]{\cal U}_{\rm final}U
=12\,U^*B_3B_3^*U,\qquad \operatorname{ran}U={\cal K}_2.} \tag{2}
\]

This direct face must still pay the Schur cost of L233's fifth
off-diagonal block against L232's positive quartic range.  After that
cost, the effective sixth face obeys the uniform bound

\[
\boxed{
[c^6]{\cal U}_{\rm eff}
\succeq \frac{215}{22}\,
U^*B_3B_3^*U\succeq0.}                             \tag{3}
\]

Thus the prepared upper metric passes the complete third ordered
transfer flag, including every partial rank change.  If the third
transfer also vanishes, (3) correctly passes control to the next
flag.  This is not an arbitrary-grade induction or a proof of the
full Crouzeix conjecture.

## 2. Exact polynomial certificate

The polynomial \({\cal C}_6\) has 53 integer, half-integer, and
quarter-integer terms.  Its coefficient \(\ell^1\)-norm is \(181\),
which proves (1).  The complete exact table is stored in

`experiments/repeated_crabb_canonical_sextic_certificate.json`.

The same file stores a 139-term Hermitian polynomial \(Z_6\) and 28
right-ideal multiplier terms.  The certificate contains rational
numbers only; no fitted decimal is used by the proof checker.

Let \({\cal F}_6\) be the sixth variable-Stein forcing obtained after
all three earlier factor columns, and put

\[
{\cal P}_6=C_6V^*.
\]

Exact partial-isometry word reduction gives

\[
\boxed{
{\cal F}_6+{\cal P}_6+{\cal P}_6^*
=Z_6-S^*Z_6S.}                                     \tag{4}
\]

At order six, the upper-complement Schur term is

\[
\begin{aligned}
{\cal S}_6={}&M_2RM_4F+M_4RM_2F+M_3RM_3F\\
&+M_2RM_2RM_2F,\qquad
R=I-F-\frac23E.
\end{aligned}
\]

Write

\[
G_j=ES^jF.
\]

The endpoint part of the certificate is the exact right-ideal
factorization

\[
\boxed{
FZ_6F+F{\cal S}_6+3G_3^*G_3
=A_1G_1+A_2G_2+(A_1G_1)^*+(A_2G_2)^*,}             \tag{5}
\]

where the two explicit polynomials \(A_1,A_2\) are reconstructed
from the 28 stored terms.

If \(x\in{\cal K}_2\), then

\[
G_jWx=ES^jWx=VB_j^*x=0,\qquad j=1,2.
\]

Physical congruence changes the sign and multiplies by four:

\[
[c^6]{\cal U}_{\rm final}
=-4W^*(FZ_6F+F{\cal S}_6)W.
\]

Compressing (5) therefore proves (2).

The exact checker leaves zero words in (4), zero words in (5), zero
words in \(Z_6-Z_6^*\), and zero words in \(E{\cal P}_6\).

## 3. The fifth-order Schur cost

L233 factors the fifth endpoint through \(G_1,G_2\).  Its two left
multipliers are

\[
\begin{aligned}
L_1&=F\left\{\frac72(S^*)^4S^3(S^*)^2
              +\frac52(S^*)^2S^3(S^*)^4\right\},\\
L_2&=F(S^*)^2S^4(S^*)^2.
\end{aligned}
\]

Two further exact word identities give, modulo terms ending in
\(G_1\),

\[
\begin{aligned}
EL_1^*F&=YG_3+\{\hbox{a polynomial}\}G_1,\\
EL_2^*F&=\{\hbox{a polynomial}\}G_1,
\end{aligned}                                      \tag{6}
\]

where

\[
Y=\frac72\{I-S^2(S^*)^2+S^*S^3(S^*)^2\}.
\]

The right-defect compression of \(Y\) is

\[
\boxed{
V^*YV=\frac72R_1,\qquad R_1=B_1^*B_1.}             \tag{7}
\]

Equations (6)--(7) and L233's physical factor \(-4\) imply that, for
\(x\in{\cal K}_2\),

\[
\boxed{
[c^5]{\cal U}_{\rm final}x
=-14B_1R_1B_3^*x.}                                 \tag{8}
\]

The \(L_2\) channel vanishes on this flag.

## 4. Completing the square

L232 gives the exact quartic face

\[
[c^4]{\cal U}_{\rm final}
=32B_1B_1^*+12B_2B_2^*
+56B_1R_1B_1^*.
\]

For \(x\in{\cal K}_2\), let

\[
z=B_3^*x,\qquad u=B_1^*y,\qquad v=B_2^*y,
\]

where \(y\) is an arbitrary first correction to the flag vector.
Using (2) and (8), the order-six quadratic form is bounded below by

\[
\begin{aligned}
&12\|z\|^2
+32\|u\|^2+56\langle R_1u,u\rangle
+12\|v\|^2
-28\operatorname{Re}\langle u,R_1z\rangle\\
&\quad\ge
\left\langle
\left[
12I-\frac{49}{2}
R_1(4I+7R_1)^{-1}R_1
\right]z,z
\right\rangle.                                     \tag{9}
\end{aligned}
\]

No pseudoinverse occurs in (9); \(4I+7R_1\) is uniformly invertible.
Since \(B_1\) is a contraction,

\[
0\preceq R_1\preceq I.
\]

The scalar function \(r^2/(4+7r)\) increases on \([0,1]\), so

\[
12I-\frac{49}{2}R_1(4I+7R_1)^{-1}R_1
\succeq
\left(12-\frac{49}{22}\right)I
=\frac{215}{22}I.                                  \tag{10}
\]

Equations (9)--(10) prove (3) uniformly through rank changes.

The leading \(Q\) in (1) also gives \(V^*C_6=0\).  L194's exact
metric chart therefore inserts this prepared row while retaining
exact lower tightness and positive Stein slack.

## 5. Rejected tempting representative

The first sparse sixth-column search imposed only the global
relations \(B_1=B_2=0\).  It passed 1,584 moderate numerical flag
tests and had the correct complete-delay trace, but its partial face
was

\[
B_3(12I-20B_1^*B_1)B_3^*.
\]

Valid rank-two Schur flags at parameter \(0.8\) made its effective
face negative, with a deterministic minimum \(-0.0352898\).  That
representative was rejected.  This is why (5), a true right-ideal
partial-flag identity, is essential; complete-delay reduction alone
is insufficient.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_canonical_sextic_preimage.py \
  --output \
  experiments/repeated_crabb_canonical_sextic_preimage_s70224.jsonl
```

The checker reconstructs every rational identity above, then audits:

1. six unstructured colligations;
2. sixteen rank-changing chains;
3. three complete double delays; and
4. three high-amplitude rank-two flags designed from the rejected
   candidate's failure mode.

All 28 records pass.  Every exact residual is zero.  The largest
direct sixth-face error is below \(9.7\times10^{-12}\), and the
smallest computed analytic middle weight is \(10.8166\), safely
above \(215/22\).

The tracked hashes are:

- certificate:
  `6e394d621637a569a8ab1ebac656e850911e1b99afb671c9f62e4aa9d302f425`;
- dataset:
  `4c8e763f56cf4422e2593b1eb93c73ba459bd4b8030397f1d612bf35b71b43a2`.

## 7. Next gate

L230--L234 now prepare the upper flag through transfer grade three.
The next useful step is not another isolated low-order expansion.
Return to L228's arbitrary-grade delayed anticommutator recursion and
use the explicit pattern learned here:

1. cancel each odd flag coefficient by a bounded right-ideal
   polynomial column;
2. make the next even direct face a positive transfer Gram; and
3. budget the preceding odd Schur cost by a uniform functional
   calculus bound such as (10).

An induction must control the growth of the polynomial-column norms
and preserve convergence in L194's analytic chart.
