# The exact real `p=7` sixth-face certificate (L179, 2026-07-24)

## 1. Result and scope

On A126's real `L=6` slice, let

\[
 z=(a_1,a_2,a_3,a_4,a_5)
\]

and define the exact base sextic by

\[
 [s^6]\delta=-{32\over225}P_6(z).
\]

The complex cubic response becomes the real cubic

\[
\begin{aligned}
C_6={}&6a_1a_2a_4-5a_1a_3a_5+2a_2a_3a_4\\
     &-2a_3a_4^2+5a_3a_5^2-6a_4^2a_5.
\end{aligned}                                       \tag{1}
\]

The exact inequality required by A126's true-normal Schur complement
is

\[
\boxed{
 P_6(z)-{1089\over290}C_6(z)^2\ge0
 \qquad(z\in\mathbb R^5).
}                                                    \tag{2}
\]

This note proves (2) by an exact rational rank-seven Gram certificate.
There is no floating-point or solver premise in the proof.

Equation (2) closes the complete **real** sixth-order face in the first
active dimension `p=7`.  It does not prove the complex `p=7` face, the
arbitrary-size block inequality, the higher-order kernel lift, or the
nonlinear tube.  The stronger sharp-looking inequality
`P_6>=9C_6^2` remains numerical.

## 2. Equality-ideal coordinates

The real phase-palindromic equality set is the union of

\[
\begin{aligned}
{\cal E}_+&=\{a_5=a_1,\ a_4=a_2\},\\
{\cal E}_-&=\{a_5=-a_1,\ a_4=-a_2,\ a_3=0\}.
\end{aligned}                                       \tag{3}
\]

Let

\[
{\cal I}_3=
\left(I({\cal E}_+)\cap I({\cal E}_-)\right)_3 .
\]

This homogeneous cubic space has dimension 21.  The checker constructs
an explicit binomial basis

\[
 q=(q_0,\ldots,q_{20})^T
\]

and verifies by direct substitution that every `q_j` vanishes on both
branches in (3).  For example,

\[
\begin{aligned}
q_0&=-a_1^2a_4+a_1a_2a_5,\\
q_1&=-a_1a_2a_3+a_1a_3a_4,\\
q_2&=-a_1^2a_3+a_1a_3a_5,
\end{aligned}
\]

with all 21 cubics listed literally in the regenerating script.  This
basis performs the necessary facial reduction before the Gram step;
a raw 35-monomial SDP is singular and obscures the exact certificate.

## 3. Sparse rank-seven certificate

Let the pivot set be

\[
 I=(0,1,2,3,7,8,10).
\]

Define the `21 x 7` rational matrix `F` by the following nonzero rows;
unlisted rows are zero:

\[
\begin{array}{c|rrrrrrr}
0&100&0&0&0&0&0&0\\
1&0&72&0&0&0&1089/29&0\\
2&0&0&20655/58&0&0&0&4563/29\\
3&0&0&0&100&0&0&-100\\
4&100&216&0&0&0&3267/29&0\\
6&0&0&0&100&0&0&-100\\
7&0&0&0&0&8262/145&0&108\\
8&0&1089/29&0&0&0&450&0\\
10&0&0&4563/29&-100&108&0&424\\
11&-100&0&0&0&0&0&0\\
12&60&0&0&0&0&0&0\\
13&0&0&0&150&0&0&-150\\
14&0&0&1521/29&0&-3042/145&0&0\\
15&180&0&0&0&0&0&0\\
16&0&0&-7605/58&0&1521/29&0&0\\
18&0&0&4563/29&0&-9126/145&0&0 .
\end{array}                                         \tag{4}
\]

Its pivot core `S=F[I,:]` is the symmetric matrix

\[
S=\begin{bmatrix}
100&0&0&0&0&0&0\\
0&72&0&0&0&1089/29&0\\
0&0&20655/58&0&0&0&4563/29\\
0&0&0&100&0&0&-100\\
0&0&0&0&8262/145&0&108\\
0&1089/29&0&0&0&450&0\\
0&0&4563/29&-100&108&0&424
\end{bmatrix}.                                      \tag{5}
\]

Put

\[
\boxed{Q=FS^{-1}F^T.}                               \tag{6}
\]

The leading principal minors of `S` are

\[
\begin{gathered}
100,\quad7200,\quad{74358000\over29},\quad
{7435800000\over29},\\
{12286915920000\over841},\quad
{4447604001941190000\over707281},\\
{6419933240166733716000\over20511149}.
\end{gathered}                                       \tag{7}
\]

They are strictly positive, so Sylvester's criterion gives `S>0`.
Equation (6) therefore gives

\[
Q\succeq0,\qquad \operatorname{rank}Q=7.             \tag{8}
\]

The pivot identity `Q[:,I]=F` is checked exactly and provides a compact
independent audit of the low-rank reconstruction.

## 4. Polynomial identity

Exact rational expansion gives

\[
\boxed{
 q(z)^TQq(z)
 =P_6(z)-{1089\over290}C_6(z)^2.
}                                                    \tag{9}
\]

Equations (8)--(9) prove (2).  The regeneration builds `P_6`
independently from the recentered endpoint recurrence, rather than
copying its 68 coefficients into the certificate.

This also explains why the earlier floating Gram programs were
ill-conditioned: the valid Gram matrix has rank seven inside a
21-dimensional equality-ideal basis, itself sitting inside the
35-dimensional cubic monomial space.

## 5. Consequence for the general attack

Three formerly uncertain pieces of A126 are now exact:

1. L178 closes the sharp last-two-coefficient edge in every size.
2. L179 closes every real direction in the first active size.
3. L180 closes every complex direction in the first active size by a
   separate ten-square tensor certificate.

The distinction remains important: this real SOS does not extend by
formal replacement of squares with moduli, because the cubic response
contains conjugate Pluecker sectors.  L180 instead uses separate
global-phase blocks on

\[
z\otimes(z\wedge J\overline z),
\]

and verifies them independently.  The next algebraic target is the
arbitrary-size anti-diagonal/interval-flux block theorem.

## 6. Exact regeneration

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_real_sixth_certificate.py \
  --output \
  experiments/crabb_full_disk_real_sixth_certificate_s70224.jsonl
```

The script verifies the equality-branch ideal, the seven positive
principal minors, the exact rank/reconstruction, and every coefficient
of (9).  A clean second run is required to reproduce the JSON
byte-for-byte.
