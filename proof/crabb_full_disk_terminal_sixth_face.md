# The terminal full-disk sixth face (L178, 2026-07-24)

## 1. Result and scope

Put `p=L+1`, assume `L>=6`, and index the `L-1` nonconstant
Toeplitz coefficients from zero.  Let

\[
 k=L-3,\qquad z_k=a,\qquad z_{k+1}=b,
\]

with every earlier coefficient zero.  Along L176's recentered exact
disk path, the complete sixth-order true-normal Schur face is strictly
negative unless `a=0` or `b=0`.

More precisely, its completed-square gain divided by the canonical
sixth-order base deficit is

\[
\boxed{
 R_k=
 {6(4k-1)^2\over
  6(4k-1)^2+169k(k-1)(k-2)}
 <1.}                                                \tag{1}
\]

Thus the edge which makes A126's stronger real `p=7` inequality sharp
is safe for the actual Schur problem in every size.  In fact `R_k`
decreases to zero.

This proves only the two-terminal-coefficient face.  It does not prove
A126's full sixth-order inequality, classify its complete kernel, or
close the nonlinear full circular-range tube.

## 2. The sparse L176 correction

Write

\[
 A=|a|^2,\qquad B=|b|^2,\qquad
 D=(k-1)A+(k+1)B.
\]

Only three intrinsic Pluecker coordinates enter L176's correction.
Substitution in its interval pulses gives the following tridiagonal
Hermitian matrix.  Its diagonal is

\[
\begin{array}{c|c}
j&B(W)_{jj}\\ \hline
0,L-1&2D/L\\
1,L-2&2(k-1)A/L-4B/L\\
2,\ldots,L-3&-8A/L-4B/L ,
\end{array}                                         \tag{2}
\]

and its upper first diagonal is

\[
 B(W)_{j,j+1}=
 \begin{cases}
  2k\overline a b/(k+2),&j=0,L-2,\\
  -4\overline a b/(k+2),&1\le j\le L-3.
 \end{cases}                                        \tag{3}
\]

Equations (2)--(3) also provide a short independent check on all
index and conjugation conventions.

## 3. Exact endpoint recurrence

Use

\[
 H(s)=I/2+sZ(z)+s^2B(W)
\]

in L122's coefficient-gauge disk model.  Expand `K^{-1}`, the
rank-one Stein equation, and the two simple generalized endpoints as
in L122 equations (15)--(17).  Sparse path selection gives

\[
\begin{array}{c|ccc}
&[s^0]&[s^2]&[s^6]\\ \hline
\lambda_-&1/2&2D/L&0\\
\lambda_+&2&8D/L&
 -128k^2A^2B/(k+2)^2 .
\end{array}                                         \tag{4}
\]

All endpoint coefficients of degrees one, three, four, and five
vanish.  The lower formula is also immediate from
`M e_0=H_{00}K e_0`.  For the upper row, (2)--(3) leave only the two
terminal long jumps, the first-diagonal correction, and their
reflections; collecting the length-two and length-three return paths
gives the last entry of (4).

Consequently, for
\(\delta=\lambda_+-4\lambda_-\), the condition-square deficit is

\[
\boxed{
 D_{6,L}(a,b)
 =-2[s^6]\delta
 ={256k^2\over(k+2)^2}|a|^4|b|^2.}                  \tag{5}
\]

The exact series engine verifies (5) symbolically with unrestricted
complex `a,b` through `L=8`, symbolically on the real slice through
`L=10`, and on a fixed Gaussian-rational complex ray through `L=15`.

## 4. Highest-mode cubic response

For a general direction, put

\[
 S_t=\sum_{0\le i<t-i}(t-2i)W_{i,t-i}.
\]

The order-three characteristic/reversed-Horner and inverse-Riemann
recurrences give, for the highest possible active mode `k=L-3`,

\[
 G_{L,k}^{(3)}
 =-{32(4k-1)\over L^2(k+1)(k+2)}
 \left\{{k+1\over2}z_kS_k
       -{k+2\over2}z_{k+1}S_{k-1}\right\}.           \tag{6}
\]

On the terminal direction,

\[
 S_k=-ka\overline b,\qquad S_{k-1}=0,
\]

and every lower mode vanishes.  Therefore

\[
\boxed{
 G_{L,k}^{(3)}
 ={16k(4k-1)\over L^2(k+2)}a^2\overline b.}          \tag{7}
\]

The full exact characteristic/Riemann engine checks (6) on generic
complex directions through `L=10` and checks the terminal
specialization (7), including the absence of other modes, through
`L=15`.

## 5. Exact Schur margin

Here the L173 remainder is `r=L-k=3`.  Its exact curvature (L173
equation (13)) becomes

\[
 b_{L,k}
 ={6(4k-1)^2+169k(k-1)(k-2)\over24L^4}.             \tag{8}
\]

The two real polarizations combine to the gain

\[
 {\cal G}_{L,k}={|G_{L,k}^{(3)}|^2\over4b_{L,k}}.
\]

Dividing this expression by (5) gives exactly (1).  The unused
fraction is

\[
\boxed{
 1-R_k=
 {169k(k-1)(k-2)\over
  6(4k-1)^2+169k(k-1)(k-2)}>0.}                     \tag{9}
\]

This is the same positive rank-one null-lift term that made L173's
apex Schur face strict.  The first values are

\[
 R_3={121\over290},\quad
 R_4={225\over901},\quad
 R_5={361\over2051},\quad
 R_6={529\over3909}.
\]

Thus A126's `P_6=9|C|^2` equality at `k=3` is not equality of the
actual Schur face: L65's curvature retains the strict `169` margin.

## 6. Consequence for the general attack

The terminal edge does not obstruct the sixth-order theorem.  It
instead identifies the correct proof mechanism:

1. a flux part matched to the cubic response;
2. a positive null-lift part with the characteristic factor `169`;
3. a base deficit which must dominate the flux block before the
   null-lift makes the Schur complement strict.

The next useful calculation is therefore the same decomposition for
all active anti-diagonals of
`z tensor (z wedge J conjugate(z))`, not a larger monomial SOS for
the isolated `p=7` polynomial.

## 7. Exact regeneration

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_terminal_sixth_exact.py \
  --minimum-length 6 --maximum-length 15 \
  --output \
  experiments/crabb_full_disk_terminal_sixth_exact_s70224.jsonl
```

The checker compares the endpoint recurrence, cubic
characteristic/Riemann response, L173 curvature completion, and the
closed ratio independently in exact SymPy arithmetic.
