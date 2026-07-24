# The full-disk Plücker correction isometry (L176, 2026-07-24)

## 1. Result

Put \(p=L+1\), \(n=L-1\), and index the nonconstant Toeplitz
coordinate as

\[
 z=(z_0,\ldots,z_{n-1}).
\]

Let

\[
 W=z\wedge J\overline z,\qquad
 W_{ij}=z_i\overline{z_{n-1-j}}
        -z_j\overline{z_{n-1-i}}.                    \tag{1}
\]

This note constructs an explicit traceless Hermitian matrix \(B(W)\)
such that the tangent of L122's full general-H disk chart,

\[
 E(W)=D_{I/2}X[B(W)],
\]

obeys the exact L65 curvature identity

\[
\boxed{\langle E(W),{\cal C}_pE(W)\rangle
       =32\{\|z\|^4-|z^TJz|^2\}
       =32\|W\|_{\wedge^2}^2.}                       \tag{2}
\]

Here \({\cal C}_p=-e_p\) is L65's positive curvature.  Thus L175's
conjectured full-disk recentering has the exact required energy in
every size.

Equation (2) does **not** yet prove L175's ambient-response identity
\(g_2=2{\cal C}_pE(W)\).  That matching is the remaining half of the
tight full-H face.

## 2. Explicit correction

The relation \(J\overline WJ=-W\) lets us use the intrinsic entries

\[
 0\le i<j\le n-1,\qquad i+j\le n-1.                  \tag{3}
\]

For \(t=1,\ldots,n-1\), put \(r=n-1-t\).  For every
\(0\le i<t-i\), define the length-\(t+2\) interval pulse

\[
 b_i^{(t)}(a)
 =2\,{\bf1}_{\{i+1\le a\le t-i\}}
  -{2(t-2i)\over t+2},
 \qquad 0\le a\le t+1.                               \tag{4}
\]

The second term removes the mean.  Define the upper offset-\(r\)
entries of \(B(W)\) by

\[
\boxed{
 B(W)_{a,a+r}
 =\sum_{0\le i<t-i}
   b_i^{(t)}(a)\,\overline{W_{i,t-i}},
 \qquad 0\le a\le t+1.}                              \tag{5}
\]

Complete (5) Hermitianly.  When \(r=0\), the entries in (5) are real
by (1).  The mean removal makes the diagonal trace zero and fixes the
irrelevant scale gauge.

For example, the first three pulse matrices are

\[
\begin{aligned}
t=1:&\quad {1\over3}(-2,4,-2),\\
t=2:&\quad (-1,1,1,-1),\\
t=3:&\quad {1\over5}
\begin{bmatrix}
-6&-2\\4&-2\\4&8\\4&-2\\-6&-2
\end{bmatrix}.                                      \tag{6}
\end{aligned}
\]

These are the exact rational patterns first exposed by L175's
full-H stationarity solve; no fitted coefficient is used below.

## 3. Differentiated disk chart

Fix \(t<n-1\), so \(r=n-1-t\ge1\), and set \(m=t+2\).
An offset-\(r\) Hermitian perturbation \(B\) contributes only to
L65's paired circle mode \(k=r\).

Differentiate

\[
X(H)=2K^{-1/2}\widehat HRK^{-1/2},\qquad
K=\widehat H+R^*\widehat HR
\]

at \(H=I/2\), including the two endpoint divided differences of
\(K^{-1/2}\).  Substitute the result in L65 equations (12)--(14).
For the offset vector \(b=(B_{a,a+r})_{a=0}^{m-1}\), the reduced
vector is

\[
{\cal D}_mb,
\]

where

\[
\begin{aligned}
({\cal D}_mb)_0&=\sqrt2(b_0-b_1),\\
({\cal D}_mb)_a&=-b_{a-1}+2b_a-b_{a+1},
                    &&1\le a\le m-2,\\
({\cal D}_mb)_{m-1}&=\sqrt2(b_{m-1}-b_{m-2}).         \tag{7}
\end{aligned}
\]

This is a direct entrywise differentiation; in particular the
endpoint square-root gauge is present.

Let \(\chi_i\) be the indicator of

\[
 I_i=\{i+1,\ldots,m-2-i\}.
\]

Since \({\cal D}_m{\bf1}=0\), equation (4) gives the flux vector

\[
 s_i={\cal D}_mb_i^{(t)}
     ={\cal D}_m(2\chi_i).                            \tag{8}
\]

If \(R_m\) is L65 equation (15)'s singular path matrix, direct
endpoint substitution gives the useful equivalent formula

\[
\boxed{s_i=16R_m\chi_i.}                              \tag{9}
\]

## 4. Exact paired-mode isometry

Let

\[
q_0=(1/\sqrt2,1,\ldots,1,1/\sqrt2)^T
\]

span \(\ker R_m\), and let \(q_{m,k}\) be L65's normalizing vector.
Put

\[
x_i=16\left(
\chi_i-
{q_{m,k}^T\chi_i\over q_{m,k}^Tq_0}q_0
\right).                                             \tag{10}
\]

Then \(q_{m,k}^Tx_i=0\) and, by (9),

\[
R_mx_i=s_i.                                          \tag{11}
\]

For \(k\ge3\), the rank-one term in L65 equation (16) therefore
annihilates \(x_i\), so

\[
{\cal K}_{m,k}s_i=x_i.
\]

For \(k=1,2\), the same equation is exactly the constrained inverse
selected by L65's \(\gamma\to\infty\) limit.

The fluxes satisfy

\[
s_j^Tq_0=0,\qquad s_j^T\chi_i=4\delta_{ij}.            \tag{12}
\]

Equations (10)--(12) now give

\[
\boxed{s_j^T{\cal K}_{m,k}s_i=64\delta_{ij}.}          \tag{13}
\]

Thus every nonterminal intrinsic complex Plücker coordinate contributes
\(64|W_{i,t-i}|^2\), with no cross terms.

## 5. Terminal diagonal mode

When \(t=n-1\), equation (5) is a real diagonal perturbation
\(d=(d_0,\ldots,d_{L-1})\).  The exact disk chart is a weighted shift.
Near the Crabb weights, its optimal similarity square is the square
of the product of its weights: the diagonal contraction metric gives
the upper bound, and the polynomial \(z^L\) gives the matching lower
bound.

Expanding that product for \(H=I/2+\tau\operatorname{diag}d\) gives
the L65 curvature

\[
\boxed{
\langle E(d),{\cal C}_pE(d)\rangle
=4\sum_{a=0}^{L-2}(d_{a+1}-d_a)^2.}                  \tag{14}
\]

The pulse columns (4) have first differences supported on distinct
pairs of edges, each with squared norm eight.  Hence their Gram
matrix under (14) is exactly

\[
32I.                                                  \tag{15}
\]

## 6. Summation

For \(t<n-1\), every intrinsic entry in (3) has a reflected partner
of equal magnitude in the full exterior square.  Thus its contribution
to \(\|W\|_{\wedge^2}^2\) is \(2|W_{i,t-i}|^2\), matching
\(64|W_{i,t-i}|^2\) in (13).

On the terminal anti-diagonal, the intrinsic entries are real and
unpaired.  Equations (14)--(15) give
\(32W_{i,n-1-i}^2\).  Summing all anti-diagonals proves (2).

## 7. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_full_disk_correction_isometry.py \
  --minimum-length 3 --maximum-length 30 \
  --output experiments/crabb_full_disk_correction_isometry_s70224.jsonl
```

The checker uses exact SymPy arithmetic.  It verifies (7)--(13) and
(15) in every length through 30; a second run is byte-identical.
