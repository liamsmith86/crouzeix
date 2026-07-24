# Leading endpoint transfer and inner tangent (L168, 2026-07-23)

## 1. Result

Let `C=C_(L+1)` be the Crabb weighted shift,

\[
C_{j,j+1}=w_j,\qquad
w_0=w_{L-1}=\sqrt2,\quad w_j=1\ (1\le j\le L-2).
\]

Fix

\[
2\le k\le\lfloor L/2\rfloor,\qquad
m=L+2-k.
\]

Let `N_m` be L115's real Riesz representative of circular normal
mode `m`, and let `E` be its first inverse-Riemann pullback at `C`.
For the endpoint resolvent

\[
F(z)=e_0^*(zI-C)^{-1}e_L=2z^{-(L+1)},               \tag{1}
\]

the logarithmic derivative is

\[
\boxed{
D\log F[E]
=kz^{-m}+\beta_{L,k}z^m,
}                                                     \tag{2}
\]

where

\[
\beta_{L,2}=0,\qquad
\beta_{L,3}={L-9\over4L},                            \tag{3}
\]

and, for `k>=4`,

\[
\beta_{L,k}
=k-3-{2(k-2)^2\over L}.                             \tag{4}
\]

Let `gamma_(L,k)` be L163's exact optimized normal/defect-column
coefficient.  The symmetrized endpoint response of defect coordinate
`e_m` is

\[
\delta_{L,k}z^m,\qquad
\delta_{L,2}=1,\quad \delta_{L,k}=\sqrt2\ (k\ge3).
                                                               \tag{5}
\]

Then

\[
\boxed{\beta_{L,k}+k=\delta_{L,k}\gamma_{L,k}}       \tag{6}
\]

and hence

\[
\boxed{
D\log F[E]-\gamma_{L,k}D\log F[\hbox{sym }e_m]
=k(z^{-m}-z^m).
}                                                     \tag{7}
\]

The right side is an anti-self-reciprocal Laurent polynomial.  On
the unit circle it is purely imaginary, so its real mean is zero.
Thus the first normalized normal, after its exact optimized defect
response is removed, is already a logarithmic-inner tangent.

This is the first explicit scalar bridge between L163's
normal/defect coefficient and L149's inner-function mechanism.

## 2. First inverse-Riemann coefficient

The entries of `N_m` lie on two diagonals:

\[
j-i=m+1\quad\hbox{and}\quad i-j=m-1,                 \tag{8}
\]

with entry `r_i r_j`, where the endpoint weights are
`r_0=r_L=1/sqrt(2)` and all interior `r_j=1`.

The first inverse-Riemann recurrence gives

\[
E=N_m-\eta_{L,k}C^{m+1}.                             \tag{9}
\]

For `k<=2`, the positive diagonal in (8) is absent and
`C^(m+1)=0`, so `eta_(L,k)=0`.  For `k>=3`, the support-vector
pairing has the following two contributions:

\[
A_-={k-1\over L},\qquad
A_+=
\begin{cases}
1/(4L),&k=3,\\
(k-3)/L,&k\ge4.
\end{cases}                                         \tag{10}
\]

The first term counts the `k` negative-diagonal entries with the two
endpoint half-weights.  The second counts the `k-2`
positive-diagonal entries; at `k=3` the sole entry joins both
endpoints and has weight `1/2`.  Substitution in the first Schwarz
coefficient gives

\[
\eta_{L,k}=A_-+A_+
=
\begin{cases}
0,&k\le2,\\
9/(4L),&k=3,\\
2(k-2)/L,&k\ge4.
\end{cases}                                         \tag{11}
\]

This is a direct specialization of the inverse-Riemann recurrence,
not a fitted coefficient.

## 3. Endpoint path count

Write

\[
R(z)=(zI-C)^{-1}=\sum_{j=0}^L C^jz^{-j-1}.
\]

Then

\[
DF[E]=e_0^*RER e_L.                                 \tag{12}
\]

Every nonzero summand in (12) is a directed weighted path from the
right endpoint to the left endpoint with one insertion of `E`.

The negative diagonal in (8) gives exactly `k` paths.  After division
by (1), every endpoint weight cancels, and their total is

\[
kz^{-m}.                                            \tag{13}
\]

The positive diagonal gives

\[
\nu_kz^m,\qquad
\nu_3={1\over4},\quad \nu_k=k-3\ (k\ge4).            \tag{14}
\]

For `k=3`, the sole inserted entry carries both endpoint half-weights.
For `k>=4`, the two boundary paths contribute one half each and the
`k-4` interior paths contribute one each.

Finally, an insertion of `C^(m+1)` can occur in `k-2` positions in
the endpoint path, so the Schwarz term in (9) contributes

\[
-\eta_{L,k}(k-2)z^m.                                \tag{15}
\]

Equations (11), (14), and (15) give

\[
\beta_{L,k}=\nu_k-\eta_{L,k}(k-2),
\]

which is exactly (3)--(4), proving (2).

## 4. Characteristic and defect checks

The same path geometry recovers L167:

\[
D\det(zI-C)[E]=-kz^{k-1}.                            \tag{16}
\]

Only the `k` entries on the negative diagonal close a trace with a
power of `C`; the polynomial correction in (9) has zero trace.

For the one-sided defect transfer,

\[
{e_m^*R(z)e_L\over e_0^*R(z)e_L}
=
\begin{cases}
\frac12z^m,&k=2,\\
\frac1{\sqrt2}z^m,&k\ge3.
\end{cases}                                         \tag{17}
\]

Reciprocal reversal (L164) supplies the matching opposite endpoint
leg, doubling (17) and proving (5).

For `k=2`, L163 gives `gamma_(L,2)=2`, so (6) is `2=2`.
For `k=3`,

\[
\beta_{L,3}+3={13L-9\over4L}
=\sqrt2\,{\sqrt2(13L-9)\over8L}.
\]

For `k>=4`,

\[
\begin{aligned}
\beta_{L,k}+k
&=2k-3-{2(k-2)^2\over L}\\
&=\sqrt2\left[
\sqrt2\left({2k-3\over2}-{(k-2)^2\over L}\right)
\right].
\end{aligned}
\]

These are precisely (6) with L163's three cases for `gamma`.

## 5. Scope and remaining gate

Equation (7) proves the desired inner-tangent cancellation only for
the **first strong coefficient** at the Crabb point.  It does not yet
prove L163.  The physical reflected path produces mixed
inverse-Riemann, endpoint, characteristic, and optimized-defect
coefficients at every weight from `d+1` through `2d`.

The next theorem must lift (7) covariantly through that complete
marked expansion.  A viable route is:

1. express the total normal variation through the moving companion
   transfer `q^*(zI-T)^(-1)e_L`;
2. use L164 to pair its two endpoint/defect legs;
3. use L162 for the zero-reflection equality transport;
4. prove the remaining one-reflection logarithmic derivative stays
   anti-self-reciprocal, so L149's real-mean-zero argument applies.

Grade one is correctly excluded: `m=L+1` is not a defect coordinate,
and its leading characteristic derivative is the nonzero constant
`-1`.  There is no matched positive mode, leaving the observed
nonzero condition derivative.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_leading_endpoint_transfer.py \
  --minimum-length 4 --maximum-length 14 \
  --output \
  experiments/crabb_leading_endpoint_transfer_s70223.jsonl
```

The checker constructs `C`, `N_m`, `E`, and the finite nilpotent
resolvent directly.  It independently verifies (2), (6), (7), and
(16) for every admissible grade through length fourteen, including
the grade-one discriminator.
