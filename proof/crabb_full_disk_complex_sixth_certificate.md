# The exact complex `p=7` sixth-face certificate (L180, 2026-07-24)

## 1. Result and scope

On A126's complex `L=6` slice, write

\[
 [s^6]\delta=-{32\over225}P_6(z).
\]

With \(W=z\wedge J\overline z\), put

\[
 C_6=6z_3W_{03}-5z_4W_{02}+2z_3W_{12}.
\]

The complete sixth-order true-normal Schur face is nonpositive exactly
when

\[
\boxed{
 P_6(z)-{1089\over290}|C_6(z)|^2\ge0.
}                                                    \tag{1}
\]

This note proves (1) for every \(z\in\mathbb C^5\) by an exact rational
ten-square certificate.  It closes the first active complex size
`p=7`.  It does not prove the stronger conjectural constant `9`, the
arbitrary-size sixth-order block inequality, the eighth-order kernel
lift, or the nonlinear full circular-range tube.

An independent upstream checker propagates five generic complex
variables through the characteristic/Riemann recurrence.  It verifies
that all quadratic normal responses vanish, mode three is the only
cubic normal mode among all ten polarizations, and its response is
exactly \(-22C_6/45\).  Thus (1) accounts for the complete complex
sixth-order Schur face, not merely a proposed mode-three subface.

## 2. Pluecker-tensor coordinates

For \(0\le i<j\le4\), define

\[
 W_{ij}=z_i\overline z_{4-j}-z_j\overline z_{4-i}.
\]

Select the following twenty coordinates

\[
 q_r=-z_aW_{ij},
\]

where the triples \((a,i,j)\), in order \(r=0,\ldots,19\), are

```text
(0,0,1) (0,0,2) (0,0,3) (0,1,2) (0,1,3)
(0,1,4) (0,2,4) (1,0,1) (1,0,2) (1,0,3)
(1,0,4) (1,1,2) (1,1,4) (1,2,3) (2,0,2)
(2,0,3) (2,0,4) (2,1,2) (2,1,3) (3,0,3).
```

Write \(x_r=\operatorname{Re}q_r\) and
\(y_r=\operatorname{Im}q_r\).  These coordinates are a sparse facial
reduction of the forty-dimensional span of
\(z\otimes(z\wedge J\overline z)\); no assertion that the twenty
complex coordinates are independent is needed.

## 3. Ten-square identity

Define

\[
\begin{aligned}
R_1={}&10x_0+10x_5-10x_{10}+18x_{15}+6x_{17},\\
R_2={}&x_7+x_4-x_9+\tfrac32x_{14},\\
R_3={}&x_1+\tfrac65x_9+\tfrac25x_{11},\\
R_4={}&3x_2+x_3,\\
R_5={}&x_8,\\
R_6={}&x_6-\tfrac25x_{13}-x_{16}
       +\tfrac25x_{18}+\tfrac65x_{19},
\end{aligned}
\]

and

\[
\begin{aligned}
I_1={}&10y_0-10y_5+10y_{10}-18y_{15}-6y_{17},\\
I_2={}&y_7-y_4+y_9-\tfrac32y_{14},\\
I_3={}&y_1-\tfrac65y_9-\tfrac25y_{11},\\
I_4={}&y_6-\tfrac25y_{13}-y_{16}
       +\tfrac25y_{18}+\tfrac65y_{19}.
\end{aligned}
\]

Exact rational expansion gives

\[
\boxed{
\begin{aligned}
P_6-{1089\over290}|C_6|^2
={}&R_1^2+100R_2^2+225R_3^2+72R_4^2+450R_5^2\\
 &+{7605\over58}R_6^2\\
 &+I_1^2+100I_2^2+225I_3^2+{7605\over58}I_4^2.
\end{aligned}
}                                                    \tag{2}
\]

Every coefficient in (2) is positive.  Hence (2) proves (1) without a
floating-point or semidefinite-solver premise.  The real and imaginary
Gram matrices have ranks six and four respectively.

## 4. How the certificate was reduced

A circle-grade invariant tensor Gram model initially has 230
coefficients but only coefficient rank 91.  Numerical facial reduction
located an exact sixty-coefficient support.  Solving its coefficient
equations over the rationals leaves eighteen affine parameters.

Positive-semidefinite completion splits into four independent grade
blocks.  Its facial equations reduce the eighteen parameters to a
one-parameter rational family \(0\le\lambda\le225\).  The endpoint
\(\lambda=0\) is (2), has Gram rank ten, and is substantially cleaner
than rationalizing a raw numerical Gram matrix.

This discovery history is not part of the proof.  The certificate
checker constructs \(P_6\) independently from the endpoint recurrence,
builds the ten displayed rank-one factors directly, and checks all 137
polarized sextic coefficients.

## 5. Consequence and next target

L179's real certificate was not being mistaken for a complex one:
the naive Hermitian lift of its rank-seven Gram matrix is false.  L180
instead retains the missing global-phase sectors through separate real
and imaginary tensor blocks.

The first-active-size obstruction is now closed over its full complex
domain.  The live algebraic task remains A126's arbitrary-size block:
derive \(D_{6,L}\) and every \(G_{L,k}^{(3)}\) in common
anti-diagonal/interval-flux coordinates, then prove the resulting
Schur blocks positive while retaining L173's null-lift curvature.
Only after that should the sixth-face kernel and its eighth-order
fallback be classified.

## 6. Exact regeneration

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_cubic_response.py \
  --output \
  experiments/crabb_full_disk_complex_cubic_response_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_sixth_certificate.py \
  --output \
  experiments/crabb_full_disk_complex_sixth_certificate_s70224.jsonl
```

The first script verifies the complete generic complex response.  The
second verifies conjugation invariance, positive rank-one factors, the
two Gram ranks, and every coefficient of (2).  Clean second runs must
reproduce both JSON files byte-for-byte.
