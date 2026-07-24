# The exact complex `p=8` sixth-face certificate (L181, 2026-07-24)

## 1. Result and scope

Let `L=7` in A126's recentered full-disk path and write
\(\delta_6=[s^6]\delta\).  The complete cubic true-normal response has
two complex modes:

\[
\begin{aligned}
G_3={176\over147}\big(
 &3z_5W_{02}-4z_3W_{04}-2z_3W_{13}
 -2\overline z_0\,\overline W_{01}\big),\\
G_4={}&G^{(3)}_{7,4},
\end{aligned}                                      \tag{1}
\]

where \(G_4\) is L178's highest-mode anti-diagonal formula.  L173 gives

\[
b_{7,3}={773\over9604},\qquad
b_{7,4}={901\over9604}.                            \tag{2}
\]

The required sixth-order Schur residual is

\[
\boxed{
{\cal R}_7(z)
=-2\delta_6(z)
-{|G_3(z)|^2\over4b_{7,3}}
-{|G_4(z)|^2\over4b_{7,4}}.
}                                                   \tag{3}
\]

This note proves

\[
\boxed{{\cal R}_7(z)\ge0\qquad(z\in\mathbb C^6)}
                                                            \tag{4}
\]

by an exact rational fifteen-square identity.  It closes the complete
complex sixth-order face in the second active size `p=8`.

This is not the general Crouzeix conjecture, a theorem for arbitrary
`8 x 8` matrices, the arbitrary-size sixth block, the higher-order
kernel lift, or the nonlinear full circular-range tube.

## 2. Exact upstream response audit

A direct characteristic expansion in six generic complex variables
was abandoned after more than an hour and 9 GB of symbolic
intermediates.  It is unnecessary.

Every degree-two response is a homogeneous quadratic, and every
degree-three response is a homogeneous cubic, in the twelve real and
imaginary Toeplitz coordinates.  For \(N=12\), use the sparse set

\[
\begin{aligned}
{\cal U}={}&\{e_i\}\\
&{}\cup\{e_i+e_j,e_i-e_j:i<j\}\\
&{}\cup\{e_i+e_j+e_k:i<j<k\}.
\end{aligned}                                      \tag{5}
\]

It has

\[
12+2{12\choose2}+{12\choose3}=364
={14\choose3}
\]

points.  The exact degree-three evaluation matrix on \({\cal U}\) has
rank `364`; its degree-two submatrix has rank `78={13 choose 2}`.
Equivalently, the familiar polarization argument first kills pure
cubes, then both two-variable cubic coefficients, then every
three-variable coefficient.

The checker evaluates the full characteristic/Riemann recurrence on
all 364 points and all twelve normal polarizations.  All 4,368
polarization evaluations are exact.  They prove:

1. every quadratic response vanishes;
2. modes three and four are the only cubic modes;
3. both formulas in (1) hold coefficientwise over the full complex
   slice.

Thus (3) contains the complete sixth-order Schur gain.

## 3. Pluecker-tensor coordinates

Put

\[
W_{ij}=z_i\overline z_{5-j}-z_j\overline z_{5-i},
\qquad q_{aij}=-z_aW_{ij},
\]

and write \(x_{aij}=\operatorname{Re}q_{aij}\),
\(y_{aij}=\operatorname{Im}q_{aij}\).

Use five ordered coordinate groups:

```text
Q0:
(0,0,4) (0,1,3) (1,0,3) (1,1,2) (2,0,2)

Q1:
(0,0,3) (0,1,2) (0,1,4) (1,0,2)
(1,0,4) (1,1,3) (2,0,3) (2,1,2)

Q2:
(0,0,2) (0,1,5) (0,2,4) (1,0,1) (1,0,5)
(1,2,3) (2,0,4) (2,1,3) (3,0,3)

Q3:
(0,0,1) (0,2,5) (0,3,4) (1,1,5)
(1,2,4) (2,0,5) (3,0,4) (3,1,3)

Q4:
(0,3,5) (1,2,5) (1,3,4) (2,1,5)
(3,0,5) (3,1,4) (4,0,4)
```

For a coefficient vector \(v\), let \(v\cdot x(Q_j)\) mean the linear
combination in the displayed order, and similarly for \(y\).

## 4. Fifteen-square identity

The nine real-coordinate factors \((\lambda,j,v)\) are

\[
\begin{array}{c|c|l}
\lambda&j&v\\ \hline
2048/9&0&(1,\frac12,0,0,0)\\
512/9&0&(0,1,-1,0,\frac32)\\
4608/25&0&(0,0,1,\frac13,0)\\
2304/25&1&(1,\frac13,0,0,\frac{10}9,\frac59,0,0)\\
256/9&1&(0,0,1,\frac32,-1,0,\frac95,\frac35)\\
64&2&(1,\frac23,0,0,-\frac23,0,\frac43,\frac23,0)\\
64&2&(0,0,1,\frac23,0,-\frac25,-1,\frac25,\frac65)\\
73984/6957&3&(1,\frac32,0,0,0,-\frac32,2,1)\\
1557504/22525&4&(1,\frac13,-\frac59,-\frac13,-1,\frac59,\frac{10}9).
\end{array}                                        \tag{6}
\]

The six imaginary-coordinate factors are

\[
\begin{array}{c|c|l}
\lambda&j&v\\ \hline
2304/25&1&(1,\frac13,0,0,-\frac{10}9,-\frac59,0,0)\\
256/9&1&(0,0,1,-\frac32,-1,0,\frac95,\frac35)\\
64&2&(1,-\frac23,0,0,\frac23,0,-\frac43,-\frac23,0)\\
64&2&(0,0,1,-\frac23,0,-\frac25,-1,\frac25,\frac65)\\
73984/6957&3&(1,-\frac32,0,0,0,\frac32,-2,-1)\\
1557504/22525&4&(1,\frac13,-\frac59,-\frac13,-1,\frac59,\frac{10}9).
\end{array}                                        \tag{7}
\]

Exact rational expansion gives

\[
\boxed{
{\cal R}_7
=\sum_{(\lambda,j,v)\in(6)}
 \lambda\,[v\cdot x(Q_j)]^2
+\sum_{(\lambda,j,v)\in(7)}
 \lambda\,[v\cdot y(Q_j)]^2.
}                                                   \tag{8}
\]

Every scale in (6)--(7) is positive.  Hence (8) proves (4).  The real
and imaginary Gram ranks are nine and six.  The checker reconstructs
\(\delta_6\) independently from the endpoint recurrence and matches
all 299 polarized sextic coefficients; there is no floating solver or
fitted coefficient premise in the proof.

## 5. Structural consequence

L180 and L181 now close the first two active complex sizes.  Their
factorizations share the same organization:

1. tensor coordinates \(z\otimes(z\wedge J\overline z)\);
2. independent circle-grade blocks;
3. rational low-rank flux factors;
4. separate real and imaginary sectors where global phase requires
   them.

The square count grows from `10` at `L=6` to `15` at `L=7`, and several
tail vectors in (6)--(7) literally extend L180's vectors.  This is
strong evidence for an all-size interval-flux LDL recurrence, but the
two certificates alone are not an induction proof.

L182 subsequently derives the arbitrary-\(L\) cubic rows
\(G_{L,k}^{(3)}\) in exactly these anti-diagonal flux coordinates.
The remaining target is the base sextic: put it in the same blocks,
then prove that the displayed low-rank elimination propagates while
retaining L173's positive null-lift curvature.  Do not run another raw
fixed-size response reconstruction.

## 6. Exact regeneration

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_L7_cubic_response.py \
  --workers 8 \
  --output \
  experiments/crabb_full_disk_complex_L7_cubic_response_s70224.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_disk_complex_L7_sixth_certificate.py \
  --output \
  experiments/crabb_full_disk_complex_L7_sixth_certificate_s70224.jsonl
```

The response checker verifies the two exact evaluation ranks and all
4,368 normal-polarization evaluations.  The certificate checker
verifies conjugation invariance, all positive factors, both Gram ranks,
and every coefficient of (8).  Clean second runs must reproduce both
JSON files byte-for-byte.
