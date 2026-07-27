# The Gau--Wu lower Schur flag has only \(n+1\) generators

> **Status and scope.**  The commutator recurrence and its
> \(n+1\)-generator consequence are exact at every finite
> nondegenerate Gau--Wu model with distinct nonzero interior zeros.
> They reduce L345's observed upper-triangular dual response to
> \(n+1\) scalar Euler identities.  Those identities themselves
> remain open, so this note does not yet prove phase covariance.

## 1. Differentiated commutation

Work in the upper-triangular Gau--Wu/Takenaka basis.  Put

\[
 L=n-1,\qquad
 S=S_\phi,\qquad
 Y_0=f(S)=e_0e_L^*.
 \tag{1}
\]

In Gau--Wu's diagonally similar matrix coordinates the last member
of (1) is \(2e_0e_L^*\); that harmless factor does not affect the
argument.  The diagonal of \(S\) is

\[
 \lambda_0=0,\quad
 \lambda_1=b_1,\ldots,\lambda_{L-1}=b_{L-1},
 \quad\lambda_L=0,                                \tag{2}
\]

where the interior zeros are nonzero and pairwise distinct.

For any simultaneous first operator/inner-function direction, write

\[
 T_\varepsilon=S+\varepsilon C+O(\varepsilon^2),
 \qquad
 Y_\varepsilon=f_\varepsilon(T_\varepsilon)
 =Y_0+\varepsilon Y_1+O(\varepsilon^2).
 \tag{3}
\]

Functional calculus gives
\(T_\varepsilon Y_\varepsilon
 =Y_\varepsilon T_\varepsilon\).  Differentiating yields

\[
 \boxed{
 [S,Y_1]=Y_0C-CY_0.}                              \tag{4}
\]

The right side has support only in the first row and final column.
It therefore has zero strict lower triangle:

\[
 [S,Y_1]_{ij}=0\qquad(i>j).                       \tag{5}
\]

This conclusion is independent of the Riemann acceleration, the
zero acceleration, and the Hessian sign.

## 2. Southwest recurrence

For an arbitrary matrix \(Y\), upper triangularity of \(S\) gives,
whenever \(i>j\),

\[
\begin{aligned}
 [S,Y]_{ij}
 ={}&(\lambda_i-\lambda_j)Y_{ij}\\
 &+\sum_{k=i+1}^{L}S_{ik}Y_{kj}
 -\sum_{k=0}^{j-1}Y_{ik}S_{kj}.                  \tag{6}
\end{aligned}
\]

Every entry on the second line is strictly southwest of
\((i,j)\).  More precisely, its index

\[
 d(i,j)=(L-i)+j                                  \tag{7}
\]

is smaller.  Start at \(d=0\), the southwest corner
\((L,0)\), and induct upwards.  Except at that corner,

\[
 \lambda_i-\lambda_j\ne0\qquad(i>j),              \tag{8}
\]

by (2).  Equations (5)--(8) therefore determine every strict lower
entry uniquely from the single scalar \(Y_{L0}\).

The diagonal entries do not occur in (6).  Consequently:

\[
\boxed{
 \text{the complete lower triangle of }Y_1
 \text{ is determined by }
 \bigl(\operatorname {diag}Y_1,(Y_1)_{L0}\bigr).}
 \tag{9}
\]

Equivalently,

\[
\boxed{
 \operatorname {tril}Y_1=0
 \quad\Longleftrightarrow\quad
 \operatorname {diag}Y_1=0
 \ \text{and}\ (Y_1)_{L0}=0.}                    \tag{10}
\]

Thus the lower-response map has complex rank at most \(n+1\).
The audit below finds the \(n\) diagonal evaluations and corner
independent on every tested full joint quotient, but that additional
rank equality is not needed for (10) and is not used as a theorem.

## 3. Consequence for L345

Let

\[
 W_\phi
 =J_\phi^{-1}S_y^*+iJ_\phi^{-1}S_x^*
 \tag{11}
\]

be L345's conditional complex dual lift.  A293 observed
\(\operatorname {tril}Y_1(W_\phi)=0\).  By (10), the exact missing
Euler content is only

\[
\boxed{
 \operatorname {diag}Y_1(W_\phi)=0,\qquad
 (Y_1(W_\phi))_{L0}=0.}                           \tag{12}
\]

The diagonal equations are the moving-root/interpolation part of
the joint response.  The corner is the single cyclic closing
condition.  Prove them in L344's fixed Toeplitz coordinates; do not
expand or estimate the remaining lower entries separately, because
(6) then closes them automatically.

At an interior-zero collision, some denominators in (8) vanish and
the recurrence becomes confluent.  Such ramified charts need their
own divided-derivative version; no collision-uniform claim is made
here.

## 4. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_lower_flag_recurrence.py
```

The checker differentiates commutation independently on every real
joint basis direction.  In two models per dimension \(3,\ldots,8\),
the lower-response rank and the diagonal-plus-corner generator rank
are both exactly \(n+1\).  The maximum commutator residual is
\(2.28\cdot10^{-15}\); the largest lower-row-space reconstruction
residual is \(4.92\cdot10^{-9}\).

The dataset SHA-256 is
`537bb48cebeb8c79ff3ba84df566f97a6fe01b5efeb1cd56ce81f70ece7476f3`.
