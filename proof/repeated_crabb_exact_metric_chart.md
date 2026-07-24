# Exact analytic metric chart near every repeated Crabb block

## 1. Result (L194, 2026-07-24)

Fix a Crabb length \(L\ge2\) and copy multiplicity \(m\ge1\).  Let

\[
A_0=C_{L+1}\otimes I_m,\qquad
P_0=\operatorname{diag}(1,2,\ldots,2,4)\otimes I_m. \tag{1}
\]

There are neighbourhoods of \(A_0\), zero, and \(P_0\) with the
following property.  For every operator \(T\) near \(A_0\), every
complex level-zero/range metric row

\[
B=(P_{01}\ \cdots\ P_{0L})\in M_{m,Lm},             \tag{2}
\]

and every Hermitian prescribed Stein Schur complement
\(R\in M_{Lm}\), there is a unique Hermitian range block

\[
C=(P_{ij})_{1\le i,j\le L}                          \tag{3}
\]

near \(\operatorname{diag}(2,\ldots,2,4)\otimes I_m\)
such that:

1. \(P-I\) has zero level-zero Schur complement;
2. \(P-T^*PT\) has prescribed range Schur complement \(R\).

The resulting metric \(P=P(T,B,R)\) is real analytic.  It satisfies

\[
P\succeq I,\qquad
P-T^*PT\succeq0\Longleftrightarrow R\succeq0.       \tag{4}
\]

The remaining upper constraint \(P\preceq4I\) is equivalent to one
analytic \(m\times m\) final-level endpoint inequality.

L194 generalizes L105 from repeated \(C_3\) to every fixed Crabb
length.  It removes convergence of forced metric coefficients as a
possible obstruction.  It does not choose the free row \(B\), prove
the endpoint sign, or establish a repeated-block neighbourhood.

## 2. Exact lower tightening

Partition level zero from levels \(1,\ldots,L\).  Given \(B,C\), set

\[
P_{00}=I+B(C-I)^{-1}B^*,\qquad
P=\begin{bmatrix}P_{00}&B\\B^*&C\end{bmatrix}.      \tag{5}
\]

Since \(C-I\succ0\) nearby, Schur complementation gives

\[
P-I\succeq0,\qquad (P-I)/(C-I)=0.                  \tag{6}
\]

Thus the complete level-zero/range row remains free, while the lower
constraint is exact to all orders.

## 3. Prescribed Stein Schur complement

Put

\[
S(T,B,C)=P-T^*PT.                                   \tag{7}
\]

At the base, \(S_{00}=I_m\) and the kernel consists of all physical
levels \(1,\ldots,L\).  Define

\[
{\cal F}(T,B,C)
=S_{KK}-S_{K0}S_{00}^{-1}S_{0K}.                   \tag{8}
\]

The chart equation is

\[
{\cal F}(T,B,C)=R.                                  \tag{9}
\]

It remains only to prove that the derivative in \(C\) is invertible.

## 4. All-length triangular Jacobian

Let the Crabb link weights be

\[
a_0=a_{L-1}=\sqrt2,\qquad
a_j=1\quad(1\le j\le L-2).                         \tag{10}
\]

For a Hermitian range direction
\(X=(X_{ij})_{1\le i,j\le L}\), direct multiplication gives

\[
\boxed{
(D_C{\cal F}\,X)_{ij}
=X_{ij}-a_{i-1}a_{j-1}X_{i-1,j-1},}                \tag{11}
\]

where every block with a zero physical index is defined to be zero.
The Schur-complement derivative contributes nothing because
\(S_{0K}=0\) at the base.

Equation (11) is triangular along every northwest/southeast block
diagonal.  If \(Y=D_C{\cal F}\,X\), its inverse is the finite
recurrence

\[
X_{ij}=Y_{ij}+a_{i-1}a_{j-1}X_{i-1,j-1}.           \tag{12}
\]

It starts on the first block row or column, where \(X=Y\).  Hence
\(D_C{\cal F}\) is an invertible real-linear map on the
\((Lm)^2\)-dimensional Hermitian range space, in every \(L,m\).

The finite-dimensional real-analytic implicit-function theorem gives
the unique solution

\[
C=C(T,B,R).                                        \tag{13}
\]

Since \(S_{00}\) remains positive, the Schur criterion gives (4).

## 5. One exact upper endpoint

Let \(J=\{0,\ldots,L-1\}\) and split the final physical level from
\(4I-P\).  Its \(J\times J\) block is positive near

\[
\operatorname{diag}(3,2,\ldots,2)\otimes I_m.
\]

Therefore

\[
P\preceq4I\Longleftrightarrow {\cal E}(T,B,R)\preceq0, \tag{14}
\]

where

\[
\boxed{
{\cal E}
=P_{LL}-4I+
P_{LJ}(4I-P_{JJ})^{-1}P_{JL}.}                     \tag{15}
\]

Thus every lower and Stein coefficient is an analytic function of
the operator, the free row, and the prescribed slack.  Only one
copy-space endpoint remains.

## 6. Role in the repeated-block attack

L193 gives the exact repeated equality manifold

\[
H=\operatorname{Toep}(Z_0,\ldots,Z_{L-1})^{-1}.
\]

L194 supplies the compatible ambient metric chart for perturbations
away from it.  The next tasks are:

1. express L193's normalized rank-\(m\) Hardy metric in L194
   coordinates;
2. derive the leading operator-valued upper endpoint in the first
   inverse-Toeplitz residual;
3. use its kernel as the next copy-space flag; and
4. add true circular normals and elliptic support branches without
   spending the same endpoint curvature twice.

The analytic chart guarantees convergence for each fixed \(L,m\), but
it does not make these sign and flag steps automatic.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_exact_metric_chart.py \
  --output \
  experiments/repeated_crabb_exact_metric_chart_s70224.jsonl
```

The checker:

1. builds the full real Jacobian of (11) on every Hermitian basis and
   verifies rank \((Lm)^2\);
2. applies (12) to every basis vector and recovers it exactly;
3. solves deterministic nearby nonlinear chart equations with
   prescribed positive Stein slack; and
4. verifies the exact lower Schur complement, recovered Stein slack,
   lower/Stein positivity, and equivalence of (14).

The standard dataset covers lengths two through five at copy
multiplicities one and two.  The proof itself is (5)--(15).
