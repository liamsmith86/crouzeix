# The endpoint response is an exact quantum-Markov coboundary

## 1. Result (L280, 2026-07-25)

Let \(S\) be a finite pure partial isometry with equal orthogonal
defect frames

\[
I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0,
\]

and transfer coefficients

\[
B_k=W^*(S^*)^kV,\qquad k\geq1.
\]

Write

\[
\Phi(K)=\sum_{k\geq1}B_kKB_k^*,\qquad
\Phi^*(H)=\sum_{k\geq1}B_k^*HB_k.                 \tag{1}
\]

Matrix innerness makes \(\Phi\) bistochastic:

\[
\Phi(I)=I,\qquad \Phi^*(I)=I.
\]

For a copy matrix \(X\in M_m\), define the polarized L212 column

\[
\boxed{
{\cal C}_k(X)=Q\left\{
 S^kWX+
 \sum_{j=1}^{k-1}(S^*)^{k-j}VX^*B_j
\right\},\qquad Q=I-VV^*.}                        \tag{2}
\]

Then L204's physical endpoint response obeys the exact identity

\[
\boxed{
{\cal M}_T(P_{\rm met}^{1/2}{\cal C}_k(X))
=8\left\{
\operatorname {sym}(B_kX^*)
-\Phi\bigl(\operatorname {sym}(X^*B_k)\bigr)
\right\}.}                                        \tag{3}
\]

Here

\[
P_{\rm met}=2I-VV^*+2WW^*,\qquad
\operatorname {sym}(Z)=\frac{Z+Z^*}{2}.
\]

Equation (3) polarizes L212: setting \(X=B_k\) and multiplying the
column by \(-7/2\) recovers exactly

\[
28\{\Phi(B_k^*B_k)-B_kB_k^*\}.
\]

The decisive specialization is \(X_k=HB_k\), with \(H=H^*\).
Normal convergence permits summation over all grades and gives

\[
\boxed{
\begin{aligned}
{\cal C}(H)&=\sum_{k\geq1}{\cal C}_k(HB_k),\\
{\cal M}_T(P_{\rm met}^{1/2}{\cal C}(H))
&=8(I-\Phi\Phi^*)H.
\end{aligned}}                                    \tag{4}
\]

Thus the entire endpoint response range has one canonical
copy-space representation:

\[
\boxed{\operatorname {ran}{\cal M}_T
=\operatorname {ran}(I-\Phi\Phi^*).}              \tag{5}
\]

No kernel projection, selected singular vector, or pseudoinverse
occurs.  At a repeated monomial Crabb block,
\(\Phi\Phi^*=I\) and \({\cal C}(H)=0\), as required.

L280 does not yet finish bounded selection.  It converts that problem
to the sharply smaller **Markov Poisson inequality**

\[
\boxed{\text{find }H=H^*\text{ with }
E+8(I-\Phi\Phi^*)H\prec0.}                        \tag{6}
\]

On a proper copy flag \(U\), (6) is read after compression by
\(U^*(\cdot)U\).  Since (4) exhausts the full response range, its
compressions exhaust L222's compressed response range as well.

For L222's physical face, L279 proves the strict separator condition
for (6).  What remains is to control the state column
\({\cal C}(H)\) uniformly when the spectral gap of
\(\Phi\Phi^*\) closes.

## 2. Polarized endpoint calculation

For \(Y=Y^*\), let

\[
{\cal H}_Y-S{\cal H}_YS^*=WYW^*,\qquad
A_Y=V^*{\cal H}_YV,\qquad
R_Y=Q{\cal H}_YV.
\]

L208 proves, order-safely,

\[
YB_k-B_kA_Y
=W^*(S^*)^kR_Y+
\sum_{j=1}^{k-1}B_jR_Y^*(S^*)^{k-j}V.             \tag{7}
\]

Pair the right side with \(X\).  The state terms in (7) are exactly
the pairing of \(R_Y\) with (2).  L208's physical normalization gives

\[
\left\langle
Y,{\cal M}_T(P_{\rm met}^{1/2}{\cal C}_k(X))
\right\rangle
=8\operatorname {Re}\operatorname {tr}
X^*(YB_k-B_kA_Y).                                 \tag{8}
\]

Cyclicity and Hermitian symmetrization give

\[
\begin{aligned}
\operatorname {Re}\operatorname {tr}X^*YB_k
&=\operatorname {tr}
Y\operatorname {sym}(B_kX^*),\\
\operatorname {Re}\operatorname {tr}X^*B_kA_Y
&=\operatorname {tr}
A_Y\operatorname {sym}(X^*B_k).
\end{aligned}
\]

Stein adjointness says

\[
\operatorname {tr}(A_YK)
=\operatorname {tr}\{Y\Phi(K)\}.
\]

Substitution in (8), for every Hermitian test \(Y\), proves (3).

## 3. The Markov specialization

Put \(X=HB_k\) in (3) and sum.  The left term is

\[
\begin{aligned}
\sum_k\operatorname {sym}(B_kB_k^*H)
&=\frac12\{\Phi(I)H+H\Phi(I)\}\\
&=H.
\end{aligned}
\]

The right input is already Hermitian:

\[
\sum_k\operatorname {sym}(B_k^*HB_k)=\Phi^*(H).
\]

This proves (4).  The series in (2) converges normally: in finite
dimension the pure colligation has spectral radius below one, so its
powers and transfer coefficients decay exponentially.

## 4. Fixed points are exactly the endpoint cokernel

On the real Hilbert space of Hermitian matrices, \(\Phi\Phi^*\) is a
positive self-adjoint contraction.  Its Dirichlet form has the exact
ordered expansion

\[
\boxed{
\begin{aligned}
\langle H,(I-\Phi\Phi^*)H\rangle_{\rm HS}
&=\|H\|_F^2-\|\Phi^*(H)\|_F^2\\
&=\sum_{k\geq1}
\|HB_k-B_k\Phi^*(H)\|_F^2.
\end{aligned}}                                    \tag{9}
\]

To verify the second line, expand the squares, use
\(\sum B_k^*B_k=I\), and take the trace.  Every term is nonnegative.
Consequently

\[
\Phi\Phi^*(H)=H
\quad\Longleftrightarrow\quad
HB_k=B_k\Phi^*(H)\quad(k\geq1).                   \tag{10}
\]

L206 proves that (10) is exactly the self-adjoint commutant condition
for \(H\in\ker{\cal M}_T^*\).  Hence

\[
\ker(I-\Phi\Phi^*)=\ker{\cal M}_T^*.
\]

Taking orthogonal complements in finite dimension proves (5).

The semidefinite alternative may therefore be restated without state
variables:

\[
E+8(I-\Phi\Phi^*)H\prec0
\]

is solvable exactly when

\[
\operatorname {tr}(YE)<0
\]

for every nonzero \(Y\succeq0\) fixed by \(\Phi\Phi^*\).  L279 proves
that condition at every first active repeated flag.

## 5. Selection consequence and remaining debt

Formula (4) separates two issues that were previously mixed.

1. **Range and analyticity are closed.**  Every Markov coboundary has
   the explicit analytic state preimage \({\cal C}(H)\), and (5)
   shows these coboundaries exhaust the whole response range.
2. **Only a quantitative Poisson estimate remains.**  As a repeated
   block becomes reducible, the nonzero spectrum of
   \(I-\Phi\Phi^*\) may approach zero.  L279 controls the physical
   face on the fixed-point blocks but does not by itself bound a
   Poisson solution or its state preimage.

The next attack should use (9), not a Moore--Penrose inverse.  Along
L197/L220's analytic flag, compare the first physical endpoint energy
with the first nonzero Dirichlet energy in (9).  It is sufficient to
show that a solution of (6) can be chosen with
\(\|{\cal C}(H)\|\) bounded (and vanishing at the Crabb apex);
\(\|H\|\) itself need not remain bounded.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_markov_response.py \
  --output \
  experiments/repeated_crabb_markov_response_s70225.jsonl
```

The checker independently verifies:

1. the polarized identity (3) for arbitrary complex \(X\);
2. both bistochastic identities;
3. the summed Markov response (4);
4. the Dirichlet/intertwining identity (9); and
5. vanishing of \({\cal C}(H)\) on a gauged repeated monomial block.

The records use three unstructured noncommuting partial isometries and
one repeated monomial block.  The tracked data file has SHA-256

```text
8017e368a9a84643fb32d3bab70b7810ec3797f0ed6e112d3aa29f4ed91e49a8
```
