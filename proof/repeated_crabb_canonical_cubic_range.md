# The canonical cubic obstruction has no pointwise cokernel

## 1. Result (L229, 2026-07-24)

Retain A172's canonical repaired metric.  Let

\[
\widehat T(c)=S+cA_1+O(c^2),\qquad
{\cal K}(c)=c^2K_2+c^3K_3+O(c^4),
\]

where \({\cal K}\) is the right-defect Schur residual of the
boundary Stein slack.  L227--L228 give

\[
K_2=E_1F+FE_1,\qquad E_1=S^*ES.
\]

The first canonical repair coefficient is especially simple:

\[
\boxed{X_2=-K_2.}                                  \tag{1}
\]

Indeed, \(S^*K_2S=0\), so \(X_2-S^*X_2S=-K_2\).

At cubic order put

\[
\boxed{
L_3=-K_3-A_1^*K_2S-S^*K_2A_1.}                   \tag{2}
\]

Then

\[
\boxed{
X_3={\cal G}_S(L_3),\qquad
[c^3]{\cal U}_{\rm can}
=-4W^*{\cal G}_S(L_3)W,}                          \tag{3}
\]

where \({\cal U}_{\rm can}\) is the physical upper-gap Schur
endpoint.  The state Schur square does not enter before order four.

The new exact fact is that the cubic endpoint in (3) lies in L204's
free-column response range at every fixed equality colligation:

\[
\boxed{
[c^3]{\cal U}_{\rm can}\in\operatorname {ran}{\cal M}_T.}       \tag{4}
\]

Thus A172's cubic obstruction has no pointwise Fredholm obstruction.
It can always be cancelled at a fixed colligation.  L229 does
**not** yet give a bounded analytic choice through the repeated
rank jump.

## 2. Exact ten-word forcing

L227's first operator coefficient is

\[
A_1=(I+F)S^*(I+E)-S^3.
\]

Expand the boundary slack and its right-defect Schur complement
through order three, insert (2), and reduce only by

\[
SS^*S=S,\qquad S^*SS^*=S^*,\qquad EF=FE=0.
\]

The result is the Hermitian ten-word polynomial

\[
\boxed{
\begin{aligned}
L_3={}&(S^*)^3S+S^*S^3
-(S^*)^4S^2-(S^*)^2S^4\\
&+(S^*)^4S^3S^*+S(S^*)^3S^4\\
&-4(S^*)^3S^3(S^*)^2-4S^2(S^*)^3S^3\\
&+3(S^*)^3S^3(S^*)^3S
+3S^*S^3(S^*)^3S^3 .
\end{aligned}}                                    \tag{5}
\]

This formula is universal: it uses no transfer rank assumption and
no copy-space commutation.

## 3. Cyclic cancellation

Every self-adjoint commutant \(H_0\) of \(S\) also commutes with
\(S^*\).  Hence cyclic trace reduction may move \(H_0\) through
every word in (5).

The ten words fall into four reduced cyclic classes:

\[
\begin{array}{c|c|c}
\text{words in (5)}&\text{cyclic reduction}&
\text{total coefficient}\\ \hline
(S^*)^3S,\ -(S^*)^4S^2 &(S^*)^2&1-1\\
S^*S^3,\ -(S^*)^2S^4&S^2&1-1\\
(S^*)^4S^3S^*,\
(S^*)^3S^3(S^*)^3S,\
(S^*)^3S^3(S^*)^2
&(S^*)^5S^3&1+3-4\\
S(S^*)^3S^4,\
S^*S^3(S^*)^3S^3,\
S^2(S^*)^3S^3
&(S^*)^3S^5&1+3-4 .
\end{array}                                       \tag{6}
\]

Therefore

\[
\boxed{\operatorname {tr}(H_0L_3)=0}              \tag{7}
\]

for every self-adjoint commutant \(H_0\).  Taking \(H_0=I\) gives
the observed scalar identity \(\operatorname {tr}L_3=0\).

## 4. Fredholm range proof

Let \(Y=Y^*\) annihilate L204's physical endpoint map, and let

\[
H_Y=\sum_{n\ge0}S^nWYW^*(S^*)^n
\]

be its balanced observability Gramian.  L206 proves exactly that

\[
H_YS=SH_Y.
\]

Stein adjointness and (3) give

\[
\begin{aligned}
\left\langle
Y,[c^3]{\cal U}_{\rm can}
\right\rangle
&=-4\operatorname {tr}
\left\{YW^*{\cal G}_S(L_3)W\right\}\\
&=-4\operatorname {tr}(H_YL_3)=0,                 \tag{8}
\end{aligned}
\]

where the last equality is (7).  Hence the cubic endpoint is
orthogonal to the complete cokernel of \({\cal M}_T\).
Finite-dimensional Fredholm gives (4).

This proof also covers reducible colligations: no generic
surjectivity or scalar-only cokernel assumption is used.

## 5. Evidence for removable analytic divisibility

Pointwise range membership is weaker than an analytic selection near
the repeated monomial, because \({\cal M}_T\) vanishes there.  The
tracked scaled Schur chains probe this remaining issue.

For multiplicity four, as the nonterminal Schur parameters are
scaled by

\[
\lambda=0.5,\ 0.2,\ 0.1,\ 0.05,
\]

the smallest nonzero singular value of \({\cal M}_T\) is
\(\Theta(\lambda)\), the cubic target is
\(\Theta(\lambda^3)\), and the minimum-coordinate preimage has

\[
\frac{\|C_3\|}{\lambda^2}
=0.0202533,\ 0.0202361,\ 0.0202202,\ 0.0202104.    \tag{9}
\]

Multiplicity five gives the same bounded pattern, with ratios
between \(0.01059\) and \(0.01092\).  Every range residual is below
\(8.4\times10^{-15}\).

Equation (9) is strong evidence that the cubic target is analytically
divisible by the collapsing endpoint map, but a numerical
pseudoinverse is not a proof.  The next target is an explicit
polynomial column, or a one-variable Smith/Schur argument proving
the \(O(\lambda^2)\) selection on every analytic path.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_canonical_cubic_selection.py \
  --output \
  experiments/repeated_crabb_canonical_cubic_selection_s70224.jsonl
```

The exact word audit regenerates (5) and reduces all cyclic trace
classes to zero.  The eight numerical records cover multiplicities
four and five and four scales approaching the repeated apex.  The
tracked SHA-256 is
`b6b6196c974adae1f78fa6d65eb7df8e3720197a9978aeea65fcbbae3b297e24`.
