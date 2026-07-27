# The Gau--Wu Hessian is a symmetric endpoint-defect flux

> **Scope.**  The endpoint identity below is exact and
> dimension-free.  The positivity of its symmetric Stein flux has
> survived complete finite-model tests but is not yet proved.  Thus
> this note sharpens the arbitrary-degree gate; it does not close it
> or establish a local theorem beyond L333.

## 1. Exact endpoint identity (L336/A283, 2026-07-26)

Let \(S=S(zf)\) be the compressed shift in L335 and use the
orthogonal decomposition

\[
 {\cal H}(zf)=\mathbb Cp\oplus{\cal M}\oplus\mathbb Cq,
 \qquad p=f,\quad q=1.
\]

The Gau--Wu similarity is

\[
 X=\operatorname{diag}(\sqrt2,I_{\cal M},1/\sqrt2),
 \qquad A=XSX^{-1}.
\]

Consider any joint second-order jet after conformal normalization
and motion of all zeros of \(f\):

\[
 Y_\varepsilon=f_\varepsilon(T_\varepsilon)
 =pq^*+\varepsilon Y_1+\varepsilon^2Y_2+o(\varepsilon^2),
 \qquad
 \operatorname{Re}\langle p,Y_1q\rangle=0.        \tag{1}
\]

Write

\[
\begin{aligned}
 r_2&=[\varepsilon^2]\,
 \langle q,(I-Y_\varepsilon^*Y_\varepsilon)q\rangle,\\
 \ell_2&=[\varepsilon^2]\,
 \langle p,(I-Y_\varepsilon Y_\varepsilon^*)p\rangle.
\end{aligned}                                     \tag{2}
\]

Then the second coefficient of the sharp scaled norm is exactly

\[
\boxed{
 [\varepsilon^2]\|XY_\varepsilon X^{-1}\|
 =-\frac12(r_2+\ell_2)
  -\frac14\left(
   |\langle q,Y_1q\rangle|^2+
   |\langle p,Y_1p\rangle|^2
  \right).}                                      \tag{3}
\]

Thus the arbitrary finite Gau--Wu Hessian is nonpositive as soon as

\[
\boxed{r_2+\ell_2\ge0.}                            \tag{4}
\]

The last two terms in (3) are an additional reserve.  Inequality
(4) is slightly stronger than the exact sign target, but every
complete test below is strictly positive.

## 2. Proof of the identity

Put

\[
\begin{array}{lll}
 a=\langle p,Y_1q\rangle,&
 b=P_{\cal M}Y_1q,&c=\langle q,Y_1q\rangle,\\
 d=\langle p,Y_1p\rangle,&
 e=P_{\cal M}Y_1^*p.&
\end{array}
\]

Simple-singular-value perturbation at
\(XY_0X^{-1}=2pq^*\) gives

\[
\begin{aligned}
 J={}&[\varepsilon^2]\|XY_\varepsilon X^{-1}\|\\
 ={}&2\operatorname{Re}\langle p,Y_2q\rangle
 +|a|^2+\frac12\|b\|^2+\frac14|c|^2
 +\frac14|d|^2+\frac12\|e\|^2.                   \tag{5}
\end{aligned}
\]

On the other hand,

\[
\begin{aligned}
 r_2&=-2\operatorname{Re}\langle p,Y_2q\rangle
      -|a|^2-\|b\|^2-|c|^2,\\
 \ell_2&=-2\operatorname{Re}\langle p,Y_2q\rangle
      -|a|^2-\|e\|^2-|d|^2.                      \tag{6}
\end{aligned}
\]

Substitution of (6) in the right side of (3) gives (5).

## 3. Why the symmetric flux is the right gate

For any orthonormal basis \(e_1,\ldots,e_{n-1}\) of
\({\cal K}_{f_\varepsilon}\), L127's Blaschke--Stein identity and
its left-handed companion give

\[
\begin{aligned}
 I-Y_\varepsilon^*Y_\varepsilon
 &=\sum_j e_j(T_\varepsilon)^*
   (I-T_\varepsilon^*T_\varepsilon)
   e_j(T_\varepsilon),\\
 I-Y_\varepsilon Y_\varepsilon^*
 &=\sum_j e_j(T_\varepsilon)
   (I-T_\varepsilon T_\varepsilon^*)
   e_j(T_\varepsilon)^* .                        \tag{7}
\end{aligned}
\]

Consequently \(r_2+\ell_2\) is a two-ended scalar Stein flux, not an
unstructured remainder of the Riemann calculation.  This reformulates
L335's Hardy comparison as the positivity of the sum of the initial
and final defect fluxes in (7).

The symmetrization is load-bearing.  In dimensions \(4,\ldots,8\),
the right and left forms separately each have two negative
directions.  A proof that bounds either endpoint by itself is
therefore false.  Their sum is positive definite in every complete
test.

## 4. The two first-order residuals

Let

\[
 T_\varepsilon=S+\varepsilon C+O(\varepsilon^2)
\]

be the shift-coordinate normalized operator, let
\(h=\dot f\), and put \(\phi=zf\).  The natural residuals are

\[
\boxed{
 {\cal R}=S^*C+C^*S,\qquad
 {\cal Z}=D\phi(S)[C]+(zh)(S).}                   \tag{8}
\]

They respectively measure failure of the contraction defect and
failure of the moving characteristic equation
\(\phi_\varepsilon(T_\varepsilon)=0\).  On every tested full joint
normal space their combined real linear map has rank exactly
\(n^2\), the whole joint dimension.  Individually their ranks are
\((n-1)^2+1\) and \(2n\), with a two-dimensional overlap.

This rank statement is numerical evidence, not yet a theorem.  It
does, however, explain why the complete joint Hessian has rank
\(n^2\) and why the raw one-defect candidate rejected after L335
missed two directions.

## 5. Exact remaining task

Polarize (7), insert the inverse-Riemann jet, and prove

\[
 [\varepsilon^2]\left\{
 \langle q,(I-Y_\varepsilon^*Y_\varepsilon)q\rangle+
 \langle p,(I-Y_\varepsilon Y_\varepsilon^*)p\rangle
 \right\}\ge0.                                   \tag{9}
\]

The proof must keep the two endpoints together.  The most credible
route is to use (8) and the model decomposition
\({\cal K}_{zf}={\cal K}_f\oplus\mathbb Cf\) to factor (9) as a
positive quadratic form in \(({\cal R},{\cal Z})\).  Do not return to
the separately false endpoint signs, a pointwise bound
\(t_E\ge\delta_E^2/2\), or the incomplete raw
partial-isometry defect.

## 6. Regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/gau_wu_two_sided_endpoint_defect.py \
  --dimensions 3,4,5,6,7,8 --samples 2 \
  --angle-count 512 \
  --output \
  experiments/gau_wu_two_sided_endpoint_defect_s70224.jsonl
```

The checker constructs the complete joint form in all \(n^2\)
physical-normal/zero-velocity variables, verifies (3) against the
independent norm-jet implementation, records the separate endpoint
inertias, tests positivity of their sum, and audits the rank of (8).

The tracked 12-record standard dataset has SHA-256
`023ef06315c2524d1dcc0e818698c6ddc1ab1fea4562b3272955b93b5f802cbb`.
An additional 42 untracked holdout models in dimensions \(3,\ldots,9\)
also had positive two-sided forms; the softest holdout eigenvalue was
\(3.43\cdot10^{-11}\).  All standard models in dimensions
\(4,\ldots,8\) had two negative directions in each one-sided form.
The holdouts confirmed one-sided indefiniteness and generically the
same index.
