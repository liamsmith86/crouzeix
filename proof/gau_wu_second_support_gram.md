# The Gau--Wu second-support term is a positive-weight negative Gram

> **Scope.**  This is one exact, dimension-free component of the
> arbitrary finite Gau--Wu Hessian.  It does not sign the remaining
> first-support/Hardy/zero-motion terms.  Its purpose is to replace
> part of L334's numerical form by a proved negative reserve and to
> state the smaller inequality that remains.

## 1. Result (L335/A282, 2026-07-26)

Let \(f\) be a finite Blaschke product of degree \(n-1\), with
\(f(0)=0\), and put \(\phi(z)=zf(z)\).  On the model space
\({\cal H}(\phi)\), let

\[
 S=S(\phi),\qquad
 X=\sqrt2\,P_{\ker S}+P_{\rm mid}
       +{1\over\sqrt2}P_{\ker S^*},\qquad
 A=XSX^{-1}.                                      \tag{1}
\]

Write \(p=f\in\ker S\) and \(q=1\in\ker S^*\).  Gau--Wu give

\[
 f(A)q=2p,\qquad W(A)=\overline{\mathbb D}.        \tag{2}
\]

For a physical direction \(E\), let

\[
 H_E(\zeta)={\overline\zeta E+\zeta E^*\over2},
 \qquad |\zeta|=1,
\]

and let \(v_\zeta\) be the unit top support vector of
\(\operatorname{Re}(\overline\zeta A)\).  With the reduced resolvent

\[
 R_\zeta=
 \left(I-\operatorname{Re}(\overline\zeta A)\right)^\dagger
 \quad\hbox{on }v_\zeta^\perp,
\]

put

\[
 t_E(\zeta)=
 \left\|R_\zeta^{1/2}H_E(\zeta)v_\zeta\right\|^2. \tag{3}
\]

This is the second support coefficient.  Let \(H_{t_E}\) be its
Schwarz transform and

\[
 K_{t_E}(w)=wH_{t_E}(w).                           \tag{4}
\]

The inverse-Riemann second operator contains the term
\(-K_{t_E}(A)\).  Its complete contribution to the second
coefficient of the sharp norm is exactly

\[
 \boxed{
 -{\cal P}_\phi(E),\qquad
 {\cal P}_\phi(E)
 =2\int_{\mathbb T}
 D_f(\zeta)t_E(\zeta)\,dm(\zeta),}                \tag{5}
\]

where

\[
 D_f(\zeta)={\zeta f'(\zeta)\over f(\zeta)}
 =\sum_{a\in Z(f)}
 {1-|a|^2\over|\zeta-a|^2}>0.                     \tag{6}
\]

The polarized form is

\[
 {\cal P}_\phi(E,F)
 =2\int_{\mathbb T}D_f(\zeta)
 \operatorname{Re}\left\langle
 R_\zeta^{1/2}H_Ev_\zeta,
 R_\zeta^{1/2}H_Fv_\zeta
 \right\rangle dm.                               \tag{7}
\]

Thus this whole part of the arbitrary-degree Hessian is a
dimension-free negative Gram.

## 2. Proof

For every analytic \(h\), model functional calculus gives

\[
 h(S)q=P_{{\cal H}(\phi)}h.
\]

The endpoint scalings in (1) therefore imply the transition identity

\[
 \langle p,h(A)q\rangle
 =2\langle f,h\rangle_{H^2}.                      \tag{8}
\]

The top singular value in (2) is simple.  If the second operator
coefficient is \(Y\), its linear contribution to the norm
coefficient is

\[
 \operatorname{Re}\langle p,Df(A)[Y]q\rangle.     \tag{9}
\]

For \(Y=-K_{t_E}(A)\), the direction commutes with \(A\), so

\[
 Df(A)[K_{t_E}(A)]
 =(f'K_{t_E})(A).                                 \tag{10}
\]

Equations (8)--(10) make (9)

\[
 -2\operatorname{Re}
 \langle f,f'K_{t_E}\rangle_{H^2}.
\]

On the circle, \(\overline f=1/f\) and
\(K_{t_E}(\zeta)=\zeta H_{t_E}(\zeta)\).  Hence

\[
\begin{aligned}
\operatorname{Re}\langle f,f'K_{t_E}\rangle
 &=\operatorname{Re}\int_{\mathbb T}
 {\zeta f'(\zeta)\over f(\zeta)}
 H_{t_E}(\zeta)\,dm\\
 &=\int_{\mathbb T}D_f(\zeta)t_E(\zeta)\,dm,
\end{aligned}
\]

because \(D_f\) is real and
\(\operatorname{Re}H_{t_E}=t_E\).  This proves (5).  Formula (6) is
the standard logarithmic derivative of the Blaschke factors, and
polarization proves (7).

## 3. Exact remaining gate

Let \({\cal H}_{\rm rest}(E,u)\) denote the complete joint Hessian
after removing the term (5), where \(u\) comprises all
\(2(n-1)\) real Blaschke-zero velocities.  The arbitrary finite-model
sign is now exactly

\[
 \boxed{
 \max_u{\cal H}_{\rm rest}(E,u)
 \le {\cal P}_\phi(E).}                           \tag{11}
\]

The reserve on the right is explicit, positive, and uses the same
support reduced resolvent already present in the Riemann jet.

Numerically, the optimized remainder in (11) is indefinite.  Its
positive index is \(2(n-2)\) in every L334 sample, while all other
directions are already favorable.  This strongly suggests that only
the model space belonging to the nonzero-zero factor \(f/z\) is
dangerous.  That inertia pattern is evidence, not part of L335.

A tempting shortcut was also rejected: the linearized rank-one
partial-isometry defect from Gau--Wu's Ando factorization misses two
Hessian directions.  It cannot replace (11) without an additional
model/characteristic-function constraint.

## 4. Regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/gau_wu_second_support_gram.py \
  --dimensions 3,4,5,6,7,8 --samples 2 \
  --angle-count 512 \
  --output experiments/gau_wu_second_support_gram_s70224.jsonl
```

The checker evaluates both sides of the fully polarized identity on
every normal basis.  Twelve models pass; the maximum matrix-entry
residual is \(5.33\cdot10^{-15}\), and every sampled Gram is positive
definite.  Dataset SHA-256:
`cb068337f93287075e54239265d493e73108305bd9f98082354b01de16663e93`.
