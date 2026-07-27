# The Gau--Wu support port in fixed Toeplitz coordinates

> **Status and scope.**  The spectral-factor transport and fixed
> triangular-resolvent formula below are exact at every finite
> nondegenerate Gau--Wu equality model.  They remove L339's moving
> boundary pseudoinverse and put its port energy in L343's
> inverse-Toeplitz coordinates.  They do not sign L343's remaining
> boundary-shape Schur form.

## 1. Two factors of the same support slack

Retain L339's Ando factor
\[
 F_\zeta=UD-\zeta E,\qquad
 F_\zeta^*F_\zeta
 =I-\operatorname {Re}(\overline\zeta A),          \tag{1}
\]
and L343's polynomial support frame \({\cal C}\).  Thus
\[
 K={\cal C}^*{\cal C},\qquad
 K=\widehat H+R^*\widehat HR,\qquad
 KA_c=2\widehat HR,                                \tag{2}
\]
where \(R\) is the unweighted superdiagonal shift,
\(\widehat H=H\oplus0\), and \(H\succ0\).

Put
\[
 B_\zeta=\widehat H^{1/2}(I-\overline\zeta R).
 \tag{3}
\]
L343's support-pencil factorization gives
\[
\boxed{
 (F_\zeta{\cal C})^*(F_\zeta{\cal C})
 =B_\zeta^*B_\zeta.}                               \tag{4}
\]
Both factors have rank \(n-1\) and the same one-dimensional right
kernel.  Their polar decompositions therefore give a partial
isometry \(V_\zeta\), unitary between their ranges, such that
\[
 F_\zeta{\cal C}=V_\zeta B_\zeta.                 \tag{5}
\]
No regularity or formula for \(V_\zeta\) will be needed; only its
isometry on the two displayed ranges is used.

## 2. Transport of a minimal port response

Let \(r_\zeta\in\operatorname {ran}F_\zeta^*\), and let
\[
 \gamma_\zeta=(F_\zeta^*)^\dagger r_\zeta          \tag{6}
\]
be the minimal solution of \(F_\zeta^*\gamma_\zeta=r_\zeta\).
It lies in \(\operatorname {ran}F_\zeta\).  Equation (5) shows that
\[
 \eta_\zeta=V_\zeta^*\gamma_\zeta
\]
is the minimal solution of
\[
 B_\zeta^*\eta_\zeta={\cal C}^*r_\zeta,
 \qquad
 \|\eta_\zeta\|=\|\gamma_\zeta\|.                 \tag{7}
\]

Since
\[
 B_\zeta^*
 =(I-\zeta R^*)\widehat H^{1/2},
\]
the first factor is invertible and
\(\operatorname {ran}\widehat H^{1/2}=e_L^\perp\).
Solvability in (7) says exactly that
\[
 (I-\zeta R^*)^{-1}{\cal C}^*r_\zeta\in e_L^\perp.
 \tag{8}
\]
The unique minimal solution is consequently
\[
\boxed{
 \eta_\zeta
 =\widehat H^{\dagger/2}
  (I-\zeta R^*)^{-1}{\cal C}^*r_\zeta,\qquad
 \|\gamma_\zeta\|=\|\eta_\zeta\|.}                \tag{9}
\]
The only pseudoinverse in (9) is the fixed matrix
\(\widehat H^{\dagger/2}=H^{-1/2}\oplus0\).  Its rank and kernel do
not vary with \(\zeta\) or with the perturbation direction.

## 3. The L339 port in these coordinates

For an ambient direction \(C\), write
\[
 {\cal H}_C(\zeta)
 ={\,\overline\zeta C+\zeta C^*\over2}.
\]
Let \(y_\zeta\) be L339's unnormalized support null state and
\[
 s_C(\zeta)
 ={\langle y_\zeta,{\cal H}_C(\zeta)y_\zeta\rangle
   \over\|y_\zeta\|^2}.
\]
Then
\[
 r_C(\zeta)
 ={\cal H}_C(\zeta)y_\zeta-s_C(\zeta)y_\zeta
 \in y_\zeta^\perp=\operatorname {ran}F_\zeta^*.
 \tag{10}
\]
The component along \(y_\zeta\) is automatically discarded by
L339's Moore--Penrose inverse, so its response is (6) with
\(r_\zeta=r_C(\zeta)\).  Equations (9)--(10) turn L339's reserve into
\[
\boxed{
 {\cal P}_\phi(C)
 =\int_{\mathbb T}
 \left\|
 \widehat H^{\dagger/2}(I-\zeta R^*)^{-1}
 {\cal C}^*
 \{{\cal H}_C(\zeta)y_\zeta-s_C(\zeta)y_\zeta\}
 \right\|^2dm(\zeta).}                            \tag{11}
\]

Thus the support port is an ordinary finite Hardy norm for one fixed
nilpotent shift and one fixed inverse-Toeplitz weight.  In
particular, L339's pointwise pseudoinverse is not an obstruction to
analytic or coefficientwise manipulation.

## 4. Compatibility with the L343 shape quotient

Let
\[
 \delta(t)=\det(I-tS),\qquad
 \alpha_\zeta=\zeta^{n-1}\delta(\overline\zeta),
\qquad
 {\bf f}_\zeta=(1,\zeta,\ldots,\zeta^{n-1})^T.
\]
L343 gives
\[
 {\cal C}{\bf f}_\zeta=\alpha_\zeta y_\zeta,\qquad
 \omega_\phi=|\alpha_\zeta|^2D_f.
 \tag{12}
\]
Consequently the forcing in (11), after multiplication by
\(\alpha_\zeta\), is
\[
\boxed{
 {\cal C}^*{\cal H}_C(\zeta){\cal C}{\bf f}_\zeta
 -s_C(\zeta)K{\bf f}_\zeta.}                      \tag{13}
\]
Multiplying once more by the fixed polynomial \(\omega_\phi\) gives
\[
 \omega_\phi\alpha_\zeta{\cal C}^*r_C(\zeta)
 =\omega_\phi{\cal C}^*{\cal H}_C(\zeta)
       {\cal C}{\bf f}_\zeta
  -\sigma_C(\zeta)K{\bf f}_\zeta,                 \tag{14}
\]
where \(\sigma_C=\omega_\phi s_C\) is exactly L343's finite
polynomial coordinate, modulo the three affine modes.
Equations (11)--(14) therefore put both pieces of the live problem
in one coefficient system:

1. the remaining direction is
   \([\sigma_C]\in{\cal P}_n^{\mathbb R}/
   \omega_\phi{\cal P}_1^{\mathbb R}\), of dimension \(2n-2\);
2. its retained reserve is the fixed \(R,H\) Hardy energy (11).

The next step is to eliminate the already strict general-\(H\) disk
direction from (14) and compare the resulting endpoint output with
(11).  No moving spectral factor or rank-dependent boundary
pseudoinverse remains.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_toeplitz_port_transport.py
```

The checker reconstructs \({\cal C}\) and \(\widehat H\)
independently, verifies (4), tests the fixed-range condition (8), and
compares the original and transported responses pointwise and after
full real polarization over the physical normal space.  It also
compares the transported Gram with L339's independently computed
weighted second-support Gram.

Twelve models, two in every dimension \(3,\ldots,8\), give maximum
spectral-factor, fixed-range, relative pointwise-energy,
transport-Gram, and support-Gram residuals respectively
\[
 6.62\cdot10^{-14},\quad
 3.85\cdot10^{-15},\quad
 2.68\cdot10^{-14},\quad
 1.44\cdot10^{-13},\quad
 1.32\cdot10^{-13}.
\]
The dataset SHA-256 is
`0e800d68f9e233738a51028b72ca9ce728f1990bb55286984e8125bbb187d874`.
