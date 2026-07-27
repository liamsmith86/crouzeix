# The Gau--Wu shape quotient in conformal Fourier coordinates

> **Status and scope.**  The conformal-Fourier coordinate chart below
> is exact.  In those coordinates the remaining L343 Schur form is
> phase covariant, and hence one complex Hermitian form of size
> \(n-1\), in every complete numerical audit.  The phase covariance
> and the negative sign of that Hermitian form are not yet proved.
> This note therefore halves the credible live target but does not
> establish the local or global Crouzeix conjecture.

## 1. An exact unweighted conformal chart

Retain L343's first support velocity \(s_C(\zeta)\) and write

\[
 \omega_\phi(\zeta)s_C(\zeta)=\sigma_C(\zeta),
 \qquad
 \omega_\phi>0,
 \tag{1}
\]

where \(\sigma_C\) is a real trigonometric polynomial of degree at
most \(n\).  The weighted polynomial quotient is

\[
 {\cal S}_\phi
 ={\cal P}^{\mathbb R}_n/
   \omega_\phi{\cal P}^{\mathbb R}_1.
 \tag{2}
\]

There is a second, conformally natural coordinate map

\[
\boxed{
 {\cal A}_\phi(C)
 =2\bigl(\widehat{s_C}(2),\ldots,\widehat{s_C}(n)\bigr)
 \in\mathbb C^{n-1}.}
 \tag{3}
\]

The factor two is the Schwarz-transform convention: these are the
coefficients of \(wH_{s_C}(w)\) belonging to support frequencies
\(2,\ldots,n\).  Constants and translations have no such
coefficients, so (3) descends to (2).

The descended map is injective.  Choose the unique representative
\(\sigma\) of a quotient class whose Fourier modes \(0,\pm1\)
vanish.  This gauge is unique: if
\(\omega_\phi h\), \(h\in{\cal P}^{\mathbb R}_1\), has those three
low modes zero, then
\(\int_{\mathbb T}\omega_\phi h^2\,dm=0\), hence \(h=0\).
If the coordinates in (3) vanish, then

\[
 \int_{\mathbb T}{\sigma^2\over\omega_\phi}\,dm
 =\int_{\mathbb T}\sigma s\,dm=0.
 \tag{4}
\]

Indeed, \(\sigma\) has only modes \(2,\ldots,n\) and their
conjugates, while precisely those modes of \(s=\sigma/\omega_\phi\)
vanish.  Positivity of \(\omega_\phi\) makes (4) imply
\(\sigma=0\).  Both spaces have real dimension \(2n-2\), so

\[
\boxed{{\cal S}_\phi\simeq\mathbb C^{n-1}}
\tag{5}
\]

exactly.  This removes the condition-number distortion introduced
by multiplying the support motion by \(\omega_\phi\); it does not
alter L343's quotient.

## 2. Exact-looking phase covariance

Let \(\widehat{\cal H}_\phi\) be L343's real Schur-maximum form and
let \(a={\cal A}_\phi(C)\).  Every complete audit gives

\[
\boxed{
 \widehat{\cal H}_\phi(e^{i\tau}a)
 =\widehat{\cal H}_\phi(a)
 \quad(\tau\in\mathbb R),}
\tag{6}
\]

to numerical precision.  Equivalently, if \(J(a)=ia\), then

\[
 J^T\widehat{\cal H}_\phi J=\widehat{\cal H}_\phi.
 \tag{7}
\]

Thus the real quadratic form has no complex-symmetric
``antiholomorphic'' block and is the realification of one Hermitian
matrix:

\[
\boxed{
 \widehat{\cal H}_\phi(a)=a^*{\bf H}_\phi a,
 \qquad {\bf H}_\phi={\bf H}_\phi^*.}
\tag{8}
\]

This explains the nearly paired spectra seen in the weighted L343
coordinates.  The pairing there is not exact because multiplication
by the nonconstant real weight \(\omega_\phi\) is not complex
linear for the Fourier complex structure.  The earlier tempting
claim that the weighted-coordinate matrix itself commutes with the
obvious \(J\) is therefore false.

Equation (8) is still an observed identity, not a theorem.  The
cleanest proof debt is to derive (7) from L344's fixed Toeplitz port:
the quadrature support velocities should produce two copies of the
same Hardy energy after the strict disk variable is eliminated.
One must prove that cancellation in the exact Schur complement,
not infer it from paired eigenvalues.

## 3. Revised live matrix

If (7) is proved, the live theorem becomes only

\[
\boxed{{\bf H}_\phi\preceq0}
\tag{9}
\]

on a complex space of dimension \(n-1\), rather than an arbitrary
real form of dimension \(2n-2\).  All tested matrices are strictly
negative.  Their soft directions remain genuine: the smallest
margin can be orders of magnitude below the largest eigenvalue, so
dimension-by-dimension floating-point Sylvester tests are not a
proof route.

The next useful identity should be a positive Gram formula for
\(-{\bf H}_\phi\), or equivalently a lossless formula for its dual
Schur matrix, in the fixed \(R,H\) coordinates of L344.  The
Hermitian collapse says that no separate correlation estimate ought
to survive in the correct coordinates.

## 4. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_conformal_shape_phase.py
```

The checker builds (3) directly from the unweighted first support
velocities, without passing through L343's weighted coordinate
basis.  It verifies rank \(2n-2\), kernel dimension \((n-2)^2\),
forms the exact numerical Schur complement, and measures (7) and
the complex-symmetric block independently.

The tracked audit contains three generic models in each dimension
\(4,\ldots,8\).  The largest phase-covariance residual is
\(3.85\cdot10^{-12}\), and the largest relative
complex-symmetric block is below the same scale.  Every Hermitian
matrix is negative; the softest tested eigenvalue is
\(-3.51\cdot10^{-7}\).  These signs and (7) are evidence only.

The dataset SHA-256 is
`41e0bec1e7a37a5ac981a59821191a7ebac44969593187203ffcffa4bc1a0ad6`.
