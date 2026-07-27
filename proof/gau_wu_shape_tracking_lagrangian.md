# Shape and root tracking form the complete Gau--Wu active chart

> **Status and scope.**  The coordinate theorem in Sections 1--2 is
> exact at every finite nondegenerate Gau--Wu model.  The
> Hardy--Lagrangian identity in Section 3 is numerical evidence only.
> In particular, this note does not prove phase covariance, the
> lower-flag Euler equations, the sign of the shape form, or the
> Crouzeix conjecture.

## 1. The joint coordinates

Put \(m=n-1\).  Retain L345's joint real space \(X_\phi\), of
dimension \(n^2\), and the exact conformal shape coordinate

\[
 {\cal A}_\phi(G)
 =2(\widehat{s_G}(2),\ldots,\widehat{s_G}(n))
 \in{\mathbb C}^m.                                \tag{1}
\]

For the simple eigenvalues \(b_1,\ldots,b_{m-1}\), let \(P_j\) be
the spectral projections and let \(v_j\) be the corresponding
Blaschke-zero velocities.  Define

\[
\begin{aligned}
 \tau_0(G,v)&=(Y_1)_{00}+(Y_1)_{LL},\\
 \tau_j(G,v)&=\operatorname {tr}(P_jG)-v_j,
 \qquad 1\le j\le m-1.
\end{aligned}                                     \tag{2}
\]

Thus \(\tau_0\) is L349's endpoint trace and the remaining entries
are L348's simple-root mismatches.  Write

\[
 {\cal O}_\phi=({\cal A}_\phi,\tau):
 X_\phi\longrightarrow{\mathbb C}^{2m}.           \tag{3}
\]

## 2. Exact coordinate theorem

The restriction of \(\tau\) to pure zero motion is triangular and
invertible.  Indeed, if the operator is fixed, then at a simple root

\[
 \tau_j=-v_j,\qquad 1\le j\le m-1.                \tag{4}
\]

At the two endpoint diagonal entries, functional calculus gives
\(\dot f(0)=-f'(0)v_0\), so

\[
 \tau_0=-2f'(0)v_0.                               \tag{5}
\]

Nondegeneracy says \(f'(0)\ne0\).  Hence
\(\tau|_{\text{zero}}\) is a complex isomorphism.

L292 proves that \({\cal A}_\phi\) has real rank \(2m\) on the
physical normal slice and that its kernel is exactly L343's
\((n-2)^2\)-dimensional disk fibre \({\cal D}_\phi\).  Equations
(4)--(5) now give

\[
\boxed{
 \operatorname {rank}_{\mathbb R}{\cal O}_\phi=4m,\qquad
 \dim_{\mathbb R}\ker{\cal O}_\phi=(n-2)^2.}       \tag{6}
\]

More precisely, projection of \(\ker{\cal O}_\phi\) to the physical
space is an isomorphism onto \({\cal D}_\phi\): once a disk-fibre
direction is chosen, (4)--(5) determine the unique zero motion for
which all tracking coordinates vanish.  Thus shape plus tracking is
the complete active quotient; there is no third unrecorded boundary
sector.

## 3. The observed Hardy--Lagrangian identity

Let \(J_\phi\) be L345's nondegenerate joint scalar Hessian and split
the real shape map as \(S=(S_x;S_y)\).  Complexify the coordinate
rows in the Hardy orientation

\[
\boxed{
 {\cal P}_\phi=
 \begin{bmatrix}
  S_y+iS_x\\
  \tau
 \end{bmatrix}.}                                  \tag{7}
\]

Every complete audit gives the stronger identity

\[
\boxed{
 {\cal P}_\phi J_\phi^{-1}{\cal P}_\phi^T=0.}      \tag{8}
\]

Equation (8) is evidence, not a theorem.  Its three blocks organize
all of the open structure:

1. the shape--shape block is exactly L345/A292 phase covariance;
2. the shape--tracking block says that the conformal dual lift has
   zero endpoint/root-tracking response, hence closes L347--L350's
   lower flag; and
3. the tracking--tracking block is the reciprocal Hardy
   polarization of the same system.

The Riesz lifts \(J_\phi^{-1}{\cal P}_\phi^T\) all have strictly
upper-triangular first Blaschke image numerically.  Equivalently, if
the physical direction is first optimized over the disk fibre and
all zero motions at fixed shape \(a\), then its tracking mismatch is
purely antiholomorphic:

\[
 \tau_{\rm opt}(a)=R_\phi\overline a              \tag{9}
\]

to solver precision, with no complex-linear part.

This is a sharper target than proving the remaining diagonal and
cyclic equations separately.  The next proof attempt should derive
(8) from L342's two endpoint residual squares and L338's Szegő
zero metric.  The two endpoint squares must remain coupled: A299
already disproves detached fixed kernel-Gram comparisons.

## 4. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_shape_tracking_lagrangian.py
```

The checker verifies (4)--(6), forms all three blocks of (8)
independently, checks (9), and evaluates the lower triangle of every
Riesz first image.  It uses two models in every dimension
\(3,\ldots,8\).  Forward isotropy residuals are scaled against the
nonzero sesquilinear pairing because the most difficult model has an
ill-conditioned joint Hessian; the backward linear-solve residual is
recorded separately.  The maximum total and individual-block
isotropy residuals are respectively \(6.60\cdot10^{-7}\) and
\(1.79\cdot10^{-6}\); the maximum antiholomorphic-map and lower-flag
residuals are \(1.07\cdot10^{-10}\) and \(1.45\cdot10^{-7}\).
The largest backward residual is \(5.41\cdot10^{-17}\).

The dataset SHA-256 is
`7f00bc8167dbbf27fc96a64b3bf7cddc28c8e4c9aeb66dbb80021b0c4d35d0ac`.
