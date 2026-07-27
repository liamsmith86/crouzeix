# Dual Hardy polarization of the Gau--Wu shape form

> **Status and scope.**  The dual-isotropy reformulation of A292's
> phase covariance is exact finite linear algebra.  The isotropy
> itself, and the associated upper-triangular response law, are
> numerical observations.  This note identifies a concrete
> Euler--Lagrange/Hardy projection to prove; it does not yet promote
> the Hermitian collapse or sign the shape form.

## 1. The joint scalar Hessian

Let \(X_\phi\) be the real space consisting of L342's physical
normal direction and all first Blaschke-zero velocities.  Its real
dimension is

\[
 \dim_{\mathbb R}X_\phi
 =\{(n-1)^2+1\}+2(n-1)=n^2.
 \tag{1}
\]

Let \(J_\phi\) be the joint scalar Hessian on this space, before
the zero velocities and L343 disk fibre are eliminated.  L342
identifies its physical Schur complement with one fourth of the
optimized similarity Hessian.  Extend A292's conformal shape map by
zero on the zero-velocity variables:

\[
 S_\phi:X_\phi\longrightarrow
 \mathbb R^{2m},\qquad
 m=n-1,
 \tag{2}
\]

where the coordinates are ordered as
\((\operatorname {Re}a,\operatorname {Im}a)\) for

\[
 a=2(\widehat s(2),\ldots,\widehat s(n)).
 \tag{3}
\]

Suppose \(J_\phi\) is nondegenerate.  The inverse form induced on
the shape dual is

\[
 K_\phi=S_\phi J_\phi^{-1}S_\phi^T.
 \tag{4}
\]

The nuisance restriction is the strict zero/disk-fibre block from
L338 and L343, so the constrained Schur form is \(K_\phi^{-1}\)
whenever \(J_\phi\) is nondegenerate.  Thus phase covariance may be
checked before either the disk or zero Schur complement is
explicitly formed.  This criterion does not itself prove
nondegeneracy of the remaining shape form.

## 2. Exact dual-isotropy equivalence

Write

\[
 S_\phi=\begin{bmatrix}S_x\\S_y\end{bmatrix},
 \qquad
 V_x=J_\phi^{-1}S_x^T,\quad
 V_y=J_\phi^{-1}S_y^T,
 \tag{5}
\]

and form the complexified dual lift

\[
 \boxed{W_\phi=V_y+iV_x.}
 \tag{6}
\]

If
\[
 K_\phi=
 \begin{bmatrix}
  K_{xx}&K_{xy}\\K_{yx}&K_{yy}
 \end{bmatrix},
\]
then direct expansion gives

\[
 \boxed{
 W_\phi^TJ_\phi W_\phi
 =K_{yy}-K_{xx}+i(K_{yx}+K_{xy}).}
 \tag{7}
\]

Since \(K_\phi\) is real symmetric, the following are equivalent:

1. \(K_\phi\) commutes with
   \({\cal J}=\left[\begin{smallmatrix}0&-I\\I&0\end{smallmatrix}\right]\);
2. \(K_{xx}=K_{yy}\) and \(K_{xy}=-K_{xy}^T\);
3. the complex dual range in (6) is totally isotropic:

   \[
   \boxed{W_\phi^TJ_\phi W_\phi=0.}               \tag{8}
   \]

The same phase identity then holds for \(K_\phi^{-1}\), so (8) is
exactly A292's missing Hermitian-collapse theorem.  It is not a
consequence of a paired numerical spectrum: it is one explicit
bilinear identity in the uneliminated \(n^2\)-variable Hessian.

## 3. The observed Hardy polarization

For \(v\in X_\phi^{\mathbb C}\), let \(Y_1(v)\) be the complex-linear
extension of the first Blaschke-image response in L336.  Every
complete audit gives

\[
 \boxed{\operatorname {tril}Y_1(W_\phi)=0,}       \tag{9}
\]

including the diagonal.  Thus the dual lift selected by the
conformal positive-frequency source has a purely strict
upper-triangular first image.  This is the finite Schur-flag version
of an analytic/Hardy polarization.

The stronger numerical diagnostic is that the right and left
endpoint forms of L336 are separately bilinearly isotropic on this
complex range:

\[
 W_\phi^TR_\phi W_\phi=0,\qquad
 W_\phi^TL_\phi W_\phi=0.                         \tag{10}
\]

Their real signs remain indefinite, as L336 already warned; (10) is
a complex bilinear statement and does not revive either false
one-sided positivity claim.

Equations (9)--(10) are not yet proved.  The next exact task is:

1. write the Euler equation
   \(J_\phi W_\phi=S_y^T+iS_x^T\) in L344's fixed
   inverse-Toeplitz port coordinates;
2. use the positive-frequency support source and the upper
   Takenaka/Schur flag to prove (9);
3. polarize L336's exact two endpoint identities on that flag and
   derive (10), hence (8).

This route proves phase covariance without assuming the desired
negative sign.  Only after (8) is closed should the retained
Hermitian form be signed.

## 4. Falsified shortcut

The inverse negative Hermitian matrices have one very large soft
rank-one component in many models.  That dominance can make a
universal-diagonal-plus-rank-one ansatz appear accurate in ordinary
relative norm.  Removing the leading component exposes additional
indefinite model-dependent corrections.  Therefore no such
rank-one formula is used or recorded as a conjecture.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_conformal_dual_polarization.py
```

The tracked audit contains three models in every dimension
\(4,\ldots,8\).  It forms (4)--(8) from the independent joint scalar
Hessian, evaluates (9), and also checks (10).  The most ill
conditioned joint form has condition number
\(1.61\cdot10^{11}\), while the largest backward solve residual is
\(9.05\cdot10^{-17}\).  Its corresponding forward residuals are the
largest in the set:

\[
\begin{array}{c|c}
\text{phase covariance}&4.29\cdot10^{-6}\\
\text{sharp isotropy}&2.30\cdot10^{-6}\\
\text{right/left isotropy}&3.60\cdot10^{-6}\\
\text{lower-triangular response}&2.60\cdot10^{-7}.
\end{array}
\]

The better-conditioned models satisfy the same identities many
orders more accurately.  These figures are evidence only; A292's
better-conditioned direct phase audit remains the primary numerical
test.

The dataset SHA-256 is
`d02f3734280a2a5594fba11cdcc1514ad6b760c6b3fb112c1b76bf5b733ce8db`.
