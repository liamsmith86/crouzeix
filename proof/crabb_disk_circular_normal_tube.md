# The disk/circular-normal tube (L174, 2026-07-24)

## 1. Theorem

Fix `p=L+1`.  In L122's exact Toeplitz disk chart write

\[
z=(z_1,\ldots,z_{L-1}),\qquad
{\cal Q}(z)=\|z\|^4-|z^TJz|^2,
\]

and

\[
\widehat H_z=\widehat{\left(\frac12I+Z(z)\right)},\qquad
K_z=\widehat H_z+S^*\widehat H_zS,\qquad
X_z=2K_z^{-1/2}\widehat H_zSK_z^{-1/2}.
\]

Here \(S\) is the unweighted shift, \(Z(z)\) is the Hermitian
Toeplitz matrix with upper diagonals \(z_k\), and the hat appends one
zero row and column.  Thus \(X_z\) is the physical-gauge disk matrix.
Let

\[
{\cal E}=\{z:z=\omega J\overline z,\ |\omega|=1\}.
\]

Choose the \(2p-4\) real coefficient-gauge Crabb support-Riesz
representatives \(R_1,\ldots,R_{2p-4}\) in Fourier modes
\(3,\ldots,p\), and transport them to physical gauge by

\[
N_\nu(z)=K_z^{1/2}R_\nu K_z^{-1/2}.
\]

For the analytic rank-one envelope \(\Gamma_p\) of L118 put

\[
F(z,y)=\Gamma_p\left(X_z+\sum_\nu y_\nu N_\nu(z)\right),
\qquad y\in\mathbb R^{2p-4}.                         \tag{1}
\]

The Riemann normalization is part of \(\Gamma_p\).  Thus
\(D_yF(z,0)\) automatically uses the pulled directions
\(N_\nu(z)-h_{N_\nu(z)}(X_z)\).  This gauge transport is the one
audited in L173; omitting it gives an inconsistent pairing of the
coefficient gradient with L65's physical Hessian.

There are a neighbourhood of the origin, an analytic map
\(y_*(z)\), and constants \(a_L,b_L>0\) such that

\[
\boxed{
F(z,y)\le
-a_L{\cal Q}(z)-b_L\|y-y_*(z)\|^2\le0.
}                                                     \tag{2}
\]

Moreover

\[
\boxed{y_*(u)=0\qquad(u\in{\cal E}).}                \tag{3}
\]

Consequently the locally optimized rank-one Stein metric proves the
complete `2` bound throughout this Toeplitz-disk/circular-normal
slice tube.  This includes the singular Crabb apex and every nearby
positive phase-palindromic equality anchor in the slice.

The theorem does not yet include the non-Toeplitz base coordinates
of the full Lewis--Overton circular-range manifold, L115's one
complex elliptic soft coordinate, or L149's compact reflected
variables.  Those are the remaining full-manifold/marked merger.

## 2. Analytic chart and normal concavity

L118 constructs \(\Gamma_p\) by analytically optimizing the
rank-one Stein defect.  Near the Crabb point, the support eigenvalue
is simple with a uniform angular gap, so the numerical-range Riemann
normalization and (1) are real analytic.

At the origin,

\[
\frac12D_y^2F(0,0)
\]

is the restriction of L65's second variation to the true
circular-normal quotient.  L115 identifies the displayed support
representatives as a basis of that quotient, and L65 is negative
definite there.  Hence, after shrinking the chart,

\[
D_y^2F(z,y)\preceq-2b_LI.                            \tag{4}
\]

The analytic implicit-function theorem gives a unique critical graph

\[
D_yF(z,y_*(z))=0.                                   \tag{5}
\]

Put

\[
H(z)=F(z,y_*(z)).
\]

Taylor's theorem and (4) give the uniform fibre estimate

\[
\boxed{
F(z,y)\le H(z)-b_L\|y-y_*(z)\|^2.
}                                                     \tag{6}
\]

Thus it remains only to prove \(H(z)\le-a_L{\cal Q}(z)\).

## 3. Exact ridge contact

For \(u\in{\cal E}\), L123 gives

\[
t_*(X_u)=4
\]

and an explicit rank-one Stein metric of condition square four.
The optimized envelope is squeezed between this feasible value and
the sharp finite-Blaschke lower value, so

\[
F(u,0)=0.                                            \tag{7}
\]

L162 proves ambient stationarity of this same local rank-one branch
after Riemann normalization.  Uniqueness of L118's defect optimizer
identifies the explicit equality metric with the analytic optimized
branch.  In particular

\[
D_yF(u,0)=0.                                         \tag{8}
\]

Equations (4)--(5) make the critical point in each normal fibre
unique, so (8) proves (3).  The envelope theorem and the remaining
directions in L162 also give

\[
\boxed{
H(u)=0,\qquad DH(u)=0
\quad(u\in{\cal E}).
}                                                     \tag{9}
\]

These are exact identities on the whole positive ridge, not merely
apex Taylor cancellations.

## 4. The maximized quartic jet

On the disk fibre \(y=0\), the optimized certificate is no larger
than L122's canonical feasible metric.  The latter has fourth-order
face

\[
[F(z,0)]_4\le-32{\cal Q}(z).                        \tag{10}
\]

This inequality is enough; optimizing the defect could only make the
pure disk face more negative.  The mixed Hessian between a Toeplitz
disk tangent and a circular normal vanishes at the apex because the
former lies in L65's kernel.  Equivalently, L173's exact response
starts quadratically:

\[
D_yF(z,0)=G_2(z)+O(\|z\|^3),\qquad y_*(z)=O(\|z\|^2).
                                                               \tag{11}
\]

The upper/prepared-dual gradient difference is \(O({\cal Q})\) by
L157 and the nonnegative-gap estimate recorded in A111, so it does
not change \(G_2\).  L65 supplies the negative normal Hessian
\(B_0\).  Analytically maximizing the degree-four part of
(10)--(11) therefore gives

\[
H_4(z)
\le-32{\cal Q}(z)+\frac14G_2(z)^TB_0^{-1}G_2(z).     \tag{12}
\]

Here \(B_0\) denotes the positive matrix of the negative normal
curvature.  L173 proves that every Schur eigenvalue is strictly below
the threshold `512`.  Since

\[
{\cal Q}(z)=4\|\mathfrak p(z)\|^2,
\]

compactness of the unit Plücker sphere gives a fixed \(c_L>0\) with

\[
\boxed{H_4(z)\le-c_L{\cal Q}(z).}                   \tag{13}
\]

Modes with no quadratic response only add negative curvature and do
not affect (13).

## 5. A singular-cone blow-up lemma

We use the following extension of L152/L155.

**Lemma.**  Suppose a real-analytic germ \(h\) satisfies

\[
h|_{\cal E}=0,\qquad Dh|_{\cal E}=0,
\]

and its fourth homogeneous term obeys

\[
h_4(z)\le-c{\cal Q}(z).
\]

Then, on a sufficiently small ball,

\[
h(z)\le-a{\cal Q}(z)
\]

for some \(a>0\).

**Proof.**  L155's determinantal-ideal argument makes all homogeneous
terms below degree four zero.  Use L124's best-phase splitting

\[
z=u+v,\qquad Tu=u,\quad Tv=-v,\quad
\|v\|\le\|u\|,\qquad
{\cal Q}(z)=4\|u\|^2\|v\|^2.                        \tag{14}
\]

If \(\|v\|\ge\delta\|u\|\), then \({\cal Q}\) controls
\(\|z\|^4\).  The fifth-order Taylor remainder is absorbed by (13)
after shrinking the ball.

If \(\|v\|<\delta\|u\|\), Taylor normally to the exact ridge:

\[
h(u+v)
=\frac12D^2h(u)[v,v]
+O(\|u\|\|v\|^3+\|v\|^4).                           \tag{15}
\]

Apply the fourth-order hypothesis to \(u+tv\), divide by \(t^2\),
and let \(t\to0\).  Equation (14) gives

\[
\frac12D^2h_4(u)[v,v]
\le-4c\|u\|^2\|v\|^2.                               \tag{16}
\]

Terms of degree at least five change the left side by
\(O(\|u\|^3\|v\|^2)\).  First shrink \(\|u\|\), then choose
\(\delta\) so the cubic-normal term in (15) is a small fraction of
\(\|u\|^2\|v\|^2\).  Equations (14)--(16) prove the lemma. \(\square\)

Apply the lemma to \(h=H\), using (9) and (13).  This gives

\[
H(z)\le-a_L{\cal Q}(z).
\]

Combining this with (6) proves (2).

## 6. Scope and regeneration

The proof uses finite-dimensional analyticity and compactness for
each fixed `p`; it does not assert dimension-uniform constants.
It also does not confuse the true circular normals with all L65
strong directions.  Directions tangent to the larger
Lewis--Overton disk manifold must be charted as base variables;
L115 is what selects the \(2p-4\) normal fibre used here.  Thus L174
is a rigorous transverse slice milestone, not yet the full
circular-range tubular theorem requested in L173 Section 6.

The load-bearing algebraic inputs regenerate with

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_circular_normal_quadratic_exact.py \
  --minimum-length 3 --maximum-length 7 \
  --output \
  experiments/crabb_circular_normal_quadratic_exact_s70223.jsonl

PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_circular_normal_plucker_schur.py \
  --minimum-length 3 --maximum-length 12 \
  --resolution 1024 \
  --output \
  experiments/crabb_circular_normal_plucker_schur_s70223.jsonl
```

The nonlinear Schur probe
`experiments/crabb_circular_normal_schur_probe.py` is an independent
finite-scale guard.  Its 36 persisted cases have positive residual
and a maximum completed-gain/dual-deficit ratio about `0.0908` after
the coefficient-to-physical gauge audit.
