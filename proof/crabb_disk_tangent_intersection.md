# Circular-range tangent reduction at every Crabb block (2026-07-22)

## 1. Result

Let `A=C_p`, `p>=3`, and let `e_p` be the nonpositive L65 second
variation of the complete similarity square.  Near `A`, denote by
`\mathcal C_p` the matrices whose numerical ranges are circular disks
with arbitrary center and radius.  Then

\[
 \boxed{\dim_{\mathbb R}
   ((\ker e_p\cap T_A{\cal C}_p)/{\cal O}_p)=2p-4,}       \tag{1}
\]

where `O_p` is L66's affine-unitary tangent space.  Equivalently,

\[
 \boxed{\dim_{\mathbb R}
   \operatorname{image}(\ker e_p\longrightarrow
          M_p/T_A{\cal C}_p)=2.}                         \tag{2}
\]

Thus, after using the exact complete-`2` theorem on circular numerical
ranges, the apparent `2p-2`-dimensional higher-order quotient in L66
has only **one complex soft normal direction**.  It is the mode-two
direction represented by `A*`.  All remaining residual equality
directions are tangent to the circular-range manifold.

This is a structural reduction, not yet a neighbourhood theorem:
one still needs a tubular certificate uniform in the two soft normal
coordinates.  The central curve is

\[
 A+cA^*,\qquad
 W(A+cA^*)=\{z+c\bar z:|z|\le1\},                       \tag{3}
\]

so the surviving normal is precisely the infinitesimal elliptic
deformation of the disk.

## 2. Imported local disk-manifold theorem

Lewis--Overton (SIMAX 2020, Theorems 6.4, 7.4, and 8.7) prove that at
every nonzero scalar superdiagonal matrix the centered disk matrices
form an analytic manifold of real codimension `2p`.  At the Crabb
block, their tangent condition says that the first support variation
is constant in the normal angle.

Adjoining complex translations allows the first Fourier harmonic as
the moving center.  Hence the arbitrary-center circular-range
manifold is analytic of codimension `2p-2`.

Put

\[
 d_0=d_{p-1}=1/\sqrt2,\qquad d_j=1\quad(1\le j\le p-2).
\]

The normalized top support vector is L65 equation (8).  Direct
substitution shows that, up to a nonzero common real factor, its
mode-`k` support coefficient is

\[
 \sigma_k(E)=
 \sum_{m-j-1=k}d_jd_mE_{jm}
 +\sum_{m-j-1=-k}d_jd_m\overline{E_{jm}}.               \tag{4}
\]

Therefore

\[
 T_A{\cal C}_p=\{E:\sigma_k(E)=0,\quad 2\le k\le p\}.    \tag{5}
\]

These are `p-1` independent complex equations, as also follows from
the imported codimension.

Every member of `C_p` satisfies the complete Crouzeix bound by the
Berger--Okubo--Ando similarity theorem: after affine normalization it
has numerical radius one and is similar to a contraction with
condition number at most two.

## 3. Intersection with the L65 modes

Use L65's paired variables `u,v`, reduced vector `t`, and matrix
`K_{r,k}`.  For every paired mode `k>=2`, equation (4) factors as

\[
 \overline{\sigma_k(E)}=h_r^Tt,\qquad
 h_r=(1/\sqrt2,1,\ldots,1,1/\sqrt2)^T,                 \tag{6}
\]

with the evident one-coordinate normalization when `r=1`.  This is
checked simply by multiplying the aggregator in L65 (13)--(14):
its endpoint `sqrt(2)` coefficients exactly restore the endpoint
weights in (4).

Now apply L65--L66 grade by grade.

* Mode one is unrestricted by (5), because it moves the center.  Its
  one-complex-dimensional quotient is tangent to `C_p`.
* In mode two, L65 equality is `t=lambda q_{r,2}` (with the scalar
  limiting interpretation at `r=1`).  Since `h_r^Tq_{r,2}>0`,
  equation (5) forces `lambda=0`.  Of L66's two complex quotient
  coordinates, one remains tangent and one is normal.  For `p=3`,
  the same statement says that the whole exceptional unpaired
  quotient is normal.
* In every paired mode `k>=3`, L65 equality is `t=0`, which implies
  (5) automatically by (6).  The one complex quotient direction is
  tangent.
* In the unpaired mode `p-1`, L65's equality condition is exactly
  the last unitary-orbit direction and satisfies (5).  At `p=3`,
  (5) removes its extra quotient coordinate.  The bottom mode is
  already zero in the L65 kernel.

Consequently the tangent quotient contains one complex mode-one
direction, one of the two complex mode-two directions when `p>=4`,
and one complex direction in each mode `3,...,p-2`.  Its real
dimension is

\[
 2+2+2(p-4)=2p-4
\]

for `p>=4`, while it is two for `p=3`.  This proves (1).  Subtracting
from L66's total quotient dimension `2p-2` proves (2).

Finally, `A*` is a soft mode-two representative.  In L65's notation
its subdiagonal weights give `t=q_{r,2}`: the two endpoint sums are
`sqrt(2)+sqrt(2)=2sqrt(2)` and every interior sum is one.  Formula
(4) is nonzero, so this representative is normal.  Equation (3)
follows directly from

\[
 \langle(A+cA^*)x,x\rangle=z+c\bar z,\qquad
 z=\langle Ax,x\rangle,
\]

and `W(A)` being the unit disk.

## 4. Independent regeneration

Run

```bash
.venv/bin/python -u experiments/crabb_disk_tangent_intersection.py \
  --output experiments/crabb_disk_tangent_intersection_s70222.jsonl
```

The checker independently reconstructs the full L65 quadratic form,
builds (4) on all real matrix coordinates, and verifies for
`p=3,...,7`:

\[
\begin{array}{c|c}
\text{quantity}&\text{dimension}\\ \hline
\operatorname{rank}e_p&p(p-2)\\
\operatorname{codim}T_A{\cal C}_p&2p-2\\
\ker e_p&p(p+2)\\
\ker e_p\cap T_A{\cal C}_p&p^2+2p-2\\
{\cal O}_p&p^2+2\\
(\ker e_p\cap T_A{\cal C}_p)/{\cal O}_p&2p-4.
\end{array}
\]

It also verifies (6) for every paired mode and checks that all
affine-unitary generators satisfy the circular tangent equations.

## 5. Next target

Do not compute all `2p-2` residual jets.  Use a tubular chart around
`C_p`, invoke the exact disk theorem on the anchor, and retain only:

1. coercive L65 quadratic normal directions;
2. the single complex elliptic normal represented by `A*`.

L117 now covers the pure elliptic family (3) in every size and proves
the exact squared similarity constant

\[
 \frac{k(c^{2p-2})}{c^{p-1}}<4,                         \tag{7}
\]

by a Jacobi/DCT rank-one Stein metric, thereby closing the all-size
gap explicitly recorded on page 46 of Kenan Li's 2021 thesis.

The remaining dimension-independent problem is uniform transverse
absorption.  The margin in (7) above the disk limit is of order
`c^(2p-2)`, so a generic continuity argument with a fixed positive
axis gap is insufficient.  Differentiate or Schur-complement L117's
exact metric and rank-one defect in the L65 coercive directions,
using the exact circular-range theorem at `c=0` as the second anchor.
