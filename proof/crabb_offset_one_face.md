# Universal offset-one elliptic Stein face (2026-07-23)

## 1. Result

Let `L>=2`, `p=L+1`, and take the real phase-one disk-equality
direction supported at coefficient offsets `1` and `L-1` (only one
copy when `L=2`).  Along A84's exact ellipse curve, L118's locally
optimized rank-one envelope satisfies

\[
\boxed{
[a^2]\Gamma_{L,1}(a,c)=-64c^2+O(c^3).
}                                                     \tag{1}
\]

This is an all-size coefficient calculation, not extrapolation from a
finite grid.  It proves the offset-one upper face needed by L132.
Consequently, if `L=qk`, the real phase-one grade-`k` direction has
the rank-one upper face

\[
-64a^2c^{2k}                                         \tag{2}
\]

after degree-`k` descent to the size-`q+1` offset-one family.

L134 subsequently polarizes the imaginary phase and proves
`-64|u|²c²` for an arbitrary complex one-pair coefficient.  Equation
(1) by itself is still not the full A84 theorem: distinct-grade
polarization, nondivisor localization, and the uniform analytic
remainder required by L124 remain open.

## 2. The explicit defect jet

Use L123's coefficient coordinates and write

\[
\begin{aligned}
S(a,c)&=C+aE+cJ(C+aE)J,\\
T(a,c)&=\phi_c(S(a,c)),\\
K(a)&=K_0+aK_1.
\end{aligned}                                        \tag{3}
\]

Let `d_0(c)` be L117's exact rank-one defect on the elliptic axis.
For a defect `d=d_0+a x`, let `Q_L(x,c)` be the coefficient of `a^2`
in the generalized endpoint condition ratio of its Stein Gramian.
This is an exact quadratic polynomial in `x`.
An `a^2` term in the defect would pair with the first defect
derivative of the condition at the axis optimizer, which is zero; it
therefore cannot change `Q_L`.

At `c=0`, the disk-equality minimizer is

\[
x^{(0)}=2\sum_{j\in\{1,L-1\}}e_j,                    \tag{4}
\]

where the set convention prevents doubling when `L=2`.  Direct
coefficient-coordinate perturbation gives

\[
Q_L(x^{(0)}+\delta,0)
=4\sum_{j=1}^{L-1}\delta_j^2+\frac83\delta_L^2.       \tag{5}
\]

Thus the defect Hessian is a unit over the power-series ring.  The
linear elliptic stationarity equation is especially sparse:

\[
\begin{array}{c|c|c}
L &[c]\nabla_xQ_L(x^{(0)},c)
  &\frac12\partial_{x_3}^2Q_L(x^{(0)},0)\\ \hline
2&0&\text{coordinate absent}\\
3&(128/3)e_3&8/3\\
L\ge4&64e_3&4 .
\end{array}                                          \tag{6}
\]

All unlisted gradient coordinates vanish.  In every case the
first-order stationarity correction is therefore

\[
\boxed{
\widehat x(c)=x^{(0)}-8c\,e_3
O(c^2),}                                             \tag{7}
\]

with the correction omitted when `L=2` and combined with (4) when
coordinates collide.

For clarity, (6) is a finite shift calculation, not a numerical
optimizer statement.  Through the required order the elliptic map is

\[
\phi_c(w)=(1+2c^2)w-cw^3+c^2w^5+O(c^3).              \tag{8}
\]

Substitute (8) in (3), solve the coefficient of `a^2` in

\[
M-T^*MT=dd^*,
\]

and use

\[
({\cal I}-C^*(\cdot)C)^{-1}F
=\sum_{n=0}^{L}(C^*)^nFC^n.                          \tag{9}
\]

At elliptic order one, only the disk seeds in (4) and coordinate
three occur.  Formula (9) gives exactly the three rows of (6).
This also shows directly that

\[
\nabla_xQ_L(\widehat x(c),c)=O(c^2).                 \tag{10}
\]

The analytic minimizer selected by L118 consequently obeys
`x_*(c)=widehat x(c)+O(c^2)`.  Since `Q_L` is quadratic with an
invertible Hessian, replacing `x_*` by `widehat x` changes its value
only by `O(c^4)`.  It is therefore enough to calculate (1) with the
explicit defect jet (7).

## 3. Endpoint forcing calculation

Write

\[
M=M_0+aM_1+a^2M_2+O(a^3)
\]

for the Stein Gramian produced by (7), and put

\[
F^{(2)}=[c^2]\{M_2-C^*M_2C\}.                        \tag{11}
\]

For `L>=7`, direct substitution of (8) and (7) gives the separated-end
diagonal

\[
\boxed{
\operatorname{diag}F^{(2)}
=(48,-96,-32,16,0,\ldots,0,16,-32,64).
}                                                     \tag{12}
\]

Nothing is hidden in the middle zeros.  At order `c^2`, every term of
(8) is a word of degree at most five.  After the disk-flat
cancellations and the correction in (7), the surviving diagonal
forcing lies within three sites of the lower end or two sites of the
upper end.  Those neighbourhoods are disjoint for `L>=7`; increasing
`L` merely inserts another zero in (12).  This is the all-size
stabilization step.

The short chains contain boundary collisions.  Direct use of the same
recurrence gives:

\[
\begin{array}{c|l}
L&\operatorname{diag}F^{(2)}\\ \hline
2&(48,-128,64)\\
3&(48,-80,-96,112)\\
4&(48,-96,-16,-16,64)\\
5&(48,-96,-32,32,-32,64)\\
6&(48,-96,32,16,16,-64,32).
\end{array}                                          \tag{13}
\]

For the Crabb shift, (9) reads at the two endpoints

\[
(M_2)_{00}=F^{(2)}_{00},\qquad
(M_2)_{LL}=4F^{(2)}_{00}+\sum_{j=1}^LF^{(2)}_{jj}.    \tag{14}
\]

Both (12) and every row of (13) therefore give

\[
\boxed{
[c^2](M_2)_{00}=48,\qquad
[c^2](M_2)_{LL}=128.}                                \tag{15}
\]

The first-amplitude endpoint rows also satisfy

\[
\begin{aligned}
(M_1)_{0j}-2(K_1)_{0j}&=O(c^2),\\
(M_1)_{Lj}-8(K_1)_{Lj}&=O(c^2),\\
(M_1)_{00}=(M_1)_{LL}&=O(c^2).
\end{aligned}                                        \tag{16}
\]

Hence generalized-eigenvector couplings and products of linear
eigenvalue shifts begin at `c^4` and do not affect the coefficient
under audit.  Since `(K_0)_{00}=(K_0)_{LL}=1/2` and the axis
generalized endpoints are `2` and `8`, (15) gives quadratic endpoint
shifts `96c^2` and `256c^2`.  Finally,

\[
[a^2c^2]\frac{\lambda_+}{\lambda_-}
=\frac{256}{2}-\frac{8\cdot96}{2^2}
=128-192=-64,                                       \tag{17}
\]

which proves (1).

## 4. Exact regeneration and scope

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_offset_one_face.py \
  --output experiments/crabb_offset_one_face_s70223.jsonl
```

The checker independently inserts (7) into the guarded exact
ellipse/Stein engine.  It verifies (12)--(17) through size 16, including
all short-chain collisions and nine instances of the stabilized
pattern.  It also reconstructs the exact defect gradient and verifies
(10) through size eight, covering every short collision and the first
separated case.  The finite grid is a regression for the coefficient
proof; the all-size step is the support separation leading to (12),
not the grid itself.

After L134's complex-phase closure, the remaining A84 order is:

1. polarize different coefficient grades;
2. prove unequal-residue localization for nondivisor grades; and
3. prove a dimension-uniform analytic remainder above the Newton face.
