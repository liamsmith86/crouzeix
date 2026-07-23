# The uniform disk-normal tube (L152, 2026-07-23)

## 1. Result

Let `p=L+1`, let `z=(z_1,...,z_(L-1))` be the Hermitian Toeplitz
coordinates in L122's exact disk chart, and let

\[
{\cal Q}(z)=\|z\|^4-|z^TJz|^2.                       \tag{1}
\]

For the coefficient-gauge disk operator

\[
K=H+R^*HR,\qquad A=2K^{-1}HR,\qquad q=He_0,
\]

let `M` be the rank-one Stein Gramian

\[
M-A^*MA=qq^*
\]

and let

\[
U_D(z)=\kappa_K(M)
\]

be its generalized condition square.  There are constants
`epsilon_L,a_L>0` such that

\[
\boxed{
U_D(z)-4\le-a_L{\cal Q}(z)\le0
\qquad(\|z\|<\epsilon_L).
}                                                     \tag{2}
\]

Equality in (2) holds exactly on L123's phase-palindromic cone

\[
{\cal E}=\{z:z=\omega J\overline z,\ |\omega|=1\}.   \tag{3}
\]

Thus the canonical rank-one metric, not merely the unknown optimum,
gives a uniform disk-normal certificate through the singular Crabb
apex.

## 2. Analyticity and exact stationarity on the equality cone

Near the origin the endpoint generalized eigenvalues of `(M,K)` are
simple.  The lower one is exactly `1/2` by L122 equation (14), so
`U_D` is real analytic.

L123 proves that on every nearby point of (3), the same canonical
metric has generalized spectrum

\[
\frac12,1,\ldots,1,2.
\]

Hence

\[
U_D(u)=4\qquad(u\in{\cal E}).                         \tag{4}
\]

It is important that (4) also has zero first normal derivative.  Fix
an equality point `u` and its characteristic finite Blaschke product
`B_u`.  Every nearby point in the Toeplitz chart still has numerical
range equal to the unit disk.  The disk spectral-set theorem and the
Stein upper certificate give

\[
\|B_u(A(z))\|^2\le t_*(A(z))\le U_D(z).              \tag{5}
\]

At `z=u`, all three quantities in (5) equal four.  The left member has
a local maximum there because it is at most four throughout the disk
chart.  The difference between the right and left members is
nonnegative and vanishes at `u`.  Therefore

\[
\boxed{DU_D(u)=0\qquad(u\in{\cal E}).}                \tag{6}
\]

This upper/lower touching argument is what prevents a linear
disk-normal term.  Merely knowing (4) would not be enough at the
singular union of phase branches.

## 3. The apex jet

L122 computes the all-size fourth-order jet of this same metric:

\[
\boxed{
U_D(z)-4=-32{\cal Q}(z)+O(\|z\|^5).
}                                                     \tag{7}
\]

L124 supplies the best-phase splitting.  After choosing an antiunitary
involution `T`, write

\[
z=u+v,\qquad Tu=u,\quad Tv=-v,\quad\|v\|\le\|u\|.
\]

Then

\[
\boxed{
{\cal Q}(z)=4\|u\|^2\|v\|^2.
}                                                     \tag{8}
\]

The vector `u` lies on one equality branch and `v` is its real
orthogonal normal.  Formula (8) is exact, not just a tangent
approximation.

## 4. Uniform blow-up at the singular apex

Put `F=U_D-4`.  We split the punctured neighbourhood into two regimes.

### 4.1 The normal is visible

Fix a small `delta>0`.  If

\[
\|v\|\ge\delta\|u\|,
\]

then (8) is bounded below by a positive multiple of `||z||^4`, with
the multiple depending only on `delta`.  The fifth-order remainder in
(7) is therefore smaller than `16Q(z)` once `||z||` is sufficiently
small.  Hence

\[
F(z)\le-16{\cal Q}(z)                                \tag{9}
\]

in this regime.

### 4.2 The normal is much smaller than the equality amplitude

Suppose instead that `||v||<delta||u||`.  On the fixed phase branch,
(4) and (6) let us Taylor expand only in the normal:

\[
F(u+v)=\frac12D_v^2F(u)[v,v]+R(u,v).                 \tag{10}
\]

Polarizing the apex jet (7)--(8) gives, uniformly over unit equality
directions,

\[
\frac12D_v^2F(u)[v,v]
=-128\|u\|^2\|v\|^2
+O(\|u\|^3\|v\|^2).                                 \tag{11}
\]

The normal Taylor remainder has at least three normal/equality legs.
Because `F` vanishes identically on both the `u` and `v` phase
branches and its first normal jet vanishes by (6), its potentially
lowest mixed part is bounded by

\[
|R(u,v)|
\le C_L\{\|u\|^3\|v\|^2
          +\|u\|\|v\|^3\}.                          \tag{12}
\]

Choose `delta` so the second term in (12) is at most, say,
`16||u||^2||v||^2`, then shrink `||u||` so the first term has the same
bound.  Equations (10)--(12) yield

\[
F(u+v)\le-96\|u\|^2\|v\|^2
=-24{\cal Q}(z).                                    \tag{13}
\]

The phase parameter ranges over a compact circle and the unit
equality directions form a compact set.  All choices above can
therefore be made uniformly.  Combining (9) and (13) proves (2).

## 5. Why this is the correct disk anchor

The ordinary expansion (7) alone does not prove (2): along paths that
approach (3) unusually quickly, a generic fifth-order remainder could
dominate `Q`.  The missing facts are precisely:

1. exact vanishing on every nonlinear equality branch, and
2. first normal stationarity on every such branch.

Those facts force the remainder estimates (10)--(12) and make the
blow-up uniform through the apex.  This is the disk-normal analogue
of L149--L150's reflected equality tube.

L151 and L152 now give two of the three disk-flat anchors:

\[
-64\sum_j|z_j|^2c^{2(L-j)},\qquad -a_L{\cal Q}(z).
\]

The remaining merger must add L117's `-16c^(2L)` and prove that the
mixed analytic remainder has at least two marked legs.  No further
pure disk-normal calculation is needed.

## 6. Regeneration

The exact inputs are regenerated by

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_toeplitz_quartic.py \
  --output experiments/crabb_disk_toeplitz_quartic_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_equality.py \
  --output experiments/crabb_palindromic_equality_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_normal_form.py \
  --output experiments/crabb_palindromic_normal_form_s70223.jsonl
```

They verify the exact disk chart, the fourth-order jet (7), the
all-size equality metric, and the normal identity (8).  The compact
blow-up in Section 4 is analytic and introduces no numerical premise.
