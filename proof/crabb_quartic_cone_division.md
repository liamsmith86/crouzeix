# Uniform division by the disk equality quartic (L155, 2026-07-23)

## 1. The analytic lemma

Let

\[
{\cal C}z=J\overline z,\qquad
{\cal E}=\{z:z=\omega{\cal C}z,\ |\omega|=1\},
\]

and

\[
{\cal Q}(z)=\|z\|^4-|z^TJz|^2.
\]

Let `F` be a real-valued real-analytic germ at the origin of
`C^(L-1)`.  Suppose:

1. `F(z)=O(||z||^4)`;
2. `F(u)=0` for every `u in E`; and
3. `DF(u)=0` for every `u in E`.

Then there are `epsilon_L,C_L>0` such that

\[
\boxed{
|F(z)|\le C_L{\cal Q}(z)
\qquad(\|z\|<\epsilon_L).
}                                                     \tag{1}
\]

The result is equally valid for a complex-valued analytic germ after
applying it to the real and imaginary parts.

This is a uniform analytic division estimate.  It does not assert
literal divisibility by the single scalar polynomial `Q`; the
codimension of the equality cone is generally larger than one.
Algebraically, the fourth-order face may be any quadratic form in the
`2 x 2` minors of the two-row matrix
`[z; J conjugate(z)]`.

## 2. Best-phase normal coordinates

L124 gives, for every `z`, an antiunitary involution `T` and a
real-orthogonal splitting

\[
z=u+v,\qquad Tu=u,\quad Tv=-v,\quad\|v\|\le\|u\|,
                                                               \tag{2}
\]

such that `u in E` and

\[
\boxed{
{\cal Q}(z)=4\|u\|^2\|v\|^2.
}                                                     \tag{3}
\]

The phase in (2) need not vary analytically with `z`.  The proof below
uses the splitting only pointwise, so this causes no coordinate
problem.

## 3. Fourth-order derivative bounds

The first hypothesis says that every Taylor coefficient of `F` of
degree less than four vanishes at the origin.  After shrinking a fixed
ball, Taylor's theorem therefore gives

\[
\boxed{
\|D^jF(x)\|\le C\|x\|^{4-j},
\qquad j=0,1,2,3,
}                                                     \tag{4}
\]

for all `x` in that ball.

For completeness, write

\[
F(x)=\sum_{d\ge4}F_d(x)
\]

as its convergent homogeneous expansion.  On a smaller ball, the
series obtained by differentiating `j<=3` times converges normally.
Each term is bounded by
`C_d||x||^(d-j)`.  Factoring `||x||^(4-j)` and summing the remaining
geometric majorant proves (4).

## 4. Taylor expansion normal to the cone

Use (2).  The second and third hypotheses give

\[
F(u)=0,\qquad DF(u)[v]=0.
\]

Taylor only in the real normal direction:

\[
F(u+v)
=\frac12D^2F(u)[v,v]
+\frac16D^3F(u+\theta v)[v,v,v]                    \tag{5}
\]

for some `0<theta<1`.  Equations (2) and (4) imply

\[
\begin{aligned}
|D^2F(u)[v,v]|
&\le C\|u\|^2\|v\|^2,\\
|D^3F(u+\theta v)[v,v,v]|
&\le C(\|u\|+\|v\|)\|v\|^3\\
&\le2C\|u\|\|v\|^3
\le2C\|u\|^2\|v\|^2.
\end{aligned}                                       \tag{6}
\]

Combining (3), (5), and (6) proves (1).  No separate angular
compactness or Łojasiewicz exponent is needed.

## 5. Application to reflected coefficients

Let `F_k(z)` be the coefficient of one independent reflected Hardy
variable of grade `k` in any real-analytic Crabb certificate.  If

\[
\boxed{
F_k=O(\|z\|^4),\qquad
F_k|_{\cal E}=0,\qquad
DF_k|_{\cal E}=0,
}                                                     \tag{7}
\]

then L155 gives

\[
|F_k(z)|\le C_{L,k}{\cal Q}(z).                     \tag{8}
\]

There are only finitely many grades in a fixed size.  On the physical
elliptic family, the corresponding variable contains `c^k`.
Consequently

\[
\left|\sum_k c^kF_k(z)\right|
\le C_L|c|{\cal Q}(z)                               \tag{9}
\]

for `|c|<1`.  Shrinking `|c|` absorbs (9) into any fixed fraction of
L152's negative disk-normal term `-a_L Q`.

Thus the corrected disk-flat merger does **not** need exact
one-reflection stationarity at a general disk point.  It needs only
the finite algebraic jets (7).  L148/L149 already give the middle
  condition.
  `experiments/crabb_disk_one_reflection_jets.py` evaluates the full
  root-free recurrence over an exact truncated-series ring.  It
  verifies the last condition at deterministic equality/normal data
  and the first condition at the apex through length ten.  Those
  finite records are adversarial evidence; the two assertions remain
  to be proved all-size.

## 6. Scope

L155 is purely analytic and dimension-by-dimension.  Its constant may
depend on `L`, which is sufficient for the local fixed-size Crabb
classification.  It does not by itself address a dimension-uniform
global Crouzeix constant.

The lemma also applies to primal/dual model gaps or optimized-envelope
coefficients.  Once (7) is established for whichever explicit
certificate is used, no further singular-cone remainder argument is
required.

Regenerate the finite jet audit with

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_one_reflection_jets.py \
  --output experiments/crabb_disk_one_reflection_jets_s70223.jsonl
```
