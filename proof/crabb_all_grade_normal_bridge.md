# All-grade circular-normal bridge: reduced endpoint gate

**Status (2026-07-23): PARTIAL, not a proof of L163.**

## 1. Candidate and proved reduction

Let

\[
2\le k\le\lfloor L/2\rfloor,\qquad
d=k+1,\qquad m=L+2-k.
\]

On the physical grade-`k` reflected path, insert the character-eligible
coercive circular normal `N_m` with weight `d`.  L163 asks whether

\[
\partial_s[\epsilon^{2d}]\Gamma_{L+1}=0.
                                                        \tag{1}
\]

The marked-sector split below rigorously removes the zero-reflection
sector and reduces the one-reflection sector to one associated dual
endpoint functional.  The ordinary singular-value Hessian part of
that functional is proved to vanish in every size.  What remains open
is one finite Hardy/projection identity for the direct mixed
endpoint/cofactor term.  Exact grades two through four and floating
grades through seven satisfy it, but the present notes do not yet
prove it all-size.

## 2. Marked weight separation

Use L149's convergent marked algebra

\[
{\cal A}
=\mathbb R\{u,\overline u,c,r,w,\overline w\},
\qquad
r=c^L,\quad w_j=c^ju_j.                             \tag{2}
\]

For the single grade-`k` amplitude, the first negative Hardy leg is

\[
w_k=c^ku_k,
\]

of physical weight

\[
\operatorname{wt}(w_k)=k+1=d.                       \tag{3}
\]

The reflected partner has weight `L-k+1`, strictly larger unless it
is the central fold; the fold is the same real marked coordinate.
The endpoint leg `r` has weight `L` and can meet the target
`2d` with a weight-`d` strong insertion only in the excluded
grade-one collision.  Two negative Hardy legs have weight at least
`2d`, so after multiplication by the strong insertion their weight is
at least `3d>2d`.

Consequently every term in (1) belongs to exactly one of:

1. the zero-reflection sector;
2. the coefficient linear in `w_k` (or its conjugate/central fold)
   and linear in the first strong coefficient.

There is no two-reflection contribution at weight `2d`.
Convergence of L149's algebra makes this a coefficient identity, not
a formal scaling heuristic.

## 3. The zero-reflection sector

Set the reflected ideal `(r,w,conjugate(w))` to zero.  L149 identifies
the remaining realization with L123's exact disk-equality companion,
up to its invertible prepared outer factor.

L162 proves at **every** such equality anchor that the normalized
rank-one Stein certificate is stationary under every ambient matrix
perturbation:

\[
DU(A)[E-h_E(A)]=0.                                  \tag{4}
\]

The strong normal series, including all of its higher
inverse-Riemann coefficients, is just one analytic ambient
perturbation depending on the benign positive variables in (2).
Substitution in (4) gives an analytic identity in those variables.
Therefore every zero-reflection coefficient of (1) vanishes.

This is where all terms nonlinear in the unreflected equality
amplitude go.  They must not be mixed into the compact
one-reflection calculation.

## 4. Positivity transfers the one-reflection row

Let `U` be L118's optimized Stein envelope and `R` the sharp prepared
Blaschke norm square.  L166 proves, by the PSD Hessian of `U-R`, that
their leading one-reflection/strong rows agree:

\[
\boxed{
[\,\epsilon^{2d}\,]\partial_{w_k}\partial_s U
=
[\,\epsilon^{2d}\,]\partial_{w_k}\partial_s R.
}                                                     \tag{5}
\]

The same statement holds after complex polarization and at a central
fold.  Hence it remains only to calculate the right side.

L145's model-space complement is the compatibility needed here: at
the compact face it constructs the rank-one Stein defect directly
from the prepared Blaschke top singular pair and agrees with the
optimized one-grade defect.  L166's PSD polarization then transfers
the scalar strong row.  Thus the endpoint response used below computes
the common coefficient in (5), rather than importing an unrelated
upper-certificate direction.

By (3), `w_k` and the first strong insertion already have total
weight `2d`.  Any positive-coordinate coefficient or any later
strong coefficient would raise the weight.  Thus the right side of
(5) is evaluated entirely at the Crabb colligation and uses only the
first normalized normal `E_d`.

This observation is the step missing from the earlier
coefficientwise attack: later strong jets are real, but they occur
only in the zero-reflection sector already killed by (4).

## 5. The associated one-reflection coefficient

At the Crabb point the endpoint transfer is

\[
F(z)=e_0^*(zI-C)^{-1}e_L=2z^{-(L+1)}.
\]

L168 proves that the first normalized normal response, after the
matching optimized symmetrized-defect response is removed, is

\[
\boxed{
D\log F
-\gamma_{L,k}D\log F[\operatorname{sym}e_m]
=k(z^{-m}-z^m).
}                                                     \tag{6}
\]

The right side is purely imaginary on the unit circle.

There is an independent polynomial-side formulation.  L169's
characteristic colligation is inner, and L170's first Schur iterate
is the prepared degree-`L` Blaschke factor.  Its grade-`k` tangent is

\[
{\dot B\over B}
=-k(z^{-m}-z^m),\qquad k\ge2.                       \tag{7}
\]

Equations (6)--(7) show that the moving endpoint/defect leg and the
moving prepared numerator are the two sides of the same
linearized colligation.  This is precisely the total endpoint
variation used in L149 equation (7).

The ordinary top-singular-value Hessian contributes no hidden
quadratic term.  L143's first prepared reflected evaluation tangent is

\[
V_k
=-4e_0e_{L-k}^*
+4\sum_{j=k}^Le_je_{j-k}^*
+4e_ke_{L-2k}^*.                                   \tag{8a}
\]

The last term includes the central-fold multiplicity.  Let

\[
W_k=D(C^L)[E_d]
=\sum_{j=0}^{L-1}C^jE_dC^{L-1-j}.
\]

The same weighted-path count as in L168 gives

\[
W_2=\operatorname{diag}(1,2,\ldots,2,1),
\]

and, for `k>=3`, `W_k` lies on superdiagonal `k-2`, with

\[
(W_k)_{j,j+k-2}
=
\begin{cases}
(k-1)\sqrt2,&j=0\ \hbox{or}\ j=L-k+2,\\
k,&\hbox{otherwise}.
\end{cases}                                         \tag{8b}
\]

For `M_0=C^L=2e_0e_L^*`, the reflected Gram tangent
`M_0^*V_k+V_k^*M_0` couples the top right endpoint `e_L` only to
`e_(L-k)`.  The strong Gram tangent
`M_0^*W_k+W_k^*M_0` couples it only to `e_(k-2)`.
These indices cannot coincide:

\[
L-k\ge k>k-2.                                       \tag{8c}
\]

Also `V_ke_L=0`.  Therefore both the direct quadratic term
`e_L^*(V_k^*W_k+W_k^*V_k)e_L` and the simple-eigenvector coupling
through the spectral gap four are zero.  The only mixed contribution
left is the direct `ws` endpoint/cofactor derivative.

L156 rewrites L149's one-reflection row as an analytic finite
coefficient functional of the endpoint-resolvent residual

\[
{\cal R}(z)
=F_{\rm raw}(z)-\Theta(1/z).                        \tag{8}
\]

At the equality point this residual is zero (the harmless endpoint
normalization factor is the same on both sides).  L156 also observes
that differentiating the coefficient functional creates no extra
term when the residual itself vanishes.

L168 gives the raw endpoint logarithmic derivative

\[
\dot F_{\rm raw}/F_{\rm raw}
=kz^{-m}+\beta_{L,k}z^m,                            \tag{9}
\]

whereas L169--L170 give

\[
{D_s\Theta(1/z)\over\Theta(1/z)}
=k(z^{-m}-z^m).                                    \tag{10}
\]

Their difference is the single positive relative mode

\[
(\beta_{L,k}+k)z^m
=\delta_{L,k}\gamma_{L,k}z^m.                      \tag{11}
\]

On the primal side this is exactly the optimized symmetrized-defect
response in L168.  The inner part has zero real boundary mean.
However, L156 proves only that the norming row is a finite analytic
coefficient functional of the residual; it does **not** state its
value on the mode in (11).  Therefore the following kernel identity
is still required:

\[
\boxed{
{\mathfrak F}_{L,k}\!\left[
  (\beta_{L,k}+k)z^m
\right]=0,\qquad
2\le k\le\lfloor L/2\rfloor ,
}                                                     \tag{12}
\]

where `F_(L,k)` is the direct mixed endpoint/cofactor functional
obtained by differentiating L149's one-reflection row in the strong
normal.  A proof of (12) must include the dependence of the prepared
numerator, reversed denominator, operator realization, and endpoint
basis.  Merely observing that `z^m` is a positive relative mode is
not enough.

L149 gives the related logarithmic formula

\[
{1\over4}D_{w_k}R
=2\operatorname{Re}
 {1\over2\pi i}\int_{\mathbb T}
 {\dot\Theta(z)\over\Theta(z)}\,{dz\over z}.         \tag{13}
\]

Its derivation includes numerator preparation, reversed denominator,
operator realization, and outer-coordinate conjugacy for a pure
one-reflection tangent on the equality family.  The missing task is
to differentiate that complete formula in the strong direction and
show that its only extra term is exactly the kernel term (12).  The
inner component then vanishes because

\[
\operatorname{Re}\{k(z^{-m}-z^m)\}=0
\quad(|z|=1).
\]

L170 supplies the analytic Schur continuation needed to perform this
differentiation without selecting a root.  It does not, by itself,
prove (12).

## 6. Grade-one discriminator

For `k=1`, the eligible mode would be `m=L+1`, outside the defect
coordinates.  L167 gives the constant characteristic variation
`dot d=-1`; equivalently the positive mode lies at the first terminal
degree `m=L+1` and is no longer removed by the degree-`L` outer
projection.  L170 correspondingly gives

\[
\dot B=0.
\]

Thus there is no degree-`L` inner tangent matching the raw endpoint
normal.  Equation (6) is unavailable, and the leftover endpoint term
is exactly the nonzero L160/L165 grade-one cross.  Any proof of (12)
must preserve this discriminator.

## 7. Current conclusion and audits

Equation (4) kills the zero-reflection part of (1).  Equation (5)
transfers its sole one-reflection part to the dual, and
(8a)--(8c) kill the ordinary quadratic singular-Hessian terms.
The direct endpoint/cofactor functional (12) is the only remaining
all-size gate in this reduction.

The result agrees with four independent finite mechanisms:

1. grade two is termwise zero (L161/L165);
2. grade three has nonzero adjoint terms but its `a^3c` and `ac^3`
   monomials cancel separately;
3. grade four has six nonzero adjoint terms with zero total;
4. floating complete recurrences through grade seven find the same
   vanishing.

These finite records are adversarial evidence, not a substitute for
(12).  In particular, the exact grade-four run has six nonzero
pairings whose total is zero, so a termwise-support shortcut is false.
L163 and the subsequent uniform tubular merger both remain open.

## 8. Exact regeneration

The sparse associated dual face is regenerated by

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_dual_normal_associated_face.py \
  --minimum-length 4 --maximum-length 14 \
  --output \
  experiments/crabb_dual_normal_associated_face_s70223.jsonl
```

The endpoint, characteristic, and Schur tangents used in
(6)--(11) are regenerated independently by

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_leading_endpoint_transfer.py \
  --minimum-length 4 --maximum-length 14 \
  --output \
  experiments/crabb_leading_endpoint_transfer_s70223.jsonl
```

An additional exact length-eight grade-four run with independent
ellipse and reflected amplitudes gives six nonzero degree pairings
whose total adjoint logarithmic derivative is exactly zero.  This
explicitly falsifies the stronger termwise-vanishing shortcut; only
the structured total cancellation is asserted.
