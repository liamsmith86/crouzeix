# The disk-flat marked merger: obstruction and corrected target (A100, 2026-07-23)

## 1. Status

The first proposed merger theorem was **false as written**.  It
claimed that L145's orbit-complement Stein metric over every nearby
Toeplitz disk point had condition square exactly equal to the norm
square of the characteristic Blaschke product.  That equality holds
on L123's phase-palindromic equality cone and through the compact
faces used by L145/L151, but not at a general disk point.

This note records the exact obstruction and the corrected route.  No
complete disk-flat theorem is claimed here.

## 2. What complementarity really gives

Let `A=A(z)` and `K=K(z)` be L122's coefficient-gauge disk model.  If

\[
\det(\xi I-A)=\xi g_z(\xi),\qquad B_z=g_z/g_z^\sharp,
\]

then Cayley--Hamilton and `ker A=Ce_0` give
`rank B_z(A)=1`.  Let `x,y` be its top right and left generalized
singular vectors and write `B_z(A)x=s y`.

Choose L145's defect

\[
q_z\perp D_z(A)^{-1}
\operatorname{span}\{x,Ax,\ldots,A^{L-1}x\},
\]

and let `P_z` be its rank-one Stein Gramian.  The model-kernel
identity proves

\[
\{P_z-B_z(A)^*P_zB_z(A)\}x=0.
\]

Consequently `x` is an endpoint eigenvector of `(P_z,K)` and

\[
x^*P_zx=s^2y^*P_zy.
\]

This proves only

\[
\kappa_K(P_z)\ge s^2=\|B_z(A)\|_K^2.
\]

Equality would additionally require `y` to be a bottom generalized
eigenvector (and no more extreme interior level).  Rank one and a
scalar Schur complement do not imply that missing assertion.

## 3. Exact rational counterexample

Take `L=3` and the real Toeplitz coefficients

\[
z_1=1/20,\qquad z_2=1/30.
\]

Construct `A,K,B_z,q_z,P_z` exactly over the rationals as above.  The
Blaschke evaluation has rank one, and the model kernel annihilates the
top line exactly.  Nevertheless, for the fixed lower singular line
`y=e_0`, its generalized Rayleigh level is

\[
\mu={2018612479358599772211200\over
4485522806668106801361},
\]

while

\[
P_ze_0-\mu Ke_0
=
\begin{pmatrix}
0\\
-1594407585875686400/4485522806668106801361\\
-6360918118107644000/13456568420004320404083\\
0
\end{pmatrix}
\ne0.
\]

Thus `e_0` is not a generalized eigenvector.  Its Rayleigh quotient
is strictly above the least generalized eigenvalue, so the exact
complementary relation implies

\[
\boxed{\kappa_K(P_z)>\|B_z(A)\|_K^2.}
\]

The nonzero Blaschke norm square in this example is

\[
{1602877538304\over400725714841}.
\]

This disproves the former proposed identity without relying on
floating-point separation.

## 4. Why the failed identity looked exact

The gap is extremely high order at the Crabb apex.  Binary64 samples
of amplitude at most `0.07` separate only around `1e-7` or less.
At amplitudes `0.12,...,0.22`, deterministic general complex samples
show a reproducibly positive gap, growing to visible size, while
`B_z(A)` remains rank one to roundoff.

The earlier checker tested only the small-amplitude regime and used a
tolerance large enough to hide precisely the gap under investigation.
The corrected checker includes the rational obstruction and treats
the finite grid only as supporting scale information.

## 5. The disk model gap is quadratic in `Q`

The exact counterexample does not affect the leading disk-normal
face.  In fact L156 implies the stronger estimate

\[
\boxed{
0\le \Delta(z,0)
=\kappa_K(P_z)-\|B_z(A_z)\|_K^2
\le C_L{\cal Q}(z)^2.
}                                                     \tag{9}
\]

### 5.1 The endpoint residual controls the model orbit

To distinguish the two defects, write

\[
q^{\rm can}=Ke_0=He_0
\]

for L122's canonical disk defect, and write `d` for the normalized
orbit-complement defect.  Put

\[
\chi(\xi)=\xi g(\xi),\quad D(\xi)=g^\sharp(\xi),\quad
T(\xi)=(q^{\rm can})^*
\operatorname{adj}(\xi I-A)e_L,\quad S=T-D.          \tag{10}
\]

L156 says that every coefficient of `S` has value and full gradient
zero on the phase-palindromic equality cone.  L155 therefore gives

\[
\|S\|_{\rm coeff}\le C_L{\cal Q}(z).                 \tag{11}
\]

Let `E` be the unique polynomial of degree at most `L` satisfying

\[
E(\xi)D(\xi)\equiv1\pmod{\chi(\xi)}.
\]

Thus `E(A)=D(A)^{-1}`.  For any polynomial `p`, Laurent expansion of
the resolvent gives the elementary moment formula

\[
(q^{\rm can})^*p(A)e_L
=[\xi^{-1}]_\infty\,{p(\xi)T(\xi)\over\chi(\xi)}.    \tag{12}
\]

Apply (12) with `p=xi^jE`, `0<=j<L`.  Since
`ED=1+chi W` for a polynomial `W`, neither
`xi^j/chi` nor `xi^jW` has a `xi^(-1)` coefficient.  Hence

\[
\boxed{
(q^{\rm can})^*D(A)^{-1}A^je_L
=[\xi^{-1}]_\infty\,{\xi^jE(\xi)S(\xi)\over\chi(\xi)}
=O({\cal Q}(z)).
}                                                     \tag{13}
\]

This makes precise the way L149's endpoint identity generates the
orbit complement: away from equality, its entire failure is a finite
linear functional of `S`.

Let

\[
{\cal O}
=\begin{bmatrix}
D(A)^{-1}e_L&D(A)^{-1}Ae_L&\cdots&
D(A)^{-1}A^{L-1}e_L
\end{bmatrix}.
\]

L123's observability determinant says that `O` has rank `L` on the
equality cone, hence throughout a fixed neighbourhood.  Normalize
the unique vector in `ker O^*` by `d_0=1/2`.  The square system

\[
{\cal O}^*d=0,\qquad e_0^*d=1/2
\]

has a uniformly bounded analytic inverse.  Since
`(q^can)_0=1/2`, equation (13) gives

\[
\boxed{\|d-q^{\rm can}\|\le C_L{\cal Q}(z).}          \tag{14}
\]

### 5.2 The lower endpoint squares the defect error

Let `P` solve

\[
P-A^*PA=dd^*.
\]

The disk model always has `Ae_0=0`, so the normalization in (14)
gives

\[
Pe_0=\frac12d,\qquad
{e_0^*Pe_0\over e_0^*Ke_0}=\frac12,\qquad
(P-\tfrac12K)e_0=\frac12(d-q^{\rm can})=O({\cal Q}).
                                                               \tag{15}
\]

On the equality cone, L123 gives the generalized spectrum
`{1/2,1,...,1,2}`.  The lower level is simple and uniformly separated.
The standard Temple residual bound applied to (15) therefore gives

\[
\lambda_{\min}(P,K)=\frac12-O({\cal Q}^2).            \tag{16}
\]

There is also no hidden singular-vector motion.  For every nearby
Toeplitz disk point,

\[
B(A)=\gamma e_0e_L^*K.                               \tag{17}
\]

Indeed `AB(A)=B(A)A=0`,
`ker A=Ce_0`, and `ker A^*=CKe_L`.  Thus the right and left
generalized singular lines are always `e_L` and `e_0`, and

\[
s^2:=\|B(A)\|_K^2={|\gamma|^2\over4}.
\]

The exact model-kernel complementarity on `e_L`, together with
`e_0^*Pe_0=1/4`, gives

\[
Pe_L={s^2\over2}Ke_L.                                \tag{18}
\]

This remains the largest generalized level by the same uniform
spectral separation.  Combining (16)--(18),

\[
\kappa_K(P)
={{s^2}/2\over 1/2-O({\cal Q}^2)}
=s^2+O({\cal Q}^2).
\]

Nonnegativity of the model gap supplies the lower inequality in (9).
This proves L157 and explains why the false exact identity was so
difficult to distinguish numerically: its error is forced above the
entire quartic disk face.

## 6. Positivity closes the mixed model gap

Work in L149's convergent marked algebra, with the real reflected
coordinate vector

\[
\eta=(r,w,\overline w),\qquad
\|\eta\|^2=|r|^2+\sum_k|w_k|^2.
\]

The variables are independent in a small real marked polydisk; the
physical family is their analytic pullback.  The sharp numerator and
model-orbit constructions remain genuine on this real slice, so the
model-kernel inequality below is nonnegative before pullback, not only
formally coefficientwise.  Let

\[
R(z,\eta)=\|B_{z,\eta}(T_{z,\eta})\|^2,\qquad
U(z,\eta)=\kappa(P_{z,\eta}),\qquad
\Delta=U-R\ge0.                                      \tag{19}
\]

The zero-reflection outer factor is an invertible coordinate
conjugacy, so `R(z,0)` and `U(z,0)` are the disk quantities in
Section 5.  L152 and L157 give

\[
R(z,0)-4\le-a{\cal Q}(z),\qquad
0\le\Delta(z,0)\le C{\cal Q}(z)^2.                  \tag{20}
\]

L154 controls the complete one-reflection row:

\[
\|D_\eta R(z,0)\|\le C{\cal Q}(z).                   \tag{21}
\]

Finally, L144's compact dual face is negative definite.  Analytic
continuity from the Crabb apex and Taylor remainder absorption give,
after shrinking the marked polydisk,

\[
R(z,\eta)-4
\le-a{\cal Q}+C{\cal Q}\|\eta\|-b\|\eta\|^2.         \tag{22}
\]

It remains to show that the nonnegative model gap cannot reintroduce
either mixed face.  The following elementary observation avoids a
new coefficient calculation.  If a nonnegative `C^2` function
`f(t)` obeys `f(0)=a_0` and `|f''|\le M`, then, whenever the minimizing
test step stays in the fixed interval,

\[
|f'(0)|^2\le2Ma_0.                                   \tag{23}
\]

Indeed use Taylor's upper bound at `t=-f'(0)/M`.
Apply (23) on every real line in the `eta` variables.  Equation (20)
then yields

\[
\|D_\eta\Delta(z,0)\|\le C{\cal Q}(z).               \tag{24}
\]

L145--L146 say that the complete compact quadratic face of `U` and
`R` agrees at the apex, including all mixed grades.  Thus
`D_\eta^2\Delta(0,0)=0`; continuity makes its norm at most an
arbitrarily small `epsilon` on a sufficiently small base
neighbourhood.  A second Taylor expansion gives

\[
0\le\Delta(z,\eta)
\le C{\cal Q}^2+C{\cal Q}\|\eta\|
 +\epsilon\|\eta\|^2.                                \tag{25}
\]

The cubic `eta` remainder has been included in the last term after
shrinking.  Since `Q`, `||eta||`, and the base radius all tend to zero,
the first two errors in (25) are absorbed by `-aQ`, and the last by
`-b||eta||^2` in (22).  Therefore

\[
\boxed{
U(z,\eta)-4
\le-a_1{\cal Q}(z)-b_1
\left(|r|^2+\sum_k|w_k|^2\right)\le0
}                                                     \tag{26}
\]

for fixed positive `a_1,b_1` in a sufficiently small marked
neighbourhood.  This proves L158: the disk-normal and reflected faces
merge without exact all-disk complementarity and without a separate
mixed-face recurrence.

## 7. Corrected merger architecture

The false exact identity is stronger than the local theorem needs.
The viable replacement is a graded positive-gap argument.

Let

\[
R(z,\eta)=\|B_{z,\eta}(T_{z,\eta})\|^2,\qquad
U(z,\eta)=\kappa(P_{z,\eta}),\qquad
\Delta=U-R\ge0,
\]

where `eta` denotes the independent reflected Hardy variables
`r=c^L` and `w_k=c^ka_k(z)`.

The required statements are:

1. **Absorbable dual one-leg sector.**
   **Closed by L154--L156.**  Exact general-disk stationarity is false
   (A101), but the endpoint first-jet theorem and analytic cone
   division prove that a coefficient of reflected grade `k` is
   bounded by `C_L|c|^kQ(z)`.  It is absorbed by L152's `-a_LQ`.
2. **Disk gap is above the normal face.**
   **Closed by L157.**  The endpoint residual controls the normalized
   orbit defect by `d-q_can=O(Q)`.  A uniform endpoint spectral gap
   squares that residual and gives `0<=Delta(z,0)<=C_LQ(z)^2`.
3. **Reflected gap is above the compact face.**
   **Closed in L158.**  L145--L146 make the complete compact Hessian
   of `Delta` zero at the apex.
4. **PSD zero-face polarization.**
   **Closed in L158 without an additional recurrence.**  Positivity
   and L157 imply `|D_eta Delta(z,0)|<=C_LQ(z)`;
   continuity of the zero compact Hessian then gives (25), which is
   absorbed by the two negative faces.
5. **Triangular marked coordinates.**
   Prove, rather than infer from an ordinary inverse-function
   theorem, the Rees estimate
   \[
   |r|^2+\sum|w_k|^2
   \asymp c^{2L}+\sum_j|z_j|^2c^{2(L-j)}.
   \]

Items one through four now give the complete disk-flat tube (26) in
the independent marked chart.  The remaining coordinate comparison
is useful for the final pullback and for absorbing L118's strong
variables; it no longer hides a primal/dual equality assumption.

## 8. Regeneration

Run

```bash
.venv/bin/python -u experiments/crabb_disk_model_complement.py \
  --output experiments/crabb_disk_model_complement_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_one_reflection_jets.py \
  --output experiments/crabb_disk_one_reflection_jets_s70223.jsonl
```

The first JSON line is the exact rational counterexample.  Remaining
lines are finite-amplitude complex grids through the requested
length.  The exact line is the proof artifact; the grid is only an
adversarial regression.  The second checker verifies the full real
first jet of (14), in both real and imaginary coefficient directions,
through length ten.
