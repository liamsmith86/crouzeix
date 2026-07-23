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

## 5. Corrected merger architecture

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

1. **Dual one-leg stationarity.**
   Prove `R(z,eta)-R(z,0)` has no term linear in `eta` for every
   nearby disk point.  L149 proves this on the equality cone; its
   logarithmic-resolvent derivation must be checked without using the
   equality normalization.
2. **Disk gap is above the normal face.**
   Prove
   `Delta(z,0)=o(Q(z))`, ideally by showing that the characteristic
   Blaschke norm and the orbit-complement metric both have L122's
   quartic `-32Q`.
3. **Reflected gap is above the compact face.**
   L145/L151 already show that the `eta^2` compact diagonal of
   `Delta` is zero.
4. **PSD zero-face polarization.**
   Since `Delta>=0`, zero diagonal on both the disk-normal and
   reflected compact faces should force the mixed face to vanish.
   Only terms strictly above the combined face then remain and are
   absorbed after shrinking the marked polydisk.
5. **Triangular marked coordinates.**
   Prove, rather than infer from an ordinary inverse-function
   theorem, the Rees estimate
   \[
   |r|^2+\sum|w_k|^2
   \asymp c^{2L}+\sum_j|z_j|^2c^{2(L-j)}.
   \]

If these five points close, combine the resulting reflected descent
with L152's canonical disk-normal estimate.  This would prove the
disk-flat tube without any false all-disk primal/dual equality.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/crabb_disk_model_complement.py \
  --output experiments/crabb_disk_model_complement_s70223.jsonl
```

The first JSON line is the exact rational counterexample.  Remaining
lines are finite-amplitude complex grids through the requested
length.  The exact line is the proof artifact; the grid is only an
adversarial regression.
