# The L70 weighted center is an exact disk-matrix curve (2026-07-22)

## 1. The analytic curve

Let `e` be real and put

\[
D(e,R)=
\begin{bmatrix}
-e/3&\sqrt2&e\\
(3\sqrt2/64)e^2R&2e/3&\sqrt2\\
-(9/64)e^3R&(3\sqrt2/64)e^2R&-e/3
\end{bmatrix}.                                                \tag{1}
\]

Define

\[
F(e,R)=81e^4R^2+(1152e^2-4096)R+4096.                         \tag{2}
\]

There is one branch of `F=0` analytic at `e=0` with `R(0)=1`:

\[
R(e)=\frac{64\{32-9e^2-8\sqrt{16-9e^2}\}}{81e^4}
=1+\frac9{32}e^2+\frac{405}{4096}e^4
+\frac{5103}{131072}e^6+O(e^8).                              \tag{3}
\]

The apparent singularity in (3) is removable.  In the L70 chart this is

\[
u(e)=\frac{3\sqrt2}{64}R(e),\qquad
v(e)=-\frac9{64}R(e),\qquad s(e)=0.                           \tag{4}
\]

In particular, the order-two correction in (3) gives exactly
`(27 sqrt(2)/2048,-81/2048)`, the center shift independently extracted from
the order-eight gradient in the higher-jet certificate.

## 2. Characteristic and Kippenhahn factorizations

The discriminant of the characteristic polynomial of (1) is

\[
-\frac9{2^{22}}Re^6(3R+64e^2-128)F(e,R).                     \tag{5}
\]

On the branch (2), let `c(e)` denote the repeated eigenvalue.  Exact division
gives

\[
\det(\lambda I-D(e,R))=(\lambda-c)^2(\lambda+2c).             \tag{6}
\]

Now form the homogeneous Kippenhahn polynomial

\[
p(x,y,z)=\det\left(zI+x\operatorname{Re}D+y\operatorname{Im}D\right).
\]

Writing `r(e)^2` for an explicit positive rational expression in `(e,R)`, the
same exact reduction modulo `F` yields

\[
\boxed{p(x,y,z)=(z-2cx)\bigl((z+cx)^2-r(e)^2(x^2+y^2)\bigr).} \tag{7}
\]

The quadratic factor is the dual equation of the circle centered at `c` with
radius `r(e)`.  The linear factor is the point `-2c`.  At `e=0`, `c=0` and
`r=1`; hence, by continuity, the point remains strictly inside the circle for
all sufficiently small `e`.  Kippenhahn's theorem therefore gives

\[
\boxed{W(D(e,R(e)))=\{z:|z-c(e)|\le r(e)\}.}                 \tag{8}
\]

Thus the repeated eigenvalue in (6) is exactly the center of the numerical-
range disk.

## 3. Complete similarity bound on the center

Normalize

\[
T_e=\frac{D(e,R(e))-c(e)I}{r(e)}.
\]

Equation (8) says `w(T_e)=1`.  The classical Berger--Okubo--Ando theorem says
that every numerical-radius contraction is similar to a contraction through
a similarity of condition number at most two.  Equivalently, the unit disk is
a complete `2`-spectral set for `T_e`.  Since the Riemann map of (8) is the
displayed affine normalization up to a disk automorphism, L21 gives

\[
\boxed{t_*(\phi_e(D(e,R(e))))\le4.}                           \tag{9}
\]

This is an exact theorem for the entire analytic center curve, not a finite-
order inference.  It explains why the rank-one formal certificate repeatedly
recentered and remained flat.

## 4. Scope

L70 proves strict leading feasible-certificate descent in every bounded
weighted direction away from (4); (9) closes (4) itself.  These two statements
do not automatically prove a full punctured neighbourhood: a sequence can
approach the disk curve faster than the leading weighted scale.  A uniform
normal estimate, or an analytic weighted Morse--Bott argument for the feasible
certificate, is still required before making that claim.

The disk theorem is classical, not new.  The apparently new content is the
exact identification of L70's hidden center with the disk-matrix locus in the
chosen affine-unitary slice.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/p3_crabb_disk_center.py
```

The script derives (2), (5)--(7), the double eigenvalue, and the series (3)
using exact SymPy polynomial remainders.  It makes no floating-point sign or
factorization decision.
