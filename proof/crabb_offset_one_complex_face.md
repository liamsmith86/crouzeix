# Complex polarization of the universal offset-one face (2026-07-23)

## 1. Result

Let `L>=3` and let the one-pair phase-palindromic coefficient be

\[
u_1=u,\qquad u_{L-1}=\overline u .
\]

Along the corresponding exact disk-equality/ellipse path, L118's
optimized rank-one amplitude Hessian satisfies

\[
\boxed{
[a^2]\Gamma_{L,1}(u;c)
=-64|u|^2c^2+O(c^3).}                                \tag{1}
\]

For `L=2`, phase palindromy forces `u` to be real and L133 already
applies.  Thus (1) is the complete one-grade offset-one coefficient in
every size.

L132 also complexifies: its Dickson reduction and critical-factor
metric lift use no reality assumption on the outer defect.  Hence for
every divisor grade `k|L`,

\[
\boxed{
[a^2]\Gamma_{L,k}(u;c)
\leq -64|u|^2c^{2k}+O(c^{3k})}                       \tag{2}
\]

as a rank-one upper face.  Equality with the full similarity optimum
is still asserted only in the cases already covered by L130; the
upper certificate is what the local Crouzeix merger needs.

This closes complex phase for one divisor grade.  It does not yet
polarize two distinct grades, treat a nondivisor grade, or prove a
uniform remainder.

## 2. Pure-imaginary exact calculation

L133 proves the real direction `u=1`.  Take now `u=i`, so the reversed
coefficient is `-i`.  In L123's companion coordinates write the
coefficient tangent as

\[
E=iF,\qquad F\ \hbox{real}.
\]

Since the ellipse parameter and its Riemann map are real,

\[
S(a,c)=C+cJCJ+ia(F-cJFJ).                            \tag{3}
\]

The coordinate-Gramian tangent is purely imaginary and Hermitian
(its imaginary part is real skew-symmetric).  Conjugation sends
`a` to `-a`; uniqueness of L118's analytic defect branch therefore
makes its first amplitude derivative purely imaginary.

The disk seed and the first elliptic correction are

\[
\boxed{
x_i(c)
=i\{-2e_1+2e_{L-1}+8ce_3\}+O(c^2).}                 \tag{4}
\]

Coordinates are combined when they collide.  Formula (4) is exactly
the general jet

\[
2\overline u\,e_1+2u\,e_{L-1}
-8c\overline u\,e_3
\ O(c^2)                                             \tag{5}
\]

at `u=i`.

Use the same map expansion as L133,

\[
\phi_c(w)=(1+2c^2)w-cw^3+c^2w^5+O(c^3),
\]

but retain conjugate transposes in the Stein equation.  Representing a
Gaussian-rational matrix by its rational real/imaginary pair gives an
exact calculation.  The defect-gradient residual at (4) is `O(c^2)`,
so L118's invertible defect Hessian again implies that replacing the
true optimizer by (4) changes the amplitude Hessian only at `O(c^4)`.

For `L>=7`, the real diagonal of the second-amplitude `c^2` Stein
forcing is the same separated pattern as in the real phase:

\[
(48,-96,-32,16,0,\ldots,0,16,-32,64).                \tag{6}
\]

The short complex chains give

\[
\begin{array}{c|l}
L&\operatorname{diag}F_i^{(2)}\\ \hline
3&(48,-80,-32,48)\\
4&(48,-96,-16,-16,64)\\
5&(48,-96,-32,32,-32,64)\\
6&(48,-96,-96,16,16,0,96).
\end{array}                                          \tag{7}
\]

In every row,

\[
F_{00}^{(2)}=48,\qquad
4F_{00}^{(2)}+\sum_{j=1}^L F_{jj}^{(2)}=128.          \tag{8}
\]

The complex endpoint couplings
`(M_1)_{0j}-2(K_1)_{0j}` and
`(M_1)_{Lj}-8(K_1)_{Lj}` have valuation at least two.
The generalized endpoint calculation is consequently unchanged:

\[
[a^2c^2]\frac{\lambda_+}{\lambda_-}
=128-4(48)=-64.                                     \tag{9}
\]

This proves (1) for `u=i`.

## 3. Conjugation kills the mixed phase

For fixed `L`, let

\[
h_L(u)=[a^2c^2]\Gamma_{L,1}(u;c).
\]

The amplitude Hessian is a real quadratic form in
`u=x+iy`.  Entrywise conjugation sends the equality/ellipse matrix for
`u` to the matrix for `conjugate(u)` and preserves every Stein
condition number.  Hence

\[
h_L(x+iy)=h_L(x-iy).                                 \tag{10}
\]

The mixed `xy` coefficient is odd under (10), so it vanishes.  L133
gives `h_L(1)=-64`; Section 2 gives `h_L(i)=-64`.  Therefore

\[
h_L(x+iy)=-64(x^2+y^2),
\]

which is (1).

## 4. Why divisor descent also complexifies

In L132, replace the real coefficient pair `(1,1)` by
`(u,conjugate(u))`.

1. The Dickson recurrence is linear in these two endpoint
   insertions.  The same first-reversal pairing kills every cross
   residue and leaves the size-`q+1` first-offset pair
   `(u,conjugate(u))`.
2. The Hermitian Toeplitz coordinate Gramian has upper coefficient
   `u` and lower coefficient `conjugate(u)`; its residue-zero
   compression is exactly the corresponding outer Gramian.
3. The critical-factor trace identity is complex-linear in the model
   function.  In the reconstruction
   `q=F(T)^{-*}Ov`, the outer rank-one defect `v` may be complex;
   `sum gamma_jw_j=0` and the metric cross cancellation are unchanged.
4. At the apex the inner metric levels are still two between outer
   endpoints one and four, so the same continuity argument preserves
   the condition number.

Thus L132's upper-envelope lift holds for arbitrary one-pair phase,
and substituting `r=c^k` in (1) proves (2).

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_offset_one_complex_face.py \
  --output experiments/crabb_offset_one_complex_face_s70223.jsonl
```

The checker keeps rational real and imaginary parts separately.  It
verifies the complex Stein equation, Hermitian generalized endpoint
formula, defect stationarity through size seven, every short collision,
and the stabilized forcing pattern through size 12.  No complex
floating-point arithmetic enters the audit.
