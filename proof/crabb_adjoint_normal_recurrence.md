# Adjoint reduction of the circular-normal recurrence

## 1. Result (L165, 2026-07-23)

Let

\[
P-T^*PT=qq^*
\]

be L118's locally optimized rank-one Stein metric, with simple
endpoint eigenpairs

\[
Pv_-=\lambda_-v_-,
\qquad
Pv_+=\lambda_+v_+,
\qquad
\|v_\pm\|=1.
\]

Define

\[
G={v_+v_+^*\over\lambda_+}
  -{v_-v_-^*\over\lambda_-}                          \tag{1}
\]

and let `Z` be the unique solution of the adjoint Stein equation

\[
\boxed{Z-TZT^*=G.}                                   \tag{2}
\]

Then

\[
\boxed{Zq=0}                                         \tag{3}
\]

and, for every operator perturbation `E`,

\[
\boxed{
D\log\kappa(P)[E]
=2\operatorname{Re}\operatorname{tr}(ZT^*PE)
=2\operatorname{Re}\langle PTZ,E\rangle_F .
}                                                     \tag{4}
\]

Thus L163 does not require two differentiated endpoint recurrences.
It is one coefficient of the adjoint matrix

\[
{\cal M}=PTZ.                                        \tag{5}
\]

On a real persymmetric path, L164 additionally gives

\[
\boxed{{\cal M}=J{\cal M}^TJ,}                       \tag{6}
\]

so only one reversed half of its support is independent.

Equations (2)--(6) are exact in every size.  They reduce the
all-grade proof gate but do not themselves prove the vanishing
claimed in L163.

## 2. Derivation of the adjoint formula

For fixed `q`, differentiating the primal Stein equation gives

\[
\dot P-T^*\dot PT
=E^*PT+T^*PE.                                       \tag{7}
\]

The differential of the logarithmic condition number with respect
to the metric is

\[
D_P\log\kappa[H]=\operatorname{tr}(GH).              \tag{8}
\]

The Stein maps

\[
{\cal L}(H)=H-T^*HT,\qquad
{\cal L}^*(Z)=Z-TZT^*
\]

are adjoint for the Frobenius pairing.  Combining (2), (7), and (8)
gives

\[
\begin{aligned}
D\log\kappa[E]
 &=\operatorname{tr}
   Z(E^*PT+T^*PE)\\
 &=2\operatorname{Re}\operatorname{tr}(ZT^*PE),
\end{aligned}
\]

which is (4).

## 3. Why every defect derivative disappears

If the defect also varies by `h`, the right side of (7) acquires

\[
hq^*+qh^*.
\]

Its contribution to (8) is

\[
2\operatorname{Re}h^*Zq.                            \tag{9}
\]

L118 stationarity permits every `h` whose first coordinate is zero,
so all coordinates of `Zq` except possibly the first vanish.
Scaling `q` scales `P` but does not change its condition number.
Taking `h=q` in (9) therefore gives `q^*Zq=0`; since `q_0=1`, the
remaining first coordinate also vanishes.  This proves (3).

This is the adjoint explanation of the normal/defect coefficient
seen in L163.  The defect response is needed to construct the
stationary base metric, but no differentiated defect jet appears in
the final normal derivative.

## 4. Endpoint commutator law

The adjoint gradient has an additional exact displacement identity.
Put

\[
{\cal M}=PTZ.
\]

The primal Stein equation and (3) give

\[
T^*{\cal M}
=T^*PTZ=(P-qq^*)Z=PZ.                               \tag{10}
\]

The adjoint Stein equation gives

\[
{\cal M}T^*
=PTZT^*=P(Z-G).
\]

Subtracting,

\[
\boxed{
T^*{\cal M}-{\cal M}T^*
=PG
=v_+v_+^*-v_-v_-^*.
}                                                     \tag{11}
\]

Thus the interior operator gradient is a commutator flux whose
discrete divergence is supported only on the two endpoint
eigendirections.  This is stronger than persymmetry alone and gives
the natural source of the telescoping in the grade-three audit.

## 5. The weighted L163 coefficient

For reflected grade `k>=1`, put

\[
d=k+1.
\]

Let `T_epsilon` be the flat normalized path, let
`dot T_epsilon` be the normalized derivative of the physical
circular-normal insertion, and expand

\[
{\cal M}_\epsilon=\sum_{j\ge0}\epsilon^j{\cal M}_j,
\qquad
\dot T_\epsilon
=\sum_{j\ge d}\epsilon^jE_j.                         \tag{12}
\]

Equation (4) gives the exact target

\[
\boxed{
\partial_s[\epsilon^{2d}]\log\kappa
=2\operatorname{Re}\sum_{j=0}^{d}
  \langle{\cal M}_j,E_{2d-j}\rangle_F .
}                                                     \tag{13}
\]

This is the promised covariant recurrence.  Compared with the direct
endpoint calculation:

1. the adjoint state is needed only through degree `d`, not `2d`;
2. every defect variation is removed by (3);
3. reciprocal reversal removes the second endpoint through L164; and
4. circle grading can be applied directly to the supports in each
   Frobenius pairing.

The direct inverse-Riemann normal derivative `E_d,...,E_(2d)` is
still present.  Therefore (13) is a strict reduction, not a proof by
notation.

## 6. Persymmetric gradient

Rank-one inverse-Stein duality works before imposing
persymmetry.  It gives an equality of optimized envelopes

\[
\Gamma(T)=\Gamma(JT^TJ).                             \tag{14}
\]

At a fixed point `T=JT^TJ`, differentiation of (14) shows that the
operator gradient is fixed by the same involution.  Since (4)
identifies that gradient with `2PTZ`, equation (6) follows.

This argument avoids trying to differentiate the matrix square roots
in L164.  It also explains the reversed support pairs seen in the
exact series.

## 7. Exact checks and the grade-one discriminator

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_adjoint_normal_gradient.py \
  --output experiments/crabb_adjoint_normal_gradient_s70223.jsonl
```

The checker independently:

1. constructs the physical disk/ellipse/normal path;
2. solves the optimized defect jets;
3. lifts the two endpoint eigenvectors only to form (1);
4. solves (2) as one Stein recurrence;
5. verifies `Zq=0` through the face degree;
6. verifies the endpoint commutator law (11);
7. evaluates every summand in (13); and
8. decomposes the normal variation into commutator and companion
   slices and verifies their pairing balance; and
9. compares (4) with a direct endpoint differentiation.

The last comparison is repeated for grades one and two.  The
grade-three record deliberately stops at the half-order adjoint
calculation; the independent order-eight endpoint equality is already
regenerated by `crabb_grade_three_normal_selection.py`.

The grade-one size-four case gives

\[
D\log\kappa=-{14\over3},
\qquad
D\kappa=-{56\over3},
\]

and has three nonzero terms in (13).  Thus the checker detects L160's
exception rather than being hard-wired to return zero.

For the grade-two size-five face, every individual pairing in (13)
is exactly zero.  This is stronger than cancellation of their sum
and agrees with L161.

Grade three is the first genuinely covariant case.  In size seven,
the five pairings for `j=0,...,4` are

\[
\begin{aligned}
&{ -3249+2500\sqrt2\over12},\quad
 -{ -1481+1224\sqrt2\over12},\quad
 { -173+10\sqrt2\over4},\\
& -{ -1057+514\sqrt2\over12},\quad
 -{ -205+132\sqrt2\over2}.
\end{aligned}                                        \tag{15}
\]

Every term is nonzero, but their exact sum is zero.  This independently
recovers the grade-three selection using an adjoint state through
degree four, rather than two endpoint states through degree eight.
It also rules out a termwise support-disjointness proof in general:
the all-grade completion must exhibit a true telescoping identity.

The half-order solver also makes the first central grade-four case
practical.  In size nine its six adjoint pairings are all nonzero and
sum exactly to zero.  The companion split gives

\[
B_4=2(-607+439\sqrt2),\qquad
H_4=-2(-607+439\sqrt2).                              \tag{16}
\]

This is a further exact case, not an all-grade extrapolation.

## 8. Companion normal form of the telescope

There is a canonical way to separate coordinate motion from
characteristic-polynomial motion.  Write `C=T_0` for the Crabb shift.
Every matrix `R` has a unique decomposition

\[
R=[C,X]+H,                                           \tag{17}
\]

where the first row of `X` is zero and `H` is supported on the bottom
row.  If `a_i=C_(i,i+1)`, the decomposition is the elementary
recurrence

\[
\begin{aligned}
X_{i+1,j}
 &={R_{ij}+{\bf1}_{j>0}a_{j-1}X_{i,j-1}\over a_i},
 &&0\le i<L,\\
H_{L,j}
 &=R_{L,j}+{\bf1}_{j>0}a_{j-1}X_{L,j-1}.            \tag{18}
\end{aligned}
\]

Applying (17) recursively to the whole series gives

\[
\dot T=[T,X]+H.                                     \tag{19}
\]

The commutator law (11) then turns (13) into

\[
{1\over2}D\log\kappa[\dot T]
=\operatorname{Re}\operatorname{tr}
 \{(v_+v_+^*-v_-v_-^*)X\}
+\operatorname{Re}\langle{\cal M},H\rangle_F.       \tag{20}
\]

Thus the desired telescope has two sharply identified pieces:

1. an endpoint boundary term caused by the change of companion basis;
2. a bottom-row term carrying the change of characteristic
   polynomial.

The exact checker verifies (19)--(20) coefficientwise.  In grade two,
both pieces vanish separately.  In grade three their target
coefficients are nonzero:

\[
\begin{aligned}
B_3&=-{-1060+683\sqrt2\over12},\\
H_3&= { -1060+683\sqrt2\over12},
\end{aligned}                                       \tag{21}
\]

and hence cancel exactly.  Grade one is the required discriminator:

\[
B_1={-13+9\sqrt2\over3},\qquad
H_1=2-3\sqrt2,\qquad B_1+H_1=-{7\over3}.             \tag{22}
\]

Keeping the equality amplitude `a` and ellipse amplitude `c`
independent sharpens (21) to

\[
\begin{aligned}
B_3(a,c)
&=-{ac\over12}
 \{(-744+531\sqrt2)a^2+(-316+152\sqrt2)c^2\},\\
H_3(a,c)&=-B_3(a,c).                                \tag{23}
\end{aligned}
\]

Thus the `a^3c` zero-reflection transport and the `ac^3`
one-reflection term cancel separately; the zero at `a=c=1` is not an
accidental cancellation between those two mechanisms.  This exact
bivariate split supports treating the first monomial through L162 and
the second through L166 plus L149.

### The leading characteristic monomial

There is nevertheless a useful all-size statement at the first strong
weight.  Let `C=C_p`, let `E_d` be the leading normalized variation,
and write

\[
\chi_C(z)=\det(zI-C)=z^{L+1}.
\]

Then

\[
\boxed{D\chi_C[E_d](z)=-kz^{k-1}.}                  \tag{24}
\]

Indeed

\[
\operatorname{adj}(zI-C)
=\sum_{r=0}^L z^{L-r}C^r.
\]

The inverse-Riemann correction is a polynomial in `C` with zero
constant term, so it contributes no trace.  In the Riesz
representative `N_(L+2-k)`, only its subdiagonal of offset
`L+1-k` can close a trace with a power of `C`.  There are exactly `k`
such entries.  The endpoint factors `sqrt(2)` in `C` cancel the
`1/sqrt(2)` endpoint factors in `N`, so each contributes one.  Jacobi
differentiation gives (24).

Equivalently, the leading companion slice in (18) is

\[
H_d=
\begin{cases}
\frac12e_Le_0^T,&k=1,\\[1mm]
\frac{k}{\sqrt2}e_Le_{k-1}^T,&k\ge2.
\end{cases}                                         \tag{25}
\]

The right-end product of Crabb weights is `2` in the first line and
`sqrt(2)` in the second, so (25) is equivalent to (24).

Thus the leading grade-`k` characteristic slice retains a factor `z`
exactly when `k>=2`; grade one has the constant variation `-1`.  This
is a genuine discriminator on L166's associated compact face and may
be used in a leading inner/companion proof.

### A full fixed-zero shortcut is false

Grade two happens to satisfy

\[
\partial_s\det T_{\epsilon,s}=0
\quad\hbox{through weight six}.
\]

This does not extend to higher grades.  In the central grade-three
case the exact determinant derivative has nonzero coefficients

\[
[\epsilon^5]=-4,\qquad
[\epsilon^7]=-2(-17+8\sqrt2),\qquad
[\epsilon^8]=-4(-17+12\sqrt2).                      \tag{26}
\]

The condition derivative nevertheless vanishes at weight eight.
Therefore the boundary/characteristic balance cannot be proved by
asserting that the normalized normal variation preserves a zero
eigenvalue or a fixed factor of the characteristic polynomial.
The fixed-degree inner argument, if successful, must include the
moving characteristic factor rather than freeze it.

Regenerate this adversarial check with

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_companion_fixed_zero_adversary.py \
  --output \
  experiments/crabb_companion_fixed_zero_adversary_s70223.jsonl
```

Consequently endpoint flux alone is not the missing proof.  The
remaining theorem must identify the bottom-row characteristic term
with the negative endpoint-basis term for every grade `k>=2`.
Equations (21)--(23) strongly point toward the companion/inner
transfer mechanisms of L149 and L162, but no all-grade identification
with those mechanisms is claimed here.

L168 subsequently identifies the first such balance explicitly:
normal minus optimized symmetrized-defect response is
`k(z^(-m)-z^m)`.  L169 then proves that the complete optimized
rank-one Stein system has an exact scalar inner colligation transfer,
whose feedthrough contains the later characteristic motion detected
in (26).  L171 removes the zero-reflection sector and the ordinary
quadratic dual terms; L172 closes the direct differentiated endpoint
functional by Cauchy-residue selection.  The convolution below remains
an independent physical-coordinate formulation.

## 9. Remaining all-grade statement

**Resolved downstream by L172.**  The formulation below records the
exact physical-coordinate version of the endpoint/characteristic
balance.

L163 is now equivalent to

\[
\boxed{
\sum_{j=0}^{k+1}
\langle{\cal M}_j,E_{2k+2-j}\rangle_F=0,
\qquad k\ge2,
}                                                     \tag{27}
\]

for the sole character-eligible circular normal
`m=L+2-k`.

A coefficientwise completion would derive the circle-graded support of
`M_j` from (2)--(3), transport the direct normal coefficients
`E_(2d-j)` through the inverse-Riemann recurrence, and prove the
boundary/characteristic balance in (19).  It no longer needs to
propagate two endpoint eigenvectors or a free final defect jet through
weight `2d`.

L171 instead splits the convolution in L149's marked algebra.  L162
kills its zero-reflection sector.  In the one-reflection sector only
the first strong coefficient can meet the target weight; its
quadratic singular-Hessian supports are disjoint.  The remaining
direct endpoint term is L156/L168's differentiated kernel identity.
L172 writes it as a Cauchy residue and proves it zero, establishing
(27) without expanding every `M_j`; grade one remains outside exactly
as required.
