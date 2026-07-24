# The full equality elliptic lift and single-Crabb chart merger

## 1. Results (L190--L192, 2026-07-24)

Fix \(p=L+1\) and one single \(p\)-square Crabb block.  Use L122's
general Hermitian disk chart

\[
 K=H+S^*HS,\qquad A=2K^{-1}HS,\qquad X=K^{1/2}AK^{-1/2},
 \tag{1}
\]

and L187's finite Hardy residual \(\Psi(H)\).  Let

\[
 {\cal E}=\{H:\Psi(H)=0\}
 \tag{2}
\]

after fixing the irrelevant positive scale of \(H\).

The conclusions are:

1. **L190 (full-equality elliptic lift).**  The nonconstant
   characteristic coefficients are real-analytic local coordinates
   on \({\cal E}\).  If

   \[
   \det(\xi I-A)=\xi g_u(\xi),\qquad
   g_u(\xi)=\xi^L+2\sum_{j=1}^{L-1}u_j\xi^j,
   \tag{3}
   \]

   put

   \[
   r=c^L,\qquad w_j=c^ju_j\quad(1\le j<L),\qquad
   {\cal N}=|r|^2+\sum_{j=1}^{L-1}|w_j|^2.
   \tag{4}
   \]

   The reflected model-complement construction extends to every
   \(H\in{\cal E}\), not only the old phase-palindromic Toeplitz
   section.  In a fixed marked neighbourhood it gives

   \[
   \boxed{U_{\rm eq}(u,c)-4\le-a_L{\cal N}}
   \tag{5}
   \]

   for a fixed \(a_L>0\).  Thus every nonzero elliptic perturbation
   of every nearby full-Hardy equality anchor is strict.

2. **L191 (residual/reflection orthogonality).**  At the Crabb
   apex the complete first-residual/reflected initial form of the
   optimized rank-one Stein envelope is

   \[
   \boxed{
   -4\|F\|_F^2-16|\rho|^2
   -64\sum_{j=1}^{L-1}|\omega_j|^2.}
   \tag{6}
   \]

   Here \(F\) is the first nonzero coefficient of \(\Psi(H)\), while
   \(\rho,\omega_j\) are the first coefficients of \(r,w_j\) at the
   same valuation.  In particular there is no \(F\)-\(r\) or
   \(F\)-\(w_j\) mixed row.  Equation (6) is for the prepared
   rank-one **upper** branch; it does not revive A101's false claim
   that the raw characteristic lower function is reflection
   stationary at every disk point.

3. **L192 (complete single-Crabb neighbourhood).**  L189's
   general-disk/circular-normal tube, L190--L191, L120, L160, and
   L163 merge by the associated-graded quadratic form and analytic
   curve selection.  Consequently L118's rank-one Stein certificate
   has condition square at most four in a full operator
   neighbourhood of every fixed single Crabb block:

   \[
   \boxed{t_*(\phi_A(A))\le4.}
   \tag{7}
   \]

   Equality in this neighbourhood occurs exactly on the centered
   circular-range equality anchors \({\cal E}\), modulo affine,
   scalar, and unitary symmetries.

L192 is a local theorem for one fixed block size.  It does not yet
merge repeated Crabb blocks, normal-block collisions, or prove the
global conjecture.

## 2. Characteristic coordinates on the full equality manifold

Two identities hold throughout the disk chart:

\[
 Ae_0=0,\qquad Ae_1=2e_0.                            \tag{8}
\]

Thus the characteristic polynomial is divisible by \(\xi^2\), and
the factor \(g\) in (3) has zero constant coefficient.

At the Crabb point, let \(H=I/2+sZ(z)\) be an arbitrary Hermitian
Toeplitz tangent.  L151's Jacobi calculation gives

\[
 [s]\det(\xi I-A)
 =2\sum_{j=1}^{L-1}\overline {z_j}\,
       \xi^{L-j+1}.                                  \tag{9}
\]

L187 says that the tangent of \({\cal E}\) is exactly the Hermitian
Toeplitz space, modulo its real diagonal scale.  Equation (9) is a
real-linear isomorphism from its \(L-1\) complex physical
coordinates to \((u_1,\ldots,u_{L-1})\).  The analytic inverse
function theorem therefore makes \(u\) a local coordinate system on
\({\cal E}\).

This also resolves the apparent dimension discrepancy in the old
phase-palindromic description.  That family imposed reciprocal
relations among the \(u_j\); the full Hardy equality manifold allows
all \(L-1\) complex characteristic coordinates independently.

## 3. Unrestricted Faber reflection

Let \(P_0=2\), \(P_1=\xi\), and

\[
 P_m(\xi)=\xi P_{m-1}(\xi)-cP_{m-2}(\xi).
 \tag{10}
\]

On the Joukowski boundary \(\xi=\zeta+c/\zeta\),

\[
 P_m(\zeta+c/\zeta)=\zeta^m+c^m\zeta^{-m}.
 \tag{11}
\]

No reciprocal or phase-palindromic relation is used here.  Applying
the Faber transform to (3) gives the exact identity

\[
 \boxed{
 ({\cal F}_cg_u)(\zeta+c/\zeta)
 =g_u(\zeta)+r\zeta^{-L}
  +2\sum_{j=1}^{L-1}w_j\zeta^{-j}.}
 \tag{12}
\]

The positive Hardy factor contains every equality-amplitude
insertion.  Every negative leg lies in the convergent marked ideal

\[
 {\mathfrak r}=(r,w_1,\ldots,w_{L-1},
                   \overline r,\overline w_1,\ldots,\overline w_{L-1}).
 \tag{13}
\]

All finite preparation, polynomial evaluation, orbit-complement,
and Stein operations use sums, products, and inverses of units.
They therefore extend analytically to

\[
 \mathbb R\{u,\overline u,c,{\mathfrak r}\}.
 \tag{14}
\]

## 4. Zero and one reflected legs

For \(H\in{\cal E}\), L187 supplies the canonical Stein metric

\[
 M-A^*MA=qq^*,\qquad q=He_0,
 \tag{15}
\]

with simple generalized endpoints \(h,4h\).  The characteristic
Blaschke product \(B=g/g^\sharp\) is the scalar transfer of the
associated defect-one unitary colligation.  It maps the two endpoint
defect lines and has norm exactly two.

Differentiate one independent negative Hardy leg in (12), including
the total preparation, operator, denominator, and endpoint-basis
motion.  The top norming functional is the invariant version of
L149 equation (7):

\[
 {1\over4}D_\delta\|B(T)\|^2
 =2\operatorname {Re}{1\over2\pi i}
   \int_{\mathbb T}{\dot\Theta(\zeta)\over\Theta(\zeta)}
             {d\zeta\over\zeta}.                    \tag{16}
\]

Here \(\Theta\) is the scalar inner transfer of the complete
colligation, not the raw coefficient numerator alone.  Since
\(\Theta_\delta\) is a differentiable family of finite inner
functions of fixed degree,

\[
 \operatorname {Re}{\dot\Theta\over\Theta}=0
 \quad\hbox{on }\mathbb T.                           \tag{17}
\]

Hence the sharp prepared norm has no one-reflection term at every
full equality anchor:

\[
 R_{\rm eq}-4\in{\mathfrak r}^2.                    \tag{18}
\]

Use L145's orbit-complement vector based on the analytic top right
line of \(B(T)\), and call its condition square \(U_{\rm eq}\).
The model-kernel relation gives

\[
 \Delta=U_{\rm eq}-R_{\rm eq}\ge0.
 \tag{19}
\]

On \({\mathfrak r}=0\), both quantities equal four for every
\(u\).  Therefore the real-analytic nonnegative gap has zero value
and gradient on that divisor:

\[
 \Delta\in{\mathfrak r}^2.                           \tag{20}
\]

This is exactly L149--L150's argument, with the phase-palindromic
companion replaced by L187's invariant defect-one colligation.

## 5. The compact face and L190

At the Crabb apex, different negative Fourier powers in (12) are
orthogonal.  L142--L146's two-reflection calculation therefore gives

\[
 U_{\rm eq}-4
 =-16|r|^2-64\sum_{j=1}^{L-1}|w_j|^2+R_3,
 \tag{21}
\]

where every monomial of \(R_3\) has at least two marked legs and an
additional factor from
\((u,\overline u,c,{\mathfrak r})\).

The calculation is unchanged for a central fold: its real and
imaginary polarizations are taken before the physical phase slice.
On the former phase-palindromic section, the high member of a
reflected pair simply had a larger \(c\)-valuation; it was not a
missing compact coordinate.

Convergence of (14) gives

\[
 |R_3|\le C\delta{\cal N}
 \tag{22}
\]

on a sufficiently small marked polydisk.  Shrinking until
\(C\delta\le8\) proves (5), for example with \(a_L=8\).

## 6. The residual/reflection initial form

The only new mixed block is between L188's first Hardy residual and
one reflected leg.  Work at the Crabb apex and assign them the same
bookkeeping degree \(d\).  Let

\[
 \Psi(H(s))=s^dF+O(s^{d+1}),\quad
 r=s^d\rho+O(s^{d+1}),\quad
 w_j=s^d\omega_j+O(s^{d+1}).                         \tag{23}
\]

Separate the optimized condition coefficient by reflection count.

* The zero-reflection part is L188:

  \[
  [s^{2d}]\Gamma_{\rm zero}=-4\|F\|_F^2.
  \tag{24}
  \]

* The two-reflection part is (21):

  \[
  [s^{2d}]\Gamma_{\rm two}
  =-16|\rho|^2-64\sum_j|\omega_j|^2.
  \tag{25}
  \]

* The one-reflection part is zero.  This is the full-Gram version of
  L156's endpoint first-jet calculation.  Write the raw cleared
  endpoint residual as

  \[
  {\cal C}_H(\xi)
  =\det(\xi K-2HS)\,
    e_0^*K(\xi K-2HS)^{-1}Ke_L
   -\det(K-2\xi S^*H),                               \tag{26}
  \]

  At \(H_0=I/2\), both pencils are triangular.  Their finite
  geometric inverses give \({\cal C}_{H_0}=0\).  The complete raw
  gradient has only the two entries

  \[
  (\nabla_H{\cal C}_{H_0})_{00}
  =(\nabla_H{\cal C}_{H_0})_{L-1,L-1}
  =\xi^{-(L+1)}
  \tag{27}
  \]

  in resolvent normalization.  Every other entry is zero.  The total
  norming functional also differentiates the two defect-line norms.
  Those are exactly the two terms in (27).  Thus the
  endpoint-normalized residual \(\widetilde{\cal C}\) used by L156
  satisfies

  \[
  \widetilde{\cal C}_{H_0}=0,\qquad
  D\widetilde{\cal C}_{H_0}[E]=0                    \tag{28}
  \]

  for every Hermitian Gram-block variation \(E\), without imposing
  an extra diagonal constraint.

  L149/L156 express each total one-reflection norming row as a finite
  coefficient functional of \(\widetilde{\cal C}_H\).  Because both
  the value and full normalized first jet vanish in (28), applying any
  reflected monomial functional gives

  \[
  \boxed{[s^{2d}]\Gamma_{\rm one}=0.}                \tag{29}
  \]

  The model-complement gap cannot add a mixed row: it is
  nonnegative, and L145--L146 make its pure reflected Hessian zero at
  the apex, so positivity of every \(2\times2\) Hessian minor forces
  all residual/reflected entries to vanish.

  This endpoint-normalized formulation is important.  Freezing the
  characteristic numerator or the endpoint lines omits the two
  cancellations in (28) and reproduces A101's false off-equality
  stationarity target.

Equations (24)--(29) prove (6).  They also show why the result is
order independent: \(d\) only locates the first nonzero coefficient;
the colligation calculation is the same quadratic polarization.

For the pure axis leg \(\rho\), L121's stronger flat endpoint
selection gives the same zero mixed row directly.  Thus (29) also
covers the \(u=0\), \(r=c^L\) corner without dividing by an equality
amplitude.

## 7. Adding the true circular normals

Let \(y\) be L173's transported true-circular-normal coordinate.
At one common valuation, the complete quadratic face decomposes as
follows.

1. **\(F\)-\(y\).**  L188 gives the weighted anti-diagonal response
   and L173's flux/null curvature proves a strict Schur complement.
   Only normal modes \(3,\ldots,L-3\) are active.
2. **\(w_1\)-\(y\).**  L160 shows that only the bottom normal
   \(Y_{L0}\) couples.  Completing it leaves the strict residual

   \[
   {32(L-1)(2L^2-L+3)\over
     L(L^2+36L-13)}|w_1|^2.
   \tag{30}
   \]

   This bottom mode is disjoint from L188's active response modes.
3. **\(w_j\)-\(y\), \(j\ge2\).**  L163 proves every eligible leading
   row is exactly zero.
4. **\(r\)-\(y\).**  L120's all-size elliptic-axis strong
   cancellation puts the first possible row strictly above the
   \(2L\) axis face.
5. **\(F\)-\((r,w)\).**  L191 gives zero.

Thus no normal-curvature square is spent twice.  The joint
associated-graded form is strictly negative whenever at least one of

\[
 F,\quad (r,w),\quad y
 \tag{31}
\]

is nonzero.

## 8. Analytic curve selection and L192

The centered radius-one disk matrices form an analytic manifold near
the Crabb block.  L122's general \(H\) chart is a local unitary
cross-section.  L115 supplies a normal splitting into the one
complex elliptic soft coordinate and the \(2p-4\) real coercive
circular normals.  Hence, after affine and unitary normalization,
every nearby operator is represented in the chart

\[
 (H,c,y)\longmapsto
 X(H)+cX(H)^*+N_H(y),                                \tag{32}
\]

followed by its normalized Riemann pullback.

All entries of L118's optimized rank-one envelope are real analytic
in each fixed phase chart.  The reflected lift (14) is convergent,
so it remains valid when the physical relations (4) are substituted.

Assume positive envelope values approach the Crabb point.
Semianalytic curve selection gives a real-analytic positive arc.
Recenter \(H\) on the analytic equality graph from L187 and let

\[
 m=\operatorname {ord}\Psi(H(s)),\quad
 n=\operatorname {ord}(r(s),w(s)),\quad
 q=\operatorname {ord}y(s),
 \tag{33}
\]

with order infinity for an identically zero block.  At the smallest
finite order among \(m,n,q\):

* a single active block has its proved negative pure face;
* two active blocks have the corresponding strict completed face;
* three active blocks have the strictly negative decomposition in
  Section 7.

Every case contradicts positivity.  If all three blocks vanish
identically, then \(c=0\), \(H\in{\cal E}\), and \(y=0\), so the arc
is exact equality rather than positive.

No positive arc exists.  Curve selection proves the full local sign
and hence L192.

## 9. Regeneration and audit boundaries

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_full_equality_elliptic_merger.py \
  --output \
  experiments/crabb_full_equality_elliptic_merger_s70224.jsonl
```

The checker:

1. constructs exact non-palindromic equality jets through order four,
   verifies (9), verifies the zero constant coefficient, and checks
   the unrestricted Faber identity (12);
2. regenerates the exact optimized face (6) for reflected grade one
   through length six and grade two in the first active size,
   including all Riemann and defect jets; and
3. solves genuinely non-Toeplitz nonlinear equality points through
   length six.  Their elliptic descent is phase independent to
   numerical error and its leading coefficient agrees with
   \(64|u_1|^2\).

An additional exact adversarial run over every real zero-Toeplitz
residual basis vector in lengths three and four gives the same zero
mixed coefficient.  Imaginary polarizations follow from the complex
colligation polarization in (26)--(27), not from treating a real-only
sparse checker as a proof.

The theorem is fixed-dimension and local.  It does not assert uniform
constants as \(p\to\infty\), does not settle repeated-block metric
flags, and does not identify a global neighbourhood covering
unrelated extremal types.
