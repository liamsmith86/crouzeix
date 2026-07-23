# The raw disk-flat elliptic face (L151, 2026-07-23)

## 1. Status

Let
`Q_L(z;c)=[a^2]Gamma_L(az,c)` be the quadratic coefficient of L118's
optimized rank-one envelope at the Crabb axis in a general Toeplitz
disk-flat direction `z=(z_1,...,z_(L-1))`.  Put

\[
{\cal E}_c(z)=\sum_{j=1}^{L-1}|z_j|^2c^{2(L-j)}.
\]

Then the complete raw principal face is

\[
\boxed{
Q_L(z;c)=-64{\cal E}_c(z)
+\text{terms strictly above the diagonal face}.
}                                                     \tag{1}
\]

More precisely, with
`D_c=diag(c^(L-1),...,c)`, the Hermitian Hessian matrix satisfies

\[
D_c^{-1}Q_L(c)D_c^{-1}\longrightarrow-64I
\qquad(c\longrightarrow0).                           \tag{1a}
\]

Consequently, for every fixed `L`, it is negative definite for all
sufficiently small nonzero `c`, with for example
`Q_L(z;c)<=-32 E_c(z)`.

The proof has three parts: the raw Faber endpoint lemma, a corrected
raw Blaschke factor and its sparse preparation, and L145's canonical
model-complement upper metric.

## 2. Full coefficient gauge

Let `p=L+1`, let `R` be the forward shift, and put

\[
 H(a)=H_0+aZ_j(\zeta),\qquad H_0=\frac12I_L\oplus0,
\]

where `Z_j` is Hermitian Toeplitz with upper `j`th diagonal equal to
`\zeta`.  In coefficient coordinates,

\[
\begin{aligned}
K(a)&=H(a)+R^*H(a)R,\\
S(a,c)&=2K(a)^{-1}\{H(a)R+cR^*H(a)\}.               \tag{2}
\end{aligned}
\]

Write `S_0=S(0,c)` and `Y_j=\partial_aS(0,c)`.  Direct multiplication
by the diagonal matrix `K(0)` gives, with `k=L-j`,

\[
\boxed{
Y_j=
2\zeta e_0(e_{j+1}^*-ce_{j-1}^*)
-2\overline\zeta e_L(e_{k+1}^*-ce_{k-1}^*).
}                                                     \tag{3}
\]

Indices in (3) lie between zero and `L` because `1<=j,k<L`.
At `c=0`,

\[
\dot A_j
=2\zeta e_0e_{j+1}^*
-2\overline\zeta e_Le_{k+1}^*.                      \tag{4}
\]

This is why the phase-palindromic companion shortcut cannot be used
for a raw coefficient: it would add the corresponding two terms from
the reversed offset.

## 3. Raw characteristic grade

The Crabb matrix `C=S_0|_(c=0)` is nilpotent of order `L+1`.
Jacobi's determinant formula and

\[
(\xi I-C)^{-1}=\sum_{m=0}^{L}\frac{C^m}{\xi^{m+1}}
\]

give

\[
\boxed{
\left.\partial_a\det(\xi I-S(a,0))\right|_{a=0}
=2\overline\zeta\,\xi^{k+1}.
}                                                     \tag{5}
\]

Indeed, the first term in (4) has zero trace against every power of
`C`.  The second term contributes only for `m=j-1`, when
`C^(j-1)e_L=e_(k+1)`.  Thus the raw offset `j`, not its smaller index,
is a characteristic perturbation of reflected grade `k=L-j`.

## 4. Exact Faber endpoint lemma

Let

\[
P_0=2,\qquad P_1(\xi)=\xi,\qquad
P_m(\xi)=\xi P_{m-1}(\xi)-cP_{m-2}(\xi)              \tag{6}
\]

be the Dickson polynomials.  Then

\[
\boxed{
\begin{aligned}
e_0^*\{DP_L(S_0)[Y_j]
       +2\overline\zeta P_k(S_0)\}
 &=4\overline\zeta e_k^*,\\
e_L^*\{DP_L(S_0)[Y_j]
       +2\overline\zeta P_k(S_0)\}
 &=4\overline\zeta c^k e_j^* .
\end{aligned}}                                       \tag{7}
\]

This is a polynomial identity in `c`.

To prove it, set `D_0=0`, `D_1=Y_j`, and differentiate (6):

\[
D_m=Y_jP_{m-1}(S_0)+S_0D_{m-1}-cD_{m-2}.            \tag{8}
\]

Use L131's endpoint paths

\[
e_0^*P_m(S_0)=2e_m^*,\qquad
e_L^*P_m(S_0)=2c^me_{L-m}^*.                         \tag{9}
\]

For `1<=s<L`, define the two path rows

\[
\begin{aligned}
A_{m,s}
&=\begin{cases}
e_{s+m}^*,&s+m\le L,\\
c^{s+m-L}e_{2L-s-m}^*,&s+m>L,
\end{cases}\\
B_{m,s}
&=\begin{cases}
c^me_{s-m}^*,&m\le s,\\
c^se_{m-s}^*,&m>s.
\end{cases}                                         \tag{10}
\end{aligned}
\]

Substitution of (3) into (8), simultaneously at the two endpoints,
gives the stronger identities

\[
\boxed{
\frac12e_0^*D_m=\zeta(A_{m,j}-B_{m,j}),\qquad
\frac12e_L^*D_m
=-\overline\zeta(A_{m,k}-B_{m,k}).
}                                                     \tag{11}
\]

They hold at `m=1` by (3).  The recurrence (8), together with (9),
moves the forward path in `A` and the reverse path in `B` by one;
the second clauses in (10) are exactly the two endpoint folds.
This proves (11) by induction.  At `m=L`, `A_(L,s)=B_(L,s)` for every
`s`, so both rows of `D_L` vanish.  Adding the two explicit
`P_k(S_0)` rows from (9) proves (7).  No root selection or asymptotic
argument enters.

The negative Hardy row in (7) has squared norm

\[
16|\zeta|^2c^{2k}.                                  \tag{12}
\]

The endpoint rows determine the grade but not, by themselves, the
sharp scalar factor.  The raw factor is derived next.

## 5. Corrected raw Blaschke factor

The factor that simultaneously works for `j<k`, `j=k`, and `j>k` is

\[
\boxed{
G_{j,a,c}=P_L+2a\overline\zeta(1+c^j)P_k,
\qquad k=L-j.
}                                                     \tag{13}
\]

The correction is `c^jP_k`, not `c^kP_j`.  On the low-reflected side
`j>k` it lies above the face and changes nothing; at a central
collision it is essential; on the high-reflected side `j<k` it
cancels every premature folded loss.  The sum of the two raw factors
for offsets `j` and `k` reproduces L140's equality-pair correction,
up to a term strictly above the smaller face.

Compose (13) with `Psi_c`, prepare its monic degree-`L` numerator, and
call the corresponding finite Blaschke product `B=N/N^sharp`.
Evaluate it on the full operator (2) and write

\[
B(T)=B_0+aB_1+a^2B_2+O(a^3).
\]

Weierstrass preparation is triangular in amplitude degree, scalar
power, and `c`-degree.  Applying its coefficient recurrence to (3)
and (13) gives, for `0<=d<k`,

\[
\boxed{
[c^d]B_1
=4E_{d,k-d}
-4{\bf1}_{d\ge j}E_{d-j,k-d+j},
}                                                     \tag{14}
\]

and at the terminal grade

\[
\boxed{
[c^k]B_1
=4\sum_{r=k}^{L}E_{r,r-k}
-4{\bf1}_{k\ge j}E_{k-j,j}.
}                                                     \tag{15}
\]

Here `E_(r,s)=e_re_s^*`.  To verify the induction, the first term in
(14) advances one step southeast under the Dickson recurrence.  Once
`d=j`, the `c^jP_k` correction starts the identical path with opposite
sign.  The two paths telescope until (15); the remaining translated
shift is the reflected endpoint row (7).  These are all terms of
weight at most `k`, because every nonleading coefficient of `Psi_c`
pairs with its inverse coefficient from `phi_c`, exactly as in L143.

In particular,

\[
\operatorname{val}_c(B_1e_L)>2k,\qquad
\operatorname{val}_c\{(H_1-4K_1)e_L\}>2k,            \tag{16}
\]

where `H_1` is the first amplitude derivative of `B^*KB`.  Continuing
the same triangular recurrence one more amplitude degree gives

\[
\boxed{
[c^d]e_0^*B_2e_L=0\ (d<2k),\qquad
[c^{2k}]e_0^*B_2e_L=-16.
}                                                     \tag{17}
\]

There is only one two-insertion path to the top row.  If it folds
before grade `k`, its copy coming from `c^jP_k` has the opposite sign;
the sole unpaired path is the terminal reflected path, whose two
endpoint weights multiply to `-16`.  This also proves (17) when
`j=k`; the two coincident paths combine before the terminal step.

At the axis,
`B_0e_L=2e_0` and `K_(0,00)=K_(0,LL)=1/2`.  Equations
(15)--(17) show that the first singular-vector Schur correction is
strictly above `2k`.  The two `B_0`--`B_2` cross terms in the top Gram
entry give `-32c^(2k)`; division by `K_(0,LL)=1/2` gives

\[
\boxed{
\|B(T)\|_K^2
=4-64|\zeta|^2a^2c^{2k}
+o(a^2c^{2k}).
}                                                     \tag{18}
\]

## 6. The raw model-complement lift

Write `B=N/D`, let `x` be its simple top right singular vector, and
choose L145's canonical defect

\[
q\perp D(T)^{-1}
\operatorname{span}\{x,Tx,\ldots,T^{L-1}x\}.         \tag{19}
\]

The orbit is triangular and has rank `L` at the axis, so (19) defines
an analytic line.  The model-kernel identity gives exact
complementarity on `x`:

\[
x^*P_T(q)x=\|B(T)\|^2y^*P_T(q)y.                    \tag{20}
\]

For the coefficient assertion, differentiate the `L` equations
(19).  Equations (14)--(16) say that every forward orbit leg below
grade `k` is paired with its shifted negative copy.  Solving the
triangular system therefore makes the two endpoint generalized
eigenvector residuals have valuation greater than `k`.  Their Schur
squares have valuation greater than `2k`.  At grade `2k`, (20) and
(17) leave the same terminal Rayleigh coefficient on the two extreme
eigenvalues.  Hence

\[
\boxed{
[a^2c^d]\{\kappa(P_T(q))-\|B(T)\|^2\}=0
\quad(0\le d\le2k).
}                                                     \tag{21}
\]

This is the raw-gauge extension of L145 equation (44).  It supplies an
explicit feasible rank-one Stein metric attaining (18), so the
optimized envelope has the same one-grade face.

## 7. Complex phases and distinct raw grades

Equations (3), (5), and (13) are complex linear in
`conjugate(zeta)`.  Polarizing the preparation over `zeta` and its
conjugate replaces `zeta^2` by `|zeta|^2`; nonzero circle-grade
monomials cannot reach the endpoint functional.  Thus (18)--(21) are
phase isotropic on the principal face.

For two raw offsets with distinct reflected grades `k` and `l`, L125
makes every mixed coefficient divisible by `c^(k+l)`.  At that grade,
the two translated shifts in (15) occupy distinct Fourier diagonals.
Their top Gram product is zero.  If one path folds, its negative copy
from the corresponding `c^jP_k` term cancels it before the endpoint.
Therefore the corrected dual face has no mixed coefficient.

Finally, the optimized primal Hessian minus the dual Hessian is
positive semidefinite for every fixed `c`.  Scale by
`diag(c^k)` over all raw grades.  Equation (21) makes every diagonal
entry of the limiting difference zero; positivity makes every
off-diagonal entry zero, exactly as in L146.  This proves (1)--(1a).

## 8. Exact finite audit

`experiments/crabb_disk_flat_elliptic_face.py` starts from (2), expands
`K(a)^(-1)` before applying the exact ellipse map, and eliminates the
rank-one defect tangent over rational truncated series.  It does not
use the phase-palindromic companion shortcut.

The persisted grid contains:

* 15 single raw offsets in lengths `2,...,8`, with reflected grades
  through three and several offsets on both sides of the midpoint;
* coefficient `-64` at every `c^(2k)` and exact zero below it;
* three polarized pairs in length five, all exactly zero through
  their first possible mixed face; and
* five phase-palindromic bridge cases where the full gauge and the
  independently derived companion gauge agree coefficientwise.

`experiments/crabb_raw_faber_endpoint.py` independently regenerates
(5) and (7) for both real and imaginary coefficient phases.

`experiments/crabb_raw_faber_blaschke.py` independently prepares
(13), verifies the complete sparse recurrence (14)--(15), the
quadratic entry (17), the dual coefficient (18), the model
orthogonality equations (19), and primal/dual equality (21).  Its
cases include every offset on both sides of the midpoint through
length five and an additional length-six high-reflected case.

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_flat_elliptic_face.py \
  --output experiments/crabb_disk_flat_elliptic_face_s70223.jsonl

.venv/bin/python -u experiments/crabb_raw_faber_endpoint.py \
  --output experiments/crabb_raw_faber_endpoint_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_raw_faber_blaschke.py \
  --output experiments/crabb_raw_faber_blaschke_s70223.jsonl
```

The three engines use different algebraic representations: direct
optimized Stein elimination, raw Dickson endpoint recurrence, and
prepared Blaschke/model complement.  Their agreement guards against a
companion-gauge or endpoint-normalization mistake.
