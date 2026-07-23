# The raw disk-flat elliptic face (2026-07-23)

## 1. Status

This note separates one new all-size endpoint lemma from the remaining
metric statement.

* **Proved below (raw Faber endpoint lemma).** A general Toeplitz disk
  coefficient at offset `j` has reflected ellipse grade `k=L-j`.
  Its characteristic and two Faber endpoint rows are exactly the same
  grade-`k` data that drove L131--L146.
* **Exact finite finding, not yet a lemma.** L118's optimized rank-one
  amplitude Hessian appears to satisfy

  \[
  Q_{L,j}(c)=-64|z_j|^2c^{2(L-j)}
  +o(c^{2(L-j)}).                                    \tag{1}
  \]

  Distinct raw offsets have no mixed coefficient on their first
  possible face.

The remaining proof debt is narrow: extend L145's model-complement
recurrence from the phase-palindromic companion gauge to the raw
coefficient gauge.  The endpoint lemma identifies the same reflected
row, but endpoint agreement alone must not be silently promoted to
metric complementary slackness.

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

Thus L140's unique Hardy correction for this raw direction is
`2a conjugate(zeta)c^kP_j`.  The scalar reflected data are exactly
those that would produce the coefficient `-4` times (10), namely
`-64|\zeta|^2c^(2k)`.

## 5. What remains before (1) is proved

The last sentence of Section 4 is not yet an upper certificate.
For a proof of (1), one must do the following in the raw gauge.

1. Prepare the corrected inner factor

   \[
   P_L+2a\overline\zeta\{P_k+c^kP_j\}
   \]

   and prove that all non-endpoint terms in (3) are inactive below
   reflected weight `k`.
2. Apply L145's orbit-complement construction to this factor and show
   that its condition gap above the inner norm has no coefficient
   through `a^2c^(2k)`.
3. Polarize distinct raw grades.  L125 forces a mixed grade
   `(k,l)` to begin no earlier than `c^(k+l)`; the endpoint rows in
   (7) are orthogonal, but a terminal fold in the metric recurrence
   still has to be checked.

These are coefficientwise triangular recurrences, not a new
optimization problem.  L139 supplies the positive defect Hessian and
L145 supplies the canonical candidate defect; the missing audit is
that the raw-gauge extra rows cannot reach the top generalized
singular pair on the principal face.

## 6. Exact finite audit

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

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_flat_elliptic_face.py \
  --output experiments/crabb_disk_flat_elliptic_face_s70223.jsonl

.venv/bin/python -u experiments/crabb_raw_faber_endpoint.py \
  --output experiments/crabb_raw_faber_endpoint_s70223.jsonl
```

The finite grids guard the algebra and indexing.  They do not replace
the model-complement recurrence in Section 5.
