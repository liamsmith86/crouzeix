# The Toeplitz disk chart and noncoercive quartic (2026-07-23)

## 1. Result and course correction

Put `p=L+1>=3` and let `R` be the `p x p` unweighted
superdiagonal shift.  There is an exact normalized disk chart through
the Crabb matrix with the following form.  For a positive Hermitian
`L x L` matrix `H`, extend it by one zero row and column to `Hhat` and
set

\[
 K=\widehat H+R^*\widehat H R,\qquad
 X(H)=2K^{-1/2}\widehat H R K^{-1/2}.                 \tag{1}
\]

Then

\[
 \boxed{W(X(H))=\overline{\mathbb D}.}                \tag{2}
\]

At `H=I/2`, equation (1) gives `X(H)=C_p`.  The Hermitian Toeplitz
curves

\[
 H(s)=\frac12I+sZ(z),\qquad
 Z(z)_{j,j+k}=z_k,\quad 1\leq k<L,                    \tag{3}
\]

give an exact `2L-2=2p-4` real-dimensional disk submanifold tangent to
all of L115's disk-flat quotient.

Along (3), an explicit rank-one Stein metric has condition-number
expansion

\[
 \boxed{
 \kappa(P(s))
 =4-32s^4{\cal Q}_L(z)+O(s^5),}                       \tag{4}
\]

where

\[
 \boxed{
 {\cal Q}_L(z)
 =\|z\|^4-\left|\sum_{k=1}^{L-1}z_kz_{L-k}\right|^2
 =\|z\|^4-|z^TJz|^2\geq0.}                            \tag{5}
\]

The inequality is just Cauchy--Schwarz.  Its equality set is large:

\[
 {\cal Q}_L(z)=0
 \quad\Longleftrightarrow\quad
 z=\omega J\overline z
 \quad\hbox{for some }|\omega|=1.                     \tag{6}
\]

Therefore the coercive quartic proposed after L121 is **false**.
Already for `p=4`,

\[
 {\cal Q}_3(z_1,z_2)
 =(|z_1|^2-|z_2|^2)^2.                                \tag{7}
\]

Moreover, optimizing the defect vector does not repair (7): the
quartic of L118's optimized rank-one envelope on this exact disk chart
is also

\[
 -32(|z_1|^2-|z_2|^2)^2.                              \tag{8}
\]

Thus the remaining merger must be stratified along (6), rather than
closed by a fictitious bound `-a||d||^4`.

This note proves the disk chart and quartic identities.  It does not
prove the mixed elliptic tube along the new null stratum.  Subsequent
L123 proves that every positive point of (6) is an exact equality
point for `t_*`; see `proof/crabb_palindromic_equality.md`.

## 2. Exact disk factorization

The matrix `K` in (1) is positive definite.  Indeed,

\[
 x^*Kx=x^*\widehat Hx+(Rx)^*\widehat H(Rx),
\]

and the first term kills only the last coordinate while the second
kills only the first.

For `|w|=1`,

\[
\begin{aligned}
K-\{\overline w\,\widehat HR+wR^*\widehat H\}
&=(I-wR^*)\widehat H(I-\overline wR)\\
&\succeq0.                                             \tag{9}
\end{aligned}
\]

Congruence by `K^{-1/2}` proves

\[
 \operatorname{Re}(\overline wX(H))\preceq I.          \tag{10}
\]

Let

\[
 f(w)=(1,w,\ldots,w^L)^T.
\]

Since `(I-\overline wR)f(w)=w^Le_L` and
`\widehat He_L=0`, the right side of (9) annihilates `f(w)`.
Consequently the top eigenvalue in (10) equals one for every `w`.
The support function is identically one, proving (2).

At `H=I/2`, `K=diag(1/2,1,...,1,1/2)`, and direct substitution in
(1) gives the Crabb weights `sqrt(2),1,...,1,sqrt(2)`.

This is the polynomial-support-vector form of Crouzeix's disk
parametrization used by Lewis--Overton.  Here it is retained as a
concrete matrix chart, rather than only as a tangent-manifold theorem.

## 3. Why Toeplitz coordinates are exactly disk-flat

General Hermitian `H`, modulo its irrelevant positive scale, has
`L^2-1=p^2-2p` real parameters.  This is the dimension of the
centered, radius-one disk manifold modulo unitary similarity near
`C_p`, by Lewis--Overton's codimension theorem.  This parameter count
is consistent with (1) being a local unitary cross-section; the
smaller Toeplitz assertion needed below is checked directly rather
than inferred from dimension alone.

Differentiate (1) at `H=I/2`.  For a Toeplitz `Z(z)`, substitution in
L65 equations (12)--(20) gives:

1. every bottom-mode square is zero;
2. the reduced vector in each paired grade lies in the displayed
   low-mode kernel; and
3. the grade-zero part is only the exact disk/scaling orbit.

Hence the L65 quadratic form vanishes.  The `L-1` complex Toeplitz
parameters are independent modulo the affine-unitary orbit, so their
dimension is `2L-2`.  L115 proves that the entire disk-flat quotient
has exactly this dimension.  Therefore (3) supplies all of it.

The nonlinear curve (3), not a least-squares correction of a linear
tangent, is the appropriate disk anchor for the mixed calculation.

## 4. Rank-one metric in coefficient coordinates

Let

\[
 G=K^{1/2},\qquad
 A=G^{-1}X(H)G=2K^{-1}\widehat HR.                    \tag{11}
\]

Set

\[
 q=\widehat He_0
\]

and let `M` be the unique Stein Gramian

\[
 M-A^*MA=qq^*.                                        \tag{12}
\]

Then

\[
 P=G^{-*}MG^{-1}
\]

is a feasible physical contraction metric for `X(H)`.  Its ordinary
eigenvalues are the generalized eigenvalues of `(M,K)`.

At `s=0`, those generalized eigenvalues are

\[
 \frac12,\quad
 \underbrace{1,\ldots,1}_{p-2\ {\rm copies}},\quad
 2.                                                    \tag{13}
\]

The endpoints are simple.  Hence ordinary analytic generalized
eigenvalue perturbation applies.

One useful exact identity persists away from zero.  Since `Ae_0=0`,
the first row of (12) gives

\[
 M_{0j}=q_0\overline{q_j}
       =\frac12K_{0j}.                                 \tag{14}
\]

Thus `1/2` remains an exact generalized eigenvalue.  The lower
endpoint has no Taylor correction.

## 5. The fourth-order calculation

Write

\[
\begin{aligned}
K(s)^{-1}&=\sum_{r\ge0}s^rN_r,
&N_0&=K_0^{-1},
&N_r&=-K_0^{-1}K_1N_{r-1},\\
A(s)&=\sum_{r\ge0}s^rA_r,
&M(s)&=\sum_{r\ge0}s^rM_r.
\end{aligned}                                         \tag{15}
\]

At each order, (12) reduces to

\[
 M_r-A_0^*M_rA_0=F_r,                                 \tag{16}
\]

where `F_r` depends only on lower coefficients.  Because `A_0` is the
nilpotent Crabb coefficient shift,

\[
 M_r=\sum_{j=0}^{L}(A_0^*)^jF_rA_0^j.                 \tag{17}
\]

Equations (15)--(17) are a finite exact recurrence.

Let `lambda_-(s),lambda_+(s)` be the endpoint generalized
eigenvalues.  Substitution of the Toeplitz entries in (15)--(17)
first gives the following endpoint-column identities, valid for every
`L`:

\[
\begin{array}{c|c}
r&\text{endpoint identity}\\ \hline
1&M_1e_L=2K_1e_L\\
2&M_2e_L=0\\
3&M_3e_L=0\\
4&e_L^*M_4e_L=-8{\cal Q}_L(z).
\end{array}                                           \tag{18}
\]

Together with (14), simple generalized-eigenvalue perturbation gives

\[
\begin{aligned}
\lambda_-(s)&=\frac12,\\
\lambda_+(s)&=2-16s^4
 \left\{\left(\sum_{k=1}^{L-1}|z_k|^2\right)^2
 -\left|\sum_{k=1}^{L-1}z_kz_{L-k}\right|^2\right\}
 +O(s^5).                                             \tag{19}
\end{aligned}
\]

Here `K_0e_L=e_L/2`, which changes the `-8` in (18) to
the `-16` in (19).  All coefficients of orders one through three
cancel.  The interior
generalized eigenvalues stay separated from the endpoints, so (18)
and (19) imply (4).

Formula (19) can also be read as a two-row Gram determinant.  The two
rows are `z` and `J conjugate(z)`, which have equal norm.  Their Gram
determinant is precisely (5), proving both its sign and the equality
condition (6).

## 6. Exact failure of coercivity after defect optimization

It is not enough to observe that one convenient metric has a null
quartic: another defect vector could conceivably make L118's envelope
strict there.  In size four this possibility can be eliminated
exactly.

Any analytic minimizing defect vector must agree with
`q=Hhat e_0` to first order; otherwise L118's positive defect-vector
Hessian creates a positive second-order term.  Write its next
coefficient as

\[
 q(s)=\widehat H(s)e_0+s^2(0,u_1,u_2,u_3)^T+O(s^3).
\]

The complete fourth-order condition-number coefficient is

\[
\begin{aligned}
-32(|z_1|^2-|z_2|^2)^2
+16(|u_1|^2+|u_2|^2)+\frac{32}{3}|u_3|^2.             \tag{20}
\end{aligned}
\]

It is minimized uniquely at `u=0`, proving (8).  In particular every
direction with `|z_1|=|z_2|` is genuinely fourth-order flat for the
analytic envelope used by L118.

## 7. The new equality-stratum target

Condition (6) says that the Toeplitz coefficient vector is
phase-conjugate-palindromic.  Exact finite computations at
noninfinitesimal rational points on this stratum give the generalized
metric spectrum (13) unchanged for `p=3,...,8`.  Binary64 SDP and
rank-one optimization also return `t_*=4` to their accuracy.

These observations suggested an exact disk equality family.  L123
subsequently proves it: the metric factors as
`M=K-qq*+2rr*`, and the matching characteristic-polynomial Blaschke
product has norm two.  The union over phase is stratified and singular
at the Crabb point.  The remaining task is therefore to compute and
control the elliptic/non-palindromic normal uniformly through that
singular apex, then merge it with L117/L121.

The former plan of proving a positive-definite quartic and applying
Young globally must not be retried.

## 8. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_disk_toeplitz_quartic.py \
  --output experiments/crabb_disk_toeplitz_quartic_s70223.jsonl
```

The checker:

1. verifies the exact disk factorization (9) on the full Toeplitz
   chart;
2. constructs (15)--(17) symbolically;
3. proves (18)--(19) with every `z_k` symbolic for `p=3,...,6`;
4. checks exact rational directions through `p=10`;
5. verifies the optimized size-four formula (20); and
6. records exact noninfinitesimal palindromic spectra through `p=8`,
   explicitly labeled finite evidence rather than an all-size proof.
