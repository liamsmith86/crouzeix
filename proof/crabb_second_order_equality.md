# Equality-space quotient for the Crabb second variation (2026-07-22)

## 1. Why the raw kernel is misleading

L65 proves that the arbitrary-size single-Crabb second variation is negative semidefinite with
real rank `p(p-2)`.  Its raw kernel therefore has dimension

\[
 2p^2-p(p-2)=p(p+2). \tag{1}
\]

Most of that space consists of motions on which the Riemann-pulled similarity square is exactly
unchanged, not merely flat to second order.  These directions must be removed before any
higher-order calculation.

## 2. Exact affine-unitary orbit

Let `A=C_p`.  The following paths preserve `t_*(phi(A))=4` exactly:

1. unitary conjugation `A -> U^*AU`;
2. complex translation `A -> A+beta I`;
3. nonzero complex scaling/rotation `A -> alpha A`.

Rotation is already contained in the unitary orbit of a weighted shift, so only real scaling
adds a new tangent direction.  The infinitesimal space is

\[
 \mathcal O_p={[A,K]+\beta I+sA: K^*=-K,\ \beta\in\mathbb C,\ s\in\mathbb R}. \tag{2}
\]

The skew-Hermitian commutant of the irreducible shift `A` consists only of scalar multiples of
`iI`.  Hence

\[
 \dim_{\mathbb R}\mathcal O_p=(p^2-1)+2+1=p^2+2. \tag{3}
\]

Every vector in (2) lies in the L65 kernel because it is tangent to an exact constant-value
path.  Quotienting (1) by (3) leaves only

\[
 \boxed{\dim_{\mathbb R}(\ker e_p/\mathcal O_p)=2p-2.} \tag{4}
\]

## 3. Grade-by-grade equality conditions

Use the circle grades and notation `(k,r,t,K_rk)` from L65.

- **Grade zero.** Equality consists of arbitrary infinitesimal weight phases together with one
  common real scaling of the Crabb weights.  This is exactly the diagonal-unitary orbit plus
  scaling; it leaves no quotient direction.
- **Paired mode `k=1`.** `K_r1` has kernel `span(q_r1)`.  Its equality space has dimension
  `2p+2`; the grade-one unitary orbit and translation use `2p` dimensions, leaving two real
  quotient directions.
- **Paired mode `k=2` (`p>=4`).** `K_r2` again has kernel `span(q_r2)`.  Equality has dimension
  `2p`; the unitary orbit uses `2(p-2)`, leaving four real quotient directions.
- **Paired modes `3<=k<=p-2`.** `K_rk` is positive definite, so equality is exactly `t=0`.
  This has dimension `2(p-k+1)`; the grade-`k` unitary orbit uses `2(p-k)`, leaving two real
  quotient directions in each mode.
- **Unpaired mode `k=p-1`.** Equality forces the two entries to sum to zero, exactly the last
  unitary-orbit direction.  The bottom-left mode is forced to vanish.

For `p=3`, the unpaired `k=2` coefficient vanishes and contributes the remaining two quotient
directions.  Thus (4) holds in every dimension.

Equivalently, after exact symmetries are removed, the higher-order problem consists of one
complex direction in mode one, two complex directions in mode two, and one complex direction
in each mode `3,...,p-2` (with the evident `p=3` exception).  This is the correct small target
for third/fourth order.

## 4. Numerical linear-algebra audit

`experiments/crabb_second_order_equality.py` independently reconstructs the entire L62 real
quadratic form and the generator space (2).  For `p=3,...,8` it checks:

- quadratic rank `p(p-2)`;
- orbit rank `p^2+2`;
- quotient dimension `2p-2`;
- `Q_p O_p=0` to solver/Fourier precision.

This audit is not used to infer (3)--(4); it catches indexing, realification, and accidental
overlap errors in the analytic count.

## 5. Next target

L67 chooses an explicit `p=3` mode-one representative and proves strict descent at order six.
L68 derives the exact invariant quartic `−4|w|⁴−31|z|²|w|²/8` on every ray with nonzero mode two.
Thus every fixed canonical straight residual `p=3` ray descends; see
`proof/p3_crabb_higher_order.md`.  A weighted blow-up where `w=O(epsilon*z)` and an exact local
affine-unitary slice with uniform coupling to the negative second-order directions remain before
this is a punctured-neighbourhood theorem.  Higher sizes then require the same calculation on
the grade representatives classified above.
