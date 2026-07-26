# Lower-neutral preparation does not make the quintic a global response

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

## 1. Result (L309, 2026-07-26)

Retain the half-scale repeated-elliptic preparation

\[
\theta=\frac12
\]

from L303--L304, but include the direct grade-two L298 direction at
order four, as required by the all-series recurrence.  Perform the
preparations in their valid triangular order:

1. insert the grade-one L298 direction at order two;
2. cancel its cubic endpoint response;
3. insert the direct grade-two L298 direction at order four;
4. add L303's parallel quartic lower neutralizer and realize the
   remaining quartic response; and
5. at order five, add the parallel direction which restores the
   complete lower Schur endpoint.

Even after all five steps, the complete upper quintic need not be a
global Markov response.

An exact scalar-copy certificate is L308's stable rational
colligation

\[
S=
\begin{bmatrix}
0&3/5&-4/5\\
0&4/5&3/5\\
0&0&0
\end{bmatrix},
\qquad V=e_1,\qquad W=e_3.                       \tag{1}
\]

Its first two transfers are

\[
B_1=-\frac45,\qquad B_2=\frac9{25}.              \tag{2}
\]

The cubic upper coefficient is exactly zero.  After the quartic
lower neutralizer, the complete lower and physical upper quartics
are

\[
[c^4]\Delta_-=-\frac{81}{1250},\qquad
[c^4]\Delta_+=-\frac{2606}{625}.                 \tag{3}
\]

These are exactly L303's all-series values

\[
\begin{aligned}
[c^4]\Delta_-&=-\theta B_2^*B_2,\\
[c^4]\Delta_+
&=4\theta B_2B_2^*
-(8\theta+4\theta^2)B_1B_1^*\\
&\quad -(8\theta-4\theta^2)(B_1B_1^*)^2 .
\end{aligned}                                    \tag{4}
\]

Before the fifth lower neutralization,

\[
[c^5]\Delta_-=\frac{3384}{3125}.                 \tag{5}
\]

Add the parallel factor column

\[
C_5=-\frac12V\,\frac{3384}{3125}                \tag{6}
\]

and solve its exact moving Stein recurrence.  The complete endpoints
then satisfy

\[
\boxed{
[c^5]\Delta_-=0,\qquad
[c^5]\Delta_+=\frac{66384}{15625}>0.}            \tag{7}
\]

For a scalar copy the balanced endpoint-response map is identically
zero: every response has trace zero, and its output is a scalar.
Therefore the nonzero number in (7) cannot be hidden in a response.

Hence

\[
\boxed{\text{prepared odd coefficients are not globally
response-null in general.}}                      \tag{8}
\]

This is a second route correction after L308.  L308 rejected odd
parity for the raw one-shot copy energy; L309 rejects it even after
the first triangular preparations and complete lower
neutralization.

L309 does **not** obstruct the finite transfer-flag strategy.  In
(1), \(B_1\ne0\), so the scalar endpoint has no surviving first
flag.  The quintic in (7) is entirely supported on an already active
transfer range (and has the favorable sign for \(c>0\)).  The correct
arbitrary-grade target is therefore weaker:

> after lower neutralization, every odd retained class must be a
> bounded hereditary factor through the cumulative active transfer
> row, while only its surviving-flag compression must be a response.

At even orders the retained class must likewise be a bounded
prior-transfer cost plus the new direct Gram, with the former charged
to earlier margins.  One must not impose a global odd-response
condition on L306--L307's prepared recurrence.

## 2. Exact construction

Let

\[
A(c)=\sum_{j=0}^5c^jA_j,\qquad
D(c)=\sum_{j=0}^5c^jD_j
\]

be the exact canonical operator and Hermitian-gauge defect factor.
The checker evaluates the repository's rational word formulas for
all \(A_j,D_j\) on (1).  It independently reconstructs the raw metric
from

\[
M-A^*MA=DD^*
\]

by exact rational Stein solves.

Let \(F_1,F_2\) be L298's exact retightening frame columns.  Start the
correction factor with

\[
C_2=\theta F_1,\qquad C_4=\theta F_2.
\]

For a chosen correction series \(C(c)\), the metric change \(X(c)\)
is determined coefficientwise by

\[
X-A^*XA
=DC^*+CD^*+CC^*.                                \tag{9}
\]

At order four add

\[
\frac12VJ_\theta,\qquad
J_\theta=
\left(4\theta+\frac{\theta^2}{2}\right)|B_1|^2
+\left(\frac{9\theta}{2}-\frac{3\theta^2}{4}\right)|B_1|^4 .
\]

Exact formal Schur complementation at \(V\) and \(W\) gives
(3)--(5).  Add (6), solve (9) once more, and the same exact formal
shorting gives (7).

No response selection, pseudoinverse, rank decision, floating
factorization, or inferred limiting value enters this certificate.

## 3. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_prepared_quintic_obstruction.py \
  --output \
  experiments/repeated_crabb_prepared_quintic_obstruction_s70226.jsonl
```

The first audit uses exact SymPy rationals, exact word formulas,
coefficientwise rational Stein solves, and exact formal Schur
inverses.  A separate floating reconstruction uses the numerical
canonical factor, L298 directions, SciPy/Numpy Schur routines, and
independent Stein solves.  It returns

\[
4.248576000000
\]

with error below \(2.1\times10^{-14}\) from (7).  The tracked dataset
SHA-256 is

```text
e7aa72f5f20dd92ac9b79af64b15106bb792abf0925ef6c3ca7845619cf4761a
```
