# A polynomial quartic preparation makes the canonical endpoint positive

> **Gauge update (2026-07-25).**  L284 subtracts a globally
> endpoint-null polynomial gauge from \(C_4\), producing an
> endpoint-equivalent quartic column that vanishes whenever \(B_1=0\).
> Later prepared orders must be recomputed in that normalized gauge.

## 1. Result (L232, 2026-07-24)

Retain L231's cubic-prepared balanced equality colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
Q=I-E,
\]

and transfer coefficients

\[
B_j=W^*(S^*)^jV,\qquad R_1=B_1^*B_1.
\]

There is a universal bounded polynomial fourth defect-factor column
\(C_4\), given explicitly below, such that adding \(c^4C_4\) to
L231's exact positive Stein factor changes the physical quartic
upper-gap endpoint to

\[
\boxed{
\begin{aligned}
[c^4]{\cal U}_{\rm final}
={}&12B_2B_2^*+32B_1B_1^*\\
&+56B_1(B_1^*B_1)B_1^*\succeq0.
\end{aligned}}                                    \tag{1}
\]

The column satisfies \(V^*C_4=0\), so its fourth-order lower endpoint
is zero.  It contains no inverse, singular vector, flag projection,
or pseudoinverse.  In particular it is real analytic through every
transfer-rank jump and obeys the dimension-free bound

\[
\boxed{\|C_4\|\le\frac{27}{2}.}                   \tag{2}
\]

Thus L231's quartic matrix-redistribution debt is closed, not merely
its scalar trace obstruction.  Formula (1) is still only one prepared
order of the metric: later mixed odd terms, higher even flags, and
the lower flag beyond quartic remain to be controlled.

## 2. Explicit ten-term column

Use words in \(S,S^*\) in the displayed order and define

\[
\boxed{
\begin{aligned}
C_4=\frac12Q\{&
10S^2(S^*)^2-(S^*)^2S^4(S^*)^2+S^2(S^*)^6
 -S^6(S^*)^2\\
&-7(S^*)^5S^3(S^*)^2
+3(S^*)^4S^3(S^*)^3\\
&+(S^*)^3S^3(S^*)^4+(S^*)^2S^3(S^*)^5\\
&-S(S^*)^3S^4(S^*)^2
+S^2(S^*)^3S^3(S^*)^2\}V .
\end{aligned}}                                    \tag{3}
\]

The braces have coefficient \(\ell^1\)-norm \(27\).  Since
\(\|Q\|,\|S\|,\|V\|\le1\), (2) follows immediately.  The leading
\(Q\) also proves \(V^*C_4=0\).

The column need not vanish on the complete-delay stratum.  There it
is an endpoint-null gauge: its fourth upper and lower responses
vanish exactly.  This is harmless at the present order but means
that every later coefficient must be recomputed in this gauge.

## 3. Exact finite certificate

Let \({\cal F}_4\) be L231's complete quartic Stein forcing and let
\(M_2=[c^2]\widetilde P\).  The constant upper-complement inverse is

\[
R=I-F-\frac23E.
\]

Write

\[
{\cal P}_4=C_4V^*.
\]

Direct reduction of (3) gives a 20-word polynomial for
\({\cal P}_4\).  There is an explicit Hermitian polynomial \(Z_4\),
with 29 Hermitian terms and 47 reduced words, for which

\[
\boxed{
{\cal F}_4+{\cal P}_4+{\cal P}_4^*
=Z_4-S^*Z_4S.}                                    \tag{4}
\]

The same witness obeys the endpoint identity

\[
\boxed{
\begin{aligned}
FZ_4F+FM_2RM_2F
={}&-8FS^*ESF-3F(S^*)^2ES^2F\\
&-14FS^*(ESFS^*E)SF .
\end{aligned}}                                    \tag{5}
\]

Both (4) and (5) are identities in the rational
partial-isometry word algebra.  They use only

\[
SS^*S=S,\qquad S^*SS^*=S^*,\qquad EF=FE=0.
\]

For auditability, the 29 exact rational terms of \(Z_4\) are
returned verbatim by
`quartic_coboundary_witness()` in the checker.  Substitution leaves
zero words in (4), zero words in (5), and zero words in
\(Z_4-Z_4^*\).  This is a finite symbolic certificate rather than a
floating regression.

## 4. From the certificate to the positive endpoint

Add \(c^4C_4\) to L231's exact Stein factor and complete the metric
with the stable variable Stein inverse.  At order four the new metric
coefficient is

\[
{\cal G}_S({\cal F}_4+{\cal P}_4+{\cal P}_4^*).
\]

Equation (4) and Stein telescoping make this coefficient exactly
\(Z_4\).  The balanced upper-gap Schur coefficient is its negative
compression, minus the already present second-order Schur square:

\[
[c^4]{\cal U}_{\rm bal}
=-FZ_4F-FM_2RM_2F.                                \tag{6}
\]

Compress (5) by \(W\).  The three state words become, in order,

\[
B_1B_1^*,\qquad B_2B_2^*,\qquad
B_1(B_1^*B_1)B_1^* .
\]

Physical congruence multiplies the left-defect endpoint by four.
Equations (5)--(6) therefore give precisely (1), including all three
coefficients and the sign.

Since \(SV=0\), every perpendicular free column \(C\) satisfies

\[
V^*{\cal G}_S(VC^*+CV^*)V=0.                     \tag{7}
\]

Thus (3) does not change the fourth lower endpoint.

## 5. Scope at rank changes and delays

The right side of (1) is a sum of oriented left Grams.  Its kernel is

\[
\ker B_1^*\cap\ker B_2^*.
\]

Hence:

1. no partial rank change can make the prepared quartic face
   negative;
2. on \(B_1=0\), it reduces to the exact delayed face
   \(12B_2B_2^*\);
3. on \(B_1=B_2=0\), it vanishes and correctly passes control to the
   next ordered flag; and
4. no commutant or reducibility exception remains at this order.

The nonvanishing endpoint-null value of \(C_4\) on a delayed
colligation is a gauge choice, not a divisibility singularity:
\(\|C_4\|\) remains uniformly bounded by (2).  It can influence
orders five and above, so those orders must use the prepared factor,
not L227's original canonical series.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_canonical_quartic_preimage.py \
  --output \
  experiments/repeated_crabb_canonical_quartic_preimage_s70224.jsonl
```

The checker first verifies (4)--(5) and Hermiticity with exact
rational word reduction.  It then reconstructs the full physical
quartic endpoint and the polynomial response on:

1. eight unstructured partial isometries;
2. sixteen rank-changing chains at multiplicities three through six
   and scales down to \(0.05\); and
3. three unstructured complete-delay colligations.

All 27 records pass.  The largest direct endpoint discrepancy is
below \(1.3\times10^{-12}\).  Every lower endpoint motion is below
\(2.9\times10^{-15}\).  On complete delays the column remains
bounded while its endpoint response is zero to roundoff, confirming
the endpoint-null gauge mechanism independently of the exact proof.

The tracked SHA-256 is
`062af4708687e52844be182a4daaa129551be172e90a6302d373c073497d03ce`.

## 7. Next gate

Re-expand the twice-prepared exact factor through orders five and
six.  First determine whether the quintic face is a removable
mixed-flag coboundary, as at cubic order.  On
\(B_1=B_2=0\), compare the sextic face with L215/L228's grade-three
result.  The present quartic gauge is nonzero there, so the
comparison must be recomputed rather than inferred from the
unprepared canonical repair.
