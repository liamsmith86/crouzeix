# Exact Schur transport does not close the later circular-normal flag

> **Scope guard.**  This note disproves one proposed proof shortcut.
> It is not a counterexample to Crouzeix's conjecture, to the repeated
> Crabb neighbourhood theorem, or to L197/L199.

## 1. Result (L319, 2026-07-26)

L199's positive first disk/circular-normal face, its reducing kernel,
the one-order gain in the cross block, and the exact Schur
congruence do **not** by themselves prove that every later Schur
quotient is positive or is a freshly recentered raw L199 face.

More precisely, there are polynomial Hermitian germs
\({\cal S}_{\rm disk}(s)\) and \({\cal S}_{\rm joint}(s)\) with:

1. the same positive semidefinite first face and the same reducing
   first kernel;
2. the same active and cross blocks, with the cross gaining one
   order;
3. a positive disk-only Schur quotient;
4. a joint correction invisible through five orders; but
5. a strictly negative joint Schur quotient.

Thus the sentence

> the old normal/Hardy transport is precisely the cross block, so
> paying \(C^*A^{-1}C\) leaves a fresh raw face

requires an endpoint-specific identity.  It is not a consequence of
Schur algebra, first-face positivity, or kernel reduction.

The live repeated merger must derive the complete later
circular-normal critical graph (including its kernel-block term)
from L194's prepared endpoint.  L197's disk positivity and L199's
first raw absorption remain valid.

## 2. Exact two-by-two certificate

For a real parameter \(s\), put

\[
A=s^2,\qquad C=s^3,
\]

and

\[
\begin{aligned}
{\cal S}_{\rm disk}(s)
 &=
\begin{bmatrix}
s^2&s^3\\
s^3&s^4+s^6
\end{bmatrix},\\
{\cal S}_{\rm joint}(s)
 &=
\begin{bmatrix}
s^2&s^3\\
s^3&s^4-s^6
\end{bmatrix}.
\end{aligned}                                      \tag{1}
\]

Both have first nonzero coefficient

\[
[s^2]{\cal S}
=
\begin{bmatrix}1&0\\0&0\end{bmatrix}\succeq0.       \tag{2}
\]

Its range is the first coordinate and its reducing kernel is the
second coordinate.  The cross is \(C=s^3\), exactly one order later
than \(A=s^2\), so \(A^{-1}C=s\) is analytic.  The transported square
is

\[
C^*A^{-1}C=s^4.                                    \tag{3}
\]

Exact Schur complementation gives

\[
\begin{aligned}
{\cal S}_{\rm disk}/A
 &=(s^4+s^6)-s^4=s^6>0,\\
{\cal S}_{\rm joint}/A
 &=(s^4-s^6)-s^4=-s^6<0
\end{aligned}                                      \tag{4}
\]

for \(s\ne0\).  Equivalently,

\[
\det{\cal S}_{\rm disk}=s^8,\qquad
\det{\cal S}_{\rm joint}=-s^8.                     \tag{5}
\]

The joint correction is only

\[
{\cal S}_{\rm joint}-{\cal S}_{\rm disk}
=
\begin{bmatrix}0&0\\0&-2s^6\end{bmatrix}.          \tag{6}
\]

It changes neither the first face, its kernel, the active block, the
cross block, nor the exact transported square.  It changes precisely
the later kernel-block coefficient which the proposed shortcut had
not derived.

## 3. Consequence for the live proof

The missing theorem must use the actual repeated endpoint, not an
abstract analytic PSD flag.  A sufficient replacement would be one
of the following.

1. Derive the full later normal response after L197
   orthogonalization and prove that its kernel-block contribution is
   L199-absorbable.
2. Prove an exact prepared-endpoint Gram/quotient representation in
   which every normal transport term is visibly included in one
   positive square.
3. For the scalar conjecture only, use L205's scalar-channel
   rigidity to bypass channel-free complete-similarity anchors, with
   a quantitative leakage estimate near the channel strata.

Until one of these is established, do not label the later
circular-normal flag or the L318 merger as proved.

## 4. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_later_circular_flag_obstruction.py \
  --output \
  experiments/repeated_crabb_later_circular_flag_obstruction_s70224.jsonl
```

The checker reconstructs (1), performs both exact Schur
complements, verifies (2)--(6), and tests the signs at three positive
rational parameter values.  The symbolic identities, not the
floating sign samples, prove L319.

The tracked dataset regenerates byte for byte with SHA-256

```text
d272005024eceb7d34da8212871a3af27c6778dc01629bba48c85f9897971952
```
