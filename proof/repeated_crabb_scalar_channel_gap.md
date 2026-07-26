# Quantitative scalar-channel gap on repeated disk equality anchors

> **Campaign scope.**  This is a scalar Crouzeix estimate on L193's
> repeated disk equality manifold.  It does not assert the stronger
> complete-similarity bound for circular-normal perturbations.

## 1. Result (L320, 2026-07-26)

Let \(T\) be one L193 inverse-block-Toeplitz disk anchor.  In L201's
balanced coordinates let

\[
C=P^{1/2}TP^{-1/2},
\qquad
P=2I-VV^*+2WW^*,                                  \tag{1}
\]

where \(C\) is a contraction and \(V,W:\mathbb C^m\to{\cal H}\)
are its right and left defect frames.  Write its matrix-inner
transfer as

\[
B(z)=\sum_{n\ge1}B_nz^n,
\qquad B_n=W^*(C^*)^nV.                            \tag{2}
\]

Define the scalar-channel score

\[
\boxed{
\sigma(B)=
\sup_{\|u\|=\|v\|=1}
\sum_{n\ge1}|u^*B_nv|^2\le1.}                     \tag{3}
\]

Then every scalar Schur function \(f\) satisfies

\[
\boxed{
\|f(T)\|^2
\le4-\bigl(1-\sqrt{\sigma(B)}\bigr)^2.}            \tag{4}
\]

By L205, \(\sigma(B)=1\) exactly when \(B\) contains a constant
one-dimensional scalar inner channel

\[
B(z)v=g(z)u.                                       \tag{5}
\]

Consequently every channel-free L193 anchor has an explicit strict
scalar Crouzeix gap.  On any compact equality-anchor set with
\(\sigma(B)\le1-\eta\), the gap is uniform:

\[
\boxed{
\|f(T)\|^2
\le4-\bigl(1-\sqrt{1-\eta}\bigr)^2<4.}             \tag{6}
\]

This supplies the quantitative reserve missing from L205.  It
changes the scalar campaign architecture:

1. channel-free noncommuting equality anchors are strict and have
   ordinary open tubes by continuity;
2. only the scalar-channel strata can support a sharp scalar
   sequence;
3. L205 splits every exact channel into a single-copy full-Hardy
   equality block plus a complementary transfer; and
4. the remaining local work is a finite channel-stratum induction
   combining this leakage reserve with L192, L199, and L318.

L320 does not by itself complete that induction or the repeated
neighbourhood theorem.

## 2. An abstract endpoint-gap lemma

The proof uses only the three eigenvalues of (1).  Let \(F\) be any
contraction on \({\cal H}\), and put

\[
s=\|V^*FW\|\le1.
\]

Then

\[
\boxed{
\|P^{-1/2}FP^{1/2}\|^2
\le4-(1-s)^2.}                                     \tag{7}
\]

To prove this, take a unit vector \(x\), and set

\[
y=P^{1/2}x,\qquad z=Fy,\qquad
\|P^{-1/2}z\|^2=4-\varepsilon.                     \tag{8}
\]

The exact loss decomposition is

\[
\begin{aligned}
\varepsilon
={}&\langle(4I-P)x,x\rangle\\
&+\{\|y\|^2-\|z\|^2\}
+\langle(I-P^{-1})z,z\rangle.                     \tag{9}
\end{aligned}
\]

All three terms are nonnegative.  Let

\[
\alpha=\|(I-WW^*)x\|,\qquad
q=\|(I-VV^*)z\|.
\]

The spectra of \(P\) on \(W^\perp\), and of \(P^{-1}\) on
\(V^\perp\), give

\[
\alpha\le\sqrt{\varepsilon/2},
\qquad q\le\sqrt{2\varepsilon}.                   \tag{10}
\]

Also

\[
\|V^*z\|^2
\ge\|P^{-1/2}z\|^2-\frac12q^2
\ge4-2\varepsilon.                                \tag{11}
\]

Write

\[
x=\beta Wu+x_\perp,\qquad
\beta=\sqrt{1-\alpha^2},\qquad \|u\|=1.
\]

Since \(\|P^{1/2}x_\perp\|\le\sqrt2\alpha\),

\[
\begin{aligned}
\sqrt{4-2\varepsilon}
&\le\|V^*z\|\\
&\le2\beta\|V^*FWu\|+\sqrt2\alpha
\le2\beta s+\sqrt{\varepsilon}.                  \tag{12}
\end{aligned}
\]

For \(0\le\varepsilon\le1\), (12) and \(\beta\le1\) imply

\[
s\ge\sqrt{1-\varepsilon/2}-\frac{\sqrt\varepsilon}{2}.
\]

Therefore

\[
\begin{aligned}
1-s
&\le1-\sqrt{1-\varepsilon/2}
  +\frac{\sqrt\varepsilon}{2}\\
&\le\sqrt\varepsilon.                             \tag{13}
\end{aligned}
\]

If \(\varepsilon\ge1\), (13)'s squared conclusion is automatic
because \(0\le s\le1\).  Thus
\(\varepsilon\ge(1-s)^2\) in all cases, proving (7).

## 3. The transfer bounds the endpoint corner

Let

\[
f(z)=\sum_{n\ge0}a_nz^n,\qquad \|f\|_{\mathbb D}\le1.
\]

Von Neumann's inequality gives \(F=f(C)\) with \(\|F\|\le1\).
Moreover \(V^*W=0\), and L201 gives

\[
V^*FW=\sum_{n\ge1}a_nB_n^*.                       \tag{14}
\]

For unit \(u,v\), scalar Hardy Parseval yields

\[
\begin{aligned}
|v^*V^*FWu|
&=\left|\sum_{n\ge1}
 a_n\overline{u^*B_nv}\right|\\
&\le
\left(\sum_{n\ge1}|u^*B_nv|^2\right)^{1/2}
\le\sqrt{\sigma(B)}.                              \tag{15}
\end{aligned}
\]

Taking the supremum in \(u,v\) gives

\[
\|V^*FW\|\le\sqrt{\sigma(B)}.                     \tag{16}
\]

Finally,

\[
f(T)=P^{-1/2}f(C)P^{1/2}.
\]

Insert (16) into (7) to prove (4).

## 4. Sharp set and continuity

L201's matrix Parseval identity gives \(\sigma(B)\le1\).
L205 proves that equality throughout the Cauchy--Schwarz chain in
(15) is equivalent to the scalar channel (5).  Hence the zero set
of the reserve in (4) is exactly the scalar equality-channel locus,
not the whole noncommuting L193 equality manifold.

On a compact finite-dimensional anchor set, the transfer state
spectral radius stays uniformly below one, so the Taylor sum in (3)
converges uniformly.  Hence \(\sigma(B)\) is the maximum of a
continuous function on two unit spheres and is continuous in the
anchor.  Equation (6) follows on compact channel-free sets.
Ordinary continuity of the numerical range, Riemann functional
calculus, and finite-dimensional operator norm then gives a strict
scalar tube around each such compact set.  No corresponding
complete-similarity claim is made.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_scalar_channel_gap.py \
  --output \
  experiments/repeated_crabb_scalar_channel_gap_s70224.jsonl
```

The checker:

1. tests (7) on deterministic random contractions with one- through
   three-dimensional defect frames;
2. reconstructs noncommuting L193 anchors and their matrix-inner
   transfer coefficients;
3. verifies (14)--(16) on deterministic scalar Schur polynomials;
4. confirms score one at the repeated apex and a strict product-state
   score on the standard noncommuting anchors.

The alternating product-state maximization is only a diagnostic for
\(\sigma(B)\).  Equations (7) and (14)--(16), not that optimizer,
prove L320.

The tracked dataset regenerates byte for byte with SHA-256

```text
973f669a8ae58492c65b73492fbb8b3c5651764ba6aefa0db94adb323c05c40b
```
