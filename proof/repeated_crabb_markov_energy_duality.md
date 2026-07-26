# Bounded Markov selection is one positive-test flux inequality

## 1. Result (L282, 2026-07-25)

Let

\[
K=I-\Phi\Phi^*
\]

act on the real Hilbert space of Hermitian copy matrices, with the
Hilbert--Schmidt inner product.  Thus \(K\succeq0\) by L280.  Let
\(U:\mathbb C^d\to\mathbb C^m\) be a copy-flag isometry and put
\(\widehat Y=UYU^*\).

For a Hermitian physical face \(E\), assume L222's strict pointwise
feasibility condition.  Define its least Markov Dirichlet cost by

\[
\mathcal Q(E;K,U)=
\inf\left\{
\langle H,KH\rangle:
U^*(E+8KH)U\preceq0,\ H=H^*
\right\}.                                         \tag{1}
\]

Then the exact convex dual is

\[
\boxed{
\mathcal Q(E;K,U)
=\sup_{Y\succeq0}
\left\{
\langle Y,U^*EU\rangle
-16\langle\widehat Y,K\widehat Y\rangle
\right\}.}                                        \tag{2}
\]

Scaling each nonzero \(Y\) also gives

\[
\boxed{
\mathcal Q(E;K,U)
=\frac1{64}
\sup_{\substack{Y\succeq0\\
\langle\widehat Y,K\widehat Y\rangle>0}}
\frac{\langle Y,U^*EU\rangle_+^2}
{\langle\widehat Y,K\widehat Y\rangle}.}           \tag{3}
\]

Directions with zero denominator have negative numerator by L222 and
therefore contribute zero to the supremum.

Consequently the entire bounded-selection problem is equivalent to
one quantitative strengthening of L279.  It is enough to prove,
uniformly on every active flag,

\[
\boxed{
\langle Y,U^*EU\rangle_+
\leq \gamma\,\|U^*B_k\|_F
\sqrt{\langle\widehat Y,K\widehat Y\rangle}
\quad(Y\succeq0).}                                \tag{4}
\]

Indeed, (3) then gives

\[
\mathcal Q(E;K,U)
\leq\frac{\gamma^2}{64}\|U^*B_k\|_F^2.            \tag{5}
\]

By L281, the balanced physical state column realizing the minimizing
response has squared norm \(2\mathcal Q\).  Thus (4) is precisely the
gap-free estimate needed for analytic remainder domination.

L282 does **not** prove (4).  It replaces a near-singular Poisson
solve by a dual off-commutant flux estimate.  L279 already proves the
zero-denominator endpoint of (4): if
\(\langle\widehat Y,K\widehat Y\rangle=0\), then
\(\langle Y,U^*EU\rangle<0\).

## 2. Reduction to an ordinary Hilbert variable

Put

\[
x=K^{1/2}H.
\]

Only the component of \(H\) orthogonal to \(\ker K\) affects either
the response or its cost.  Hence (1) is equivalently

\[
\inf\left\{
\|x\|^2:
U^*\{E+8K^{1/2}x\}U\preceq0,\quad
x\in\operatorname {ran}K
\right\}.                                         \tag{6}
\]

For a dual matrix \(Y\succeq0\), the Lagrangian is

\[
\|x\|^2+\langle Y,U^*EU\rangle
+8\langle K^{1/2}\widehat Y,x\rangle.              \tag{7}
\]

Since \(K^{1/2}\widehat Y\in\operatorname {ran}K\), its infimum over
\(x\) occurs at

\[
x=-4K^{1/2}\widehat Y
\]

and equals

\[
\langle Y,U^*EU\rangle
-16\langle\widehat Y,K\widehat Y\rangle.           \tag{8}
\]

L222 supplies a strictly feasible point, so finite-dimensional Slater
duality proves (2).  This argument remains valid when \(K\) is
singular and uses no pseudoinverse.

## 3. Homogeneous ratio

Fix \(Y\succeq0\), and abbreviate

\[
e=\langle Y,U^*EU\rangle,\qquad
d=\langle\widehat Y,K\widehat Y\rangle.
\]

Replacing \(Y\) by \(tY\), \(t\geq0\), turns (8) into

\[
te-16t^2d.                                        \tag{9}
\]

If \(e\leq0\), its maximum is zero.  If \(e>0\), strict feasibility
forces \(d>0\), and the maximum is attained at
\(t=e/(32d)\), with value \(e^2/(64d)\).  Taking the supremum over
positive rays proves (3).

## 4. Physical interpretation

L281 identifies

\[
\langle\widehat Y,K\widehat Y\rangle
=2\|R_{\widehat Y}\|_F^2.                          \tag{10}
\]

Thus (4) can be written as

\[
\langle Y,U^*EU\rangle_+
\leq\sqrt2\,\gamma\,\|U^*B_k\|_F
\|R_{\widehat Y}\|_F.                              \tag{11}
\]

The right side is the product of:

1. the first active physical transfer amplitude; and
2. the exact observability defect measuring how far \(Y\) is from a
   genuine reducing-copy obstruction.

This is the natural quantitative continuation of L279.  Its
zero-defect case is already strict and negative.  The next proof
should extend L279's weighted trace calculation off the commutant and
bound every resulting error term by the product in (11), rather than
constructing or estimating a Poisson inverse.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_markov_energy_duality.py \
  --output \
  experiments/repeated_crabb_markov_energy_duality_s70225.jsonl
```

The deterministic audit uses the two-copy depolarizing Markov
Laplacian \(K(H)=H-\operatorname {tr}(H)I/2\).  For diagonal faces
\(E=\operatorname {diag}(a,b)\), \(a>0\) and \(a+b<0\), both sides of
(2)--(3) equal \(a^2/32\).  The checker compares the closed form with
independent primal and dual semidefinite programs.  The exact
argument above, not the floating audit, proves L282.  The tracked
dataset has SHA-256

```text
559e4a7ee92179c2cca7085e9942433158bbb700902988bb3a865d7451369a33
```
