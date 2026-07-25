# The canonical Stein repair fails the partial-flag metric bound

## 1. Status (A172, 2026-07-24)

Retain L219--L228's balanced equality colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

the boundary metric \(P_{\rm bl}(c)\), and its elliptic Stein slack

\[
H(c)=P_{\rm bl}(c)-\widehat T(c)^*
                         P_{\rm bl}(c)\widehat T(c).
\]

L227 defines

\[
\begin{aligned}
G(c)&=V^*H(c)V,\\
{\cal K}(c)&=H(c)-H(c)VG(c)^{-1}V^*H(c),\\
X_{\rm can}(c)&=-{\cal G}_{\widehat T(c)}({\cal K}(c)).
\end{aligned}                                      \tag{1}
\]

The repaired slack is exactly positive:

\[
(P_{\rm bl}+X_{\rm can})
-\widehat T^*(P_{\rm bl}+X_{\rm can})\widehat T
=HVG^{-1}V^*H\succeq0.                             \tag{2}
\]

However, the canonical repaired metric does **not** in general obey
the required upper bound through a rank-changing partial flag:

\[
\boxed{
P^{1/2}(P_{\rm bl}+X_{\rm can})P^{1/2}
\npreceq4I.}                                       \tag{3}
\]

This is a falsification of a proposed metric, not of Crouzeix's
conjecture and not of L227's contraction identity.

The obstruction is already cubic.  Let \({\cal U}_{\rm can}(c)\) be
the physical upper-gap Schur complement after eliminating the state
complement of \(W\mathbb C^m\).  Its second coefficient is L227's
favorable face

\[
[c^2]{\cal U}_{\rm can}=12B_1B_1^*.               \tag{4}
\]

If \(B_1\) has rank one and \(U\) is an isometry onto
\(\ker B_1^*\), the tracked noncommuting examples have

\[
\boxed{
U^*[c^3]{\cal U}_{\rm can}U
\quad\hbox{indefinite}.}                           \tag{5}
\]

For example, at multiplicity four its eigenvalues are

\[
(-0.0119640699,\ 0.0000766396,\ 0.0090452341).     \tag{6}
\]

Since (4) has no active-to-kernel cross row, the copy-space Schur
cross-square starts only at order \(c^4\).  Thus a negative
eigenvalue in (5) implies

\[
\lambda_{\min}{\cal U}_{\rm can}(c)<0
\]

for every sufficiently small positive \(c\).  This proves (3) for
each audited colligation up to the numerical construction error.

## 2. Why this is a local obstruction

The examples use L224's lossless Schur realizations with

\[
\Gamma_0=0,\qquad
\Gamma_j=\lambda\Gamma_j^{(0)}
\quad(1\le j<m),
\]

and a fixed terminal unitary.  Hence they converge to the repeated
monomial colligation as \(\lambda\to0\).

For multiplicities four and five, (5) remains indefinite at

\[
\lambda=1,\ 0.5,\ 0.2,\ 0.1.
\]

At multiplicity four the negative eigenvalue is respectively

\[
-1.1964\,10^{-2},\quad
-1.4746\,10^{-3},\quad
-9.2061\,10^{-5},\quad
-1.1382\,10^{-5}.                                  \tag{7}
\]

The approximately cubic scaling in \(\lambda\) is consistent with a
mixed first/second Schur-parameter term.  Therefore shrinking the
equality neighbourhood does not make the unmodified canonical
repair a certificate: after fixing any nonzero \(\lambda\), taking
\(c\) sufficiently small exposes the cubic face before the next
even transfer Gram.

## 3. Compatibility with the complete-delay results

There is no contradiction with L228.  On the exact delayed stratum
\(B_1=0\), every coefficient below \(c^4\), including the cubic
one, vanishes.  L228's candidate then predicts the favorable
grade-two even face.

The obstruction occurs while **approaching** that stratum:
\(B_1\) is nonzero but singular, its \(c^2\) Gram is removed on
\(\ker B_1^*\), and the mixed cubic term survives there.  Thus:

1. proving L228's complete-delay anticommutator remains valuable;
2. complete-delay covariance alone cannot prove the partial-flag
   metric sandwich for the canonical repair; and
3. the odd mixed flags must be cancelled before using the delayed
   even Grams.

The full cubic coefficient has trace zero in every tracked case,
to at worst \(1.9\times10^{-14}\).  This agrees with L204's fact
that homogeneous free-row endpoint responses are trace-free.  It
strongly suggests that the cubic obstruction should be removed by
an analytic defect-frame/free-row correction rather than estimated
against the even boundary budget.

## 4. Corrected next target

The next constructive step should be:

1. derive \( [c^3]{\cal U}_{\rm can}\) in ordered transfer words;
2. prove its compression pairs to zero with every L222 reducing
   commutant separator;
3. exhibit a bounded analytic preimage under L204's endpoint map,
   preferably a polynomial polarization of L208--L212; and
4. insert that correction through L194's exact lower-tight/Stein
   chart before advancing to the \(c^4\) flag.

This is a finite mixed-flag preparation problem.  It is narrower
than constructing the entire one-image series, but it cannot be
skipped by the exact contraction repair.

## 5. Regeneration and numerical guard

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_canonical_repair_flag_obstruction.py \
  --output \
  experiments/repeated_crabb_canonical_repair_flag_obstruction_s70224.jsonl
```

The eight records cover multiplicities four and five and four scales
approaching the repeated apex.  They regenerate the direct ellipse
series, boundary slack, canonical Schur residual, variable Stein
repair, state endpoint Schur complement, and the cubic copy-kernel
compression.  The tracked SHA-256 is
`a3ae62c284c85e79fcae48dc52b1a352362bba6b8ea4fdd7b5de22f2cb56e0e2`.

During this audit, SciPy's automatic real bilinear Lyapunov branch
failed its defining equation on a nonnormal \(12\times12\) rational
partial isometry.  The shared Stein solver now checks its residual
and retries in complex arithmetic.  This numerical issue does not
cause (5): the tracked Schur realizations are complex, and every
repair coefficient satisfies its Stein equation to below
\(4.3\times10^{-15}\).
