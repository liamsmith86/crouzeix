# Gau--Wu equality contains non-Crabb sharp disk strata

> **Campaign scope.**  This is an exact course correction for the
> post-L329 globalization plan.  The family below is classical
> Gau--Wu structure, not a new equality theorem.  It proves that the
> known scalar sharp set is strictly larger than the repeated-Crabb
> collision family.

## 1. Result (L330, 2026-07-26)

For real \(a\), \(0<|a|<1\), put

\[
 r_a=\sqrt{2(1-a^2)},\qquad
 G_a=
 \begin{bmatrix}
 0&r_a&-2a\\
 0&a&r_a\\
 0&0&0
 \end{bmatrix},
 \tag{1}
\]

and let

\[
 f_a(z)=z\,{z-a\over1-az}.                        \tag{2}
\]

Then

\[
 \boxed{W(G_a)=\overline{\mathbb D},\qquad
 f_a(G_a)=2E_{13},\qquad \|f_a(G_a)\|=2.}         \tag{3}
\]

Since \(f_a\) is a degree-two finite Blaschke product,
\(\max_{\mathbb D}|f_a|=1\).  Thus every \(G_a\) is an exact scalar
Crouzeix equality matrix.

For \(a\ne0\),

\[
 \sigma(G_a)=\{0,0,a\}.                           \tag{4}
\]

An affine image of a Crabb block has only one eigenvalue.  Therefore
\(G_a\) is not affine-similar, hence not affine-unitarily equivalent,
to \(C_3\).  Direct sums of such blocks similarly give exact sharp
disk strata outside the repeated nilpotent family.

Consequently L329 covers the zero-collision corner \(a=0\) and a
neighbourhood of it, but it cannot by itself cover the global known
equality set.  The post-L329 statement

> every sharp sequence is covered by the repeated-Crabb charts

is false before any numerical search is run.

## 2. Elementary verification

Let \(\zeta=e^{i\theta}\) and

\[
 H_\theta=\operatorname{Re}(\zeta^{-1}G_a).
\]

Direct determinant expansion gives

\[
 \det(\lambda I-H_\theta)
 =(\lambda-1)(\lambda+1)
 \left(\lambda-\frac a2(\zeta+\zeta^{-1})\right).
 \tag{5}
\]

The three support eigenvalues are therefore

\[
 1,\quad -1,\quad a\cos\theta.                    \tag{6}
\]

Because \(|a|<1\), the largest is exactly \(1\) for every
\(\theta\).  The support function of \(W(G_a)\) is identically one,
which proves the disk identity in (3).

The rational functional calculus is also exact:

\[
\begin{aligned}
 f_a(G_a)
 &=G_a(G_a-aI)(I-aG_a)^{-1}\\
 &=2E_{13}.                                       \tag{7}
\end{aligned}
\]

Equations (5)--(7) prove all assertions without relying on a
floating numerical-range reconstruction.

## 3. Relation to Gau--Wu

This is Corollary 3 of Gau--Wu (2009) specialized to the Blaschke
zeros \(0,a\).  Their general theorem gives, for every finite
Blaschke product \(f\) with \(f(0)=0\), an irreducible disk model

\[
 X_\phi S(\phi)X_\phi^{-1},\qquad \phi=zf,
\]

which attains \(\|f(T)\|=2\).  Arbitrary nonzero Blaschke zeros
therefore produce many higher-dimensional nonnilpotent sharp disk
models as well.

The original contribution here is only campaign bookkeeping:
L326 had already recorded the classical model theorem, but the first
post-L329 globalization wording accidentally treated L192/L329's
monomial collision charts as though they covered the whole Gau--Wu
family.  L330 corrects that inference.

## 4. Revised frontier

The next local target is not another repeated-Crabb coefficient.  It
is:

1. prove, or falsify, a scalar neighbourhood theorem at an arbitrary
   fixed Gau--Wu disk model;
2. start with the explicit \(3\times3\) family \(G_a\), where the
   support spectrum (6) has the uniform gap \(1-|a|\); and
3. only after the full disk-model equality manifold is covered,
   ask whether every global sharp sequence approaches that manifold.

The distinction between these two questions is essential.  An
equality-stratum classification cannot globalize L329 if known
equality points themselves lie outside its local charts.

## 5. Audit

Run

```bash
.venv/bin/python -u experiments/gau_wu_noncrabb_sharp_stratum.py \
  --output experiments/gau_wu_noncrabb_sharp_stratum_s70224.jsonl
```

The checker verifies (5), (7), the norm, and the non-Crabb spectral
certificate exactly for several rational \(a\)'s, while also checking
the symbolic identities with \(a\) left free.
