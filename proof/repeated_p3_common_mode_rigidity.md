# Common-mode rigidity on the flat copy stratum (2026-07-22)

## 1. Unique outer Laurent modes

L90 gives, up to a scalar multiple of the identity,

\[
 Q_{Z,w}(q)=Q_Z(q)+\ell_w(q)Z+overline{\ell_w(q)}Z^*,                 \tag{1}
\]

with

\[
 \ell_w(q)=\frac{\sqrt2}{16}(q^{-1}w-q^3\bar w).                     \tag{2}
\]

The quadratic term `Q_Z` has only Fourier modes `0,+2,-2`.  Hence the
non-scalar coefficients at the unique outer modes are

\[
 [q^3]Q_{Z,w}=-\frac{\sqrt2}{16}\bar w Z,
 \qquad
 [q^{-3}]Q_{Z,w}=-\frac{\sqrt2}{16}wZ^*.                              \tag{3}
\]

## 2. Every common top vector reduces `Z`

Assume `w!=0` and a nonzero vector `x` is an eigenvector of `Q_{Z,w}(q)` for
every unit `q`.  Project the Laurent identity

\[
 Q_{Z,w}(q)x=\lambda(q)x
\]

onto `x^perp`.  Every Fourier coefficient must vanish separately.  Equations
(3) give

\[
 P_{x^\perp}Zx=0,
 \qquad P_{x^\perp}Z^*x=0.                            \tag{4}
\]

Thus

\[
 Zx=\alpha x,qquad Z^*x=\bar\alpha x,                \tag{5}
\]

so `span{x}` is a reducing subspace of `Z`.  In particular:

\[
 \boxed{w\ne0\text{ and `Z` unitarily irreducible }
 \Longrightarrow\text{ no common top branch.}}        \tag{6}

By L87, the matrix-Jensen term is then strictly negative on every candidate
certificate-maximizing vector.

## 3. Equality blocks are already single-block cases

Decompose `Z` into its minimal reducing blocks.  Relation (5) says that any
block carrying an L90 common top vector must be one-dimensional.  On that
block the full repeated matrix is simply

\[
 C_3+\bar zX_0+zY_0+W(w),                              \tag{7}
\]

which lies in L73's single-`C3` neighbourhood.  Its numerical range is a
complete `2`-spectral set.  The numerical range of the full direct sum contains
that block's numerical range, so domain monotonicity preserves its bound.

Any higher-dimensional irreducible reducing block has strict second-order
matrix-Jensen margin by (6).  Therefore the common `w` mode introduces no new
irreducible higher-order equality mechanism: equality is carried only by
single-block summands already controlled by L73.

This is a fixed-stratum reduction, not yet a uniform estimate as `w->0`; that
weighted transition returns to the `w=0` nonnormal copy-matrix problem.

## 4. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_copy_matrix.py
```

The exact L90 checker proves (1)--(3).  The implication (3)--(5) is coefficient
comparison in a finite Laurent polynomial and uses no numerical threshold.
