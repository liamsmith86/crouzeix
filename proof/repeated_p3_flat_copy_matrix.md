# Flat copy-matrix reduction and the normal stratum (2026-07-22)

## 1. One matrix encodes every flat graph and loop

On L87's higher-order center, generator one and the common strong mode vanish.
Let `X_0,Y_0` be L74's generator-zero pair.  Since

\[
 X_0=-\frac34Y_2,
 \qquad X_2=-\frac43Y_0,                               \tag{1}
\]

every generator-zero/generator-two edge and every complex generator-zero loop
is encoded by one arbitrary copy-space matrix `Z`:

\[
 \boxed{E_Z=Z^*\otimes X_0+Z\otimes Y_0.}              \tag{2}
\]

Indeed, for an oriented edge `a<b`, set

\[
 Z_{ab}=\alpha_0^{ab},
 \qquad Z_{ba}=-\frac43\overline{\alpha_2^{ab}},       \tag{3}
\]

and on a loop set `Z_aa=gamma_a`.  Equations (1)--(3) reproduce both physical
blocks of every L74 edge and L84's diagonal loop.  Conversely, these entries
recover `Z`, so (2) is a bijective parameterization of the complete flat
graph/loop quotient.

The remaining common flat single-block mode is

\[
 W(w)=w(E_{10}+E_{21}),
\]

applied identically to every copy.

## 2. Effective support formula

The first top-support compression of (2) vanishes identically.  Exact reduced-
resolvent multiplication gives the compact second effective matrix

\[
 \boxed{
 Q_Z(q)=\frac5{128}(ZZ^*+Z^*Z)
 -\frac3{128}\{q^2Z^2+q^{-2}(Z^*)^2\}.}               \tag{4}
\]

Formula (4) replaces all orientation-dependent edge bookkeeping by ordinary
copy-matrix algebra.  It also explains several earlier cases:

- a two-copy off-diagonal `Z` has both `ZZ*+Z*Z` and `Z^2` scalar, which is
  L75's second-order pair flatness;
- a generator-zero star is a rank-one row/column configuration in `Z`;
- a length-two directed path has `Z^2!=0`, producing the oscillating second-
  support coupling seen in L83.

## 3. Exact complete-2 theorem on the normal stratum

Suppose `Z` is normal and `Z,w` are sufficiently small.  Choose a copy unitary
`U` with

\[
 U^*ZU=\operatorname{diag}(z_1,\ldots,z_m).            \tag{5}
\]

Conjugating by `U tensor I_3` turns

\[
 I_m\otimes C_3+E_Z+I_m\otimes W(w)
\]

into the direct sum of

\[
 B_j=C_3+\bar z_jX_0+z_jY_0+W(w).                     \tag{6}
\]

Each `B_j` lies in L73's full single-`C3` neighbourhood, so `W(B_j)` is a
complete `2`-spectral set for `B_j`.  Let

\[
 \Omega=W(\operatorname{diag}(B_1,\ldots,B_m))
        =\operatorname{conv}\bigcup_jW(B_j).           \tag{7}
\]

For every matrix-valued rational function `F` analytic on `Omega`, domain
monotonicity gives

\[
 \|F(B_j)\|\le2\sup_{W(B_j)}\|F\|
             \le2\sup_\Omega\|F\|.
\]

Taking the maximum over the direct-sum blocks proves

\[
 \boxed{\Omega\text{ is a complete `2`-spectral set for the full repeated
 matrix whenever `Z` is normal and small.}}            \tag{8}
\]

This is an exact local theorem on an arbitrary-multiplicity repeated-block
stratum, not merely a Taylor-sign statement.

## 4. Nonnormal core (subsequently closed raywise by L93)

Every `Z` is unitarily triangularizable.  In Schur form its strictly upper
part contributes generator-zero edges only; the diagonal is already covered
by (8).  Thus the unresolved L87 center is now a nonnormal upper-triangular
copy matrix, with its second support governed by (4).

The pure two-copy nonnormal block is L77's strict cubic model, and rank-two
stars are closed at second order by L80.  L93 subsequently avoids a
size-by-size Schur classification: its metric-kernel flag is strict unless it
isolates a reducing scalar-support block, and a Clifford reduction shows that
every irreducible terminal block has size at most two.  Thus all fixed
nonnormal directions are now closed.  The remaining problem is uniformity as
the flag/Jensen and two-copy cubic margins collapse together.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_flat_copy_matrix.py
```

The checker uses a fully symbolic `2 x 2` matrix with independent formal
adjoints.  It verifies the zero first compression and proves (4) entry by
entry.  The proof of (4) is the dimension-free Kronecker expansion; the
symbolic matrix checks every one of its four scalar tensor coefficients.
