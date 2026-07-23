# Uniform arbitrary-copy metric-flag tube (2026-07-22)

## 1. Statement and scope

Fix a finite copy multiplicity `m`.  L93 proves that every fixed direction
in the repeated-`C3` flat core is controlled, but its hierarchical metric
choice was not uniform when the common-top and leakage ranks changed.
L110--L112 close the terminal blocks on which that hierarchy can stabilize.
This note proves the missing gluing statement:

\[
\boxed{\text{The complete flat copy core has a complete-`2`
neighbourhood at }I_m\otimes C_3.}                    \tag{1}
\]

The neighbourhood may depend on `m`; no dimension-free radius is asserted.
The later task is to glue (1) to L86's non-flat and losing-space variables.

## 2. Stable flag decomposition

For normalized flat data `(Z,w)`, start with the common-top space `K_0` of
L93 and iterate

\[
 K_{j+1}=\{x\in K_j:Zx,Z^*x\in K_j\}.                \tag{2}
\]

The chain has length at most `m`.  On `K_j minus K_{j+1}`, L93's metric
variation has compression

\[
 -5\sqrt2\,P_jZ(I-P_j)Z^*P_j
 -\frac{15\sqrt2}{4}\,P_jZ^*(I-P_j)ZP_j.             \tag{3}
\]

Its kernel is exactly `K_{j+1}`.  The stable space `K_infty` reduces `Z`.
On it, every vector is a common top vector, so the support is scalar and

\[
 Z^2=\alpha I,\qquad ZZ^*+Z^*Z=\beta I.              \tag{4}
\]

Minimal reducing pieces of (4) have dimension one or two.  L111 closes
repeated nonnormal two-dimensional pieces and L112 closes their collision
with arbitrary-multiplicity normal pieces.  If `w!=0`, L91 already forces
every common top vector to reduce `Z`, so the same stable decomposition
applies without a further irreducible case.

Thus every zero of all leading flag forms is, after a copy unitary, a direct
sum of:

- normal one-dimensional copy blocks, covered by L88; and
- two-dimensional terminal blocks, covered uniformly by L110.

Direct sums are covered exactly by domain monotonicity.

## 3. The tangent kernel has no hidden directions

Choose a local unitary slice through one such block-diagonal zero stratum.
The common kernel of the Jensen defect and all forms (3) consists of
directions which preserve every stable reducing block, plus tangents inside
the scalar-support blocks.

L111 proves that at a repeated nonnormal block every internal relation
tangent is a unitary-orbit tangent plus the same three two-dimensional
parameters on every multiplicity coordinate.  L112 proves that at a normal
`+/-` collision the only unitary-invariant tangent is a rectangular edge
whose SVD is a direct sum of L110 pairs.  Consequently the tangent kernel
of the complete flag certificate is exactly the tangent space to the
already covered blockwise stratum.

On a complementary normal slice, the finite collection of quadratic flag
forms is therefore coercive.  For each base point there is a hierarchical
metric combination and a constant `c>0` such that

\[
 {\cal E}_{\rm lead}\preceq-cs^2\|R\|^2 I            \tag{5}
\]

on the remaining active endpoint, where `R` is distance in that slice.
Here `s` is the physical distance from the repeated base; equivalently,
after dividing out the common second-order factor `s^2`, the normalized
coercivity is `-c||R||^2`.
This is just the finite-dimensional kernel criterion (3): choose the first
flag form on `K_0 minus K_1`, then successively smaller positive weights on
`K_1 minus K_2`, and so on.  The strict complement of `K_0` retains the
canonical Jensen gap.

The normalized direction space is compact.  Away from any fixed
neighbourhood of the blockwise zero strata, finitely many such metric
choices give one uniform strict leading margin.  It remains only to absorb
remainders in their tubular neighbourhoods.

## 4. Block-sign parity removes linear transverse remainders

Let `P` be one stable reducing block, `Q=I-P`, and split a nearby flat copy
matrix as

\[
 Z=Z_{\rm bd}+R,\qquad
 Z_{\rm bd}=PZP+QZQ,\qquad
 R=PZQ+QZP.                                           \tag{6}
\]

The block-sign unitary `J=P-Q` fixes `Z_bd` and sends `R` to `-R`.
The same copy unitary conjugates the full repeated matrix.  Therefore its
numerical range, normalized Riemann map, similarity value, and an averaged
metric branch are invariant under `R -> -R`.

In particular:

- every diagonal block of the upper Schur endpoint is even in `R`;
- every cross block is odd in `R`;
- after taking the Schur complement of a strictly negative complementary
  endpoint, its contribution to the active block is quadratic in `R`.

The statement is exact and survives support crossings because it is unitary
equivalence, not branch differentiation.

Use the block-diagonal direct-sum certificate as the normal anchor.  As in
L110, compare the Riemann maps of the whole block-diagonal and actual
numerical ranges.  On a proper common-top stratum the complementary support
has a fixed normalized gap; the block-sign symmetry and the cluster Schur
complement improve the support enlargement from a merely transverse bound
to

\[
 0\le h_A-h_{\rm bd}\le Cs^2\|R\|^2.                 \tag{7}
\]

The inclusion disk map and its Stein Schur slack are consequently
`O(||R||^2)`.  Insert that slack into L105's analytic chart and subtract the
inherited block-diagonal endpoint before expanding.

Every term not included in (5) now contains two factors of `R`; local
analyticity gives

\[
 \|{\cal R}_{\rm exact}\|\le C s^3\|R\|^2,           \tag{8}
\]

where `s` is the distance to the repeated base.  After shrinking `s`,
(8) is absorbed by (5).

If the normalized complementary support gap in (7) collapses, its kernel is
promoted into `K_0`; this is precisely the next flag stratum, already
included in the finite decomposition of Section 2.  If the stable block
itself reaches a full-common-top collision, L111--L112 and L110 supply its
tube.  Hence no inverse gap is used uniformly across a rank change.

## 5. Finite-cover conclusion

Proceed by the finite flag length.  At a point with strict Jensen or flag
margin, (5) and (8) give an open complete-`2` neighbourhood.  At a zero
stratum, use the exact blockwise certificate and its L110 tube, then apply
the even transverse estimate (8).  Rank-drop points are handled by enlarging
the stable block and moving to the next member of the same finite
stratification.

For fixed `m`, the unit sphere of normalized flat data modulo the compact
copy-unitary group is compact.  A finite subcover gives one neighbourhood
of the repeated base on the entire flat core, proving (1).

The proof does not assert a uniform radius as `m` grows.  It also does not
yet include L86's generator-one, common-strong, winner--loser, or positive
first-Jensen variables; those have their own explicit negative forms and are
the next gluing step.

## 6. Dependencies and regeneration

The exact algebraic inputs regenerate from

```bash
.venv/bin/python -u experiments/repeated_p3_flat_metric_flag.py
.venv/bin/python -u experiments/repeated_p3_scalar_support_rigidity.py
.venv/bin/python -u experiments/repeated_p3_normal_collision.py
```

The first proves (3) and its exact kernel, the second proves tangent rigidity
of every repeated nonnormal full-common-top block, and the third proves the
normal-collision SVD reduction.  The parity and finite-cover arguments in
Sections 4--5 are exact finite-dimensional consequences of block-unitary
equivariance and L105's analytic chart.
