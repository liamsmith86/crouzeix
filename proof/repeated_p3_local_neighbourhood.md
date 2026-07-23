# A complete-`2` neighbourhood of every fixed repeated `C3` block

## 1. The theorem

Fix a finite multiplicity `m` and put

\[
 M_m=I_m\otimes C_3.
\]

There is a neighbourhood `U_m` of `M_m` in the complex `3m x 3m` matrices
such that every `A in U_m` satisfies

\[
\boxed{W(A)\text{ is a complete `2`-spectral set for }A.} \tag{1}
\]

Equivalently, for a normalized Riemann map `phi_A` of `W(A)`,

\[
 t_*(\phi_A(A))\le4.                                  \tag{2}
\]

The radius may depend on `m`.  This is a local theorem at repeated
three-level Crabb blocks, not the general Crouzeix conjecture.

## 2. The finite equality stratification

L61 splits first directions into two cases.  A positive support-compression
Jensen gap gives a strict negative first endpoint.  Its zero set has a
maximal winner/loser decomposition, and L86 gives the complete winner
endpoint

\[
\begin{aligned}
 {\cal E}_W={}&-16\operatorname{mean}
   \{\lambda_{\max}(Q)I-Q\}
   -8H_1^2-\frac{21}{4}|v|^2I\\
 &-\frac{25}{8}A_0A_0^*
   -8A_1A_1^*
   -\frac{50}{9}A_2A_2^*.                            \tag{3}
\end{aligned}
\]

Every summand is negative semidefinite.  L87 identifies their common
kernel exactly: on a kernel vector, the generator-one variable `H_1`, the
common strong coordinate `v`, and all winner--loser coefficients vanish,
while the remaining support family has a common top vector.  After promoting
every common losing kernel into the winner space, the only surviving
variables are precisely L87's flat copy core.

L113 proves a full neighbourhood theorem on that flat core for fixed `m`.
Its zero-margin strata reduce, by L110--L112, to direct sums of the
single-block disk critical manifolds from L71--L73, together with harmless
unitary and affine coordinates.  Denote their union by `D_m`.

Winner dimension, losing dimension, common-top dimension, and the L93 flag
dimensions take only finitely many values.  They give a finite rank
stratification of a compact normalized slice.  At a rank drop, a common
kernel is promoted to the adjacent larger stratum; no inverse gap is carried
across the change.

## 3. Normal coercivity away from the flat tangent

Fix one stratum and use a local affine-unitary slice.  Split a perturbation
coordinate into a flat tangent `x` and a complementary strong coordinate
`y`.  Formula (3), L82's strict losing lift, and L87 imply

\[
 {\cal Q}_x(y)\preceq0,\qquad
 \ker{\cal Q}_x\text{ consists exactly of reducing flat
 blocks untouched by }y.                               \tag{4}
\]

The kernel in (4) is not inferred by dimension: L87 gives the vanishing
equations, and L111--L113 prove that their solutions are exactly
unitary/blockwise flat motion.  Split off every maximal reducing flat kernel
block and take the Schur endpoint on the remaining active quotient.  On the
chosen normal slice, compactness therefore gives

\[
 \boxed{{\cal Q}_{x,\rm act}(y)
        \preceq-c\|y_{\rm act}\|^2I_{\rm act}}         \tag{5}
\]

for some `c>0` locally on the stratum.  A component of `y` not represented
in `y_act` leaves another reducing flat block and is absorbed into `x`;
there is no uncontrolled endpoint direction.  The losing sector is strictly
feasible by L82 whenever its mean gap is nonzero.  If that gap tends to zero,
its common kernel is promoted to the winner space, which is one of the
finitely many adjacent strata already included above.

Similarly, on the complement of a neighbourhood of the zero first-Jensen
face, L61's first endpoint has one uniform negative margin.  Thus the only
remaining analytic issue is a tube around `D_m`, where the first endpoint
vanishes and (5) is the first strong term.

## 4. Ambient stationarity on the critical stratum

L72 proves that the explicit condition-four certificate on every
single-block disk curve is stationary after the first conformal correction
under **every complex ambient matrix direction**.  Take direct sums of those
certificates.  Domain monotonicity contracts each summand for the convex-hull
numerical range.

As in L106--L110, use the Riemann map of the whole block-diagonal numerical
range and then compose with the inclusion disk map for the actual domain.
L107 embeds the resulting PSD Stein Schur slack, and L108 moves to the
zero-slack branch.  Common scalar/conformal terms cancel from the endpoint
jet by the general form of L109.  Consequently L72's stationarity survives
direct sums and the whole-domain normalization:

\[
 D_y{\cal E}(x,0)=0\qquad(x\in D_m).                  \tag{6}
\]

Winner--loser and reducing-block cross variables also have an exact sign
unitary which sends `y` to `-y`; their active diagonal endpoint is even.
Thus (6) does not rely on differentiating an active support branch.

## 5. Uniform strong-variable remainder

Use L105's exact analytic lower/Stein chart and subtract the inherited
critical-stratum endpoint before expanding.  L106 supplies a fixed-contour
factor for every off-stratum operator term.  Equation (6) removes all terms
linear in the strong variable.  Hence, for `x` at distance `s` from the
repeated base and a normal displacement `y`,

\[
 {\cal E}(x,y)
 ={\cal E}(x,0)+{\cal Q}_x(y)+{\cal R}(x,y),\qquad
 \|{\cal R}(x,y)\|
 \le C(s\|y\|^2+\|y\|^3).                             \tag{7}
\]

At support crossings, the same estimate follows sequentially from L101.
After the common lower-order recentering, the residual support perturbation
is quadratic in `y`; the normalized profiles are a compact continuous
family, so the radial/logarithmic proof transfers it to the Riemann map with
an `o(||y||^2)` remainder.  No smooth largest-eigenvalue selection is used.

Combine (5) and (7).  After shrinking `s` and `||y||`,

\[
 {\cal E}_{\rm act}(x,y)
 \preceq{\cal E}_{\rm act}(x,0)
 -\frac c2\|y_{\rm act}\|^2I_{\rm act}.              \tag{8}
\]

L113 gives a nonpositive endpoint on every split-off flat block and on
`E_act(x,0)`.  Cross endpoint blocks are handled by the same Schur
complements used in L86/L113; their quadratic contribution is already part
of (5), and higher terms are in (7).  Hence the full exact zero-slack chart
metric at `(x,y)` satisfies the upper bound.  Its lower and Stein constraints
hold automatically by L105.

## 6. First-Jensen and rank-change patches

It remains to justify that the preceding local estimates cover every
approach to the repeated base.

- If the normalized first-Jensen gap stays positive, L61's strict
  first-order lift applies.
- If it tends to zero, pass to the zero-Jensen stratum and use (3).
- If any term in (3) has a normalized margin, its negative second-order form
  and the strict-lift argument apply.
- If every term tends to zero, L87 puts the direction in the flat tangent,
  where L113 applies; the normal residual is then controlled by (8).
- If winner, losing, common-top, or flag rank changes, enlarge the kernel
  and move to the adjacent member of the finite stratification.  L82 and
  L113 were formulated with maximal kernels precisely for this transition.

This is also a contradiction proof: a hypothetical countersequence has a
normalized subsequence in one stratum.  The first nonzero transverse
coefficient is strict unless it is tangent to the flat stratum; choosing a
nearest flat point rules out the latter, while L113 controls sequences
remaining in the flat tube.

The normalized slice modulo the compact copy-unitary group is compact.
Finitely many first-Jensen, strong-normal, flat-tube, and rank-change patches
therefore give one neighbourhood `U_m`.  This proves (1)--(2).

## 7. Dependencies and audit

The load-bearing exact identities regenerate from

```bash
.venv/bin/python -u experiments/repeated_p3_multiwinner_gap.py
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
.venv/bin/python -u experiments/repeated_p3_flat_metric_flag.py
.venv/bin/python -u experiments/repeated_p3_slack_transfer.py
.venv/bin/python -u experiments/repeated_p3_scalar_support_rigidity.py
.venv/bin/python -u experiments/repeated_p3_normal_collision.py
```

They prove the complete quadratic form (3), its exact kernel, the flat-flag
coercivity, the joint normal/transverse endpoint jet, and the absence of
hidden multiplicity tangents.  L101, L105--L110, and the finite-dimensional
strict-lift/compactness argument supply the analytic passage from those jets
to the exact neighbourhood theorem.
