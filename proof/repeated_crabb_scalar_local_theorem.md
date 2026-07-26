# A scalar Crouzeix neighbourhood of every fixed repeated Crabb block

> **Campaign scope.**  This is a fixed-dimension **scalar** local
> theorem.  It does not prove the completely bounded conjecture, does
> not give a radius uniform in \(L\) or \(m\), and does not settle
> arbitrary matrices.  The proof assembles the previously banked
> finite flags; it introduces no new all-grade coefficient claim.

## 1. Result (L329, 2026-07-26)

Fix \(p=L+1\ge3\) and \(m\ge1\), and put

\[
 A_0=C_p\otimes I_m.
\]

There is an operator-norm neighbourhood \({\cal U}_{p,m}\) of
\(A_0\) such that every \(A\in{\cal U}_{p,m}\) satisfies the scalar
Crouzeix conjecture:

\[
 \boxed{
 \|g(A)\|
 \le2\max_{z\in W(A)}|g(z)|
 \quad\hbox{for every polynomial }g.}             \tag{1}
\]

Equivalently, if \(\Phi_A\) is a normalized Riemann map from
\(\operatorname{int}W(A)\) to \(\mathbb D\), then

\[
 \boxed{
 \|f(\Phi_A(A))\|\le2
 \quad\hbox{for every scalar Schur function }f.}   \tag{2}
\]

For \(m=1\), L192 proves the stronger complete-\(2\) statement.
For \(m>1\), (1)--(2) are scalar only.  L319 remains a valid
obstruction to upgrading this proof to a complete-similarity
neighbourhood.

## 2. Finite extremal and geometric reductions

Let \(N=pm\).  Crouzeix's finite extremal theorem (2007, Theorem
2.1) permits a maximizing scalar function for an \(N\times N\)
matrix to be chosen as a finite Blaschke product of degree at most
\(N-1\).  This is the finite Schur--Pick/Hermite reduction; it is
also restated in BGG+20 and Li's extremal-function account.  Together
with the compact unit sphere, it makes the scalar norming data a
compact finite-dimensional family.  Degree loss is included by
allowing terminal Schur parameters on the unit circle.

Apply affine normalization, a local unitary slice, and the numerical-
range Riemann pullback used in L61/L101/L115.  Any hypothetical
countersequence to (2) has a subsequence with:

1. convergent finite Schur data and norming vectors;
2. a maximal first-support winner/loser decomposition;
3. a fixed L197 disk flag after finite ramification;
4. fixed transfer and top-spectral ranks; and
5. fixed scalar product-line frames on each nonzero rank stratum.

If a support eigenvalue crosses, use L101's sequential normalized
profile estimate.  No globally differentiable largest support
eigenvector or maximizing channel is selected.

It is therefore enough to exclude one ramified real-analytic
counterarc on each of finitely many rank strata.  This is the same
finite-dimensional curve-selection use as L192 and L114, now with
the scalar Schur data included among the compact parameters.

## 3. First-support Jensen gate

Let \(J\ge0\) be L61's first support-compression Jensen gap.
If its first nonzero coefficient is positive, L61 gives a strict
negative similarity endpoint one order before any disk/normal/
elliptic response.  Since a complete similarity upper bound is
stronger than (2), that stratum cannot contain a scalar counterarc.

On \(J=0\), take the maximal common-winner copy space.  L199's
Fourier argument gives, coefficient by coefficient:

1. every true circular-normal copy coefficient is scalar on the
   winner;
2. its winner/loser cross blocks vanish; and
3. the losing compression is strict unless its common kernel is
   promoted into the winner.

Thus every scalar near-sharp branch reaches a maximal zero-Jensen
winner.  Arbitrary noncommuting normal amplitudes do not survive on
that winner; treating them by a false complete tensorization is
unnecessary.

## 4. One same-direction disk/channel reserve

Fix the actual scalar norming direction \(k\) on the top cluster and
retain

\[
 \Delta_k
 =\langle D_Xk,k\rangle+1-\sigma_X(k).             \tag{3}
\]

L325's proof before minimization gives a fixed local \(a>0\) with

\[
 4-\|f(T)x\|^2\ge a\Delta_k.                       \tag{4}
\]

L324 stacks every scalar disk/circular-normal and disk/reflected
response before completing it.  If \(d\) is the accumulated L197
loss and \(\delta\) the terminal transfer leakage, its response vector
\(R\) satisfies

\[
 \|R\|\le C(\delta+d).                             \tag{5}
\]

L327 proves on the **same** norming branch that

\[
 \delta+d\asymp\Delta_k.                           \tag{6}
\]

Let \(y\) collect the finitely many true circular-normal and marked
reflected response coordinates after L199's scalarization.  L192 and
L199 supply a retained curvature
\(\Gamma\succeq\gamma I\).  The joint part of the scalar excess is
bounded by

\[
 -a\Delta_k-\frac12y^*\Gamma y
 +2\operatorname{Re}\langle y,R\rangle.            \tag{7}
\]

One vector completion, not separate modewise estimates, gives

\[
\begin{aligned}
(7)
&\le-a\Delta_k+\frac2\gamma\|R\|^2\\
&\le-a\Delta_k+C_1\Delta_k^2
\le-\frac a2\Delta_k                               \tag{8}
\end{aligned}
\]

after shrinking.  This is the load-bearing merger.  It uses neither
L322's equality-anchor reserve nor an L197 square from another
direction, and it spends every normal curvature only once.

## 5. Retain the elliptic margin

On the pure repeated elliptic chart, L318 constructs one exact metric
\(P_{\rm ell}\) with

\[
 I\preceq P_{\rm ell}\preceq4I,\qquad
 T_{\rm ell}^*P_{\rm ell}T_{\rm ell}
 \preceq P_{\rm ell}.                              \tag{9}
\]

Its finite prepared upper endpoint satisfies

\[
 {\cal U}(c)\succeq
 \frac72\beta_L^2c^{2L}I\qquad(c\ne0),             \tag{10}
\]

after the fixed-half square completion and terminal-tail absorption.
The interior spectrum of \(4I-P_{\rm ell}\) remains uniformly
separated from zero.  Apply L325's exact three-loss identity (9) with
this prepared \(P_{\rm ell}\).  Dropping only its nonnegative
contraction and output losses gives

\[
 4-\|P_{\rm ell}^{-1/2}f(C_{\rm ell})
       P_{\rm ell}^{1/2}x\|^2
 \ge
 \langle(4P_{\rm ell}^{-1}-I)y,y\rangle.           \tag{10a}
\]

Schur congruence and the fixed interior gap make the smallest
near-sharp value of the right side uniformly comparable to
\(\lambda_{\min}{\cal U}(c)\).  Thus (10), rather than upper-endpoint
positivity by itself, is a strict scalar norm reserve on every
near-sharp input direction.

L317 puts every nonprincipal mixed transfer block one extra power of
\(c\) above its paid direct Gram.  L324 includes every corresponding
**marked** reflected/normal response in the single vector \(R\), and
L327 pays its disk/channel coefficient.  The unmarked principal
elliptic block stays inside L318's exact grouped endpoint; it is not
treated as a perturbative response.  Consequently the completion (8)
can be made while retaining a fixed fraction of (10).  L307's
\(-C^*C\) term remains inside L318's complete square; it is never
split or discarded.

At \(c=0\), (10) vanishes and (4)--(8) are the disk/channel estimate.
For \(c\ne0\), (10) and (8) are separate nonnegative margins.  No
metric from one anchor is added to a metric from another: (4) bounds
the scalar disk norm, while (9)--(10) are the exact prepared
elliptic certificate, and L324 is precisely their scalar response
comparison in the common physical chart.

## 6. The finite flag and its terminal face

Along a counterarc, apply (8) at the first active L197 quotient.
If its normalized coefficient is nonzero, (8) is strict.  If it
vanishes, L197's exact least-squares congruence descends to the next
copy kernel.  L324 retains all earlier quotient columns, so no prior
disk range reappears unpaid.  There are at most \(m\) strict disk
steps.

If the transfer leakage vanishes to the current order, L321 carries
the common eigenline to the next nonzero transfer jet.  If it
vanishes identically, the scalar channel is exact.  At an actual
scalar equality endpoint:

1. L326 applies Gau--Wu's classical theorem and splits an orthogonal
   scalar disk-model summand;
2. L328 proves that summand has dimension exactly \(p\) and lies in
   L192's single-block tube; and
3. the reducing complement is near
   \(C_p\otimes I_{m-1}\).

This gives induction on \(m\).  The base \(m=1\) is L192.  At an
exact reducing split, domain monotonicity handles the convex hull of
the block numerical ranges, and the block functional calculus norm
is their maximum.  Winner/loser cross variables have already vanished
at the zero-Jensen gate.

If \(\Delta_k=0\) but the chosen scalar function does not attain norm
two, ordinary compactness gives a strict scalar patch; Gau--Wu is
invoked only for actual equality, not from the vanishing of the
coarse defect alone.

## 7. Associated-graded remainder lemma and support crossings

The last step is an ordinary finite-dimensional lemma, but its
hypotheses must be matched to the banked endpoint results rather than
asserted from ideal membership.

On a fixed ramified stratum, use L197's analytic square coordinates
\(\eta_j\), L323's analytic off-channel frames, L199's scalar normal
coordinates, and L318's finite prepared elliptic row.  Let \(z\)
denote the joint off-channel collection.  Equations (2), (10), and
L327 give

\[
 \|z\|^2\asymp\delta+d\asymp\Delta_k.             \tag{11}
\]

Recenter all tangential variables on the channel-preserving split
model.  For the squared scalar excess

\[
 {\cal Q}(A,f,x)=\|f(\Phi_A(A))x\|^2-4\|x\|^2,
\]

the finite Taylor expansion through the first active valuation has
the form

\[
\begin{aligned}
 {\cal Q}\le {\cal Q}_{\rm split}
 &-aJ-a_0\|z\|^2-\tfrac12y^*\Gamma y
2\operatorname{Re}\langle y,R(z)\rangle\\
 &-b_0{\cal E}_{\rm ell}
o\!\left(J+\|z\|^2+\|y\|^2+{\cal E}_{\rm ell}\right),             \tag{12}
\end{aligned}
\]

where \({\cal E}_{\rm ell}\) is the actual directionwise L318 upper
margin (and is bounded below on a near-sharp top direction by a fixed
multiple of \(\lambda_{\min}{\cal U}(c)\)).  Here
\({\cal Q}_{\rm split}\le0\) by L192 and the induction hypothesis.

For clarity, the justification of every block in (12) is:

1. L61 and L199 give the \(J\) and losing-space terms and remove all
   nonscalar zero-Jensen normal coefficients.
2. L325 gives the directionwise \(-a_0\|z\|^2\) term after (11).
   The simultaneous gauge in L323--L324 removes degrees zero and one
   in \(z\), so every coefficient linear in a normal or marked
   reflected variable satisfies
   \(\|R(z)\|=O(\|z\|^2)\).
3. L192/L199 give \(\Gamma\succ0\).  All finite response coordinates
   are stacked before this one Hessian is completed.
4. L317 says that the only zero-slack elliptic quadratic blocks are
   L318's principal direct Grams.  Every other grouped block has an
   extra \(c\) and is absorbed in L318's single fixed-half block
   completion.  L318 also absorbs the complete analytic elliptic
   tail.  Thus no raw root or isolated reflected summand is being
   declared lower order.
5. After those complete blocks are removed, Taylor's theorem leaves
   only a displayed quadratic term times a vanishing chart
   coefficient.  Tangential terms with no transverse factor remain
   in \({\cal Q}_{\rm split}\).  This proves the little-\(o\) term in
   (12) on the fixed stratum.

Completing the \(y\)-block in (12), using
\(\|R(z)\|=O(\|z\|^2)\), costs \(O(\|z\|^4)\).  It is therefore paid
by a fixed fraction of the linear \(\|z\|^2\) reserve after shrinking,
exactly as in (8).  The L318 elliptic completion has already been
performed as one grouped block, so this scalar completion does not
spend its principal margin a second time.

The little-\(o\) statement is curvewise, which is all the
contradiction argument needs.  At support crossings L101 supplies the
same estimate sequentially from compact normalized profiles.  At a
disk-flag rank drop, L197 factors the new even valuation before the
graph inverse is used.  Thus no inverse spectral gap is carried across
a rank change.

Combine (8), the retained part of (10), the strict Jensen/losing
margins, and (12).  The first nonzero transverse coefficient on any
counterarc is strictly negative.  If every transverse coefficient
vanishes identically, Section 6 places the arc in an exact reducing
split controlled by L192 and the induction hypothesis.  Both cases
contradict positive scalar excess.

The normalized slice modulo the compact copy-unitary group has only
finitely many winner, disk-flag, transfer-rank, and model-split
strata.  Their local patches therefore contain one common
neighbourhood \({\cal U}_{p,m}\).  This proves (1)--(2).

## 8. What this does and does not resolve

L329 proves the scalar conjecture locally at every fixed repeated
Crabb block.  Together with L192 it covers both simple and repeated
members of this known equality family.

It does **not** prove the global Crouzeix conjecture.  A global proof
still needs either:

1. a classification showing every possible sharp sequence reduces
   to these local models; or
2. a separate uniform estimate on the compact complement of their
   neighbourhoods.

It also does not prove the completely bounded version.  L319's exact
later-Schur obstruction explains why the scalar channel argument
cannot simply be promoted to a complete endpoint.

## 9. Audit

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/repeated_crabb_scalar_local_merger.py \
  --output \
  experiments/repeated_crabb_scalar_local_merger_s70224.jsonl
```

The exact checker verifies the master completion in (7)--(8) across
multiple rational reserves, curvatures, response constants, safe
defect radii, and normal coordinates.  It audits the final algebraic
assembly only.  The geometric and analytic inputs are the cited
lemmas, not the finite computation.
