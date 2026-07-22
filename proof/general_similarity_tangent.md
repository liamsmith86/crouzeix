# First-order similarity reduction at repeated Crabb blocks (2026-07-22)

## 1. Scope and status

This note linearizes the general L21 contraction-similarity program at a repeated Crabb
equality block.  It gives a finite conformal tangent, solves the resulting semidefinite tangent
program in closed form, and proves that its value is nonpositive in every direction.  Under a
uniform first-order conformal expansion, a direct feasible-metric lift then proves that the
upper right directional derivative of `t_*` is nonpositive.

This is a **local first-order theorem**, not a neighbourhood theorem and not a proof of
Crouzeix's conjecture.  Directions in the equality case of the final Jensen inequality may
still increase at second order.

## 2. Repeated disk equality block

Let `C_p` be the `p x p` Crabb weighted shift, with weights

\[
 (\sqrt2,\sqrt2) \quad(p=3),\qquad
 (\sqrt2,1,\sqrt2) \quad(p=4).
\]

Thus `W(C_p)=D`.  The experiments use either rotated copy-major sums
`A_0=direct_sum(e^{i gamma_j} C_p)` or the level-major form `C_p tensor I_m`; these are
unitarily equivalent.  Put

\[
 A_\epsilon=A_0+\epsilon E,\qquad \epsilon\downarrow0.
\]

The perturbations have zero scalar part.  Scaling `A_epsilon` and its numerical range together
does not alter the pulled-back operator.  A different normalization of the Riemann map only
postcomposes it by a disk automorphism.  The value `t_*` is invariant under such
postcomposition: a metric making `T` a contraction also makes every disk automorphism of `T`
a contraction, and the inverse automorphism gives the reverse inequality.

## 3. Directional support function

For

\[
 H_0(\theta)=\operatorname{Re}(e^{-i\theta}A_0),\qquad
 H_E(\theta)=\operatorname{Re}(e^{-i\theta}E),
\]

the top eigenvalue of `H_0(theta)` is one, has multiplicity `m`, and is separated from the
rest.  If the columns of `V(theta)` span that eigenspace, ordinary Hermitian degenerate
perturbation theory gives the exact right derivative

\[
 s_E(\theta)=\lambda_{\max}\!\left(V(\theta)^*H_E(\theta)V(\theta)\right). \tag{1}
\]

Consequently the numerical-range support function is

\[
 h_\epsilon(\theta)=1+\epsilon s_E(\theta)+o(\epsilon). \tag{2}
\]

If the largest eigenvalue in (1) stays simple with a uniform gap, the remainder is uniform and
`s_E` is smooth.  Crossings require a nonsmooth shape-derivative treatment; the one-sided
formula (1) itself remains valid.

## 4. Finite conformal tangent

Use the disk gauge in which the conformal map `Psi_epsilon:D -> W(A_epsilon)` fixes zero and
has positive derivative.  Write

\[
 \Psi_\epsilon(w)=w+\epsilon F_E(w)+o(\epsilon).
\]

On the unit circle its normal displacement must equal (1):

\[
 \operatorname{Re}\bigl(e^{-i\theta}F_E(e^{i\theta})\bigr)=s_E(\theta).
\]

The Schwarz integral therefore gives, with
`s_hat(k)=(2 pi)^{-1} integral s_E(theta)e^{-ik theta} dtheta`,

\[
 F_E(w)=\widehat s(0)w+2\sum_{k\ge1}\widehat s(k)w^{k+1}. \tag{3}
\]

Inverting the map and applying analytic functional calculus yields

\[
 \phi_\epsilon(A_\epsilon)
 =A_0+\epsilon G_E+o(\epsilon),\qquad
 G_E=E-F_E(A_0). \tag{4}
\]

Because `A_0^p=0`, (3) truncates **exactly** in (4):

\[
 F_E(A_0)=\widehat s(0)A_0
   +2\sum_{k=1}^{p-2}\widehat s(k)A_0^{k+1}. \tag{5}
\]

Thus the entire first-order conformal calculation uses only `p-1` Fourier coefficients of the
scalar function (1), even though the boundary perturbation has infinitely many modes.

## 5. Explicit equality metric and tangent SDP

For Crabb weights `w_0,...,w_{p-2}`, define

\[
 P_0=\operatorname{diag}(1,|w_0|^2,|w_0w_1|^2,\ldots).
\]

Write its scalar level entries as `p_j`, starting with `p_0=1`.

The level values are `(1,2,4)` for `p=3` and `(1,2,2,4)` for `p=4`; repeat this metric in the
same copy-major or level-major ordering as `A_0`.  Then

\[
 I\preceq P_0\preceq4I,\qquad A_0^*P_0A_0\preceq P_0, \tag{6}
\]

and `t_*(A_0)=4`.  The reverse inequality follows from `||A_0^(p-1)||=2`: any
similarity to a contraction must have condition number at least two.

Set

\[
 P=P_0+\epsilon X+o(\epsilon),\qquad t=4+\epsilon d+o(\epsilon).
\]

Let `U_1,U_4,U_D` be orthonormal bases of the kernels of
`P_0-I`, `4I-P_0`, and `D_0=P_0-A_0^*P_0A_0`.  Linearizing the three primal LMIs gives

\[
 U_1^*XU_1\succeq0,\qquad U_4^*(dI-X)U_4\succeq0, \tag{7}
\]

\[
 U_D^*\left[X-A_0^*XA_0-G_E^*P_0A_0-A_0^*P_0G_E\right]U_D\succeq0. \tag{8}
\]

Minimizing `d` subject to (7)--(8) is the finite tangent-cone SDP.  These kernel compressions
are the exact first-order tangent conditions for the PSD cone.

## 6. Closed-form dual and universal sign

Work in the level-major form `A_0=C_p tensor I_m`, put `L=p-1`, and denote the positive Crabb
weights by `w_0,...,w_{L-1}`.  The three active kernels are respectively level zero, level
`L`, and levels `1,...,L`.  If `W_1,W_4,W_D` are the embedded PSD dual variables, dual
stationarity is

\[
 W_4+A_0W_DA_0^*=W_1+W_D,\qquad \operatorname{tr}W_4=1. \tag{9}
\]

The off-diagonal block recurrences in (9) force `W_D` to be level diagonal.  All diagonal
blocks are determined by one density matrix `Y`:

\[
 (W_D)_{jj}=\frac4{p_j}Y\quad(1\le j\le L),\qquad
 W_4=Y,\quad W_1=4Y,\quad Y\succeq0,\quad\operatorname{tr}Y=1. \tag{10}
\]

Conversely every density matrix `Y` gives a dual feasible triple, so there are no missing dual
faces.  If `(G_E)_{j,j+1}` denotes an adjacent level block, strong duality gives

\[
 d_E=\lambda_{\max}\!\left[
  4\sum_{j=0}^{L-1}\frac{(G_E)_{j,j+1}+(G_E)_{j,j+1}^*}{w_j}
 \right]. \tag{11}
\]

Only the linear term `s_hat(0)A_0` in (5) reaches the first superdiagonal.  Hence, with

\[
 M_E=\sum_{j=0}^{L-1}\frac{E_{j,j+1}+E_{j,j+1}^*}{w_j},
\]

equation (11) becomes

\[
 d_E=4\lambda_{\max}(M_E)-8L\widehat s(0). \tag{12}
\]

The normalized top eigenvector `r` of `Re C_p` has endpoint entries `(2L)^{-1/2}` and
interior entries `L^{-1/2}`.  In particular

\[
 w_jr_jr_{j+1}=\frac1L. \tag{13}
\]

Taking the zeroth Fourier coefficient of the compressed perturbation in (1) and using (13)
gives

\[
 \frac1{2\pi}\int_0^{2\pi}V(\theta)^*H_E(\theta)V(\theta)\,d\theta
 =\frac1{2L}M_E. \tag{14}
\]

Convexity of the maximum eigenvalue (equivalently, test a top eigenvector of the mean against
each pointwise maximum) now yields

\[
 \widehat s(0)=\frac1{2\pi}\int\lambda_{\max}(V^*H_EV)
 \ge \frac1{2L}\lambda_{\max}(M_E). \tag{15}
\]

Combining (12)--(15) proves the exact formula

\[
 \boxed{d_E=-4J_E\le0},\qquad
 J_E:=2L\widehat s(0)-\lambda_{\max}(M_E)\ge0. \tag{16}
\]

This argument holds for every multiplicity and every perturbation, not only for the sampled
families.  It also does not require the compression maximum in (1) to be simple; simplicity is
needed only for the clean conformal expansion used in the next step.

## 7. Feasible-metric lift

The tangent SDP is strictly feasible after allowing a sufficiently large `d`: choose a scalar
level-diagonal `X` with `X_00>0`, successive inequalities
`X_jj>w_{j-1}^2 X_{j-1,j-1}`, and `d>X_LL`, then scale it to dominate the fixed `G_E` terms.
Thus an optimizer of (7)--(8) can be perturbed by an arbitrarily small multiple of a strict
tangent direction.

Assume now that (4) holds uniformly in operator norm.  For a strict tangent solution, the
active kernel blocks of all three exact primal LMIs are positive at order `epsilon`; their
inactive blocks remain uniformly positive, and the off-diagonal Schur-complement correction is
only order `epsilon^2`.  Therefore it lifts to an exact feasible metric for all sufficiently
small positive `epsilon`.  Letting the added strict perturbation tend to zero proves

\[
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(\phi_\epsilon(A_\epsilon))-4}{\epsilon}
 \le d_E=-4J_E\le0. \tag{17}
\]

If `J_E>0`, the similarity square is strictly below four for all sufficiently small positive
`epsilon`.  If `J_E=0`, (17) controls only first order.  Proving a full neighbourhood theorem
therefore requires classifying the Jensen-equality directions and treating their next
nonvanishing variation.

Equality in (15) has a precise elementary characterization.  Let
`B_E(theta)=V(theta)^*H_E(theta)V(theta)` and choose a top eigenvector `y` of its mean.  The
continuous function

\[
 \lambda_{\max}(B_E(\theta))-y^*B_E(\theta)y
\]

is nonnegative and has integral zero exactly when `J_E=0`.  Hence

\[
 J_E=0\quad\Longleftrightarrow\quad
 \text{there is one unit }y\text{ that is a top eigenvector of }B_E(\theta)
 \text{ for every }\theta. \tag{18}
\]

For a single Crabb block (`m=1`) this is automatic, so all variation begins at second order.
For repeated blocks every splitting with no common maximizing copy-vector has `J_E>0`;
equation (18), rather than all perturbations, is the exceptional set needing the next
calculation.  A dimension/genericity statement has not yet been proved.

## 8. Numerical cross-check

`experiments/general_similarity_tangent_probe.py` constructs (1)--(8).  A 65,536-angle run on
the two earlier ultralocal cross directions gave:

| `p,m` | Clarabel `d` | SCS `d` | nonlinear `(t_*(1e-5)-4)/1e-5` |
|---|---:|---:|---:|
| `3,2` | -1.5570293183 | -1.5570293187 | -1.557011 |
| `4,2` | -1.6369727114 | -1.6369726798 | -1.636959 |

A wider 16,384-angle run tested 12 cases: `p=3,4`, `m=2,3`, and one independently seeded
full, cross-copy, and operator-weight direction for each pair.  Both solvers found a negative
derivative in every case, ranging from `-1.5570293` to `-4.0337326`.  Their maximum objective
difference was `3.17e-8`, and their maximum difference from the independent closed form (16)
was `3.14e-8`; the most negative reported tangent slack was `-7.56e-9`.

Comparing each SCS prediction with the earlier full nonlinear value at `epsilon=1e-4` gives a
maximum relative slope error `1.20e-4` (maximum absolute error `4.78e-4`).  Recomputing the
two cross cases at 65,536 rather than 16,384 angles changes the objective by at most
`3.04e-11`.  These agreements validate the implementation and the first-order model; the
universal sign itself follows from (15), not from the samples.

Data:

- `experiments/general_similarity_tangent_s9173401.jsonl` (four high-resolution records);
- `experiments/general_similarity_tangent_all_s9173401.jsonl` (24 cross-solver records).

Reproduction:

```bash
.venv/bin/python -u experiments/general_similarity_tangent_probe.py \
  --output experiments/general_similarity_tangent_s9173401.jsonl
.venv/bin/python -u experiments/general_similarity_tangent_probe.py \
  --block-sizes 3 4 --multiplicities 2 3 \
  --families full cross operator_weight --support-resolution 16384 \
  --output experiments/general_similarity_tangent_all_s9173401.jsonl
```

## 9. Remaining proof problem

Equation (18) classifies first-order equality, but it does not make every equality direction an
exact symmetry.  The live local question is to quotient the common-maximizer directions by
the equality-preserving orbits (unitary similarity, scalar affine changes, disk automorphisms,
and repeated-block motion), then compute second order on any transverse remainder.  The
single-block case necessarily lies entirely in this second-order problem.  The separate
regularity debt is to establish (4) for nonsmooth compression-eigenvalue crossings rather than
assuming the uniform conformal expansion.

## 10. Literature calibration

Greenbaum--Lewis--Overton (2017) already proved nonnegative directional derivative of the
scalar Crouzeix ratio at a Crabb/monomial pair and identified second order on an active
manifold as the next local question.  Lewis--Overton (2020) proved partial smoothness of the
numerical radius at single superdiagonal disk matrices.  Those results motivate the present
equality classification but do not cover repeated support multiplicity or the stronger,
domain-recalibrated completely bounded similarity quantity `t_*`.  The density-matrix dual
collapse and Jensen formula (16) therefore appear new after this targeted audit; a broader
publication-level novelty search is still required.
