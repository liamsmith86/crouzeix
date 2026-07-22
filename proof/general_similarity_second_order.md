# Second-order similarity reduction at a single Crabb block (2026-07-22)

## 1. Why second order is unavoidable

L61 proves that the first-order L21 similarity-square change at a repeated Crabb block is
`-4` times a support-compression Jensen gap.  For a **single** Crabb block the compression is
scalar, so that gap vanishes in every direction.  This is also the basic unresolved face of a
repeated block after a common maximizing copy-vector has been selected.

This note derives both the finite second conformal coefficient and the finite second-order
metric SDP, then validates them numerically.  The SDP's universal sign is not proved.  Thus
this is a **proved reduction plus numerical sign evidence**, not a local-neighbourhood theorem.

## 2. Expansions and first-order complementarity

Let `A=C_p`, let `P_0` be the explicit Crabb metric from L61, and work in any smooth disk gauge
for which

\[
 T_\epsilon=\phi_\epsilon(A+\epsilon E)
 =A+\epsilon G+\epsilon^2H+o(\epsilon^2). \tag{1}
\]

Use metric and bound expansions

\[
 P_\epsilon=P_0+\epsilon X+\epsilon^2Y,
 \qquad t_\epsilon=4+\epsilon^2 e. \tag{2}
\]

Write

\[
 D_0=P_0-A^*P_0A,
\]

\[
 D_1=X-A^*XA-G^*P_0A-A^*P_0G, \tag{3}
\]

\[
 \begin{aligned}
 D_2={}&Y-A^*YA-H^*P_0A-A^*P_0H-G^*P_0G\\
       &{}-G^*XA-A^*XG .
 \end{aligned} \tag{4}
\]

At multiplicity one, the L61 dual density matrix is the scalar one.  Its embedded dual blocks
are strictly positive on all three active kernels.  Complementarity therefore forces every
first-order active compression to vanish:

\[
 X_{00}=0,\qquad X_{LL}=0,\qquad (D_1)_{K_DK_D}=0, \tag{5}
\]

where `L=p-1` and `K_D` consists of levels `1,...,L`.

## 3. Second-order PSD lemma

Suppose, in kernel/range coordinates,

\[
 M_\epsilon=
 \begin{bmatrix}0&0\\0&R\end{bmatrix}
 +\epsilon M_1+\epsilon^2M_2+o(\epsilon^2),\qquad R\succ0,
\]

and `(M_1)_KK=0`.  The second-order tangent condition is

\[
 (M_2)_{KK}-(M_1)_{KR}R^{-1}(M_1)_{RK}\succeq0, \tag{6}
\]

or equivalently the affine block LMI

\[
 \begin{bmatrix}(M_2)_{KK}&(M_1)_{KR}\\
 (M_1)_{RK}&R\end{bmatrix}\succeq0. \tag{7}
\]

This is just the Schur complement.  Strict versions lift to the exact PSD constraint for all
sufficiently small positive `epsilon`; as in L61, a limiting argument handles a non-strict
second-order optimizer.

## 4. The finite second-order SDP

For the lower metric constraint, the active kernel is level zero and the positive range block
is

\[
 R_1=(P_0-I)_{\{1,\ldots,L\}}.
\]

For the upper metric constraint, the active kernel is level `L` and

\[
 R_4=(4I-P_0)_{\{0,\ldots,L-1\}}.
\]

Finally `D_0` equals one on level zero and zero on `K_D`.  Applying (7) to the three primal
LMIs gives

\[
 \begin{bmatrix}
 Y_{00}&X_{0R_1}\\X_{R_10}&R_1
 \end{bmatrix}\succeq0, \tag{8}
\]

\[
 \begin{bmatrix}
 e-Y_{LL}&-X_{LR_4}\\-X_{R_4L}&R_4
 \end{bmatrix}\succeq0, \tag{9}
\]

\[
 \begin{bmatrix}
 (D_2)_{K_DK_D}&(D_1)_{K_D0}\\
 (D_1)_{0K_D}&1
 \end{bmatrix}\succeq0. \tag{10}
\]

Minimize `e` over Hermitian `X,Y` subject to (5) and (8)--(10).  This is an ordinary finite
SDP: the apparently quadratic Schur penalties have been retained as affine block LMIs.  Given
the conformal coefficients `G,H`, no boundary discretization remains in (8)--(10).

If (1) has a uniform `o(epsilon^2)` remainder, the same strict-lift argument as L61 gives an
upper second-order Dini bound by the optimum of (8)--(10).  Equality with the actual value
coefficient would require the reverse sensitivity argument; the computations below strongly
support it but do not prove it.

## 5. Finite analytic formula for `G,H`

Put

\[
 H_0(\theta)=\operatorname{Re}(e^{-i\theta}A),\qquad
 H_E(\theta)=\operatorname{Re}(e^{-i\theta}E),
\]

and let `v_theta` be the normalized top eigenvector of `H_0(theta)`.  The top eigenvalue is
simple and uniformly separated for a single Crabb block.  Standard Hermitian perturbation
theory gives the first two numerical-range support coefficients

\[
 s(\theta)=v_\theta^*H_E(\theta)v_\theta, \tag{11}
\]

\[
 q(\theta)=v_\theta^*H_E(\theta)
 (I-H_0(\theta))^\dagger H_E(\theta)v_\theta, \tag{12}
\]

where the pseudoinverse acts on the orthogonal complement of `v_theta`.

Write the normalized conformal map as

\[
 \Psi_\epsilon(w)=w+\epsilon F(w)+\epsilon^2K(w)+o(\epsilon^2)
\]

and on the circle set `e^{-i theta}F(e^{i theta})=s(theta)+iv(theta)`.  Differentiating the
boundary tangent shows that the first normal-angle shift is

\[
 \delta(\theta)=v(\theta)-s'(\theta). \tag{13}
\]

Expanding the support line at that shifted normal gives the second normal data

\[
 \kappa(\theta)=q(\theta)-\frac12\delta(\theta)^2. \tag{14}
\]

Consequently both conformal coefficients follow from Schwarz integrals:

\[
 F(w)=\widehat s(0)w+2\sum_{k\ge1}\widehat s(k)w^{k+1},\qquad
 K(w)=\widehat\kappa(0)w+2\sum_{k\ge1}\widehat\kappa(k)w^{k+1}. \tag{15}
\]

Series inversion and analytic functional calculus now give

\[
 G=E-F(A),\qquad
 H=-DF(A)[E]+(F'F)(A)-K(A), \tag{16}
\]

where `DF(A)[E]` is the Frechet derivative.  Nilpotence makes (16) finite: `F(A)` and `K(A)`
need degrees below `p`, while `DF(A)[E]` needs only the coefficients of `F` through degree
`2p-1`, since every term is `A^r E A^s` with `r,s<p`.  Thus (11)--(16) remove numerical
Riemann-map fitting from the mathematical reduction.

## 6. Disk-gauge check

The physical-zero gauge used to derive L61 and the centroid gauge used by the Theodorsen code
need not have the same `G`.  Their maps differ by disk automorphisms.  Infinitesimally this
means

\[
 G_{\rm fit}-G_{\rm physical}=a_0I+i\gamma A-\overline{a_0}A^2. \tag{17}
\]

The script independently constructs `G_physical` from (11)--(16) and fits the difference to
`(I,A,A^2)`.  Across the 12-case run, the largest residual in (17) is
`7.15e-9`, and the largest automorphism-coefficient constraint error is `1.23e-9`.
Because `t_*` is invariant under disk automorphisms and the whole curve is fitted in one gauge,
the second-order SDP can be evaluated in either consistent pair `(G,H)`.  The analytic and
fitted-gauge SDP values agree to the solver/map accuracy reported below.

## 7. Numerical validation

`experiments/general_similarity_second_order_probe.py` constructs (11)--(16), samples each map
at `epsilon=(-2,-1,0,1,2)h`, independently fits through degree four, and solves (8)--(10) in
both gauges with Clarabel and SCS.
The main run uses `h=.005`, boundary resolutions 2048 and 4096, and support resolution 65536.
It tests `p=3,4`, three full directions and three structured tridiagonal/operator-weight
directions at each size: 12 cases and 24 cross-solver records.

Every analytic predicted second-order coefficient is negative:

\[
 -3.94535181\ \le e\le\ -0.00284309. \tag{12}
\]

The largest Clarabel/SCS difference is `3.47e-7`.  Halving the support resolution changes the
analytic `H` by at most `1.31e-16`; the independent fitted `H` changes by at most `7.09e-8`,
and the analytic/fitted-gauge SDP values agree within `3.57e-7`.  The nonlinear quotient
`(t_*(T_h)-4)/h^2` has the same negative sign in every case and differs from the prediction by
at most `1.02e-2`, consistent with an `O(h)` remainder.
Two structured `p=3` cases are nearly flat (`e=-.00284` and `-.00427`) and are the sharpest
targets for any sign proof.

A half-step rerun on the four direction-zero cases leaves the analytic predictions unchanged
and changes the independently fitted predictions by at most `3.51e-8`.  The finite nonlinear
quotients generally move toward the predictions; the nearly flat `p=3` case is already close
to the full SDP noise scale, so its finite quotient is not used as a sign certificate.

Data:

- `experiments/general_similarity_second_order_all_s9173401.jsonl` (24 main records);
- `experiments/general_similarity_second_order_halfstep_s9173401.jsonl` (eight half-step
  records).

Reproduction:

```bash
.venv/bin/python -u experiments/general_similarity_second_order_probe.py \
  --directions 3 \
  --output experiments/general_similarity_second_order_all_s9173401.jsonl
.venv/bin/python -u experiments/general_similarity_second_order_probe.py \
  --fit-step .0025 \
  --output experiments/general_similarity_second_order_halfstep_s9173401.jsonl
```

## 8. Next analytic target

Dualize (8)--(10) after substituting the finite formulas (11)--(16), and expose its dependence
on the support data as L61 did at first order.  The immediate falsification target is the nearly
flat structured `p=3` face.  A universal `e<=0`, combined with L61, would prove second-order
local nonincrease at single and repeated Crabb blocks; zero values would still need the next
nonvanishing order or an exact equality-orbit classification.
