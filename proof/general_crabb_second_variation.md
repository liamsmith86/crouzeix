# Arbitrary-size second variation at a single Crabb block (2026-07-22)

## 1. Result and scope

Let `p=L+1>=3`, let `A=C_p` be the canonical Crabb weighted shift with weights

\[
 a_0=a_{L-1}=\sqrt2,\qquad a_j=1\quad(1\le j\le L-2),
\]

and let `P_0=diag(1,2,...,2,4)`.  For an arbitrary complex perturbation `E`, let `e_p(E)`
denote the optimum of the finite L62 second-order metric program.  Then

\[
 \boxed{e_p(E)\le0\qquad(E\in\mathbb C^{p\times p}).} \tag{1}
\]

Consequently, under the uniform conformal expansion established in L62,

\[
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(\phi_\epsilon(A+\epsilon E))-4}{\epsilon^2}
 \le e_p(E)\le0. \tag{2}
\]

This proves the single-block second-order sign in every dimension.  It does **not** prove a
full neighbourhood theorem: the quadratic form has a large kernel, and repeated-block
compression crossings require separate regularity.

## 2. General metric elimination

Put

\[
 Q=G^*P_0A+A^*P_0G.
\]

First-order complementarity leaves the first row

\[
 x_j=X_{0j},\qquad 1\le j\le L,
\]

free and determines the rest recursively:

\[
 X_{jk}=a_{j-1}a_{k-1}X_{j-1,k-1}+Q_{jk},qquad j,k\ge1. \tag{3}
\]

The analytic conformal tangent makes the terminal consistency condition `X_LL=0` hold.
The lower and upper metric Schur complements give

\[
 Y_{00}\ge\sum_{j=1}^{L-1}|x_j|^2+\frac13|x_L|^2, \tag{4}
\]

\[
 e-Y_{LL}\ge\frac13|X_{L0}|^2+rac12\sum_{j=1}^{L-1}|X_{Lj}|^2. \tag{5}
\]

Let `R` be the five forcing terms subtracted from `Y-A^*YA` in L62 equation (4), let
`d=(D_1)_{1:L,0}`, and put `C_D=R_{1:L,1:L}+dd^*`.  Equality in the Stein Schur complement can
be solved entry by entry.  Its diagonal recurrence is

\[
 Y_{LL}=4Y_{00}+2\sum_{j=1}^{L-1}(C_D)_{jj}+(C_D)_{LL}. \tag{6}
\]

Equations (3)--(6) reduce the SDP to an unconstrained quadratic in `x`.  Its quadratic part is

\[
 8\sum_{j=1}^{L-1}|x_j|^2+\frac83|x_L|^2. \tag{7}
\]

Thus the metric minimizer is unique and elementary.  The remaining issue is the sign after
inserting the conformal coefficients `G,H`.

## 3. Exact path Green function for the conformal coefficient

At normal angle `theta`, put `z=e^{i theta}`.  The normalized top support eigenvector is

\[
 v_\theta=\frac1{\sqrt L}
 \left(\frac1{\sqrt2},z,z^2,\ldots,z^{L-1},\frac{z^L}{\sqrt2}\right)^T. \tag{8}
\]

Let `V=diag(v_0(0),...,v_0(L))`, and let `Delta` be the forward incidence matrix of the path.
The support defect has the exact ground-state factorization

\[
 I-H_0(0)=\frac1{2L}V^{-1}\Delta^*\Delta V^{-1}. \tag{9}
\]

For

\[
 b_\theta=H_E(\theta)v_\theta-s(\theta)v_\theta,
\]

the right side is orthogonal to `v_theta`.  Solving the path difference equation in (9) gives
the second support coefficient without diagonalizing any trigonometric matrix:

\[
 q(\theta)=2L\sum_{\ell=0}^{L-1}
 \left|\sum_{j=0}^{\ell}\overline{(v_\theta)_j}(b_\theta)_j\right|^2. \tag{10}
\]

Formula (10), the normal correction `kappa=q-delta^2/2`, and the finite Schwarz formulas in
L62 determine `G,H` using only Laurent-polynomial arithmetic.  They are the inputs to the mode
calculation below.

## 4. Circle-mode reduction

The action

\[
 E\longmapsto e^{-i\alpha}D_\alpha^*ED_\alpha,qquad
 D_\alpha=\operatorname{diag}(1,e^{i\alpha},\ldots,e^{iL\alpha}), \tag{11}
\]

fixes `A`.  An entry `E_jm` has grade `m-j-1`.  The quadratic form `e_p` is invariant under
(11), so distinct grade pairs do not interact.

Fix `1<=k<=L-1` and put `r=L-k`.  Align the grades `+k` and `-k` by defining

\[
 u_j=\overline{E_{j,j+k+1}}\quad(0\le j<r),\qquad
 v_j=E_{j+k-1,j}\quad(0\le j\le r+1). \tag{12}
\]

For `r>=2`, define `t in C^r` by

\[
\begin{aligned}
 t_0&=u_0+\eta_kv_0+\sqrt2v_1,\\
 t_j&=u_j+v_{j+1}\qquad(1\le j\le r-2),\\
 t_{r-1}&=u_{r-1}+\sqrt2v_r+\eta_kv_{r+1},
\end{aligned} \tag{13}
\]

where `eta_1=1/sqrt(2)` and `eta_k=1` for `k>=2`.  If `r=1`, put

\[
 t_0=u_0+\sqrt2v_0+2v_1+\sqrt2v_2. \tag{14}
\]

Let `B_r` be the `r x r` canonical Crabb block, put

\[
 R_r=\frac14\left(I-\frac{B_r+B_r^*}{2}\right),\qquad
 q_{r,k}=\left(\frac{k+2}{\sqrt2},1,\ldots,1,
                    \frac{k+2}{\sqrt2}\right)^T. \tag{15}
\]

For `r>=2`, define

\[
 \mathcal K_{r,k}=
 \begin{cases}
 \displaystyle\lim_{\gamma\to\infty}
     (R_r+\gamma q_{r,k}q_{r,k}^*)^{-1},&k=1,2,\\[6pt]
 \displaystyle\left(R_r+\frac{3}{2k(k-1)(k-2)}
     q_{r,k}q_{r,k}^*\right)^{-1},&k\ge3.
 \end{cases} \tag{16}
\]

For `r=1`, set

\[
 \mathcal K_{1,k}=0\ (k=1,2),\qquad
 \mathcal K_{1,k}=\frac{k(k-1)(k-2)}{6(k+1)^2}\ (k\ge3). \tag{17}
\]

Direct substitution of (10) into the finite `G,H` formula, followed by (3)--(7), gives the
whole paired-mode contribution in the compact form

\[
 \boxed{e_p^{(k)}(E)=-t^*\mathcal K_{r,k}t\qquad(1\le k\le L-1).} \tag{18}
\]

The two unpaired bottom modes are

\[
 e_p^{(L)}(E)=-\frac{(L-1)(L-2)}{3L}
       |E_{L-1,0}+E_{L,1}|^2, \tag{19}
\]

\[
 e_p^{(L+1)}(E)=-\frac{L^2+36L-13}{6L}|E_{L0}|^2. \tag{20}
\]

For auditability, the calculation leading to (18) is only the following finite chain:
expand (10) grade by grade; take Fourier coefficients `0,...,L-1` of `kappa`; insert them in
L62 equation (16); use recurrence (3); and complete the diagonal square (7).  Multiplying the
resulting reduced block by the right side of (16) gives the identity matrix.  No spectral
root or fitted coefficient enters.

## 5. Positivity of every nonzero mode

The support matrix `(B_r+B_r^*)/2` has largest eigenvalue one, so `R_r` is positive
semidefinite with a one-dimensional kernel.  The positive vector spanning that kernel has
nonzero inner product with `q_rk`.  Therefore:

- for `k>=3`, the matrix inside the inverse in (16) is positive definite;
- for `k=1,2`, the limit in (16) exists and is positive semidefinite of rank `r-1`;
- the scalar in (17) is nonnegative;
- both coefficients in (19)--(20) are nonnegative for `L>=2`.

Hence every contribution (18)--(20) is nonpositive.

## 6. The grade-zero weighted-shift mode

The grade-zero entries are the first superdiagonal.  A perturbation supported there remains a
complex weighted shift `W_epsilon`.  Its numerical range is the centered disk of radius
`w(W_epsilon)`, so the Riemann-pulled operator is exactly

\[
 T_\epsilon=W_\epsilon/w(W_\epsilon).
\]

For normalized weights `beta_j`, the diagonal metric defined by

\[
 p_0=1,\qquad p_{j+1}=|\beta_j|^2p_j \tag{21}
\]

makes `T_epsilon` a contraction.  Near the Crabb weights the least and greatest levels remain
the two endpoints.  Moreover the numerical-radius power inequality gives

\[
 p_L=\frac{\|W_\epsilon^L\|^2}{w(W_\epsilon)^{2L}}\le4. \tag{22}
\]

Thus (21) is an exact condition-two similarity locally, not merely a second-order ansatz.  In
particular the grade-zero contribution to the L62 optimum is nonpositive.  Together with
(18)--(20) and circle-mode orthogonality, this proves (1).

## 7. Independent regeneration

`experiments/general_crabb_second_order_modes.py` independently:

1. reconstructs `G,H` from sampled support perturbations;
2. performs the general scalar metric elimination (3)--(7);
3. compares every nonzero mode with (13)--(20);
4. reconstructs full real quadratic matrices at low sizes and checks their ranks and signs.

The deterministic run tests every mode for `p=3,...,30`.  The largest absolute
formula/elimination residual is `2.34e-9` (relative residual `3.03e-13`); all low-size full
forms are negative semidefinite, with rank

\[
 p(p-2). \tag{23}
\]

The exact all-coordinate certificates L63--L64 separately prove the `p=3,4` specializations
with zero polynomial residuals.

## 8. Remaining local debt

Equation (23) leaves a real equality space of dimension `p(p+2)`.  Before going to third order,
quotient it by disk automorphisms, unitary similarity, scaling, and tangent directions to the
disk-matrix manifold.  The residual equality modes are the only ones capable of defeating a
full neighbourhood maximum.  For repeated Crabb blocks, L61 handles every positive Jensen-gap
direction; common-maximizer directions still require a matrix-valued version of the present
single-block calculation.
