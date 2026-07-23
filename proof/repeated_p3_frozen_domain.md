# Frozen-domain transverse functional calculus (2026-07-22)

## 1. Setting

Let `A` be sufficiently close to `I_m tensor C_3`.  In a fixed copy
decomposition write

\[
 A=N+E,\qquad N=\operatorname{diag}(B_1,\ldots,B_m),
                                                               \tag{1}
\]

where `E` has zero diagonal copy blocks.  Put

\[
 \Omega=\operatorname{int}W(A)
\]

and let `f=phi_A:Omega->D` be any normalized Riemann map.  The point is to
compare `f(A)` and `f(N)` while keeping this **same** possibly nonsmoothly
chosen function `f`.

## 2. Exact normal certificate for the perturbed domain

Every `B_j` is a compression of `A`, so

\[
 W(B_j)\subset W(A).                                    \tag{2}
\]

Let `psi_j:D->int W(B_j)` be a normalized Riemann map and put

\[
 S_j=\psi_j^{-1}(B_j),\qquad g_j=f\circ\psi_j.          \tag{3}
\]

L73 gives a positive metric `P_j`, depending only on `B_j`, such that

\[
 I\preceq P_j\preceq4I,\qquad
 S_j^*P_jS_j\preceq P_j.                               \tag{4}
\]

By (2), `g_j` is a holomorphic self-map of the disk.  Regard `S_j` as a
contraction in the `P_j` inner product.  Von Neumann's inequality and
holomorphic functional calculus give

\[
 g_j(S_j)^*P_jg_j(S_j)\preceq P_j.                     \tag{5}
\]

Composition in the functional calculus gives

\[
 g_j(S_j)=f(B_j).                                      \tag{6}
\]

Consequently

\[
 \boxed{
 P_N=\operatorname{diag}(P_1,\ldots,P_m),\quad
 I\preceq P_N\preceq4I,\quad
 f(N)^*P_Nf(N)\preceq P_N.}                            \tag{7}
\]

This is an exact certificate for the block-diagonal normal comparison
operator using the **actual perturbed domain** `W(A)`.  Neither `P_N` nor
the proof differentiates `A->phi_A`.

## 3. Uniform transverse Cauchy expansion

After shrinking the neighbourhood once, support-function perturbation gives

\[
 \frac34\mathbb D\subset W(A),                         \tag{8}
\]

while the spectra of both `A` and `N` lie in `|z|<1/4`.  Thus the fixed
circle

\[
 \Gamma=\{|z|=1/2\}                                   \tag{9}
\]

lies in every `Omega`, surrounds both spectra, and stays a fixed positive
distance from them.  In particular,

\[
 \sup_{z\in\Gamma}\|(z-N)^{-1}\|\le C_0               \tag{10}
\]

uniformly.  Since `|f(z)|<=1` on `Omega`, Dunford calculus and the Neumann
resolvent series give, whenever `C_0||E||<1`,

\[
\begin{aligned}
 f(N+E)-f(N)
 =\sum_{k\ge1}{1\over2\pi i}\int_\Gamma
 f(z)\{(z-N)^{-1}E\}^k(z-N)^{-1}\,dz.                 \tag{11}
\end{aligned}
\]

The series converges absolutely in operator norm.  Its `k`th term is
homogeneous of degree `k` in `E`, and

\[
 \left\|\text{term}_k\right\|
 \le C_1(C_0\|E\|)^k.                                 \tag{12}
\]

In particular,

\[
 \boxed{
 \|f(A)-f(N)\|\le C\|E\|,\qquad
 \left\|f(A)-f(N)-\sum_{k=1}^p{\cal D}_k(f,N)[E^k]
 \right\|\le C_p\|E\|^{p+1},}                         \tag{13}
\]

with constants uniform over all nearby matrices and all their normalized
Riemann maps.

The same proof applies to any finite collection of Fréchet coefficients and
to their insertion into the analytic metric chart of L105.  The Riemann map
acts only as a bounded parameter on the fixed contour.

## 4. Consequence for the terminal tube

Equations (7) and (13) solve the main regularity problem behind L101:

- the exact normal comparison remains completely contractive for the actual
  enlarged numerical-range domain;
- every off-normal operator and metric remainder has an explicit transverse
  factor;
- all higher transverse coefficients have a uniform geometric majorant.

Combined with L105, this removes the possibility that the forced lower or
Stein metric coefficients diverge because of support-branch crossings.
If `g_j` is not an automorphism, disk functional calculus can create Stein
slack, so (7) need not lie in the zero-slack version of L105's chart.  L107
resolves this compatibility exactly by retaining the nonnegative Stein
Schur complement as a parameter; no tightening operation is required.

One quantitative step remains before a neighbourhood theorem.  In the
two-copy terminal chart, write `r` for the normal split and `delta` for the
Schur edge.  The desired reduced endpoint estimate is

\[
 {\cal E}_{\rm red}\preceq
 16\{\overline Q_3-\operatorname{mean}\lambda_{\max}(Q_3)I\}
 +R,\qquad
 \|R\|\le\omega(r+\delta)\,
 \delta(r^2+\delta^2),                                \tag{14}
\]

where `omega(s)->0` uniformly.  The stronger bound
`C delta(r+delta)^3` suffices.

For `r` bounded below relative to `delta`, L104's off-diagonal transfer has
an inverse of size `O(1/r)` and should center the traceless endpoint by a
quantitative implicit-function argument.  For `r` comparable to or smaller
than `delta`, one must not divide by `r`; L99/L100's pure-edge cubic gap is
uniform and absorbs the bound in (14).  Proving (14) is now the precise
two-regime analytic gate.

## 5. Scope

This is an exact functional-calculus lemma, not the terminal-tube theorem.
It does not prove (14), select the free endpoint-centering block uniformly,
or lift through L93's collapsing metric flag.  Its role is to remove the
false obstacle that nonsmooth dependence of the Riemann map on `A` prevents
transverse factorization: no such differentiation is needed.
