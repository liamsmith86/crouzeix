# The analytic rank-one envelope and weighted transverse gate (2026-07-22)

## 1. Result and scope

Let `C=C_p`, `p=L+1>=3`, and let

\[
 T(A)=\phi_A(A)
\]

for any locally analytic normalization of the Riemann map near `C`.  There is
a real-analytic scalar certificate `Gamma_p(A)` on a full neighbourhood of
`C` such that

\[
 t_*(T(A))-4\leq\Gamma_p(A),\qquad \Gamma_p(C)=0.        \tag{1}
\]

It has three load-bearing properties:

1. its Hessian in the matrix variable is exactly L65's nonpositive
   second variation `e_p`;
2. on the complete elliptic axis,

   \[
   \Gamma_p(C+cC^*)={k(c^{2L})\over c^L}-4
   =-16c^{2L}+O(c^{4L});                                \tag{2}
   \]

3. after analytic maximization in L65's coercive directions, the remaining
   variables are exactly L66's `2p-2` real quotient coordinates.  L115
   splits these into `2p-4` disk-tangent coordinates and the one complex
   elliptic normal.

Thus L117 is not an isolated one-dimensional certificate: it lies on the
same unique analytic metric branch whose Hessian L65 already controls.
The remaining local theorem is a weighted mixed-term problem on this
finite normal form.  This note does not claim that mixed problem is closed.

## 2. Rank-one Stein coordinates

For a strictly stable matrix `T` and a vector `q`, define

\[
 P(T,q)=\sum_{n\ge0}(T^*)^nqq^*T^n.                    \tag{3}
\]

It is the unique solution of

\[
 P-T^*PT=qq^*\succeq0.                                 \tag{4}
\]

Whenever `P>0`,

\[
 t_*(T)\leq\kappa(P):={\lambda_{\max}P\over\lambda_{\min}P}. \tag{5}
\]

At `T=C`, take `q=e_0`.  Then

\[
 P(C,e_0)=P_0=\operatorname{diag}(1,2,\ldots,2,4).      \tag{6}
\]

Fix the irrelevant scale by writing

\[
 q=e_0+(0,x_1,\ldots,x_L)^T.
\]

The Stein recurrence is

\[
 P_{jk}=a_{j-1}a_{k-1}P_{j-1,k-1}+q_j\bar q_k.         \tag{7}
\]

Since the endpoint eigenvalues one and four in (6) are simple, ordinary
second-order eigenvalue perturbation gives

\[
 \boxed{\kappa(P(C,q))-4
 =8\sum_{j=1}^{L-1}|x_j|^2+\frac83|x_L|^2+O(\|x\|^3).} \tag{8}
\]

For completeness, the coefficient eight has three contributions.  For
`j<L`, the `LL` diagonal recurrence contributes `2|x_j|^2`, its coupling
to the corresponding interior level contributes another `2|x_j|^2`,
and the downward shift of the simple lower eigenvalue contributes
`4|x_j|^2`.  At `j=L`, the three terms are `1`, `1/3`, and `4/3`,
giving `8/3`.

Hence the defect-vector Hessian is positive definite in every dimension.

## 3. The analytic envelope

Near a single Crabb block, the top support eigenvalue is simple with a
uniform angular gap and its boundary is strictly convex and analytic.
The near-circle Fourier implicit-function argument used in L73 therefore
makes `T(A)` real analytic for every fixed `p`.

Apply the analytic implicit-function theorem to the critical equations

\[
 \nabla_x\kappa(P(T(A),e_0+(0,x)))=0.                   \tag{9}
\]

Equation (8) gives a unique analytic solution `x=x_*(A)` near zero.  Put

\[
 \Gamma_p(A)=\kappa(P(T(A),q_*(A)))-4.                 \tag{10}
\]

Equations (4)--(5) prove (1).  This construction is independent of the
local Riemann-map gauge: a disk automorphism transforms a rank-one Stein
defect by the standard resolvent congruence, preserves `P` and its
condition number, and the inverse automorphism gives the converse.

## 4. Identification with L65

The first- and second-order expansion of (4) is precisely L62's metric
program with equality in its Stein Schur complement.  L65's elimination
shows that the minimizing second-order program always takes this equality:
its remaining free first row has the positive quadratic part

\[
 8\sum_{j=1}^{L-1}|x_j|^2+\frac83|x_L|^2.              \tag{11}
\]

Conversely, a nearby rank-one positive defect has a unique analytic
Cholesky vector after fixing its first coordinate, so every equality
solution of the Schur recurrence is represented by (4).  Minimizing (11)
is exactly differentiating through (9).  Therefore

\[
 \boxed{\frac12D^2\Gamma_p(C)[E,E]=e_p(E).}             \tag{12}
\]

In particular its coercive rank is `p(p-2)`, with the mode kernels and
constants already proved in L65.

## 5. Exact contact with the L117 axis

L117's diagonal metric `P_c` obeys

\[
 P_c-T_c^*P_cT_c=q_cq_c^*,\qquad
 \kappa(P_c)={k(c^{2L})\over c^L}.                      \tag{13}
\]

L116 gives the same number as a lower bound for every contraction
similarity.  Thus `P_c` is globally optimal.  Since `P_c->P_0` and its
defect vector tends to `e_0`, it lies on the unique local branch (9).
This proves (2), not just an upper estimate.

## 6. Analytic splitting and the exact remaining variables

Choose a real-linear complement to the affine-unitary orbit and split it
as

\[
 E=y+z,
\]

where `y` is L65's negative eigenspace and `z` is its kernel quotient.
The Hessian in `y` is negative definite.  The analytic implicit-function
theorem therefore gives a unique local maximizer

\[
 y=y_*(z),\qquad
 H_p(z)=\Gamma_p(C+y_*(z)+z),\qquad
 \Gamma_p(C+y+z)\leq H_p(z).                            \tag{14}
\]

L66 gives `dim_R z=2p-2`.  L115 then identifies a linear splitting

\[
 z=(d,c),\qquad d\in\mathbb R^{2p-4},\quad c\in\mathbb C, \tag{15}
\]

where `d` is tangent to the exact circular-range manifold and `cC^*` is
the sole normal soft direction.  The disk theorem controls the nonlinear
anchors associated with `d`; L117 controls the physical axis `d=y=0`.
What remains is to compare these two anchors after the strong maximization
in (14).

## 7. A sufficient weighted-gradient theorem

First freeze the disk-flat variables `d=0`.  Let `g_c` be the gradient of
`Gamma_p` at `C+cC^*` projected onto the strong space.  Uniform negative
definiteness and (2) show:

\[
 \boxed{\|g_c\|=o(c^L)\quad\Longrightarrow\quad
 H_p(0,c)<0\quad(0<|c|\ll1).}                           \tag{16}
\]

Indeed, analytic strong maximization gives

\[
 y_*(0,c)=O(\|g_c\|),\qquad
 H_p(0,c)=\Gamma_p(C+cC^*)+O(\|g_c\|^2).
\]

The correction is `o(c^(2L))`, while (2) is
`-16c^(2L)+o(c^(2L))`.  Away from the maximizing graph, the negative
quadratic in `y-y_*` supplies the required absorption.

This criterion avoids any sharp constant comparison.  It is the immediate
load-bearing target before treating the additional disk-flat variables.

## 8. Exact low-size weighted evidence

The sparse checker forms

\[
 C+cC^*+s c^L Y
\]

and differentiates (10) in `s`, holding L117's optimizing defect vector as
permitted by the envelope theorem.  It derives the Riemann map by a
reduced-support-resolvent recurrence, then solves the linearized Stein
equation exactly.

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_transverse_weighted_gradient.py \
  --output experiments/crabb_transverse_weighted_gradient_s70222.jsonl
```

For every real matrix direction in sizes `p=3,4,5`, and for four additional
strong directions in sizes `p=4,5,6`, every coefficient through the
dangerous degree `2L` vanishes exactly.  The first two exposed terms are

\[
\begin{array}{c|c|c}
p&Y&D\Gamma(C+cC^*)[c^LY]\\ \hline
4&E_{30}&-16c^7+O(c^8),\\
5&(E_{30}+E_{41})/\sqrt2&-32c^{10}+O(c^{11}).
\end{array}                                             \tag{17}
\]

Thus their unscaled gradients are `-16c^4` and `-32c^6`, agreeing with
independent conformal-map finite differences and satisfying (16).
The exact cancellations also hold through degree ten for tested grade-four
and grade-six directions at `p=6`.

These checks strongly suggest the all-size path-length rule
`g_c=o(c^L)`, but a finite list is not its proof.  The next task is to
extract that cancellation directly from the sparse support/Stein
recurrences (most likely as an endpoint path-length selection rule), then
add the disk-flat variables in (15).
