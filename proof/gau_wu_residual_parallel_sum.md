# The Gau--Wu endpoint gap is one exact parallel-sum square

> **Status and scope.**  The parallel-sum identity and optimal
> tracking formula below are exact at every finite nondegenerate
> Gau--Wu model.  They use L354's exact residual normal form and do
> not assume its numerically observed physical shape factorization.
> The common physical Hessian and its sign remain open.

## 1. Conjugate endpoint weights

Retain L354's antiunitary endpoint restriction
\[
 \Gamma:p^\perp\longrightarrow q^\perp.
 \tag{1}
\]
The two weights in L342's exact square gap are
\[
\begin{aligned}
 D_-&=\operatorname {diag}(1,\ldots,1,3)
      &&\text{on }p^\perp,\\
 D_+&=\operatorname {diag}(3,1,\ldots,1)
      &&\text{on }q^\perp.
\end{aligned}                                      \tag{2}
\]
Canonical endpoint exchange gives
\[
\boxed{\Gamma^*D_+\Gamma=D_-.}                    \tag{3}
\]
Equivalently, (3) is the exact reason that L352's metric and
zero-motion residual graphs are orthogonal for the weighted endpoint
norm.

For fixed physical forcing and zero motion, varying L342's free
metric boundary row changes the residual pair by
\[
 (u,-\Gamma u),\qquad u\in p^\perp.               \tag{4}
\]
The pure metric map onto \(u\) is bijective: its domain and codomain
have the same real dimension, and L342's strictly positive
homogeneous metric form makes its kernel zero.

## 2. Parallel-sum elimination

Let \(a={\cal A}_\phi(C)\), let \(\tau\) be L351's tracking
coordinate, and retain L354's exact relation
\[
 r_+=-\Gamma r_-+
 {\cal T}_\phi\overline\tau+\Delta_\phi(C).        \tag{5}
\]
Put
\[
 w={\cal T}_\phi\overline\tau+\Delta_\phi(C)
 \in q^\perp.                                     \tag{6}
\]
For fixed \(C,\tau\), (4)--(6) turn L342's complete residual square
into
\[
 \|u\|_{D_-}^2+\|w-\Gamma u\|_{D_+}^2.            \tag{7}
\]
Set \(y=\Gamma u\).  Equation (3) gives
\[
\begin{aligned}
 \|u\|_{D_-}^2+\|w-\Gamma u\|_{D_+}^2
 &=\|y\|_{D_+}^2+\|w-y\|_{D_+}^2\\
 &=2\left\|y-\frac w2\right\|_{D_+}^2
   +\frac12\|w\|_{D_+}^2.                         \tag{8}
\end{aligned}
\]
Surjectivity in (4) permits \(y=w/2\).  Therefore
\[
\boxed{
 \min_x\{
 {\cal F}_C(x)-4{\cal J}_C(\tau)\}
 =\frac12\|
 {\cal T}_\phi\overline\tau+\Delta_\phi(C)
 \|_{D_+}^2.}                                     \tag{9}
\]

If
\[
 e_\phi(C)=\min_x{\cal F}_C(x)
 \tag{10}
\]
is L342's common physical similarity Hessian, then (9) is
equivalently the complete joint scalar normal form
\[
\boxed{
 e_\phi(C)-4{\cal J}_C(\tau)
 =\frac12\|
 {\cal T}_\phi\overline\tau+\Delta_\phi(C)
 \|_{D_+}^2.}                                     \tag{11}
\]
No semidefinite variable, moving metric, or opaque real Schur
complement remains in (11).

## 3. Exact optimal tracking

L354 defines \({\cal T}_\phi\) by composing L351's invertible
pure-zero tracking map with L338's map
\[
 k\longmapsto-2zk.
\]
Multiplication by \(z\) maps \({\cal K}_f\) isomorphically onto
\(q^\perp=z{\cal K}_f\).  Thus
\({\cal T}_\phi\) is a complex isomorphism.

Equation (11) has the unique zero-motion maximizer
\[
\boxed{
 \overline{\tau_{\rm opt}(C)}
 =-{\cal T}_\phi^{-1}\Delta_\phi(C).}              \tag{12}
\]
The endpoint component of L354's corrector is
\(\overline{\tau_0}\), while
\(\Delta_\phi(C)\in{\cal M}\).  Hence (12) also gives
\[
 \boxed{(\tau_{\rm opt})_0=0,}                    \tag{13}
\]
recovering the endpoint-trace Euler equation directly from the
residual geometry.

At each fixed physical direction \(C\), equation (12) makes the
zero-optimized tracking law
\(\tau_{\rm opt}(C)=R_\phi\overline{{\cal A}_\phi(C)}\)
for a fixed complex-linear \(R_\phi\) exactly equivalent to L354's
remaining physical statement
\[
 \Delta_\phi(C)={\cal B}_\phi a
 \quad\text{with }{\cal B}_\phi\text{ complex linear}. \tag{14}
\]
Indeed, (12) and invertibility of \({\cal T}_\phi\) give both
directions of this fixed-physical equivalence.  In particular, (14)
implies the antiholomorphic tracking law observed in L351 after
also optimizing the physical disk fibre at fixed shape.  That
fully optimized law alone does **not** imply (14): it constrains
\(\Delta_\phi\) only on the selected physical lift, not on every
direction in \(\ker{\cal A}_\phi\).  Disk-fibre annihilation and
complex linearity of the descended response remain separate proof
debts.  See `audit_l352_l355_20260905.md`, Section 2, for an exact
quadratic countermodel to the stronger converse.

## 4. What still remains

Equation (11) completely resolves the metric/zero portion of the
joint Hessian.  It does not determine the integration constant
\(e_\phi(C)\).  The live phase proof still requires both:

1. prove (14) from L344's fixed Hardy port; and
2. prove that \(e_\phi(C)\), after L343's disk elimination, is
   type \((1,1)\) in \(a\).

Only after both steps is the common shape Hessian Hermitian.  Its
negative sign then remains to be represented by a positive Gram.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_residual_parallel_sum.py
```

The checker independently rebuilds L342's complete residual square,
verifies (3), eliminates the metric block directly, and compares the
result with the right side of (9).  It then compares (11) with the
independently assembled scalar and similarity Hessians and verifies
the optimizer equation in (12).

Two models in every dimension \(3,\ldots,8\) give weight,
parallel-sum, scalar-gap, and optimizer-equation residuals within the
bounds \(1.60\cdot10^{-14}\), \(4.23\cdot10^{-15}\),
\(4.46\cdot10^{-15}\), and \(6.46\cdot10^{-14}\), respectively.

The dataset SHA-256 is
`2c99ea653b2b5ddf671db159f115d4e07ac42ca12b1209df429154e5a2797e1c`.
