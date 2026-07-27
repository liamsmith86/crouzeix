# The Gau--Wu endpoint residual has a signed complex structure

> **Status and scope.**  The covariance below is exact at every
> finite nondegenerate Gau--Wu model.  It identifies a universal
> complex structure in L342's two endpoint squares.  It does not yet
> prove that L351's physical shape/tracking quotient intertwines this
> structure, so phase covariance, the lower-flag equations, and the
> Hessian sign remain open.

## 1. The residual map at zero physical forcing

Retain L342's metric boundary coordinate
\(x\in{\mathbb C}^{m}\), \(m=n-1\), and write an inner tangent as
L338's

\[
 h_k=\phi k-J_fk,\qquad
 k=\sum_j{\overline v_j\over1-\overline a_jz}.     \tag{1}
\]

At zero physical direction, L342's endpoint residual map is

\[
 {\cal R}_0(x,v)=(r_-(x,v),r_+(x,v))
 \in{\mathbb C}^{m}\oplus{\mathbb C}^{m}.          \tag{2}
\]

Define the signed endpoint phase

\[
 {\cal I}_{\rm end}(u_-,u_+)=(-iu_-,iu_+).         \tag{3}
\]

Then

\[
\boxed{
 {\cal R}_0(ix,-iv)
 ={\cal I}_{\rm end}{\cal R}_0(x,v).}              \tag{4}
\]

Thus the metric boundary and zero-motion variables enter with
opposite complex orientations, while the two endpoint residuals also
carry opposite orientations.

## 2. Metric-boundary proof

For pure metric motion the first Blaschke image is zero, so

\[
 r_-=-2a_-,\qquad r_+=-2a_+.                      \tag{5}
\]

The lower endpoint motion is

\[
 a_-=-(P_0-I)_{\Pi\Pi}^{-1}Z_{\Pi0}.              \tag{6}
\]

Since \(x=e_0^*Z\Pi\), replacing \(x\) by \(ix\) replaces the
adjoint column \(Z_{\Pi0}\) by \(-iZ_{\Pi0}\).  Hence

\[
 a_-(ix)=-ia_-(x).                                \tag{7}
\]

After conjugation by the fixed Gau--Wu diagonal, L342's homogeneous
extension \(Z\) is a selfadjoint truncated Toeplitz operator.
Canonical model conjugation exchanges the endpoint columns and is
antilinear.  It therefore sends the factor \(-i\) at the lower
column to \(+i\) at the upper column.  The fixed diagonal operator
in

\[
 a_+=(4I-P_0)^{-1}Z_{\{0,\ldots,L-1\},L}          \tag{8}
\]

is complex linear, so it preserves that scalar factor.  Consequently

\[
 a_+(ix)=ia_+(x).                                 \tag{9}
\]

Equations (5), (7), and (9) prove (4) for pure metric motion.

## 3. Zero-motion proof

Replacing the zero velocity \(v\) by \(-iv\) replaces the kernel
coordinate in (1) by \(ik\).  Since \(J_f\) is antilinear, L338's
endpoint identities give

\[
\begin{aligned}
 h_{ik}(S)q&=-J_f(ik)=iJ_fk=-i\,h_k(S)q,\\
 h_{ik}(S)^*p&=-z(ik)=i\{-zk\}=i\,h_k(S)^*p.
\end{aligned}                                     \tag{10}
\]

The fixed Gau--Wu endpoint scaling is real and does not change these
phases.  Hence the lower zero residual is multiplied by \(-i\) and
the upper zero residual by \(+i\), proving (4) for pure zero motion.
Linearity proves the full identity.

## 4. Consequences and the remaining bridge

The two weights in L342's exact square gap are real diagonal:

\[
 D_-=\operatorname {diag}(1,\ldots,1,3),\qquad
 D_+=\operatorname {diag}(3,1,\ldots,1).          \tag{11}
\]

They commute with (3).  Therefore

\[
 \|r_-\|_{D_-}^2+\|r_+\|_{D_+}^2                 \tag{12}
\]

is invariant under the signed phase in (4).  Moreover, L342's
absence of an \(x\)-\(v\) cross term says that the metric and
zero-motion residual ranges are orthogonal for this weighted norm;
bijectivity makes them complementary half-dimensional ranges.

This is the universal residual-side complex structure sought after
L351.  It reduces the **mixed matching part** of the remaining phase
theorem to one precise intertwining statement:

> after eliminating L343's disk fibre, determine whether the physical
> forcing in L342 identifies L351's phase
> \((a,\tau)\mapsto(ia,-i\tau)\) with
> \({\cal I}_{\rm end}\).

There is a separate load-bearing term.  The identity
\({\cal F}_C(x)-4{\cal J}_C(h)=\|{\cal R}(C,x,h)\|^2\)
determines all terms involving \(x\) or \(h\), but it leaves a common
purely physical quadratic ``integration constant'' in
\({\cal F}\) and \(4{\cal J}\).  Even a proof of the commutative
diagram would only make the residual-square gap phase covariant; it
would not by itself make that common physical term Hermitian.
Therefore the full L351 identity still requires both:

1. the residual-forcing intertwiner for the mixed
   shape/tracking blocks; and
2. a type-\((1,1)\) formula for the common pure physical coefficient,
   most naturally from L344's fixed port.

Together those would give
\({\cal P}_\phi J_\phi^{-1}{\cal P}_\phi^T=0\).
Equation (4) proves neither remaining item by itself.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_residual_complex_covariance.py
```

The checker constructs L342's complete residual map independently,
tests both halves of (4), verifies invariance of (11), and checks
weighted orthogonality and bijectivity of the metric/zero ranges.
Two models in every dimension \(3,\ldots,8\) give metric- and
zero-covariance residuals at most \(3.50\cdot10^{-16}\) and
\(2.19\cdot10^{-16}\); weighted orthogonality is below
\(2.55\cdot10^{-16}\).  The smallest singular value of the complete
residual map is \(3.89\cdot10^{-4}\) in the most ill-conditioned
model.

The dataset SHA-256 is
`4b300f489b23ad51830014139640c1d403fc598f4a191bdb16d0751d8dc0e471`.
