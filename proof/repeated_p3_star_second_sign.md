# Arbitrary-multiplicity star sign at repeated `C3` (2026-07-22)

## 1. Star-coupled repeated block

Let

\[
 A_0=I_m\otimes C_3,
\]

choose the L61 common maximizing copy vector `y=e_0`, and put
`r=m-1`.  Consider a perturbation with zero diagonal and
orthogonal-to-orthogonal copy blocks.  For every `xi_k in y^perp`, its two
selected/cross blocks are the L74 representatives with coefficients

\[
 (\alpha_{0k},\alpha_{1k},\alpha_{2k})\in\mathbb C^3.
\]

Collect coefficients of each generator type into vectors

\[
 a_j=(\alpha_{j1},\ldots,\alpha_{jr})^T\in\mathbb C^r,
 \qquad j=0,1,2.                                      \tag{1}
\]

Then the upper second-order Dini change of the L21 similarity square is
nonpositive for every `m` and every triple (1).  This extends L76 from one
orthogonal copy to simultaneous star coupling with arbitrarily many copies.

## 2. Effective support mean

Let `Q(q)` be the `m x m` second effective top-support matrix,

\[
 Q(q)=P_{\rm top}V(q)R(q)V(q)P_{\rm top},             \tag{2}
\]

and put

\[
 \kappa(q)=\lambda_{\max}Q(q),\qquad
 k_0=\frac1{2\pi}\int_0^{2\pi}\kappa(e^{i\theta})\,d\theta. \tag{3}
\]

The first support compression vanishes by L74.  Exact Laurent averaging in
(2) gives

\[
 \overline Q:=\frac1{2\pi}\int Q(e^{i\theta})\,d\theta
 =\begin{bmatrix}\mu&0\\0&G\end{bmatrix},            \tag{4}
\]

where

\[
\begin{aligned}
G={}&\frac5{128}\bar a_0a_0^T
    +\frac14\bar a_1a_1^T
    +\frac5{72}\bar a_2a_2^T,\\
\mu={}&\operatorname{tr}G.
                                                               \tag{5}
\end{aligned}
\]

The selected-copy entry is the sum of the pairwise means, while the
`y^perp` block is their polarized Gram matrix.  There are no mean cross
entries between those two sectors.

Unlike the multiplicity-two case, `kappa` need not be `pi`-periodic when
`m>2`; its first Fourier coefficient can be nonzero.  The second inverse-map
term evaluated at the repeated Crabb block therefore has the general form

\[
 K_2(A_0)=k_0A_0+2\widehat\kappa(1)A_0^2.             \tag{6}
\]

The `A_0^2` term must be retained in the metric calculation.  It cancels
exactly from the final endpoint below.

## 3. Full second metric and endpoint identity

For the cross pair indexed by `k`, place the L76 first metric block

\[
 Z(a^{(k)})=\begin{bmatrix}
0&-3\sqrt2\alpha_{0k}/8&0\\
\alpha_{2k}/\sqrt2&-2\sqrt2\alpha_{1k}
                         &3\sqrt2\alpha_{0k}/4\\
0&-\sqrt2\alpha_{2k}&0
\end{bmatrix}                                         \tag{7}
\]

between `y` and `xi_k`, and use its adjoint in the reverse block.  Their sum
is a Hermitian first metric tangent whose three active compressions vanish.

Insert the lower and Stein Schur penalties, then propagate the second metric
through the two positive Crabb levels exactly as in L76.  The upper copy-space
endpoint is

\[
 \boxed{
 {\cal E}=16(\overline Q-k_0I_m)-8B,
 \qquad
 B=\begin{bmatrix}
 \|a_1\|^2&0\\0&\bar a_1a_1^T
 \end{bmatrix}\succeq0.}                              \tag{8}
\]

Identity (8) includes every mixed product between distinct orthogonal copies.
The contribution from `widehat kappa(1)` in (6) is exactly zero at the upper
endpoint; omitting it before checking this cancellation would be invalid.

## 4. Sign

Pointwise Hermitian order gives

\[
 Q(q)\preceq\kappa(q)I_m.
\]

Averaging and using (3) yields

\[
 \overline Q\preceq k_0I_m.                            \tag{9}
\]

Both terms in (8) are therefore negative semidefinite.  Choosing the scalar
metric endpoint to be `lambda_max(cal E)` and using the strict-lift/limiting
argument of L62 proves

\[
 \boxed{
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(T_\epsilon)-4}{\epsilon^2}
 \le\lambda_{\max}({\cal E})\le0.}                    \tag{10}
\]

This is a matrix Jensen argument at the **effective support level**, not an
assumption that pairwise certificates can simply be added.

## 5. Exact equality classification

The star structure makes `Q(q)` block diagonal:

\[
 Q(q)=\operatorname{diag}(s(q),R(q)),\qquad R(q)\succeq0. \tag{11}
\]

Exact Laurent algebra gives the trace-reversal identity

\[
 s(q)=\operatorname{tr}R(-q).                          \tag{12}
\]

If `a_1=0`, the odd part in (12) vanishes and

\[
 s(q)=\operatorname{tr}R(q)\ge\lambda_{\max}R(q).     \tag{13}
\]

Thus the selected copy vector is a common top eigenvector of `Q(q)`, belongs
to the kernel of `B`, and makes (8) equal to zero.  Hence
`lambda_max(cal E)=0` whenever `a_1=0`.

Conversely suppose `a_1` is nonzero and `cal E` has a null vector.  Since both
terms in (8) are NSD, write that unit vector as `(0,z)` with
`a_1^Tz=0`; it must also be a common top eigenvector of `Q(q)`.  Apply this at
`q` and `-q`.  From (12),

\[
 z^*R(q)z\ge\operatorname{tr}R(-q),\qquad
 z^*R(-q)z\ge\operatorname{tr}R(q).                   \tag{14}
\]

Each left side is at most the trace of its own PSD matrix.  Adding (14)
forces equality everywhere.  Therefore `R(q)` has rank one with fixed range
`span{z}`.  Its mean `G` also has that range, so (5) forces every `a_j` to be
proportional to `conj(z)`.  In particular `a_1=c conj(z)`, while
`a_1^Tz=c||z||^2=0`; hence `a_1=0`, a contradiction.  We have proved

\[
 \boxed{\lambda_{\max}({\cal E})=0\iff a_1=0}          \tag{15}
\]

for nonzero or zero star data alike.

There is one further split inside the explicit endpoint equality (15).  A
fixed orthogonal-copy vector joins
the selected top eigenspace exactly when

\[
 \dim\operatorname{span}\{a_0,a_2\}\le1.             \tag{16}
\]

Indeed, equality in the trace bound (13) is equivalent to fixed rank-one
range for `R(q)`, and its mean range is the span of `conj(a_0),conj(a_2)`.
The rank-one case (16) is unitarily a single cross pair and is closed cubically
by L77.  Genuine rank-two data require at least two orthogonal copies, and the
selected support branch is their only common top branch.

This does **not** say that the optimal second-order similarity coefficient
vanishes.  The metric tangent (7) is only one feasible choice.  A nonlinear
Riemann-map/SDP check on the canonical orthogonal data
`a_0=(1,0), a_2=(0,1)` gives quotients `(t_*-4)/epsilon^2` tending to about
`-.566`, not zero.  Thus additional homogeneous first-metric tangents improve
(8), and the full multi-copy second-order metric optimization remains active.
The earlier inference that rank-two data first move at fourth order was false.

## 6. Scope and regeneration

The theorem covers arbitrary multiplicity but only the pure star quotient:
there are no diagonal perturbation blocks and no blocks internal to
`y^perp`.  Those pieces, their mixtures with (1), and the common-maximizer
order inequality remain before a repeated-block neighbourhood theorem.

Run

```bash
.venv/bin/python -u experiments/repeated_p3_star_second_sign.py
```

The checker uses two independent symbolic orthogonal copies.  It reconstructs
the complete support mean (4)--(5), keeps an arbitrary complex first conformal
mode in (6), builds the full second metric, and proves (8) entry by entry.
It also proves (12) and the flat trace identity (13) exactly.
Every entry for general `r` is a sum of diagonal terms or a polarization of
one pair of coefficient vectors, so the two-symbolic-copy identity proves the
dimension-independent formula.
