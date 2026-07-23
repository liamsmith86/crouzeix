# Exact equality on the phase-palindromic disk stratum (2026-07-23)

## 1. The theorem

Put `p=L+1>=3`.  Let `R` be the unweighted `p x p`
superdiagonal shift.  Let `H` be a positive Hermitian Toeplitz
`L x L` matrix, extended by a final zero row and column, and normalize
its diagonal to `1/2`.  Write

\[
 H_{j,j+k}=z_k,\qquad 1\leq k<L,
\]

and suppose

\[
 z_k=\omega\overline {z_{L-k}}
 \quad(1\leq k<L),\qquad |\omega|=1.                 \tag{1}
\]

As in L122, set

\[
 K=H+R^*HR,\qquad X=2K^{-1/2}HRK^{-1/2}.             \tag{2}
\]

Then `W(X)` is the closed unit disk and

\[
 \boxed{
 \sup_{\|f\|_{\mathbb D}\leq1}\|f(X)\|=2,\qquad
 t_*(X)=4.}                                           \tag{3}
\]

Thus L122's quartic null cone is not hiding a higher-order descent:
it is an exact all-size equality family for both the scalar Crouzeix
norm and the optimal similarity square.  The union over `omega` is
stratified and singular at the Crabb point, where all phase branches
meet; it is not asserted to be one smooth manifold there.

The proof gives both sides explicitly:

* a positive rank-one Stein metric of condition number four; and
* a finite Blaschke product whose value at `X` has norm two.

This settles the disk anchor, not the elliptic normal to this equality
family.  That mixed normal remains the local-merger gate.

## 2. Removing the phase

Let

\[
 D_\theta=\operatorname{diag}(1,e^{i\theta},\ldots,e^{iL\theta}).
\]

Replacing `H` by `D_theta^* H D_theta` sends
`z_k` to `e^{ik theta}z_k` and sends `omega` in (1) to
`omega e^{iL theta}`.  Choose `theta` so that the latter is one.
Under this change,

\[
\begin{aligned}
K'&=D_\theta^*KD_\theta,\\
A'&=e^{-i\theta}D_\theta^*AD_\theta,\qquad
A:=2K^{-1}HR.
\end{aligned}                                         \tag{4}
\]

Scalar rotation and unitary similarity change neither side of (3).
It therefore suffices to prove the theorem under

\[
 z_k=\overline {z_{L-k}}.                             \tag{5}
\]

All formulas below use this gauge.  The invariant formulas in
Sections 4--5 use

\[
 q=He_0,\qquad r=J\overline q,                        \tag{6}
\]

where `J` reverses the `p` coordinates.

## 3. The inverse disappears: an explicit companion form

Condition (5) implies

\[
q-r=\frac12(e_0-e_L),\qquad
Ke_0=q,\qquad Ke_L=r.                                 \tag{7}
\]

The apparently rational matrix `A=2K^{-1}HR` is actually

\[
\boxed{
\begin{aligned}
Ae_0&=0,\\
Ae_1&=2e_0,\\
Ae_j&=e_{j-1}+2z_{j-1}(e_0-e_L),
       \qquad 2\leq j\leq L .
\end{aligned}}                                        \tag{8}
\]

To verify (8), multiply its right side by `K`.  The `j=1` column gives
`2Ke_0=2He_0`.  For `j>=2`, Toeplitz structure gives

\[
 Ke_{j-1}+z_{j-1}(e_0-e_L)=2He_{j-1}.                \tag{9}
\]

The interior entries in (9) are simply two identical Toeplitz
entries.  At the first boundary the missing copy is `z_(j-1)`;
at the last boundary (5) identifies the extra copy with the same
number.  Equations (7) and (9) prove `K A=2HR`, hence (8).

In particular,

\[
\ker A=\mathbb Ce_0,\qquad
\ker A^*=\math Cr,\qquad \operatorname{rank}A=p-1.    \tag{10}
\]

The second assertion follows directly from

\[
 A^*r=A^*Ke_L=(KA)^*e_L=2R^*He_L=0.                 \tag{11}
\]

## 4. The exact Stein metric

Define

\[
 \boxed{M=K-qq^*+2rr^*.}                             \tag{12}
\]

Then

\[
 \boxed{M-A^*MA=qq^*.}                               \tag{13}
\]

Here is a direct all-size verification.  First,

\[
 A^*q=2R^*q=:u,\qquad A^*r=0.                        \tag{14}
\]

Using (8) and `KA=2HR`, entrywise Toeplitz cancellation gives

\[
 K-A^*KA=2qq^*-2rr^*-uu^*.                           \tag{15}
\]

For completeness, put `q_0=1/2`, `q_L=0`.  In the gauge (5),
`r_i=q_i` for `1<=i<L`, `r_0=0`, `r_L=1/2`, and
`u_i=2q_(i-1)` for `i>=1`.  The independent entries of the left side
of (15) reduce to

\[
\begin{array}{c|c}
\text{indices}&(K-A^*KA)_{ij}\\ \hline
i=0&2q_0\overline {q_j}\\
1\leq i,j<L&-4q_{i-1}\overline {q_{j-1}}\\
1\leq i<L,\ j=L&
 -q_i-4q_{i-1}\overline {q_{L-1}}\\
i=j=L&-\frac12-4|q_{L-1}|^2 .
\end{array}                                           \tag{16}
\]

These are exactly the corresponding entries of the right side of
(15).  Substitution of (14)--(15) into (12) proves (13).

The generalized spectrum of `(M,K)` is immediate.  Equations (7)
give

\[
\begin{aligned}
q^*K^{-1}q&=\frac12,&
r^*K^{-1}r&=\frac12,&
q^*K^{-1}r&=0.
\end{aligned}                                         \tag{17}
\]

Consequently

\[
 \operatorname{spec}(K^{-1/2}MK^{-1/2})
 =\left\{\frac12,\,
   \underbrace{1,\ldots,1}_{p-2},\,2\right\}.          \tag{18}
\]

Thus `M` is positive, (13) is a contraction certificate, and its
condition number relative to `K` is exactly four.  Explicitly, the
physical metric is

\[
 P=K^{-1/2}MK^{-1/2},\qquad P-X^*PX\succeq0,
\]

and its ordinary condition number is four.  Von Neumann's inequality
at every matrix level therefore proves `t_*(X)<=4`.

## 5. A matching finite Blaschke product

Let `xi` denote the scalar polynomial variable.  The companion form
(8) gives

\[
 \det(\xi I-A)=\xi g(\xi),\qquad
 g(\xi)=\xi^L+2\sum_{k=1}^{L-1}z_k\xi^k.             \tag{19}
\]

Define the reversed polynomial

\[
 g^\sharp(\xi)
 =\xi^L\overline {g(1/\overline\xi)}
 =1+2\sum_{k=1}^{L-1}z_k\xi^k,                       \tag{20}
\]

where the last equality uses (5).

All roots of `g` lie strictly inside the disk.  Indeed, (13) makes
`A` a contraction in the `M` metric.  Equality for a unimodular
eigenvalue would force its eigenvector to be orthogonal to
`q,A^*q,...,(A^*)^Lq`.  But direct row reduction using (8) gives

\[
 \det
 \begin{bmatrix}
 q^*\\ q^*A\\ \vdots\\ q^*A^L
 \end{bmatrix}
 =2^{L-1}\det H\ne0.                                  \tag{21}
\]

Here the determinant on the right is that of the original positive
`L x L` Toeplitz block, not of its singular `p x p` extension.
Hence there is no unitary eigenvalue.  It follows that

\[
 B(\xi)=\frac{g(\xi)}{g^\sharp(\xi)}                 \tag{22}
\]

is a finite Blaschke product and `||B||_D=1`.

Cayley--Hamilton and (10) now determine `B(A)` almost completely:

\[
 AB(A)=0,\qquad B(A)A=0.
\]

Therefore `B(A)=gamma e_0r^*` for a scalar `gamma`.  Since
`B'(0)=2z_1`, `Ae_1=2e_0`, and `A^je_1=0` for `j>=2`,

\[
 B(A)e_1=4z_1e_0.
\]

Also `r^*e_1=z_1`, so `gamma=4` whenever `z_1` is nonzero.
Continuity proves the same formula at `z_1=0`:

\[
 \boxed{B(A)=4e_0r^*.}                               \tag{23}
\]

The `K`-operator norm of this rank-one matrix is

\[
\begin{aligned}
\|B(A)\|_K
&=4\sqrt{(e_0^*Ke_0)(r^*K^{-1}r)}\\
&=4\sqrt{\frac12\frac12}=2.                           \tag{24}
\end{aligned}
\]

Since `X=K^(1/2) A K^(-1/2)`, this is exactly
`||B(X)||=2`.  Equations (22)--(24) prove the lower bound in (3);
Section 4 proves the matching upper bound.

For the original phase in (1), diagonal gauging transports (23) to
the same invariant rank-one formula after the corresponding
normalization of the Blaschke product.  In any event, its norm is
unchanged.

## 6. Consequence for the Crabb merger

L122's fourth-order Gram determinant vanishes precisely in the
tangent directions to (1).  The present theorem identifies the
nonlinear object behind that degeneracy:

\[
\{\text{positive phase-palindromic Toeplitz }H\}
\quad\longmapsto\quad
\{\text{exact disk matrices with }t_*=4\}.            \tag{25}
\]

Therefore no pure disk-flat expansion can be coercive on the whole
non-orbit disk-flat quotient at the Crabb point.  The remaining local
task is a normal-form theorem near the stratified equality family
(25): resolve its fixed-phase branches, then show that the elliptic
coordinate and the non-palindromic disk normal produce strict descent
with controlled mixed remainders uniformly through the singular apex.

## 7. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_equality.py \
  --output experiments/crabb_palindromic_equality_s70223.jsonl
```

The checker verifies symbolically through `p=6` and on exact rational
directions through `p=10`:

1. the companion identity `KA=2HR`;
2. the Stein factorization (12)--(15);
3. the three `K`-inner products (17);
4. the characteristic polynomial (19);
5. the observability determinant (21); and
6. the polynomial identity
   `g(A)=4e_0r^*g^sharp(A)`, which is (23) without a matrix inverse.

It also checks nongauged exact examples with `omega=-1` through
`p=5`, using the deliberately slower original inverse formula for
`A`.
