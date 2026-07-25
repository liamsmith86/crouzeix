# Every analytic chain port has one common tail channel

## 1. Result (L251, 2026-07-25)

Retain L237's complete delay-chain decomposition

\[
\Xi=
\begin{bmatrix}
L_r&k_rW^*\\
We_{r-1}^*&\Xi_-
\end{bmatrix}
\quad\hbox{on}\quad
(\mathbb C^m)^r\oplus\mathcal H_- ,
\]

and write

\[
H_r(z)=zI-L_r,\qquad
D_r(z)=zI-\Xi_--\gamma_r(z,c)WW^* .
\]

Let \(\Delta_0,\ldots,\Delta_r\) be L237's continuants.  For every
function \(f\) analytic on a neighbourhood of \(\sigma(\Xi)\), define

\[
\boxed{
\mathcal Q_{r,j}^{\,f}
=\frac1{2\pi i}\oint_\Gamma
 f(z)\frac{\Delta_j(z,c)}{\Delta_r(z,c)}
 D_r(z)^{-1}\,dz,\qquad 0\le j<r .}               \tag{1}
\]

Then the two oppositely oriented direct-map blocks use the **same**
tail operator:

\[
\boxed{
\begin{aligned}
e_j^*f(\Xi)J
 &=\mu_jc^{\,r-j}W^*\mathcal Q_{r,j}^{\,f},\\
J^*f(\Xi)e_j
 &=\mathcal Q_{r,j}^{\,f}W,
\end{aligned}
\qquad
\mu_j=
\begin{cases}
2,&j=0,\\
1,&1\le j<r.
\end{cases}}                                      \tag{2}
\]

Here \(J:\mathcal H_-\to(\mathbb C^m)^r\oplus\mathcal H_-\) is the
tail inclusion.  In particular, at the remote physical endpoint,

\[
\boxed{
e_0^*f(\Xi)J=2c^rW^*\mathcal Q_{r,0}^{\,f},
\qquad
J^*f(\Xi)e_0=\mathcal Q_{r,0}^{\,f}W.}             \tag{3}
\]

Taking \(f=\phi_c\) proves that the factor two in L245 survives the
**complete** theta/Riemann functional calculus.  It is not a
resolvent-only heuristic and it cannot be detached from the common
operator \(\mathcal Q_{r,0}^{\phi_c}\).

Equation (2) also sharpens the bookkeeping in A194.  In every scalar
trace or determinant cycle, tail-to-chain and chain-to-tail crossings
occur in pairs, although a cycle may enter and leave through different
chain rows.  A one-sided term from \(\mathcal B^\sharp(x)\) therefore
cannot be assigned an independent scalar endpoint cost: it remains
inside a closed product of \(\mathcal Q\)-channels and any intervening
chain blocks.  The only exceptional port scalar is the physical
remote multiplier \(2c^r\).

L251 does **not** yet prove

\[
[c^{2k}]\log\mathcal V=4\|B_k\|_F^2.
\]

The remaining content is now a closed-channel trace identity: insert
(2) into L250's output determinant, use L244 on the common
\(\mathcal Q\)-background, and reduce the resulting paired tail
cycles by matrix-inner autocorrelation.  L249 remains the guardrail:
changing the terminal multiplier while leaving
\(\mathcal Q_{r,0}^{\phi_c}\) and the metric frozen destroys the
physical balance.

## 2. Proof

L245 gives the full block inverse

\[
(zI-\Xi)^{-1}=
\begin{bmatrix}
H_r^{-1}
+(H_r^{-1}k_r)W^*D_r^{-1}W(e_{r-1}^*H_r^{-1})
&
(H_r^{-1}k_r)W^*D_r^{-1}\\
D_r^{-1}W(e_{r-1}^*H_r^{-1})
&
D_r^{-1}
\end{bmatrix}.                                    \tag{4}
\]

Its two chain Green columns are

\[
e_{r-1}^*H_r^{-1}e_j=\frac{\Delta_j}{\Delta_r},
                                                               \tag{5}
\]

and

\[
e_j^*H_r^{-1}k_r
=
\begin{cases}
2c^r/\Delta_r,&j=0,\\
c^{r-j}\Delta_j/\Delta_r,&1\le j<r.
\end{cases}                                      \tag{6}
\]

Compress the upper-right and lower-left blocks of (4) to the
\(j\)-th chain row.  Equations (5)--(6) give

\[
\begin{aligned}
e_j^*(zI-\Xi)^{-1}J
 &=\mu_jc^{r-j}\frac{\Delta_j}{\Delta_r}
   W^*D_r^{-1},\\
J^*(zI-\Xi)^{-1}e_j
 &=D_r^{-1}W\frac{\Delta_j}{\Delta_r}.
\end{aligned}                                    \tag{7}
\]

The continuant quotient is scalar, so it may be moved through
\(W,W^*\), and \(D_r^{-1}\) without changing the order of any copy
matrix.  Apply the holomorphic functional calculus

\[
f(\Xi)=\frac1{2\pi i}\oint_\Gamma
f(z)(zI-\Xi)^{-1}\,dz .
\]

Substitution of (7) gives exactly (1)--(2).  The case \(j=0\) is
(3).  No truncation, commutativity assumption, or delayed-transfer
valuation is used.

## 3. Closed-cycle consequence

Relative to chain \(\oplus\) tail, every term of a trace-log
expansion is a closed block path.  Each visit from the tail to the
chain must therefore be matched by a later departure, possibly
through a different chain row.  Equation (2) packages a typical
paired crossing as

\[
\mu_jc^{r-j}
\operatorname {tr}\!\left(
\cdots\mathcal Q_{r,i}^{\,g}W
\mathcal C_{i,j}
W^*\mathcal Q_{r,j}^{\,f}\cdots
\right),                                         \tag{8}
\]

where \(\mathcal C_{i,j}\) denotes the intervening chain block product
and the displayed copy order is inherited from the original cycle.
Thus the apparent zero \(c\)-valuation of
\(\mathcal B^\sharp(x)\) in an isolated resolvent block is not the
valuation of a scalar volume contribution.  The scalar calculation
must retain the paired channel (8).

This consequence is deliberately weaker than the still-open
autocorrelation evaluation.  It rules out detached one-sided
bookkeeping but does not restore A192's symmetric endpoint valuation
or assert that every paired cycle vanishes.

## 4. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_analytic_chain_port_balance.py \
  --output \
  experiments/repeated_crabb_analytic_chain_port_balance_s70224.jsonl
```

The checker builds generic noncommuting tails at delays one through
six.  For a fixed polynomial functional calculus it independently:

1. forms \(f(\Xi)\) by matrix powers;
2. reconstructs every \(\mathcal Q_{r,j}^{\,f}\) by Cauchy
   quadrature from the retained Schur resolvent; and
3. verifies both identities in (2), including the unique factor two
   at \(j=0\).

The equations above, not the floating quadrature, prove L251.
The tracked dataset regenerates byte for byte with SHA-256
`11277f02a8f9270eb64b0bf1ab467f3db8ede97bf4fd16ee97419fd8f8a7f445`.
