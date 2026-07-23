# Normal scalar-support collisions split into two-copy pairs (2026-07-22)

## 1. The collision stratum

The remaining boundary of L111 is the normal case in L93's scalar-support
relations

\[
 Z^2=\alpha I,\qquad ZZ^*+Z^*Z=\beta I.              \tag{1}
\]

After a phase rotation and scaling, a non-scalar normal solution is

\[
 Z_0=sJ,\qquad
 J=\begin{bmatrix}I_p&0\\0&-I_q\end{bmatrix},
 \qquad s>0.                                         \tag{2}
\]

The multiplicities `p,q` may be arbitrary.  This note proves that every
first-order equality-preserving perturbation of (2), modulo a copy unitary,
is a rectangular direct sum of the two-copy terminal chart already closed
by L110.

## 2. Exact relation tangent

Write

\[
 X=\begin{bmatrix}A&B\\C&D\end{bmatrix}.
\]

The linearized square relation is

\[
 Z_0X+XZ_0=\dot\alpha I.
\]

Its diagonal blocks and vanishing cross blocks give

\[
 A=aI_p,\qquad D=-aI_q,                              \tag{3}
\]

for one complex scalar `a`; the rectangular blocks `B,C` are unrestricted.
The linearized anticommutator relation is

\[
 XZ_0+Z_0X^*+X^*Z_0+Z_0X=\dot\beta I.               \tag{4}
\]

Its off-diagonal blocks cancel identically.  With (3), both diagonal blocks
equal `4s Re(a)`, so (4) imposes no further condition and fixes
`dot beta=4s Re(a)`.

Thus the full real tangent dimension is

\[
 2+4pq.                                              \tag{5}
\]

The two scalar dimensions move the common `+/-` normal eigenvalues; all
apparent multiplicity freedom is in `B,C`.

## 3. Unitary gauge and singular-value splitting

For a rectangular `L:p -> q`, use the skew-Hermitian generator

\[
 K=\begin{bmatrix}0&L\\-L^*&0\end{bmatrix}.
\]

Its orbit tangent is

\[
 [K,Z_0]
 =\begin{bmatrix}0&-2sL\\-2sL^*&0\end{bmatrix}.      \tag{6}
\]

Choose `L=C^*/(2s)`.  Then

\[
 X+[K,Z_0]
 =\begin{bmatrix}
 aI_p&B-C^*\\0&-aI_q
 \end{bmatrix}.                                      \tag{7}
\]

Hence the only unitary-invariant coupling is the arbitrary rectangular
matrix

\[
 R=B-C^*.                                            \tag{8}
\]

Take a singular-value decomposition `R=U Sigma V^*`.  The block-diagonal
copy unitary `diag(U,V)` fixes `Z_0` and turns (7) into an orthogonal direct
sum of

\[
 \begin{bmatrix}s+ta&t\sigma_j\\0&-s-ta\end{bmatrix},
                                                               \tag{9}
\]

one for every nonzero singular value, together with unmatched
one-dimensional normal coordinates.  Every matrix in (9) is precisely
L100/L110's trace-zero two-copy Schur chart; the scalar `a` is its normal
weighted tangent.  The unmatched coordinates remain on L88's exact normal
stratum.

## 4. Uniform consequence

The off-diagonal unitary orbit in (6) has real dimension `2pq`; the quotient
has the expected `2pq` real singular-edge coordinates plus the two common
normal parameters.  More importantly, (7)--(9) are exact linear algebra,
not a dimension-only argument.

Therefore a weighted sequence approaching a normal full-common-top block
has only two equality-compatible behaviours:

1. it stays in the normal stratum, covered exactly by L88; or
2. after a moving copy unitary, its first transverse coefficient is a
   direct sum of independent L110 two-copy terminal pairs.

Any failure of (3)--(4) produces a non-scalar first variation of L93's
effective support and hence a strict Jensen/metric-flag coefficient.
Together with L111, this closes every **full common-top** multiplicity
collision.  The remaining arbitrary-copy gate is a proper common-top space
whose descending L93 flag changes rank, followed by the L86 losing-space
lift.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_normal_collision.py
```

For every `1<=p,q<=4`, the checker constructs the full complex tangent and
proves by exact real-linear rank that

\[
 \dim T_{\rm full}=4pq+2,\quad
 \dim T_{\rm hom}=4pq,\quad
 \dim T_{\rm orbit}=2pq.
\]

It also keeps independent symbolic entries in a `2 x 3` rectangular pair
and proves (7) entry by entry.  The final SVD is the standard exact
rectangular singular-value theorem and is dimension-independent.
