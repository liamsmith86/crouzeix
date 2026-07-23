# Complete fixed-direction sign at repeated `C3` (2026-07-22)

## 1. Maximal winner/loser decomposition

Let `A_0=I_m tensor C_3` and consider a first-order direction on L61's
zero-Jensen face.  After the common first boundary motion is removed, the top
copy-space compression has the form

\[
 B(q)=\begin{bmatrix}0&0\\0&B_-(q)\end{bmatrix}\preceq0.                 \tag{1}
\]

Choose the winner space `W` maximally, so

\[
 \bigcap_{|q|=1}\ker B_-(q)=\{0\}.                                     \tag{2}
\]

The zero off-diagonal block in (1) is forced by negative semidefiniteness,
not assumed: a Hermitian matrix with a zero principal block and a nonzero
cross block cannot be semidefinite.  Hence every winner--loser physical block
satisfies L74's cross equation.

Continuity and (2) give

\[
 \overline B_-:=\operatorname{mean}B_-(q)\prec0.                        \tag{3}
\]

By L61's tangent duality and Slater argument, the losing lower, upper, and
Stein constraints can therefore all be made strict at first order while the
global bound coefficient remains zero.  Arbitrary losing diagonal and
losing--losing blocks are absorbed by that strict slack.

## 2. Winner--loser coefficient matrices

For every winner--loser pair, use L74's three complex quotient coefficients.
Collect generator `j` into a matrix

\[
 A_j:\mathcal L\longrightarrow\mathcal W,
 \qquad j=0,1,2.                                                       \tag{4}
\]

The two-winner/one-loser radial model is enough to polarize every quadratic
winner endpoint.  Give both winners radial derivative `+dC_3`, the loser
`-dC_3`, and use L81's optimized cross metric on both edges.  Exact
second-metric elimination yields

\[
 \begin{aligned}
 \mathcal E_{W/L}={}&16(\overline Q_{W/L}-k_0I_{\mathcal W})
 -\frac{25}{8}A_0A_0^*-8A_1A_1^*\\
 &-\frac{50}{9}A_2A_2^*,                                               \tag{5}
 \end{aligned}
\]

where `Q_W/L(q)` is the effective second support through the losing copy and
`k_0` is the conformal second mean.  The coefficient is independent of the
positive radial gap.  Every entry in general multiplicity is a polarization
of the audited two-winner/one-loser paths; additional losers sum the Gram
matrices in (5).

## 3. Add the complete winner sector

Inside `W`, L85 gives all common motions, relative loops, and winner--winner
edges.  The winner-internal and winner--loser quadratic pieces have disjoint
intermediate copy sectors, so all mixed products vanish in the winner
endpoint.  The conformal scalar is inserted only once, using the top
eigenvalue of the **total** effective support

\[
 Q(q)=Q_W(q)+Q_{W/L}(q).                                                \tag{6}
\]

Combining L85 with (5) therefore gives the complete winner endpoint

\[
 \boxed{
 \begin{aligned}
 \mathcal E_W={}&16(\overline Q-\operatorname{mean}
       \lambda_{\max}(Q)I)-8H_1^2-\frac{21}{4}|v|^2I\\
 &-\frac{25}{8}A_0A_0^*-8A_1A_1^*
  -\frac{50}{9}A_2A_2^*\preceq0.
 \end{aligned}}                                                       \tag{7}
\]

The first term is nonpositive by pointwise matrix order.  Every remaining
term is the negative of a positive semidefinite matrix.

The radial calculation proves the universal cross block in (7); replacing
the radial losing tangent by the strict tangent supplied by (3) does not
change the winner endpoint.  It only makes the losing constraints inactive
at second order.  This is precisely L82's arbitrary-gap argument, now
polarized over the whole winner space.

## 4. Fixed-direction consequence

L61 already gives strict negative first-order change whenever its Jensen gap
is positive.  On the zero-Jensen face, (7) proves nonpositive upper
second-order change.  Thus, under the uniform one-sided conformal expansion
used by L61--L62:

\[
 \boxed{\text{every fixed direction at every repeated `C3` block is
 controlled through its first two potentially active orders.}}        \tag{8}
\]

This remains an upper-certificate statement.  It is not a full repeated-block
neighbourhood theorem: when (7) has a kernel, strict higher-order terms must
be controlled uniformly as the losing mean gap (3) tends to zero.  L67--L73
and L77/L80 handle the pure boundary models, but their simultaneous weighted
transition remains.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_multiwinner_gap.py
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
```

The first checker keeps two independent winner--loser coefficient triples, an
independent winner--winner triple, an arbitrary positive gap, and an arbitrary
first conformal mode.  It verifies the strict losing compression, proves (5)
entry by entry, and independently checks the absence of mixed internal/cross
terms.  The second regenerates the complete L85 winner term.  Their block-
sector additivity is the exact finite decomposition in section 3, not a
numerical extrapolation.
