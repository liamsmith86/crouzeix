# The all-size Crabb descent gradient and pure elliptic tube (2026-07-23)

## 1. Result

Let `p=L+1>=3`, `C=C_p`, and

\[
 A_c=C+cC^*,\qquad 0<c<1.
\]

Let `Gamma_p` be L118's analytic rank-one upper envelope and put

\[
 P(z)=P_{L,c}(z)
 =2c^{L/2}T_L\!\left({z\over2\sqrt c}\right),\qquad
 r=c^L.                                                 \tag{1}
\]

For `Pi` the compression onto `span{e_0,e_L}`, define

\[
 Q_c=\Pi P(A_c)\Pi
 =\begin{pmatrix}0&2\\2r&0\end{pmatrix}.               \tag{2}
\]

Write

\[
 \tau(r)={k(r^2)\over r}.
\]

Then, for every complex matrix direction `Y`,

\[
 \boxed{
 D\Gamma_p(A_c)[Y]
 ={\tau'(r)\over2}\operatorname{Re}\left\{
 (DP(A_c)[Y])_{L0}-r(DP(A_c)[Y])_{0L}\right\}.}         \tag{3}
\]

This is an exact all-size identity, not only a leading expansion.  Moreover,
as `c->0`,

\[
 \boxed{\|D\Gamma_p(A_c)\|
 =O\!\left(c^{L+\lfloor L/2\rfloor}\right)
 =O(c^{L+1}).}                                         \tag{4}
\]

Consequently L118's weighted-gradient gate holds:

\[
 \|\operatorname{Proj}_{\rm strong}
 \nabla\Gamma_p(C+cC^*)\|=o(c^L).                      \tag{5}
\]

After analytic maximization in L65's coercive directions,

\[
 H_p(0,c)<0\qquad(0<|c|\ll1).                          \tag{6}
\]

Thus the full strong transverse tube over the **pure elliptic face** is
closed in every size.  L115's additional disk-flat variables and their mixed
merger with `c` remain; this note does not claim the complete neighbourhood
of `C_p`.

## 2. A polynomial conditional expectation

L119 proves the proper-map identity

\[
 B_c\circ\phi_c=\phi_r\circ P,\qquad
 P(\zeta+c/\zeta)=\zeta^L+r/\zeta^L.                  \tag{7}
\]

For a scalar function `h` analytic near the closed ellipse `E_c`, define
the normalized trace over the `L` polynomial fibres,

\[
 ({\cal E}h)(w)={1\over L}\sum_{P(z)=w}h(z),            \tag{8}
\]

where roots are counted with multiplicity.  Symmetry of the roots makes
this analytic across the critical values.

The Crabb/Dickson conditional-expectation identity is

\[
 \boxed{\Pi h(A_c)\Pi=({\cal E}h)(Q_c).}                \tag{9}
\]

One direct proof starts on the boundary.  Put

\[
 z_j=\zeta_j+c/\zeta_j,\qquad
 \zeta_j^L=\xi,\qquad 0\leq j<L.
\]

Averaging Laurent monomials over the `L` roots kills every exponent not
divisible by `L`.  The survivors are precisely polynomials in
`\xi+r/\xi=P(z_j)`.  Applying the same root-of-unity filter to powers of
`A_c` leaves their `0,L` compression.  Hence (9) holds for polynomials,
and uniform approximation gives the analytic statement.

At the two critical values, the multiplicities in (8) give endpoint
weights one and interior weights two.  These are exactly twice the
DCT-I weights `e_j^2`, which is the spectral version of the same proof.

## 3. The differential quadrature identity

For `|\zeta|=1`, the normalized top support vector of `A_c` is

\[
 u(\zeta)_m={d_m\zeta^m\over\sqrt L},\qquad
 d_0=d_L={1\over\sqrt2},\quad d_m=1\ (0<m<L).          \tag{10}
\]

The key simplification is that this vector is independent of `c` in the
exterior ellipse coordinate.  For `\xi=\zeta^L`, the corresponding top
support vector of the outer block (2) is

\[
 v(\xi)={e_0+\xi e_L\over\sqrt2}.                      \tag{11}
\]

Differentiate the Dickson recurrence for `P`.  A root-of-unity filter,
or direct substitution of (10), gives

\[
 \boxed{
 v(\xi)^*\Pi DP(A_c)[Y]\Pi v(\xi)
 ={1\over L}\sum_{\zeta_j^L=\xi}
 P'(z_j)\,u(\zeta_j)^*Yu(\zeta_j).}                    \tag{12}
\]

This is an algebraic identity for arbitrary complex `Y`.  To see the
filter explicitly, use

\[
 P'(z_j)
 =L\,{\zeta_j^L-r\zeta_j^{-L}\over
          \zeta_j-c\zeta_j^{-1}}                       \tag{13}
\]

and insert (10); averaging retains exactly the four outer-compression
coefficients on the left.

Geometrically, `u(\zeta)^*Yu(\zeta)` is a valid complex representative
of the first numerical-boundary displacement of `W(A_c+sY)`: its normal
component is the standard support derivative.  Multiplication by
`P'(z_j)` maps that displacement to the descended boundary.  Tangential
choices remain tangential under the conformal boundary map and do not
affect support functions.  Therefore (12) says:

> The first support variation of `W(Q_c+s Pi DP[Y] Pi)` is the average,
> over the `L` sheets, of the first variation obtained by mapping
> `W(A_c+sY)` through `P`.

## 4. Shape derivatives commute with the fibre trace

Choose analytic Riemann gauges and set

\[
 f_s=B_c\circ\phi_{A_c+sY},\qquad f_0=\phi_r\circ P.
\]

Let a dot denote the derivative at zero.  At corresponding boundary
points `z_j(s)`, the identity `|f_s(z_j(s))|=1` gives

\[
 \operatorname{Re}\left[
 \overline{\phi_r(w)}
 \{\dot f(z_j)+\phi_r'(w)P'(z_j)\dot z_j\}\right]=0,
 \qquad w=P(z_j).                                      \tag{14}
\]

Average (14) over the fibres.  By (12), the averaged mapped displacement
is exactly the numerical-boundary displacement of the two-dimensional
matrix

\[
 Q_c+s\,\Pi DP(A_c)[Y]\Pi.
\]

It follows that `E dot f` obeys the boundary condition for the first
shape derivative of the Riemann map of that `2 x 2` numerical range.
The analytic solution is unique after choosing the same disk gauge.
Indeed, the difference of two solutions divided by `phi_r` has purely
imaginary boundary values and at most one simple pole, at the zero of
`phi_r`; it is exactly an infinitesimal disk automorphism.  Such a term in
another gauge is harmless because the similarity optimum and L119's
derivative are gauge invariant.

Now differentiate the matrix function.  The fixed scalar part factors:

\[
 f_0(A)=\phi_r(P(A)).
\]

Since `P(A_c)` is block diagonal after pairing reversal coordinates,
the `0,L` compression of its Fréchet derivative depends only on
`\Pi DP[Y]\Pi`.  For the shape part, (9) gives

\[
 \Pi\dot f(A_c)\Pi=({\cal E}\dot f)(Q_c).
\]

Thus the endpoint derivative in L119 is exactly the endpoint derivative
of the normalized two-dimensional elliptic pullback along
`\Pi DP[Y]\Pi`.  This proves the descent chain rule

\[
 D\Gamma_p(A_c)[Y]
 =D\Gamma_2(Q_c)[\Pi DP(A_c)[Y]\Pi].                  \tag{15}
\]

For the matrix in (2), translation, scaling, and unitary directions do
not change the elliptic parameter to first order.  If `Z` is its
perturbation, then

\[
 \dot r={1\over2}\operatorname{Re}(Z_{L0}-rZ_{0L}).    \tag{16}
\]

Since the exact two-dimensional value is `tau(r)`, equations (15)--(16)
prove (3).

## 5. Finite-path valuation

The polynomial (1) is the Dickson polynomial characterized by

\[
 P_0(A)=2I,\qquad P_1(A)=A,\qquad
 P_n(A)=AP_{n-1}(A)-cP_{n-2}(A).                      \tag{17}
\]

Equivalently,

\[
 P_L(A)=\sum_{j=0}^{\lfloor L/2\rfloor}
 (-1)^j{L\over L-j}\binom{L-j}{j}c^jA^{L-2j}.         \tag{18}
\]

Consequently every term in `DP_L(A_c)[Y]` has the form

\[
 c^j A_c^aYA_c^b,\qquad a+b=L-2j-1.                  \tag{19}
\]

In a path contributing to a row of `e_L^*A_c^a`, at least
`ceil(a/2)` steps use the subdiagonal `cC^*`; otherwise a length-`a`
path cannot end at the upper index `L`.  Likewise every path in
`A_c^be_0` uses at least `ceil(b/2)` such steps.  Hence the `L,0`
entry of (19) is divisible by

\[
 c^{j+\lceil a/2\rceil+\lceil b/2\rceil}
 \quad\hbox{and therefore by}\quad
 c^{\lfloor L/2\rfloor}.                              \tag{20}
\]

The `0,L` entry has no negative powers of `c`, and its coefficient in
(3) is multiplied by `r=c^L`.  Thus

\[
 (DP[Y])_{L0}-r(DP[Y])_{0L}
 =O(c^{\lfloor L/2\rfloor})\|Y\|.                    \tag{21}
\]

Finally Jacobi's product gives

\[
 \tau(r)=4-16r^2+O(r^4),\qquad
 {\tau'(r)\over2}=-16r+O(r^3).                        \tag{22}
\]

Combining (3), (21), and `r=c^L` proves (4)--(5).

## 6. Strict sign after strong maximization

L118 gives

\[
 \Gamma_p(A_c)=\tau(c^L)-4
 =-16c^{2L}+O(c^{4L}).                                \tag{23}
\]

The gain from optimizing L65's strong variables is the square of the
strong gradient:

\[
 O(\|\operatorname{Proj}_{\rm strong}\nabla\Gamma_p(A_c)\|^2)
 =O(c^{2L+2\lfloor L/2\rfloor})
 =o(c^{2L}).                                           \tag{24}
\]

Equations (23)--(24) prove (6).  By L118's negative-definite analytic
splitting, all nearby strong directions away from the maximizing graph
only decrease the envelope further.

The remaining local obstruction is therefore no longer the coercive
transverse space.  It is solely the weighted merger of the elliptic
coordinate with L115's `2p-4` disk-flat coordinates.

## 7. Independent exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_descent_gradient.py \
  --output experiments/crabb_descent_gradient_s70223.jsonl
```

For `p=3,...,10`, the checker uses symbolic Dickson recurrences to
regenerate:

1. the reversal identity (11);
2. the conditional expectation (9) on monomials through degree `3L`;
3. the arbitrary-direction differential quadrature (12), with all `p^2`
   entries independent; and
4. the exact first power `c^floor(L/2)` in (21).

It then compares (3) with the independent support/Riemann/Stein engine
on three exposed coefficients:

\[
 -64c^5\quad(p=3\hbox{ axis}),\qquad
 -16c^7\quad(p=4\hbox{ bottom}),\qquad
 -32c^{10}\quad(p=5\hbox{ penultimate}).
\]

All identities are exact SymPy equalities; no floating-point conformal
map or SDP enters this audit.
