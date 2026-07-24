# The circular-normal completed square (L160, 2026-07-23)

## 1. Result and corrected geometry

Put `p=L+1>=3` and let `C=C_p`.  The coercive space in L118 must not
be optimized as one undifferentiated block over L158's Toeplitz disk
chart.  Part of that space is tangent to the much larger analytic
manifold of matrices whose numerical ranges are disks.  Maximizing in
those tangent directions cancels L152's Toeplitz quartic to first
order.  This is expected and harmless: the Berger--Okubo--Ando theorem
already gives the complete constant two on the whole disk manifold.

Use instead L115's tubular splitting:

1. an exact circular-range anchor;
2. the one complex soft normal represented by `C*`; and
3. the `2p-4` real coercive circular normals.

On L149--L158's reflected exceptional divisor, let

\[
 w_k=c^ku_k,\qquad 1\leq k\leq\lfloor L/2\rfloor ,
 \qquad r=c^L.                                      \tag{1}
\]

In the offset-one sector, the only compact quadratic coupling between
`w_1` and a coercive circular normal is with the bottom mode

\[
 y=Y_{L0}.                                          \tag{2}
\]

The joint face of L118's optimized rank-one envelope is

\[
\boxed{
 -64|w_1|^2
 +b_L\operatorname{Re}(\overline {w_1}y)
 +e_p(Y),\qquad
 b_L=-{8(5L-1)\over L}.
}                                                    \tag{3}
\]

For `L=2`, the weighted path also contains the central reflected
coordinate at the same order and the actual flat coefficient is
`-80`; discarding its extra `-16` gives the uniform upper face (3).
L65 gives

\[
 e_p(Y)\leq-B_L|y|^2,\qquad
 B_L={L^2+36L-13\over6L},                            \tag{4}
\]

with the remaining circular-normal modes contributing additional
nonpositive squares.  Completing (3)--(4) leaves

\[
\boxed{
 -R_L|w_1|^2,\qquad
 R_L=
 {32(L-1)(2L^2-L+3)\over
  L(L^2+36L-13)}>0.
}                                                    \tag{5}
\]

Thus allowing every true coercive circular normal does **not** destroy
the first reflected face.

This closes the offset-one circular-normal Schur complement.  It does
not by itself prove the full single-Crabb neighbourhood theorem.  The
higher-grade circular-normal cross faces still require an all-size
selection audit, and one still needs a uniform analytic tubular lift
over the complete disk manifold rather than only the Toeplitz
cross-section used by L158.

## 2. Why the earlier full-strong absorption target is false

Let `T_C C_disk` denote L115's circular-range tangent and let `K`
be L65's equality kernel.  The correct dimension count is in the
normal quotient, not in a Euclidean intersection with the range of
the Hessian:

\[
\begin{aligned}
\dim_{\mathbb R}(M_p/T_C{\cal C}_{\rm disk})&=2p-2,\\
\dim_{\mathbb R}\operatorname{image}
 (K\longrightarrow M_p/T_C{\cal C}_{\rm disk})&=2,\\
\dim_{\mathbb R}N_{\rm coercive}&=2p-4.              \tag{6}
\end{aligned}
\]

The two-dimensional image is the elliptic soft normal represented by
`C*`; the last line is any complementary normal section.  One must
not replace (6) by
`dim(range(-D^2 Gamma) intersect T_C C_disk)=(p-2)^2`:
that tempting identity already fails at `p=6`.

Some of L118's strong directions supply curvature corrections between
different disk charts.  A finite-difference probe on a
non-palindromic Toeplitz disk path finds

\[
 {(\hbox{full-strong completed-square gain})
  \over(\hbox{Toeplitz quartic deficit})}\longrightarrow1. \tag{7}
\]

In contrast, the gradient projected onto `N_coercive` is zero to the
precision of the independent Riemann/Stein engine.  Equation (7)
falsifies the proposed estimate that all L118 strong variables could
be absorbed by L158's prepared norm.  It does not threaten the
conjecture; it says that the exact disk manifold is the necessary
anchor, exactly as L115 anticipated.

## 3. Offset-one Fourier selection

Under L65's circle action,

\[
Y\longmapsto e^{-i\theta}D_\theta^*YD_\theta,
\]

the entry `Y_jm` has character `m-j-1`.  Insert the two offset-one
equality characters and the first ellipse reflection in L62's Fourier
recurrence.  At the first reflected weight, their only surviving
coercive circular-normal character is the unique bottom entry
`Y_L0`.  Consequently the real pairing in (3) is invariant and every
other coercive normal has zero offset-one cross coefficient.
Conjugation and circle rotation turn the real calculation below into
the complex invariant `Re(conjugate(w_1)y)`.

Finite grade-two probes show an additional power of `c`, but that is
evidence for the next selection theorem, not promoted here to an
all-grade result.

## 4. Exact bottom cross coefficient

It remains to calculate the scalar `b_L`.  Give the phase-one
offset-one equality amplitude and the ellipse parameter the same
bookkeeping weight `epsilon`.  Let `X(epsilon)` be L122's exact
physical disk chart and put

\[
 A_{\epsilon,s}
 =X(\epsilon)+\epsilon X(\epsilon)^*
  +s\epsilon^2e_Le_0^*.                              \tag{8}
\]

Let `T_{\epsilon,s}` be its normalized Riemann pullback.  At `s=0`,
solve L118's rank-one defect stationarity through its first two jets.
The envelope theorem permits that defect to be held fixed when
differentiating in `s`.  Direct expansion gives

\[
\boxed{
 \partial_s[\epsilon^4]\,
 \kappa(P(T_{\epsilon,s},q_*(\epsilon)))
 =-40+{8\over L}.
}                                                    \tag{9}
\]

Here is an inverse-free derivation of (9).  Write

\[
 H_\epsilon=H_0+\epsilon H_1,\qquad
 K_\epsilon=H_\epsilon+R^*H_\epsilon R,\qquad
 G_\epsilon=K_\epsilon^{-1/2}.
\]

If `G_j=[epsilon^j]G_epsilon`, the identity
`G_epsilon K_epsilon G_epsilon=I` gives, entrywise,

\[
(G_j)_{ab}
=-{1\over\sqrt{(K_0)_{aa}}+\sqrt{(K_0)_{bb}}}
 \sum_{\substack{i+\ell+m=j\\i,m<j}}
 (G_iK_\ell G_m)_{ab}.                              \tag{10}
\]

Thus the physical coefficients of
`X=2GHRG` require no matrix square root or inverse.  Insert them in
L62's support Green formula and the finite Schwarz recurrence for the
Riemann map.  The bottom insertion has

\[
v_\theta^*e_Le_0^*v_\theta={e^{-iL\theta}\over2L}.
\]

All Fourier modes except the bottom character vanish.  For `L>=4`,
support separation makes the two endpoint derivatives stabilize
except for the explicit bottom support weight.  Direct substitution
in the Stein and simple-eigenpair recurrences gives

\[
\begin{aligned}
\partial_s[\epsilon^4]\lambda_-
 &=16+8\sqrt2,\\
\partial_s[\epsilon^4]\lambda_+
 &=24+{8\over L}+32\sqrt2.                           \tag{11}
\end{aligned}
\]

At total weight four, the recurrence sees only bounded endpoint
neighbourhoods and the explicitly spanning endpoint entries in
`H_1` and the bottom insertion.  Substituting (10) reduces every term
to the prefix and suffix sums from L65 equation (10).  After their
cancellations, the lower sum is `16+8sqrt(2)` and the upper sum is
`24+32sqrt(2)+8/L`; the sole surviving nonconstant length dependence
is inherited from the bottom support factor `1/(2L)`.  Increasing `L`
therefore only inserts zero middle rows.  This is the all-size
stabilization step.

The strong insertion has no endpoint coefficient below total weight
four.  Hence differentiation of the ratio gives

\[
\boxed{
 [\epsilon^4s](\lambda_+/\lambda_-)
 =\partial_s[\epsilon^4](\lambda_+-4\lambda_-)
 =-{8(5L-1)\over L}.                                \tag{12}
}
\]

The short chains use the same recurrence before endpoint separation:

\[
\begin{array}{c|c|c|c}
L&\partial_s[\epsilon^4]\lambda_-
 &\partial_s[\epsilon^4]\lambda_+
 &\partial_s[\epsilon^4](\lambda_+/\lambda_-)\\ \hline
2&8(\sqrt2+2)&4(7+8\sqrt2)&-36\\
3&8(1+\sqrt2)&16(-1+6\sqrt2)/3&-112/3.
\end{array}                                         \tag{13}
\]

Equations (11)--(13) prove (9).  The endpoint weights in those sums are
`(1/sqrt(2),1,...,1,1/sqrt(2))`.

The coefficient of `s^2` in this weighted face is the ordinary L65
second variation of `e_Le_0^*`, namely `-B_Ls^2`.  L134 supplies the
flat coefficient `-64` (or the stronger `-80` when `L=2`).  Therefore

\[
\begin{aligned}
-64|w_1|^2
 +b_L\operatorname{Re}(\overline {w_1}y)-B_L|y|^2
\leq
-\left(64-{b_L^2\over4B_L}\right)|w_1|^2,
\end{aligned}                                       \tag{14}
\]

and simplifying the parenthesis gives (5).

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_circular_normal_face.py \
  --output experiments/crabb_circular_normal_face_s70223.jsonl
```

The checker independently:

1. constructs `K^(-1/2)` from (10);
2. reconstructs the physical disk/ellipse path (8);
3. derives the Riemann pullback from the support recurrence;
4. solves both rank-one defect jets exactly;
5. verifies the two endpoint derivatives and their cancellation to
   (9);
6. checks every real coercive normal character through `p=5`; and
7. compares the completed residual with (5).

The default exact audit covers `p=3,...,6`; longer exact runs through
`p=8` give the same rational formula.  The checker is a regeneration
of the displayed all-size endpoint calculation, not the reason for
extrapolating it.
