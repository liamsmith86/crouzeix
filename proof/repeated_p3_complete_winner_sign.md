# Complete tied-winner second-order sign at repeated `C3` (2026-07-22)

## 1. Complete first-order normal form

Let `A_0=I_m tensor C_3` and restrict to a copy sector on which the first
top-support compression is scalar for every boundary direction.  L74 and L84
already classify, modulo block-unitary orbit, every off-diagonal edge and every
relative diagonal zero-support loop.

It remains to identify a diagonal motion shared by all copies.  Put that
common block in L69's affine-unitary normal slice

\[
 D(z,w,s,v)=
 \begin{bmatrix}
 -z/3&s&\bar z\\ w&2z/3&-s\\ v&w&-z/3
 \end{bmatrix}.                                      \tag{1}
\]

The `z` part equals `D_0(bar z)+(5z/12)I`; the scalar term is a global
translation and `D_0` is already an L84 loop.  The `s` part is `D_1(-s)` and
is also an L84 loop.  Thus the only new common directions are

\[
 W(w)=\begin{bmatrix}0&0&0\\w&0&0\\0&w&0\end{bmatrix},
 \qquad V(v)=vE_{20}.                                 \tag{2}
\]

Consequently L74, L84, and (2) exhaust the tied-winner first-order quotient;
no general diagonal block remains hidden.

## 2. First conformal motion

All graph edges and L84 loops have zero first support compression.  The common
part (2) has scalar first support

\[
 s(q)=\frac{\sqrt2}{4}(wq^{-2}+\bar wq^2)
      +\frac18(vq^{-3}+\bar vq^3).                    \tag{3}
\]

The first inverse Riemann-map coefficient is therefore

\[
 F(\zeta)=\frac{\bar w}{\sqrt2}\zeta^3
           +\frac{\bar v}{4}\zeta^4.                 \tag{4}
\]

Although `F(C_3)=0`, its Frechet derivative on the full perturbation `E` must
be retained:

\[
 \begin{aligned}
 DF(C_3)[E]={}&\frac{\bar w}{\sqrt2}
 (EC_3^2+C_3EC_3+C_3^2E)\\
 &+\frac{\bar v}{4}(C_3EC_3^2+C_3^2EC_3).             \tag{5}
 \end{aligned}
\]

Dropping (5) gives the wrong second coefficient.

## 3. Normal-angle correction

Let `Q(q)` be the complete second effective support matrix built from the
physical first perturbation, including all common, loop, and graph paths.  If
`delta(q)` is the first normal-angle shift, the second boundary datum is

\[
 \kappa(q)=\lambda_{\max}Q(q)-\frac12\delta(q)^2.      \tag{6}
\]

The two Fourier modes in (3) are orthogonal, and direct differentiation gives

\[
 \operatorname{mean}\delta^2
 =\frac94|w|^2+\frac12|v|^2.                          \tag{7}
\]

Thus, writing `lambda_0=mean lambda_max Q`,

\[
 \widehat\kappa(0)=\lambda_0-rac98|w|^2-rac14|v|^2. \tag{8}
\]

## 4. Complete endpoint identity

Use L83's pairwise edge metric and L84's diagonal loop metric.  Put every real
generator-one loop coefficient, including the common `s` motion from (1), on
the diagonal of the Hermitian copy matrix `H_1`.  Propagate the complete
second metric with (5) and the Schwarz coefficient (6), retaining its first
Fourier mode.

Before substituting (8), exact elimination gives

\[
 \mathcal E=16\overline Q-16\widehat\kappa(0)I-8H_1^2
             -18|w|^2I-\frac{37}{4}|v|^2I.            \tag{9}
\]

Using (7)--(8), the flat `w` contribution cancels exactly and (9) becomes

\[
 \boxed{
 \mathcal E=16(\overline Q-\lambda_0I)-8H_1^2
             -\frac{21}{4}|v|^2I\preceq0.}            \tag{10}
\]

Indeed, `Q(q)<=lambda_max(Q(q))I` pointwise, `H_1^2>=0`, and the last term is
nonpositive.  Therefore every fixed tied-winner first-order direction at a
repeated `C3` block has nonpositive upper second-order similarity change.

The coefficient `21/4` agrees independently with L63's single-block curvature
in the `v` direction.  Setting `w=v=0` recovers L84 exactly.

## 5. Equality and remaining scope

Equality in this explicit certificate requires

\[
 v=0,\qquad H_1x=0,\qquad
 Q(q)x=\lambda_{\max}Q(q)x\quad\text{a.e.}             \tag{11}
\]

This is a complete **second-order sign**, not yet a repeated-block
neighbourhood theorem.  Pure common flat motions are governed at higher order
by L67--L73; pure pair/star graph flats by L77/L80.  Their simultaneous
weighted equality transitions, and uniform coupling to first-order losing
sectors from L82, remain.

## 6. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_common_maximizer.py
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
```

The first checker proves completeness of the edge and loop quotients.  The
second retains independent symbolic `w,v`, three complex and three real loops,
and all nine complex edges on a full triangle.  It reconstructs `Qbar`, inserts
(5), uses (7)--(8), and verifies (10) entry by entry.  As in L83--L84, every
quadratic term is a path on at most three copy vertices, so the triangle
contains every general-multiplicity path type.
