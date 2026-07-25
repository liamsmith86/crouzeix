# Every complete delay has one scalar continued-fraction reflection

## 1. Result (L237, 2026-07-24)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

and assume a complete delay of length \(r\geq1\):

\[
B_1=\cdots=B_r=0,\qquad B_j=W^*(S^*)^jV.           \tag{1}
\]

Put \(W_j=S^jW\), remove

\[
W_0\mathbb C^m\oplus\cdots\oplus W_{r-1}\mathbb C^m,
\]

and let \(J:{\cal H}_r\to{\cal H}\) include the retained space.  As
usual,

\[
S_r=J^*SJ,\qquad V_r=J^*V,\qquad W_r=J^*S^rW,
\]

\[
E_r=V_rV_r^*,\qquad F_r=W_rW_r^*.
\]

L235's one-delay arrow has the following exact arbitrary-delay
version.  Let \(N_r\) be the forward block shift on
\((\mathbb C^m)^r\), let \(F_0\) project onto its first block, and put

\[
L_r(c)=N_r+c(I+F_0)N_r^*.                           \tag{2}
\]

The balanced ellipse pencil

\[
\Xi_S(c)=S+c(I+F)S^*(I+E)
\]

has the block form

\[
\boxed{
\Xi_S(c)=
\begin{bmatrix}
L_r(c)&c(I+F_0)e_{r-1}W_r^*\\
W_re_{r-1}^*&\Xi_{r,-}(c)
\end{bmatrix},}                                    \tag{3}
\]

where

\[
\Xi_{r,-}(c)=S_r+cS_r^*(I+E_r).                    \tag{4}
\]

Define scalar continuants

\[
\Delta_0=1,\qquad \Delta_1=z,\qquad
\Delta_2=z^2-2c,
\]

\[
\Delta_j=z\Delta_{j-1}-c\Delta_{j-2}\quad(j\geq3),
\]

and

\[
\gamma_1(z,c)=\frac{2c}{z},\qquad
\gamma_r(z,c)=c\frac{\Delta_{r-1}(z,c)}
                         {\Delta_r(z,c)}
\quad(r\geq2).                                     \tag{5}
\]

Then the entire removed chain has the exact retained resolvent

\[
\boxed{
J^*(zI-\Xi_S(c))^{-1}J
=\{zI-\Xi_{r,-}(c)-\gamma_r(z,c)F_r\}^{-1}.}        \tag{6}
\]

Thus arbitrary delay growth is not a growing noncommutative word
problem: all left-terminal excursions are resummed by the one scalar
continued fraction \(\gamma_r\).

There is a sharper associated-graded fact.  Let

\[
\eta(z,c)=\frac{z-\sqrt{z^2-4c}}2
=\frac cz\,C\!\left(\frac c{z^2}\right),            \tag{7}
\]

where \(C(t)=\sum_{j\geq0}C_jt^j\) is the Catalan series.  Then

\[
\boxed{
\gamma_r(z,c)-\eta(z,c)
=\frac{c^r}{z^{2r-1}}+O(c^{r+1}).}                 \tag{8}
\]

The coefficient in (8) is exactly one for every \(r\).  Hence the
effect of the doubled original endpoint first reaches the promoted
tail after exactly \(r\) reflected steps, with no coefficient growth.

Finally, L221's retained boundary metric is

\[
\boxed{
\begin{aligned}
J^*P_{\rm bl}^S(c)J
={}&I-\sum_{n\geq1}\frac{q^n}{1+q^n}
   (S_r^*)^nE_rS_r^n\\
&+q^r\sum_{n\geq0}q^nS_r^nF_r(S_r^*)^n,
\qquad q=c^2 .
\end{aligned}}                                     \tag{9}
\]

Equations (6), (8), and (9) give an arbitrary-grade replacement for
the expanding cubic-through-sextic word lists.

L237 does **not** yet prove L228's slack anticommutator.  The remaining
content is to show that the first left reflection in (8), paired with
the first right-boundary reflection and the Schur square, is exactly

\[
E_1F_r+F_rE_1,\qquad E_1=S^*ES,                    \tag{10}
\]

at order \(c^{2r+2}\).  This is now a two-boundary coefficient
calculation with a unit left-reflection coefficient, rather than an
unbounded terminal-crossing expansion.

## 2. Delay-chain block form

Under (1), L209 makes

\[
W_0,\ldots,W_r
\]

isometric, mutually orthogonal frames.  Relative to the removed
chain followed by the retained space,

\[
S=
\begin{bmatrix}
N_r&0\\
W_re_{r-1}^*&S_r
\end{bmatrix},\qquad
E=\begin{bmatrix}0&0\\0&E_r\end{bmatrix},
\qquad
F=\begin{bmatrix}F_0&0\\0&0\end{bmatrix}.          \tag{11}
\]

Taking the adjoint in (11) and multiplying by the two endpoint
weights gives

\[
(I+F)S^*(I+E)=
\begin{bmatrix}
(I+F_0)N_r^*&(I+F_0)e_{r-1}W_r^*\\
0&S_r^*(I+E_r)
\end{bmatrix}.
\]

This proves (3).

## 3. Eliminate the clean chain exactly

The matrix \(zI-L_r(c)\) is scalar block tridiagonal.  Its first
reverse edge has weight \(2c\); all later reverse edges have weight
\(c\).  Expanding its determinant along the last row gives precisely
the recurrence for \(\Delta_r\).

For \(r\geq2\), Cramer's rule yields

\[
e_{r-1}^*(zI-L_r(c))^{-1}e_{r-1}
=\frac{\Delta_{r-1}}{\Delta_r}.                    \tag{12}
\]

For \(r=1\), the factor \(I+F_0=2I\) sits in the coupling block
instead.  In both cases the Schur complement of \(zI-L_r(c)\) in
\(zI-\Xi_S(c)\) is

\[
zI-\Xi_{r,-}(c)-\gamma_r(z,c)F_r.
\]

The retained block of the inverse is the inverse of that Schur
complement, proving (5)--(6).  Holomorphic functional calculus may
now be applied to (6) on a sufficiently large fixed contour.

## 4. The endpoint reflection has coefficient one

Set \(t=c/z^2\) and write

\[
\gamma_r(z,c)=\frac cz\,h_r(t).
\]

The last-corner Schur complement gives the scalar iteration

\[
h_1(t)=2,\qquad
h_r(t)=\frac1{1-th_{r-1}(t)}\quad(r\geq2).          \tag{13}
\]

The Catalan series is the normalized fixed point

\[
C(t)=\frac1{1-tC(t)},\qquad C(0)=1.                \tag{14}
\]

Subtracting (14) from (13) gives the exact error recursion

\[
h_r-C
=\frac{t(h_{r-1}-C)}
       {(1-th_{r-1})(1-tC)}.                       \tag{15}
\]

Since

\[
h_1-C=1+O(t),
\]

induction in (15) proves

\[
h_r(t)-C(t)=t^{r-1}+O(t^r).
\]

Multiplication by \(c/z\) proves (8), including its coefficient.
This is an all-\(r\) scalar proof, not a fitted continuation from the
finite-grade jets.

## 5. Boundary metric and remaining gate

Formula (9) is L221 with all \(r\) removed layers at once.  The
right-defect orbits compress without a shift.  The left orbits
\(W_0,\ldots,W_{r-1}\) disappear, while

\[
S^{r+n}W=JS_r^nW_r,
\]

which supplies the factor \(q^r\).

The comparison object \(\eta\) in (7) is the unreflected half-line
self-energy.  Equation (8) therefore isolates the first left
boundary image.  A proof of L228 should now:

1. insert the full L125 theta/ODE scalar coefficients, not merely
   their Newton bottom edge;
2. subtract the half-line zero/one-reflection cancellation;
3. keep the single coefficient in (8);
4. pair it with the right defect \(E_r\); and
5. include the right-defect Schur square.

Only the fifth step can change the raw two-boundary cross into the
anticommutator (10).  No further delay-chain expansion is needed.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_multidelay_terminal_resolvent.py \
  --output \
  experiments/repeated_crabb_multidelay_terminal_resolvent_s70224.jsonl
```

The six deterministic records use delays one through six, defect
multiplicities one through three, and noncommuting tails.  They audit
(3), (5)--(6), (9), every earlier transfer zero, and the exact Catalan
prefix in (8).  The tracked SHA-256 is
`78718693beb33f65a82c5f255a8f998ed7599edb4a80cf1b6ae7625adf87bf65`.
