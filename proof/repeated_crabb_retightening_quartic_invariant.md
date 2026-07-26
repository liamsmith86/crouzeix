# The moving quartic has an exact two-ended transfer signature

## 1. Result (L302, 2026-07-26)

Retain L298--L301's balanced pure partial-isometry colligation

\[
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\]

and put

\[
B=W^*S^*V.                                       \tag{1}
\]

Let

\[
\begin{aligned}
A(c)&=S+cA_1+c^2A_2+O(c^3),\\
D(c)&=V+cD_1+c^2D_2+O(c^3)
\end{aligned}                                    \tag{2}
\]

be the moving operator and canonical defect frame from L227.  Insert
L298's grade-one fixed-base retightening direction
\(({\mathscr R},{\mathscr F})\).  L299's first remaining defect is

\[
{\cal N}_3
=-\{A_1^*{\mathscr R}S+S^*{\mathscr R}A_1
     +D_1{\mathscr F}^*+{\mathscr F}D_1^*\}.      \tag{3}
\]

Use the convenient metric-only cubic cancellation

\[
X_3=-{\cal G}_S({\cal N}_3),\qquad C_3=0.         \tag{4}
\]

This representative need not be used in the final metric.  L301
proves that changing it changes the next endpoint only by a bounded
Markov response.  Its lower corner is selection-independent as well:
for L301's successor

\[
{\cal T}(X,C)
=-\{A_1^*XS+S^*XA_1+D_1C^*+CD_1^*\},
\]

one has \(V^*{\cal T}(X,C)V=0\), using
\(SV=0\), \(V^*S^*=0\), and
\(D_1=-2(S^*)^2V\).  Hence both endpoint classes computed below are
independent of the chosen cubic preimage.

The quartic defect left after (4) is

\[
\begin{aligned}
{\cal N}_4=-\{&
A_2^*{\mathscr R}S+S^*{\mathscr R}A_2
+A_1^*{\mathscr R}A_1\\
&+D_2{\mathscr F}^*+{\mathscr F}D_2^*
+{\mathscr F}{\mathscr F}^*\\
&+A_1^*X_3S+S^*X_3A_1\}.                         \tag{5}
\end{aligned}
\]

Set

\[
Q_4=W^*{\cal G}_S({\cal N}_4)W.                  \tag{6}
\]

For every L206 endpoint-cokernel test \(Y=Y^*\), whose observability
Gramian commutes with \(S\),

\[
\boxed{
\langle Y,Q_4\rangle
=\frac{15}{4}\operatorname {tr}(YBB^*)
 +\frac12\operatorname {tr}\{Y(BB^*)^2\}.}        \tag{7}
\]

The same quartic state forcing has the exact lower corner

\[
\boxed{
K_4:=V^*{\cal N}_4V
=6B^*B+\frac94(B^*B)^2\succeq0.}                 \tag{8}
\]

Put

\[
P_4=\frac{15}{4}BB^*+\frac12(BB^*)^2.            \tag{9}
\]

Then

\[
\boxed{Q_4-P_4\in\operatorname {ran}{\cal M}_S.}\tag{10}
\]

There is an important two-ended consequence.  The metric-only
quartic cancellation \(-{\cal G}_S({\cal N}_4)\) has lower motion
\(-K_4\), so it **cannot** be inserted into a lower-tight metric
unchanged.  The compulsory parallel correction

\[
{\cal G}_S(VK_4V^*)                               \tag{11}
\]

neutralizes that lower corner.  Its upper commutant class is

\[
6BB^*+\frac94(BB^*)^2.
\]

Combining this with the negative of (9) leaves the exact
lower-neutral upper class

\[
\boxed{
R_4=\frac94BB^*+\frac74(BB^*)^2\succeq0.}         \tag{12}
\]

This is a **cost**, not a favorable negative face.  Crucially, it is
a pure prior-\(B\) factor:

\[
U^*R_4U=0\qquad\text{whenever }U^*B=0.            \tag{13}
\]

Thus the first nonlinear even correction has no pointwise obstruction
on the next surviving upper flag after its lower endpoint is restored.
On the already active \(B\)-range it is a higher-order cost beneath
the strict quadratic margin.

L302 closes the quartic **two-ended pointwise homology class**.  It
does not yet:

1. supply a uniformly bounded analytic column realizing (10) and
   (12) through rank changes; or
2. combine this Stein-response signature with every nonlinear
   moving Schur-graph term in the complete metric endpoint.

Those are the next ordered-flag/Smith and mixed-graph gates.  In
particular, (7) alone must not be advertised as a completed
lower-tight quartic correction.

## 2. The grade-one direction is finite

L298 gives

\[
\begin{aligned}
{\mathscr F}
&=-\frac12VB^*B-\frac12(I-E)SWB\\
&=-\frac12SWB,                                    \tag{14}
\end{aligned}
\]

because \(VB^*=ESW\).  Therefore

\[
{\mathscr R}
={\cal G}_S(V{\mathscr F}^*+{\mathscr F}V^*).     \tag{15}
\]

In fact this Green sum has only two nonzero terms.  If

\[
H=V{\mathscr F}^*+{\mathscr F}V^*,
\]

then exact partial-isometry reduction gives

\[
(S^*)^2HS^2=0,\qquad
\boxed{{\mathscr R}=H+S^*HS.}                    \tag{16}
\]

Hence every part of (5) is a finite word polynomial except the
metric-only Green representative \(X_3\).

## 3. The cubic Green term telescopes in the trace

Put

\[
K=SA_1^*+A_1S^*.                                 \tag{17}
\]

The first elliptic tangent satisfies the exact identity

\[
K=2\{S^2+(S^*)^2-S^3S^*-S(S^*)^3\}.             \tag{18}
\]

Therefore

\[
Y_2=2\{S^2+(S^*)^2\}
\]

solves the dual Stein equation

\[
\boxed{Y_2-SY_2S^*=K.}                           \tag{19}
\]

Let \(H_Y\) be L206's observability Gramian.  If
\(H_YS=SH_Y\), it commutes with every coefficient in (2).  Cyclicity,
Stein adjointness, (4), and (19) give

\[
\begin{aligned}
&\operatorname {tr}H_Y
 \{-A_1^*X_3S-S^*X_3A_1\}\\
&\qquad
=\operatorname {tr}{\cal G}_S({\cal N}_3)H_YK\\
&\qquad
=\operatorname {tr}{\cal N}_3H_YY_2.             \tag{20}
\end{aligned}
\]

Thus the only apparently nonlocal term in (5) becomes a finite
polynomial inside every relevant pairing.

## 4. Exact cyclic-word certificate

Substitute (14)--(16), the exact canonical \(D_1,D_2\), and (20)
into \(\operatorname {tr}{\cal N}_4\).  Reduce using only

\[
SS^*S=S,\qquad S^*SS^*=S^*,\qquad EF=FE=0.       \tag{21}
\]

After cyclic trace collection, both the quartic expression and the
state lift of (9) have exactly the following five classes:

\[
\begin{array}{c|r}
\text{cyclic word}&\text{coefficient}\\ \hline
I&17/4\\
S^*S&-17/4\\
(S^*)^2S^2&-17/4\\
(S^*)^3S^3&15/4\\
\{(S^*)^3S^3\}^2&1/2 .
\end{array}                                      \tag{22}
\]

Their difference has zero exact cyclic classes.  Since purity gives

\[
\operatorname {tr}\{W^*{\cal G}_S(Z)W\}
=\operatorname {tr}Z,                            \tag{23}
\]

the scalar case of (7) follows:

\[
\boxed{
\operatorname {tr}Q_4
=\frac{15}{4}\|B\|_F^2
 +\frac12\operatorname {tr}\{(B^*B)^2\}.}         \tag{24}
\]

The lower corner needs no cyclic reduction and no Green sum.  Since
\(SE=0\) and \(ES^*=0\), the \(A_2\) and \(X_3\) cross terms in
(5) vanish between two copies of \(E\).  Direct word reduction gives

\[
\boxed{
E{\cal N}_4E
=6\,ESFS^*E+\frac94(ESFS^*E)^2
=VK_4V^*.}                                       \tag{25}
\]

This proves (8) and exposes why the attractive upper sign in (7)
cannot be used without paying attention to the lower endpoint.

## 5. Weighted commutant identity

Take the spectral resolution

\[
H_Y=\sum_\lambda\lambda\Pi_\lambda.               \tag{26}
\]

L206 proves that every \(\Pi_\lambda\) reduces \(S,S^*,E,F\) and
both defect spaces.  The canonical coefficients, finite formulas
(14)--(16), and both Stein inverses in (4)--(6) therefore restrict
to each reducing block.  Apply (24) on every block, multiply by
\(\lambda\), and sum.  The transfer intertwining

\[
YB=BA,\qquad H_YV=VA
\]

then turns the two block sums into

\[
\frac{15}{4}\operatorname {tr}(YBB^*)
\quad\text{and}\quad
\frac12\operatorname {tr}\{Y(BB^*)^2\}.
\]

This proves (7).  Both terms are nonnegative for \(Y\succeq0\), and
their sum vanishes exactly on the reducing support where \(B=0\).

Finally, L206 identifies these \(Y\)'s with
\((\operatorname {ran}{\cal M}_S)^\perp\).  Equation (7) says
\(Q_4-P_4\) annihilates that orthogonal complement, proving (10) by
finite-dimensional Fredholm duality.

For the lower neutralization, Stein adjointness and \(H_YV=VA\)
give

\[
\begin{aligned}
\langle Y,W^*{\cal G}_S(VK_4V^*)W\rangle
&=\operatorname {tr}(AK_4)\\
&=6\operatorname {tr}(YBB^*)
  +\frac94\operatorname {tr}\{Y(BB^*)^2\}.       \tag{27}
\end{aligned}
\]

Subtract (7) from (27).  The result is the pairing with \(R_4\) in
(12).  Hence the lower-neutral correction

\[
-{\cal G}_S({\cal N}_4)+{\cal G}_S(VK_4V^*)
\]

is endpoint-equivalent, after a perpendicular Markov response, to
the upper motion \(R_4\).  Equation (13) is immediate from
\(U^*B=0\).

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_retightening_quartic_invariant.py \
  --output \
  experiments/repeated_crabb_retightening_quartic_invariant_s70226.jsonl
```

The checker reconstructs (3)--(8) independently from the moving
canonical metric and defect factor.  Exact rational-word arithmetic
checks (16), Hermiticity of (3), the telescope (19), the lower
corner (25), and the zero-class residual in (22).  The 22 numerical
records include
unstructured colligations, rank-changing chains, complete first
delays, and nonscalar weighted reducing direct sums.  The tracked
dataset has SHA-256

```text
c8e3b4edfadcb95d8d313f20e3d1056865983df8f51ac254d1d9974d524a6816
```
