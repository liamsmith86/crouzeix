# Independent audit of L17 and L21 (2026-07-21)

This note re-derives the two load-bearing lemmas without using the confidence
labels in their original proof notes.  The purpose is adversarial: identify
the exact hypotheses and check the two steps most likely to hide a reversed
inequality or an unjustified symmetry reduction.

## 1. L17: Schwarzian comparison and the positive-solution step

Let $g$ be real $C^3$ on $[-D,D]$, with $g'>0$, $g(0)=0$, and
$g'(0)=1$.  Put $Q=Sg/2$.  The standard quotient representation is

\[
 y_1={g\over\sqrt{g'}},\qquad y_2={1\over\sqrt{g'}},\qquad
 y_j''+Qy_j=0. \tag{A1}
\]

The normalizations give $y_1=s$ and $y_2=c+ps$, where
$s(0)=0,s'(0)=1,c(0)=1,c'(0)=0$.  Therefore the unknown $p$ cancels:

\[
 {1\over g(D)}-{1\over g(-D)}
 ={c(D)\over s(D)}-{c(-D)\over s(-D)}. \tag{A2}
\]

The point requiring an audit is the variational interpretation of (A2).
For $h\in H^1_0(-D,D)$ with $h(0)=1$, minimize

\[
 E_Q[h]=\int_{-D}^D(h'^2-Qh^2). \tag{A3}
\]

The positive solution $y_2=1/\sqrt{g'}$ gives the exact ground-state
identity

\[
 E_Q[h]=\int_{-D}^D y_2^2\left((h/y_2)'\right)^2\ge0. \tag{A4}
\]

The boundary term vanishes because $h(\pm D)=0$.  Thus the form is
positive and the affine minimization problem is well posed.  Its minimizer
solves $h''+Qh=0$ separately on the two half intervals.  Integration by
parts gives

\[
 \min E_Q=h'(0-)-h'(0+)
 ={c(D)\over s(D)}-{c(-D)\over s(-D)}, \tag{A5}
\]

including the signs in (A2).  The positive solution also supplies
disconjugacy, so the denominators in (A2)--(A5) do not vanish.

If $Q\ge-1$, then, for every admissible $h$,

\[
 E_Q[h]\le E_{-1}[h]. \tag{A6}
\]

The direction is correct: a larger potential $Q$ makes the integrand
$h'^2-Qh^2$ smaller.  Taking minima preserves this direction.  The two
constant-potential half-interval minimizers give
$\min E_{-1}=2\coth D$, so (A2) is at most $2\coth D$.  Multiplication by
$\tanh D/2$ is exactly D2.

For the original map $g=G\circ\tanh$, the chain rule gives

\[
 Sg=(SG\circ\tanh)\operatorname{sech}^4-2. \tag{A7}
\]

Hence $SG\ge0$ implies $Q=Sg/2\ge-1$.  No symmetry of $G$ or of the two
nodes is used.  The audit finds no gap in L17.  Its precise regularity
hypotheses are $G\in C^3$ and $G'>0$ on the closed node interval; the
degenerate endpoints in EL4 still require the stated continuity limit.

## 2. L21: SDP duality and parity averaging

For $r(T)<1$, consider

\[
 \min t:\quad P-I\succeq0,\quad tI-P\succeq0,\quad
 P-T^*PT\succeq0. \tag{A8}
\]

Slater is genuine.  The Lyapunov series
$P_0=\sum_{j\ge0}(T^*)^jT^j$ satisfies
$P_0-T^*P_0T=I$ and $P_0\succeq I$.  Taking $\lambda P_0$ with
$\lambda>1$, then $t>\lambda\|P_0\|$, makes all three constraints strict.

Pairing the constraints with $X,Y,Z\succeq0$ and minimizing the Lagrangian
over $P,t$ gives, with no sign ambiguity,

\[
 \max\operatorname{tr}X:\quad
 X-Y+Z-TZT^*=0,\qquad\operatorname{tr}Y=1. \tag{A9}
\]

For $D=Z-TZT^*$, (A9) says $X=Y-D\succeq0$.  The exact auxiliary problem is

\[
 \min\{\operatorname{tr}Y:Y\succeq0,\ Y\succeq D\}
 =\operatorname{tr}D_+. \tag{A10}
\]

The upper bound is attained by $Y=D_+$.  For the lower bound, compress any
admissible $Y$ to the positive spectral subspace of $D$ and use $Y\succeq0$.
Scaling the ray $Z\mapsto\lambda Z$ now gives

\[
 t_*(T)=\max\left(1,\sup_{Z\succeq0,Z\ne0}
 {\operatorname{tr}(Z-TZT^*)_-\over
  \operatorname{tr}(Z-TZT^*)_+}\right). \tag{A11}
\]

The denominator cannot vanish on a nonzero ray: $D\preceq0$ would imply
$Z\preceq TZT^*\preceq\cdots\to0$ by strict stability.

Finally suppose $JTJ=-T$.  Averaging a primal metric with $JPJ$ is valid.
For the dual, one must **not** average the nonlinear ratio directly.  Average
the full feasible triple $(X,Y,Z)$ in (A9).  Its traces and objective are
unchanged and the averaged $Z$ commutes with $J$.  If its objective is
$d=1+n-p>1$, where $p=\operatorname{tr}D_+$ and
$n=\operatorname{tr}D_-$, feasibility gives $p\le1$ and $n>p$.  The best
rescaling of that invariant ray has value

\[
 {n\over p}-(1+n-p)=(1-p)\left({n\over p}-1\right)\ge0. \tag{A12}
\]

Thus an invariant ratio witness exists at least at the original violating
level.  This justifies the parity restriction without a convexity claim
about the ratio.  The audit finds no gap in L21.

## 3. L29: clean regeneration of all exact certificates

The three certificate entry points cited by L29 were rerun from a clean
detached shell, with no saved symbolic expansion files supplied:

1. `slice_odd_centered_completion_certificate.py` rebuilt the 309,479-term
   upper and 565,425-term lower deep series and exited exact after 85 minutes;
2. `slice_odd_small_edge_certificate.py` regenerated both residual series and
   all rational theta checks, with worst correction ratio below `0.907895`;
3. `slice_odd_compact_certificate.py` closed its bridge, ridge, and direct
   ranges with no pending box (15,267, 15,047, and 17,594 splits).

The wrapper ended with `PASS 2026-07-21T15:10:52-07:00`.  This verifies
reproducibility of the computer-assisted part; it is independent of the
separate algebraic review already recorded in
`proof/slice_odd_block_reduction.md`.

## 4. Audit verdict

Both lemmas survive independent derivation.  L17's essential hidden-looking
step is secured by the explicit ground-state identity (A4), and L21's
nonlinear symmetry issue is secured only by averaging the linear dual triple
and then rescaling as in (A12).  These arguments should be retained whenever
the lemmas are reused; shorter claims that the ratio itself averages are not
valid substitutes.  L29's three regenerating checkers also pass cleanly, so
the audit debt identified for all three load-bearing results is discharged.
