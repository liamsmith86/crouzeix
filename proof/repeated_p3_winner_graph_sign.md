# Tied-winner graph sign for repeated `C3` (2026-07-22)

## 1. Complete winner graph

Let `A_0=I_m tensor C_3` and suppose the entire copy space is tied at first
support order after the common boundary motion and orbit directions have been
removed.  For every unordered copy pair `a<b`, place one L74 quotient pair in
the blocks `(E_ba,E_ab)`:

\[
 E_{ba}=\sum_{j=0}^2\bar\alpha_j^{ab}X_j,
 \qquad
 E_{ab}=\sum_{j=0}^2\alpha_j^{ab}Y_j.                 \tag{1}
\]

All diagonal copy blocks are zero in this canonical quotient.  Equation (3)
of L74 implies that the complete first top-support compression vanishes, not
merely one selected row and column.

Define the Hermitian copy-space matrix `H_1` by

\[
 (H_1)_{ab}=\alpha_1^{ab},\quad
 (H_1)_{ba}=\bar\alpha_1^{ab},\quad (H_1)_{aa}=0.      \tag{2}
\]

Let `Q(q)` be the second effective support matrix on the `m`-dimensional top
copy space,

\[
 Q(q)=P_{\rm top}V(q)R(q)V(q)P_{\rm top},qquad
 \kappa(q)=\lambda_{\max}Q(q),                        \tag{3}
\]

and put

\[
 \overline Q=\frac1{2\pi}\int Q(e^{i\theta})\,d\theta,
 \qquad k_0=\frac1{2\pi}\int\kappa(e^{i\theta})\,d\theta. \tag{4}
\]

## 2. Pairwise first metric

For each edge `a<b`, put the L76 metric block

\[
 Z(\alpha^{ab})=\begin{bmatrix}
0&-3\sqrt2\alpha_0^{ab}/8&0\\
\alpha_2^{ab}/\sqrt2&-2\sqrt2\alpha_1^{ab}
                         &3\sqrt2\alpha_0^{ab}/4\\
0&-\sqrt2\alpha_2^{ab}&0
\end{bmatrix}                                         \tag{5}
\]

in copy block `(a,b)` and its adjoint in `(b,a)`.  All first active lower,
upper, and Stein compressions vanish.  Insert the complete lower and Stein
Schur penalties and propagate the second metric through both positive Crabb
levels.

The second Schwarz map can have a nonzero first Fourier mode, so its evaluation
at `A_0` contains an `A_0^2` term.  As in L78, retain an arbitrary complex
coefficient for that term during the calculation rather than assuming parity.

## 3. Endpoint identity and sign

Exact edge-path collection gives

\[
 \boxed{
 {\cal E}=16(\overline Q-k_0I_m)-8H_1^2.}             \tag{6}
\]

The first conformal Fourier mode cancels identically.  Since

\[
 Q(q)\preceq\kappa(q)I_m
\]

pointwise, averaging gives `overline Q <= k_0 I_m`.  Also `H_1^2` is PSD.
Therefore

\[
 {\cal E}\preceq0,
 \qquad
 \limsup_{\epsilon\downarrow0}
 \frac{t_*(T_\epsilon)-4}{\epsilon^2}
 \le\lambda_{\max}({\cal E})\le0.                    \tag{7}
\]

This proves nonpositive second-order change for the full tied-winner quotient,
including arbitrary cycles and all winner--winner edges.  L78 is the special
case in which the graph is a star.

## 4. Orientation warning

It is tempting to define analogous Hermitian matrices `H_0,H_2` and write

\[
 \overline Q\stackrel{?}{=}
 \frac5{128}H_0^2+\frac14H_1^2+\frac5{72}H_2^2.
\]

That identity is **false on graphs with cycles**.  For a three-copy triangle,
oriented two-edge paths contribute mixed generator-0/generator-2 terms.  They
cancel on a star's diagonal/Gram structure but not around the complete graph.
The proof of (6) uses the actual mean (4), so it does not rely on this false
shortcut.

This orientation dependence also means that arbitrary copy-unitary changes
do not preserve the zero-diagonal pairwise gauge term by term; diagonal orbit
pieces must be removed again after changing copy basis.

## 5. Equality and remaining scope

Equality in this particular endpoint requires simultaneously

\[
 H_1x=0,
 \qquad Q(q)x=\kappa(q)x\quad\text{for almost every }q. \tag{8}
\]

As L79 demonstrated for a star, equality of this explicit certificate need
not be equality of the optimized second-order metric.  Random complete
three-copy triangles tested numerically have strict negative quadratic
coefficients, but a full graph-equality optimization has not yet been proved.

The theorem also assumes zero diagonal winner blocks.  L84 subsequently
classifies every relative diagonal support-kernel residual and proves that
adding all of them preserves (6).  Coupling the resulting graph/loop data to a
common single-block motion, and obtaining a uniform estimate as losing sectors
merge into the winner graph, remain.

## 6. Regeneration and dimension independence

Run

```bash
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
```

The checker uses a symbolic three-copy triangle with all nine complex edge
coefficients independent.  It reconstructs `Q(q)`, its exact Laurent mean,
the full pairwise first metric, every second Schur penalty, and proves (6)
entry by entry while retaining an arbitrary first conformal mode.

Every quadratic endpoint entry in general multiplicity is a sum over
two-edge paths and hence involves at most three copy vertices.  The symbolic
triangle contains every such oriented path type; relabeling and summing proves
the dimension-independent identity.
