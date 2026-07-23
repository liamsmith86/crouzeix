# Tied-winner diagonal quotient and graph sign (2026-07-22)

## 1. Diagonal support-kernel quotient

For one diagonal copy block `D`, let

\[
 \ell_D(q)=r(q)^*\{q^{-1}D+qD^*\}r(q),
 \qquad r(q)=\tfrac12(q^{-1},\sqrt2,q)^T.              \tag{1}
\]

The condition `ell_D(q)=0` for every unit `q` is the diagonal analogue of
L74's cross equation.  Its real Fourier matrix has rank seven on the 18 real
entries of `D`, so its kernel has dimension 11.  Infinitesimal within-copy
unitary commutators `[C_3,K]`, `K*=-K`, have rank eight inside that kernel.

Use the L74 generators `X_j,Y_j`.  The orthogonal three-real-dimensional
quotient is

\[
 D_0(\gamma)=\bar\gamma X_0+\gamma Y_0,
 \qquad D_1(s)=sX_1\quad(\gamma\in\mathbb C,\ s\in\mathbb R),              \tag{2}
\]

because `X_1=Y_1`.  The real and imaginary parts of `D_0` and the real
direction `D_1` satisfy (1), are orthogonal to every unitary commutator, and
together with the rank-eight orbit have rank 11.  Thus (2) is the complete
diagonal support-kernel quotient, not an ansatz.

## 2. Add the diagonal loops to the complete winner graph

Keep every off-diagonal L74 edge from L83 and put

\[
 E_{aa}=D_0(\gamma_a)+D_1(s_a)                         \tag{3}
\]

on each diagonal copy block.  Define the Hermitian copy-space matrix `H_1`
by

\[
 (H_1)_{ab}=\alpha_1^{ab},\quad
 (H_1)_{ba}=\bar\alpha_1^{ab},\quad
 (H_1)_{aa}=s_a.                                      \tag{4}
\]

For an edge, retain L76's metric block `Z(alpha)`.  For the loop (3), use the
Hermitian diagonal metric block

\[
 Z(\gamma_a,s_a/2,0)+Z(\gamma_a,s_a/2,0)^*.           \tag{5}
\]

All lower, upper, and Stein first active compressions still vanish.

Let `Q(q)`, `Qbar`, and `k_0=mean lambda_max Q(q)` denote the complete second
effective support data, now including every edge--edge, edge--loop, and
loop--loop path.  Exact propagation of the lower and contraction Schur
penalties through both Crabb levels gives the unchanged identity

\[
 \boxed{{\cal E}=16(\overline Q-k_0I)-8H_1^2.}         \tag{6}
\]

The first conformal Fourier mode again cancels.  Since
`Q(q)<=lambda_max(Q(q))I` pointwise and `H_1^2>=0`, (6) proves

\[
 \boxed{{\cal E}\preceq0}.                            \tag{7}
\]

Thus all canonical diagonal **relative** motions can be added to the entire
tied-winner graph without creating a positive second-order coefficient.  The
real loop `s_a` simply fills the diagonal of the same generator-one matrix
which already controlled the off-diagonal graph.

## 3. Equality and exact remaining scope

Equality of this explicit endpoint requires

\[
 H_1x=0,\qquad Q(q)x=\lambda_{\max}(Q(q))x
 \quad\hbox{for almost every }q.                      \tag{8}
\]

Generic cyclic and long-path examples are strict.  Star/pair configurations
can satisfy (8); L77 and L80 close their pure fixed rays, but a classification
of every graph equality and its weighted mergers remains.

There is also a distinct common-motion term.  If the first support compression
is a scalar `h(q)I` rather than zero, subtracting one reference diagonal block
leaves exactly the quotient (2), but the common term `I_m tensor D` remains.
Equation (7) does **not** yet control mixed terms between that shared single-
block motion and the relative graph/loop data.  L85 subsequently retains and
closes those mixed terms; they cannot be silently discarded as gauge.

## 4. Regeneration

Run

```bash
.venv/bin/python -u experiments/repeated_p3_common_maximizer.py
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
```

The first checker proves the ranks `7/8/3/11` and exact orthogonality for the
diagonal quotient.  The second uses a fully symbolic three-copy triangle,
three independent complex `D_0` loops, three real `D_1` loops, and all nine
complex edge coefficients.  It reconstructs the actual Laurent mean of
`Q(q)` and verifies (6) entry by entry.  Every quadratic term is a path on at
most three copy vertices, so the triangle audit contains all general-
multiplicity path types.
