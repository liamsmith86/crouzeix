# Hereditary transfer flags survive exact nonlinear Schur transport

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

> **Route update.**  L291 subsequently proves that the hereditary
> state form below is sufficient but stronger than necessary.  Every
> admissible positive-interior transport is a bounded graph gauge plus
> its literal endpoint residual, and flag zero is equivalent
> pointwise to hereditary factorization of that endpoint residual.
> The live A178 target is therefore the bounded endpoint factor, not
> state-ideal membership.

## 1. Result (L290, 2026-07-26)

Let \(S\) be a matrix-valued partial-isometry colligation with initial
and final defect frames

\[
V:\mathbb C^m\longrightarrow{\cal H},\qquad
W:\mathbb C^m\longrightarrow{\cal H},
\]

and projections \(E=VV^*\), \(F=WW^*\).  Write

\[
B_j=W^*(S^*)^jV,\qquad D_j=ES^jF.                \tag{1}
\]

Let \(H,\widetilde H\) be Hermitian state matrices or formal series
whose interior blocks relative to \(F\) are invertible, and let
\(J_H,J_{\widetilde H}\) be their endpoint graph columns from L289.
Assume that their difference has the finite hereditary form

\[
\boxed{
\Delta H=\widetilde H-H
=\sum_{j=1}^r\{A_jD_j+D_j^*A_j^*\}.}             \tag{2}
\]

The multipliers \(A_j\) may themselves be matrices or formal series.
Then the complete nonlinear Schur-endpoint change factors exactly as

\[
\boxed{
\begin{aligned}
{\cal U}(\widetilde H)-{\cal U}(H)
=\sum_{j=1}^r\{&
  (J_{\widetilde H}^*A_jV)B_j^*\\
 &+B_j(J_H^*A_jV)^*\}.
\end{aligned}}                                   \tag{3}
\]

Consequently, if \(U\) is any endpoint column satisfying

\[
B_1^*U=\cdots=B_r^*U=0,                          \tag{4}
\]

then

\[
\boxed{
U^*\{{\cal U}(\widetilde H)-{\cal U}(H)\}U=0.}    \tag{5}
\]

Thus the state-side hereditary delay channels \(D_j=ES^jF\), with
the generator exposed on the graph-facing side, are carried by the
full nonlinear Schur shorting into the endpoint flag ideal generated
by \(B_j\).  This remains true for arbitrary motion of both graph
columns and automatically includes every nonlinear Schur cross.

This is an exact sufficient invariant for A178.  It does **not** show
that the physical all-grade transport automatically has form (2);
constructing such representatives, or proving (3)'s endpoint
factorization directly for their mixed-graph pairing, remains open.

## 2. Graph and transfer identities

Regard \(W\) as the fixed endpoint inclusion and put
\(P=I-WW^*\).  The graph column is

\[
J_H=W-P(PHP)^{-1}PHW.                            \tag{6}
\]

It obeys

\[
W^*J_H=I,\qquad FJ_H=W,\qquad PHJ_H=0.           \tag{7}
\]

The same identities hold for \(J_{\widetilde H}\).  Therefore every
delay channel has the graph-independent endpoint value

\[
\boxed{
D_jJ_H=ES^jFJ_H
=ES^jW
=VB_j^*.}                                        \tag{8}
\]

Equation (8), rather than raw membership in a two-sided word ideal,
is the structural reason the transfer flag survives shorting.

## 3. Proof of the hereditary factorization

L289 gives the exact two-graph identity

\[
{\cal U}(\widetilde H)-{\cal U}(H)
=J_{\widetilde H}^*\Delta HJ_H.                  \tag{9}
\]

Insert one summand from (2).  Equation (8) gives

\[
\begin{aligned}
J_{\widetilde H}^*A_jD_jJ_H
  &=(J_{\widetilde H}^*A_jV)B_j^*,\\
J_{\widetilde H}^*D_j^*A_j^*J_H
  &=B_j(J_H^*A_jV)^*.
\end{aligned}                                    \tag{10}
\]

Summing (10) proves (3).  Multiplying it on the left by \(U^*\) and
on the right by \(U\), then using (4), proves (5).

The proof is coefficientwise valid for formal series when the two
interior constant coefficients are invertible.  In particular, the
multiplier coefficients may depend on all earlier preparation
grades; no assumption that the graph stays fixed is present.

The same formula gives the quantitative bound

\[
\boxed{
\|\Delta{\cal U}\|
\le\sum_{j=1}^r
\bigl(
 \|J_{\widetilde H}^*A_jV\|
 +\|J_H^*A_jV\|
\bigr)\|B_j\|.}                                  \tag{11}
\]

Hence L289's locally bounded analytic graph columns do not add a new
growth loss once a bounded hereditary state representative is known.

## 4. Scope: this is not the false raw ideal claim

L285 disproves automatic delay-ideal preservation for the raw
transport forcing.  There is no contradiction.  A general word
\(XD_jY\) need not expose \(D_j\) against a graph column, because the
arbitrary right multiplier \(Y\) destroys (8).  Two-sided membership
alone therefore need not imply flag annihilation after Schur
shorting.

The invariant in (2) is deliberately one-sided and hereditary:
every nonadjoint term ends in \(D_j\), and Hermitian symmetry supplies
the adjoint term.  Equivalently, the relevant object is the
mixed-graph endpoint pairing, not a raw state-polynomial quotient.

L286 and L288 are finite examples of the resulting endpoint
behavior.  L290 does not retroactively assert that each of their
stored metric witnesses literally has form (2).  It says that future
recurrence steps may be certified in either of two precise ways:

1. construct a bounded state representative of form (2); or
2. use L289 and factor the complete mixed-graph pairing directly as
   the right side of (3).

Any proposed induction must also preserve the physical direct Gram
and preceding Schur budgets.  Vanishing on a surviving flag is only
the transport part of that certificate.

## 5. Consequence for the finite-jet campaign

For a fixed repeated length \(L\), L289 requires only grades through
\(2L\).  L290 identifies the correct sufficient invariant at each of
those finitely many stages:

\[
[c^n]\Delta H
\quad\hbox{hereditarily generated by}\quad
D_1,\ldots,D_r
\]

implies that its complete nonlinear endpoint motion vanishes on
\(\bigcap_{j\le r}\ker B_j^*\).

The remaining A178 theorem is still substantial:

1. derive an arbitrary-grade finite recurrence whose relevant state
   changes have this hereditary form, or whose mixed-graph pairings
   factor as in (3);
2. bound the chosen multipliers and all retained even faces;
3. prove positive margins through the terminal grade \(L\).

L290 replaces the vague phrase “right-ideal preservation” by an
exact graph-compatible condition.  It is an elementary consequence
of L289 and (8), not a claim of a new general Schur-complement
theorem.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_graph_flag_ideal.py \
  --output \
  experiments/repeated_crabb_graph_flag_ideal_s70226.jsonl
```

The checker uses generic positive old states and hereditary
perturbations.  It verifies:

1. both graph equations and endpoint normalizations;
2. L289's exact two-graph identity;
3. the factorization (3);
4. compression to the surviving flag; and
5. the norm bound (11).

The deterministic suite contains one unstructured case and
rank-changing cases of defect multiplicities \(3,4,5\), with
nontrivial surviving flags in the latter three.  The tracked dataset
has SHA-256

```text
c0bd78fd90f4d1c7260cff29fb3a4c4d7bf8b10eb3497d6365bbc8061fef60a8
```
