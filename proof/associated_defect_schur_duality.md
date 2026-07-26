# Initial and final defect Schur faces are associated-grade dual

> **Closure note (2026-07-25).**  L283 applies this operator
> conjugacy to L278's intact face and proves L228's full delayed
> matrix anticommutator, not merely its trace.

## 1. Result (L246, 2026-07-24)

Let

\[
B(\varepsilon)=S+O(\varepsilon)
\]

be a square analytic matrix series, where \(S\) is a partial
isometry with equal defect projections

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad Q=I-E,\qquad P=I-F.
\]

Define its initial and final defect operators

\[
D_R=I-B^*B,\qquad D_L=I-BB^*,
\]

and take their Schur complements away from the nonzero defect
pivots:

\[
\begin{aligned}
K_R&=QD_RQ-QD_RE(ED_RE)^{-1}ED_RQ,\\
K_L&=PD_LP-PD_LF(FD_LF)^{-1}FD_LP.
\end{aligned}                                      \tag{1}
\]

If

\[
K_R=O(\varepsilon^d),
\]

then

\[
\boxed{
K_L=O(\varepsilon^d),\qquad
[\varepsilon^d]K_L
=S[\varepsilon^d]K_RS^*.}                         \tag{2}
\]

The converse holds with \(S^*\).  In particular the first nonzero
initial and final Schur faces have the same trace, rank, inertia, and
nonzero spectrum.

The statement survives an analytic metric similarity.  Let
\(R(\varepsilon)=I+O(\varepsilon)>0\), put

\[
\begin{aligned}
H_R&=R-A^*RA,\\
H_L&=R^{-1}-AR^{-1}A^*,
\end{aligned}                                      \tag{3}
\]

and suppose \(A(0)=S\).  If the fixed-\(E\) Schur residual of \(H_R\)
first appears in degree \(d\), then the fixed-\(F\) residual of
\(H_L\) first appears in the same degree and their leading faces
again obey (2).

Applied to L225's boundary metric, this converts the unresolved
right-defect face into an equivalent **left-defect dual-Stein face**.
It does not prove L228 or the candidate delayed trace law by itself.
Its value is that the left endpoint is the remote clean-chain
endpoint in L245, so the remaining state-lift cancellation can be
audited from either end without changing the first face.

For L228's proposed face the conjugate has a particularly simple
form.  With \(E_1=S^*ES\) and
\(F_j=S^jF(S^*)^j\),

\[
\begin{aligned}
S(E_1F_{k-1})S^*
&=SS^*ES^kF(S^*)^k\\
&=(I-F)EF_k=EF_k.
\end{aligned}
\]

Taking adjoints gives

\[
\boxed{
S(E_1F_{k-1}+F_{k-1}E_1)S^*
=EF_k+F_kE.}                                      \tag{4}
\]

Consequently L228 is equivalent, at its first nonzero grade, to the
left-final-defect target

\[
\boxed{[c^{2k}]K_L=EF_k+F_kE.}                    \tag{5}
\]

This is the active Hardy cell \((0,k)\), rather than \((1,k-1)\).
It aligns the unknown face with the zero-th coisometric boundary row
in L244.

## 2. Exact graph proof

Write the right graph injection

\[
G_R
=Q-E(ED_RE)^{-1}ED_RQ.                            \tag{6}
\]

It obeys

\[
ED_RG_R=0,\qquad
QD_RG_R=K_R,
\]

and hence, as a map from \(Q\mathcal H\) into \(\mathcal H\),

\[
\boxed{D_RG_R=K_R.}                               \tag{7}
\]

Likewise

\[
G_L
=P-F(FD_LF)^{-1}FD_LP
\]

satisfies

\[
\boxed{D_LG_L=K_L.}                               \tag{8}
\]

The elementary defect intertwining identity is

\[
\boxed{D_LB=BD_R.}                                \tag{9}
\]

Multiply (9) on the right by \(G_R\) and on the left by \(G_L^*\).
Using (7)--(8) and self-adjointness gives

\[
K_L\,U=V\,K_R,                                    \tag{10}
\]

where

\[
U=PB G_R:Q\mathcal H\to P\mathcal H,\qquad
V=G_L^*BQ:Q\mathcal H\to P\mathcal H.
\]

At \(\varepsilon=0\),

\[
G_R=Q,\qquad G_L=P,\qquad U=V=S|_{Q\mathcal H}.
\]

The last map is unitary from \(Q\mathcal H\) onto \(P\mathcal H\).
Thus \(U\) is invertible as a formal series, and (10) gives the exact
quotient identity

\[
K_L=V K_R U^{-1}.                                 \tag{11}
\]

If \(K_R=O(\varepsilon^d)\), taking its first coefficient in (11)
proves

\[
[\varepsilon^d]K_L
=S[\varepsilon^d]K_RS^*.
\]

Applying the same argument to \(B^*\) proves the converse.

## 3. Metric congruence

Set

\[
\widetilde B=R^{1/2}AR^{-1/2}.
\]

Then

\[
\begin{aligned}
H_R&=R^{1/2}(I-\widetilde B^*\widetilde B)R^{1/2},\\
H_L&=R^{-1/2}(I-\widetilde B\widetilde B^*)R^{-1/2}.
\end{aligned}                                      \tag{12}
\]

Congruence transports the two pivot spaces to analytic graph
perturbations of \(E\mathcal H\) and \(F\mathcal H\).  A Schur
complement is the quadratic form induced on the corresponding
quotient.  If that quotient form vanishes through degree \(d-1\),
changing its analytic graph coordinates by \(I+O(\varepsilon)\)
cannot alter its degree-\(d\) coefficient.  Since \(R(0)=I\),
equation (2) for \(\widetilde B\) therefore gives the same equation
for the fixed-frame Schur residuals of (3).

This also shows why lower-order metric square-root coefficients do
not have to be expanded to use the lemma.

## 4. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/associated_defect_schur_duality.py \
  --output \
  experiments/associated_defect_schur_duality_s70224.jsonl
```

The checker constructs unitary left/right analytic paths through
generic partial isometries, adds a first non-partial-isometric
perturbation independently in degrees one through six, and applies
an unrelated positive analytic metric congruence.  It verifies all
earlier vanishings, the initial/final face conjugacy, and trace
equality.  The algebra above, not the floating audit, proves the
result.  The tracked dataset SHA-256 is
`80117de3080c1c1986827d2438c58232807227ac1ad0b385680b211fcf977283`.
