# The Gau--Wu boundary loss is a first-jet quantity

> **Status and scope.**  The cancellation and reduction below are
> exact at every finite nondegenerate Gau--Wu model.  They remove the
> second-support curvature from L340's loss, but they do not prove
> that loss is positive or contractive.  The sign of the boundary
> shape form and the Crouzeix conjecture remain open.

## 1. The inverse-Riemann second jet

Let \(A=XSX^{-1}\) be a Gau--Wu disk model.  For a raw physical
direction \(E\), let \(s_E\) and \(t_E\) be its first and second
support coefficients.  Put

\[
 K_s(w)=wH_{s_E}(w),\qquad K_t(w)=wH_{t_E}(w),
\]

and let \(a_E=H_{\mathbb T}s_E-\partial_\theta s_E\) be the
first angular reparametrization of the boundary correspondence.
Write \(K_{a^2/2}=wH_{a_E^2/2}(w)\).

The inverse-Riemann normalized jet is

\[
 T_\varepsilon=A+\varepsilon G+\varepsilon^2H+o(\varepsilon^2),
\]

where

\[
\begin{aligned}
 G&=E-K_s(A),\\
 H&=-DK_s(A)[E]+K_s'(A)K_s(A)
     +K_{a^2/2}(A)-K_t(A).
\end{aligned}                                      \tag{1}
\]

Because \(E=G+K_s(A)\), and because \(K_s(A)\) commutes with \(A\),

\[
 DK_s(A)[K_s(A)]=K_s'(A)K_s(A).
\]

Thus the apparent second-order expression collapses exactly to

\[
\boxed{
 H=-DK_s(A)[G]+K_{a^2/2}(A)-K_t(A).}             \tag{2}
\]

The part preceding \(-K_t(A)\) depends only on the first jet
\((E,s_E,G)\).

## 2. Paired-frame expansion

Move to shift coordinates and write

\[
 C=X^{-1}GX,\quad
 D^{(1)}=-DK_s(S)[C]+K_{a^2/2}(S),\quad
 D=D^{(1)}-K_t(S).
\]

Retain L341's frames \(R_\varepsilon,L_\varepsilon\), with

\[
 R_0=I-pp^*,\qquad L_0=I-qq^*,
\]

and first responses \(R_1,L_1\).  Define

\[
\begin{aligned}
 \Delta_R^{(1)}&=-(S^*C+C^*S),\\
 \Delta_L^{(1)}&=-(SC^*+CS^*),\\
 \Delta_R^{(2)}(Z)&=-(S^*Z+Z^*S+C^*C),\\
 \Delta_L^{(2)}(Z)&=-(SZ^*+ZS^*+CC^*).
\end{aligned}                                      \tag{3}
\]

Since the moving frames retain rank \(n-1\), L341 gives

\[
 p^*R_2p=\|R_1p\|^2,\qquad
 q^*L_2q=\|L_1q\|^2.
\]

Consequently the optimized two-sided endpoint form has the exact
decomposition

\[
\boxed{
 {\cal H}_{\rm phys}(E)
 ={\cal B}(E)+\Xi(E),}                             \tag{4}
\]

where

\[
\begin{aligned}
 {\cal B}(E)
 &=\operatorname {tr}\{\Delta_R^{(2)}(D)R_0\}
  +\operatorname {tr}\{\Delta_L^{(2)}(D)L_0\},\\
 \Xi(E)
 &=\operatorname {tr}\{\Delta_R^{(1)}R_1\}
   +\|R_1p\|^2\\
 &\quad+\operatorname {tr}\{\Delta_L^{(1)}L_1\}
   +\|L_1q\|^2.                                   \tag{5}
\end{aligned}
\]

Here \(R_1,L_1\) use L338's exactly optimizing zero tangent.  All
traces in (5) are real.

## 3. Exact cancellation of the second support

The \(-K_t(S)\) part of \({\cal B}\) contributes precisely

\[
\boxed{
 {\cal B}_t(E)=2{\cal P}_\phi(E),}                 \tag{6}
\]

where \({\cal P}_\phi\) is L335's positive second-support Gram.
One proof is to combine L335 with L336: \(-K_t(A)\) contributes
\(-{\cal P}_\phi\) to the sharp norm Hessian, while L336 multiplies
the two-sided endpoint form by \(-1/2\), and the first-image penalty
has no second-support term.

Put

\[
 {\cal B}_1(E)
 =\operatorname {tr}\{\Delta_R^{(2)}(D^{(1)})R_0\}
  +\operatorname {tr}\{\Delta_L^{(2)}(D^{(1)})L_0\}.
\]

Equations (4)--(6) give

\[
 {\cal H}_{\rm phys}
 =2{\cal P}_\phi+{\cal B}_1+\Xi.                  \tag{7}
\]

Therefore L340's loss satisfies the first-jet identity

\[
\boxed{
 {\cal G}_\phi(E)
 =2\|E\|_F^2-{\cal B}_1(E)-\Xi(E).}               \tag{8}
\]

Every occurrence of \(t_E\), \(K_t\), and the moving support reduced
resolvent has cancelled.  Using \(Sp=0\) and \(S^*q=0\), the base
term can also be written

\[
\boxed{
 {\cal B}_1
 =4\operatorname {Re}\operatorname {tr}
   \{S^*(DK_s(S)[C]-K_{a^2/2}(S))\}
  -2\|C\|_F^2+\|Cp\|^2+\|C^*q\|^2.}              \tag{9}
\]

Thus (8) is determined entirely by the first conformal support
velocity, the first normalized operator direction, and the paired
first frame responses.

## 4. Revised proof target

L340 numerically observes

\[
 0\preceq{\cal G}_\phi
 \preceq2I+2{\cal P}_\phi.                        \tag{10}
\]

The upper inequality is exactly L338's remaining endpoint-flux
positivity.  Equation (8) separates its two jobs cleanly:

1. factor the first-jet form \(2I-{\cal B}_1-\Xi\) while keeping the
   two endpoint frames coupled; and
2. prove that this first-order output is dominated by the raw input
   energy plus L339's support-port energy.

Do not attempt to factor the two frame terms separately, and do not
carry \(K_t\) into the loss factorization: (6) proves that doing so
duplicates a term which has already cancelled.

The still more aggressive shortcut \({\cal G}_\phi=-\Xi\) is also
false.  It discards the first conformal-curvature reserve
\(2I-{\cal B}_1\).  In the standard models, \(-\Xi\) has minimum
eigenvalue between \(-2.95\) and \(-4\), whereas
\(2I-{\cal B}_1\) is positive definite.  The cancellation between
these two first-jet pieces is load-bearing.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u experiments/gau_wu_first_jet_loss.py \
  --output experiments/gau_wu_first_jet_loss_s70224.jsonl
```

The checker independently reconstructs the support jets, verifies
(2), solves both L341 first-frame Stein equations, polarizes
(4), (6), and (8), and reproduces L340's rank law.  The numerical
rank and positivity are retained only as evidence; the jet
reduction and cancellation are the exact claims.  Five models in
dimensions \(4,\ldots,8\) give second-jet and support-cancellation
errors below \(8.28\cdot10^{-16}\), and endpoint/loss reconstruction
errors below \(7.70\cdot10^{-12}\).  The dataset SHA-256 is
`15507d65397ecf1e0845c5d7831eb5ebabfd802b16b45e8234e592802546a78e`.
