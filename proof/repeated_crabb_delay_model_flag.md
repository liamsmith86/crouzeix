# Exact delay covariance of the Schur model and boundary metric

## 1. Result (L221, 2026-07-24)

Retain the finite pure partial isometry and transfer

\[
\begin{gathered}
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\\
B(z)=W^*(I-zS^*)^{-1}V=\sum_{n\ge1}B_nz^n.
\end{gathered}                                      \tag{1}
\]

Assume a complete delay of length \(r\ge1\):

\[
B_1=\cdots=B_r=0.                                  \tag{2}
\]

Let \(S_r,V_r,W_r\) be L209's deflated colligation after removing

\[
W,SW,\ldots,S^{r-1}W,
\]

and let

\[
\widetilde B(z)
=W_r^*(I-zS_r^*)^{-1}V_r
=\sum_{n\ge1}B_{n+r}z^n.                           \tag{3}
\]

Then the transfer, model kernel, L220 feature flag, right-defect
orbits, and L219 boundary metric all obey exact delay covariance.

First,

\[
\boxed{B(z)=z^r\widetilde B(z).}                   \tag{4}
\]

For the right model kernel

\[
K_B(z,w)=\frac{I-B(z)^*B(w)}{1-\overline zw},
\]

equation (4) gives the orthogonal splitting

\[
\boxed{
K_B(z,w)
=\sum_{j=0}^{r-1}(\overline zw)^jI
(\overline zw)^rK_{\widetilde B}(z,w).}            \tag{5}
\]

If \(\Phi_j\) are L220's ordered Schur features of \(B\), and
\(\widetilde\Phi_j\) are those of the tail \(\widetilde B\), then

\[
\boxed{
\Phi_j(z)=z^jI\quad(0\le j<r),\qquad
\Phi_{r+j}(z)=z^r\widetilde\Phi_j(z).}              \tag{6}
\]

Thus complete delay removal is an exact orthogonal tail operation in
the canonical Schur flag, not merely an associated-graded statement.

Let \(J:{\cal H}_r\to{\cal H}\) include the retained state space.
Every forward right-defect orbit has the exact orthogonal expansion

\[
\boxed{
(S^*)^nV
=J(S_r^*)^nV_r
 \sum_{j=0}^{r-1}S^jWB_{n+j}
\qquad(n\ge0).}                                    \tag{7}
\]

Finally, with \(q=c^2\), \(a_n=q^n/(1+q^n)\), and L219's metric

\[
P_{\rm bl}^{S}(c)
=I-\sum_{n\ge1}a_n(S^*)^nES^n
 \sum_{n\ge1}q^nS^nF(S^*)^n,
\]

its retained compression is

\[
\boxed{
\begin{aligned}
J^*P_{\rm bl}^{S}(c)J
={}&I-\sum_{n\ge1}a_n
  (S_r^*)^nV_rV_r^*S_r^n\\
&+q^r\sum_{n\ge0}q^n
  S_r^nW_rW_r^*(S_r^*)^n.
\end{aligned}}                                     \tag{8}
\]

Equations (5)--(8) identify every raw cross block and every metric
weight created by delay removal.  They do not yet prove that the
prepared physical elliptic endpoint is a model-kernel norm.  That
single identification remains the gate between L220's orthogonal
feature Gram and the desired \(-16\) physical Gram.

## 2. Transfer and model-kernel splitting

L210 gives

\[
\widetilde B_n=B_{n+r}
\]

under a complete delay.  Assumption (2) therefore proves (4)
coefficientwise.

Put \(t=\overline zw\).  Then

\[
\begin{aligned}
I-B(z)^*B(w)
&=I-t^r\widetilde B(z)^*\widetilde B(w)\\
&=(1-t^r)I
 t^r\{I-\widetilde B(z)^*\widetilde B(w)\}.
\end{aligned}
\]

Since

\[
\frac{1-t^r}{1-t}=\sum_{j=0}^{r-1}t^j,
\]

division by \(1-t\) proves (5).  A sum of positive kernels is an
orthogonal direct sum of their reproducing-kernel spaces, so

\[
{\cal K}_B
={\cal K}_{z^rI_m}\oplus z^r{\cal K}_{\widetilde B}. \tag{9}
\]

This proof uses no choice of state basis and no commutation of copy
matrices.

## 3. Compatibility with the ordered Schur flag

In L218's Schur recursion, a zero parameter gives

\[
F_j(z)=zF_{j+1}(z),\qquad A_j(z)=I.
\]

Equation (4) means that the first \(r\) steps have zero parameters.
L220 defines

\[
\Phi_j(z)=z^jA_j(z)A_{j-1}(z)\cdots A_0(z).
\]

For \(j<r\), every factor is the identity, proving the first formula
in (6).  For \(j=r+\ell\), split the ordered product after \(A_r\):

\[
\begin{aligned}
\Phi_{r+\ell}(z)
&=z^{r+\ell}A_{r+\ell}(z)\cdots A_r(z)\\
&=z^r\widetilde\Phi_\ell(z).
\end{aligned}
\]

This proves the second formula in (6), with the multiplication order
unchanged.  Substitution in L220's feature sum independently
regenerates (5).

## 4. State-orbit splitting

L209 proves that

\[
W,SW,\ldots,S^{r-1}W
\]

are mutually orthonormal and orthogonal to the retained space.
For

\[
g_n=(S^*)^nV,
\]

the coefficient of \(g_n\) on the \(j\)-th removed layer is

\[
(S^jW)^*g_n
=W^*(S^*)^{n+j}V
=B_{n+j}.                                         \tag{10}
\]

The retained component is

\[
J^*g_n
=(S_r^*)^nV_r,                                    \tag{11}
\]

because \(S^nJ=JS_r^n\) and \(V=JV_r\).  Combining the orthogonal
components (10)--(11) proves (7).

Formula (7) explains the early future-transfer images in the
one-image frame.  A raw forward orbit already contains all removed
components \(S^jWB_{n+j}\); Schur orthogonalization removes them as
one triangular block rather than grade by grade.

## 5. Boundary-metric compression

The right-orbit part of (8) follows from (11):

\[
J^*(S^*)^nES^nJ
=(S_r^*)^nV_rV_r^*S_r^n.                           \tag{12}
\]

For the left orbits, \(S^nW\) lies in the removed line when
\(0\le n<r\).  For \(n=r+\ell\),

\[
S^{r+\ell}W
=JS_r^\ell W_r.                                   \tag{13}
\]

Hence

\[
\begin{aligned}
J^*\left\{\sum_{n\ge1}q^nS^nF(S^*)^n\right\}J
&=\sum_{\ell\ge0}q^{r+\ell}
  S_r^\ell W_rW_r^*(S_r^*)^\ell,
\end{aligned}
\]

which is the second line of (8).  Norm convergence follows from
purity, as in L219.

## 6. Consequence for the physical gate

L216 already proves that every **linear** retained endpoint response
is unchanged by the same delay removal.  L221 adds the two nonlinear
structures that were not explicit there:

1. the model kernel and its Schur-whitened features split
   orthogonally by (5)--(6);
2. every lower-order boundary-metric block and cross row is fixed by
   (7)--(8).

Therefore no further coordinate or multiplication-order ambiguity
remains in delay deflation.  What is still unproved is the
load-bearing interface

\[
\text{physical prepared two-reflection quotient}
=-16\times
\text{Schur-feature Gram}.                          \tag{14}
\]

Once (14) is proved at one active feature layer, (5)--(8) and L216
transport it to every delayed layer.  L221 deliberately does not use
(14) in its proof.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_delay_model_flag.py \
  --output \
  experiments/repeated_crabb_delay_model_flag_s70224.jsonl
```

The deterministic records use fully delayed unstructured
colligations for the state-orbit and boundary-metric identities, and
independent noncommuting strict Schur tails for the kernel/feature
identities.  Delay lengths one through four and defect multiplicities
one and two are covered.  They verify (4)--(8), including the ordered
L220 feature shift.  The equations above, not the floating audit,
prove L221.  The tracked dataset SHA-256 is
`a47c7e965cb5d48a69bd10c943856e454b478fe0b0cf91118037ccf55dd59b10`.
