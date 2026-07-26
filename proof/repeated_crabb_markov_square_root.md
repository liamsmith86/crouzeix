# The Markov response has an exact observability square root

## 1. Result (L281, 2026-07-25)

Retain L280's pure partial-isometry colligation, bistochastic transfer
channel

\[
\Phi(K)=\sum_{k\geq1}B_kKB_k^*,
\]

and endpoint response \({\cal M}_T\).  For \(H=H^*\), let

\[
\begin{aligned}
{\cal H}_H-S{\cal H}_HS^*&=WHW^*,\\
A_H&=V^*{\cal H}_HV=\Phi^*(H),\\
R_H&=(I-VV^*){\cal H}_HV.
\end{aligned}                                      \tag{1}
\]

Then L280's Markov coboundary has the much shorter exact state
preimage

\[
\boxed{
{\cal M}_T\left(2P_{\rm met}^{1/2}R_H\right)
=8(I-\Phi\Phi^*)H.}                               \tag{2}
\]

More importantly, its state-column norm is exactly the Markov
Dirichlet energy:

\[
\boxed{
\left\|2R_H\right\|_F^2
=2\left\langle H,(I-\Phi\Phi^*)H\right\rangle_{\rm HS}.}       \tag{3}
\]

The polarized form is

\[
\boxed{
\left\langle Y,(I-\Phi\Phi^*)H\right\rangle_{\rm HS}
=2\operatorname {Re}\operatorname {tr}(R_Y^*R_H).}             \tag{4}
\]

Thus a closing Markov spectral gap does not itself make the physical
state preimage singular.  The remaining selection target is exactly:

> solve L280's Poisson inequality with bounded Dirichlet energy.

There is no longer a separate state-synthesis or analytic-preimage
estimate.  The auxiliary copy variable \(H\) may diverge; the metric
chart sees only \(2R_H\), whose squared norm is (3).

L281 still does not prove that the physical faces admit a uniformly
bounded-energy Poisson solution.  L282 subsequently replaces that
selection debt by an equivalent quantitative off-commutant flux
inequality.

## 2. The defect column and transfer adjoint

The Stein series gives

\[
{\cal H}_H
=\sum_{n\geq0}S^nWHW^*(S^*)^n.
\]

Consequently

\[
A_H
=\sum_{n\geq1}B_n^*HB_n
=\Phi^*(H),                                       \tag{5}
\]

where the \(n=0\) cell vanishes because \(W^*V=0\).
As in L206,

\[
{\cal H}_HW=WH,\qquad
{\cal H}_HV=VA_H+R_H.                              \tag{6}
\]

The second identity is the orthogonal decomposition into the right
defect and its complement.

## 3. Exact energy balance

Multiply the Stein equation in (1) on the right by \(S\).
Since \(W^*S=0\),

\[
{\cal H}_HS=S{\cal H}_H(I-VV^*).
\]

Therefore

\[
S^*{\cal H}_HS
=(I-VV^*){\cal H}_H-R_HV^*.                       \tag{7}
\]

Multiply the Stein equation by \({\cal H}_H\), take the trace, and
use \({\cal H}_HW=WH\):

\[
\|H\|_F^2
=\operatorname {tr}{\cal H}_H^2
-\operatorname {tr}(S^*{\cal H}_HS{\cal H}_H).    \tag{8}
\]

Insert (7).  From (6),

\[
\begin{aligned}
\operatorname {tr}\{VV^*{\cal H}_H^2\}
&=\|A_H\|_F^2+\|R_H\|_F^2,\\
\operatorname {tr}(R_HV^*{\cal H}_H)
&=\|R_H\|_F^2.
\end{aligned}
\]

Hence (8) becomes

\[
\boxed{
\|H\|_F^2-\|A_H\|_F^2=2\|R_H\|_F^2.}             \tag{9}
\]

On the other hand, bistochasticity and (5) give

\[
\begin{aligned}
\left\langle H,(I-\Phi\Phi^*)H\right\rangle
&=\|H\|_F^2-\|\Phi^*(H)\|_F^2\\
&=\|H\|_F^2-\|A_H\|_F^2.
\end{aligned}
\]

Together with (9), this proves (3).  Real polarization proves (4).

## 4. Endpoint response

L204/L208's physical normalization says that every balanced
perpendicular column \(C\) satisfies

\[
\left\langle
Y,{\cal M}_T(P_{\rm met}^{1/2}C)
\right\rangle
=8\operatorname {Re}\operatorname {tr}(C^*R_Y).   \tag{10}
\]

Set \(C=2R_H\).  Equations (4) and (10) give

\[
\begin{aligned}
\left\langle
Y,{\cal M}_T(2P_{\rm met}^{1/2}R_H)
\right\rangle
&=16\operatorname {Re}\operatorname {tr}(R_H^*R_Y)\\
&=8\left\langle Y,(I-\Phi\Phi^*)H\right\rangle .
\end{aligned}
\]

This holds for every Hermitian \(Y\), proving (2).

At a repeated monomial apex, \(I-\Phi\Phi^*=0\), so (3) forces
\(R_H=0\) for every \(H\).  The preimage therefore vanishes without
choosing a kernel frame.

## 5. Consequence for bounded selection

For a proper flag \(U\), seek

\[
U^*\{E+8(I-\Phi\Phi^*)H\}U\prec0.
\]

L222+L279 prove pointwise feasibility.  L281 says the corresponding
free metric row has balanced norm

\[
\sqrt{
2\langle H,(I-\Phi\Phi^*)H\rangle}.
\]

Thus it is sufficient—and now necessary for this canonical
square-root choice—to prove a flagwise energy estimate

\[
\langle H,(I-\Phi\Phi^*)H\rangle
\leq C\,{\cal E}_{\rm active},
\]

where \({\cal E}_{\rm active}\) is the first transfer energy on that
flag.  L279 identifies the separator margin with exactly that energy.
The next step should derive the inequality from L197/L220's ordered
feature flag or falsify it numerically on scaled rank chains.

## 6. Regeneration

The checker is shared with L280:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_markov_response.py \
  --output \
  experiments/repeated_crabb_markov_response_s70225.jsonl
```

In addition to L280's audits, it verifies (2), (3), and the collapse
of the polarized all-grade column to \(2R_H\).  The tracked dataset
has SHA-256

```text
18a8c075ceb72d50871d7ff2d1ae443e8e5732dc2ce0e3f6957dc41901dc7fda
```
