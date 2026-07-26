# The canonical transfer flux already has the gap-free energy bound

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

> **Scope correction (L295).**  The flux estimate below is exact,
> but the raw unretightened L227/L285 upper face is not bounded above
> by L212's lower-tight target: at the monomial apex their faces are
> `-12I` and `-16I`.  L297 shows the difference is exactly the
> favorable lower face and its retightening response; the effective
> lower-tight comparison remains open.  L296 supplies the estimate
> for arbitrary polarized L280 responses in the raw branch.

## 1. Result (L294, 2026-07-26)

Retain L280's bistochastic transfer channel

\[
\Phi(K)=\sum_{n\ge1}B_nKB_n^*,\qquad
\Phi^*(Y)=\sum_{n\ge1}B_n^*YB_n,
\]

and Markov Laplacian

\[
K_{\rm M}=I-\Phi\Phi^*.
\]

For each transfer grade \(k\), define L212's canonical endpoint flux

\[
\boxed{
G_k=\Phi(B_k^*B_k)-B_kB_k^*.}                    \tag{1}
\]

Then every Hermitian copy-space test \(Y\) satisfies the exact
gap-free estimate

\[
\boxed{
|\langle Y,G_k\rangle_{\rm HS}|
\le
\|B_k\|_F
\sqrt{\langle Y,K_{\rm M}Y\rangle_{\rm HS}}.}      \tag{2}
\]

More precisely, if \(A_Y=\Phi^*(Y)\), then

\[
\boxed{
\langle Y,G_k\rangle
=-\operatorname {Re}\operatorname {tr}
B_k^*(YB_k-B_kA_Y),}                              \tag{3}
\]

while L280 gives

\[
\boxed{
\langle Y,K_{\rm M}Y\rangle
=\sum_{n\ge1}\|YB_n-B_nA_Y\|_F^2.}                \tag{4}
\]

Thus (2) is ordinary Cauchy--Schwarz applied to one summand of the
exact physical observability defect.  No spectral gap, inverse,
rank projection, compactness argument, or pseudoinverse occurs.

The consequence for L282 is immediate.  If a first physical endpoint
face on a flag has the canonical form

\[
\boxed{
E_k^{\rm can}
=12B_kB_k^*-28\Phi(B_k^*B_k)
=-16B_kB_k^*-28G_k,}                              \tag{5}
\]

then every \(Y\succeq0\) obeys

\[
\boxed{
\langle Y,E_k^{\rm can}\rangle_+
\le
28\|B_k\|_F
\sqrt{\langle Y,K_{\rm M}Y\rangle}.}              \tag{6}
\]

Moreover L212's explicit polynomial column has endpoint response
\(28G_k\), so it changes (5) exactly to

\[
\boxed{-16B_kB_k^*\preceq0.}                      \tag{7}
\]

Equations (2) and (6) close the **quantitative Markov-energy
inequality** for the canonical channel-minus-Gram target.  They do
not prove that the fully transported physical endpoint has form (5),
nor do they turn a rank-jumping flag projection into one analytic
polynomial column.  The physical identification was the conditional
equation (6) of L212 and raw simultaneous superposition is false by
A179/L220.

For the L212 exact-axis branch, the live gate is:

> prove that its all-grade first physical flag face equals (5), or is
> bounded above by it modulo already favorable terms.

That statement uses a lower-tight normalization.  L295 disproves
comparing the raw L227/L285 upper face with (5), while L297 shows the
lower-retightened effective face is apex-aligned with it.  If one
retains the raw repair normalization, one must instead separate
L283's favorable
\(-12B_kB_k^*\) face and identify only the mixed transported
remainder as a flux controlled below.  On a moving partial flag, the
compressed form supplies L282's gap-free energy certificate without
inserting the flag projection into the metric.

## 2. Exact pairing identity

Put

\[
A_Y=\Phi^*(Y).
\]

Channel adjointness gives

\[
\begin{aligned}
\langle Y,\Phi(B_k^*B_k)\rangle
&=\langle A_Y,B_k^*B_k\rangle\\
&=\operatorname {tr}(B_k^*B_kA_Y).
\end{aligned}
\]

Also

\[
\langle Y,B_kB_k^*\rangle
=\operatorname {tr}(B_k^*YB_k).
\]

Subtracting proves

\[
\langle Y,G_k\rangle
=-\operatorname {tr}B_k^*(YB_k-B_kA_Y).
\]

The left side is real because \(Y\) and \(G_k\) are Hermitian, so
inserting the real part gives (3).

## 3. Cauchy--Schwarz and the physical face

Equation (3) gives

\[
|\langle Y,G_k\rangle|
\le\|B_k\|_F\|YB_k-B_kA_Y\|_F.                   \tag{8}
\]

The \(k\)-th squared norm on the right is one nonnegative summand of
L280's exact identity (4).  This proves (2).

For \(Y\succeq0\),

\[
\langle Y,-16B_kB_k^*\rangle\le0.
\]

Using (5), dropping this favorable term, and applying (2) gives (6).
Finally, L212 proves that its universal polynomial column \(C_k\)
satisfies

\[
{\cal M}_T(C_k)=28G_k.
\]

Adding that response to (5) proves (7).

For any isometry \(U\), put \(P=UU^*\),
\(\widehat Y=UYU^*\), and define the flagged flux

\[
G_{k,P}
=P\{\Phi(B_k^*PB_k)-B_kB_k^*\}P.                 \tag{8a}
\]

The same calculation gives

\[
\begin{aligned}
\langle \widehat Y,G_{k,P}\rangle
=-\operatorname {Re}\operatorname {tr}
(PB_k)^*P(\widehat YB_k-B_kA_{\widehat Y}).
\end{aligned}                                     \tag{8b}
\]

Therefore

\[
\boxed{
|\langle \widehat Y,G_{k,P}\rangle|
\le
\|PB_k\|_F
\sqrt{\langle\widehat Y,K_{\rm M}\widehat Y\rangle}.}          \tag{8c}
\]

In particular, the flagged canonical target

\[
\begin{aligned}
E_{k,P}^{\rm can}
&=12P B_kB_k^*P
  -28P\Phi(B_k^*PB_k)P\\
&=-16P B_kB_k^*P-28G_{k,P}
\end{aligned}                                                   \tag{8d}
\]

obeys

\[
\boxed{
\langle \widehat Y,E_{k,P}^{\rm can}\rangle_+
\le
28\|PB_k\|_F
\sqrt{\langle\widehat Y,K_{\rm M}\widehat Y\rangle}.}           \tag{8e}
\]

The inequality itself holds for every \(P\).  When \(P\) is the
grade-\(k\) transfer flag, so that \(PB_j=0\) for \(j<k\), L208's
pointwise flag column has exactly this compressed endpoint response.
Its literal use of \(P\) is not an analytic selection through rank
changes.  Inequality (8c), however, is precisely L282's
projection-free energy criterion: L281 realizes the resulting bounded
energy by the physical observability column.  Thus L294 does not
smuggle a discontinuous flag projection into the metric.

## 4. Weighted form

For nonnegative summable weights \(r_k\), put

\[
G(r)=\Phi\left(\sum_kr_kB_k^*B_k\right)
-\sum_kr_kB_kB_k^*.
\]

Summing (3) and applying Cauchy--Schwarz across the transfer grades
gives

\[
\boxed{
|\langle Y,G(r)\rangle|
\le
\left(\sum_kr_k\|B_k\|_F^2\right)^{1/2}
\left(\sum_kr_k\|YB_k-B_kA_Y\|_F^2\right)^{1/2}.} \tag{9}
\]

In particular, if \(0\le r_k\le r_{\max}\), the second factor is at
most

\[
\sqrt{r_{\max}}\,
\sqrt{\langle Y,K_{\rm M}Y\rangle}.
\]

This is the exact energy estimate for L212's weighted all-grade
candidate.  A179 says that inserting its raw columns simultaneously
does not produce the required physical endpoint series; (9) does not
revive that false claim.

The same proof gives the branch-neutral corollary

\[
E_{a,b}=-aB_kB_k^*-bG_k,\qquad a\ge0,\ b\in\mathbb R,
\]

\[
\boxed{
\langle Y,E_{a,b}\rangle_+
\le |b|\,\|B_k\|_F
\sqrt{\langle Y,K_{\rm M}Y\rangle}.}              \tag{10}
\]

L295 shows that this general form, rather than comparison of the raw
repair upper face with (5), is the correct possible interface to
L285/L289 when its favorable lower budget is retained.  L297 removes
the apex obstruction to the effective normalization, and L298
subsequently closes the desired complete-delay lower-tight face by
using one seventh of this preimage together with a negative orbit
Gram.  The nonlinear moving-series transport remains open.

## 5. Relation to L293

L293 expands the right side of (2) in L218's matrix-Schur
coordinates.  If

\[
\Gamma_j=\varepsilon\Delta_j+O(\varepsilon^2),
\]

then

\[
\langle Y,K_{\rm M}Y\rangle
=2\varepsilon^2\sum_j\|[Y,\Delta_j]\|_F^2
+O(\varepsilon^3).
\]

Thus L294 supplies the exact nonlinear numerator pairing whose first
Schur jet is controlled by L293's commutator square.  The two routes
are now literally the same Cauchy--Schwarz estimate.  L295 separates
the raw repair upper face from L212's lower-tight target, and L297
identifies their apex gap as lower retightening.  The remaining
physical work is either to prove the effective lower-tight base (5),
or retain the raw favorable lower budget and derive an
\(E_{a,b}\)-type formula for its mixed remainder.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_canonical_flux_inequality.py \
  --output \
  experiments/repeated_crabb_canonical_flux_inequality_s70226.jsonl
```

The exact SymPy checker uses both noncommuting matrix-unit
depolarizing channels and noncommuting random-unitary-type
bistochastic channels in multiplicities two through four.  It
verifies:

1. both unital identities;
2. the exact pairing identity (3);
3. L280's Dirichlet identity (4);
4. nonnegativity of the Cauchy--Schwarz slack in (2);
5. the flagged pairing, Cauchy--Schwarz, target, and correction
   identities (8a)--(8e);
6. the canonical correction (7); and
7. (6) and (8e) in cases where the uncorrected canonical candidate
   has a genuinely positive pairing.

Every residual and inequality is exact; no floating tolerance is
used.
The tracked dataset has SHA-256

```text
1ea9886ad59dc1c0facfa93f43cc17b00a67db57c67a512cceca9621b5fa8110
```
