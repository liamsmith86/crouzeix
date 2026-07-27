# The reduced Gau--Wu physical loss has boundary-size rank

> **Status.**  Every statement in this note about positivity and rank
> is numerical evidence, not a theorem.  The exact ingredients
> \(H_{\rm phys}\) and \({\cal P}_\phi\) are supplied by L338 and
> L339.  The purpose of the experiment is to identify the smallest
> plausible factorization target before opening its proof.

## 1. The reduced physical loss (L340/A287, 2026-07-26)

After L338 completes every zero-velocity square, write

\[
 H_{\rm phys}(C)
 =Q(C,0)-2\|\kappa_C\|^2.                         \tag{1}
\]

This is the exact physical Schur remainder whose nonnegativity would
close L336.  Let \({\cal P}_\phi\) be L335/L339's exact support-port
Gram and define the loss form

\[
\boxed{
 {\cal G}_\phi
 =2I+2{\cal P}_\phi-H_{\rm phys}.}                \tag{2}
\]

In twelve complete generic models, two in each dimension
\(4,\ldots,9\), the following exact-looking pattern holds:

\[
\boxed{
 {\cal G}_\phi\succeq0,\qquad
 \operatorname{rank}{\cal G}_\phi=6n-14,\qquad
 \dim\ker{\cal G}_\phi=(n-4)^2.}                  \tag{3}
\]

The physical normal dimension is

\[
 (n-1)^2+1=(6n-14)+(n-4)^2,                      \tag{4}
\]

so (3) accounts for the entire quotient with no unexplained
numerical directions.  Null singular values are at roundoff while
the smallest active singular values remain macroscopically
separated.

## 2. Why this matters

If (3) is proved, the physical remainder has the exact organization

\[
 H_{\rm phys}
 =2I+2{\cal P}_\phi-{\cal L}_\phi^*{\cal L}_\phi \tag{5}
\]

for a response map of real rank \(6n-14\).  The arbitrary
quadratic-size normal problem then reduces to a boundary response
whose size grows only linearly with \(n\).  On the
\((n-4)^2\)-dimensional kernel, the desired sign is automatic:

\[
 H_{\rm phys}=2I+2{\cal P}_\phi\succ0.            \tag{6}
\]

The remaining theorem would be the sharp port contraction

\[
 \|{\cal L}_\phi C\|^2
 \leq2\|C\|^2+2\|\gamma_C\|_{L^2(\mathbb T)}^2,   \tag{7}
\]

where L339 identifies
\({\cal P}_\phi(C)=\|\gamma_C\|_{L^2}^2\).

Equation (7), not the numerical factorization (3), is the actual sign
needed for L336.

## 3. Proof target and falsification rule

The next derivation should start from L127's paired Stein identities
and L339's unitary completion \(U=S+qp^*\).  It should identify an
explicit response \({\cal L}_\phi C\) and prove both:

1. its Gram is exactly (2); and
2. the lossless energy estimate (7).

Any candidate must reproduce the rank \(6n-14\) and the nullity
\((n-4)^2\) before a proof route is opened.  A factor with full
quadratic rank is structurally stale, even if it matches a few small
dimensions.

## 4. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_physical_loss_rank.py \
  --output experiments/gau_wu_physical_loss_rank_s70224.jsonl
```

The checker rebuilds the complete L336 endpoint form, performs
L338's exact zero Schur completion, inserts L339's support Gram, and
audits the positivity and rank in (3).  The 12-record dataset has
SHA-256
`5611d1090336303e28f7a423e91114ad53174e3350c62f4350acfa9f127e5572`.
