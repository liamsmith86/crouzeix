# The L212 target does not dominate the canonical-repair branch

## 1. Result (L295, 2026-07-26)

L294 proves a gap-free energy estimate for L212's canonical target

\[
E_k^{\rm L212}
=12B_kB_k^*-28\Phi(B_k^*B_k)
=-16B_kB_k^*-28G_k.                              \tag{1}
\]

That target must not be identified with, or used as an upper bound
for, the L227/L230--L289 canonical-repair preparation.  Already at a
repeated monomial apex of length \(L\),

\[
\boxed{
E_L^{\rm repair}=-12I,\qquad
E_L^{\rm L212}=-16I.}                            \tag{2}
\]

Consequently

\[
\boxed{
E_L^{\rm repair}-E_L^{\rm L212}=4I\succ0,}       \tag{3}
\]

so the proposed comparison

\[
E_L^{\rm repair}\preceq E_L^{\rm L212}
\]

is false.  The difference cannot be hidden among favorable negative
Grams or prior-flag factors: at the apex every prior transfer row is
zero and the difference is strictly positive on the terminal flag.

This is a scope obstruction, not a failure of L294.  It separates two
valid metric constructions:

1. L227/L283's canonical contraction repair has terminal upper face
   \(-12I\), which is already strictly favorable;
2. L214--L215's exact-axis preparation has L212's base and reaches
   \(-16I\) after the L212 correction.

The live route must therefore do one of the following.

1. Prove the all-grade physical base identity for the
   L214--L215/L212 exact-axis construction itself; or
2. stay in the L227/L285 canonical-repair branch and express its
   mixed transported remainder as a controlled transfer flux, rather
   than compare its complete face with (1).

The second option can still use L294.  More generally, for
\(a\geq0\) and real \(b\),

\[
E_{a,b}=-aB_kB_k^*-bG_k
\]

satisfies, for every \(Y\succeq0\),

\[
\boxed{
\langle Y,E_{a,b}\rangle_+
\le |b|\,\|B_k\|_F
\sqrt{\langle Y,(I-\Phi\Phi^*)Y\rangle}.}         \tag{4}
\]

Thus the useful content of L294 is the flux estimate, not the
incorrect cross-branch comparison.

## 2. Exact apex calculation

At a repeated monomial apex,

\[
B_n=0\quad(n<L),\qquad B_L=U
\]

for a unitary copy matrix \(U\).  Hence

\[
\Phi(K)=UKU^*,\qquad
B_L^*B_L=B_LB_L^*=I,
\]

and therefore

\[
G_L=\Phi(I)-I=0.
\]

Equation (1) gives \(E_L^{\rm L212}=-16I\).
Independently, L283 proves that L227's canonical repair has upper
face \(-12B_LB_L^*=-12I\).  Subtraction proves (3).

The same calculation shows why L212's response cannot repair the
mismatch inside the canonical-repair branch: its response is
\(28G_L=0\) at the apex.

Finally, (4) follows directly from L294:

\[
\langle Y,-aB_kB_k^*\rangle\le0
\]

and

\[
|\langle Y,G_k\rangle|
\le\|B_k\|_F
\sqrt{\langle Y,(I-\Phi\Phi^*)Y\rangle}.
\]

## 3. Consequence for the campaign

Do not continue the instruction

> show the L285/L289 complete physical endpoint is no larger than
> L212's canonical target.

It fails at the most symmetric equality point.  Instead, first
separate the complete negative Gram already supplied by L283 from
the mixed graph-transport remainder.  The next load-bearing
question is whether that remainder has an exact full or weighted
\(G_k\)-flux representation, or an equally direct pairing with
L280's observability defects.  Such a representation would invoke
(4) without demanding the false \(-12I\preceq-16I\) comparison.

## 4. Exact regeneration

Run

```bash
.venv/bin/python -u \
  experiments/repeated_crabb_canonical_target_scope_obstruction.py \
  --output \
  experiments/repeated_crabb_canonical_target_scope_obstruction_s70226.jsonl
```

The exact SymPy checker uses terminal unitaries in multiplicities one
through four.  It verifies both bistochastic identities, \(G_L=0\),
the two endpoint faces in (2), the strict \(4I\) comparison
violation, and the vanishing L212 response.  No floating tolerance is
used.  The tracked dataset has SHA-256

```text
d983482024c5fd0b998ea0028edafcd56d8336c01d6a4fe4aa493ca10ba0264d
```
