# Removing the moving feedthrough by one Schur step (L170, 2026-07-23)

## 1. Exact construction

Let `theta` be L169's scalar characteristic inner function of degree
`p`, and put

\[
d=\theta(0).
\]

Define its first Schur iterate

\[
\boxed{
{\cal S}\theta(z)
={\theta(z)-d\over z(1-\overline d\,\theta(z))}.
}                                                     \tag{1}
\]

Then `S theta` is a finite Blaschke product of degree `p-1`.
Moreover it depends real-analytically on every analytic matrix/metric
parameter near the Crabb point, where `d=0`.

This gives the canonical way to continue the prepared degree-`L`
Blaschke certificate when the constant characteristic term moves.
It requires no eigenvalue selection and remains analytic when the
eigenvalues split nonanalytically.

## 2. Proof

The disk automorphism

\[
b_d(w)={w-d\over1-\overline d\,w}
\]

maps the unit circle to itself and satisfies
`b_d(theta(0))=0`.  Hence `b_d composed theta` is inner and divisible
by `z`.  Division by `z` removes one Blaschke factor, proving that
(1) is inner of degree `p-1`.

The numerator in (1) has zero constant term by definition, so its
division by `z` is coefficientwise analytic.  The denominator is a
unit near `d=0`.  Complex conjugation of `d` makes the dependence
real-analytic, which is exactly the regularity used by the real
weighted campaign.

At a disk-equality anchor, L169 has

\[
\theta(z)=zB(z),\qquad d=0.
\]

Thus

\[
{\cal S}\theta=B,                                   \tag{2}
\]

so (1) is not a new certificate: it is the intrinsic continuation of
the prepared L123/L149 Blaschke factor.

## 3. Linearization at the Crabb point

At the size-`p` Crabb block,

\[
\theta_0(z)=z^p,\qquad
({\cal S}\theta_0)(z)=z^{p-1}.
\]

For a real analytic variation, differentiation of (1) at `d=0`
gives

\[
\boxed{
\dot B(z)
={\dot\theta(z)-\dot d\over z}
+\dot d\,z^{2p-1}.
}                                                     \tag{3}
\]

For L167's leading grade-`k` eligible normal, put
`m=p+1-k`.  The determinant quotient in L169 gives

\[
\dot\theta(z)
=-kz^{k-1}+kz^{2p+1-k}.                             \tag{4}
\]

### Grade one

Here

\[
\dot\theta=-1+z^{2p},\qquad \dot d=-1.
\]

The two terms in (3) cancel:

\[
\boxed{\dot B=0.}                                   \tag{5}
\]

Thus the moving characteristic feedthrough does not create a
degree-`p-1` Blaschke tangent that could balance the raw endpoint
normal.  This is the intrinsic grade-one exception.

### Every grade at least two

Now `dot d=0`, so (3)--(4) give

\[
\boxed{
\dot B(z)=-kz^{k-2}+kz^{2p-k}.
}                                                     \tag{6}
\]

Relative to `B_0=z^(p-1)`,

\[
{\dot B\over B_0}
=-k(z^{-m}-z^m),                                    \tag{7}
\]

which is a logarithmic-inner tangent.  This is the polynomial-side
version of L168's endpoint identity.  The two derivations are
independent: L168 counts weighted endpoint paths, while (3)--(7) use
only L167's characteristic derivative and the Schur algorithm.

## 4. Nonlinear significance

A115 found that the grade-three characteristic feedthrough becomes
nonzero at later mixed weights.  Equation (1), unlike division by
`z`, remains analytic there and incorporates those terms
automatically.  It is therefore the correct degree-reduction map for
the prospective L149/L169 bridge.

L170 by itself does not prove L163.  At this stage the remaining
statement was:

> In the one-reflection marked sector, identify the frozen prepared
> Blaschke cross transferred by L166 with the real-mean variation
> induced by the Schur iterate (1) of L169's optimized colligation.

Equations (5)--(7) prove the required grade discriminator and remove
root tracking from that task.  They do not by themselves identify
the norming functional.

L171 subsequently separates reflection count and proves that sparse
singular-Hessian support removes every ordinary quadratic term.  It
leaves one direct endpoint/cofactor kernel identity: L156's complete
differentiated norming functional must be evaluated on the
optimized-defect mode left after the inner tangent is removed.

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_leading_endpoint_transfer.py \
  --minimum-length 4 --maximum-length 14 \
  --output \
  experiments/crabb_leading_endpoint_transfer_s70223.jsonl
```

In addition to L168's endpoint identities, the checker constructs
the characteristic determinant quotient, differentiates (1), and
verifies (5)--(7) for every grade through length fourteen.
