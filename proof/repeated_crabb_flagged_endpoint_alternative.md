# A semidefinite alternative for every repeated endpoint flag

> **Closure note (2026-07-25).**  L279 proves the all-grade delayed
> trace \(\operatorname {tr}E_{2k,\mathrm{eff}}
> =-16\|B_k\|_F^2\), closing this note's pointwise range
> obstruction.  L280 subsequently identifies the whole response range
> with the coboundaries of \(I-\Phi\Phi^*\) and supplies an explicit
> analytic state preimage.  A uniform Poisson/Dirichlet bound through
> closing Markov gaps remains open.

## 1. Result (L222, 2026-07-24)

Retain the physical repeated-equality data

\[
P-T^*PT=VV^*,\qquad
P=2I-VV^*+2WW^*,
\]

and L204's real-linear endpoint response

\[
{\cal M}_T(C)
=W^*{\cal G}_T(VC^*+CV^*)W,\qquad V^*C=0,          \tag{1}
\]

where

\[
{\cal G}_T(H)-T^*{\cal G}_T(H)T=H.
\]

Let \(U:\mathbb C^d\to\mathbb C^m\) be an isometry selecting any
left-copy flag, and define

\[
{\cal M}_{T,U}(C)=U^*{\cal M}_T(C)U.                \tag{2}
\]

For a proposed Hermitian endpoint face \(E\in{\rm Herm}_d\), the
following are equivalent:

\[
\boxed{
\begin{array}{ll}
\text{(i)}&
\text{there is a }C,\ V^*C=0,\text{ such that }
E+{\cal M}_{T,U}(C)\prec0;\\[2mm]
\text{(ii)}&
\operatorname {tr}(YE)<0\text{ for every nonzero }Y\succeq0\\
&\text{whose observability Gramian }H_{UYU^*}
\text{ commutes with the balanced partial isometry }S.
\end{array}}                                       \tag{3}
\]

Here

\[
\begin{aligned}
S&=P^{1/2}TP^{-1/2},\\
H_{\widehat Y}
&=\sum_{n\ge0}S^nW\widehat YW^*(S^*)^n,\qquad
\widehat Y=UYU^* .                                  \tag{4}
\end{aligned}
\]

Equivalently, the only obstructions to making a flagged face
strictly negative are positive reducing-copy weights:

\[
\boxed{
H_{\widehat Y}S=SH_{\widehat Y}.}                  \tag{5}
\]

This is a strict semidefinite alternative, not a numerical
surjectivity assertion.  It has two consequences for the repeated
elliptic campaign.

1. On an irreducible proper flag, the compressed response is
   typically onto, so **any** Hermitian face can be repaired.
2. If it is not onto, every positive separator in (3) comes from a
   genuinely reducing colligation summand.  On that summand the
   earlier transfer rows vanish whenever
   \(U^*B_1=\cdots=U^*B_{k-1}=0\).  The unresolved matrix-valued
   **pointwise range obstruction** therefore reduces to the scalar
   trace of the fully delayed physical face on reducing summands.

L222 removes the need to identify every partial-flag face with the
exact matrix \(-16B_kB_k^*\) for pointwise feasibility.  It does
**not** prove the remaining fully delayed trace sign, nor does the
separation theorem supply a bounded analytic correction as the
commutant dimension changes.  The all-grade trace sign (grades two
and three are L214--L215) and a bounded pathwise/Schur selection are
the two remaining steps needed to turn this alternative into a local
metric theorem.

## 2. The compressed adjoint

Give every Hermitian matrix space the real trace inner product.
For \(Y=Y^*\in M_d\), cyclicity gives

\[
\begin{aligned}
\langle Y,{\cal M}_{T,U}(C)\rangle
&=\langle UYU^*,{\cal M}_T(C)\rangle.               \tag{6}
\end{aligned}
\]

Consequently

\[
\boxed{
{\cal M}_{T,U}^*(Y)
={\cal M}_T^*(UYU^*).}                              \tag{7}
\]

L204 computes the full adjoint.  If \({\cal Z}_{\widehat Y}\) solves

\[
{\cal Z}_{\widehat Y}
-T{\cal Z}_{\widehat Y}T^*
=W\widehat YW^*,
\]

then, up to an irrelevant positive scalar normalization,

\[
{\cal M}_T^*(\widehat Y)
=(I-VV^*){\cal Z}_{\widehat Y}V.                   \tag{8}
\]

Thus

\[
Y\perp\operatorname {ran}{\cal M}_{T,U}
\quad\Longleftrightarrow\quad
(I-VV^*){\cal Z}_{UYU^*}V=0.                       \tag{9}
\]

L206 proves, with no genericity assumption, that (9) is equivalent
to (5).  It also gives Hermitian \(A\) such that

\[
\boxed{
H_{\widehat Y}W=W\widehat Y,\qquad
H_{\widehat Y}V=VA,\qquad
\widehat YB_n=B_nA\quad(n\ge1).}                   \tag{10}
\]

Therefore the annihilator in (3) is exactly the self-adjoint
commutant localized to the selected left flag.

## 3. Semidefinite separation

Put

\[
{\cal R}=\operatorname {ran}{\cal M}_{T,U}
\subset{\rm Herm}_d.
\]

This is a finite-dimensional real linear subspace and hence closed.
The affine space \(E+{\cal R}\) meets the open negative cone if and
only if it cannot be strongly separated from that cone.

Suppose first that it does not meet the cone.  The separating
hyperplane theorem supplies a nonzero Hermitian \(Y\) satisfying

\[
\langle Y,R\rangle=0\quad(R\in{\cal R}),\qquad
\langle Y,N\rangle\le0\quad(N\preceq0),\qquad
\langle Y,E\rangle\ge0.                            \tag{11}
\]

The middle condition is equivalent to \(Y\succeq0\).  Thus the
separator is already a nonzero positive matrix in
\({\cal R}^{\perp}\), and it obeys

\[
\operatorname {tr}(YE)\ge0.                        \tag{12}
\]

Conversely, (12) is incompatible with
\(E+R\prec0\), because a nonzero positive matrix pairs strictly
negatively with every negative-definite matrix.  Hence

\[
(E+{\cal R})\cap\{N:N\prec0\}\ne\varnothing
\quad\Longleftrightarrow\quad
\operatorname {tr}(YE)<0
\]

for every nonzero \(Y\succeq0\) in \({\cal R}^{\perp}\).  Combining
this with (7)--(9) proves (3).

## 4. Localization to fully delayed reducing summands

Assume a flag is still open at transfer grade \(k\):

\[
U^*B_j=0\qquad(1\le j<k).                           \tag{13}
\]

Let \(Y\succeq0\) be a separator in (3), and set
\(\widehat Y=UYU^*\).  By (5), every spectral projection of
\(H_{\widehat Y}\) reduces both \(S\) and \(S^*\).  Equations (10)
show that its left and right defect spaces are the corresponding
spectral subspaces of \(\widehat Y\) and \(A\).

On every positive spectral subspace of \(\widehat Y\), (13) says
that the first \(k-1\) left transfer rows vanish.  The intertwining
relations (10) show that the matching right spectral space is
invariant as well.  Hence that reducing colligation summand has a
complete delay of length \(k-1\), in precisely the sense of
L209--L221.

It follows that the separator tests in (3) do not require a new
noncommutative partial-flag formula.  They are weighted sums of trace
tests on completely delayed reducing summands.  A strictly negative
fully delayed trace on every summand where \(B_k\ne0\) makes (3)
strict.  If \(B_k\) also vanishes, that summand advances to the next
flag.  L201's invertible terminal \(B_L\) makes this process finite.

This conclusion is pointwise.  Near a rank jump, the norm of a
chosen preimage can in principle diverge even though the affine range
meets the negative cone at every nonzero parameter.  L197's analytic
flag or L220's Schur coordinates must still be used to prove bounded
selection along approaching paths.

## 5. What remains

The next target is the gauge-invariant higher-contact trace law

\[
\operatorname {tr}E_{2k,\mathrm{eff}}
=-16\|B_k\|_F^2                                    \tag{14}
\]

on a completely delayed reducing summand, after lower tightening
and previous upper Schur elimination.  L203 proves (14) for \(k=1\);
L214--L215 prove the stronger matrix identity for \(k=2,3\).

A proof of (14) for all \(k\), together with L222, is sufficient for
the pointwise sign/range part of the repeated elliptic flag.  One
must additionally prove bounded pathwise selection through rank
changes before invoking analytic remainder domination.  Exact
equality with a prescribed matrix Gram at every partial flag is no
longer required.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_flagged_endpoint_alternative.py \
  --output \
  experiments/repeated_crabb_flagged_endpoint_alternative_s70224.jsonl
```

The deterministic audit checks full, proper irreducible, and proper
reducing flags.  It verifies the compressed-adjoint identity,
commutant localization of every numerical annihilator, full response
on irreducible proper flags, a known strictly feasible correction,
and the positive-separator obstruction on reducing flags.  The
finite-dimensional separation argument above, not the floating
audit, proves L222.  The tracked dataset SHA-256 is
`e04f7a4fdb8d97a6e843f434e882d840839c4f624721dedbb7e78c81f8383728`.
