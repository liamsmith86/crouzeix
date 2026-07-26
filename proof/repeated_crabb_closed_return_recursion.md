# One-delay covariance is trace-only and retains a deep divergence

## 1. Result (L269, 2026-07-25)

L257's associated one-delay recursion is not an operator identity,
even after L258's final-row renewal.  At the first nontrivial delayed
grade, its operator difference contains an irreducible unilateral
boundary term of depth \(k+1\).  That term survives in the operator
quotient but is paired with its shifted copy as a trace-zero Stein
divergence.

This disproves the stronger version of A213's second stop condition
which asked that **no** boundary term deeper than \(k\) survive.
The corrected goal is:

1. isolate every depth-\(k+1\) term into a delayed trace-zero
   divergence; and
2. apply L266 only to the remaining depth-at-most-\(k\) return.

The claim is exact already at grade two.  The same structure is
verified by exact rational word arithmetic through grade five, but
that finite continuation is evidence rather than an all-grade
theorem.

## 2. Exact grade-two face

Use \(s=S\), \(a=S^*\), and

\[
Q_j=a^js^j,\qquad R_j=s^ja^j.
\]

Impose the first complete delay \(B_1=0\).  Let
\(\mathfrak D_2\) be the degree-four L258 closed-return face of the
full system, and let \(\mathfrak D_1^{\rm tail}\) be the embedded
degree-two face of the independently balanced deflated tail.
Exact word reduction gives

\[
\boxed{
\begin{aligned}
\mathfrak D_2-\mathfrak D_1^{\rm tail}
={}&-2I+5Q_1-R_1-3Q_2+R_2\\
&+G_2-SG_2S^*,
\end{aligned}}                                   \tag{1}
\]

where

\[
G_2=S(S^*)^4S^4S^*.                              \tag{2}
\]

Both full and tail closed-return series have every earlier
coefficient equal to zero.  The seven words in (1) are nonzero
canonical representatives in the \(B_1=0\) operator quotient.

On the unilateral half-line, (2) reduces to the boundary projection
\((L^*)^3L^3\), so it has boundary depth three.  Since \(k=2\), this
is exactly the first depth excluded by L266's low-depth theorem.
The shifted term has the same trace but is a distinct operator word.
Thus the deep remainder cannot be removed by strengthening L266 or
by claiming that L243's contour products are all shallow.

## 3. Why the trace still transports

The radial part of (1) has zero trace.  Indeed,

\[
\operatorname {tr}R_j=\operatorname {tr}Q_j,
\]

and L252 gives, under \(B_1=0\),

\[
\operatorname {tr}Q_1=n-m,\qquad
\operatorname {tr}Q_2=n-2m.
\]

Therefore

\[
\begin{aligned}
\operatorname {tr}
(-2I+5Q_1-R_1-3Q_2+R_2)
&=-2n+4\operatorname {tr}Q_1
  -2\operatorname {tr}Q_2\\
&=0.                                             \tag{3}
\end{aligned}
\]

The delayed cyclic quotient also reduces both \(G_2\) and
\(SG_2S^*\) to the same representative \(Q_4\).  Hence

\[
\operatorname {tr}(G_2-SG_2S^*)=0.               \tag{4}
\]

Equations (3)--(4) prove

\[
\operatorname {tr}\mathfrak D_2
=\operatorname {tr}\mathfrak D_1^{\rm tail}
\]

without asserting operator covariance.  This is the exact mechanism
behind L268's equal traces and unequal ranks.

## 4. Finite all-grade diagnostic

For every audited \(2\le k\le5\), exact reduction gives

\[
\mathfrak D_k-\mathfrak D_{k-1}^{\rm tail}
=P_k(I,Q_1,\ldots,Q_k,R_1,R_2)
+G_k-SG_kS^*,                                    \tag{5}
\]

where

\[
G_k=S(S^*)^{k+2}S^{k+2}S^*.                      \tag{6}
\]

The first term in (5) is radial, the second has unilateral boundary
depth \(k+1\), and the complete cyclic difference is zero.  The
operator differences contain respectively `7,8,9,8` canonical words.
The full active faces themselves contain only `10,11,11,12` words,
far fewer than the pre-renewal mass expansions.

The repeated form (5) is a strong guide, not an induction.  An
all-grade proof must derive (5) from L243's zero/one-kernel sectors
and L251/L258's paired renewal.  L270 subsequently proved that the
exact monomial axes force the two scalar moment identities once (5)
is known, so the remaining gate is the structural decomposition
itself.  Computing more grades would not close that step.

## 5. Exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_closed_return_recursion.py \
  --maximum-grade 5 \
  --output \
  experiments/repeated_crabb_closed_return_recursion_s70225.jsonl
```

The checker now reuses L252's exact quotient engine to construct
L258's closed-return defect directly.  It independently constructs
the full and embedded-tail series, checks every earlier vanishing,
verifies the exact seven-word formula (1), removes (6), checks that
the remainder is radial, and cyclically reduces the complete
difference to zero.

The tracked dataset SHA-256 is
`0398cee7f133b13ab4bc1f332e5d8e030f8ce86e92929e547cc9b30cce48a988`.
Both modified scripts pass Ruff and `py_compile`.

## 6. Corrected live gate

The A213 proof cannot require literal boundary depth at most \(k\)
for the complete operator face.  The exact grade-two face disproves
that statement.  The viable replacement is

\[
\text{active return}
=\text{L266-shallow part}
+(Y_k-SY_kS^*)
+\text{other trace-null commutators},            \tag{7}
\]

with the delayed boundary flux of every displayed divergence proved
to vanish.

Thus the next all-grade calculation should seek the universal
recurrence for \(P_k\) and \(G_k\) in (5), not try to make the deep
terms disappear.  If (5) is proved, the one-delay scalar recursion
follows immediately and iteration ends at L256.
