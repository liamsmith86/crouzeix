# Every word of the unweighted return fan has coefficient one

## 1. Result (L275, 2026-07-25)

Retain L274's complete-delay hypotheses and unweighted L258 renewal
\(\mathfrak U_k=I-Z_{\rm ret}\).  Define

\[
\begin{aligned}
A_{k,j}&=(S^*)^jS^{k+2}(S^*)^{k+2-j},
&&0\leq j\leq k,\\
B_{k,j}&=S^{j+2}(S^*)^{k+2}S^{k-j},
&&0\leq j\leq k-1,\\
G_k&=S(S^*)^{k+2}S^{k+2}S^* .
\end{aligned}                                    \tag{1}
\]

Then the complete fan

\[
\boxed{
\mathcal F_k
=\sum_{j=0}^kA_{k,j}
+\sum_{j=0}^{k-1}B_{k,j}
+G_k}                                            \tag{2}
\]

has coefficient one word by word in the active unweighted face:

\[
\boxed{
[A_{k,j}][c^{2k}]\mathfrak U_k
=[B_{k,j}][c^{2k}]\mathfrak U_k
=[G_k][c^{2k}]\mathfrak U_k=1.}                 \tag{3}
\]

This is an arbitrary-grade induction, not an extrapolation from the
finite fan audits.

The load-bearing associated layer statement is the following
fan-restricted strengthening of L274.  Let

\[
\Delta_k^{\mathfrak U}
=[c^{2k}]\mathfrak U_k
-[c^{2k-2}]\mathfrak U_{k-1}^{\rm tail}.         \tag{4}
\]

Then

\[
\boxed{
\begin{aligned}
[A_{k,j}]\Delta_k^{\mathfrak U}&=0
&& (0\leq j<k),\\
[A_{k,k}]\Delta_k^{\mathfrak U}&=1,\\
[B_{k,j}]\Delta_k^{\mathfrak U}&=0
&& (1\leq j<k),\\
[B_{k,0}]\Delta_k^{\mathfrak U}&=1,\\
[G_k]\Delta_k^{\mathfrak U}&=1,\qquad
[SG_kS^*]\Delta_k^{\mathfrak U}=-1.
\end{aligned}}                                   \tag{5}
\]

L274 proves the four nonzero coefficients in (5).  The inherited
zeros follow from the same first-visit decomposition: after L258 has
resummed the final row, a closed path which does not create the new
outer terminal crossing is the embedded tail path with the identical
L251 analytic channel.  A path which does create that crossing has,
modulo the available lower delays, one of the four outer run patterns
in the last three lines of (5), and no inherited fan pattern.

L275 proves the **coefficients** of the full fan.  It does not prove
that the active unweighted face has no additional nonradial words and
does not prove the lower faces vanish.  L276 subsequently shows that
literal exclusion of every additional word is stronger than the
scalar problem requires: active cyclic radiality through index
\(k+2\), together with lower radiality, is enough.  In particular,
(3) must not be silently upgraded to the full support equality
\([c^{2k}]\mathfrak U_k=\hbox{radial}+\mathcal F_k\).

## 2. Fan-restricted first-visit decomposition

Write L258's renewal as

\[
Z_{\rm ret}
=PZP+\sum_{n\geq0}PZF(FZF)^nFZP.                \tag{6}
\]

Remove the first left wandering layer and put

\[
P=I-F,\qquad T=SP=S^2S^*,\qquad
F_{\rm t}=SFS^* .
\]

The independently balanced tail renewal has the same form as (6)
with \(S,I,F\) replaced by \(T,P,F_{\rm t}\).

Group every word contributing to (6) by its first and last visit to
the removed \(F\)-row.  L251 says that the two visits use one common
analytic tail channel; the physical factor two remains inside that
paired channel.  L258 has already summed every intervening \(FZF\)
loop.  There are therefore two disjoint cases.

### No new outer crossing

Let \(\Pi_k^{\rm old}\) retain precisely the words in (9) below.
The first-visit grouping gives the projected associated identity

\[
\boxed{
\Pi_k^{\rm old}[c^{2k}]\mathfrak U_k
=\iota_k\!\left(
\Pi_{k-1}^{\rm fan}
[c^{2k-2}]\mathfrak U_{k-1}^{\rm tail}
\right),}                                       \tag{7}
\]

where \(\iota_k\) is the tail-word embedding in (11).  To see (7),
take a grouped closed path which does not use the first unavailable
terminal crossing and replace every maximal retained segment by the
corresponding tail segment.  The complete-delay identities

\[
T^h=S^{h+1}S^*,\qquad
(T^*)^h=S(S^*)^{h+1}\qquad(h\geq1)               \tag{8}
\]

make this a coefficient-preserving bijection.  The metric coefficients
on such a segment agree because the active coefficient is deleted in
both systems and every lower orbit is covered by the available delay
relations.  Removing the clean layer removes exactly the paired
entry/exit contour degree \(c^2\), accounting for the coefficient
shift \(2k\mapsto2k-2\).  The common L251 channel preserves the
remaining direct-map scalar coefficient and matrix order.  The words
selected by \(\Pi_k^{\rm old}\) are exactly

\[
A_{k,0},\ldots,A_{k,k-1},
\qquad
B_{k,1},\ldots,B_{k,k-1}.                        \tag{9}
\]

Their three alternating run lengths uniquely determine the deleted
outer layer, so the correspondence has no multiple preimages; its
inverse restores that layer and inserts the same retained \(P\)
between the adjacent tail segments.  Hence it is a bijection of the
already-grouped closed paths, not merely a surjection of reduced
words.

Thus their full and embedded-tail coefficients are identical and
cancel in (4).

This bijection is deliberately stated **after** L258's Schur renewal
has summed the final-row loops.  It is false term by term in the raw
direct/immediate/later split.  For example, at grade four the three
raw full-minus-tail contributions to \(A_{4,0}\) are
\((-32,64,-32)\), not three zeros; only the complete physical sum
vanishes.  This is the same cancellation guardrail as L274.

### A new outer crossing

If the path uses the new crossing, mark its first entrance and last
exit.  Between them, L243's active filtration permits only the
zero/one nonbackground sector, and L258 has already closed every
final-row loop.  Read the unreduced outer \(S,S^*\) runs from those
two marked visits.  Modulo

\[
ES^hF=0\qquad(1\leq h<k),
\]

the only fan or one-shifted-fan normal forms are

\[
A_{k,k},\qquad B_{k,0},\qquad G_k,\qquad SG_kS^*. \tag{10}
\]

Any other placement either:

1. uses an available lower delay and returns to the no-new-crossing
   case;
2. changes an outer run and is not one of the fan words (1); or
3. contains another nonbackground selection, which L243 places above
   degree \(2k\).

This proves the zero statements in (5).  L274's exhaustive source
table supplies the four coefficients in (10), including the
direct/immediate-return cancellation for the crossed pair.  Hence
(5) is proved without a whole-series full/tail identity.  L235's
counterexample to such a whole-series identity is untouched.

## 3. Tail fan embedding

Apply (8) directly to the tail fan at grade \(k-1\).  For every
allowable index,

\[
\boxed{
\begin{aligned}
A_{k-1,j}(T)&=A_{k,j}(S),
&&0\leq j\leq k-1,\\
B_{k-1,j}(T)&=B_{k,j+1}(S),
&&0\leq j\leq k-2,\\
G_{k-1}(T)&=SG_k(S)S^* .
\end{aligned}}                                   \tag{11}
\]

For example, the first line follows by inserting

\[
(T^*)^j=S(S^*)^{j+1},\qquad
T^{k+1}=S^{k+2}S^*
\]

and repeatedly using \(S^*SS^*=S^*\) and \(SS^*S=S\).
The endpoint \(j=0\) is

\[
T^{k+1}(T^*)^{k+1}=S^{k+2}(S^*)^{k+2}.
\]

The second and third lines are the same three-letter reduction.
These are word identities in the complete-delay quotient, not trace
identities.

## 4. Induction

At grade one, direct expansion of L256's exact jet gives

\[
\begin{aligned}
[c^2]\mathfrak U_1
={}&4I-4Q_1-4R_2\\
&+A_{1,0}+A_{1,1}+B_{1,0}+G_1.                 \tag{12}
\end{aligned}
\]

Thus every word of \(\mathcal F_1\) has coefficient one.

Assume (3) at grade \(k-1\).  By definition (4),

\[
[c^{2k}]\mathfrak U_k
=\Delta_k^{\mathfrak U}
+[c^{2k-2}]\mathfrak U_{k-1}^{\rm tail}.         \tag{13}
\]

Equations (11) and the induction hypothesis put coefficient one on

\[
A_{k,0},\ldots,A_{k,k-1},
\quad
B_{k,1},\ldots,B_{k,k-1},
\quad
SG_kS^* .
\]

Equation (5) adds the two missing words \(A_{k,k},B_{k,0}\), adds
\(G_k\), and cancels the inherited \(SG_kS^*\).  The result is
exactly (2), with coefficient one on all \(2k+2\) words.  This proves
(3) for every \(k\).

## 5. Independent exact audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_full_fan_coefficients.py \
  --maximum-physical-grade 6 \
  --maximum-embedding-grade 12 \
  --output \
  experiments/repeated_crabb_full_fan_coefficients_s70225.jsonl
```

The checker independently verifies:

1. the fan-restricted layer law (5) in physical grades two through
   six;
2. the grade-one base (12);
3. coefficient one on every full fan word through grade six;
4. absence of any additional nonradial full or layer word through
   grade six; and
5. all three embedding identities (11) through grade twelve.

Item 4 is a finite falsification audit of the still-open support
theorem, not part of the arbitrary-grade proof above.

The tracked dataset regenerates with SHA-256
`b35c69ab580c8cb019b47f584c99b408c7404a9a94d975f902726b8321bf599b`.
