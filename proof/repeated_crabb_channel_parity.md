# Every doubled-Hardy closed endpoint is even in the transfer channel

## 1. Result (L316, 2026-07-26)

Retain L255's doubled-Hardy representation of the physical ellipse
pencil.  Its only colligation-dependent object is the transfer
Hankel channel

\[
{\cal H}={\cal O}_R{\cal O}_L^*,\qquad
{\cal P}=P_{\{0\}}={\cal J}_0{\cal J}_0^*,
\]

where \({\cal J}_0:\mathbb C^m\to H^2(\mathbb C^m)\) is the
row-zero inclusion.

Let

\[
{\cal J}=
\begin{bmatrix}I&0\\0&-I\end{bmatrix}             \tag{1}
\]

be the orientation-sign operator.  If
\(\mathbb X({\cal H})\) denotes L255's complete doubled pencil, then

\[
\boxed{
\mathbb X(-{\cal H})
={\cal J}\mathbb X({\cal H}){\cal J}.}            \tag{2}
\]

Consequently every analytic functional-calculus expression obeys

\[
\boxed{
f(\mathbb X(-{\cal H}))
={\cal J}f(\mathbb X({\cal H})){\cal J}.}         \tag{3}
\]

The identity persists under sums, products, adjoints, analytic
inverses, block-diagonal universal metric weights, and Schur
complements which eliminate a whole orientation block.  Hence:

1. every diagonal orientation block of a closed physical expression
   is an even series in \({\cal H},{\cal H}^*\);
2. every off-diagonal block is odd; and
3. every same-endpoint sandwich at \(V\) or \(W\) is even, after its
   exact Hardy embedding is included.

The endpoint embeddings are

\[
\boxed{
\begin{aligned}
\begin{bmatrix}{\cal O}_RV\\{\cal O}_LV\end{bmatrix}
&=
\begin{bmatrix}{\cal J}_0\\{\cal H}^*{\cal J}_0\end{bmatrix},\\
\begin{bmatrix}{\cal O}_RW\\{\cal O}_LW\end{bmatrix}
&=
\begin{bmatrix}{\cal H}{\cal J}_0\\{\cal J}_0\end{bmatrix}.
\end{aligned}}                                    \tag{4}
\]

Under \({\cal H}\mapsto-{\cal H}\), the first column in (4)
transforms by \({\cal J}\), while the second transforms by
\(-{\cal J}\).  The harmless global minus sign cancels in a
quadratic sandwich.

Therefore a nonconstant retained closed endpoint cannot contain
exactly one transfer-Hankel factor:

\[
\boxed{\text{channel degree }1\text{ is absent.}} \tag{5}
\]

Every nonresponse endpoint term is either channel-independent
(the universal elliptic-axis background) or contains at least two
physical transfer channels.  This is compatible with L309's
nonzero odd **elliptic-degree** example: odd powers of \(c\) need
not have odd transfer-channel degree.

There is also a coordinate-free gauge interpretation.  Replacing
the left defect frame by \(W\mapsto-W\) leaves \(S,E,F\), every state
metric, and every physical state operation unchanged, while

\[
B_j\mapsto-B_j,\qquad {\cal H}\mapsto-{\cal H}.
\]

Both occurrences of \(W\) in a same-endpoint sandwich cancel their
signs.  Thus (5) is not an artifact of a particular Hardy window.
L298's columns and L305/L311's physical normalized selections are
made from gauge-invariant state products, so the same parity applies
to the **grouped closed endpoint** of the prepared recurrence.  It
need not apply to a particular one-root hereditary summand chosen
during L310's decomposition; such summands may occur, but their total
odd-channel part must cancel or lie in the grouped response.

L316 closes the parity part of L315's multiplier-energy problem.  It
does not by itself prove that a two-channel term appearing at
elliptic degree \(d\) has the sharp weighted index inequality
\(j+k\le d\), nor that its coefficient is dominated by the retained
Gram margin.  That final two-channel weight/energy estimate remains
the live gate.

## 2. Pencil covariance

L255 gives

\[
\mathbb X({\cal H})
=
\begin{bmatrix}
X_{RR}&X_{RL}\\
X_{LR}&X_{LL}
\end{bmatrix},                                    \tag{6}
\]

where:

1. \(X_{RR}\) is the universal right half-line pencil plus
   \(c{\cal H}{\cal P}{\cal S}{\cal H}^*{\cal P}\);
2. \(X_{LL}\) is universal;
3. \(X_{RL}\) is linear in \({\cal H}\); and
4. \(X_{LR}\) is linear in \({\cal H}^*\).

Thus the diagonal blocks are unchanged by
\({\cal H}\mapsto-{\cal H}\), while both off-diagonal blocks change
sign.  Conjugation by (1) has exactly the same effect, proving (2).

For a polynomial \(f\), equation (3) follows by taking powers in
(2).  Analytic functional calculus follows by local uniform
polynomial approximation or by the resolvent:

\[
(zI-\mathbb X(-{\cal H}))^{-1}
={\cal J}(zI-\mathbb X({\cal H}))^{-1}{\cal J}.   \tag{7}
\]

No coefficient, orientation, or theta term is separated.  In
particular, L260's warning remains in force: the two physical
orientations are combined before parity is used.

## 3. Closure under physical operations

Call a block expression \({\cal J}\)-covariant if

\[
Y(-{\cal H})={\cal J}Y({\cal H}){\cal J}.         \tag{8}
\]

If \(Y,Z\) are covariant, so are

\[
Y+Z,\quad YZ,\quad Y^*,\quad Y^{-1}              \tag{9}
\]

whenever the inverse exists.  A universal block-diagonal weight
commutes with \({\cal J}\), so inserting L236's \(D_R,D_L\) does
not alter (8).

Write a covariant Hermitian block as

\[
Y=
\begin{bmatrix}A&B\\B^*&D\end{bmatrix}.
\]

Then \(A,D\) are even and \(B\) is odd.  If \(D\) is invertible, its
Schur complement

\[
A-BD^{-1}B^*                                     \tag{10}
\]

is even.  The same holds with the orientations reversed.  Iterated
closed returns and graph eliminations inherit the property because
they are made only from (9)--(10).

This proves parity for the full grouped operation, not merely for a
finite truncation of its Neumann expansion.

## 4. Endpoint covariance

The right Hardy analysis of \(V\) has only row zero:

\[
{\cal O}_RV={\cal J}_0.
\]

The opposite frame gives

\[
({\cal O}_LV)_j=W^*(S^*)^jV=B_j,
\]

which is exactly row \(j\) of
\({\cal H}^*{\cal J}_0\).  Likewise

\[
{\cal O}_LW={\cal J}_0,\qquad
{\cal O}_RW={\cal H}{\cal J}_0.
\]

Thus (4) holds exactly.

Let \(v_E({\cal H})\) denote either endpoint column.  It obeys

\[
v_E(-{\cal H})
=\epsilon_E{\cal J}v_E({\cal H}),\qquad
\epsilon_E\in\{1,-1\}.                            \tag{11}
\]

For every covariant expression \(Y\),

\[
\begin{aligned}
v_E(-{\cal H})^*Y(-{\cal H})v_E(-{\cal H})
&=v_E({\cal H})^*
 {\cal J}{\cal J}Y({\cal H})
 {\cal J}{\cal J}v_E({\cal H})\\
&=v_E({\cal H})^*Y({\cal H})v_E({\cal H}).
\end{aligned}                                    \tag{12}
\]

Therefore the sandwich is even.  Its Taylor/free-word expansion has
no monomial containing an odd total number of
\({\cal H},{\cal H}^*\), proving (5).

## 5. Consequence for the preparation recurrence

L315 proves that every retained root at transfer grade \(k\) has a
displayed bounded bridge presentation and appears no earlier than
elliptic degree \(2k\).  L316 adds the complementary fact that a
closed retained endpoint cannot stop after exposing only that one
channel.

Thus the dangerous alternative

\[
c^d\{KB_k^*+B_kK^*\}
\quad\text{with }K\text{ genuinely channel-independent}       \tag{13}
\]

is absent from the grouped closed endpoint.  Either (13) is a
response/gauge term which disappears from the closed same-orientation
quantity, or \(K\) contains at least one additional transfer
channel.

What remains is quantitative rather than parity-based:

1. expose that second channel with its correct index and \(c\)-weight;
2. write the result as a bounded block quadratic form in
   \([cB_1,c^2B_2,\ldots]\); and
3. show that every nonprincipal block has an extra power of \(c\),
   while the principal diagonal is exactly L283's positive direct
   Gram.

L243's reflection-power filtration and L313's paid reverse-edge
normal form are the appropriate tools for this last step.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_channel_parity.py \
  --output \
  experiments/repeated_crabb_channel_parity_s70226.jsonl
```

The checker uses finite Hardy windows with exact terminal rows
discarded.  For random noncommuting Hankel channels it independently
verifies:

1. the pencil covariance (2);
2. covariance of polynomial functional calculus and resolvents;
3. even diagonal and odd off-diagonal blocks;
4. covariance of both endpoint columns in (4);
5. evenness of their closed metric sandwiches; and
6. evenness after both orientation Schur complements.

All 12 tracked records pass with every reported sign-covariance,
endpoint-evenness, and Schur-evenness residual exactly zero in
floating arithmetic.  The dataset has SHA-256

```text
03b23df2857cfd71edf69da66d2dd7df1a8c72a3c667e9a378df8663ed7a5711
```

The sign-conjugation identities above, not the floating audit, prove
L316.
