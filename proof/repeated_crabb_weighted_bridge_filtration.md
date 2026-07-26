# The affine preparation recurrence preserves elliptic bridge valuation

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

## 1. Result (L315, 2026-07-26)

Retain the balanced pure partial-isometry colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\]

\[
B_k=W^*(S^*)^kV,\qquad
G_k=ES^kF=VB_k^*W^*,
\]

and L305's norm-closed two-sided bridge ideals

\[
{\cal I}_k
=\overline{{\cal A}G_k{\cal A}
            +{\cal A}G_k^*{\cal A}},
\qquad {\cal A}=\operatorname {Alg}(S,S^*).       \tag{1}
\]

Call a finite formal state series

\[
Z(c)=\sum_{d\ge0}c^dZ_d
\]

**elliptically bridge-filtered** if every coefficient has a finite
decomposition

\[
\boxed{
Z_d=\sum_{2k\le d}Z_{d,k},\qquad
Z_{d,k}\in{\cal I}_k.}                            \tag{2}
\]

Membership in (2) is always accompanied by a displayed bounded root
presentation

\[
Z_{d,k}
=\sum_\nu\alpha_\nu p_\nu G_kq_\nu
 \sum_\mu\beta_\mu r_\mu G_k^*t_\mu,             \tag{2a}
\]

or its norm-convergent stable-Stein completion, with the finite
coefficient/word-length sum locally bounded on the fixed
repeated-block neighbourhood.  This presentation—not abstract ideal
membership—is part of the definition.  It prevents (2) from becoming
vacuous when \(\operatorname {Alg}(S,S^*)\) is the full matrix
algebra and supplies the rank-stable \(O(\|B_k\|)\) constants below.
No factor is selected by division.

For a state-to-copy column \(C(c)\), apply (2)--(2a) to
\(C(c)V^*\).

Then the complete L285/L299/L306 affine preparation recurrence
preserves (2):

1. L298's seed pair
   \[
   X(c)=\sum_k r_kc^{2k}X_k,\qquad
   C(c)=\sum_k r_kc^{2k}C_k
   \]
   is filtered, in fact
   \[
   X_k,\ C_kV^*,\ VC_k^*\in{\cal I}_k;            \tag{3}
   \]
2. multiplication on either side by any finite coefficient of the
   raw moving operator or frame preserves every root \(k\);
3. the entire six-term successor
   \[
   \begin{aligned}
   {\cal N}(X,C)
   ={}&-\Delta^*XS-S^*X\Delta-\Delta^*X\Delta\\
      &-HC^*-CH^*-CC^*
   \end{aligned}                                  \tag{4}
   \]
   is filtered;
4. L285's compulsory lower parallel elimination, fixed-base Stein
   inversion, and L305/L311's normalized response selection all
   preserve the same filtration; and
5. therefore every finite prepared successor obtained by iterating
   those operations remains filtered.

Here

\[
\Delta=A-S=O(c),\qquad H=D-V=O(c).
\]

The linear moving terms in (4) actually gain at least one degree:
a grade-\(k\) root entering at degree \(d\ge2k\) reappears there
only at degree at least \(d+1\).

At the upper endpoint, L311 converts every displayed root in
(2a) to

\[
\boxed{
c^d\{K_{d,k}B_k^*+B_kK_{d,k}^*\}
+{\cal M}_S(\widehat C_{d,k}),}                  \tag{5}
\]

with

\[
\begin{aligned}
c^dB_k
 &=c^{d-k}(c^kB_k),\\
\|c^d\widehat C_{d,k}\|_F
 &\le C_{d,k}c^{d-k}\|c^kB_k\|_F .               \tag{6}
\end{aligned}
\]

Thus no recurrence step can expose transfer grade \(k\) before
physical elliptic degree \(2k\).  The displayed presentations make
the factors and response columns polynomial/rank-stable and
automatically satisfy L292's Smith valuation conditions along every
analytic rank-changing arc.

L315 closes the **root valuation** part of the grouped L306/L307
weighting problem.  It does not prove the sign of the hereditary
factors in (5).  The remaining preparation gate is now purely a
multiplier-energy/margin statement:

> show that the \(K_{d,k}\) produced at odd orders factor through
> already active weighted transfer rows, and that at even orders
> their square-completion cost is dominated by retained earlier
> margins after the new direct Gram is isolated.

L313 is the finite physical-word tool for those multipliers.  L307's
\(\theta(1-\theta)C^*C\) reserve and affine leftover must remain in
the final accounting.  L315 neither proves Crouzeix's conjecture nor
closes the repeated-Crabb metric sandwich.

## 2. The L298 seed has the sharp root

L298 writes

\[
\begin{aligned}
D_k
&=-\frac12Q\left\{
 S^kWB_k+
 \sum_{j=1}^{k-1}(S^*)^{k-j}VB_k^*B_j
\right\},\\
C_k&=-\frac12VB_k^*B_k+D_k,                      \tag{7}\\
X_k&=-{\cal G}_S(VB_k^*B_kV^*)
     +{\cal G}_S(VD_k^*+D_kV^*).
\end{aligned}
\]

L314 sharpens L305's cumulative statement to

\[
\begin{aligned}
C_kV^*
={}&-\frac12G_kG_k^*\\
&-\frac12Q\left\{
 S^kG_k^*
 +\sum_{j=1}^{k-1}(S^*)^{k-j}G_kG_j^*
\right\}\in{\cal I}_k.                           \tag{8}
\end{aligned}
\]

Also

\[
VB_k^*B_kV^*=G_kG_k^*\in{\cal I}_k.              \tag{9}
\]

The stable Stein sum

\[
{\cal G}_S(Y)=\sum_{n\ge0}(S^*)^nYS^n
\]

preserves every closed two-sided ideal.  More strongly, applying it
to the explicit roots in (8)--(9) gives a normally convergent rooted
presentation whose coefficient/word-length sum is controlled by a
fixed power moment of \(S\).  Equations (7)--(9) therefore prove
(3) in the stronger sense (2a).  Multiplication by \(c^{2k}\) gives
exactly (2), with no lost degree and no cumulative-root ambiguity.

## 3. Algebraic closure of the filtration

Let \(Z_{d,k}\in{\cal I}_k\).  Since \({\cal I}_k\) is a two-sided
star ideal,

\[
pZ_{d,k}q,\quad q^*Z_{d,k}^*p^*
\in{\cal I}_k                                      \tag{10}
\]

for arbitrary finite state polynomials \(p,q\).  The same remains
true for the displayed presentations: multiplication only lengthens
the prefixes/suffixes and multiplies their finite local bounds.
Norm-convergent stable Stein sums preserve the presentation by
completion.

If \(Y_{e,\ell}\in{\cal I}_\ell\), then

\[
Z_{d,k}Y_{e,\ell}\in{\cal I}_k\cap{\cal I}_\ell. \tag{11}
\]

It may be assigned to either root.  Its total degree obeys

\[
d+e\ge2k+2\ell\ge2k,\ 2\ell.                     \tag{12}
\]

Therefore finite sums, adjoints, products, and fixed-base Stein
inversion preserve (2)--(2a), including their local bounds.

The raw coefficients of \(A,D\) are finite physical polynomials:
L313 proves this for \(A\), while L314's graph factor proves it for
\(D\).  Since \(\Delta,H=O(c)\), multiplying a component
\(c^dZ_{d,k}\) by either motion adds at least one degree and leaves
its root in \({\cal I}_k\).  This proves the claim for the first five
terms of (4).  For the square,

\[
C_dC_e^*
=(C_dV^*)(VC_e^*),                               \tag{13}
\]

and (11)--(12) apply.  Hence all of (4) is filtered.

This is why neither L217's premature raw-frame images nor the
unbounded word length of a Stein expansion causes a physical weight
loss: the already paid bridge sits in \(X\) or \(C\), and a two-sided
multiplier cannot remove it.

## 4. Lower elimination and response selection

At one L285 step, let the Hermitian forcing coefficient
\(R_d\) have a decomposition (2).  Its lower corner is

\[
K_d=V^*R_dV.
\]

The parallel column

\[
C_d^\parallel=-\frac12VK_d
\]

has state lift

\[
\boxed{
C_d^\parallel V^*
=-\frac12ER_dE.}                                  \tag{14}
\]

Equation (10) shows that (14) preserves every root and degree.
The lower-neutral forcing

\[
R_d^0=R_d-VK_dV^*
\]

is therefore filtered, and so is

\[
X_d^0={\cal G}_S(R_d^0).                          \tag{15}
\]

L305 decomposes the upper endpoint of each displayed rooted
component into a copy root plus an explicit response whose lift stays in the same
\({\cal I}_k\).  L310 exposes the literal hereditary factor, and
L311 replaces its possibly nonideal standalone response by a
physical-root-specific column

\[
\widehat C_{d,k}V^*\in{\cal I}_k.                 \tag{16}
\]

Its coefficient and word-length estimate is linear in the bound of
the incoming presentation.  Using (16) for the selected
perpendicular correction and then
applying (15) proves that the complete triangular correction pair is
again filtered.  Induction proves item 5 of Section 1.

The argument does not assert that an arbitrary endpoint-null column
is filtered.  L285's counterexample to naïve ideal invariance remains
valid.  The induction uses the physical L298 seed and L311's
normalized physical-root selection.

## 5. Weighted endpoint extraction

For one finite rooted term

\[
c^dpG_kq,\qquad d\ge2k,
\]

L311 gives directly from the presentation

\[
W^*{\cal G}_S\{pG_kq+(pG_kq)^*\}W
=KB_k^*+B_kK^*+{\cal M}_S(\widehat C),
\]

with

\[
\|\widehat C\|_F
\le C_{p,q,k}\|B_k\|_F,\qquad
\widehat CV^*\in{\cal I}_k.                      \tag{17}
\]

Multiplying by \(c^d\) and extracting \(c^k\) gives (5)--(6).
For a stable Stein sum, L305's completion argument gives the same
conclusion with a finite local constant involving

\[
\Lambda_S=\sum_{n\ge0}(n+1)\|S^n\|.
\]

This constant is uniformly finite after shrinking a fixed
repeated-block neighbourhood.  Since L289 requires only the finite
jet through the terminal grade, no infinite-in-\(k\) summability is
claimed or needed.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_weighted_bridge_filtration.py \
  --output \
  experiments/repeated_crabb_weighted_bridge_filtration_s70226.jsonl
```

The checker constructs mixed L298 corrections, the actual canonical
moving operator/frame jets, and every coefficient of (4).  It tests
unstructured and rank-collapsing colligations as well as complete
delays of depths one through three.  It independently verifies:

1. the sharp seed identities (8)--(9);
2. the fixed-base Stein equations;
3. disappearance of every grade-\(k\) seed and successor when
   \(B_k=0\);
4. absence of a successor before degree \(2k+1\); and
5. bounded normalized successor size under rank collapse.

All 12 tracked records pass.  The largest premature coefficient is
\(2.38\times10^{-15}\), the largest successor attached to a vanished
root is \(2.06\times10^{-14}\), and the dataset has SHA-256

```text
e98a59d70874201afb32272a3e2408e38bda2bc00019f8fdf7cef5ca1ec4b3fb
```

The ideal and degree induction above, rather than the floating audit,
proves L315.
