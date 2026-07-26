# The complete quartic response has a gap-free relative flux bound

## 1. Result (L304, 2026-07-26)

Retain L303's balanced pure partial-isometry colligation and put

\[
E=VV^*,\qquad F=WW^*,\qquad
B=W^*S^*V,\qquad L=BB^*.                         \tag{1}
\]

Fix the partial-retightening scale

\[
\boxed{\theta=\frac12.}                          \tag{2}
\]

After L303's quartic parallel lower neutralization, let \(U_4\) be
the complete **balanced** upper coefficient, including every moving
Schur-graph cross.  L303 identifies its retained commutant class as

\[
P_4=-\frac54L-\frac34L^2.                        \tag{3}
\]

Define the response residual

\[
{\cal E}_4=U_4-P_4.                              \tag{4}
\]

Let

\[
\begin{aligned}
H_Y-SH_YS^*&=WYW^*,\\
Z_Y&=(I-E)H_YV
\end{aligned}                                    \tag{5}
\]

for an arbitrary Hermitian endpoint test \(Y\), and let
\(\Gamma_S\) be the Frobenius operator norm of the stable Stein
inverse \(X\mapsto X-S^*XS\).  Then

\[
\boxed{
|\langle Y,{\cal E}_4\rangle|
\le
\left(
\frac{177}{16}+330\Gamma_S^2+33\Gamma_S
\right)
\|B\|_F\|Z_Y\|_F.}                               \tag{6}
\]

The balanced response map

\[
{\cal M}_S(C)
=W^*{\cal G}_S(VC^*+CV^*)W,\qquad V^*C=0,        \tag{7}
\]

has adjoint \({\cal M}_S^*Y=2Z_Y\).  Hilbert quotient duality
therefore gives a perpendicular quartic column \(C_4^{(0)}\) with

\[
\boxed{
{\cal M}_S(C_4^{(0)})={\cal E}_4,\qquad
\|C_4^{(0)}\|_F
\le
\left(
\frac{177}{32}+165\Gamma_S^2+\frac{33}{2}\Gamma_S
\right)\|B\|_F.}                                 \tag{8}
\]

This is the convenient metric-only cubic representative used to
calculate L302--L303.  Replacing it by L300's admissible bounded
cubic column changes the quartic only by L301's successor response.
Including that change gives an admissible quartic column \(C_4\)
with the nonsharp explicit bound

\[
\boxed{
\|C_4\|_F
\le
\left(
\frac{177}{32}+565\Gamma_S^2+\frac{113}{2}\Gamma_S
\right)\|B\|_F.}                                 \tag{9}
\]

Thus L303's remaining quartic Smith/response gate is closed:

1. the selected column is \(O(\|B_1\|)\), including through rank
   collapse;
2. it preserves L303's lower endpoint and retained upper margin; and
3. along every L197 analytic failure arc, the uniform bound forces
   L292's Smith divisibilities and hence an analytic removable
   selection.

L304 completes the cubic-through-quartic bounded selection for one
fixed partial scale.  It is not an arbitrary-grade induction.  The
next valid step is to extract the relative-commutator recurrence
behind the certificate below, not to compute an isolated higher
grade.

## 2. Relative commutators carry both small factors

The physical first bridge is

\[
\boxed{
G:=ESF=VB^*W^*,\qquad
G^*=FS^*E=WBV^*.}                                \tag{10}
\]

Hence every word sandwich \(pGq\) or \(pG^*q\), with \(p,q\) words
in \(S,S^*\), satisfies

\[
\|pGq\|_F,\ \|pG^*q\|_F\le\|B\|_F.              \tag{11}
\]

The observability Stein equation gives the exact commutators

\[
\boxed{
[H_Y,S]=-SZ_YV^*,\qquad
[H_Y,S^*]=VZ_Y^*S^*.}                            \tag{12}
\]

For \(\ell\in\{S,S^*\}\) and
\(J=pGq\) or \(pG^*q\), cyclicity and (11)--(12) imply

\[
\boxed{
|\operatorname {tr}\{H_Y[\ell,J]\}|
\le\|B\|_F\|Z_Y\|_F.}                            \tag{13}
\]

This is the desired product estimate.  An ordinary cyclic
commutator would give only \(\|Z_Y\|\); inserting the physical bridge
inside the commutator supplies the additional \(\|B_1\|\) needed at
rank collapse.

## 3. Exact 28-term certificate

Use \(s=S\), \(a=S^*\), \(g=G\), and \(h=G^*\).  L303's exact
complete upper trace lift uses L302's finite dual-telescope
replacement for one Green term.  After subtracting (3), denote the
resulting state polynomial by \(D_{\rm eff}\).  Exact
partial-isometry word reduction proves

\[
\boxed{
D_{\rm eff}
=\sum_{\nu=1}^{28}
\alpha_\nu
[\ell_\nu,p_\nu b_\nu q_\nu],\qquad
b_\nu\in\{g,h\}.}                                \tag{14}
\]

The nonzero data are:

| \(\alpha\) | \(\ell\) | \(p\) | \(b\) | \(q\) |
|---:|:---:|:---|:---:|:---|
| \(1/8\) | \(a\) | \(1\) | \(h\) | \(ss\) |
| \(-1/2\) | \(s\) | \(1\) | \(g\) | \(aa\) |
| \(-3/8\) | \(a\) | \(ss\) | \(h\) | \(1\) |
| \(1/4\) | \(s\) | \(1\) | \(h\) | \(ssss\) |
| \(9/16\) | \(s\) | \(1\) | \(h\) | \(ssaa\) |
| \(-1/2\) | \(a\) | \(a\) | \(g\) | \(aaa\) |
| \(-1/4\) | \(s\) | \(ss\) | \(h\) | \(ss\) |
| \(-1/2\) | \(a\) | \(aa\) | \(g\) | \(aa\) |
| \(1/4\) | \(s\) | \(saaa\) | \(g\) | \(1\) |
| \(5/16\) | \(s\) | \(aass\) | \(h\) | \(1\) |
| \(1/4\) | \(a\) | \(aaaa\) | \(g\) | \(1\) |
| \(1/4\) | \(s\) | \(1\) | \(g\) | \(aassss\) |
| \(21/16\) | \(a\) | \(1\) | \(g\) | \(aasssa\) |
| \(1/4\) | \(s\) | \(a\) | \(g\) | \(aaass\) |
| \(1/4\) | \(s\) | \(ass\) | \(h\) | \(sss\) |
| \(1/4\) | \(a\) | \(aaa\) | \(g\) | \(aas\) |
| \(-1/4\) | \(a\) | \(saaa\) | \(g\) | \(aa\) |
| \(1/4\) | \(s\) | \(asss\) | \(h\) | \(ss\) |
| \(1/4\) | \(a\) | \(aaaass\) | \(h\) | \(1\) |
| \(1/4\) | \(a\) | \(1\) | \(h\) | \(sssaaass\) |
| \(1/4\) | \(a\) | \(1\) | \(h\) | \(ssaaasss\) |
| \(-11/16\) | \(s\) | \(1\) | \(g\) | \(aasssaaa\) |
| \(1/4\) | \(s\) | \(1\) | \(g\) | \(aaasssaa\) |
| \(1/4\) | \(a\) | \(saaa\) | \(g\) | \(aaas\) |
| \(1/4\) | \(a\) | \(saaasss\) | \(h\) | \(s\) |
| \(3/8\) | \(a\) | \(sssaaass\) | \(h\) | \(1\) |
| \(1/4\) | \(a\) | \(ssaaasss\) | \(h\) | \(1\) |
| \(-25/16\) | \(s\) | \(aaasssaa\) | \(g\) | \(1\) |

The coefficient norm is

\[
\boxed{\sum_{\nu=1}^{28}|\alpha_\nu|=\frac{177}{16}.}          \tag{15}
\]

Equations (13)--(15) give

\[
|\operatorname {tr}(H_YD_{\rm eff})|
\le\frac{177}{16}\|B\|_F\|Z_Y\|_F.              \tag{16}
\]

The certificate was found by exact rational sparse linear algebra,
but (14) is independently regenerated from its displayed terms; the
optimization is not part of the proof.

## 4. The dual telescope remains gap-free off the cokernel

Let \({\cal N}_3\) be L299's cubic forcing and put

\[
R_3={\cal G}_S({\cal N}_3).
\]

Write \(A_1\) for the first moving operator tangent and

\[
Y_2=2\{S^2+(S^*)^2\},\qquad
Y_2-SY_2S^*=SA_1^*+A_1S^*.                     \tag{17}
\]

L302 replaces the nonlocal quartic term

\[
A_1^*R_3S+S^*R_3A_1
\]

by \({\cal N}_3Y_2\) inside cokernel pairings.  For an arbitrary
\(Y\), the exact error is

\[
\begin{aligned}
\varepsilon_Y={}&
\operatorname {tr}R_3
\{[S,H_Y]A_1^*+[A_1,H_Y]S^*\}\\
&+\operatorname {tr}{\cal N}_3
 {\cal G}_S^\vee([S,H_Y]Y_2S^*),                \tag{18}
\end{aligned}
\]

where \({\cal G}_S^\vee\) is the dual Stein inverse.  Indeed,

\[
H_YY_2-SH_YY_2S^*
=H_Y(SA_1^*+A_1S^*)-[S,H_Y]Y_2S^*,
\]

and Stein adjointness gives (18).

The product form of \(A_1\) gives \(\|A_1\|\le5\).  Its exact reduced
word expansion has commutator length norm \(24\), so

\[
\|[A_1,H_Y]\|_F\le24\|Z_Y\|_F,\qquad
\|Y_2\|\le4.
\]

Consequently

\[
\boxed{
|\varepsilon_Y|
\le33\Gamma_S\|{\cal N}_3\|_F\|Z_Y\|_F.}        \tag{19}
\]

L298--L300 give

\[
\|{\cal N}_3\|_F
\le(20\Gamma_S+2)\|B\|_F.                       \tag{20}
\]

Only half of this term occurs at \(\theta=1/2\).  Combining
(16), (19), and (20) proves (6).

## 5. Restore the admissible cubic representative

The convenient metric-only cubic has norm

\[
\|X_3^{(0)}\|_F
\le\Gamma_S(20\Gamma_S+2)\|B\|_F.               \tag{21}
\]

L300 supplies an admissible perpendicular column \(C_3\) with

\[
\|C_3\|_F\le20\Gamma_S\|B\|_F.
\]

The compatible **difference direction** from the metric-only
representative is

\[
\Delta X_3={\cal G}_S(VC_3^*+C_3V^*),
\]

and obeys

\[
\|\Delta X_3\|_F
\le2\Gamma_S\|C_3\|_F
\le40\Gamma_S^2\|B\|_F.                         \tag{22}
\]

The admissible cubic metric is \(X_3^{(0)}+\Delta X_3\).  The pair
\((\Delta X_3,C_3)\) is exactly the fixed-base Stein difference
between the two representatives.  L301 transports it one order
forward with a quartic response column of norm at most

\[
\begin{aligned}
10\|\Delta X_3\|_F+2\|C_3\|_F
\le(400\Gamma_S^2+40\Gamma_S)\|B\|_F.           \tag{23}
\end{aligned}
\]

Adding (23) to (8) proves (9).  L301 also proves that this successor
has zero lower corner, so L303's lower neutralization is unchanged.

On a fixed repeated-block neighbourhood, \(\Gamma_S\) is uniformly
bounded.  Therefore (9) is uniform and vanishes with \(B_1\).
Along any analytic rank-changing arc, the response equation has a
locally bounded solution.  L292's Smith criterion removes every
apparent division, producing an analytic solution.  L197 curve
selection then rules out a hidden local rank stratum with blow-up.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_retightening_quartic_flux.py \
  --output \
  experiments/repeated_crabb_retightening_quartic_flux_s70226.jsonl
```

The checker reconstructs (14) exactly from the 28 displayed
commutators and obtains zero residual words.  It independently
reconstructs the complete L303 Schur coefficient, the nonlocal
telescope error, both flux inequalities, and a minimum-norm response
column on unstructured, rank-changing, and completely delayed
colligations, including reducible direct sums with a nontrivial
endpoint cokernel.  All 19 records pass.  The tracked dataset has
SHA-256

```text
881bd89b39171ab2603148d54092f605683ece63d65bb61e5394e0e7347f725c
```
