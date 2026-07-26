# The complete retightening difference is safe through quartic order

## 1. Result (L303, 2026-07-26)

Retain L302's balanced pure partial-isometry colligation and write

\[
B=B_1=W^*S^*V,\qquad A=B^*B,\qquad L=BB^*.       \tag{1}
\]

Let \(M_{\rm raw}(c)\) be L227's exact canonical repair metric.  Fix

\[
0<\theta\le1,                                    \tag{2}
\]

add \(\theta\) times L298's grade-one retightening direction at order
two, and use the correspondingly scaled L302 cubic and quartic Stein
cancellations.  Denote
the complete lower and physical upper Schur-gap differences from the
raw metric by

\[
\Delta_-(c),\qquad \Delta_+(c).                   \tag{3}
\]

These are complete Schur endpoints, not merely defect
compressions.  Through the quartic,

\[
\begin{aligned}
[c^2]\Delta_-&=-\theta A,&
[c^3]\Delta_-&=0,\\
[c^2]\Delta_+&=4\theta L,&
[c^3]\Delta_+&\in\operatorname {ran}{\cal M}_S .
                                                               \tag{4}
\end{aligned}
\]

L300 supplies a gap-free bounded column canceling the cubic response.
At order four, the exact lower matrix is

\[
\boxed{
[c^4]\Delta_-
=-\left(4\theta+\frac{\theta^2}{2}\right)A
 -\left(\frac{9\theta}{2}
         -\frac{3\theta^2}{4}\right)A^2.}        \tag{5}
\]

For every L206 cokernel test \(Y=Y^*\),

\[
\boxed{
\langle Y,[c^4]\Delta_+\rangle
=(8\theta-2\theta^2)\operatorname {tr}(YL)
 +(10\theta+\theta^2)\operatorname {tr}(YL^2).}  \tag{6}
\]

Now put

\[
J_\theta
=\left(4\theta+\frac{\theta^2}{2}\right)A
 +\left(\frac{9\theta}{2}
         -\frac{3\theta^2}{4}\right)A^2
\succeq0                                        \tag{7}
\]

and add at order four the analytic parallel Stein direction

\[
{\cal G}_S(VJ_\theta V^*),\qquad
C_{\parallel}=\frac12VJ_\theta.                  \tag{8}
\]

It restores the lower endpoint exactly.  Its physical upper-gap
effect has a factor \(-4\).  Consequently the lower-neutral complete
quartic class is

\[
\boxed{
[c^4]\Delta_-\equiv0,\qquad
[c^4]\Delta_+
\equiv
-(8\theta+4\theta^2)L
-(8\theta-4\theta^2)L^2
\pmod{\operatorname {ran}{\cal M}_S}.}           \tag{9}
\]

The upper term in (8) is a cost, but only on the already active
first-transfer range.  If \(U^*B=0\), then

\[
U^*\{(8\theta+4\theta^2)L
 +(8\theta-4\theta^2)L^2\}U=0.                  \tag{10}
\]

Moreover the retained quadratic margin dominates it for small \(c\):

\[
\begin{aligned}
&4\theta c^2L
-c^4\{(8\theta+4\theta^2)L
      +(8\theta-4\theta^2)L^2\}\\
&\quad
=4\theta c^2L^{1/2}
\{I-c^2((2+\theta)I+(2-\theta)L)\}L^{1/2}\\
&\quad\succeq0\qquad(|c|\le1/2),                 \tag{11}
\end{aligned}
\]

because \(0\preceq L\preceq I\).

For the recommended choice \(0<\theta<1\), the raw leading lower face
\(+A\) retains the strict margin

\[
(1-\theta)c^2A.                                  \tag{12}
\]

Thus L298's partial lower retightening, L300's cubic cancellation, L302's
two-ended quartic correction, and **all complete Schur-graph crosses**
are pointwise compatible through order four:

1. the lower endpoint retains a positive active-range margin and
   advances unchanged on \(\ker B\);
2. the next surviving upper flag is unchanged; and
3. the previously active upper range retains a positive margin.

The endpoint choice \(\theta=1\) remains a useful lower-tight
associated-graded limit, but (12) then vanishes; L303 does not claim
that the raw higher lower faces stay positive in that limit.  A fixed
partial scale such as \(\theta=1/2\) avoids that unnecessary debt.

L303 closes the complete quartic **pointwise quotient and margin**
for \(0<\theta<1\).
By itself it does not provide a uniformly bounded analytic
perpendicular column realizing the response congruences through rank
changes.  L304 subsequently closes that L281/L292 flux/Smith gate at
the fixed robust choice \(\theta=1/2\).  Promoting the mechanism to an
arbitrary-grade recurrence remains open.

## 2. Only the second metric cross enters

Let

\[
M_{\rm raw}(c)=I+c^2M_2+O(c^3).                  \tag{13}
\]

The first metric coefficient vanishes.  L227's exact construction
and its canonical slack factor give the finite polynomial

\[
\boxed{
M_2=P_{{\rm bl},2}-{\cal K}_2,}                  \tag{14}
\]

where \(P_{{\rm bl},2}\) is L219's second boundary coefficient and
\({\cal K}_2\) is the lifted slack Schur residual.

L298 changes the second metric coefficient by
\({\mathscr R}\).  Since the constant cross blocks and first metric
coefficients vanish at both equality endpoints, the complete Schur
shorting first sees this change quadratically at order four.  No
cubic metric coefficient enters that square.

For polynomials \(C,D\) and an interior inverse \(G\), put

\[
{\cal X}(C,D;G)
=CGD^*+DGC^*+DGD^*.                              \tag{15}
\]

This is exactly the change from \(CGC^*\) to
\((C+D)G(C+D)^*\).

## 3. Exact lower endpoint

Put \(Q=S^*S=I-E\).  The constant lower-gap interior inverse is

\[
G_-=2Q-\frac23F.                                 \tag{16}
\]

Its raw and retightening second cross blocks are

\[
C_-=EM_2Q,\qquad D_-=\theta E{\mathscr R}Q.      \tag{17}
\]

L302's quartic metric coefficient is
\(-{\cal G}_S({\cal N}_4)\).  Since \(SV=0\), its direct lower
compression is

\[
-V^*{\cal N}_{4,\theta}V
=-6\theta A
 -\left(\frac{5\theta}{2}
         -\frac{\theta^2}{4}\right)A^2.          \tag{18}
\]

Therefore the complete lower difference is the state lift of

\[
-V\left\{
6\theta A
+\left(\frac{5\theta}{2}
       -\frac{\theta^2}{4}\right)A^2
\right\}V^*
-{\cal X}(C_-,D_-;G_-).                          \tag{19}
\]

Substitute L302's two-term formula for \({\mathscr R}\), insert
(14), and reduce only by

\[
SS^*S=S,\qquad S^*SS^*=S^*,\qquad EF=FE=0.
\]

The exact residual after subtracting \(-VJ_\theta V^*\) has zero
words.  This proves (5).  Adding (8) changes the quartic lower
endpoint by \(+J_\theta\); it cannot enter any quartic Schur square
because it first appears at order four.  Hence the first part of
(9) follows.

## 4. Exact upper commutant class

The constant upper-gap interior inverse is

\[
G_+=I-F-\frac23E.                                \tag{20}
\]

The corresponding second cross blocks are

\[
C_+=-FM_2(I-F),\qquad
D_+=-\theta F{\mathscr R}(I-F).                  \tag{21}
\]

The scaled L302 calculation proves that the direct balanced quartic
response has commutant class

\[
\left(4\theta-\frac{\theta^2}{4}\right)L
+\frac{\theta}{2}L^2.                            \tag{22}
\]

The full balanced upper-gap difference is that direct response minus
the cross change:

\[
Q_{4,\theta}-{\cal X}(C_+,D_+;G_+).              \tag{23}
\]

Exact cyclic word collection reduces (21) to

\[
\boxed{
\left(2\theta-\frac{\theta^2}{2}\right)L
+\left(\frac{5\theta}{2}
       +\frac{\theta^2}{4}\right)L^2}            \tag{24}
\]

against every reducing-commutant test.  The physical left-defect
congruence multiplies (24) by four, proving (6).

The parallel direction (8) changes the balanced upper gap by

\[
-W^*{\cal G}_S(VJ_\theta V^*)W.
\]

On every L206 reducing block its weighted pairing is

\[
-\left(4\theta+\frac{\theta^2}{2}\right)
 \operatorname {tr}(YL)
-\left(\frac{9\theta}{2}-\frac{3\theta^2}{4}\right)
 \operatorname {tr}(YL^2).                       \tag{25}
\]

Subtract (25)'s positive polynomial from (24), then multiply by four:

\[
4\left\{
\left(2\theta-\frac{\theta^2}{2}\right)L
+\left(\frac{5\theta}{2}+\frac{\theta^2}{4}\right)L^2
-\left(4\theta+\frac{\theta^2}{2}\right)L
-\left(\frac{9\theta}{2}-\frac{3\theta^2}{4}\right)L^2
\right\}
=-(8\theta+4\theta^2)L
 -(8\theta-4\theta^2)L^2.
\]

This proves the upper identity in (9) by L206 Fredholm duality.

## 5. Independence of the cubic representative

The metric-only cubic in L302 is only a convenient calculator.  L300
replaces it by a bounded perpendicular column that cancels the
complete cubic response.  If two such cubic choices differ by a
fixed-base Stein direction \((X,C)\), L301 proves that their first
moving successor has upper endpoint in
\(\operatorname {ran}{\cal M}_S\).  It also has zero lower corner:

\[
V^*{\cal T}(X,C)V=0.                              \tag{26}
\]

The cubic coefficient cannot enter the quartic Schur square, because
the constant and linear cross coefficients vanish.  Therefore (5),
(6), and (9) are independent of the selected bounded cubic
representative.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_retightening_full_quartic.py \
  --output \
  experiments/repeated_crabb_retightening_full_quartic_s70226.jsonl
```

The exact audit reconstructs both constant interior inverses, the
second raw and retightening crosses, the scaled L302 direct response,
and both Schur-square changes.  Each residual is a polynomial of
degree at most two in \(\theta\); exact checks at
\(\theta=-1,1/2,1\) leave zero lower words and zero upper cyclic
classes, proving the coefficient law.  The 22
numerical records independently reconstruct the full raw, prepared,
and lower-neutral Schur series on unstructured, rank-changing,
completely delayed, and nonscalar reducing direct sums.  The tracked
dataset has SHA-256

```text
7e8cd50fbcb95d99866089c89a0d9b0b684bc32403fa153ae3c822628d8c41d0
```
