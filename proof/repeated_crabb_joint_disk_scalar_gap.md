# A joint disk-endpoint and scalar-channel gap

> **Campaign scope.**  This is an abstract endpoint-angle theorem
> applied to L193's actual canonical disk metric.  It measures disk
> upper loss and scalar-channel loss on the same norming direction.
> It is not yet the full repeated-Crabb neighbourhood theorem.

## 1. Result (L325, 2026-07-26)

Let \(T\) be any sufficiently near repeated block-disk point in
L193's chart, let \(P\) be its normalized canonical Hardy metric, and
put

\[
 C=P^{1/2}TP^{-1/2}.
\]

Then

\[
 I\preceq P\preceq4I,\qquad
 I-C^*C=VV^*,                                     \tag{1}
\]

where \(V\) is an isometry and \(C\) is pure.  The lower eigenspace of
\(P\) is exactly \(\operatorname{ran}V\).  Let \(X\) be an isometry
onto the \(m\)-dimensional spectral cluster of \(P\) near four.  After
one fixed neighbourhood shrink there are constants
\(\beta_-,\beta_+>0\) such that

\[
\begin{aligned}
 P|_{(\operatorname{ran}V)^\perp}&\succeq
 (1+\beta_-)I,\\
 P|_{(\operatorname{ran}X)^\perp}&\preceq
 (4-\beta_+)I.                                    \tag{2}
\end{aligned}
\]

Define the top-cluster disk defect

\[
 D_X=X^*(4P^{-1}-I)X\succeq0.                     \tag{3}
\]

Purity and (1) give the observability Parseval identity

\[
 \sum_{n\ge0}(C^*)^nVV^*C^n=I.                   \tag{4}
\]

For a unit \(k\in\operatorname{ran}X\), define its best constant
scalar Hardy output-line score

\[
 \sigma_X(k)=
 \sup_{\|v\|=1}
 \sum_{n\ge0}
 |\langle C^nk,Vv\rangle|^2
 \le1,                                            \tag{5}
\]

and the **joint scalar disk defect**

\[
 \boxed{
 \Delta(P,C)=
 \min_{\substack{k\in\operatorname{ran}X\\\|k\|=1}}
 \left\{
 \langle(4P^{-1}-I)k,k\rangle+
 1-\sigma_X(k)
 \right\}.}                                       \tag{6}
\]

There is a locally uniform \(c_{L,m}>0\) such that every scalar Schur
function \(f\) obeys

\[
 \boxed{
 \|P^{-1/2}f(C)P^{1/2}\|^2
 \le4-c_{L,m}\Delta(P,C).}                        \tag{7}
\]

The two summands in (6) are evaluated on the **same** top-cluster
direction.  This corrects the invalid shortcut of adding L322's
equality-anchor channel reserve to an independently chosen L197 disk
compression.

On L193's equality manifold, \(\operatorname{ran}X\) is L201's
output-defect space, \(D_X=0\), and (6) reduces to

\[
 \Delta=1-\sigma(B).                              \tag{8}
\]

Thus (7) extends the qualitative linear content of L322 to the whole
canonical disk tube while automatically incorporating every L197
flag layer.  L322 remains sharper on the equality manifold.

## 2. Exact three-loss identity

Fix a unit physical input \(x\), put

\[
 y=P^{1/2}x,\qquad F=f(C),\qquad z=Fy.
\]

Then

\[
\begin{aligned}
 4-\|P^{-1/2}FP^{1/2}x\|^2
={}&
 \langle(4P^{-1}-I)y,y\rangle\\
 &+\{\|y\|^2-\|Fy\|^2\}\\
 &+\langle(I-P^{-1})Fy,Fy\rangle.                 \tag{9}
\end{aligned}
\]

Every term is nonnegative.  This is the exact condition-four
similarity chain with its three losses retained separately.

Let

\[
 q=\frac{\beta_-}{1+\beta_-},\qquad
 \mu=\frac{\beta_+}{4-\beta_+}.                   \tag{10}
\]

Equations (2) give

\[
 I-P^{-1}\succeq
 q(I-VV^*),\qquad
 (4P^{-1}-I)|_{(\operatorname{ran}X)^\perp}
 \succeq\mu I.                                    \tag{11}
\]

If \(a=\|X^*y\|\), \(b=\|(I-XX^*)y\|\), and
\(R^2=a^2+b^2=\|y\|^2\), choose a unit
\(k\in\operatorname{ran}X\) in the direction of \(XX^*y\) when
\(a\ne0\).  Put

\[
 d=\langle(4P^{-1}-I)k,k\rangle,\qquad
 s=\sqrt{\sigma_X(k)}.                            \tag{12}
\]

The first term of (9) is at least

\[
 da^2+\mu b^2.                                    \tag{13}
\]

## 3. The contraction-row angle

Take a unit \(v\) in the direction of \(V^*Fy\).  Set

\[
 t=|\langle Fk,Vv\rangle|\le s.                   \tag{14}
\]

Since \(F^*Vv\) has norm at most one, its component perpendicular to
\(k\) has norm at most \(\sqrt{1-t^2}\).  Therefore

\[
 \|V^*Fy\|
 \le ta+\sqrt{1-t^2}\,b.                          \tag{15}
\]

The last two terms in (9) satisfy

\[
\begin{aligned}
 \|y\|^2-\|Fy\|^2+
 \langle(I-P^{-1})Fy,Fy\rangle
 &\ge q\{\|y\|^2-\|V^*Fy\|^2\}\\
 &\ge q\{R^2-(ta+\sqrt{1-t^2}b)^2\}.              \tag{16}
\end{aligned}
\]

The first inequality uses \(\|Fy\|\le\|y\|\) after splitting the
output into \(\operatorname{ran}V\oplus
(\operatorname{ran}V)^\perp\).

Spend half of the complement loss in (13).  The remaining angle form
in \((a,b)\) is

\[
 \frac\mu2b^2+
 q\{a^2+b^2-(ta+\sqrt{1-t^2}b)^2\}.               \tag{17}
\]

Its \(2\times2\) matrix has

\[
 \operatorname{tr}=q+\frac\mu2,\qquad
 \det=\frac{q\mu}{2}(1-t^2).                      \tag{18}
\]

Consequently

\[
 (17)\ge
 \frac{q\mu}{2q+\mu}(1-t^2)R^2
 \ge
 \frac{q\mu}{2q+\mu}(1-s^2)R^2.                  \tag{19}
\]

No triangle estimate is used; the contraction-row
\(\sqrt{1-t^2}\) is the load-bearing linear-gap mechanism.

## 4. Add the disk face on the same direction

Shrink so that \(0\preceq D_X\preceq d_{\max}I\) for a fixed
\(d_{\max}>0\).  The unspent terms

\[
 da^2+\frac\mu2b^2
\]

obey

\[
 da^2+\frac\mu2b^2
\ge
 \min\left\{1,\frac{\mu}{2d_{\max}}\right\}
 dR^2.                                            \tag{20}
\]

Combine (19)--(20), and use \(R^2=\langle Px,x\rangle\ge1\).
With

\[
 c_{L,m}=
 \min\left\{
 \frac{q\mu}{2q+\mu},
 1,\frac{\mu}{2d_{\max}}
 \right\}>0,                                      \tag{21}
\]

equations (9)--(20) give

\[
4-\|P^{-1/2}FP^{1/2}x\|^2
 \ge c_{L,m}\{d+1-\sigma_X(k)\}
\ge c_{L,m}\Delta(P,C).                          \tag{22}
\]

The first inequality in (22) is the directionwise statement

\[
4-\|P^{-1/2}FP^{1/2}x\|^2
\ge c_{L,m}\Delta_k,\qquad
\Delta_k=d+1-\sigma_X(k).                         \tag{22a}
\]

L327 retains this \(\Delta_k\) and proves that it is quantitatively
equivalent to L324's accumulated response coordinates on the same
scalar branch.  One must not minimize over a different \(k\) before
paying a direction-attached response.

If \(a=0\), (13) supplies the fixed loss
\(\mu R^2\); decreasing \(c_{L,m}\) if necessary gives the same
conclusion.  Taking the supremum over unit \(x\) proves (7).

## 5. Relation with the L197 flag

Let \(G=4I-P\succeq0\).  The top spectral graph over the final copy
level and L194's final Schur graph are related by analytic invertible
maps.  Moreover,

\[
 4P^{-1}-I=P^{-1}(4I-P).                           \tag{23}
\]

Therefore \(D_X\) is analytically congruent, up to uniformly positive
left/right weights, to L197's positive endpoint
\(S_{\rm disk}=-{\cal E}\).  Along every ramified analytic arc,

\[
 \langle D_Xk,k\rangle
 \asymp\sum_j\|\eta_j\|^2,                         \tag{24}
\]

with exactly L324's accumulated quotient columns.  Hence (6) cannot
lose a later disk flag layer or spend an earlier one twice.

## 6. Terminal equality rigidity

Suppose equality is approached in (7) with
\(\Delta=0\).  Then some unit \(k\in\operatorname{ran}X\) has

\[
 Pk=4k,\qquad \sigma_X(k)=1.                      \tag{25}
\]

Equation (4) and equality in the output-line score imply, for a unit
\(v\),

\[
 V^*C^nk=b_n v,\qquad
 \sum_n|b_n|^2=1.                                 \tag{26}
\]

If a scalar Schur function also attains norm two, every inequality in
(9), (14)--(16), and coefficient Cauchy--Schwarz is equality.  Its
coefficient sequence is proportional to \((b_n)\), has \(H^2\) norm
one, and is therefore scalar inner.  In particular,

\[
 f(C)k=Vv.                                        \tag{27}
\]

The observability model from (4) intertwines \(C\) with the backward
shift.  Equations (26)--(27) generate a scalar backward-shift cyclic
subspace.  On L193's full equality manifold, L205's characteristic
kernel argument proves directly that this is a reducing scalar
full-Hardy summand.

For a **partial** L197 top endpoint away from that manifold, the
reducing conclusion is now supplied by L326 rather than by silently
strengthening (26)--(27).  Gau--Wu's 2009 equality theorem for
numerical contractions says that scalar norm equality two forces an
orthogonal reducing summand

\[
 X_\phi S(\phi)X_\phi^{-1},\qquad \phi=zf,
\]

with \(f\) inner and \(f(0)=0\).  Their companion theorem gives this
summand numerical range equal to the disk.  L326's local commutant
argument shows that its dimension is a positive multiple of
\(L+1\), so the orthogonal complement is a strictly
smaller-multiplicity repeated-Crabb problem.

The distinction still matters: \(\sigma_X(k)=1\) alone can describe
an \(H^2\) output-line state.  Scalar innerness and the reducing split
are asserted only after norm equality forces (27) and Gau--Wu applies.
Also, the general L193 disk point need not be nilpotent, so the model
inner function is not silently replaced by a monomial.

## 7. Audit

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/repeated_crabb_joint_disk_scalar_gap.py \
  --output \
  experiments/repeated_crabb_joint_disk_scalar_gap_s70224.jsonl
```

The exact checker verifies the trace/determinant calculation
(17)--(19), the disk allocation (20), and the final constant (21) on
Pythagorean rational angle/input grids and multiple spectral gaps.
The proof is (9)--(22), not the finite audit.

The standard dataset contains 225 records and has SHA-256

```text
bd7e63c2ffaa8ba1bdbe5e6c61cfae1efbb722c107e36eef92b9dfb07cd0c967
```
