# Canonical normal form of the Gau--Wu endpoint residuals

> **Status and scope.**  The elimination of the free metric row and
> all Blaschke-zero velocities below is exact at every finite
> nondegenerate Gau--Wu model.  It reduces L352's mixed intertwiner to
> one purely physical map with only \(n-2\) complex output
> coordinates.  The final factorization of that map through the
> conformal shape is numerical evidence only.  The common pure
> physical term in L342 also remains open, so this note proves neither
> phase covariance nor the Hessian sign.

## 1. The two endpoint spaces

Retain
\[
 {\cal K}_{zf}=\mathbb Cp\oplus{\cal M}\oplus\mathbb Cq,
 \qquad \dim{\cal M}=n-2,
\]
and L350's canonical conjugation \({\cal C}\), which exchanges
\(p\) and \(q\).  Its restriction
\[
 \Gamma:p^\perp\longrightarrow q^\perp
 \tag{1}
\]
is antiunitary.  Write L342's two endpoint residuals as
\[
 r_-\in p^\perp,\qquad r_+\in q^\perp.
 \tag{2}
\]
Define their conjugated sum
\[
 d=r_++\Gamma r_-.
 \tag{3}
\]

## 2. The metric graph cancels exactly

At zero physical and zero motion, the first Blaschke image is zero.
The free boundary row \(x\) produces the selfadjoint
truncated-Toeplitz metric extension from L342.  Canonical
conjugation exchanges its two active endpoint columns.  Substitution
in the two active-kernel motion formulas gives
\[
 a_+(x)=-\Gamma a_-(x).
 \tag{4}
\]
Since \(r_\pm=-2a_\pm\) in this sector,
\[
\boxed{r_+(x)=-\Gamma r_-(x),\qquad d(x)=0.}
\tag{5}
\]
Thus the complete metric residual range is the antiunitary graph in
(5).  This is stronger than recording only its signed phase in L352.

## 3. The zero graph is corrected by tracking

For a pure moving-inner direction, use L338's unique coordinate
\[
 h_k=(zf)k-J_fk,\qquad k\in{\cal K}_f.
 \tag{6}
\]
L338 gives
\[
 r_-=-J_fk,\qquad r_+=-zk.                       \tag{7}
\]
The two canonical conjugations obey
\[
 {\cal C}(J_fk)=zk.                              \tag{8}
\]
Consequently the pure-zero residual range is the opposite graph,
\[
 r_+=\Gamma r_-,\qquad d=-2zk.                  \tag{9}
\]

Let \(\tau\in\mathbb C^{n-1}\) be L351's endpoint/simple-root
tracking coordinate.  Its restriction to pure zero motion is the
complex isomorphism
\[
 \tau_0=-2f'(0)v_0,\qquad
 \tau_j=-v_j\quad(1\le j\le n-2).                \tag{10}
\]
Let \(k_\tau\) be the L338 kernel coordinate of the unique pure-zero
motion having tracking value \(\tau\).  Since \(k_\tau\) is
complex-linear in \(\overline\tau\), (9) defines a fixed
complex-linear map
\[
 {\cal T}_\phi:\overline{\mathbb C^{n-1}}\to q^\perp,
 \qquad
 {\cal T}_\phi\overline\tau=-2zk_\tau,            \tag{11}
\]
with the harmless fixed Gau--Wu diagonal included when ordinary
matrix coordinates are used.

For a general physical/metric/zero first direction, define
\[
\boxed{
 \Delta_\phi
 =r_++\Gamma r_- -{\cal T}_\phi\overline\tau.}
\tag{12}
\]
Linearity, (5), and (9)--(11) prove that \(\Delta_\phi\) is
independent of both the free metric row and every zero velocity.  It
is a purely physical response.

There is one further exact cancellation.  The \(p\)-component of
\(d\) is the conjugate endpoint trace:
\[
 \langle p,d\rangle=\overline{\tau_0}.           \tag{13}
\]
The same component of (11) is
\(\overline{\tau_0}\), by (7), (10), and endpoint exchange.  Hence
\[
\boxed{\Delta_\phi(C)\in{\cal M}.}               \tag{14}
\]
The mixed residual problem has therefore fallen from two
\((n-1)\)-component endpoint vectors to one physical
\((n-2)\)-component vector.

## 4. Exact remaining mixed-intertwiner criterion

Let
\[
 a={\cal A}_\phi(C)
 =2(\widehat{s_C}(2),\ldots,\widehat{s_C}(n))
 \in\mathbb C^{n-1}
 \tag{15}
\]
be A292's exact conformal shape coordinate.  Equation (12) rewrites
the full endpoint relation as
\[
\boxed{
 r_+=-\Gamma r_-+
 {\cal T}_\phi\overline\tau+\Delta_\phi(C).}
\tag{16}
\]

The physical part of L352's desired signed intertwiner is now
equivalent to the single statement
\[
\boxed{
 \Delta_\phi(C)={\cal B}_\phi a
 \quad\hbox{for a complex-linear }
 {\cal B}_\phi:\mathbb C^{n-1}\to{\cal M}.}
\tag{17}
\]
Indeed, under
\[
 (a,\tau,r_-,r_+)\longmapsto
 (ia,-i\tau,-ir_-,ir_+),
\]
antilinearity of \(\Gamma\) and complex linearity of
\({\cal T}_\phi,{\cal B}_\phi\) make every term in (16) acquire the
same factor \(i\).  Conversely, signed invariance of the residual
range and surjectivity of L351's active chart force the induced
physical map to be complex linear in \(a\).

Equation (17) has two concrete pieces:

1. \(\Delta_\phi\) must vanish on L343's disk fibre
   \(\ker{\cal A}_\phi\);
2. the descended real-linear map on the shape quotient must have no
   \(\overline a\) part.

These are still proof debts.  They should be derived from L344's
fixed Hardy port, not from numerical row-space inclusion.

Even after (17), L352's second debt remains: the common pure
physical coefficient shared by L342's scalar and similarity forms
must separately be shown to be type \((1,1)\).

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_canonical_residual_normal_form.py
```

The checker constructs \(\Gamma\) independently from L350's
canonical conjugation, verifies the metric graph (5), constructs
\({\cal T}_\phi\) only from the pure-zero tracking restriction, and
checks (12)--(14) on the complete physical/metric/zero basis.

Two models in every dimension \(3,\ldots,8\) give metric-graph,
zero-normal-form, endpoint-trace, and corrected-homogeneous
residuals below \(9.0\cdot10^{-15}\).

As evidence only, the checker also fits (17).  Its relative
factorization and disk-kernel residuals are below
\(3.9\cdot10^{-15}\), and the real rank of the physical response is
exactly \(2n-4\) in all twelve models.  That machine-precision
pattern is the next proof target, not a promoted theorem.

The dataset SHA-256 is
`9d3cadeb5deeb464b707c5fc8171f80b1bbe5296360d265a0b96cacd7fe46f02`.
