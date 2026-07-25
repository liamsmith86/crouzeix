# The active right-half-line metric may be restored for free

## 1. Result (L264, 2026-07-25)

Let \(S\) be a partial isometry with

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad P=I-F.
\]

Let \(A(c)=S+O(c)\).  Suppose two formal metrics \(R\) and
\(\widehat R\) both isolate the final row through degree \(d\):

\[
R=F\oplus PRP+O(c^{d+1}),\qquad
\widehat R=F\oplus P\widehat RP+O(c^{d+1}),       \tag{1}
\]

agree below degree \(d\), and satisfy

\[
\widehat R-R=c^dX+O(c^{d+1}),\qquad X=PXP.        \tag{2}
\]

For either metric define L250's balanced initial-defect mass

\[
\mu_R
=-m+\operatorname {tr}
\left\{
(R-A^*FA)^{-1}(R-A^*RA)
\right\}.                                       \tag{3}
\]

Then the exact first response is

\[
\boxed{
[c^d](\mu_{\widehat R}-\mu_R)
=\operatorname {tr}\{(F-E)X\}
=-\operatorname {tr}(EX).}                       \tag{4}
\]

Consequently the active scalar volume is unchanged whenever

\[
\boxed{\operatorname {tr}(EX)=0.}                \tag{5}
\]

The convenient sufficient condition \(EX=XE=0\) will be used below.

Equation (4) is also an important guardrail: an arbitrary retained
metric coefficient is **not** free.  Its initial-defect compression
changes the active trace.  This agrees with L247 rather than
contradicting it.

## 2. Port-isolated right-half-line gauge

Write L219's metric as

\[
R=I+R_{\rm right}+R_{\rm left},
\]

where

\[
\begin{aligned}
R_{\rm right}
&=-\sum_{j\geq1}\frac{q^j}{1+q^j}(S^*)^jES^j,\\
R_{\rm left}
&=\sum_{j\geq1}q^jS^jF(S^*)^j,\qquad q=c^2.
\end{aligned}                                    \tag{6}
\]

Assume complete delay

\[
B_1=\cdots=B_{k-1}=0
\]

and put

\[
X_k=[q^k]R,\qquad X_k^R=[q^k]R_{\rm right}.
\]

Starting from L250's whole-coefficient deletion

\[
R^\circ=R-q^kX_k,
\]

define

\[
\boxed{
R^\triangleright
=R^\circ+q^kPX_k^RP.}                             \tag{7}
\]

Then

\[
\boxed{
[c^{2k}]\mu_{R^\triangleright}
=[c^{2k}]\mu_{R^\circ},}                          \tag{8}
\]

and the final row/cross blocks remain isolated.  Moreover the
complete active coefficient of the **right-half-line metric** is
restored on the retained block.  The active left-orbit coefficient
remains deleted.

Thus A194 may be evaluated in the port-isolated right-half-line
gauge (7).  L240/L244's \(D_R\) is present through the active degree,
while the transfer-dependent left-boundary face that carries the
initial-defect mass has not been silently restored.

L264 does not evaluate the remaining physical numerator.

## 3. Proof of the response formula

Put

\[
N_R=R-A^*FA,\qquad H_R=R-A^*RA.
\]

Because \(FS=0\),

\[
N_R(0)=I,\qquad H_R(0)=E.                         \tag{9}
\]

The two inputs agree below degree \(d\), so their difference at
degree \(d\) is the linear response at \((A,R)=(S,I)\).  From (2),

\[
\delta N=X,\qquad
\delta H=X-S^*XS.
\]

Using \(\delta(N^{-1})=-X\),

\[
\begin{aligned}
[c^d](\mu_{\widehat R}-\mu_R)
&=\operatorname {tr}\{-XE+X-S^*XS\}\\
&=\operatorname {tr}\{X(I-E)-XSS^*\}\\
&=\operatorname {tr}\{X(S^*S-SS^*)\}\\
&=\operatorname {tr}\{(F-E)X\}.                  \tag{10}
\end{aligned}
\]

Since \(FX=0\), equation (10) is (4).  This proves (5).

## 4. The restored right coefficient misses the initial defect

The coefficient \(X_k^R\) is a linear combination of

\[
E_j=(S^*)^jES^j,\qquad j\mid k,\quad j\geq1.
\]

Purity gives the orthogonal right wandering resolution

\[
I=\sum_{j\geq0}(S^*)^jES^j.
\]

Hence

\[
EE_j=E_jE=0\qquad(j\geq1).                        \tag{11}
\]

Because \(EF=0\), one also has \(EP=PE=E\).  Therefore

\[
E(PX_k^RP)=(PX_k^RP)E=0.                         \tag{12}
\]

Apply (5) with \(d=2k\) and \(X=PX_k^RP\) to obtain (8).
L250 already proves that \(R^\circ\) isolates the final row; adding
a \(PXP\) coefficient preserves that isolation.

## 5. Guardrail

The left coefficient

\[
[q^k]R_{\rm left}=S^kF(S^*)^k
\]

is not gauge-free.  Its initial-defect compression has trace

\[
\operatorname {tr}\{E S^kF(S^*)^k\}
=\|B_k\|_F^2.
\]

Restoring it would change (3) at the active degree by
\(-\|B_k\|_F^2\).  Likewise the deleted final corner/cross row cannot
be restored.  The valid gauge is exactly (7): restore the universal
right-half-line metric, not the complete physical retained metric.
