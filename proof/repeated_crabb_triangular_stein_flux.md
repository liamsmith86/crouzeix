# The radial volume target is one triangular Stein boundary flux

## 1. Result (L265, 2026-07-25)

Retain L252's pure partial isometry

\[
I-S^*S=E,\qquad I-SS^*=F,
\]

and put

\[
E_j=(S^*)^jES^j,\qquad Q_j=(S^*)^jS^j.
\]

For \(k\geq1\), define the finite triangular right-orbit potential

\[
\boxed{
X_k=\sum_{j=0}^{k}(k+1-j)E_j.}                    \tag{1}
\]

Then L252's three-term radial operator is exactly its Stein
divergence:

\[
\boxed{
Q_{k+2}-(k+2)Q_1+(k+1)I
=X_k-S^*X_kS.}                                    \tag{2}
\]

Consequently

\[
\boxed{
\operatorname {tr}
\{Q_{k+2}-(k+2)Q_1+(k+1)I\}
=\operatorname {tr}(FX_k)
=\sum_{j=1}^{k}(k+1-j)\|B_j\|_F^2,}              \tag{3}
\]

where \(B_j=W^*(S^*)^jV\) for defect frames
\(E=VV^*\), \(F=WW^*\).

Under complete delay,

\[
B_1=\cdots=B_{k-1}=0,
\]

only the last orbit in (1) reaches the final defect, and

\[
\boxed{
\operatorname {tr}(FX_k)=\|B_k\|_F^2.}            \tag{4}
\]

Thus A194's remaining cyclic congruence is equivalently the discrete
flux law

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathfrak D_{\rm ret}
=4\operatorname {tr}(X_k-S^*X_kS).}              \tag{5}
\]

L265 does not prove (5).  It supplies a more structural target for
the L242--L245/L258 calculation: after L264 restores the complete
right-half-line metric, the full physical numerator should reduce,
modulo trace-null terms, to four times the finite divergence of the
triangular right-orbit potential (1).

This is compatible with the indefinite active residual.  Equation
(5) is a trace conservation law, not a positive operator-Gram
factorization.

## 2. Operator telescope

Since \(S^*S=I-E\),

\[
\begin{aligned}
Q_{j+1}
&=(S^*)^jS^*SS^j\\
&=(S^*)^j(I-E)S^j\\
&=Q_j-E_j.
\end{aligned}                                    \tag{6}
\]

Therefore

\[
Q_{k+2}-I=-\sum_{j=0}^{k+1}E_j.                  \tag{7}
\]

Also \(Q_1=I-E\), so

\[
\begin{aligned}
Q_{k+2}-(k+2)Q_1+(k+1)I
&=Q_{k+2}-I+(k+2)E\\
&=(k+1)E_0-\sum_{j=1}^{k+1}E_j.                  \tag{8}
\end{aligned}
\]

On the other hand,

\[
S^*X_kS
=\sum_{j=1}^{k+1}(k+2-j)E_j.
\]

Subtracting this from (1) gives exactly (8), proving (2).

## 3. Boundary trace

Trace cyclicity and \(SS^*=I-F\) give

\[
\operatorname {tr}(X_k-S^*X_kS)
=\operatorname {tr}(FX_k).                       \tag{9}
\]

For every \(j\geq1\),

\[
\begin{aligned}
\operatorname {tr}(FE_j)
&=\operatorname {tr}
\{W^*(S^*)^jVV^*S^jW\}\\
&=\|V^*S^jW\|_F^2\\
&=\|B_j\|_F^2,
\end{aligned}                                    \tag{10}
\]

while \(FE_0=0\).  Insert (1) into (9) and use (10) to prove (3).
Complete delay leaves only \(j=k\), whose triangular coefficient is
one, proving (4).

## 4. Interface with the live proof

L253--L261 express the same scalar as a Toeplitz-window leakage and
then as the norm of the complete future leakage row.  L265 adds the
state-side conservation form:

\[
\text{Hardy leakage energy}
=\text{right-orbit Stein flux}.
\]

This is the natural form after L264.  The right-half-line metric is
now complete through the active degree, and L244/L262 make its bulk
coisometric.  A proof of A194 should therefore seek an exact
coefficientwise divergence identity for the full L258 numerator:

\[
[c^{2k}]\mathfrak D_{\rm ret}
=4(X_k-S^*X_kS)+\text{trace-null terms}.
\]

The universal relative jet in L242 limits that verification to the
three coefficients \(D_0,D_1,D_2\); L245/L251 must supply the omitted
chain flux before the divergence is taken.
