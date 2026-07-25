# The first new Schur face ignores analytic motion of the defect graph

## 1. Result (L263, 2026-07-25)

Let \(J:\mathbb C^m\to{\cal H}\) be an isometry, put

\[
F=JJ^*,\qquad P=I-F,
\]

and let \(d(c):\mathbb C^m\to{\cal H}\) be a formal analytic column
such that

\[
d(0)=J,\qquad a(c):=J^*d(c)\ \hbox{is invertible}.
\]

Suppose a Hermitian formal series has the form

\[
H(c)=d(c)d(c)^*+c^s\Delta(c),\qquad s\geq1.       \tag{1}
\]

Let \(K(c)\) be its Schur residual away from the fixed \(F\)-row:

\[
K=PHP-PHF(FHF)^{-1}FHP.                           \tag{2}
\]

Define the graph quotient

\[
\boxed{
T(c)=P-Pd(c)a(c)^{-1}J^*.}                        \tag{3}
\]

Then

\[
\boxed{
K(c)=c^sT(c)\Delta(c)T(c)^*+O(c^{2s}).}           \tag{4}
\]

In particular \(T(0)=P\), so the first new face is

\[
\boxed{
[c^s]K=P\Delta(0)P,\qquad
\operatorname {tr}[c^s]K
=\operatorname {tr}\{P\Delta(0)P\}.}             \tag{5}
\]

Thus analytic movement of the rank-\(m\) background defect graph
cannot insert a hidden theta weight into its first nonzero Schur
face.

There is an equivalent Wold-coordinate statement.  If L262's formal
unitary \({\cal W}_c^*=I+O(c)\) is used to write

\[
\Delta(c)={\cal W}_c^*\Lambda(c){\cal W}_c,
\]

then

\[
\boxed{
\operatorname {tr}[c^s]K
=\operatorname {tr}\{P\Lambda(0)P\}.}            \tag{6}
\]

This closes the fixed-port qualification left implicit in L262:
the Wold unitary need not preserve \(F\) exactly, because its motion
is one order too late to affect a first new face.

L263 does **not** identify the physical perturbation
\(\Lambda(0)\).  For A194 one must still prove, using L243,
L245/L251, and L258, that after the word-free graph update and the
full two-orientation cancellation its retained trace is

\[
4\operatorname {tr}(P_k{\cal L}_B^2P_k).
\]

The point of L263 is that this remaining problem is now purely the
physical channel numerator.  Neither the moving defect column, the
fixed Schur port, nor the half-line metric can change its leading
scalar coefficient.

## 2. Differential Schur formula

Write

\[
d=\begin{bmatrix}a\\b\end{bmatrix},
\qquad
\Delta=
\begin{bmatrix}
\Delta_{FF}&\Delta_{FP}\\
\Delta_{PF}&\Delta_{PP}
\end{bmatrix}
\]

relative to \(F{\cal H}\oplus P{\cal H}\), using \(J\) to identify
\(F{\cal H}\) with \(\mathbb C^m\).  The blocks of the background
Gram \(dd^*\) are

\[
aa^*,\qquad ab^*,\qquad ba^*,\qquad bb^*.
\]

Put \(\varepsilon=c^s\).  The standard inverse expansion gives

\[
\begin{aligned}
(aa^*+\varepsilon\Delta_{FF})^{-1}
={}&(aa^*)^{-1}\\
&-\varepsilon(aa^*)^{-1}\Delta_{FF}(aa^*)^{-1}
+O(\varepsilon^2).
\end{aligned}
\]

Substitution in (2) cancels the entire background Gram and leaves

\[
\begin{aligned}
\varepsilon^{-1}K
={}&\Delta_{PP}
-ba^{-1}\Delta_{FP}
-\Delta_{PF}(a^{-1})^*b^*\\
&+ba^{-1}\Delta_{FF}(a^{-1})^*b^*
+O(\varepsilon).
\end{aligned}                                    \tag{7}
\]

The right side of (7) is exactly \(T\Delta T^*\), proving (4).
Because \(d(0)=J\), one has \(b(0)=0\), \(a(0)=I_m\), and hence
\(T(0)=P\).  Extracting the coefficient of \(c^s\) proves (5).

## 3. Wold transport

L262 gives

\[
{\cal W}_c^*=I+O(c),\qquad {\cal W}_c=I+O(c).
\]

Consequently

\[
T(c){\cal W}_c^*
=P+O(c).
\]

Insert
\(\Delta={\cal W}_c^*\Lambda{\cal W}_c\) into (4).  At degree
\(s\), every nonconstant coefficient of either exterior factor
would add at least one further power of \(c\).  Therefore

\[
[c^s]K=P\Lambda(0)P.
\]

Taking the finite stabilized trace proves (6).  No bounded
boundary functional calculus or rank-continuity argument is used.

## 4. Exact scope at the live gate

For L240's background,

\[
d(0)=J_0,\qquad J_0^*d(c)=\vartheta_3(c^2)^{-1}I_m,
\]

so every hypothesis above holds.  L258 already guarantees that the
physical final-row calculation is a Schur quotient of this type.
Accordingly the remaining proof may work entirely in the constant
Hardy coordinates when extracting the first genuinely new face.

What remains unproved is not a metric or graph issue:

1. isolate the physical \(\Delta\) after summing both ellipse
   orientations;
2. use L243 to retain only the zero/one inverse-kernel sectors;
3. show that its retained trace is four times the complete L259
   leakage-row norm; and
4. invoke L261 only after that complete row has been formed.

An operator Gram claim is neither required nor expected: finite
audits show that the active residual matrix can be indefinite.  L263
is only a first-face trace/Schur transport identity.
