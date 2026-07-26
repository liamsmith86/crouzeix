# The complete rooted quotient is one copy-space energy defect

## 1. Result (L306, 2026-07-26)

Retain L305's balanced stable partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\]

its transfer channel \(\Phi\), and its perpendicular response map
\({\cal M}_S\).  Define the dual initial-copy closure

\[
\boxed{
\Psi_S(X)
=V^*{\cal G}_S^\vee(X)V,\qquad
{\cal G}_S^\vee(X)=\sum_{n\ge0}S^nX(S^*)^n.}      \tag{1}
\]

Let \(A(c),D(c)\) be the moving raw operator and defect frame in
L299, and let \(X(c)=X(c)^*\), \(C(c)\) be a prepared metric/frame
correction satisfying the fixed-base equation

\[
\boxed{X-S^*XS=VC^*+CV^*.}                       \tag{2}
\]

Put

\[
\Delta=A-S,\qquad H=D-V.                         \tag{3}
\]

L299's complete correction defect is

\[
\begin{aligned}
{\cal N}={}&-\Delta^*XS-S^*X\Delta-\Delta^*X\Delta\\
&-HC^*-CH^*-CC^*.
\end{aligned}                                     \tag{4}
\]

Its complete upper Stein endpoint has the exact quotient form

\[
\boxed{
W^*{\cal G}_S({\cal N})W
=\Phi({\mathfrak q})+{\cal M}_S(C_{\rm rsp}),}    \tag{5}
\]

where the retained copy-space energy matrix is

\[
\boxed{
\begin{aligned}
{\mathfrak q}={}&
-\Psi_S\!\left(
 \Delta^*XS+S^*X\Delta+\Delta^*X\Delta
 \right)\\
&-\left(C^*H+H^*C+C^*C\right).
\end{aligned}}                                    \tag{6}
\]

For every fixed repeated-block neighbourhood and finite prepared
jet, \(C_{\rm rsp}\) has a rank-stable bound linear in the active
bridge amplitudes.  With the physical bridge-word presentations from
L298 and L305, its state lift can be chosen in the same cumulative
bridge module.

Thus the all-grade retained-root problem is no longer a growing
state-word or mixed-graph calculation.  It is the sign and
factorization of the single finite copy matrix (6).

There is an equivalent unnormalized form.  Before subtracting the
fixed-base equation, put

\[
\begin{aligned}
{\mathfrak q}_{\rm full}
={}&\Psi_S(X-A^*XA)\\
&-\{C^*D+D^*C+C^*C\},\\
{\mathfrak q}_{\rm base}
={}&\Psi_S(X-S^*XS)\\
&-\{C^*V+V^*C\}.
\end{aligned}                                    \tag{7}
\]

Then

\[
\boxed{{\mathfrak q}
={\mathfrak q}_{\rm full}-{\mathfrak q}_{\rm base},} \tag{8}
\]

and expansion of \(A=S+\Delta\), \(D=V+H\) gives (6) directly.

L306 is a structural all-series reduction.  It does not prove that
(6) is favorable.  The remaining partial-retightening gate is now
precise:

1. its odd prepared coefficients must be response-null; and
2. its even quotient, after the compulsory lower neutralization,
   must be a two-ended transfer Gram dominated by an earlier retained
   margin.

L300 and L303--L304 prove those statements only at cubic and quartic
order.  No isolated higher coefficient is promoted here.

## 2. Canonical copy closure of a bridge-ideal state

Let \(Z=Z^*\) be a finite bridge-ideal word polynomial.  L305 gives a
rooted retained channel plus an ideal-valued response.  The retained
channel may be replaced, without a rank-changing loss, by the
coordinate-free closure (1):

\[
\boxed{
W^*{\cal G}_S(Z)W
=\Phi(\Psi_S(Z))+{\cal M}_S(K_Z).}                \tag{9}
\]

The column \(K_Z\) is perpendicular and obeys

\[
\|K_Z\|_F
\le c_Z(S)\sum_j\|B_j\|_F,                       \tag{10}
\]

where the finite constant depends on the displayed word
presentation and on a stable power moment of \(S\), but not on the
ranks of the \(B_j\).

To see the mechanism, let

\[
H_Y-SH_YS^*=WYW^*,\qquad
A_Y=V^*H_YV,\qquad
Z_Y=(I-E)H_YV.                                   \tag{11}
\]

Stability and \(I-S^*S=E\) give the Wold resolution

\[
\boxed{I=\sum_{n\ge0}(S^*)^nES^n.}               \tag{12}
\]

Insert (12) in \(\operatorname{tr}(H_YZ)\).  Move \(H_Y\) through the
power words and use

\[
[H_Y,S]=-SZ_YV^*,\qquad
[H_Y,S^*]=VZ_Y^*S^*.                             \tag{13}
\]

The part with no commutator is

\[
\operatorname{tr}\{A_Y\Psi_S(Z)\}
=\operatorname{tr}\{Y\Phi(\Psi_S(Z))\}.           \tag{14}
\]

Every error term contains one \(Z_Y\).  Because every displayed
summand of \(Z\) contains a physical bridge
\(ES^jF=VB_j^*W^*\) or its adjoint, the same term also contains
\(B_j\).  Therefore

\[
\left|
\left\langle Y,\,
W^*{\cal G}_S(Z)W-\Phi(\Psi_S(Z))
\right\rangle
\right|
\le2c_Z(S)
\left(\sum_j\|B_j\|_F\right)\|Z_Y\|_F.           \tag{15}
\]

Since \({\cal M}_S^*Y=2Z_Y\), Hilbert quotient duality gives (9)--(10).
Expanding the commutators in (13) and using L305 termwise chooses the
column with its state lift in the original bridge module.

For a norm-convergent stable Stein sum, the same proof applies by
completion.  One possible local bound uses

\[
\Lambda_S=\sum_{n\ge0}(n+1)\|S^n\|<\infty.       \tag{16}
\]

This is uniformly finite after shrinking a fixed repeated-block
neighbourhood.  No Markov spectral-gap inverse occurs.

## 3. Covariant frame products close in copy space

The response cokernel is characterized by \(Z_Y=0\).  In that case
(13) says \(H_Y\) commutes with \(S,S^*\).  Equation (11) also gives
the transfer intertwinings

\[
\boxed{
YB_j=B_jA_Y,\qquad
A_YB_j^*=B_j^*Y.}                                \tag{17}
\]

The physical moving frame and every L298/L305 prepared column are
covariant:

\[
\boxed{H_YD=DA_Y,\qquad H_YC=CA_Y
\quad(Z_Y=0).}                                   \tag{18}
\]

For the raw frame this follows coefficientwise from its Hermitian
factor gauge.  Its constant coefficient is \(V\), so (18) starts
with \(H_YV=VA_Y\).  The raw operator, metric, and slack coefficients
are built equivariantly from \(S,S^*\) and their stable Stein
closures, so they commute with \(H_Y\).  If (18) holds for all earlier
factor coefficients, the next Hermitian factor remainder \(R_n\)
therefore commutes with \(H_Y\).  Also \(E\), \(Q=I-E\), and \(H_Y\)
commute, while
\(A_Y(V^*R_nV)=(V^*R_nV)A_Y\).  Hence both terms in the gauge formula
\[
Q R_nV+\frac12V(V^*R_nV)
\]
intertwine \(H_Y\) with \(A_Y\), proving the induction.  For
L298/L305 columns, (18) follows directly from (17) and their displayed
word formulas.

Consequently, for any two covariant columns \(D,C\),

\[
\begin{aligned}
\operatorname{tr}(H_YDC^*)
 &=\operatorname{tr}(A_YC^*D),\\
\operatorname{tr}(H_YCD^*)
 &=\operatorname{tr}(A_YD^*C),\\
\operatorname{tr}(H_YCC^*)
 &=\operatorname{tr}(A_YC^*C).                   \tag{19}
\end{aligned}
\]

Equations (9), (14), and (19) show that the endpoint of

\[
X-A^*XA-DC^*-CD^*-CC^*
\]

has retained representative
\(\Phi({\mathfrak q}_{\rm full})\).  Subtracting the fixed-base
representative gives (5)--(8).

Away from the exact cokernel, the failures of (18) are again finite
word commutators controlled by \(Z_Y\), and every correction factor
contains an active bridge.  L305 therefore supplies the same
gap-free \(O(\|B_j\|)\) response bound rather than merely pointwise
cokernel orthogonality.

## 4. Why (6) is the correct next object

Formula (6) separates three effects without expanding a grade:

1. the two linear moving-operator crosses;
2. the quadratic moving-operator cross; and
3. the complete frame energy
   \(C^*H+H^*C+C^*C\).

The negative square \(-C^*C\) is already favorable.  The only
possible indefinite content is the paired motion of the raw operator
and frame.  At cubic order L300 proves that paired motion is
response-null.  At quartic order L302--L304 plus L303's lower
neutralization reduce it to a polynomial in \(B_1B_1^*\) supported
on the earlier active range.

An all-grade proof should now work directly with (6), using the
theta/ODE covariance of \(A,D\) and the transfer intertwinings (17).
It should not:

1. enumerate another state-word coefficient;
2. assume the two ellipse orientations cancel separately (L260);
3. discard the favorable square \(-C^*C\); or
4. replace the bounded response by a rank-changing pseudoinverse.

## 5. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_copy_energy_quotient.py \
  --output \
  experiments/repeated_crabb_copy_energy_quotient_s70226.jsonl
```

The checker constructs the actual degree-four canonical moving
operator and Hermitian-gauge defect frame.  It tests simultaneous
retightening grades \(1,2,3\) on unstructured colligations, grades
\(1,2\) through three rank-collapse scales, and the nonvacuous grades
\(2,3\) after a complete first delay.  It independently verifies the
fixed-base equation, the six-term defect, both forms (6)--(8), and
that the residual endpoint lies in the perpendicular response range.
All 16 tracked records pass.  A further 100 untracked mixed-grade
stress cases passed with maximum quotient reconstruction error
\(1.20\times10^{-18}\) and response synthesis error
\(1.74\times10^{-18}\).  A separate multiplicity-three reducing
colligation with a nonscalar response cokernel tested (18) directly:
the raw-frame and prepared-column covariance residuals were
\(1.67\times10^{-14}\) and \(6.36\times10^{-18}\), respectively.

The tracked dataset regenerates byte-identically with SHA-256

```text
cc538f27b3c5395a90def31a062bad12bf2449589cecf0823b139764f1275d65
```

These floating audits test the complete nonlinear assembly and rank
collapse; the proof is (9)--(19).
