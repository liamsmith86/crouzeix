# Cyclic support through index \(k+2\) is enough for delay covariance

## 1. Result (L276, 2026-07-25)

Retain L250/L258's completely delayed, active-metric-deleted
closed-return defect

\[
\mathfrak D_k(c)=I-R_{P,k}Z_{{\rm ret},k}
\]

and the independently balanced deflated tail
\(\mathfrak D_{k-1}^{\rm tail}\), with \(k\geq2\).  Put

\[
\Delta_k=[c^{2k}]\mathfrak D_k
-[c^{2k-2}]\mathfrak D_{k-1}^{\rm tail}.          \tag{1}
\]

The following two support statements are sufficient for the complete
one-delay scalar recursion.

1. Every odd lower face vanishes, and for every \(h<k\) the even face
   \([c^{2h}]\mathfrak D_k\) is a two-sided radial polynomial
   \[
   \sum_{j=0}^{J_h}a_{h,j}Q_j
   +\sum_{j=1}^{J_h}b_{h,j}R_j,
   \qquad
   Q_j=(S^*)^jS^j,\quad R_j=S^j(S^*)^j,           \tag{2}
   \]
   for some finite \(J_h\).  No uniform sharp index bound is needed
   for a fixed formal coefficient.
2. In the universal partial-isometry/delay quotient, modulo trace
   cyclicity (and without Cayley--Hamilton), the active associated
   difference (1) is a radial polynomial of index at most \(k+2\):
   \[
   \Delta_k\ \overset{\rm tr}{\sim}\
   P_k(Q)=\sum_{j=0}^{k+2}p_{k,j}Q_j.             \tag{3}
   \]

Then every lower face in (2) is the **zero operator**, and

\[
\boxed{\operatorname {tr}\Delta_k=0.}             \tag{4}
\]

Consequently

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathfrak D_k
=4\|B_k\|_F^2,}                                  \tag{5}
\]

and all earlier coefficients of \(\mathfrak D_k\) vanish.

This is a strictly weaker stop condition than L271's exact
full-support formula.  Words outside L275's proved fan are harmless
if they cyclically radialize without exceeding index \(k+2\).
Therefore the campaign need not prove literal absence of every extra
operator word in order to close the scalar volume identity.

L276 proves the implication (2)--(5).  It does **not** itself prove
the two support premises.  L277 subsequently proves every lower odd
face vanishes in arbitrary grade, and L278 proves every lower even
face vanishes as well.  Exact word arithmetic verifies the remaining
active-cyclic statement through the existing finite range, but its
arbitrary-grade theorem is now the sole gate on this route.

## 2. Exact axes kill every lower radial face

Fix \(h<k\).  Take a unilateral shift of length
\(N>\max\{2J_h,k+1\}\), or an \(m\)-fold inflation of it, with its first
transfer cell beyond grade \(k\).  L117/L241's exact axis metric has
zero normalized defect through degree \(2k\).  The edge deletion is
at degree \(2k\), so it does not alter the lower coefficient
\(2h\).  Hence the polynomial (2) vanishes on every such shift.

On a length-\(N\) shift, the diagonal projections

\[
I,\ Q_1,\ldots,Q_{J_h},\
R_1,\ldots,R_{J_h}                                \tag{6}
\]

are linearly independent when \(N>2J_h\).  Indeed, successive
diagonal differences from the right endpoint isolate the \(Q_j\)
coefficients, successive differences from the left endpoint isolate
the \(R_j\) coefficients, and an interior entry isolates the identity
coefficient.  Thus every coefficient in (2) is zero.

Together with the assumed odd-face cancellation, this proves

\[
[c^d]\mathfrak D_k=0\qquad(d<2k),                 \tag{7}
\]

The operator conclusion in (7), rather than mere trace vanishing, is
important: it makes the active log-determinant coefficient linear,
as required by L258.

## 3. Three scalar tests force the active trace

Write

\[
p_k(z)=\sum_{j=0}^{k+2}p_{k,j}z^j.
\]

L252 gives, under
\(B_1=\cdots=B_{k-1}=0\),

\[
\operatorname {tr}Q_j=
\begin{cases}
n-jm,&0\leq j\leq k+1,\\
n-(k+2)m+\|B_k\|_F^2,&j=k+2.
\end{cases}                                      \tag{8}
\]

Therefore (3) has the universal trace form

\[
\boxed{
\operatorname {tr}\Delta_k
=p_k(1)n-p_k'(1)m
+p_{k,k+2}\|B_k\|_F^2.}                          \tag{9}
\]

The three scalar coefficients in (9) are forced by exact cases.

### Long axes

Use arbitrarily long inflated monomial axes whose first transfer cell
lies beyond \(k\).  Both the full and deflated associated faces are
zero by L241, and L247's deletion response is zero because \(B_k=0\).
Thus (9) vanishes for arbitrarily many values of \(n/m\), forcing

\[
p_k(1)=0,\qquad p_k'(1)=0.                       \tag{10}
\]

### Active monomial

Now use the length-\(k\) monomial channel
\(B(z)=z^kI_m\).  L241 gives total delayed Schur trace \(2m\) in
every grade, while L247's metric contribution is \(-2m\); hence the
edge-deleted face is \(4m\).  Its one-layer deflated tail is the
length-\((k-1)\) monomial channel and has the same edge-deleted face
\(4m\) (with L256 supplying the final grade-one case).  Therefore

\[
\operatorname {tr}\Delta_k=0.
\]

Since \(\|B_k\|_F^2=m\), equations (9)--(10) force

\[
p_{k,k+2}=0.                                     \tag{11}
\]

Substitution of (10)--(11) into (9) proves (4) for every completely
delayed partial isometry, with arbitrary defect multiplicity and
arbitrary noncommuting future transfer coefficients.

## 4. Iteration and the volume coefficient

Equation (4) is the associated recursion

\[
\operatorname {tr}[c^{2k}]\mathfrak D_k
=\operatorname {tr}[c^{2k-2}]
\mathfrak D_{k-1}^{\rm tail}.                    \tag{12}
\]

Iterate (12) through the clean delay flag.  The final tail has active
transfer \(B_k\) at relative grade one, and L256 proves

\[
\operatorname {tr}[c^2]\mathfrak D_1
=4\|B_k\|_F^2.
\]

Together with the operator lower vanishings (7), L258 gives

\[
[c^{2k}]\log\det(I+\mathfrak D_k)
=\operatorname {tr}[c^{2k}]\mathfrak D_k,
\]

which proves (5).

## 5. New live support target

After L277--L278, the remaining proof should no longer seek the stronger
statement

\[
[c^{2k}]\mathfrak U_k
=\text{radial}+\mathcal F_k
\]

unless that exact identity falls out for free.  It is enough to prove
the sole still-open statement:

1. every active associated word, including any word outside L275's
   fan, cyclically radializes with index at most \(k+2\).

L277--L278 prove every lower face vanishes.  L275 already proves
coefficient one on the complete fan, and every
fan word cyclically reduces to \(Q_{k+2}\).  L273 supplies the exact
metric convolution.  Thus only cyclic control of possible **extra**
active words remains; their individual coefficients are irrelevant.
