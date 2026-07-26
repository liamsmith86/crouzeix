# The fixed-base retightening leaves one compact nonlinear defect

## 1. Result (L299, 2026-07-26)

Let an analytic raw contraction metric and defect frame satisfy

\[
M(c)-A(c)^*M(c)A(c)=D(c)D(c)^*,\qquad
A(0)=S,\quad D(0)=V.                              \tag{1}
\]

Let \(R(c),F(c)\) be L298's weighted balanced metric/frame
correction.  At the fixed base colligation,

\[
R-S^*RS=VF^*+FV^*.                                \tag{2}
\]

Put

\[
\Delta=A-S,\qquad H=D-V.                          \tag{3}
\]

If one inserts L298 directly into the moving pair, the exact
remaining Stein defect is

\[
\begin{aligned}
{\cal N}(c)
={}&(M+R)-A^*(M+R)A-(D+F)(D+F)^*\\
={}&-\Delta^*RS-S^*R\Delta-\Delta^*R\Delta\\
 &-HF^*-FH^*-FF^*.                               \tag{4}
\end{aligned}
\]

This is the entire nonlinear debt.  There are no hidden inverse,
Schur, or endpoint terms in (4).  Since

\[
\Delta,H=O(c),\qquad R,F=O(c^2),
\]

one has

\[
\boxed{{\cal N}(c)=O(c^3).}                       \tag{5}
\]

Thus L298 correctly leaves every first completely delayed even face
unchanged.  However, (4) is **not** automatically a prior-flag factor.
There is an exact rational rank-one example in which its cubic
endpoint response has nonzero compression to \(\ker B_1^*\).
Consequently:

1. L298's fixed-base weighted sum must not be promoted directly to a
   moving nonlinear metric;
2. the effective route still needs one active mixed odd homology
   correction; and
3. that correction should be sought using L280/L296 or L289/L292,
   not by another complete-delay calculation.

This is a precise stop on the no-extra-correction shortcut, not a
failure of L298 or of the lower-tight route.

## 2. Exact defect identity

Subtract the raw equality (1) from the proposed corrected equality.
The residual is

\[
R-A^*RA-DF^*-FD^*-FF^*.                           \tag{6}
\]

Using \(A=S+\Delta\) gives

\[
A^*RA
=S^*RS+\Delta^*RS+S^*R\Delta+\Delta^*R\Delta.     \tag{7}
\]

Using \(D=V+H\) and then (2) in (6)--(7) cancels the two
fixed-base linear terms \(VF^*+FV^*\), leaving exactly (4).
The order statement (5) follows immediately.

This derivation is all-series and arbitrary-dimensional.  It replaces
the proposed next isolated-grade expansion by one six-term formula.

## 3. Cubic coefficient

Write

\[
\begin{aligned}
A(c)&=S+cA_1+O(c^2),&
D(c)&=V+cD_1+O(c^2),\\
R(c)&=c^2R_1+O(c^4),&
F(c)&=c^2F_1+O(c^4).
\end{aligned}                                      \tag{8}
\]

Then

\[
\boxed{
[c^3]{\cal N}
=-\{A_1^*R_1S+S^*R_1A_1+D_1F_1^*+F_1D_1^*\}.}   \tag{9}
\]

For L227's boundary repair,

\[
\begin{aligned}
A_1&=(I+WW^*)S^*(I+VV^*)-S^3,\\
H_1&=-(A_1^*S+S^*A_1),\\
D_1&=H_1V-\frac12V(V^*H_1V).
\end{aligned}                                      \tag{10}
\]

Since \(SV=0\), one has \(V^*H_1V=0\).  Equations (9)--(10)
and L298's \(V^*F_1=-B_1^*B_1/2\) then give

\[
\boxed{V^*[c^3]{\cal N}V=0.}                     \tag{11}
\]

So the first nonlinear debt is already lower-tight.  Its obstruction
is entirely at the upper endpoint.

## 4. Exact rational partial-flag obstruction

Take state dimension six and copy dimension two.  Let

\[
V=[e_1\ e_2],\qquad W=[e_5\ e_6],
\]

let the initial domain frame be \([e_3\ e_4\ e_5\ e_6]\), and let
the final range frame be \([e_1\ e_2\ e_3\ e_4]\).  Between them use
the rational orthogonal bridge

\[
U=
\begin{bmatrix}
3/5&-4/5&0&0\\
4/13&3/13&-12/13&0\\
48/65&36/65&5/13&0\\
0&0&0&1
\end{bmatrix}.                                    \tag{12}
\]

The resulting partial isometry \(S\) satisfies

\[
\det(\lambda I-S)=\lambda^5(\lambda-48/65),
\]

so it is spectrally stable, and

\[
B_1=W^*S^*V=
\begin{bmatrix}
0&-12/13\\
0&0
\end{bmatrix},\qquad
\ker B_1^*=\operatorname {span}\{e_2\}.            \tag{13}
\]

Construct \(R_1,F_1\) from L298 and \({\cal N}_3\) from (9).
Exact rational solution of the Stein equation gives

\[
\begin{aligned}
Q_3
&:=W^*{\cal G}_S({\cal N}_3)W\\
&=
\begin{bmatrix}
-\frac{2985984}{54865681}&
-\frac{139470336}{3566269265}\\
-\frac{139470336}{3566269265}&
\frac{2985984}{54865681}
\end{bmatrix}.                                    \tag{14}
\end{aligned}
\]

It is trace zero and indefinite, but on the surviving first flag

\[
\boxed{
e_2^*Q_3e_2
=\frac{2985984}{54865681}>0.}                     \tag{15}
\]

Thus \(Q_3\) is not a prior-\(B_1\) endpoint factor: such a factor
would have zero compression to \(\ker B_1^*\).  The cubic response
must be actively repaired.

Equation (15) does not rule out that repair.  Indeed the trace-zero
and indefinite structure is consistent with a homogeneous endpoint
response.  What remains is to derive a bounded L296 polarization
column (or L292 valuation certificate) for (9).

## 5. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_moving_retightening_defect.py \
  --output \
  experiments/repeated_crabb_moving_retightening_defect_s70226.jsonl
```

The checker reconstructs (12)--(15) with exact SymPy rational
arithmetic.  Nine additional noncommuting rank-chain cases in
multiplicities two through four verify the zero lower corner,
trace-zero indefinite upper response, and nonzero surviving-flag
compression while approaching repeated monomial strata.  The tracked
dataset has SHA-256

```text
caea7e9741a255ba6319e573050be51cd951889f8755bc286156885cae9622a9
```
