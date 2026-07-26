# Endpoint-null gauge transport breaks at the next even coefficient

## 1. Result (L287, 2026-07-25)

L286 makes L284's quartic and quintic delay normalizations
simultaneous through degree five.  Apply L285's recurrence once more,
include L284's individually normalized sextic column, and perform the
compulsory lower-endpoint elimination.

The resulting sixth metric forcing \(H_6^0\) has

\[
EH_6^0E=0,
\]

but in general its Stein response has a nonzero upper trace:

\[
\boxed{
\operatorname{tr}\!\left(
W^*\mathcal G_S(H_6^0)W
\right)\ne0.}                                     \tag{1}
\]

Every remaining lower-preserving factor correction is a
perpendicular column, up to a forcing-null skew parallel part.
L204 proves that every perpendicular endpoint response has trace
zero.  Therefore no such correction can make the sixth metric
coefficient globally endpoint-null.

This is an exact obstruction, not a failed finite search.  On the
rational scalar-defect partial isometry in Section 4,

\[
\boxed{
W^*\mathcal G_S(H_6^0)W
=\frac{151040}{177147}>0.}                        \tag{2}
\]

Consequently L284--L286's convenient convention

> every transported metric coefficient has zero compression at both
> defect endpoints

cannot be the arbitrary-grade induction.

L287 does **not** disprove bounded repeated-elliptic selection.
At even order the nonzero metric trace must be retained together with
the direct transfer Gram and the nonlinear endpoint Schur terms.
The correct recursion should cancel removable odd faces but control,
rather than erase, the physical even face.  This is consistent with
L230--L234's observed odd-cancellation/even-Gram pattern.

## 2. The transported sextic forcing

Let the changes already fixed by L286 be

\[
\begin{aligned}
\Delta D_4&=-N_4, &
\Delta M_4&=-Z_4,\\
\Delta D_5&=-N_5+K_5, &
\Delta M_5&=Z_5^{\rm tr}.
\end{aligned}                                     \tag{3}
\]

Use L285's triangular recurrence at order six and denote its raw
forcing by \(R_6\).  In factor-lift notation it contains

\[
\begin{aligned}
R_6={}&
\Delta D_4D_2^*+D_2\Delta D_4^*
+\Delta D_5D_1^*+D_1\Delta D_5^*\\
&+\sum_{\substack{a+b+d=6\\d\in\{4,5\}}}
T_a^*\Delta M_dT_b .                              \tag{4}
\end{aligned}
\]

There is no quadratic \(\Delta D_i\Delta D_j^*\) term yet because
the first changed factor degree is four.

Let \(\mathcal P_6=N_6V^*\) be L284's endpoint-null sextic gauge and
put

\[
H_6=R_6-\mathcal P_6-\mathcal P_6^*.              \tag{5}
\]

L285's parallel lower correction is

\[
\mathcal P_6^\parallel=-\frac12EH_6E.             \tag{6}
\]

Its forcing contribution is \(-EH_6E\), so

\[
\boxed{H_6^0=H_6-EH_6E.}                          \tag{7}
\]

Exact word reduction gives

\[
\begin{array}{c|c|c}
\text{polynomial}&\text{reduced words}&\ell^1\\ \hline
H_6^0&260&756\\
\mathcal P_6^\parallel&40&48.
\end{array}                                       \tag{8}
\]

Moreover,

\[
\boxed{
EH_6^0E=0,\qquad
\mathcal P_6^\parallel\in\mathcal I_3,}            \tag{9}
\]

where \(\mathcal I_3\) is the triple-delay ideal.  Thus the lower
parallel correction itself does not spoil L284's desired triple-delay
normalization.  The obstruction lies at the upper endpoint.

## 3. Why a nonzero trace is decisive

After (6), every further factor column that keeps the lower endpoint
zero can be written as

\[
C=VA+C_\perp,\qquad V^*C_\perp=0,
\]

where the Hermitian part of \(A\) is already fixed by the lower
equation.  A skew-Hermitian \(A\) contributes no Hermitian forcing.
The effective remaining freedom is therefore \(C_\perp\).

L204's endpoint map is

\[
\mathcal M_S(C_\perp)
=W^*\mathcal G_S(VC_\perp^*+C_\perp V^*)W.
\]

The identity copy matrix is always in its adjoint cokernel, so

\[
\boxed{\operatorname{tr}\mathcal M_S(C_\perp)=0.} \tag{10}
\]

For defect multiplicity one, (10) says
\(\mathcal M_S(C_\perp)=0\) identically.  Hence any scalar example
with (1) disproves an endpoint-null completion.

## 4. Exact rational counterexample

Take

\[
S=
\begin{bmatrix}
0&7/9&-4/9&-4/9\\
0&-4/9&1/9&-8/9\\
0&-4/9&-8/9&1/9\\
0&0&0&0
\end{bmatrix},
\qquad
V=e_1,\qquad W=e_4.                               \tag{11}
\]

Direct rational multiplication gives

\[
S^*S=I-VV^*,\qquad SS^*=I-WW^*,\qquad V^*W=0.
\]

Its characteristic polynomial is

\[
\lambda^2(3\lambda+2)^2/9,
\]

so it is pure and spectrally stable, with spectral radius \(2/3\).
Its first transfer coefficient is nonzero:

\[
B_1=W^*S^*V=V^*SW=-\frac49.                      \tag{12}
\]

Evaluate the 260-word rational polynomial \(H_6^0\) from (7) on
(11), and solve the \(4\times4\) Stein equation exactly over
\(\mathbb Q\):

\[
X-S^*XS=H_6^0(S,S^*).
\]

The exact residual is zero and

\[
W^*XW=\frac{151040}{177147},
\]

which proves (2).  No floating approximation or rank decision occurs
in this certificate.

## 5. Corrected arbitrary-grade target

Do not demand a globally endpoint-null \(Z_r\) at every transported
order.  L287 proves that this convention already fails at the second
successor.

The next useful object is the **complete effective even face**:

1. the nonzero endpoint compression of
   \(\mathcal G_S(H_{2k}^0)\);
2. the direct transfer Gram supplied by L283/L234;
3. the Schur cross-squares generated by earlier odd coefficients; and
4. the endpoint response of any selected perpendicular correction.

The induction must factor or bound their sum on the surviving partial
flag.  At odd successors, an endpoint-null homological cancellation
like L286 may still be the correct step.  At even successors, the
trace-carrying face is physical data and must be budgeted, not gauged
away.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_sextic_transport_trace_obstruction.py \
  --output \
  experiments/repeated_crabb_sextic_transport_trace_obstruction_s70225.jsonl
```

The checker reconstructs (3)--(9) in exact word algebra, verifies the
triple-delay divisibility of the parallel correction, evaluates (11)
with SymPy rationals, and solves the Stein equation exactly.  The
tracked dataset has SHA-256

```text
372921b75a883e3066e1d910082e5cb805d4af7b2e37d2999a9863daa3f5bb1c
```
