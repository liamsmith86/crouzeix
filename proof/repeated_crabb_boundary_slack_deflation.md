# The delayed boundary-slack face and its grade-one trace

## 1. Status (L225a and candidate L225b, 2026-07-24)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\qquad B_j=W^*(S^*)^jV,
\]

and the repeated equality metric

\[
P=2I-E+2F.
\]

With \(q=c^2\), let

\[
R(c)=P_{\rm bl}(c)
=I-\sum_{j\ge1}\frac{q^j}{1+q^j}(S^*)^jES^j
 \sum_{j\ge1}q^jS^jF(S^*)^j.                       \tag{1}
\]

Suppose

\[
B_1=\cdots=B_{k-1}=0.                               \tag{2}
\]

Two parts of the first delayed face are now exact.

1. L223 gives the physical upper metric face

   \[
   -4c^{2k}B_kB_k^*.
   \]

2. The matching lower metric face is

   \[
   \boxed{c^{2k}B_k^*B_k.}                         \tag{3}
   \]

There is also a sharply isolated remaining identity.  Let

\[
\widehat T(c)
=P^{1/2}\phi_c(T+cT^*)P^{-1/2},
\qquad T=P^{-1/2}SP^{1/2},
\]

and form the balanced Stein slack

\[
H(c)=R(c)-\widehat T(c)^*R(c)\widehat T(c).          \tag{4}
\]

Put \(Q=I-E\) and take its Schur complement away from the right
defect:

\[
\boxed{
\mathcal K_S(c)
=QH(c)Q-QH(c)V\{V^*H(c)V\}^{-1}V^*H(c)Q.}          \tag{5}
\]

For grade one,

\[
\boxed{
\operatorname {tr}[c^2]\mathcal K_S(c)
=2\|B_1\|_F^2.}                                    \tag{6}
\]

Equation (6) is proved below.  On every fully delayed grade tested,
the stronger covariance

\[
\boxed{
[c^{2k}]\mathcal K_S(c)
=J\bigl([c^2]\mathcal K_{S_{k-1}}(c)\bigr)J^*}     \tag{7}
\]

holds to binary64 precision.  Here \(S_{k-1}\) is L209's colligation
after removing \(W,\ldots,S^{k-2}W\), and \(J\) includes its retained
state space.  Equation (7) is still a **candidate identity**, not a
proved lemma.

If (7), or only its trace, is proved, then (3), L223, and (6) imply
the all-grade effective trace law

\[
\boxed{
\operatorname {tr}E_{2k,\rm eff}
=-16\|B_k\|_F^2,}                                  \tag{8}
\]

which is precisely the remaining pointwise sign required by L222.
Thus the matrix-Gram bridge has been reduced to one state-space
Schur-defect covariance.

## 2. Exact lower boundary face

L219 gives

\[
\begin{aligned}
R(c)-P^{-1}
={}&\sum_{j\ge1}\left(\frac12-\frac{q^j}{1+q^j}\right)
  (S^*)^jES^j\\
&+\frac14F+\sum_{j\ge1}q^jS^jF(S^*)^j.             \tag{9}
\end{aligned}
\]

The constant kernel of (9) is \(V\mathbb C^m\).  Since \(SV=0\),

\[
V^*(S^*)^jES^jV=0.
\]

Also \(V^*FV=0\), while

\[
\begin{aligned}
V^*S^jF(S^*)^jV
&=(W^*(S^*)^jV)^*(W^*(S^*)^jV)\\
&=B_j^*B_j.                                        \tag{10}
\end{aligned}
\]

Therefore (2) gives

\[
V^*\{R(c)-P^{-1}\}V
=q^kB_k^*B_k+O(q^{k+1}).                           \tag{11}
\]

The cross row also begins at \(q^k\), so its Schur square begins at
\(q^{2k}=O(q^{k+1})\).  Physical congruence contributes no factor at
the lower defect because \(P^{1/2}V=V\).  This proves (3), including
its right-Gram orientation.

## 3. The grade-one slack Schur trace

Write

\[
\widehat T(c)=S+c\dot S+c^2\ddot S+O(c^3).
\]

The Riemann coefficients give

\[
\begin{aligned}
\dot S&=J-S^3,\qquad J=PS^*P^{-1},\\
\ddot S
&=2S-(S^2J+SJS+JS^2)+S^5.                         \tag{12}
\end{aligned}
\]

From (1),

\[
R(c)=I+c^2R_2+O(c^4),
\qquad R_2=SFS^*-S^*ES.                            \tag{13}
\]

Expand (4):

\[
\begin{aligned}
H_0&=E,\\
H_1&=-(\dot S^*S+S^*\dot S),\\
H_2&=R_2-\ddot S^*S-S^*\ddot S-\dot S^*\dot S
     -S^*R_2S.                                     \tag{14}
\end{aligned}
\]

Because \(V^*H_0V=I\), the second coefficient of (5) is

\[
[c^2]\mathcal K_S
=QH_2Q-AA^*,
\qquad A=QH_1V=-S^*\dot SV.                        \tag{15}
\]

Take the scalar trace.  Since \(S=SQ\) and \(S^*=QS^*\),

\[
\begin{aligned}
\operatorname {tr}[c^2]\mathcal K_S
={}&\operatorname {tr}\{R_2(F-E)\}\\
&-2\operatorname {Re}\operatorname {tr}(S^*\ddot S)
-\|\dot SQ\|_F^2-\|A\|_F^2.                        \tag{16}
\end{aligned}
\]

L203's order-preserving partial-isometry reduction gives, with

\[
t=\|W^*S^*V\|_F^2=\|B_1\|_F^2,
\]

\[
\begin{aligned}
\|A\|_F^2+\|\dot SQ\|_F^2
&=\|\dot S\|_F^2-16t,\\
\|\dot S\|_F^2
+2\operatorname {Re}\operatorname {tr}(S^*\ddot S)
&=12t.                                             \tag{17}
\end{aligned}
\]

It remains only to insert the boundary-metric term.  The relations

\[
FS=0,\quad S^*F=0,\quad SE=0,\quad ES^*=0
\]

give

\[
\begin{aligned}
\operatorname {tr}(R_2F)&=-t,\\
\operatorname {tr}(R_2E)&= t.
\end{aligned}                                      \tag{18}
\]

Hence

\[
\operatorname {tr}\{R_2(F-E)\}=-2t.
\]

Substitute this and (17) in (16):

\[
-2t-12t+16t=2t.
\]

This proves (6) without commuting copy matrices.

## 4. Why the trace \(2\) would finish the pointwise sign

At the first face, replacing the boundary metric by \(R+c^{2k}X\)
changes the balanced Stein slack by

\[
X-S^*XS.
\]

To remove the complementary Schur defect, its \(Q\)-corner must
contain

\[
-K_k,\qquad K_k=[c^{2k}]\mathcal K_S.              \tag{19}
\]

Terms \(VC^*+CV^*\) change the rank-\(m\) factor but not the following
trace identity.  The two orbit resolutions imply, for
\(\mathcal L_S(X)=X-S^*XS=Y\),

\[
\begin{aligned}
\operatorname {tr}(W^*XW)&=\operatorname {tr}Y,\\
\operatorname {tr}(V^*XV)&=\operatorname {tr}(V^*YV).
\end{aligned}
\]

After physical congruence,

\[
\boxed{
\Delta_{\rm up}-4\Delta_{\rm low}
=-4\operatorname {tr}K_k.}                         \tag{20}
\]

The boundary metric starts with

\[
\operatorname {tr}E_{\rm up}=-4\|B_k\|_F^2,\qquad
\operatorname {tr}E_{\rm low}=+\|B_k\|_F^2.
\]

Retightening the lower endpoint sets
\(\Delta_{\rm low}=-\|B_k\|_F^2\).  If
\(\operatorname {tr}K_k=2\|B_k\|_F^2\), equation (20) gives

\[
\Delta_{\rm up}=-12\|B_k\|_F^2.
\]

Adding the original upper face yields (8).

This argument determines the trace of every admissible repair; it
does not construct a bounded repair through rank jumps.  L222 still
requires that separate selection step after the trace sign closes.

## 5. Evidence for the deflation covariance

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_boundary_slack_deflation.py \
  --output \
  experiments/repeated_crabb_boundary_slack_deflation_s70224.jsonl
```

The deterministic audit covers fully delayed unstructured
colligations through grade five.  It checks:

1. every earlier coefficient of (5) vanishes;
2. the leading matrix coefficient equals the embedded deflated
   grade-one coefficient in (7);
3. its trace is \(2\|B_k\|_F^2\);
4. both boundary faces have the orientations in L223 and (3); and
5. the reconstructed total trace is \(-16\|B_k\|_F^2\).

Only Sections 2--4 are proved.  The finite-grade audit is evidence
for (7), whose all-grade state-word or Schur-model proof remains the
next task.  The tracked dataset SHA-256 is
`e60e26045252d1ba6270b77f94032168eebbad21aa5d2908f34bb890fb9fdf38`.
