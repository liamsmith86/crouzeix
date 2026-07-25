# The retained first-reflection jet is universal through the slack face

## 1. Result (L242, 2026-07-24)

Retain a pure partial isometry

\[
I-S^*S=E,\qquad I-SS^*=F,\qquad EF=0,
\]

and put

\[
J=S^*(I+E),\qquad
\Xi_-(c)=S+cJ.
\]

For a complete removed delay of length \(r\), L237 gives the retained
direct-map block

\[
\mathcal A_r(c)
=\frac1{2\pi i}\oint
\phi_c(z)
\{zI-\Xi_-(c)-\gamma_r(z,c)F\}^{-1}\,dz.          \tag{1}
\]

Let

\[
\eta(z,c)=\frac{z-\sqrt{z^2-4c}}2
\]

be the half-line self-energy, and define \(\mathcal A_\infty\) by
replacing \(\gamma_r\) with \(\eta\) in (1).

For every \(r\geq2\), the complete difference through the degree
needed by L228 is

\[
\boxed{
\mathcal A_r-\mathcal A_\infty
=(-1)^rc^{2r}
\{D_0+cD_1+c^2D_2\}
+O(c^{2r+3}),}                                    \tag{2}
\]

where the following are identities in the free word algebra on
\(S,J,F\):

\[
D_0=\sum_{a+b=1}S^aFS^b,                          \tag{3}
\]

\[
D_1=JF+FJ-\sum_{a+b=3}S^aFS^b,                   \tag{4}
\]

and

\[
\boxed{
\begin{aligned}
D_2={}&3D_0+\sum_{a+b=5}S^aFS^b\\
&-\sum_{a+b+d=2}
\{S^aJS^bFS^d+S^aFS^bJS^d\}\\
&-2\sum_{a+b+d=1}S^aFS^bFS^d.
\end{aligned}}                                    \tag{5}
\]

No coefficient in (3)--(5) depends on \(r\).  After imposing the
partial-isometry relations, they simplify further to

\[
\boxed{
\begin{aligned}
D_0={}&SF,\\
D_1={}&FJ-S^3F,\\
D_2={}&SF+S^5F
-\sum_{a+b=2}S^aJS^bF
-\sum_{a+d=2}S^aFJS^d.
\end{aligned}}                                    \tag{6}
\]

Thus L239's exact contour collapse has a fixed state-word jet in the
entire target window; no new delay-dependent direct-map polynomial
appears at higher grade.

L242 is **not** L228.  The retained block (1) is only one part of the
Stein product.  Amplitude entering the removed chain also contributes
to \(A^*P_{\rm bl}A\); omitting it leaves the false half-line residual
found during this calculation.  L240 closes that missing flux
exactly.  Equations (2)--(6) should be used only inside the complete
L238--L240 state-lifted calculation.

## 2. Scalar coefficients in the target window

Write

\[
\phi_c(w)=\sum_{n\geq0}a_n(c)w^{2n+1}.
\]

Extracting the first two nonzero diagonals from L125's theta/ODE
recurrence gives, for every \(n\geq1\),

\[
\boxed{
a_n(c)=(-1)^n
\{c^n+(2n+1)c^{n+2}\}
+O(c^{n+4}).}                                     \tag{7}
\]

The exceptional linear coefficient is

\[
a_0(c)=1+2c^2+O(c^4),
\]

which is the grade-one correction already shown to be load-bearing
in L235.  For completeness, (7) is an all-\(n\) identity rather than
a pattern read from finitely many coefficients.  In L125's Newton
equation use

\[
\frac{k}{c}=4-16c^2+O(c^4),\qquad
\alpha^2=1-8c^2+O(c^4),
\]

and write its inverse-map blow-up as \(H=H_0+c^2H_1+\cdots\).
Parameterize

\[
x=\frac{y}{(1+y)^2},\qquad H_0(x)=1+y.
\]

The unique next triangular coefficient satisfies

\[
H_1\!\left(\frac{y}{(1+y)^2}\right)
=-\frac{(2+y+y^2)(1+y)}{1-y}.
\]

Series reversion therefore gives the two-diagonal generating identity

\[
\phi_c(w)
=\frac{w}{1+cw^2}
+c^2w\frac{2+cw^2+c^2w^4}{(1+cw^2)^2}
+\text{terms four weighted degrees later}.
\]

Expanding the last display gives the exceptional \(a_0\) and (7) for
every \(n\geq1\).

L237's scalar error recursion gives one more required coefficient:

\[
\boxed{
\gamma_r-\eta
=\frac{c^r}{z^{2r-1}}
+(2r-2)\frac{c^{r+1}}{z^{2r+1}}
+O(c^{r+2})
\qquad(r\geq2).}                                  \tag{8}
\]

Finally,

\[
\eta=\frac cz+\frac{c^2}{z^3}+O(c^3).             \tag{9}
\]

Equations (7)--(9) are the only scalar data entering (2).

## 3. Laurent-residue count

Linearize the resolvent in

\[
\delta_r=\gamma_r-\eta:
\]

\[
R_{\gamma_r}-R_\eta
=R_\eta\delta_rFR_\eta
+R_\eta\delta_rFR_\eta\delta_rFR_\eta+\cdots .    \tag{10}
\]

Two \(\delta_r\) insertions first contribute at degree \(c^{4r}\).
For \(r\geq2\),

\[
4r>2r+2,
\]

so only the term linear in \(\delta_r\) can enter (2).  This is a
proved filtration for the retained direct block; it is not a license
to linearize the later endpoint inverse in the full Stein problem.

Expand each \(R_\eta\) in the three letters

\[
S,\qquad cJ,\qquad
\eta F=\left(\frac cz+\frac{c^2}{z^3}+\cdots\right)F.
\]

The residue condition leaves the following exhaustive table.

\[
\begin{array}{c|c|c}
\text{relative degree}&\text{source}&\text{word family}\\ \hline
0&a_r,\ \delta_r^{(0)}&D_0\\
1&a_r\text{ with one }J;\ a_{r+1}&D_1\\
2&a_r^{(n+2)};\ \delta_r^{(1)};
   a_{r+2};\ a_{r+1}\text{ with one }J;
   a_{r+1}\text{ with one }\eta F&D_2.
\end{array}                                       \tag{11}
\]

At relative degree two, the correction in \(a_r\) contributes
\((2r+1)D_0\), while the second term of (8), paired with
\(a_{r+1}\), contributes \(-(2r-2)D_0\).  Their sum is the universal
\(3D_0\) in (5).  A single background \(\eta F\) can lie on either
side of the distinguished reflection and gives the factor \(2\) in
the last line of (5).  The remaining residue placements give the
two middle sums.  This proves (2)--(5).

## 4. Partial-isometry collapse

Use

\[
FS=0,\qquad S^*F=0,\qquad EF=0.
\]

Then (3) leaves only \(SF\).  In (4),

\[
JF=0,\qquad
\sum_{a+b=3}S^aFS^b=S^3F,
\]

which gives the second line of (6).  In (5), every word after the
last \(F\) vanishes, and the double-\(F\) sum reduces to \(SF\).
Combining its coefficient \(-2\) with \(3D_0\) gives the first term
of the last line of (6).  This proves the simplified form.

## 5. Independent exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_reflected_direct_jet.py \
  --output \
  experiments/repeated_crabb_reflected_direct_jet_s70224.jsonl
```

The checker:

1. generates the scalar direct coefficients from L125's exact
   theta/ODE recurrence;
2. generates \(\gamma_r-\eta\) from the continued-fraction recursion;
3. enumerates every Laurent-resolvent placement through relative
   degree two in the unreduced free-word algebra; and
4. compares exactly with (3)--(5).

Delays two through five all pass with rational residual zero.  A
separate symbolic-delay enumeration verifies that every polynomial
in \(r\geq3\) cancels to (3)--(5) before the finite records are
written.  The three universal polynomials have respectively
\(2,6,23\) unreduced words.  The tracked SHA-256 is
`d9b408c36052d76c970436847f791131db22f83718fb4d212456667adada6824`.
