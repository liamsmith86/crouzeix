# The all-size elliptic defect-Hessian edge (L139, 2026-07-23)

## 1. Statement

Fix an L117 elliptic Crabb axis of length `L` and use the coefficient
defect scale of L118, so the zeroth defect coordinate is fixed and the
free coordinates are

\[
x=(x_1,\ldots,x_L).
\]

Let `H_c` be the real quadratic matrix of the pure-defect part of the
rank-one condition envelope.  Then:

1. `H_c` is exactly block diagonal by the parity of the coordinate
   index;
2. its diagonal Newton constants are

   \[
   (H_c)_{ii}=4+o(1)\quad(i<L),\qquad
   (H_c)_{LL}=\frac83+o(1);
   \]

3. if `1<=i<j<=L` and `j-i=2r`, then

   \[
   (H_c)_{ij}=
   \begin{cases}
   8c^r+O(c^{r+1}),&j<L,\\[2mm]
   \dfrac{16}{3}c^r+O(c^{r+1}),&j=L.
   \end{cases}                                      \tag{1}
   \]

In fact the error gains two powers in the nonterminal Fourier leg;
the weaker displayed remainder is enough for the factorization.

Consequently, if

\[
H_c={\cal L}_c{\cal D}_c{\cal L}_c^*
\]

is its exact unit-lower LDL* factorization, then

\[
\begin{aligned}
({\cal D}_c)_{ii}&=4+o(1) &&(i<L),\\
({\cal D}_c)_{LL}&=\frac83+o(1),\\
({\cal L}_c)_{ji}&=
\begin{cases}
2c^{(j-i)/2}+O(c^{(j-i)/2+1}),&j<L,\\[1mm]
\dfrac43c^{(L-i)/2}+O(c^{(L-i)/2+1}),&j=L,
\end{cases}
\end{aligned}                                      \tag{2}
\]

for same-parity `i<j`, while the opposite-parity entries vanish.
Thus the interior outer factor has associated series

\[
\boxed{
1+2c{\cal S}^2+2c^2{\cal S}^4+\cdots
=\frac{1+c{\cal S}^2}{1-c{\cal S}^2}.
}                                                   \tag{3}
\]

This proves the Hessian-factor half of A98.  It does not yet insert the
operator/Faber linear forcing or prove the completed square.

## 2. Exact DCT modal reduction

Use L138 with

\[
P_0=\operatorname{diag}(p_0,\ldots,p_L),\qquad
p_0=1,\quad p_m\longrightarrow2\ (0<m<L),\quad
p_L\longrightarrow4.
\]

For a homogeneous defect multiplier `s`, write its sampled cosine
expansion as

\[
s(\theta_n)=\sum_{t=0}^L g_t\cos(t\theta_n),\qquad
\theta_n=\frac{\pi n}{L}.                            \tag{4}
\]

Put `H_s=U diag(s)U^*` and `B=D^{-1}H_sD`.  DCT-I
orthogonality and the endpoint identity

\[
\cos(L\theta_n)\cos(m\theta_n)
=\cos((L-m)\theta_n)
\]

give, for `0<t<L`,

\[
\begin{array}{ll}
B_{0t}=c^{t/2}g_t/\sqrt2,&
B_{t0}=c^{-t/2}g_t/\sqrt2,\\[1mm]
B_{L,L-t}=c^{-t/2}g_t/\sqrt2,&
B_{L-t,L}=c^{t/2}g_t/\sqrt2.                         \tag{5}
\end{array}
\]

For `t=L`, the same formulas hold without the factors `1/sqrt(2)`.
The mode `g_0` is the irrelevant common defect-scale direction.

L138's endpoint Schur formula uses only the zeroth and `L`th rows of
`B`.  Equations (5) therefore prove an exact fact: **different DCT
modes never couple in the pure-defect condition Hessian**.  Over
complex coordinates each mode is a real `2 x 2` phase block, but
there are still no cross terms between distinct `g_t`.

For example, on the real part of a nonterminal mode its exact scalar
weight is

\[
\begin{aligned}
w_t(c)=\frac12\bigg[
&p_{L-t}c^{-t}
+\frac{(p_{L-t}c^{-t/2}+p_Lc^{t/2})^2}
       {p_L-p_{L-t}}\\
&-p_Lp_tc^t
-p_L\frac{(p_tc^{t/2}+c^{-t/2})^2}
             {1-p_t}
\bigg].                                             \tag{6}
\end{aligned}
\]

Since the axis levels tend to `(1,2,...,2,4)`,

\[
w_t(c)=4c^{-t}+o(c^{-t})\quad(0<t<L).                \tag{7}
\]

For the endpoint mode, the missing `1/2` in (5) gives instead

\[
w_L(c)=\frac83c^{-L}+o(c^{-L}).                      \tag{8}
\]

The leading terms (7)--(8) are phase isotropic.  Phase-dependent
pieces contain both the small and large entries in (5), so they are
strictly above this edge.

## 3. Jacobi multiplication is the LDL outer factor

Let `q=c^2` be the Jacobi nome.  The defect-to-homogeneous change in
L138 is

\[
s(\theta_n)
=a_c(\theta_n)\sum_{i=1}^Lc^{i/2}x_i\cos(i\theta_n),
                                                            \tag{9}
\]

up to a common scalar tending to one, where

\[
a_c(\theta)=\frac{\operatorname{dn}(2K\theta/\pi)}{k'}.
\]

The standard Jacobi Fourier series is

\[
a_c(\theta)=a_0(c)\left[
1+4\sum_{\nu\ge1}
\frac{c^{2\nu}}{1+c^{4\nu}}\cos(2\nu\theta)
\right].                                             \tag{10}
\]

Sampling aliases in (10) are higher order for every frequency
strictly below `L`.  Product-to-sum now gives, in the associated
Newton filtration,

\[
\boxed{
c^{-j/2}g_j
=x_j+
2\sum_{\substack{i<j\\j-i\ {\rm even}}}
c^{(j-i)/2}x_i
+\text{strictly higher terms}.
}                                                   \tag{11}
\]

The coefficient two is not fitted: it is one half of the Fourier
coefficient four in (10).  All Fourier frequencies are even, so
(11) also proves exact parity separation.

For completeness, consider `i<j` with `j-i=2r`.  In the modal sum,
the unique contribution at weight `c^r` is mode `t=j`:

\[
\begin{aligned}
c^{-j/2}g_j^{(j)}&=x_j+\cdots,\\
c^{-j/2}g_j^{(i)}&=2c^rx_i+\cdots.
\end{aligned}
\]

Every `t>j` gains at least `2(t-j)/2`, every `t<j` gains through the
downward Fourier leg, and every finite DCT alias uses a larger Jacobi
frequency.  Combining this unique term with (7)--(8) gives a cross
coefficient `16c^r` in the quadratic form for `j<L`, hence
`(H_c)_{ij}=8c^r`; at `j=L` it gives `32c^r/3`, hence
`(H_c)_{iL}=16c^r/3`.  This proves (1).

## 4. LDL edge

In the unit-lower LDL recursion,

\[
({\cal L}_c)_{ji}
=\frac{(H_c)_{ji}
-\sum_{k<i}({\cal L}_c)_{jk}({\cal D}_c)_{kk}
             \overline{({\cal L}_c)_{ik}}}
       {({\cal D}_c)_{ii}}.                           \tag{12}
\]

For same parity, every summand with `k<i` has valuation

\[
\frac{j-k}{2}+\frac{i-k}{2}
=\frac{j-i}{2}+(i-k),
\]

strictly above the edge.  Divide (1) by the diagonal constant four.
This gives coefficient two in every interior row and coefficient
`(16/3)/4=4/3` in the terminal row, proving (2)--(3).

## 5. Independent regression

`experiments/crabb_defect_hessian_factor.py` reconstructs the complete
pure-defect Hessian from the original Stein recurrence over exact
rational truncated series, without using (4)--(12).  Its persisted
length-eight audit checks:

- nine interior and three terminal same-parity edges;
- constants `4` and `8/3`;
- coefficients `2` and `4/3`;
- absence of every lower and opposite-parity term.

That finite regeneration is a guard against normalization mistakes;
the DCT/Jacobi argument above is the all-size proof.
