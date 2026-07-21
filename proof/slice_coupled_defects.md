# Stein-defect parameterization of the coupled slice obstruction

The two modal norm bounds close every dual ray supported on one parity block.
This note eliminates the remaining contraction LMIs and classifies the possible
rank faces of a genuinely coupled optimum.  The result does not yet prove
`L20`; it reduces the live part to two positive \(2\times2\) defect matrices,
with rank at most one on every dual-active block.

## 1. Exact defect coordinates

Retain the modal notation of `slice_similarity_duality.md`:

\[
 G_o=U^TD^2U,\qquad G_e=V^TD^2V,\qquad
 \Sigma=\operatorname{diag}(\tau_1,\tau_2),
 \quad 0<\tau_2<\tau_1<1. \tag{1}
\]

For a modal metric pair \((H_o,H_e)\), define its two Stein defects by

\[
 Q_o=H_o-c\Sigma H_e\Sigma,
 \qquad
 Q_e=H_e-c^{-1}\Sigma H_o\Sigma. \tag{2}
\]

The two contraction inequalities are exactly \(Q_o,Q_e\succeq0\).
Eliminating first \(H_e\), then \(H_o\), gives

\[
 \begin{aligned}
 H_o-\Sigma^2H_o\Sigma^2&=Q_o+c\Sigma Q_e\Sigma,\\
 H_e-\Sigma^2H_e\Sigma^2&=Q_e+c^{-1}\Sigma Q_o\Sigma.
 \end{aligned} \tag{3}
\]

Let \(\mathcal K\) denote Schur multiplication by the two-node Szegő
kernel

\[
 K_{ij}=\frac1{1-\tau_i^2\tau_j^2}. \tag{4}
\]

Because every denominator is positive, (3) has the unique solution

\[
 \boxed{\begin{aligned}
 H_o&=\mathcal K(Q_o+c\Sigma Q_e\Sigma),\\
 H_e&=\mathcal K(Q_e+c^{-1}\Sigma Q_o\Sigma).
 \end{aligned}} \tag{5}
\]

Conversely, (5) constructed from any \(Q_o,Q_e\succeq0\) has precisely the
defects (2).  Thus (5) parameterizes **every** contraction metric; no
contraction inequality remains to be checked.  The kernel matrix in (4) is
positive definite, so positive semidefiniteness of the resulting metrics
follows from the Schur product theorem.  The sandwich constraints below
enforce positive definiteness.

Consequently the complete slice target is equivalent to finding nonzero
\(Q_o,Q_e\succeq0\) for which

\[
 G_o\preceq H_o\preceq4G_o,
 \qquad
 G_e\preceq H_e\preceq4G_e, \tag{6}
\]

with \(H_o,H_e\) given explicitly by (5).  This is a two-node Pick-kernel
form of `L20`.

## 2. KKT rank-face classification

Let \(P\) be an optimal primal metric and let
\(Z=\operatorname{diag}(Z_o,Z_e)\) be an optimal parity-reduced dual witness.
Strong duality was proved in the modal note.  Complementary slackness for the
contraction constraint gives

\[
 \operatorname{tr}(Z_o\widehat Q_o)=0,
 \qquad
 \operatorname{tr}(Z_e\widehat Q_e)=0, \tag{7}
\]

where \(\widehat Q_o,\widehat Q_e\) are the physical-coordinate contraction
defects.  They are congruent to (2), so their ranks agree with those of
\(Q_o,Q_e\).  Since both factors in each trace are positive semidefinite,

\[
 \operatorname{rank}Z_i+\operatorname{rank}Q_i\le2
 \quad(i=o,e). \tag{8}
\]

Both dual blocks cannot be positive definite.  If they were, (8) would force
\(Q_o=Q_e=0\), hence \(P=T^*PT\).  Iterating and using \(r(T)<1\) would give
\(P=0\), contradicting \(P\succeq I\).

For a potentially sharp witness there are three cases, up to swapping parity:

1. **one-block:** \(Z_o=0\) (a maximizer may be chosen with
   \(\operatorname{rank}Z_e=1\));
2. **rank-one/rank-one:** both dual blocks have rank one, so both defects have
   rank at most one;
3. **rank-one/full:** one dual block has rank one and the other rank two, so
   the corresponding defects have ranks at most one and zero.

The first face is already closed.  If \(Z_o=0\), then

\[
 D_o=-BZ_eB^*,\qquad D_e=Z_e,
\]

and its trace ratio is at most \(\|B\|^2\le4\); the other parity uses
\(\|C\|^2\le4\).  A maximizing one-block witness may be chosen rank one.

Hence any hypothetical optimum above four must be on face 2 or 3.  On both
faces the metric in (5) has

\[
 Q_o=aa^T,\qquad Q_e=bb^T, \tag{9}
\]

where either vector may be zero.  Up to a common scale, (9) has only three
real parameters: two projective angles and one relative magnitude.  This is
the exact remaining coupled search space.

## 3. Explicit rank-one formula

Writing \(a=(a_1,a_2)^T\) and \(b=(b_1,b_2)^T\), formula (5) becomes

\[
 (H_o)_{ij}=\frac{a_ia_j+c\tau_i\tau_jb_ib_j}
 {1-\tau_i^2\tau_j^2},\qquad
 (H_e)_{ij}=\frac{b_ib_j+c^{-1}\tau_i\tau_ja_ia_j}
 {1-\tau_i^2\tau_j^2}. \tag{10}
\]

The physical metric blocks are recovered by

\[
 P_o=(DU)^{-T}H_o(DU)^{-1},\qquad
 P_e=(DV)^{-T}H_e(DV)^{-1}. \tag{11}
\]

Thus the coupled theorem has no hidden operator variables: it asks whether a
KKT minimizer on a coupled face can have condition number above four.  The
default SDP cases reproduce their optima from (10) to numerical precision on
the coupled faces; `experiments/slice_coupled_defects.py` audits the exact
reconstruction and the KKT rank restrictions.

## 4. Transfer reduction of the rank-one/rank-one face

The remaining face also has an exact colligation reduction.  Write

\[
 Z_o=xx^T,\qquad Z_e=yy^T,
\]

and take the spectral positive and negative factors of its two dual blocks:

\[
 \begin{aligned}
 xx^T-(By)(By)^T&=p_op_o^T-n_on_o^T,\\
 yy^T-(Cx)(Cx)^T&=p_ep_e^T-n_en_e^T.
 \end{aligned} \tag{12}
\]

Each difference of two rank-one matrices has at most one positive and one
negative eigenvalue, so this includes the semidefinite limits by allowing a
factor to vanish.  After harmless sign choices, equality of the corresponding
row Gramians gives two rotations.  Thus, for some \(a,b\in[-1,1]\), with
\(s=\sqrt{1-a^2}\) and \(t=\sqrt{1-b^2}\),

\[
 \begin{aligned}
 By&=ax+sn_o,&p_o&=-sx+an_o,\\
 Cx&=by+tn_e,&p_e&=-ty+bn_e.
 \end{aligned} \tag{13}
\]

Stack

\[
 X=\binom{x}{y},\quad P=\binom{p_o}{p_e},\quad
 N=\binom{n_o}{n_e},\quad
 A=\operatorname{diag}(aI_2,bI_2),\quad
 S=\operatorname{diag}(sI_2,tI_2). \tag{14}
\]

Equations (13) become

\[
 TX=AX+SN,\qquad P=-SX+AN. \tag{15}
\]

Multiplying the second equation by \(S\) and using \(S^2=I-A^2\) gives

\[
 SP=-(I-AT)X. \tag{16}
\]

For \(|a|,|b|<1\), \(I-AT\) is invertible: the eigenvalues of
\((AT)^2\) are \(ab\tau_i^2\), whose moduli are strictly below one.  Hence

\[
 N=\mathscr R_{a,b}(T)P,
 \qquad
 \mathscr R_{a,b}(T)
 =-S^{-1}(T-A)(I-AT)^{-1}S. \tag{17}
\]

Since the factors in (12) are spectral,

\[
 {\operatorname{tr}D_-\over\operatorname{tr}D_+}
 ={\|N\|^2\over\|P\|^2}
 \le\|\mathscr R_{a,b}(T)\|^2. \tag{18}
\]

The apparent singularity at \(a=\pm1\) or \(b=\pm1\) disappears after block
elimination.  Put \(F=BC\), \(G=CB\), and \(z=ab\).  Then exactly

\[
 \boxed{
 \mathscr R_{a,b}(T)=
 \begin{bmatrix}
 (aI-bF)(I-zF)^{-1}&-st(I-zF)^{-1}B\\
 -st(I-zG)^{-1}C&(bI-aG)(I-zG)^{-1}
 \end{bmatrix}.} \tag{19}
\]

Formula (19) extends continuously to the closed square.  Consequently the
single sufficient theorem

\[
 \boxed{\|\mathscr R_{a,b}(T)\|\le2
 \quad\text{for every }$a,b$\in[-1,1]^2} \tag{RT}
\]

closes the rank-one/rank-one face and therefore proves `L20`.

## 5. Matrix-valued Blaschke structure

The transfer is not an arbitrary two-parameter rational matrix.  Let

\[
 L=\operatorname{diag}(D,\sqrt cD),\qquad
 \mathcal U=\operatorname{diag}(U,V),\qquad
 T_0=\begin{bmatrix}0&\Sigma\\\Sigma&0\end{bmatrix}.
\]

The modal formula gives \(T=L\mathcal U T_0\mathcal U^TL^{-1}\).  Because
\(A\) and \(S\) are scalar on each parity block, they commute with both
\(L\) and \(\mathcal U\).  Therefore

\[
 \mathscr R_{a,b}(T)
 =L\mathcal U\mathscr R_{a,b}(T_0)\mathcal U^TL^{-1}. \tag{20}
\]

After grouping the two coordinates belonging to the same node,
\(\mathscr R_{a,b}(T_0)\) is the direct sum of

\[
 K_{a,b}(w)={1\over1-abw^2}
 \begin{bmatrix}
 a-bw^2&-st w\\
 -st w&b-aw^2
 \end{bmatrix},qquad w=\tau_1,\tau_2. \tag{21}
\]

This is a rational \(2\times2\) inner function.  Direct algebra gives

\[
 K_{a,b}(1/z)^TK_{a,b}(z)=I,qquad
 \det K_{a,b}(z)={ab-z^2\over1-abz^2}. \tag{22}
\]

Thus the determinant is a degree-two scalar Blaschke product, but the
interior transfer is genuinely matrix-valued.  It does, however, close the
large-$c$ range immediately.  A rational inner matrix is contractive at every
point of the disk, so (21) gives

\[
 \|\mathscr R_{a,b}(T_0)\|\le1.
\]

The two orthogonal modal rotations do not change this norm, while the
similarity in (20) has

\[
 \operatorname{cond}(L)=c^{-3/2}.
\]

Consequently

\[
 \boxed{\|\mathscr R_{a,b}(T)\|\le c^{-3/2}\le2
 \quad\text{whenever }c\ge2^{-2/3}.} \tag{22a}
\]

This uses the complete matrix-valued inner structure, not a scalar boundary
reduction, and removes every parameter $a,b,r,u$ in that nome range.

If \(a=\pm1\) or \(b=\pm1\),
(19) becomes a direct sum of a signed identity and a scalar automorphism of
\(F\) or \(G\); this boundary is already contained in the proved even sector.
At \(a=b=0\), the transfer is \(-T\), so L24 and L26 close that point.  The
new content of (RT) is the open interior of the parameter square.

`experiments/slice_rank_one_transfer.py` audits (19) and (22) with exact
rational/symbolic arithmetic.  Its deterministic floating-point search is
sharp on the two default rank-one/rank-one SDP cases: the squared transfer
norms \(3.136329347\) and \(2.491664377\) reproduce the primal-dual optima.
This sharpness is evidence only, not the proof of (RT).

## 6. Polynomial transfer defect

There is a second exact simplification which removes all resolvents from the
norm-two test.  Put $X_0=I-AT$.  Expanding both sides and using that $A$
commutes with $S=(I-A^2)^{1/2}$ gives the operator-ball identity

\[
 I-\mathscr R_{a,b}(T)^*\mathscr R_{a,b}(T)
 =S X_0^{-*}(I-T^*T)X_0^{-1}S. \tag{24}
\]

Indeed, after multiplying by $X_0^*$ and $X_0$, the identity reduces to

\[
 X_0^*S^{-2}X_0-(T^*-A)S^{-2}(T-A)=I-T^*T, \tag{25}
\]

whose mixed terms cancel and whose two remaining coefficients use
$S^{-2}(I-A^2)=I$.  With $Y=X_0^{-1}S$, (24) implies

\[
 4I-\mathscr R_{a,b}(T)^*\mathscr R_{a,b}(T)
 =Y^*\mathcal Q_{a,b}(T)Y,
 \qquad
 \mathcal Q_{a,b}(T)=I-T^*T+3X_0^*S^{-2}X_0. \tag{26}
\]

Thus (RT) is equivalent in the open parameter square to
$\mathcal Q_{a,b}(T)\succeq0$.  Multiplication by the positive scalar
$s^2t^2$ makes this a polynomial LMI.  In parity blocks it is

\[
 s^2t^2\mathcal Q_{a,b}(T)=
 \begin{bmatrix}E_o&-3J\\-3J^T&E_e\end{bmatrix}, \tag{27}
\]

where

\[
 \begin{aligned}
 E_o&=t^2(4-a^2)I+s^2(4b^2-1)C^*C,\\
 E_e&=s^2(4-b^2)I+t^2(4a^2-1)B^*B,\\
 J&=at^2B+bs^2C^*.
 \end{aligned} \tag{28}
\]

Both diagonal blocks are already positive semidefinite from L24 and L26.  For
example, if $4b^2-1\ge0$, positivity of $E_o$ is immediate.  Otherwise
$C^*C\preceq4I$, so

\[
 E_o\succeq3(a^2+4b^2-5a^2b^2)I\succeq0. \tag{29}
\]

The last scalar expression is bilinear in $a^2,b^2\in[0,1]$ and is
nonnegative at all four corners.  Swapping $a,b$ gives
$E_e\succeq3(4a^2+b^2-5a^2b^2)I\succeq0$ in its only nontrivial sign
range.

There is also a direct two-square decomposition of the entire polynomial
core.  For column vectors $x,y$, expansion of (27)--(28) gives

\[
\begin{aligned}
 \left\langle
 \begin{bmatrix}E_o&-3J\\-3J^T&E_e\end{bmatrix}
 \binom{x}{y},\binom{x}{y}\right\rangle
={}&3t^2\|x-aBy\|^2+3s^2\|y-bCx\|^2\\
 &+s^2t^2\bigl(\|x\|^2-\|Cx\|^2
                    +\|y\|^2-\|By\|^2\bigr).
                                                        \tag{29a}
\end{aligned}
\]

Equivalently, with column concatenation understood,

\[
\begin{aligned}
 s^2t^2\mathcal Q_{a,b}(T)
={}&3t^2
 \binom{I}{-aB^T}\begin{bmatrix}I&-aB\end{bmatrix}
 +3s^2
 \binom{-bC^T}{I}\begin{bmatrix}-bC&I\end{bmatrix}\\
 &+s^2t^2\operatorname{diag}(I-C^TC,I-B^TB).
                                                        \tag{29b}
\end{aligned}
\]

This identity is exact and remains meaningful on the parameter boundary.
It immediately recovers the already closed edges: if $a=\pm1$, only the
first square remains, and if $b=\pm1$, only the second remains.  At
$a=b=0$, (29b) reduces to
$\operatorname{diag}(4I-C^TC,4I-B^TB)\succeq0$ by L24 and L26.  In the
interior, it isolates the sole possible loss as the central Stein-type
block defect $\operatorname{diag}(I-C^TC,I-B^TB)$; the two residual squares
must compensate for block singular values above one.  The identity is
symbolically audited in `experiments/slice_rank_one_transfer.py`.

Therefore the remaining content is precisely the off-diagonal Schur
inequality in (27), not positivity of its diagonal pieces.  This formulation
has no elliptic resolvents and is only quartic in $a,b$; it is the preferred
starting point for using the modal relation (23).  The exact audit of (24)--
(28) is included in `experiments/slice_rank_one_transfer.py`.

### 6.1 Hyperbolic square completion

The Schur coupling in (27) has an exact second-stage factorization.  For
$q\in(-1,1)$ define

\[
 A_q={2(4-q^2)\over1-q^2},\qquad
 D_q={2(4q^2-1)\over1-q^2},\qquad
 E_q={6q\over1-q^2}. \tag{30}
\]

The key scalar identity is

\[
 A_qD_q-E_q^2=-16. \tag{31}
\]

Put

\[
 \mathcal H=
 \begin{bmatrix}
 \sqrt{A_a}I&-E_aB/\sqrt{A_a}\\
 -E_bC/\sqrt{A_b}&\sqrt{A_b}I
 \end{bmatrix},\qquad
 \Delta=\operatorname{diag}(C/\sqrt{A_b},B/\sqrt{A_a}).
\]

Expanding and using (31) gives the exact square completion

\[
 2\mathcal Q_{a,b}(T)=\mathcal H^*\mathcal H-16\Delta^*\Delta. \tag{32}
\]

Let

\[
 u={E_a\over A_a}={3a\over4-a^2},\qquad
 v={E_b\over A_b}={3b\over4-b^2}.
\]

Both lie strictly between $-1$ and $1$.  Since the spectrum of $BC$ is
$\{\tau_1^2,\tau_2^2\}$, the matrix $\mathcal H$ is invertible.  Therefore
(27) is equivalent to the small-gain inequality

\[
 \|\mathcal S_{a,b}(T)\|\le1,\qquad
 \mathcal S_{a,b}(T)=4\Delta\mathcal H^{-1}. \tag{33}
\]

Writing $F=BC$ and $G=CB$ removes the inverse of the $4\times4$ block:

\[
 \mathcal S_{a,b}(T)=4
 \begin{bmatrix}
 {C(I-uvF)^{-1}\over\sqrt{A_aA_b}}&
 {uG(I-uvG)^{-1}\over A_b}\\
 {vF(I-uvF)^{-1}\over A_a}&
 {B(I-uvG)^{-1}\over\sqrt{A_aA_b}}
 \end{bmatrix}. \tag{34}
\]

This exposes a concrete scalar sufficient inequality.  If $S_{ij}$ denote
the four $2\times2$ blocks in (34), then

\[
 \boxed{\sum_{i,j=1}^2\|S_{ij}\|^2\le1} \tag{BE}
\]

implies (33): replace every block by its operator norm to obtain a scalar
$2\times2$ comparison matrix, whose operator norm is at most its Frobenius
norm.  Thus (BE) would prove (RT).  It is stronger than necessary, but it is
numerically robust: deterministic global searches over all five live real
parameters approach one only at the already sharp $c\to0$, $a,b\to0$
boundary.  The default rank-one/rank-one cases give $0.790699$ and $0.742709$.
More strongly, the search still passes when the exact conformal ratio $r=H(p)$
is replaced by every point of its rigorous cubic envelope

\[
 s_0p+a_3p^3\le r\le s_0p+(1-s_0)p^3. \tag{35}
\]

The cubic information is essential: the looser strip $s_0p\le r\le p$ has
numerical values as large as $1.3067$.  Thus a proof of (BE) may use only the
same algebraic theta data as L29, but cannot fall back to the soft inequality
$r\le p$.  This is evidence, not a certificate.  Formulae (30)--(34) are
exact and are audited by `experiments/slice_rank_one_transfer.py`; the two
optional searches are reproduced with `--global-energy` and
`--envelope-energy`.

### 6.2 Exact small-nome boundary of the block-energy target

The sharp-looking corner of (BE) can be evaluated without floating-point
optimization.  Set $p=0$ and $a=b=0$, and write $t=\tan\theta$ for the free
left modal angle and $g=\sqrt{k/c}$.  Direct substitution in the modal block
formula gives

\[
 B={g\over\sqrt{1+t^2}}
   \begin{bmatrix}1&0\\ct&0\end{bmatrix},\qquad
 C={g\over\sqrt{1+t^2}}
   \begin{bmatrix}c&t\\0&0\end{bmatrix}. \tag{36}
\]

At $a=b=0$ the off-diagonal blocks in (34) vanish and the other two are
$C/2$ and $B/2$.  Hence the angle cancels exactly:

\[
 \sum_{i,j}\|S_{ij}\|^2
 ={1\over4}(\|B\|^2+\|C\|^2)
 ={k(1+c^2)\over4c}. \tag{37}
\]

L23's nome inequality $k/c\le4/(1+c^2)^2$ therefore proves

\[
 \sum_{i,j}\|S_{ij}\|^2\le{1\over1+c^2}<1. \tag{38}
\]

Thus the exact boundary never exceeds one, although it tends to one as
$c\downarrow0$.  An independent 80-decimal audit at
$c=10^{-2},10^{-4},10^{-8},10^{-12},10^{-20}$ finds
$(1-\mathrm{BE})/c^2\to3$ and no excess in a scaled neighbourhood.  The old
coordinate $p=p_*+c^4x$ belongs to L29's different odd-block ridge; it is not
the exact equality boundary of L44.  The audit is reproduced by
`experiments/slice_rank_one_ridge_check.py`.  The finite neighbourhood scan is
falsification evidence only; (37)--(38) are the proved boundary statement.

### 6.3 A sharp product theorem for the two parity blocks

There is a second exact consequence of the conformal coupling which is not
visible from the separate bounds in L24 and L26:

\[
 \boxed{\|B\|\,\|C\|\le2.} \tag{39}
\]

This section proves (39).  It will close both coordinate axes of (RT).
Write the common modal matrix as

\[
 M=U\operatorname{diag}(\sqrt{k},p\sqrt{k})V^T,
 \qquad \tan v=r\tan u,
\]

where $p=\tau _2/\tau _1$ and $r=H(p)$.  If $m_{ij}$ are the entries of
$M$, then

\[
 {m_{12}\over m_{21}}=-\delta,qquad
 \delta={p-r\over1-pr},qquad
 w:=m_{21}^2
 ={k\tan ^2u(1-pr)^2\over(1+\tan ^2u)(1+r^2\tan ^2u)}. \tag{40}
\]

Consequently

\[
 0\le w\le {k(1-pr)^2\over(1+r)^2}. \tag{41}
\]

The interval is exact: the maximum occurs at $\tan ^2u=1/r$, and every
intermediate value is attained.  Put

\[
 S=k(1+p^2),\quad d=kp,\quad
 \alpha=(1-c^2)(1-\delta^2/c^2),\quad
 \beta=(1-c^2)(c^{-2}-\delta^2),
\]

and $x=S-\alpha w$, $y=S+\beta w$.  Directly from
$B=c^{-1/2}DMD^{-1}$ and $C=c^{1/2}DM^TD^{-1}$,

\[
 \begin{array}{c|c|c}
 &\operatorname{tr}(X^*X)&\det X\\ \hline
 B&x/c&d/c\\
 C&cy&cd.
 \end{array} \tag{42}
\]

For a $2\times2$ matrix whose squared singular values are
$\lambda\ge\mu$, define

\[
 U(X)=\operatorname{tr}(X^*X)
       -{|\det X|^2\over\operatorname{tr}(X^*X)}.
\]

Then $U(X)-\|X\|^2=\mu^2/(\lambda+\mu)\ge0$.  Thus (39) follows from

\[
 \left(x-{d^2\over x}\right)
 \left(y-{d^2\over y}\right)\le4. \tag{43}
\]

Set $w=kW$, $x=kX$, and $y=kY$, and parameterize (41) by

\[
 W={ (1-pr)^2\over(1+r)^2}\,\omega,qquad0\le\omega\le1.
\]

After multiplication by $XY>0$, (43) is the scalar polynomial target

\[
 k^2(X^2-p^2)(Y^2-p^2)-4XY\le0. \tag{44}
\]

Only elementary nome bounds are needed.  With $q=c^2$,
$s_0=\theta _3(q)^{-2}$, and $n^2\ge2n-1$,

\[
 \theta _3(q)\le1+{2q\over1-q^2},\qquad
 s_0\ge h:={(1-q^2)^2\over(1+2q-q^2)^2}. \tag{45}
\]

The nonnegative Taylor coefficients of $H$ also give
$s_0p\le r\le p$.  Hence it is enough to use the larger interval

\[
 r=p(1-q\lambda\eta),\quad0\le\eta\le1,qquad
 \lambda={1-h\over q}={4(1+q-q^2)\over(1+2q-q^2)^2}. \tag{46}
\]

For completeness, here is the exact polynomial submitted to the finite
certificate.  Put

\[
 L=4(1+q-q^2),\quad R=(1+2q-q^2)^2,\quad R_0=R-qL\eta,
\]
\[
 G=R+pR_0,\qquad D=R-p^2R_0,
\]

and

\[
\begin{aligned}
 X_n&=(1+p^2)G^2-(1-q)(D^2-p^2qL^2\eta^2)\omega,\\
 Y_n&=q(1+p^2)G^2+(1-q)(D^2-p^2q^3L^2\eta^2)\omega,\\
 N_x&=X_n^2-p^2G^4,\qquad
 N_y=Y_n^2-p^2q^2G^4. \tag{47}
\end{aligned}
\]

Then $X=X_n/G^2$ and $qY=Y_n/G^2$.  L23 gives
$k\le4c/(1+c^2)^2$, while trivially $k\le1$.  If the product multiplying
$k^2$ in (44) is negative, (44) is immediate; otherwise either upper bound
may replace $k$.  After clearing positive factors the two sufficient
integer polynomials are

\[
 \begin{aligned}
 P_{\rm lo}&=16N_xN_y-4(1+q)^4X_nY_nG^4,
      &&0\le c\le37/125,\\
 P_{\rm hi}&=N_xN_y-4qX_nY_nG^4,
      &&59/200\le c\le1. \tag{48}
 \end{aligned}
\]

Both are nonpositive.  The exact checker
`experiments/slice_block_product_certificate.py` regenerates the 15,953-
and 15,084-term polynomials and converts them to integer Bernstein form.
For $P_{\rm hi}$, bisecting each coordinate once after restricting $c$
gives 16 boxes, all with strictly negative Bernstein coefficients.  For
$P_{\rm lo}$, four outer boxes cover

\[
 p\ge37/125\quad\hbox{or}\quad
 \omega\le51/250\quad\hbox{or}\quad\omega\ge199/250.
\]

The remaining cube is resolved without a limiting numerical margin.  Let
$d=|\omega-1/2|$ and choose the largest of $c,p,d$.  The four exact blow-up
charts are

\[
 \begin{array}{ll}
 p=cy, & \omega=1/2+c(2z-1),\\
 c=py, & \omega=1/2+p(2z-1),\\
 c=dy,\ p=dv, & \omega=1/2+d,\\
 c=dy,\ p=dv, & \omega=1/2-d,
 \end{array}\qquad y,v,z\in[0,1]. \tag{49}
\]

After scaling the leading variable by $37/125$, all 71,145; 322,245;
1,095,633; and 1,095,633 Bernstein coefficients, respectively, are
nonpositive.  Arithmetic is over Python integers; zero coefficients remain
exactly zero.  Since a polynomial on a box is a convex combination of its
Bernstein coefficients, (48), then (44), (43), and finally (39) follow.

There is an immediate transfer consequence.  When $b=0$, formula (19)
factors as

\[
 \mathscr R_{a,0}(T)=
 \operatorname{diag}(I,C)
 \begin{bmatrix}aI&-\sqrt{1-a^2}I\\-\sqrt{1-a^2}I&-aI\end{bmatrix}
 \operatorname{diag}(I,B). \tag{50}
\]

The middle factor is orthogonal.  L24, L26, and (39) imply
$\max(1,\|B\|)\max(1,\|C\|)\le2$, so
$\|\mathscr R_{a,0}(T)\|\le2$.  The factorization with $B,C$ swapped proves
the same statement when $a=0$.  Thus both full coordinate axes of (RT) are
closed, not merely their common center.

### 6.4 A sharper nome bound and the complete zero-node face

The block-energy ridge has a second boundary mechanism which is invisible in
the estimate used in L23.  Jacobi's product gives, with $q=c^2$,

\[
 {k\over4c}=\prod_{n\ge1}
 \left({1+q^{2n}\over1+q^{2n-1}}\right)^4. \tag{51}
\]

Every omitted factor is at most one.  Retaining the first two factors and
using $0\le q\le1/4$ gives the sharper estimate

\[
 \boxed{k\le {4c\over1+4c^2}\qquad(0\le c\le1/2).} \tag{52}
\]

Indeed, after clearing denominators, (52) reduces to

\[
 (1+q)^4(1+q^3)^4
 -(1+4q)(1+q^2)^4(1+q^4)^4\ge0. \tag{53}
\]

The left side is $q^2Q(q)$.  All 24 Bernstein coefficients of $Q(x/4)$
on $0\le x\le1$ are positive.  This exact rational check is regenerated by
`experiments/slice_sharp_nome_certificate.py`.

There is also a useful exact reduction on the complete zero-node face
$p=0$.  Put

\[
 u={3a\over4-a^2},\quad v={3b\over4-b^2},\qquad
 h_a={2(1-a^2)\over4-a^2},\quad
 h_b={2(1-b^2)\over4-b^2},
\]

and $K_0=k(1+c^2)/c$.  Formula (36) shows

\[
 \|B\|^2+\|C\|^2=K_0,qquad
 \|B\|\|C\|\le {K_0\over2}, \tag{54}
\]

where the second inequality is sharp at $\tan^2\theta=1$.  The four blocks
in (34) are rank one on this face, so (BE) follows from

\[
 {h_ah_bK_0+u^2h_b^2k^2+v^2h_a^2K_0^2/4
  \over(1-uvk)^2}\le1. \tag{55}
\]

Changing the sign of exactly one of $a,b$ only increases the denominator, so
it is enough to take $u,v\ge0$.  The left side increases with $k$.  For
$c\ge1/2$, substitute $k\le1$ in (55), put $c=(1+x)/2$, and clear the
positive denominator.  The result has degree $(4,4,4)$ in $(x,a,b)$, and
all 125 tensor-product Bernstein coefficients are nonnegative (16 are
exactly zero).  Consequently

\[
 \boxed{\text{(BE), hence (RT), holds when }p=0,\ c\ge1/2.} \tag{56}
\]

For $c\le1/2$, put $x=2c$ and insert (52), now
$k\le2x/(1+x^2)$, in (55).  Clearing the positive denominator produces a
degree $(4,4,4)$ polynomial $P(x,a,b)$.  Its nonnegativity on the unit cube
has the following finite exact Bernstein cover:

* on $b\ge15/16$, split at the exact ridge $a=x$ and use
  $a=xy$ below it and $a=x+(1-x)y$ above it;
* on $b\le15/16$, use the compact boxes $x\ge1/16$ and, when
  $x\le1/16$, first $1/16\le b\le15/16$, then
  $b\le1/16$, $a\ge1/16$;
* in the remaining cube $x,a,b\le1/16$, choose a largest coordinate $t$
  and use one of
  \[
   (x,a,b)={1\over16}(t,ty,tz),\quad
   {1\over16}(ty,t,tz),\quad
   {1\over16}(ty,tz,t),
  \]
  after dividing the transformed polynomial by its exact common factor
  $t^2$.

All chart maps have rational coefficients.  The two ridge charts, the
$a\ge1/16$ chart, and the three origin charts are nonnegative without
subdivision; the compact $x\ge1/16$ chart needs eight exact dyadic
bisections and the intermediate $b$ chart needs one.  Thus the cover is
finite and uses no floating-point
sign decisions.  It includes the curved zero set $b=1$, $a=x$ and closes the
degenerate origin by continuity.  The same checker regenerates every
coefficient.  Combining the two halves gives

\[
 \boxed{\text{(BE), hence (RT), holds when }p=0,\quad0<c<1.} \tag{57}
\]

The ridge $b=1$, $a=2c$ remains the correct centered coordinate for a
neighbourhood with small positive $p$; the sharp boundary is not confined to
$a=b=0$.

## 7. Remaining target

The live theorem is now (27), equivalently (RT), with the special conformal coupling

\[
 \tan v=r\tan u,\qquad r=H(p),\qquad p=\tau_2/\tau_1. \tag{23}
\]

Arbitrary rotations and node pairs do not obey the factor-four bound.  A
scalar reduction is also insufficient: numerical interior maxima can exceed
all scalar-automorphism boundary values.  A proof must retain the
Blaschke--Potapov structure (21) while using (23).  The preferred targets are
either the polynomial Schur complement (27) or the stronger four-block energy
inequality (BE).  L47--L48 remove the complete axes $a=0$ and $b=0$;
L49--L50 add the sharper small-nome estimate and close the complete $p=0$
face.  The live transfer square is genuinely two-parameter with
$abp\ne0$, and a neighbourhood of its singular boundary must retain
$b\to1$, $a-2c\to0$, $p\to0$.  Proving
(RT) completes the elliptic \(4\times4\) slice; it does not by itself settle the
general conjecture.
