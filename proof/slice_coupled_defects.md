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
 \quad\text{for every }(a,b)\in[-1,1]^2} \tag{RT}
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
interior transfer is genuinely matrix-valued.  If \(a=\pm1\) or \(b=\pm1\),
(19) becomes a direct sum of a signed identity and a scalar automorphism of
\(F\) or \(G\); this boundary is already contained in the proved even sector.
At \(a=b=0\), the transfer is \(-T\), so L24 and L26 close that point.  The
new content of (RT) is the open interior of the parameter square.

`experiments/slice_rank_one_transfer.py` audits (19) and (22) with exact
rational/symbolic arithmetic.  Its deterministic floating-point search is
sharp on the two default rank-one/rank-one SDP cases: the squared transfer
norms \(3.136329347\) and \(2.491664377\) reproduce the primal-dual optima.
This sharpness is evidence only, not the proof of (RT).

## 6. Remaining target

The live theorem is now (RT), with the special conformal coupling

\[
 \tan v=r\tan u,\qquad r=H(p),\qquad p=\tau_2/\tau_1. \tag{23}
\]

Arbitrary rotations and node pairs do not obey the factor-four bound.  A
scalar reduction is also insufficient: numerical interior maxima can exceed
all scalar-automorphism boundary values.  A proof must retain the
Blaschke--Potapov structure (21) while using (23), for example through a
four-dimensional Schur complement or a two-node matrix-valued Pick
inequality.  Proving (RT) completes the elliptic \(4\times4\) slice; it does
not by itself settle the general conjecture.
