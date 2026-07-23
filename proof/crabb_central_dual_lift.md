# The exact central Blaschke metric lift (2026-07-23)

## 1. Result

Use L126's central family

\[
 T=T_k(a,c),\qquad U=B_{k,c}(T),\qquad
 V=\operatorname{span}\{e_0,e_k,e_{2k}\},
\]

in physical coordinates.  Put `r=c^k`, and let `O:C^3 -> V` be the
outer isometry.  L126 and L129 give

\[
 U=OU_3O^*\oplus U_{\rm in},\qquad
 t_*(U)=t_*(U_3),\qquad U_3=T_1(a,r)                 \tag{1}
\]

locally at `(a,c)=(0,0)`.

There is also a local upper lift with no condition-number cost:

\[
\boxed{
 t_*(T_k(a,c))=t_*(T_1(a,c^k))
 \quad (|a|+c\hbox{ sufficiently small}).}           \tag{2}
\]

The neighbourhood in (2) may depend on `k`.  This is an exact local
identity for every fixed central degree, not a uniform-in-dimension
neighbourhood theorem.

The proof constructs the full rank-one Stein metric explicitly.  Its
key input is a finite-Blaschke trace identity whose multiplier is the
outer critical factor of `B`.

## 2. The outer critical factor

First consider any real finite Blaschke product of degree `k`,

\[
 B(z)={N(z)\over D(z)},\qquad
 N(z)=\prod_{\nu=1}^k(z-b_\nu),\qquad
 D(z)=\prod_{\nu=1}^k(1-b_\nu z),                    \tag{3}
\]

where `b_nu in (-1,1)`.  A harmless real unit phase has been omitted.
Put

\[
 W=N'D-ND'.                                          \tag{4}
\]

The critical points of a finite Blaschke product occur in
reflection pairs across the unit circle, and none lies on the circle.
Choose the factor `Q` containing the critical points outside the disk.
With roots at infinity included projectively,

\[
 W=\kappa Q Q^\#,\qquad
 Q^\#(z)=z^{k-1}\overline{Q(1/\overline z)}.          \tag{5}
\]

Thus `deg Q<=k-1`, and `Q` has no zero in the closed disk.  Define

\[
 F={Q\over D}.                                       \tag{6}
\]

Then `F` is a zero-free member of the model space `K_B`.  Recall that
every member of that space has the form

\[
 f={p\over D},\qquad \deg p<k.                       \tag{7}
\]

## 3. Constant unweighted fiber trace

For `u` away from the critical values, let

\[
 {\cal T}_u(f)
 ={1\over k}\sum_{B(t)=u}{f(t)\over F(t)}
 ={1\over k}\sum_{B(t)=u}{p(t)\over Q(t)}.           \tag{8}
\]

Then

\[
\boxed{
 {\cal T}_u(f)
 ={ \langle f,F\rangle_{H^2}\over\|F\|_{H^2}^2},
 \qquad f\in K_B,}                                   \tag{9}
\]

so the left side is independent of the fiber value `u`.

Here is a direct proof.  Set `R_u=N-uD`.  The sum in (8) is the sum of
the residues at the zeros of `R_u` of

\[
 {p(z)\over Q(z)}{R_u'(z)\over R_u(z)}\,dz.           \tag{10}
\]

Every other finite pole is an outside critical point `q` of `B`.
Since

\[
 {R_u'\over R_u}={D'\over D}+{B'\over B-u},           \tag{11}
\]

the second term has a zero of at least the same order as `Q` at `q`.
The residue there therefore comes only from the first term and is
independent of `u`.  The same argument in the coordinate `1/z`
handles a critical point at infinity.  The residue theorem proves
that (8) is independent of `u`, with multiplicities and critical
values obtained by continuation.

It remains to identify the constant.  On the unit circle,

\[
 \overline{F(z)}={zQ^\#(z)\over N(z)}.
\]

Using (5),

\[
\begin{aligned}
 \langle f,F\rangle_{H^2}
 &= {1\over2\pi i}\int_{|z|=1}
    {p(z)Q^\#(z)\over D(z)N(z)}\,dz\\
 &= {1\over\kappa\,2\pi i}\int_{|z|=1}
    {p(z)\over Q(z)}
    \left({N'(z)\over N(z)}-{D'(z)\over D(z)}\right)dz\\
 &= {1\over\kappa}\sum_{N(t)=0}{p(t)\over Q(t)}.      \tag{12}
\end{aligned}
\]

The last sum is the fiber trace at `u=0`.  Taking `f=F` in (12)
shows `kappa ||F||^2=k`, after the harmless normalization in (5).
Equation (9) follows.

This argument is basis-free.  In an orthonormal real model-space
basis write

\[
 F=\sum_{j=0}^{k-1}\gamma_jf_j,\qquad
 \beta=\|F\|^2=\sum_j\gamma_j^2.                     \tag{13}
\]

Then (9) says

\[
 {1\over k}\sum_{B(t)=u}{f_j(t)\over F(t)}
 ={\gamma_j\over\beta}.                              \tag{14}
\]

## 4. Central subcritical conditional expectation

Return to coefficient coordinates for L126.  Write

\[
 S=S_k(a,c)=S_k(0,c)+a\,u v^*,
\]

\[
 u=2(e_0-e_{2k}),\qquad
 v=e_{k+1}-ce_{k-1}.                                 \tag{15}
\]

Let `E:C^3 -> C^(2k+1)` embed the outer coordinates.  For the Dickson
polynomials

\[
 P_0=2,\qquad P_1=z,\qquad
 P_m=zP_{m-1}-cP_{m-2},
\]

one has the exact subcritical identities

\[
\boxed{
 E^*P_m(S)E=0,\qquad 1\leq m<k.}                     \tag{16}
\]

For the unperturbed `S_0=S_k(0,c)`, (16) is immediate from L126's
finite Dickson polynomial model: a nonzero grade `m<k` cannot return
to an outer residue.

The amplitude does not change this compression.  The same complete
column formulas give

\[
\begin{aligned}
 v^*P_j(S_0)E&=0 &&(0\leq j\leq k-2),\\
 E^*P_j(S_0)u&=0 &&(1\leq j\leq k-1).                \tag{17}
\end{aligned}
\]

Expand the Dickson recurrence at `S_0+a uv^*`.  Every term containing
an amplitude insertion has a rightmost insertion.  The Dickson word
between that insertion and `E` has grade at most `m-1<=k-2`, so its
outer compression contains the first zero in (17).  Terms without an
insertion have the already-zero unperturbed compression.  This proves
(16), including every power of `a`.

The polynomials `1,P_1,...,P_(k-1)` form a basis modulo `P_k`.
Consequently polynomial division and (16) give

\[
\boxed{
 E^*g(S)E
 =({\cal E}_P g)(P_k(S)|_V),\qquad
 ({\cal E}_P g)(w)={1\over k}\sum_{P_k(z)=w}g(z).}    \tag{18}
\]

It is enough to prove (18) for polynomials and then use analytic
approximation.  The same root-of-unity calculation as L120 shows
that `E_P P_m=0` for `1<=m<k`, so the polynomial quotient on both
sides is identical.

The proper-map identity

\[
 B_{k,c}\circ\phi_c=\phi_{c^k}\circ P_k              \tag{19}
\]

identifies the `P_k` fibers in (18) with the `B_(k,c)` fibers in
(8).  Passing from coefficient to physical coordinates does not
alter a scalar compression because L126's coordinate Gramian reduces
`V`.  Therefore (9) and (18) prove

\[
\boxed{
 O^*\left({f\over F}\right)(T)O
 ={\langle f,F\rangle\over\|F\|^2}I_3,
 \qquad f\in K_B.}                                   \tag{20}
\]

This is the missing root-of-unity identity.

## 5. Explicit defect lift and cross cancellation

Define the injective three-column map

\[
 J=F(T)^{-*}O.                                       \tag{21}
\]

Equation (20), in the basis (13), becomes

\[
 O^*f_j(T)^*J={\gamma_j\over\beta}I_3.               \tag{22}
\]

Let `v_3` be the rank-one Stein defect vector of the local optimal
size-three metric, and put

\[
 q=Jv_3,\qquad q_j=f_j(T)^*q.                        \tag{23}
\]

Write

\[
 q_j={\gamma_j\over\beta}Ov_3+w_j,\qquad
 w_j\in V^\perp.                                     \tag{24}
\]

Because `F=sum gamma_j f_j`,

\[
\sum_j\gamma_jq_j
 =F(T)^*F(T)^{-*}Ov_3
 =Ov_3.                                              \tag{25}
\]

The outer part of the left side is also

\[
 {1\over\beta}\sum_j\gamma_j^2Ov_3=Ov_3.
\]

Subtracting proves the formerly open identity

\[
\boxed{\sum_j\gamma_jw_j=0.}                         \tag{26}
\]

Equivalently,

\[
 \Pi_V\left(\sum_jq_jq_j^*\right)(I-\Pi_V)=0,         \tag{27}
\]

because the outer component in (24) is always collinear with `v_3`.

There is also an exact dual interpretation.  For every `3 x 3`
matrix `Z`, L127's dual transfer obeys

\[
\begin{aligned}
 {\cal R}_{B,T}(OZO^*)J
 &=\sum_j f_j(T)OZO^*f_j(T)^*J\\
 &={1\over\beta}F(T)OZ.                              \tag{28}
\end{aligned}
\]

Thus, if `Z>=0` has kernel `v_3`, its positive lift annihilates the
explicit vector `q=Jv_3`.  No numerical corank inference is needed.

## 6. Equality of the two similarity optima

L127's Gramian composition can be written

\[
 {\cal P}_T(q)=\sum_{j=0}^{k-1}{\cal P}_U(q_j).       \tag{29}
\]

The operator `U=U_3+U_in` reduces `V`.  Equation (27) therefore makes
the metric in (29) reduce `V` as well.  Its outer compression is

\[
 O^*{\cal P}_T(q)O
 ={1\over\beta}{\cal P}_{U_3}(v_3).                  \tag{30}
\]

At `(a,c)=(0,0)`, `B(z)=z^k`, `F=1`, and `q=e_0`.  The full metric is

\[
 {\cal P}_{C_{2k+1}}(e_0)
 =\operatorname{diag}(1,2,\ldots,2,4).               \tag{31}
\]

Its outer levels are `1,2,4`, while every remaining level is the
strictly interior value two.  The optimal size-three defect and all
entries in (29) vary continuously on L118's analytic branch.
Therefore, after shrinking a neighbourhood for the fixed `k`, every
inner metric level remains strictly between the two active outer
endpoints.  Equations (29)--(30) give

\[
 t_*(T)\leq t_*(U_3).                                \tag{32}
\]

Conversely, inner-function contraction gives `t_*(T)>=t_*(U)`, and
L129 gives `t_*(U)=t_*(U_3)` locally.  This proves (2).

As a corollary, the complete central amplitude Hessian in size
`2k+1` is the size-three Hessian at `r=c^k`.  In particular its first
nonzero joint term is

\[
 -64a^2c^{2k}.                                       \tag{33}
\]

## 7. Deterministic regeneration

Run

```bash
.venv/bin/python -u \
  experiments/crabb_central_dickson_descent.py \
  --output experiments/crabb_central_dickson_descent_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_central_dual_lift.py \
  --output experiments/crabb_central_dual_lift_s70223.jsonl
```

The exact checker verifies (16)--(17) symbolically through `k=12`.
The independent floating checker factors the Blaschke critical
Wronskian, constructs `F` and `J` directly, verifies (20), (26), and
(28), and regenerates (30) and the final condition on degrees two
through five.  The floating checker is a regression for the analytic
proof, not its replacement.
