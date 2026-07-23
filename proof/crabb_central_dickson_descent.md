# Central equality directions descend exactly to size three (2026-07-23)

## 1. The algebraic theorem

Let `k>=1`, put `p=2k+1`, and work in L123's coefficient coordinates.
Let `C=C_p` have first superdiagonal weight two and all remaining
superdiagonal weights one.  With `J` denoting coordinate reversal, set

\[
 S_k(c)=C+cJCJ.                                        \tag{1}
\]

The phase-one equality coefficient at the central offset `k` produces
the exact pencil

\[
 S_k(a,c)=S_k(c)+aE_k(c),
\]

\[
 E_k(c)=2(e_0-e_{2k})
       (e_{k+1}^*-ce_{k-1}^*).                        \tag{2}
\]

Let the Dickson polynomials be

\[
 P_0(z)=2,\qquad P_1(z)=z,\qquad
 P_n(z)=zP_{n-1}(z)-cP_{n-2}(z).                      \tag{3}
\]

Put `V=span{e_0,e_k,e_(2k)}` and let `Pi` be its coordinate
compression.  Then `V` reduces both the polynomial and its first
amplitude derivative:

\[
\boxed{
\begin{aligned}
\Pi P_k(S_k(c))\Pi
 &=
 \begin{pmatrix}
 0&2&0\\ c^k&0&1\\0&2c^k&0
 \end{pmatrix},\\[4pt]
\Pi DP_k(S_k(c))[E_k(c)]\Pi
 &=
 \begin{pmatrix}
 -2c^k&0&2\\0&0&0\\2c^k&0&-2
 \end{pmatrix},
\end{aligned}}                                        \tag{4}
\]

and all four cross blocks between `V` and `V^\perp` vanish.

The first matrix in (4) is exactly the size-three Crabb ellipse with
parameter

\[
 r=c^k.
\]

The second is exactly its central phase-palindromic equality tangent.
Thus the central offset-`k` amplitude is an exact size-three outer
event after one degree-`k` descent.

## 2. Exact coefficient proof

The recurrence (3) has the explicit form

\[
 P_n(z)=\sum_{j=0}^{\lfloor n/2\rfloor}
 (-1)^j{n\over n-j}{n-j\choose j}c^jz^{n-2j}.        \tag{5}
\]

There is a convenient finite polynomial model for this cancellation.
Over `Q(c)` set

\[
 q_0=1,\qquad q_j=c^{-j}P_j\quad(1\leq j<2k),\qquad
 q_{2k}={P_{2k}\over2c^{2k}},
\]

and impose the terminal relation

\[
 P_{2k+1}=cP_{2k-1}.                                  \tag{6}
\]

Multiplication by `z` in the ordered basis `(q_0,...,q_(2k))` is
exactly `S_k(c)`: the factor two at the two terminal edges comes from
`P_0=2` and the normalization of `q_(2k)`.  The Dickson product law

\[
 P_mP_n=P_{m+n}+c^nP_{m-n}\quad(m>n),\qquad
 P_n^2=P_{2n}+2c^n                                   \tag{7}
\]

and (6) now give the complete column formula

\[
P_k(S_k)e_j=
\begin{cases}
c^ke_k,&j=0,\\
c^{k-j}e_{k-j}+c^ke_{k+j},&0<j<k,\\
2e_0+2c^ke_{2k},&j=k,\\
e_{j-k}+c^{2k-j}e_{3k-j},&k<j<2k,\\
e_k,&j=2k.
\end{cases}                                           \tag{8}
\]

Although the model used powers of `c^{-1}`, (8) is polynomial in `c`,
so the identity extends to `c=0`.  In particular,

\[
\begin{aligned}
P_k(S_k)e_0&=c^ke_k,\\
P_k(S_k)e_k&=2e_0+2c^ke_{2k},\\
P_k(S_k)e_{2k}&=e_k.                                  \tag{9}
\end{aligned}
\]

The other cases of (8) contain no outer row.  Thus (9) and the full
column formula prove both cross-block assertions for `P_k(S_k)` and
give the first matrix in (4).

For the derivative, insert (2) into the differentiated recurrence

\[
\begin{aligned}
D_0&=0,\qquad D_1=E_k,\\
D_n&=E_kP_{n-1}(S_k)+S_kD_{n-1}-cD_{n-2}.             \tag{10}
\end{aligned}
\]

Using (7) before and after the one `E_k` insertion gives the following
full column formula at `n=k`:

\[
D_ke_j=
\begin{cases}
-2c^ke_j+2c^{k-j}e_{2k-j},&0\leq j<k,\\
0,&j=k,\\
2c^{2k-j}e_{2k-j}-2e_j,&k<j\leq2k.
\end{cases}                                           \tag{11}
\]

For completeness, (11) can also be checked without the polynomial
model: substitute (5) in

\[
 DP_k(S)[E]=\sum_\ell d_{k,\ell}c^\ell
 \sum_{r=0}^{k-2\ell-1}S^rES^{k-2\ell-1-r},
\]

where `d_(k,ell)` is the coefficient in (5), and pair every path at
its first reversal.  The unpaired monotone paths give exactly (11).
Equivalently, (10), the three-term product law (7), and the two
terminal relations verify (11) by induction.  Its three outer columns
are

\[
\begin{aligned}
D_ke_0&=-2c^ke_0+2c^ke_{2k},\\
D_ke_k&=0,\\
D_ke_{2k}&=2e_0-2e_{2k}.                              \tag{12}
\end{aligned}
\]

All remaining columns in (11) have only inner rows.  Hence (11)--(12)
prove the second identity and both derivative cross-block zeros in
(4).

## 3. Conformal consequence

The degree-`k` Chebyshev--Blaschke product satisfies the proper-map
identity

\[
 B_{k,c}\circ\phi_c=\phi_{c^k}\circ P_k.              \tag{9}
\]

Therefore

\[
 B_{k,c}(T_k(a,c))
 =\phi_{c^k}(P_k(S_k(a,c))).                          \tag{10}
\]

To first order in `a`, (4) makes the outer block on the right exactly
the size-three equality/ellipse path with parameter `c^k`; the
remaining fibers are reducing at that order.

This explains two otherwise surprising A85 facts:

1. the first central Hessian event occurs at `a²c^(2k)`; and
2. its coefficient is the size-three coefficient `-64`.

The exact formal audit gives the stronger finite observation

\[
 H_{5,2}(c)=H_{3,1}(c^2)\pmod {c^{14}},                \tag{11}
\]

where `H_(p,k)` denotes the optimized amplitude Hessian.  Equation
(11) is evidence for a metric-fiber self-similarity, not part of L126.

## 4. What remains for A84

L126 proves the all-size polynomial descent, but (9) alone has the
wrong direction for an upper similarity bound: a metric making `T` a
contraction also makes `B(T)` a contraction, not conversely.

The next lemma must lift the active size-three rank-one metric through
the inactive inner fibers.  At the L117 axis their generalized metric
levels lie strictly between the two active endpoints.  The expected
Schur/implicit-function argument is:

1. solve the inner-fiber Stein equations using their strict endpoint
   gaps;
2. show that their feedback into the active block costs more than
   `c^(2k)`;
3. retain the size-three endpoint coefficient `-64`; and
4. localize a noncentral offset to its central `2k+1` path window.

Distinct offset grades must then be polarized.  L126 does not yet
prove that metric-lifting statement or the uniform Newton remainder.

## 5. Exact regeneration

Run

```bash
.venv/bin/python -u \
  experiments/crabb_central_dickson_descent.py \
  --output experiments/crabb_central_dickson_descent_s70223.jsonl
```

The checker constructs the polynomial and derivative recurrences
symbolically and verifies both outer identities and all cross-block
zeros through `k=12`.
