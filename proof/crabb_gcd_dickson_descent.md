# Common-divisor Dickson descent for one Crabb grade (2026-07-23)

## 1. Exact reduction

Let

\[
L=dq,\qquad k=ds,\qquad
1\le s\le\lfloor q/2\rfloor,
\]

and take the phase-palindromic equality coefficient at grades
`k` and `L-k`.  In the coefficient coordinates of L123, write its
ellipse pencil as `S_(L,k)(a,c)`.  For the degree-`d` Dickson
polynomial

\[
P_0=2,\qquad P_1=z,\qquad
P_m=zP_{m-1}-cP_{m-2},
\]

put

\[
V_d=\operatorname{span}\{e_0,e_d,e_{2d},\ldots,e_{qd}\}.
\]

Then `V_d` and its complement reduce the complete polynomial pencil,
and

\[
\boxed{
P_d(S_{L,k}(a,c))|_{V_d}
=S_{q,s}(a,c^d).
}                                                     \tag{1}
\]

The equality includes every amplitude power; in fact the left side is
affine in `a`.  The coordinate Gramian reduces the same space and

\[
\boxed{
K_{L,k}(a)|_{V_d}=K_{q,s}(a).
}                                                     \tag{2}
\]

Finally,

\[
\boxed{
\Pi_{V_d}P_m(S_{L,k}(a,c))\Pi_{V_d}=0,
\qquad1\le m<d.
}                                                     \tag{3}
\]

Taking `d=gcd(L,k)` makes `gcd(q,s)=1`.  Thus every one-grade
localization problem reduces exactly to a coprime length-grade pair.
L132 is the special case `d=k`, for which `s=1`.

## 2. Residue proof

Use the finite Dickson polynomial model from L126 and L132.  Impose
the terminal relation

\[
P_{L+1}=cP_{L-1}
\]

and use

\[
P_mP_n=P_{m+n}+c^nP_{m-n}\quad(m>n),\qquad
P_n^2=P_{2n}+2c^n.                                  \tag{4}
\]

Multiplication by `P_d` changes a polynomial grade only by `+d` or
`-d`, including at a terminal fold.  It therefore preserves every
residue class modulo `d`.  On residue zero, divide all indices by
`d`.  The axis part becomes exactly

\[
C_{q+1}+c^dJC_{q+1}J.
\]

The equality tangent has coefficient offsets `ds` and `d(q-s)`.
After the same division its residue-zero part is precisely the
grade-`s` tangent in length `q`.  This proves the linear part of (1).

There are no higher amplitude powers.  Between two tangent insertions
in a degree-`d` Dickson word there are at most `d-2` axis factors.
The tangent row is supported one coordinate from indices divisible by
`d`, while its column vector is supported at the two terminal
endpoints.  A path of length at most `d-2` cannot return between those
supports.  Hence every word with two insertions is zero.  The same
distance argument before the first possible return gives (3).

The Toeplitz coordinate Gramian has offsets `ds` and `d(q-s)`.
It preserves residues modulo `d`, and its residue-zero compression
has offsets `s` and `q-s` with exactly the endpoint multiplicities of
the reduced chart.  This proves (2).

## 3. Rank-one metric consequence

Let `T_(L,k)(a,c)=phi_c(S_(L,k)(a,c))`.  The proper-map identity

\[
B_{d,c}\circ\phi_c=\phi_{c^d}\circ P_d
\]

and (1) show that the active descended block is the complete
`T_(q,s)(a,c^d)` family.  At the apex, every other residue class is
an unweighted size-`q` nilpotent shift with metric levels equal to two,
strictly between the active endpoint levels one and four.  It remains
inactive locally for each fixed triple `(d,q,s)`.

Equations (2)--(3) are exactly the hypotheses used in L132's outer
critical-factor reconstruction.  Repeating that proof gives

\[
\boxed{
t_*(T_{q,s}(a,c^d))
\le t_*(T_{L,k}(a,c))
\le \Gamma_{q+1}(T_{q,s}(a,c^d)),
}                                                     \tag{5}
\]

where `Gamma_(q+1)` is L118's locally optimized rank-one envelope on
the reduced family.  The reconstruction is complex-linear in the
outer defect, so (5) permits an arbitrary one-pair phase
`(u,conjugate(u))`.

As in L132, (5) is a sandwich, not an assertion that the full and
reduced similarity optima are equal away from the cases where outer
rank-one optimality is already known.

If the reduced coprime grade has upper face

\[
-64|u|^2(c^d)^{2s},
\]

then (5) transfers it without condition cost to
`-64|u|²c^(2k)`.  Consequently the remaining nondivisor diagonal
problem may be restricted to

\[
\boxed{\gcd(L,k)=1.}                                 \tag{6}
\]

## 4. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_gcd_dickson_descent.py \
  --output experiments/crabb_gcd_dickson_descent_s70223.jsonl
```

The checker audits every coprime reduced pair through
`2<=d<=5`, `2<=q<=8`.  It verifies (1)--(3), exact amplitude
linearity, both polynomial and Gramian cross-block zeros, and every
inactive apex shift using exact SymPy polynomial arithmetic.
