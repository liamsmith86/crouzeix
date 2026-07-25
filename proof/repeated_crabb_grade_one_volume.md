# The edge-deleted volume has the universal grade-one response

## 1. Result (L256, 2026-07-25)

Let \(S\) be an arbitrary finite-dimensional partial isometry with
initial and final defects

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad \operatorname{rank}E
=\operatorname{rank}F=m.
\]

Write the balanced direct map as \(A(c)\).  At grade one, delete the
full active \(c^2\) coefficient from L219's boundary metric and call
the resulting metric \(R^\circ\).  L250's normalized mass is

\[
\mu^\circ(c)
=-m+\operatorname{tr}\!\left(N(c)^{-1}H(c)\right),\qquad
N=R^\circ-A^*FA,\quad H=R^\circ-A^*R^\circ A .
\]

Then, without any scalar, commutative, or monomial restriction,

\[
\boxed{[c^2]\mu^\circ(c)=4\|B_1\|_F^2.}             \tag{1}
\]

The constant and linear coefficients vanish.  Thus L250's open
all-delay identity is proved at its relative grade-one base.  The
remaining all-grade theorem is an **associated delay-covariance**
statement: removing \(k-1\) clean layers must transport the first
active scalar coefficient to (1).  L235 already shows that this cannot
be strengthened to equality of the complete tail series.

## 2. Exact low-order direct map

Put

\[
J=(I+F)S^*(I+E).
\]

The theta/ODE coefficients in L125 give

\[
A=S+cA_1+c^2A_2+O(c^3),                            \tag{2}
\]

where

\[
\begin{aligned}
A_1&=J-S^3,\\
A_2&=2S-(S^2J+SJS+JS^2)+S^5.                     \tag{3}
\end{aligned}
\]

The term \(2S\) in \(A_2\) is load-bearing.  Omitting it would repeat
the false Newton-bottom-edge shortcut ruled out in L235.

At grade one, deleting the active metric face leaves

\[
R^\circ=I+O(c^4).                                  \tag{4}
\]

Because \(FS=0\) and \(S^*F=0\), (2)--(4) imply

\[
N=I-c^2A_1^*FA_1+O(c^3),\qquad
N^{-1}=I+c^2A_1^*FA_1+O(c^3).                     \tag{5}
\]

For \(H=\sum c^jH_j\),

\[
\begin{aligned}
H_0&=E,\\
H_1&=-(S^*A_1+A_1^*S),\\
H_2&=-(S^*A_2+A_2^*S+A_1^*A_1).                  \tag{6}
\end{aligned}
\]

Consequently the only degree-two scalar face is

\[
\operatorname{tr}H_2+
\operatorname{tr}(A_1^*FA_1E).                   \tag{7}
\]

## 3. Cyclic reduction

Let \(Q_j=(S^*)^jS^j\), and write \(\equiv_{\rm cyc}\) for equality
after taking the trace and imposing only the partial-isometry and
orthogonal-defect relations.  Direct expansion of (3), followed by
exact cyclic reduction, gives

\[
\begin{aligned}
H_1&\equiv_{\rm cyc}0,\\
H_2&\equiv_{\rm cyc}36Q_1-12Q_3-24I,\\
A_1^*FA_1E&\equiv_{\rm cyc}16Q_3-48Q_1+32I.
\end{aligned}                                     \tag{8}
\]

The two degree-two pieces are not separately the desired energy.
Their sum is

\[
4(Q_3-3Q_1+2I).                                   \tag{9}
\]

L252's exact radial telescope at \(j=1,3\) reads

\[
\operatorname{tr}Q_1=n-m,\qquad
\operatorname{tr}Q_3=n-3m+\|B_1\|_F^2.            \tag{10}
\]

Taking the trace of (9) and using (10) proves (1).  This also explains
why the raw Stein contribution and final-row whitening must remain
paired: only their cyclic sum collapses to the transfer energy.

## 4. Exact audit

The identities in (8) were regenerated with the exact rational word
engine used by L228 and L252.  Independently, L252's grade-one checker
reconstructs the complete L250 mass and verifies that its cyclic face
equals (9), with zero radial and transfer differences:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_cyclic_radial_volume.py \
  --maximum-grade 1 \
  --output /tmp/repeated_crabb_grade_one_volume.jsonl
```

The one-record audit has SHA-256

```text
87940887a9d67e82b40e31668f4c059984ff83df3130aa533b1ef029449171d6
```

No floating arithmetic enters this verification.  This proves the
universal relative base case only; it does not prove the required
associated coefficient covariance through an arbitrary clean delay.
