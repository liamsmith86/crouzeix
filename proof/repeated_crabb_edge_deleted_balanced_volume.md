# The delayed edge deletion removes every metric square root

## 1. Result (L250, 2026-07-25)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\qquad EF=0,
\]

and L219's boundary metric

\[
R(c)=I-\sum_{j\geq1}\frac{q^j}{1+q^j}(S^*)^jES^j
 \sum_{j\geq1}q^jS^jF(S^*)^j,\qquad q=c^2.       \tag{1}
\]

Assume a complete delay

\[
B_1=\cdots=B_{k-1}=0,\qquad
B_j=W^*(S^*)^jV.                                  \tag{2}
\]

Let \(X_k=[q^k]R\), delete that **entire** coefficient, and put

\[
R^\circ=R-q^kX_k,\qquad P=I-F.
\]

Then the final defect is metrically isolated through the active
order:

\[
\boxed{
R^\circ
=F\oplus PR^\circ P+O(q^{k+1}).}                  \tag{3}
\]

Let \(A(c)\) be the balanced direct ellipse map and define

\[
\widetilde A=(R^\circ)^{1/2}
 A(R^\circ)^{-1/2}.
\]

Use L248's notation

\[
M=I-\widetilde A^*F\widetilde A,\qquad
H_R=I-\widetilde A^*\widetilde A.
\]

In balanced coordinates put instead

\[
N=R^\circ-A^*FA,\qquad
H=R^\circ-A^*R^\circ A.                           \tag{4}
\]

Then every coefficient through the live face can be calculated
without a metric square root:

\[
\boxed{
[c^d]\operatorname {tr}(M^{-1}H_R)
=[c^d]\operatorname {tr}(N^{-1}H),
\qquad 0\leq d\leq2k.}                            \tag{5}
\]

Equivalently, L248's edge-deleted volume has the balanced determinant
form

\[
\boxed{
\mathcal V_{\widetilde A}
=2^p
\frac{
\det\{R^\circ-A^*(F+\tfrac12PR^\circ P)A\}
}{
\det\{R^\circ-A^*FA\}
}
+O(c^{2k+2}),}                                    \tag{6}
\]

where \(p=\operatorname {rank}P\).

By Sylvester's determinant identity, (6) also has the output-oriented
form

\[
\boxed{
\mathcal V_{\widetilde A}
=2^p
\frac{\det(I-D^\circ Z)}
     {\det(I-FZ)}
+O(c^{2k+2}),\qquad
\begin{cases}
D^\circ=F+\tfrac12PR^\circ P,\\
Z=A(R^\circ)^{-1}A^*.
\end{cases}}                                      \tag{7}
\]

This is the preferred interface with L244: its coisometric identity
is exactly an equality of the form \(A R^{-1}A^*=R^{-1}\) on the
retained half-line.  L243--L245 can therefore be inserted into \(Z\)
before the finite final port is eliminated.

L264 subsequently computes the response to a retained active metric
coefficient \(X=PXP\) as \(-\operatorname {tr}(EX)\).  Thus that
block is not arbitrary.  The right-half-line part is nevertheless
free because every positive right wandering orbit is orthogonal to
\(E\).  One may therefore restore
\(c^{2k}P([c^{2k}]R_{\rm right})P\), while leaving the active left
orbit, \(F\)-corner, and cross row deleted.  This port-isolated
right-half-line gauge is the preferred L244 interface.

Thus A194 is now the balanced-coordinate scalar identity

\[
\boxed{
[c^{2k}]
\left\{-m+\operatorname {tr}(N^{-1}H)\right\}
=4\|B_k\|_F^2,}                                   \tag{8}
\]

with all earlier coefficients zero.  Equation (8) remains open.
L250 removes the analytic square roots and places the target directly
in the already-developed right Stein slack \(H\), while retaining
the exact final-row whitening \(N^{-1}\).  It does not separate the
large terms hidden inside that inverse.

## 2. Final-port isolation

Write

\[
E_j=(S^*)^jES^j,\qquad
F_j=S^jF(S^*)^j.
\]

Since \(FS=0\),

\[
FF_j=F_jF=0\qquad(j\geq1).                        \tag{9}
\]

Moreover,

\[
\begin{aligned}
W^*E_jW&=B_jB_j^*,\\
E_jW&=(S^*)^jVB_j^*.
\end{aligned}                                     \tag{10}
\]

The scalar weight in the first sum of (1) is

\[
\frac{q^j}{1+q^j}
=\sum_{\ell\geq1}(-1)^{\ell-1}q^{j\ell}.          \tag{11}
\]

Hence the coefficient of \(q^h\) in either \(FRF-F\) or \(PRF\)
uses only \(B_j\) with \(j\mid h\).  By (2), both blocks vanish for
\(h<k\).  At \(h=k\), every proper-divisor term still vanishes and
the only possible block is part of the full matrix coefficient
\(X_k\).  Deleting \(X_k\) removes it, proving (3).  The next possible
term is \(q^{k+1}=c^{2k+2}\).

This argument is why the **whole** coefficient must be deleted.
Removing only its final corner would leave a cross row and would not
justify (5).

## 3. Balanced initial-defect mass

Replace \(R^\circ\), modulo \(O(c^{2k+2})\), by the block-diagonal
series in (3).  Functional calculus at the identity gives

\[
(R^\circ)^{1/2}F
=F(R^\circ)^{1/2}=F.                              \tag{12}
\]

Therefore

\[
\begin{aligned}
M
&=(R^\circ)^{-1/2}
 \{R^\circ-A^*FA\}
 (R^\circ)^{-1/2},\\
H_R
&=(R^\circ)^{-1/2}
 \{R^\circ-A^*R^\circ A\}
 (R^\circ)^{-1/2}.                                \tag{13}
\end{aligned}
\]

Cyclicity of the trace in (13) gives

\[
\operatorname {tr}(M^{-1}H_R)
=\operatorname {tr}(N^{-1}H),
\]

which proves (5) through degree \(2k\).  Notice that \(H\) is exactly
the balanced right Stein slack already used in L225--L228.  The new
ingredient is the row denominator \(N\), which resums the terms that
would otherwise appear as the final Schur square.

For the determinant, L248 gives

\[
\mathcal V_{\widetilde A}
=2^p
\frac{\det(I-\widetilde A^*D_F\widetilde A)}
     {\det(I-\widetilde A^*F\widetilde A)},
\qquad D_F=F+\frac12P.
\]

Equation (12) yields

\[
(R^\circ)^{1/2}D_F(R^\circ)^{1/2}
=F+\frac12PR^\circ P.
\]

Conjugate both determinant arguments by
\((R^\circ)^{1/2}\).  Their common determinant factor cancels and
gives (6).

Factoring \(R^\circ\) from each determinant in (6) and applying
\(\det(I-UV)=\det(I-VU)\) gives (7).

## 4. Consequence and guardrail

The normalized initial defect always has \(m\) eigenvalues equal to
one because L248's whitened main row has \(m\)-dimensional kernel.
The remaining small eigenvalues need not be positive; generic delayed
audits are indefinite even though their sum has coefficient
\(+4\|B_k\|_F^2\).  Therefore L250 supports a **trace-energy**
argument, not a positive leakage-Gram claim.

The safe continuation of A194 is now:

1. insert L243--L245's physical zero/one-reflection lift into
   \(A\), \(H\), and \(N\);
2. use L244's coisometry inside the product \(N^{-1}H\);
3. use matrix-inner autocorrelation to cancel the later
   \(\mathcal B^\sharp(x)\) rows; and
4. retain the doubled remote cell only after those cancellations.

## 5. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_edge_deleted_balanced_volume.py \
  --output \
  experiments/repeated_crabb_edge_deleted_balanced_volume_s70224.jsonl
```

The checker uses generic complete delays through grade five.  It
verifies the final metric block and cross-row vanishing in (3),
constructs the metric square root as a formal series, and compares
(5)--(7) coefficientwise with the literal normalized mass, both
determinant orientations, and the final-defect Schur trace.  It also
repeats the open \(+4\) coefficient as a falsification audit.  The
algebra above proves (3)--(7); the observed value in (8) is not used
as proof.

The tracked dataset SHA-256 is

```text
6d41542e06d6bc7e6e50ea91e57a37ad91f887b1ab91a40a3c2cd4079f69c4fb
```
