# The raw weighted retightening has a genuine quintic copy obstruction

## 1. Result (L308, 2026-07-26)

Retain L298's unprepared weighted correction

\[
X(c)=\sum_{k\ge1}c^{2k}R_k,\qquad
C(c)=\sum_{k\ge1}c^{2k}F_k,                      \tag{1}
\]

insert it directly into L299's moving canonical pair, and let
\({\mathfrak q}_{\rm raw}(c)\) be L306's copy successor.  L300 proves
that its cubic coefficient is a response.  That odd cancellation does
**not** persist automatically:

\[
\boxed{
[c^5]\operatorname {tr}{\mathfrak q}_{\rm raw}
\not\equiv0.}                                     \tag{2}
\]

There is an exact rational scalar-copy example.  Take

\[
S=
\begin{bmatrix}
0&3/5&-4/5\\
0&4/5& 3/5\\
0&0&0
\end{bmatrix},
\qquad V=e_1,\qquad W=e_3.                        \tag{3}
\]

Then

\[
I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad
\rho(S)=4/5,\qquad B_1=W^*S^*V=-4/5.             \tag{4}
\]

Exact partial-isometry word algebra and rational evaluation give

\[
\boxed{
[c^5]\operatorname {tr}{\mathfrak q}_{\rm raw}
=\frac{33264}{15625}>0.}                          \tag{5}
\]

For copy dimension one the homogeneous endpoint-response map is
identically zero: every response has trace zero, and its output is a
scalar.  Hence (5) is not a hidden Markov response.

L308 disproves only the **one-shot raw parity shortcut**.  It does not
contradict L300's cubic response, L301's successor transport, or the
possibility that the recursively prepared series has response-null
odd coefficients.  Instead it proves that those preparations are
load-bearing.  The live all-series object is the L307 affine
recurrence after each preceding coefficient is canceled or retained;
one may not factor the original L298 successor as though it were
already prepared.

## 2. Exact coefficient

Write the moving canonical pair as

\[
A=S+\sum_{j\ge1}c^jA_j,\qquad
D=V+\sum_{j\ge1}c^jD_j.                          \tag{6}
\]

Only the first two weighted retightening directions can enter degree
five.  Directly from L299,

\[
\begin{aligned}
{\cal N}_{5,\rm raw}=-\{&
A_1^*R_2S+S^*R_2A_1+
A_3^*R_1S+S^*R_1A_3\\
&+A_1^*R_1A_2+A_2^*R_1A_1\\
&+D_1F_2^*+F_2D_1^*
 +D_3F_1^*+F_1D_3^*\}.                           \tag{7}
\end{aligned}
\]

Purity gives

\[
\operatorname {tr}{\mathfrak q}_{\rm raw,5}
=\operatorname {tr}{\cal N}_{5,\rm raw}.         \tag{8}
\]

The only nonlocal term in (7) is the \(R_2\) cross.  If

\[
R_2-S^*R_2S=VF_2^*+F_2V^*,
\]

then L302's exact dual telescope

\[
Y_2=2(S^2+(S^*)^2),\qquad
Y_2-SY_2S^*=SA_1^*+A_1S^*
\]

replaces its trace by the finite word
\(\operatorname {tr}\{(VF_2^*+F_2V^*)Y_2\}\).
After partial-isometry reduction, (8) has eight nonzero cyclic
classes:

\[
\begin{array}{c|r}
\text{word in }a=S^*,\,s=S&\text{coefficient}\\ \hline
aa,\ ss&-1/2\\
aaaaasss,\ aaasssss&17/2\\
aaaaasssaaasss,\ aaasssaaasssss&-6\\
aaaasssaaaasss,\ aaassssaaassss&-2 .
\end{array}                                      \tag{9}
\]

Evaluating (9) on (3) gives (5) exactly.  Independently assembling the
actual matrices in (7) numerically from the rational input, including
the full Stein solutions, gives the same value to
\(4.5\times10^{-16}\).

## 3. Consequence for the recurrence

The tempting assertion

\[
\text{“all odd coefficients of the raw L298 copy energy are
responses”}
\]

is false.  The valid sequence is triangular:

1. cancel L300's cubic by an affine L307 correction;
2. carry its L301 quartic successor and L303 lower neutralization;
3. only then test the resulting quintic coefficient; and
4. formulate the arbitrary-grade recurrence for those **prepared**
successors.

This is precisely why L307's forcing \(R\) cannot be suppressed from
the all-series argument.

## 4. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_raw_quintic_obstruction.py \
  --output \
  experiments/repeated_crabb_raw_quintic_obstruction_s70226.jsonl
```

The checker constructs (7) twice: in exact reduced word algebra with
the dual telescope, and by an independent floating matrix assembly
on the rational input (3).  It also verifies (4), the exact zero
cubic cyclic class, and the agreement of the quintic coefficient with
(5).  The tracked dataset SHA-256 is

```text
0b755aec840af01085b7edf84078a3299174e025a2ef70d26360829ab8f603ec
```
