# Principal-face completed-square target (2026-07-23)

## 1. New exact grade-four audit

The guarded rational Hessian engine now reaches the first previously
untested coprime noncentral grade

\[
(L,k)=(9,4),\qquad \gcd(L,k)=1.
\]

It gives

\[
[a^2]\Gamma_{9,4}(a,c)=-64c^8+O(c^9)
\]

exactly: every lower coefficient is zero, and the two quadratic
metric endpoint coefficients are `(48,128)`, so the gauge-invariant
combination is

\[
128-4(48)=-64.                                      \tag{1}
\]

This extends the finite exact frontier from grades at most three to
grade four.  It is evidence for, not a proof of, all-size
principal-face locality.

## 2. Transported defect and the one-coordinate correction

Let `d_c` be L117's axis defect and let `U(u)` be the lower
triangular Toeplitz coefficient operator from L131.  The natural
transported disk defect tangent is

\[
x_{\rm tr}=2U(u)d_c.                                \tag{2}
\]

In `(L,k)=(9,4)`, the exact optimized defect agrees with (2) in every
coefficient below `c^4`.  At weight four there is exactly one
difference:

\[
[c^4](x_*-x_{\rm tr})=-8e_8.                        \tag{3}
\]

The disk defect Hessian on this interior coordinate is four.
Consequently the completed-square improvement from (2) to the
optimizer is exactly

\[
4|-8|^2=256.                                        \tag{4}
\]

Since the optimized face is `-64`, the transported certificate has
coefficient

\[
-64+256=192.                                        \tag{5}
\]

The earlier `(L,k)=(7,3)` audit has the same structure with correction
`+8c^3e_6`: the transported value is again `192` and optimization
subtracts `256`.

An adversarial length scan shows why this raw prefix must not be
mistaken for the theorem.  The `+8c^3e_(L-1)` difference persists for
grade three at lengths eight through ten, whereas the raw grade-four
difference is already zero at lengths ten through twelve.  Thus the
square center can move into the transported term as the terminal
geometry changes.  Only the quotient normal form below is a plausible
length-independent statement.

## 3. The sharpened all-size target

L131's reflected Faber row is

\[
r(u;c)=4\sum_j u_jc^j e_{L-j}^*,\qquad
\|r\|^2=16\sum_j|u_j|^2c^{2j}.                       \tag{6}
\]

The data and constants (3)--(5) identify the missing Schur statement
more precisely than “the kernel is diagonal.”  After removing all
strictly lower-weight transported-defect terms, there should be an
isometric endpoint-to-defect map `J_def` for which the associated
face is

\[
\boxed{
{\cal Q}_{\rm face}(u,\eta)
=4\|\eta-2J_{\rm def}r(u;c)\|^2-4\|r(u;c)\|^2 .
}                                                     \tag{7}
\]

Minimizing (7) gives

\[
-4\|r\|^2=-64\sum_j|u_j|^2c^{2j},                   \tag{8}
\]

while setting the residual correction to zero gives
`12||r||^2`, i.e. the observed `+192` one-grade transported value.
Because the reflected coordinate rows in (6) are orthogonal, (7)
would also kill every distinct-grade first-face cross term.

Equation (7) is a **conjectured normal form**, not a proved lemma.
The raw defect coordinate is gauge-dependent; `J_def` must be
derived after restoring homogeneous defect scale and quotienting its
null scaling direction.  Short grades and central collisions can
place the square center partly inside (2), so raw optimizer-prefix
comparison is not itself the invariant statement.

## 4. A spectral-factor clue in the universal defect Hessian

The quadratic part in the free defect tangent is independent of the
equality direction.  An exact LDL factorization of this universal
Hessian through length eight has the Newton edge

\[
\begin{aligned}
D_i&=4+O(c), &&i<L,\\
D_L&=\frac83+O(c),\\
{\cal L}_{i,i-2r}&=2c^r+O(c^{r+1}), &&i<L,\\
{\cal L}_{L,L-2r}&=\frac43c^r+O(c^{r+1}).
\end{aligned}                                       \tag{9}
\]

All lower and wrong-parity coefficients vanish exactly in the audit.
The interior lower factor is therefore the finite-path truncation of

\[
\frac{1+cz^2}{1-cz^2}
=1+2\sum_{r\ge1}c^rz^{2r}.                           \tag{10}
\]

This is not an arbitrary fitted series.  Equation (10) is the
analytic outer factor of the Newton edge of L117's periodized-sech
Szego weight: its Fourier coefficient at distance `r` starts
`2c^r`.  It gives a concrete route to `J_def`: conjugate the defect
quadratic by this lower factor, handle the endpoint `4/3` weight, and
then apply L135's DST-I transform.

The constants now fit without another numerical parameter.  L117's
axis defect has associated spatial edge

\[
d_{\rm edge}(x)
=1+4\sum_{r\ge1}(-x)^r
=\frac{1-3x}{1+x},\qquad x=c\,{\cal S}^2.             \tag{11}
\]

Therefore whitening the transported defect by (10) gives

\[
\frac{1+x}{1-x}\,2d_{\rm edge}(x)
=2\frac{1-3x}{1-x}
=2-4\sum_{r\ge1}x^r.                                \tag{12}
\]

At every first reflected grade the transported tail in (12) has
magnitude four, while L131's reflected Faber row has magnitude four
with the opposite normal orientation.  Their gap is eight, exactly
the correction measured in (3).  With disk Hessian weight four, the
square cost is `4*8^2=256`; the surviving negative row energy is
`-4*4^2=-64`.

Equations (10)--(12) explain both previously mysterious constants
`+192` and `-64`.  The remaining issue is rigor, not coefficient
discovery: prove the LDL edge (9) all-size and verify that reversal,
terminal weighting, and complex DST polarization implement the
claimed opposite orientation without an equal-weight alias.

The exact finite LDL audit does **not** prove (9) in every length.
An all-size proof should derive (9) directly from L117's Szego kernel
or from the Stein path recurrence, rather than extrapolate the
factorization.

## 5. Proof obligations

1. Restore the homogeneous defect variation, including its scale
   coordinate, and write the exact axis quadratic before choosing
   `x_0=0`.
2. Use L135's DCT/DST coordinates to identify the quotient map
   `J_def` and prove that it is isometric on the principal associated
   grade.
3. Show that every terminal fold not represented by (6) gains at
   least one additional power of `c`; this is the missing
   Schur/filtration commutation in A94.
4. Polarize the resulting real factorization over complex
   phase-palindromic coefficients.  L134 supplies the required
   within-grade calibration.

This route needs only the principal face.  A97 already shows that an
all-`c` central-window comparison is false, so (7) should not be
strengthened into such a comparison.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_principal_face_locality.py \
  --output experiments/crabb_principal_face_locality_s70223.jsonl
```

The checker reuses the guarded rational ellipse/Stein engine,
independently reconstructs `2U(u)d_c`, verifies agreement below the
reflected grade, and records the exact disk-Hessian cost of the first
transport correction.

The separate falsification scan

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_transport_correction_scan.py \
  --output experiments/crabb_transport_correction_scan_s70223.jsonl
```

checks that the raw correction itself is not length invariant.

The universal-Hessian factor audit is

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_defect_hessian_factor.py \
  --output experiments/crabb_defect_hessian_factor_s70223.jsonl
```
