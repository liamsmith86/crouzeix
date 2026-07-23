# Weighted transition on the `3x3` Crabb slice (2026-07-22)

## 1. Weighted chart

Use the exact local slice from L69 and write

\[
 A_\epsilon=C_3+\epsilon R_1(z)+\epsilon^2R_2(u)
 +\epsilon^3\{sY+vV\},                                      \tag{1}
\]

where `z,u,v` are complex, `s` is real, `Y_01=1`, `Y_12=-1`, and `V_20=1`.
The circle action has weights one, two, and three on `(z,u,v)`, respectively,
and fixes `s`.  For `z != 0`, scaling and the circle action reduce the leading
calculation to `z=1`.  This is the joint regime left open by L67--L69.

For the analytic rank-one Stein branch, let

\[
 \kappa(P_\epsilon)=4+\epsilon^6H(u,v,s)+O(\epsilon^7).       \tag{2}
\]

The remainder is uniform when `(u,v,s)` stays in a fixed compact set.  This is
a feasible-metric upper certificate; no assertion that this branch is the
global L21 optimizer is needed.

## 2. Exact coefficient extraction

The required weighted coefficients occur at much lower ordinary orders:

* L67 gives the constant `C=-171/4096`.
* L68 gives the mode-two quadratic coefficient `A=-31/8`.
* The ordinary fifth-order invariant linear in the mode-two coordinate is
  `B=-123 sqrt(2)/256`.  It is separated from the cubic invariant using the
  two paths `w=1,2`; the cubic coefficient is `23 sqrt(2)/8` and does not enter
  weighted order six.
* The ordinary fourth-order term linear in `v` is `E=-117/128`.
* The ordinary cubic mode-two/bottom coupling is `D=-6 sqrt(2)`.
* L63 supplies the transverse squares `-21|v|^2/4` and `-8s^2`.

Circle invariance and conjugation leave exactly the corresponding complex
invariants.  Thus, after setting `z=1`,

\[
\begin{aligned}
H(u,v,s)={}&-\frac{171}{4096}-\frac{31}{8}|u|^2
-\frac{123\sqrt2}{256}\operatorname{Re}u
-\frac{117}{128}\operatorname{Re}v\\
&-6\sqrt2\operatorname{Re}(v\bar u)
-\frac{21}{4}|v|^2-8s^2.                                  \tag{3}
\end{aligned}
\]

The sparse exact regeneration derives every coefficient from the support
characteristic equation, boundary reparameterization, inverse functional
calculus, and Stein recurrence.  It does not fit numerical conformal maps.

## 3. Complete-square collapse

Exact algebra turns (3) into

\[
\boxed{
H=-8s^2
-\frac{21}{4}\left|v+\frac{39}{448}+\frac{4\sqrt2}{7}u\right|^2
-\frac{25}{56}\left|u-\frac{3\sqrt2}{64}\right|^2.}        \tag{4}
\]

Consequently the leading weighted form is nonpositive and has the unique zero

\[
 u_0=\frac{3\sqrt2}{64},\qquad v_0=-\frac9{64},\qquad s_0=0. \tag{5}
\]

This closes the feared sign obstruction in the weighted transition: every
bounded weighted direction except (5) has a strictly negative leading
feasible certificate.  It does **not** yet prove a full punctured-neighbourhood
theorem, because the unique center (5) must be controlled at higher order.

## 4. Equality-center phenomenon

Exact continuation reveals substantial additional flatness.  On
the path (5), the optimized coefficients through order eight vanish and the
fixed-center order-ten coefficient is `-13851/4194304`.  However the order-eight
gradient shifts the center by

\[
 u=u_0+\frac{27\sqrt2}{2048}\epsilon^2+O(\epsilon^4),\qquad
 v=v_0-\frac{81}{2048}\epsilon^2+O(\epsilon^4),             \tag{6}
\]

and completing (4) adds exactly `+13851/4194304`, cancelling the fixed-center
order-ten descent.  With (6), the exact certificate remains flat through order
twelve; the `Y` direction contributes `-8s^2` at that scale.  This strongly
suggests a hidden analytic equality center for the rank-one certificate, not a
positive direction.  A structural characterization is preferable to an
indefinite sequence of higher-order expansions.

Section 4 records the audited frontier but is not used to upgrade the theorem
in Section 3.  Its higher-jet regeneration is separated from the load-bearing
leading-form certificate.

## 5. Regeneration

Run

```bash
.venv/bin/python -u experiments/p3_crabb_weighted_slice.py
.venv/bin/python -u experiments/p3_crabb_center_jet.py
```

`experiments/p3_sparse_series.py` is specialized to the canonical real `p=3`
Crabb path.  Against the older independent symbolic implementation it agrees
through order four, reproduces L67--L68, and gives the same ordinary fifth
coefficient when its Riemann output is passed to the older Stein engine.
The second command independently rebuilds the fixed-center tenth coefficient,
the eighth-order center gradient and completion, the corrected order-twelve
jet, and the reversal-transpose symmetry controlling the last real direction.
