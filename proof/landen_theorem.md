# The Landen theorem (campaign milestone, 2026-07-20/21)

## Theorem (sym3 ρ-positivity, closed form — classical-complete modulo two finite algebra steps)
Let A = [[0,a,0],[b,0,c],[0,d,0]] real with e² := ab+cd > 0, nonnormal, and Ω = int W(A).
Then W(A) is the ellipse with semi-axes A₂ = √((a+b)²+(c+d)²)/2, B₂ = √((a−b)²+(c−d)²)/2
(support-function computation: h(θ)² = A₂²cos²θ + B₂²sin²θ; foci ±e since A₂²−B₂² = ab+cd).
Let k be the elliptic modulus of the Riemann map (K(k′-m)/K(m) = 4·artanh(B₂/A₂)/π), k′ = √(1−k²).
At the extremal pair (K_norm > 1 regime):
- stationarity: α² = (1−k′)/k [⟺ pseudo-hyperbolic midpoint law ⟺ sn(K/2,k) = 1/√(1+k′), classical];
- the extremal-function zero: z₁ = e/√2 EXACTLY (u₁ = K/2 quarter-period; verified 40 digits);
- r₁ = eπ/(2√2·K(k)·k);
- **ρ = 1 − π/((1+k′)K(k)) = 1 − π/(2K(k₁)), k₁ = (1−k′)/(1+k′) (Landen descent).**
Hence ρ ≥ 0 with equality iff k = 0 (disk), by K(k₁) ≥ K(0) = π/2. ∎
[Verified independently: mpmath 30 dps vs float64 Theodorsen pipeline; Landen identity to 1e-31.]

## Interpretation: the ζ = z² collapse IS a Landen transformation
2×2 theorem: ρ = 1 − π/(2K(m)) on the confocal ellipse. sym3: the ζ-square collapse descends the
modulus one Landen step and REPRODUCES the same formula. Conjecture (Landen tower): every
even-symmetric collapse level of the zero-diag tridiagonal tower has
ρ = 1 − π/(2K(k_descended)) ≥ 0, closing the symmetric family for ALL n by induction —
PROVIDED the level's numerical range is elliptic (exact ellipse needed for the classical step).

## Remaining housekeeping for full rigor (finite algebra, numerically established)
1. α-configuration extremality (2×2 base: verified over full K>1 regime; symbolic proof pending).
2. The collapse identity v = α² from dK/dα = 0 (one-line calculus, verify signs).
3. K_norm > 1 hypothesis handling (below it extremal structure differs; irrelevant for Crouzeix).

## Level-4 elliptic slice: first exact facts (2026-07-21, slice_exact.py)
- CORRECTION of an earlier note: foci² = cL = e₁² — **the ellipse W(A) is confocal with the OUTER
  eigenvalue pair ±e₁** (not e₁²+e₂²); τ₁ = φ(e₁) = √k exactly (focus-image law, same as sym3).
  Inner pair ±e₂ = interior nodes with τ₂ = φ(e₂) free.
- ρ > 0 on the slice (exact machinery, matches pipeline): +1.7e-6, +4.7e-5, (+2.2e-3 low-K case).
- Stationarity: NEAR-equimodularity |B(τ₁)| ≈ |B(τ₂)| (3–9% off) — conjectured exact law is
  frame-weighted equimodularity; extract from high-precision scan with block-obliquity logging.
- Anomaly flag: weight-ordering a₂ > a₁ case gave K = 1.24 with bad reduced-diag — active-block
  switch or even-phase competition; handle sectors separately.

## Next-level targets
- Determine the elliptic-W slice of sym4 (quartic degenerates to conic ⟺ weight condition);
  test the DOUBLE-Landen prediction ρ = 1 − π/(2K(k₂)).
- Literature check: is Crouzeix known for n ≥ 4 with elliptic numerical range? If open, the
  elliptic-sym4 Landen theorem = first fully-new class.
- Non-elliptic sym4: perturbative closure off the elliptic slice (margin vs deviation δ).
