"""High-precision checks for the analytic EL4 Schwarzian proof.

The proof is in ``proof/el4_schwarzian_theorem.md``.  This script independently checks:

1. the Schwarzian gap against its positive Fourier series;
2. the identity between the coordinate-free D2 functional and the EL4 formula; and
3. the resulting nonnegative EL4 margins on a deterministic grid.

These checks are regression tests, not substitutes for the proof.
"""

from __future__ import annotations

from dataclasses import dataclass

from mpmath import ellipfun, ellipk, mp


sn = ellipfun("sn")
cn = ellipfun("cn")
dn = ellipfun("dn")


@dataclass(frozen=True)
class EllipticData:
    k: mp.mpf
    m: mp.mpf
    K: mp.mpf
    Kp: mp.mpf
    s: mp.mpf
    nome: mp.mpf

    @classmethod
    def from_k(cls, k: mp.mpf) -> "EllipticData":
        m = k * k
        K = ellipk(m)
        Kp = ellipk(1 - m)
        return cls(k=k, m=m, K=K, Kp=Kp, s=mp.pi / (2 * K), nome=mp.exp(-mp.pi * Kp / K))


def hyperbolic_midpoint(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """Return the real pseudo-hyperbolic midpoint of b < a in (-1, 1)."""
    return mp.tanh((mp.atanh(a) + mp.atanh(b)) / 2)


def schwarzian_gap(data: EllipticData, U: mp.mpf) -> mp.mpf:
    """Return S(sin^2(sU)) - S(k sn^2(U)); it must be nonnegative."""
    x = 2 * U
    elliptic_term = 1 / sn(x, data.m) ** 2
    trigonometric_term = data.s**2 / mp.sin(data.s * x) ** 2
    return 6 * (elliptic_term - trigonometric_term) - 2 * (1 + data.m - data.s**2)


def fourier_gap(data: EllipticData, U: mp.mpf, terms: int = 500) -> mp.mpf:
    """Positive Fourier series equal to one sixth of ``schwarzian_gap``."""
    x = 2 * U
    total = mp.mpf("0")
    for n in range(1, terms + 1):
        q2n = data.nome ** (2 * n)
        coefficient = n * q2n / (1 - q2n)
        total += coefficient * (1 - mp.cos(2 * n * data.s * x))
    return 8 * data.s**2 * total


def el4_values(data: EllipticData, U2: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """Return (coordinate-free D2 value, displayed EL4 Theta)."""
    a = data.k
    b = data.k * sn(U2, data.m) ** 2
    v = hyperbolic_midpoint(a, b)
    W = mp.ellipf(mp.asin(mp.sqrt(v / data.k)), data.m)
    delta = (a - v) / (1 - a * v)

    Ga = mp.mpf("1")
    Gb = mp.sin(data.s * U2) ** 2
    Gv = mp.sin(data.s * W) ** 2
    Gprime = (
        data.s
        * mp.sin(data.s * W)
        * mp.cos(data.s * W)
        / (data.k * sn(W, data.m) * cn(W, data.m) * dn(W, data.m))
    )
    d2_value = delta * (1 - v * v) * Gprime * (1 / (Ga - Gv) + 1 / (Gv - Gb)) / 2

    theta = (
        delta
        * data.s
        * (1 - v * v)
        * mp.sin(data.s * W)
        * mp.cos(data.s * U2) ** 2
        / (
            2
            * data.k
            * sn(W, data.m)
            * cn(W, data.m)
            * dn(W, data.m)
            * mp.cos(data.s * W)
            * (mp.sin(data.s * W) ** 2 - mp.sin(data.s * U2) ** 2)
        )
    )
    return d2_value, theta


def main() -> None:
    mp.dps = 70
    max_fourier_error = mp.mpf("0")
    max_formula_error = mp.mpf("0")
    min_schwarzian_gap = mp.inf
    min_el4_margin = mp.inf

    for k_text in ("0.05", "0.2", "0.5", "0.8", "0.95"):
        data = EllipticData.from_k(mp.mpf(k_text))
        for index in range(1, 20):
            U = data.K * index / 20
            gap = schwarzian_gap(data, U)
            series_gap = fourier_gap(data, U)
            min_schwarzian_gap = min(min_schwarzian_gap, gap)
            max_fourier_error = max(max_fourier_error, abs(gap / 6 - series_gap))

            d2_value, theta = el4_values(data, U)
            max_formula_error = max(max_formula_error, abs(d2_value - theta))
            min_el4_margin = min(min_el4_margin, 1 - theta)

    tolerance = mp.mpf("1e-55")
    assert min_schwarzian_gap > 0
    assert min_el4_margin > 0
    assert max_fourier_error < tolerance
    assert max_formula_error < tolerance

    print(f"minimum Schwarzian gap: {mp.nstr(min_schwarzian_gap, 12)}")
    print(f"minimum EL4 margin:       {mp.nstr(min_el4_margin, 12)}")
    print(f"max Fourier error:        {mp.nstr(max_fourier_error, 5)}")
    print(f"max D2/Theta error:       {mp.nstr(max_formula_error, 5)}")


if __name__ == "__main__":
    main()
