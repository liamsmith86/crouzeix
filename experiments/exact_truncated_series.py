"""Small exact-arithmetic truncated power-series utility.

The research checkers only need univariate series with rational
coefficients.  Keeping this minimal implementation separate avoids
repeating truncation, valuation, division, and square-root logic in
each formal checker.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable


@dataclass(frozen=True)
class Series:
    """A rational power series truncated at a fixed exclusive order."""

    coefficients: tuple[Fraction, ...]

    @classmethod
    def constant(cls, value: int | Fraction, order: int) -> Series:
        """Construct a constant series."""

        if order <= 0:
            raise ValueError("the truncation order must be positive")
        return cls((Fraction(value), *(Fraction(0) for _ in range(order - 1))))

    @classmethod
    def from_coefficients(
        cls,
        coefficients: Iterable[int | Fraction],
        order: int,
    ) -> Series:
        """Construct a series, padding or truncating to ``order``."""

        values = [Fraction(value) for value in coefficients]
        values = values[:order]
        values.extend(Fraction(0) for _ in range(order - len(values)))
        return cls(tuple(values))

    @classmethod
    def monomial(
        cls,
        degree: int,
        order: int,
        coefficient: int | Fraction = 1,
    ) -> Series:
        """Construct ``coefficient * c**degree``."""

        if degree < 0:
            raise ValueError("the monomial degree must be nonnegative")
        values = [Fraction(0)] * order
        if degree < order:
            values[degree] = Fraction(coefficient)
        return cls(tuple(values))

    @property
    def order(self) -> int:
        """Return the exclusive truncation order."""

        return len(self.coefficients)

    def _coerce(self, other: int | Fraction | Series) -> Series:
        if isinstance(other, Series):
            if other.order != self.order:
                raise ValueError("series truncation orders do not match")
            return other
        return Series.constant(other, self.order)

    def __add__(self, other: int | Fraction | Series) -> Series:
        other_series = self._coerce(other)
        return Series(
            tuple(
                left + right
                for left, right in zip(
                    self.coefficients,
                    other_series.coefficients,
                    strict=True,
                )
            )
        )

    __radd__ = __add__

    def __neg__(self) -> Series:
        return Series(tuple(-value for value in self.coefficients))

    def __sub__(self, other: int | Fraction | Series) -> Series:
        return self + (-self._coerce(other))

    def __rsub__(self, other: int | Fraction | Series) -> Series:
        return self._coerce(other) - self

    def __mul__(self, other: int | Fraction | Series) -> Series:
        other_series = self._coerce(other)
        values = []
        for degree in range(self.order):
            values.append(
                sum(
                    (
                        self.coefficients[index]
                        * other_series.coefficients[degree - index]
                        for index in range(degree + 1)
                    ),
                    Fraction(0),
                )
            )
        return Series(tuple(values))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> Series:
        if exponent < 0:
            raise ValueError("negative powers are not supported")
        result = Series.constant(1, self.order)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result *= base
            base *= base
            remaining //= 2
        return result

    def valuation(self) -> int:
        """Return the first nonzero degree, or ``order`` for zero."""

        return next(
            (
                degree
                for degree, coefficient in enumerate(self.coefficients)
                if coefficient
            ),
            self.order,
        )

    def inverse_unit(self) -> Series:
        """Invert a series with nonzero constant coefficient."""

        constant = self.coefficients[0]
        if not constant:
            raise ZeroDivisionError("the series is not a unit")
        values = [Fraction(0)] * self.order
        values[0] = 1 / constant
        for degree in range(1, self.order):
            values[degree] = -sum(
                (
                    self.coefficients[index] * values[degree - index]
                    for index in range(1, degree + 1)
                ),
                Fraction(0),
            ) / constant
        return Series(tuple(values))

    def __truediv__(self, other: int | Fraction | Series) -> Series:
        other_series = self._coerce(other)
        numerator_valuation = self.valuation()
        denominator_valuation = other_series.valuation()
        if denominator_valuation == self.order:
            raise ZeroDivisionError("division by the zero series")
        if numerator_valuation < denominator_valuation:
            raise ValueError("the quotient would contain negative powers")

        numerator = Series.from_coefficients(
            self.coefficients[numerator_valuation:],
            self.order,
        )
        denominator = Series.from_coefficients(
            other_series.coefficients[denominator_valuation:],
            self.order,
        )
        shift = Series.monomial(
            numerator_valuation - denominator_valuation,
            self.order,
        )
        return shift * numerator * denominator.inverse_unit()

    def square_root_unit(self) -> Series:
        """Return the normalized square root when the constant is one."""

        if self.coefficients[0] != 1:
            raise ValueError("only square roots with constant one are supported")
        values = [Fraction(0)] * self.order
        values[0] = Fraction(1)
        for degree in range(1, self.order):
            quadratic = sum(
                (
                    values[index] * values[degree - index]
                    for index in range(1, degree)
                ),
                Fraction(0),
            )
            values[degree] = (self.coefficients[degree] - quadratic) / 2
        return Series(tuple(values))

    def coefficient(self, degree: int) -> Fraction:
        """Return one coefficient, treating degrees past the cutoff as zero."""

        if degree < 0:
            raise ValueError("the coefficient degree must be nonnegative")
        if degree >= self.order:
            return Fraction(0)
        return self.coefficients[degree]

    def truncated_coefficients(self, exclusive_order: int) -> tuple[Fraction, ...]:
        """Return coefficients below ``exclusive_order``."""

        return self.coefficients[:exclusive_order]

    def __str__(self) -> str:
        terms = [
            f"{coefficient}*c^{degree}"
            for degree, coefficient in enumerate(self.coefficients)
            if coefficient
        ]
        return " + ".join(terms) if terms else "0"
