"""Minimal first-order dual numbers — **legacy** first-order smoke test.

A dual number is ``a + b·eps`` with ``eps^2 = 0``. Evaluating a real function on
``x0 + 1·eps`` yields ``f(x0) + f'(x0)·eps`` exactly, so the imaginary part is the
first derivative with no truncation or subtractive-cancellation error.

This is intentionally tiny — enough to evaluate a scalar nonlinear residual with a
single perturbed parameter and prove the hypercomplex residual → R^(1) pathway. It
is **not** the production backend and must not be grown toward OTI. The production
endgame is arbitrary-order OTILib residual evaluation (see
``residual_core/algebra/otilib_adapter.py`` and
``residual_core/core/oti_rhs_provider.py``); higher orders and multivariate
truncation belong there, not here.

Supported ops (per task): + - * / , integer power, real/imag extraction, and
mixing with plain Python/NumPy real scalars.
"""

from __future__ import annotations

import numbers
from typing import Union

Number = Union[int, float, "Dual1"]


class Dual1:
    __slots__ = ("real", "imag")

    def __init__(self, real: float, imag: float = 0.0):
        self.real = float(real)
        self.imag = float(imag)

    # ---- helpers --------------------------------------------------------
    @staticmethod
    def _coerce(x: Number) -> "Dual1":
        if isinstance(x, Dual1):
            return x
        if isinstance(x, numbers.Real):
            return Dual1(float(x), 0.0)
        return NotImplemented

    def __repr__(self) -> str:
        return "Dual1(%r, %r)" % (self.real, self.imag)

    # ---- addition / subtraction ----------------------------------------
    def __add__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return Dual1(self.real + o.real, self.imag + o.imag)

    __radd__ = __add__

    def __sub__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return Dual1(self.real - o.real, self.imag - o.imag)

    def __rsub__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return Dual1(o.real - self.real, o.imag - self.imag)

    def __neg__(self) -> "Dual1":
        return Dual1(-self.real, -self.imag)

    def __pos__(self) -> "Dual1":
        return Dual1(self.real, self.imag)

    # ---- multiplication / division -------------------------------------
    def __mul__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        # (a + b eps)(c + d eps) = ac + (ad + bc) eps   (eps^2 = 0)
        return Dual1(self.real * o.real,
                     self.real * o.imag + self.imag * o.real)

    __rmul__ = __mul__

    def __truediv__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        # (a + b eps)/(c + d eps) = a/c + (b c - a d)/c^2 eps
        c = o.real
        return Dual1(self.real / c,
                     (self.imag * c - self.real * o.imag) / (c * c))

    def __rtruediv__(self, other: Number) -> "Dual1":
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return o.__truediv__(self)

    # ---- integer power --------------------------------------------------
    def __pow__(self, n: int) -> "Dual1":
        if not isinstance(n, numbers.Integral):
            raise TypeError("Dual1 supports integer powers only, got %r" % (n,))
        n = int(n)
        if n == 0:
            return Dual1(1.0, 0.0)
        # (a + b eps)^n = a^n + n a^(n-1) b eps
        if n > 0:
            a_pow = self.real ** (n - 1)
            return Dual1(a_pow * self.real, n * a_pow * self.imag)
        # negative integer power via reciprocal
        recip = Dual1(1.0, 0.0) / self
        return recip ** (-n)

    # ---- comparisons (compare on the real part) ------------------------
    def __float__(self) -> float:
        return self.real

    def __eq__(self, other) -> bool:
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return self.real == o.real and self.imag == o.imag


# --------------------------------------------------------------------------- #
# extraction helpers (work on Dual1 or plain reals)
# --------------------------------------------------------------------------- #
def real_part(x: Number) -> float:
    return x.real if isinstance(x, Dual1) else float(x)


def imag_part(x: Number) -> float:
    return x.imag if isinstance(x, Dual1) else 0.0


def seed(value: float) -> Dual1:
    """A perturbed input value ``value + 1·eps`` (unit first-order seed)."""
    return Dual1(float(value), 1.0)
