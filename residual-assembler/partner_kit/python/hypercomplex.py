"""Minimal hypercomplex scalars for the partner kit (self-contained).

This is the *seam* where the derivative scalar type is chosen. It provides a
first-order dual number (``Dual1``) as a **legacy smoke test** — enough to extract
R^(1) with no dependencies so a partner can run the kit standalone. The
**production** backend is OTILib (arbitrary order, multi-parameter); see
``docs/otilib_partner_quickstart.md``. Do not grow ``Dual1`` toward OTI — swap in
the OTILib scalar instead.

Design intent: the kit passes a *seeded* parameter into the partner's
scalar-generic residual and reads the derivative from the imaginary part. Swapping
``Dual1`` for an OTILib scalar (higher orders, multivariate truncation) requires
only that the partner's residual stay written against a generic scalar type — the
extraction API (``real_part`` / ``imag_part`` / ``seed``) stays the same.
"""

from __future__ import annotations

import numbers
from typing import Union

Number = Union[int, float, "Dual1"]


class Dual1:
    """First-order dual number ``a + b·eps`` with ``eps^2 = 0``.

    Evaluating a real function on ``x0 + 1·eps`` gives ``f(x0) + f'(x0)·eps``.
    Supports the operations a scalar residual typically needs.
    """

    __slots__ = ("real", "imag")

    def __init__(self, real: float, imag: float = 0.0):
        self.real = float(real)
        self.imag = float(imag)

    @staticmethod
    def _c(x: Number):
        if isinstance(x, Dual1):
            return x
        if isinstance(x, numbers.Real):
            return Dual1(float(x), 0.0)
        return NotImplemented

    def __repr__(self):
        return "Dual1(%r, %r)" % (self.real, self.imag)

    def __add__(self, o):
        o = self._c(o)
        return NotImplemented if o is NotImplemented else Dual1(self.real + o.real, self.imag + o.imag)
    __radd__ = __add__

    def __sub__(self, o):
        o = self._c(o)
        return NotImplemented if o is NotImplemented else Dual1(self.real - o.real, self.imag - o.imag)

    def __rsub__(self, o):
        o = self._c(o)
        return NotImplemented if o is NotImplemented else Dual1(o.real - self.real, o.imag - self.imag)

    def __neg__(self):
        return Dual1(-self.real, -self.imag)

    def __mul__(self, o):
        o = self._c(o)
        if o is NotImplemented:
            return NotImplemented
        return Dual1(self.real * o.real, self.real * o.imag + self.imag * o.real)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        if o is NotImplemented:
            return NotImplemented
        c = o.real
        return Dual1(self.real / c, (self.imag * c - self.real * o.imag) / (c * c))

    def __rtruediv__(self, o):
        o = self._c(o)
        return NotImplemented if o is NotImplemented else o.__truediv__(self)

    def __pow__(self, n):
        if not isinstance(n, numbers.Integral):
            raise TypeError("Dual1 supports integer powers only")
        n = int(n)
        if n == 0:
            return Dual1(1.0, 0.0)
        if n > 0:
            a = self.real ** (n - 1)
            return Dual1(a * self.real, n * a * self.imag)
        return (Dual1(1.0, 0.0) / self) ** (-n)

    def __float__(self):
        return self.real


def real_part(x: Number) -> float:
    return x.real if isinstance(x, Dual1) else float(x)


def imag_part(x: Number) -> float:
    return x.imag if isinstance(x, Dual1) else 0.0


def seed(value: float) -> Dual1:
    """A first-order unit seed ``value + 1·eps`` for one perturbed parameter."""
    return Dual1(float(value), 1.0)


# ---- scalar-type seam (for future OTI/HYPAD) ------------------------------- #
# The kit selects a scalar 'algebra' by name. Only 'dual1' (order 1) is
# implemented; 'oti'/'hypad' are declared so higher orders can plug in later
# WITHOUT changing the partner's residual code or the extraction API.
SUPPORTED_ALGEBRAS = ("dual1",)
PLANNED_ALGEBRAS = ("oti", "hypad")


def make_seed(value: float, algebra: str = "dual1", order: int = 1):
    if algebra == "dual1":
        if order != 1:
            raise NotImplementedError("dual1 supports order 1 only; order=%d "
                                      "needs an OTI/HYPAD scalar." % order)
        return seed(value)
    raise NotImplementedError(
        "algebra %r not implemented in the kit yet (planned: %s). The residual "
        "must remain scalar-generic so it can be swapped in later."
        % (algebra, ", ".join(PLANNED_ALGEBRAS)))
