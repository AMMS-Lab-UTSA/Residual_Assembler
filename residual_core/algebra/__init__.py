"""Hypercomplex / automatic-differentiation algebras.

Currently only a tiny first-order dual number (``dual1.Dual1``) used to prove the
hypercomplex residual → sensitivity RHS pathway. This is deliberately minimal and
is NOT the production OTI/HYPAD library.
"""

from .dual1 import Dual1, real_part, imag_part, seed

__all__ = ["Dual1", "real_part", "imag_part", "seed"]
