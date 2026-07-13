"""Requested module name ``solid_c3d8`` — re-exports the C3D8 solid backends.

The concrete implementations live in ``solid_c3d8_small_strain`` and
``solid_c3d8_finite_strain`` (kept as-is so existing CP tests and imports do not
break). This module gives the canonical short name asked for by the project
layout and a convenient default.
"""

from .solid_c3d8_small_strain import SolidC3D8SmallStrain
from .solid_c3d8_finite_strain import SolidC3D8FiniteStrain

# Canonical default when a caller just says "solid_c3d8": finite strain, the home
# of the verified crystal-plasticity backend.
SolidC3D8 = SolidC3D8FiniteStrain

__all__ = ["SolidC3D8", "SolidC3D8SmallStrain", "SolidC3D8FiniteStrain"]
