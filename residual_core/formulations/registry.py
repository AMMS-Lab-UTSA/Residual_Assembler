"""Formulation registry — name -> Formulation instance, with capability specs.

This builds a ``core.registry.Registry`` of every formulation backend the
framework ships. The core assembler resolves ``model.element_formulation`` keys
through this registry; the model inspector and requirements engine reason over
each backend's ``.spec`` (a ``BackendSpec``) without importing any physics.

Crystal plasticity is *not* special here — ``solid_c3d8_finite_strain`` is just
one of several registered backends (truss, beam, stress-driven, UEL, shell
placeholder).
"""

from __future__ import annotations

from ..core.registry import Registry
from .solid_c3d8_small_strain import SolidC3D8SmallStrain
from .solid_c3d8_finite_strain import SolidC3D8FiniteStrain
from .stress_driven_adapter import StressDrivenC3D8
from .truss2 import Truss2
from .beam2 import Beam2
from .shell_placeholder import ShellPlaceholder
from .nonlinear_spring1 import NonlinearSpring1
from .nonlinear_bar1 import NonlinearBar1


_ALWAYS = (SolidC3D8SmallStrain, SolidC3D8FiniteStrain, StressDrivenC3D8,
           Truss2, Beam2, ShellPlaceholder, NonlinearSpring1, NonlinearBar1)


def build_formulation_registry() -> Registry:
    reg = Registry("formulation")
    for cls in _ALWAYS:
        reg.register(cls())
    # optional UEL adapter (Mode 3); registered unconfigured (no routine yet)
    try:
        from .uel_adapter import UelAdapter
        reg.register(UelAdapter())
    except Exception:
        pass
    return reg


def default_formulation_registry() -> Registry:
    return build_formulation_registry()
