"""Formulation backends + a default registry.

Registry maps a formulation key -> a Formulation instance. The assembler resolves
element bindings (model.element_formulation) through this dict; it never imports
a concrete formulation itself.
"""

from .base import Formulation, ElementResult
from .solid_c3d8_small_strain import SolidC3D8SmallStrain
from .solid_c3d8_finite_strain import SolidC3D8FiniteStrain
from .stress_driven_adapter import StressDrivenC3D8
from .truss2 import Truss2
from .beam2 import Beam2
from .shell_placeholder import ShellPlaceholder
from .nonlinear_spring1 import NonlinearSpring1
from .nonlinear_bar1 import NonlinearBar1
from .registry import build_formulation_registry, default_formulation_registry


def default_formulations():
    """Return {name: instance} for all registered backends (dict view).

    Backward-compatible with the original API. For the richer capability-aware
    view (BackendSpec per backend) use ``default_formulation_registry()``.
    """
    reg = {}
    for cls in (SolidC3D8SmallStrain, SolidC3D8FiniteStrain, StressDrivenC3D8,
                Truss2, Beam2, ShellPlaceholder, NonlinearSpring1, NonlinearBar1):
        inst = cls()
        reg[inst.name] = inst
    # optional UEL adapter (Mode 3) if present
    try:
        from .uel_adapter import UelAdapter
        u = UelAdapter()
        reg[u.name] = u
    except Exception:
        pass
    return reg


__all__ = ["Formulation", "ElementResult", "SolidC3D8SmallStrain",
           "SolidC3D8FiniteStrain", "StressDrivenC3D8", "Truss2", "Beam2",
           "ShellPlaceholder", "NonlinearSpring1", "default_formulations",
           "build_formulation_registry", "default_formulation_registry"]
