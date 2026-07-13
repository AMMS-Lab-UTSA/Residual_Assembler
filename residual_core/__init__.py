"""Residual Assembler -- a model-agnostic residual assembly framework.

The framework is formulation-agnostic *by architecture*. The core assembles

    R(y, q, a, t) = F_internal - F_external + F_constraints

by dispatching each element to a registered *formulation backend* and each
material point to a registered *material backend*. The core makes no assumption
about displacement-only DOFs, element topology, mechanics, stress/strain,
UMAT, or crystal plasticity. A formulation becomes *supported* once a backend
satisfying the formulation contract is registered and verified.

Verified example backends today:
  - stress-driven C3D8 solid (exported-field residual verification)
  - UMAT-driven C3D8 crystal plasticity (one example backend, not the core)
  - simple 2-node truss/bar (truss2)
  - simple 2-node 3D beam/frame (beam2)
  - UEL-direct residual adapter (skeleton)
  - shell backend contract (placeholder, documented, not yet runnable)

Quick start::

    from residual_core import ResidualProblem
    p = ResidualProblem.from_abaqus("model.inp")
    p.inspect()                              # elements, materials, missing data
    R = p.assemble(mode="stress-driven")     # once a field export is attached

or the ``resasm`` command line (``python -m residual_core.ui.cli``).
"""

from __future__ import annotations

from .ui.wizard import ResidualProblem
from .formulations.registry import (
    build_formulation_registry, default_formulation_registry)
from .materials.registry import (
    build_material_registry, default_material_registry)
from .core.results import (
    ResidualAssemblyResult, TangentResult, TangentSource, StateResult,
    ValidationReport)
from .core.sensitivity_package import (
    SensitivityRHSResult, SensitivityPackage, AlgebraMetadata)
from .core.rhs_provider import (
    SensitivityRHSProvider, DualNumberRHSProvider, default_rhs_provider)
from .algebra.dual1 import Dual1

__all__ = [
    "ResidualProblem",
    "build_formulation_registry", "default_formulation_registry",
    "build_material_registry", "default_material_registry",
    "ResidualAssemblyResult", "TangentResult", "TangentSource", "StateResult",
    "ValidationReport", "SensitivityRHSResult", "SensitivityPackage",
    "AlgebraMetadata", "SensitivityRHSProvider", "DualNumberRHSProvider",
    "default_rhs_provider", "Dual1",
]
