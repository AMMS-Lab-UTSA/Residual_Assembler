"""Isotropic linear-elastic material backend.

This is a real, RUNNABLE material update (unlike the crystal-plasticity UMAT,
which needs Intel ifort + Abaqus to compile). It exists so the full
material-update-driven pipeline (Mode 2) -- material update -> residual -> tangent
-> finite-difference verification -- can be exercised end-to-end offline through
the formulation-agnostic assembler. It also serves as the minimal reference
implementation of the Material contract for anyone adding a new material.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from .base import Material, MaterialBinding

# neutral shared helper (core), NOT a formulation import -> keeps materials
# independent of formulations (layering: materials depend downward on core only)
from ..core.voigt import isotropic_D


class IsotropicElastic(Material):
    """sigma = D : eps (small strain). Stateless (n_state_vars = 0).

    Constants convention (binding.constants): [E, nu]. Falls back to E=200e3,
    nu=0.3 if not supplied. Cauchy stress == PK for small strain; DDSDDE == D.
    """

    name = "isotropic_elastic"
    n_state_vars = 0
    stress_measure = "cauchy"
    tangent_measure = "ddsdde"
    kinematic_input = "small_strain"
    constitutive_kind = "stress_strain"
    input_variables = ("strain", "dstrain")
    output_variables = ("stress", "tangent")
    parameters = ("E", "nu")
    supported_formulations = ("solid_c3d8_small_strain",)
    notes = "Real runnable reference material (no Abaqus/ifort needed)."

    def _EnuD(self, binding: MaterialBinding):
        c = list(binding.constants) if binding and binding.constants else []
        E = float(c[0]) if len(c) >= 1 else 200000.0
        nu = float(c[1]) if len(c) >= 2 else 0.3
        return E, nu, isotropic_D(E, nu)

    def evaluate(self, kinematics: Dict[str, Any], state_prev: np.ndarray,
                 binding: MaterialBinding, time, dtime: float,
                 fields: Optional[dict], options: Optional[dict]):
        E, nu, D = self._EnuD(binding)
        eps = np.asarray(kinematics["strain"], dtype=float)   # (6,) engineering shear
        sigma = D @ eps
        diag = {"material": self.name, "E": E, "nu": nu}
        return sigma, D, (state_prev if state_prev is not None else np.zeros(0)), diag
