"""Material backend contract.

A *material* answers the constitutive question only: given a kinematic input at a
material (integration) point and the previous state, return the stress, the
material tangent, and the updated state. It knows NOTHING about elements, shape
functions, or assembly -- that is the formulation's job (see formulations/base.py).

This is exactly the Abaqus UMAT contract, generalized:

    evaluate(kinematics, state_prev, binding, time, dtime, fields, options)
        -> (stress, tangent, state_new, diagnostics)

so a UMAT, a Python elastic law, or a crystal-plasticity return map all present
the same face to the formulation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class MaterialBinding:
    """What a formulation needs to drive a material at an element.

    Carried in the assembler's ``properties`` slot (per element). Keeps the
    material instance together with its constants and state size so the
    formulation never hard-codes any material knowledge.
    """
    material: "Material"
    constants: List[float] = field(default_factory=list)   # PROPS
    n_state_vars: int = 0                                   # NSTATV / *Depvar
    name: str = ""
    section: Optional[dict] = None


class Material(ABC):
    """Abstract constitutive backend.

    Subclasses declare their measures so the formulation can bridge correctly:

    - ``stress_measure``  : 'cauchy' (true stress, spatial) | 'pk2' | 'pk1'
    - ``tangent_measure`` : 'ddsdde' (Abaqus objective-rate d sigma/d eps) |
                            'material' (dS/dE) | None
    - ``kinematic_input`` : what evaluate() expects in ``kinematics`` --
                            'small_strain' (needs 'strain','dstrain') or
                            'deformation_gradient' (needs 'F0','F1').

    n_state_vars mirrors Abaqus NSTATV. The formulation supplies/receives the
    per-IP state slab; the material only reads state_prev and returns state_new.
    """

    name: str = "material"
    n_state_vars: int = 0
    stress_measure: str = "cauchy"
    tangent_measure: str = "ddsdde"
    kinematic_input: str = "small_strain"

    # ---- capability declaration (see docs/adding_a_material.md) ----
    # A material is NOT required to be a UMAT. These fields let the framework
    # describe any constitutive law: continuum stress-strain, beam moment-
    # curvature, cohesive traction-separation, thermal/diffusion flux, etc.
    constitutive_kind: str = "stress_strain"   # what law family this is
    input_variables: Tuple[str, ...] = ("strain", "dstrain")
    output_variables: Tuple[str, ...] = ("stress", "tangent")
    parameters: Tuple[str, ...] = ()           # PROPS names, if known
    supported_formulations: Tuple[str, ...] = ()   # empty = any compatible
    limitations: Tuple[str, ...] = ()
    notes: str = ""

    @property
    def tangent_available(self) -> bool:
        return self.tangent_measure is not None

    @property
    def spec(self):
        """Machine-readable capability descriptor (core.registry.BackendSpec)."""
        from ..core.registry import BackendSpec
        return BackendSpec(
            backend_name=self.name,
            kind="material",
            dof_types=(),
            required_inputs=tuple(self.input_variables),
            optional_inputs=tuple(self.parameters),
            supported_modes=tuple(self.supported_formulations),
            limitations=tuple(self.limitations),
            notes=self.notes or ("kind=%s; input=%s; stress=%s; tangent=%s"
                                 % (self.constitutive_kind, self.kinematic_input,
                                    self.stress_measure, self.tangent_measure)),
        )

    @abstractmethod
    def evaluate(self, kinematics: Dict[str, Any], state_prev: np.ndarray,
                 binding: MaterialBinding, time, dtime: float,
                 fields: Optional[dict], options: Optional[dict]):
        """Return (stress_voigt(6,), tangent(6,6) or None, state_new(nstatev,),
        diagnostics: dict).

        stress_voigt is in Abaqus Voigt order (11,22,33,12,13,23) and in the
        measure named by ``stress_measure``.
        """
        raise NotImplementedError

    def init_state(self, binding: MaterialBinding, coords=None,
                   n_ip: int = 8) -> np.ndarray:
        """Initial per-IP state slab (n_ip, n_state_vars). Override to seed
        history (e.g. a UMAT's first-increment STATEV initialization)."""
        n = binding.n_state_vars or self.n_state_vars
        return np.zeros((n_ip, n), dtype=float)
