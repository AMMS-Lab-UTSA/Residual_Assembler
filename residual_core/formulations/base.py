"""Formulation backend contract.

A *formulation* owns the finite-element weak form for one element family: its
DOFs, shape functions, integration rule, the kinematics computed from the DOFs,
how it obtains the field entering the weak form (from a material update, from an
exported field, or from a user element), and how it turns that into the element
residual and tangent.

The global assembler (core/assembler.py) is formulation-AGNOSTIC: it only calls

    eval_element(element_id, element_type, coords, dofs, solution_state,
                 material_state, properties, time, dtime, fields, options)
        -> (element_residual, element_tangent, updated_state, diagnostics)

and scatters the result. All element-specific mathematics lives in subclasses of
Formulation, here under formulations/. See docs/formulation_contract.md for the
11-point contract every backend must satisfy before it is trusted.

Generic residual (documented in docs/architecture.md):

    R(y, q, a, t) = F_internal(y, q, a, t) - F_external(t) + F_constraints(y, t)

eval_element returns only the element's F_internal contribution (and its
tangent); F_external and F_constraints are added by the core assembler.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import numpy as np

from ..core.registry import BackendSpec


@dataclass
class ElementResult:
    """Convenience container; eval_element returns the 4-tuple form for the
    assembler, but backends may build this and call ``.as_tuple()``."""
    residual: np.ndarray
    tangent: Optional[np.ndarray] = None
    updated_state: Optional[np.ndarray] = None
    diagnostics: Dict[str, Any] = field(default_factory=dict)

    def as_tuple(self):
        return self.residual, self.tangent, self.updated_state, self.diagnostics


class Formulation(ABC):
    """Abstract element formulation.

    Class attributes a backend MUST set:
      name            : unique registry key
      element_types   : tuple of Abaqus element types handled (e.g. ("C3D8",))
      dof_types       : per-node DOFs, ordered (e.g. ("UX","UY","UZ"))
      sign_convention : "residual" if the returned vector is +F_internal
                        contribution to R = Fint - Fext (the standard here), or
                        "rhs" if it is Abaqus-UEL RHS = -R (must be documented).
      verification_levels : tuple of ints from docs/verification_strategy.md that
                        apply to this backend.
    """

    name: str = "formulation"
    element_types: Tuple[str, ...] = ()
    dof_types: Tuple[str, ...] = ("UX", "UY", "UZ")
    sign_convention: str = "residual"
    verification_levels: Tuple[int, ...] = ()

    # ---- capability declaration (consumed by the model inspector / CLI) ----
    # These let the framework reason about a backend without importing physics.
    supported_modes: Tuple[str, ...] = ()      # e.g. ("stress-driven","material-replay")
    required_inputs: Tuple[str, ...] = ()      # e.g. ("coords","connectivity","dofs")
    optional_inputs: Tuple[str, ...] = ()
    limitations: Tuple[str, ...] = ()
    verification_tests: Tuple[str, ...] = ()
    notes: str = ""

    # ---- richer capability declaration (audit item 3) ----------------------
    required_inputs_by_mode: Dict[str, Tuple[str, ...]] = {}   # per-mode inputs
    optional_inputs_by_mode: Dict[str, Tuple[str, ...]] = {}
    material_interface_needed: bool = False    # calls a Material.evaluate()?
    state_requirements: str = "none"           # none | history-dependent
    tangent_support: str = "none"              # analytic | fd | none
    verification_status: str = "unspecified"   # verified | simple | contract-only | skeleton

    @property
    def dofs_per_node(self) -> int:
        return len(self.dof_types)

    @property
    def spec(self) -> BackendSpec:
        """Machine-readable capability descriptor (see core/registry.py)."""
        return BackendSpec(
            backend_name=self.name,
            kind="formulation",
            supported_element_types=tuple(self.element_types),
            dof_types=tuple(self.dof_types),
            supported_modes=tuple(self.supported_modes),
            required_inputs=tuple(self.required_inputs),
            optional_inputs=tuple(self.optional_inputs),
            verification_tests=tuple(self.verification_tests),
            limitations=tuple(self.limitations),
            notes=self.notes,
            required_inputs_by_mode={k: tuple(v) for k, v in
                                     dict(self.required_inputs_by_mode).items()},
            optional_inputs_by_mode={k: tuple(v) for k, v in
                                     dict(self.optional_inputs_by_mode).items()},
            material_interface_needed=bool(self.material_interface_needed),
            state_requirements=self.state_requirements,
            tangent_support=self.tangent_support,
            verification_status=self.verification_status,
        )

    @abstractmethod
    def n_nodes(self, element_type: str) -> int:
        """Number of nodes for the given element type."""
        raise NotImplementedError

    @abstractmethod
    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        """Evaluate one element.

        Parameters
        ----------
        element_id     : int
        element_type   : str (e.g. "C3D8")
        coords         : (n_nodes, ndim) reference nodal coordinates
        dofs           : (n_nodes*dofs_per_node,) element solution (node-major)
        solution_state : dict of extra solution context (e.g. {'U_prev': ...})
        material_state : (n_ip, n_state) previous per-IP state, or None
        properties     : per-element properties (e.g. a MaterialBinding), or None
        time           : (2,) [step time, total time] at increment start
        dtime          : float increment
        fields         : dict of extra fields (e.g. {'stress_ip': (n_ip,6)})
        options        : dict of flags (e.g. {'compute_tangent': True})

        Returns
        -------
        (element_residual (ndof,),
         element_tangent  (ndof,ndof) or None,
         updated_state    (n_ip,n_state) or None,
         diagnostics      dict)
        """
        raise NotImplementedError
