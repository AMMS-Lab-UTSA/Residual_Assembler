"""Stress-driven residual adapter (Mode 1).

The field entering the weak form (here: the integration-point Cauchy stress) is
supplied EXTERNALLY -- exported from Abaqus or another solver -- rather than
produced by a material update. The adapter just assembles

    f_int,e = sum_k B^T sigma_k dv        (no material, no tangent)

This verifies the finite-element residual assembly INDEPENDENTLY of any material
model: if Abaqus' own exported stresses reproduce Abaqus' reactions through this
adapter, the assembly is correct. Applicable to any formulation whose weak-form
field is known; this implementation delegates the B-operator/assembly to the
verified C3D8 kernel and reads sigma from ``fields``.

Expected ``fields`` layout:  fields['stress_ip'][element_id] -> (n_ip, 6) Cauchy
stress in Abaqus Voigt order. ``options['config']`` selects 'finite' (current
config, the nlgeom-matching path) or 'small' (reference config).
"""

from __future__ import annotations

import numpy as np

from .base import Formulation
from . import c3d8_kernel as k


class StressDrivenC3D8(Formulation):
    name = "stress_driven_c3d8"
    element_types = ("C3D8",)
    dof_types = ("UX", "UY", "UZ")
    sign_convention = "residual"       # returns +F_internal contribution to R
    verification_levels = (0, 1, 3, 5)
    supported_modes = ("stress-driven",)
    required_inputs = ("coords", "connectivity", "dofs", "stress_ip")
    optional_inputs = ("config",)
    limitations = ("C3D8 (8-node hex) only",
                   "needs exported integration-point stress; no material update",
                   "no material tangent (stress-driven)")
    verification_tests = ("zero-field", "patch", "global-equilibrium")
    notes = "Preferred verification mode: isolates assembly from the material update."
    required_inputs_by_mode = {
        "stress-driven": ("coords", "connectivity", "dofs", "stress_ip")}
    optional_inputs_by_mode = {"stress-driven": ("config",)}
    material_interface_needed = False
    state_requirements = "none"
    tangent_support = "none"
    verification_status = "verified"

    def n_nodes(self, element_type):
        return 8

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        fields = fields or {}
        config = options.get("config", "finite")
        coords = np.asarray(coords, dtype=float)
        u_e = np.asarray(dofs, dtype=float)

        stress_map = fields.get("stress_ip")
        if stress_map is None or element_id not in stress_map:
            raise ValueError("stress_driven_c3d8 requires fields['stress_ip'][%r] "
                             "(exported integration-point stress)." % (element_id,))
        sigma_ip = np.asarray(stress_map[element_id], dtype=float).reshape(-1, 6)

        if config == "finite":
            r = k.element_internal_force_finite_strain(coords, u_e, sigma_ip)
        elif config == "small":
            r = k.element_internal_force_small_strain(coords, sigma_ip)
        else:
            raise ValueError("options['config'] must be 'finite' or 'small'")
        # no material tangent in stress-driven mode; state passes through
        return r, None, material_state, {"formulation": self.name,
                                         "config": config, "n_ip": 8}
