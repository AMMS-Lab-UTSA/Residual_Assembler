"""C3D8 small-strain solid formulation (Mode 2: material-update-driven).

Weak form (reference configuration, engineering-shear Voigt):

    f_int,e = sum_k  B0(X, xi_k)^T sigma_k  detJ0_k w_k
    K_e     = sum_k  B0^T D_k B0            detJ0_k w_k

sigma_k, D_k come from the bound Material at each integration point. This is the
small-strain sibling of solid_c3d8_finite_strain; use it for materials whose
kinematic_input is 'small_strain' (e.g. IsotropicElastic) and for the offline
Level-4 tangent verification. All numerics delegate to the verified c3d8_kernel.
"""

from __future__ import annotations

import numpy as np

from .base import Formulation
from . import c3d8_kernel as k


class SolidC3D8SmallStrain(Formulation):
    name = "solid_c3d8_small_strain"
    element_types = ("C3D8",)
    dof_types = ("UX", "UY", "UZ")
    sign_convention = "residual"       # returns +F_internal contribution
    verification_levels = (0, 1, 2, 3, 4, 5)
    supported_modes = ("material-replay",)
    required_inputs = ("coords", "connectivity", "dofs", "material")
    optional_inputs = ("dofs_prev", "material_state")
    limitations = ("small strain", "C3D8 (8-node hex) only",
                   "needs a small-strain material (kinematic_input='small_strain')")
    verification_tests = ("zero-field", "patch", "FD-tangent")
    required_inputs_by_mode = {
        "material-replay": ("coords", "connectivity", "dofs", "material",
                            "material_parameters")}
    optional_inputs_by_mode = {"material-replay": ("dofs_prev", "material_state")}
    material_interface_needed = True
    state_requirements = "material-defined (stateless for elastic)"
    tangent_support = "analytic"
    verification_status = "verified"

    def n_nodes(self, element_type):
        return 8

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        solution_state = solution_state or {}
        compute_tangent = options.get("compute_tangent", True)
        u_e = np.asarray(dofs, dtype=float)
        u_prev = np.asarray(solution_state.get("dofs_prev",
                                               np.zeros_like(u_e)), dtype=float)
        coords = np.asarray(coords, dtype=float)

        binding = properties
        material = getattr(binding, "material", None)
        if material is None:
            raise ValueError("solid_c3d8_small_strain needs a material binding "
                             "(Mode 2). For Mode 1 use stress_driven_adapter.")
        nstate = getattr(binding, "n_state_vars", 0) or material.n_state_vars
        pts, wts = k.ABAQUS_C3D8_GAUSS.points, k.ABAQUS_C3D8_GAUSS.weights

        r = np.zeros(24)
        K = np.zeros((24, 24)) if compute_tangent else None
        state_new = np.zeros((8, nstate)) if nstate else None
        for ip in range(8):
            B0, detJ0 = k.b_matrix_reference(coords, pts[ip])
            eps = B0 @ u_e
            deps = B0 @ (u_e - u_prev)
            s_prev = (material_state[ip] if material_state is not None
                      else np.zeros(nstate))
            sigma, D, s_new, _ = material.evaluate(
                {"strain": eps, "dstrain": deps}, s_prev, binding,
                time, dtime, fields, options)
            w = detJ0 * wts[ip]
            r += (B0.T @ sigma) * w
            if compute_tangent and D is not None:
                K += (B0.T @ D @ B0) * w
            if state_new is not None:
                state_new[ip] = s_new
        return r, K, state_new, {"formulation": self.name, "n_ip": 8}
