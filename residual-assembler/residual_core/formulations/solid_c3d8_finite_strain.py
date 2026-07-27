"""C3D8 finite-strain solid formulation (Mode 2: material-update-driven).

This is the home of the verified crystal-plasticity backend. The examples are
``nlgeom=YES`` and the CP UMAT returns Cauchy stress, so the internal force is
assembled in the CURRENT configuration:

    f_int,e = sum_k  B_spatial(x, xi_k)^T sigma_cauchy_k  detJ_current_k w_k
    K_e     = sum_k  B_spatial^T D_k B_spatial detJ_k w_k + K_geo_k

The bound Material (kinematic_input='deformation_gradient') receives F0/F1 per IP
and returns Cauchy stress + DDSDDE + updated STATEV. All element numerics delegate
to the verified c3d8_kernel; this class only bridges DOFs<->kinematics<->material
and packs/unpacks per-IP state.

NOTE on the tangent measure: turning the finite-strain UMAT DDSDDE into the exact
Abaqus AMATRX (objective-rate + geometric split) is the error-prone item deferred
to the Abaqus comparison (see docs/limitations.md and tests). The residual (which
only needs Cauchy stress) is exact; the tangent here uses the kernel's
conventional geometric term with the material DDSDDE.
"""

from __future__ import annotations

import numpy as np

from .base import Formulation
from . import c3d8_kernel as k


def _F_at(coords, u_e, u_prev, xi):
    """Deformation gradients F0 (prev) and F1 (current) at natural coord xi,
    using reference-config shape-function gradients: F = I + du/dX."""
    _, _, dNdX = k._b_from_coords(coords, xi)          # (8,3) d N_a / dX
    ue = u_e.reshape(8, 3)
    up = u_prev.reshape(8, 3)
    F1 = np.eye(3) + ue.T @ dNdX
    F0 = np.eye(3) + up.T @ dNdX
    return F0, F1


class SolidC3D8FiniteStrain(Formulation):
    name = "solid_c3d8_finite_strain"
    element_types = ("C3D8",)
    dof_types = ("UX", "UY", "UZ")
    sign_convention = "residual"
    verification_levels = (0, 1, 2, 3, 4, 5, 6, 7)
    supported_modes = ("material-replay",)
    required_inputs = ("coords", "connectivity", "dofs", "material",
                       "material_state", "time_increments")
    optional_inputs = ("dofs_prev",)
    limitations = ("finite strain, current-configuration assembly",
                   "C3D8 (8-node hex) only",
                   "needs a deformation-gradient material (e.g. UMAT / crystal_plasticity)",
                   "consistent tangent approximate pending Abaqus AMATRX comparison")
    verification_tests = ("zero-field", "rigid-body", "patch", "FD-tangent",
                          "CP-kernel-bitwise")
    notes = "Home of the verified crystal-plasticity backend (one example backend)."
    required_inputs_by_mode = {
        "material-replay": ("coords", "connectivity", "dofs", "material",
                            "material_parameters", "material_state",
                            "solution_history", "time_increments")}
    optional_inputs_by_mode = {"material-replay": ("dofs_prev",)}
    material_interface_needed = True
    state_requirements = "history-dependent"
    tangent_support = "analytic (approx; pending Abaqus AMATRX comparison)"
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
            raise ValueError("solid_c3d8_finite_strain needs a material binding "
                             "(Mode 2). For Mode 1 use stress_driven_adapter.")
        nstate = getattr(binding, "n_state_vars", 0) or material.n_state_vars
        pts = k.ABAQUS_C3D8_GAUSS.points

        sigma_ip = np.zeros((8, 6))
        D_ip = np.zeros((8, 6, 6))
        state_new = np.zeros((8, nstate)) if nstate else None
        for ip in range(8):
            F0, F1 = _F_at(coords, u_e, u_prev, pts[ip])
            s_prev = (material_state[ip] if material_state is not None
                      else np.zeros(nstate))
            kin = {"F0": F0, "F1": F1,
                   "element": element_id, "ip": ip + 1}
            sigma, D, s_new, _ = material.evaluate(
                kin, s_prev, binding, time, dtime, fields, options)
            sigma_ip[ip] = sigma
            D_ip[ip] = D if D is not None else 0.0
            if state_new is not None:
                state_new[ip] = s_new

        r = k.element_internal_force_finite_strain(coords, u_e, sigma_ip)
        K = None
        if compute_tangent:
            K = k.element_tangent(coords, u_e, D_ip, sigma_ip, mode="finite")
        return r, K, state_new, {"formulation": self.name, "n_ip": 8}
