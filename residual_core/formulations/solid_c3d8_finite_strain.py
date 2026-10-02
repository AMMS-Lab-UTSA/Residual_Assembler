"""Bounded C3D8 total isotropic hyperelasticity in the current configuration.

For L = delta F F^-1 = d + w and Abaqus D = tau^Jaumann/J,
delta sigma = D:d + w sigma - sigma w - sigma tr(d).
Linearizing spatial gradients and volume gives B.T c B + Kgeo with
c:d = D:d - d sigma - sigma d. See docs/evidence/recovery_finite.md.
History-dependent rotating materials are deliberately unsupported.

Integration: ``options["integration"]`` is ``"full"`` (default; 2x2x2 Gauss,
plain F and B -- what every published number of this formulation used) or
``"mean_dilatation"`` (Abaqus C3D8: the material sees Fbar with the element's
mean volume change and the force is the B-bar force with weight w0 Jbar,
residual_core.formulations.c3d8_nlgeom; reproduces Abaqus 2021 reactions on a
distorted element to 4.7e-8). With mean dilatation the tangent is the exact
linearisation through Fbar: for l = dFbar Fbar^-1 = d + w,
delta sigma = D:d - sigma tr(d) + w sigma - sigma w at every point, plus the
fixed-stress derivative of the B-bar force.
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
    limitations = ("C3D8 full integration, real positive-J geometry",
                   "isotropic total hyperelastic materials only; no history transport",
                   "requires Cauchy stress and Kirchhoff-Jaumann/J tangent")
    verification_tests = ("test_finite_strain_hyperelastic",)
    notes = "Bounded total hyperelastic path; no generic finite-strain UMAT claim."
    required_inputs_by_mode = {
        "material-replay": ("coords", "connectivity", "dofs", "material",
                            "material_parameters", "material_state",
                            "solution_history", "time_increments")}
    optional_inputs_by_mode = {"material-replay": ("dofs_prev",)}
    material_interface_needed = True
    #: the material kinematic input this element drives; with several C3D8
    #: replay backends, the element's material chooses by it
    accepted_kinematic_inputs = ("deformation_gradient",)
    state_requirements = "stateless total hyperelasticity; F0/F1 explicit"
    tangent_support = "analytic, exact weak linearization for declared measures"
    verification_status = "bounded-hyperelastic"

    def n_nodes(self, element_type):
        return 8

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        solution_state = solution_state or {}
        compute_tangent = options.get("compute_tangent", True)
        u_e = np.asarray(dofs)
        u_prev = np.asarray(solution_state.get("dofs_prev", np.zeros_like(u_e)))
        if any(value.dtype.kind not in "fi" or not np.all(np.isfinite(value))
               or value.size != 24 for value in (u_e, u_prev)):
            raise ValueError("finite-strain C3D8 requires 24 real finite DOFs; "
                             "OTI is supported for material parameters at fixed real geometry")
        coords = np.asarray(coords)
        if (element_type != "C3D8" or coords.shape != (8, 3)
            or coords.dtype.kind not in "fi" or not np.all(np.isfinite(coords))):
            raise ValueError("finite-strain formulation requires C3D8 with finite (8,3) coordinates")

        binding = properties
        material = getattr(binding, "material", None)
        if material is None:
            raise ValueError("solid_c3d8_finite_strain needs a material binding "
                             "(Mode 2). For Mode 1 use stress_driven_adapter.")
        expected = {"kinematic_input": "deformation_gradient",
                    "stress_measure": "cauchy",
                    "tangent_measure": "kirchhoff_jaumann_over_j",
                    "response_kind": "isotropic_total_hyperelastic"}
        for attribute, required in expected.items():
            actual = getattr(material, attribute, None)
            if actual != required:
                raise ValueError("solid_c3d8_finite_strain: material %s has %s=%r; "
                                 "requires %r. Only bounded isotropic total hyperelasticity "
                                 "is available; generic rotating/history UMATs are unsupported."
                                 % (material.name, attribute, actual, required))
        if (binding.n_state_vars or material.n_state_vars
                or (material_state is not None and np.size(material_state))):
            raise ValueError("solid_c3d8_finite_strain does not support history state or DROT transport")
        pts = k.ABAQUS_C3D8_GAUSS.points
        integration = options.get("integration", "full")
        if integration in ("mean_dilatation", "selective_reduced"):
            return self._eval_mean_dilatation(element_id, coords, u_e, u_prev, material,
                                              binding, time, dtime, fields, options,
                                              compute_tangent)
        if integration != "full":
            raise ValueError("solid_c3d8_finite_strain: integration must be 'full' or "
                             "'mean_dilatation', got %r" % (integration,))

        sigma_ip = []
        D_ip = []
        for ip in range(8):
            for configuration in (coords, coords + u_e.reshape(8, 3),
                                  coords + u_prev.reshape(8, 3)):
                jacobian = configuration.T @ k.shape_grad_natural(pts[ip])
                if np.linalg.det(jacobian) <= 0:
                    raise ValueError("finite-strain C3D8 requires positive reference/current/previous "
                                     "Jacobian at element %s IP %s" % (element_id, ip + 1))
            F0, F1 = _F_at(coords, u_e, u_prev, pts[ip])
            kin = {"F0": F0, "F1": F1,
                   "element": element_id, "ip": ip + 1}
            sigma, D, s_new, _ = material.evaluate(
                kin, np.zeros(0), binding, time, dtime, fields, options)
            if np.asarray(sigma).shape != (6,) or np.size(s_new):
                raise ValueError("finite-strain material must return stress (6,) and no history state")
            sigma_ip.append(sigma)
            if compute_tangent:
                if D is None or np.asarray(D).shape != (6, 6):
                    raise ValueError("finite-strain tangent requires material DDSDDE (6,6)")
                D_ip.append(k.kirchhoff_jaumann_to_spatial(D, sigma))

        r = k.element_internal_force_finite_strain(coords, u_e, sigma_ip)
        K = None
        if compute_tangent:
            K = k.element_tangent(coords, u_e, D_ip, sigma_ip, mode="finite")
        return r, K, None, {"formulation": self.name, "n_ip": 8,
                    "tangent_measure": "exact_weak_linearization",
                    "history": "none", "rotation": "global isotropic total response"}

    def _eval_mean_dilatation(self, element_id, coords, u_e, u_prev, material, binding,
                              time, dtime, fields, options, compute_tangent):
        """Mean-dilatation path. Geometry (coords, u, u_prev) is real -- refused
        otherwise in ``eval_element`` -- so Fbar, its inverse and its
        directional derivative are real. The stress and DDSDDE take the dtype
        of the material constants: with Dual1/OTI constants they are object
        arrays, never cast to float, and the force and tangent are assembled
        from them by real-coefficient linear combinations (C3D8Nlgeom.force /
        dforce_fixed_stress), so their parameter derivatives survive."""
        from .c3d8_nlgeom import C3D8Nlgeom, voigt_eng
        el = C3D8Nlgeom(np.asarray(coords, dtype=float)[None], "mean_dilatation")
        U1 = u_e.reshape(1, 8, 3)
        U0 = u_prev.reshape(1, 8, 3)
        F1 = el.Fbar(U1)[0]
        F0 = el.Fbar(U0)[0]
        sigma_ip, D_ip = [], []
        for ip in range(8):
            if np.linalg.det(F1[ip]) <= 0 or np.linalg.det(F0[ip]) <= 0:
                raise ValueError("finite-strain C3D8 requires positive Jacobians at element %s "
                                 "IP %s" % (element_id, ip + 1))
            kin = {"F0": F0[ip], "F1": F1[ip], "element": element_id, "ip": ip + 1}
            sigma, D, s_new, _ = material.evaluate(kin, np.zeros(0), binding, time, dtime,
                                                   fields, options)
            if np.asarray(sigma).shape != (6,) or np.size(s_new):
                raise ValueError("finite-strain material must return stress (6,) and no history state")
            sigma_ip.append(np.asarray(sigma))
            D_ip.append(None if D is None else np.asarray(D))
        live = any(a.dtype == object for a in sigma_ip + [D for D in D_ip if D is not None])
        number = object if live else np.result_type(*sigma_ip, float)
        sigma_ip = np.array(sigma_ip, dtype=number)
        r = el.force(U1, sigma_ip[None])[0]
        K = None
        if compute_tangent:
            if any(D is None or D.shape != (6, 6) for D in D_ip):
                raise ValueError("finite-strain tangent requires material DDSDDE (6,6)")
            E = el.unit_directions(1)
            dF = el.dFbar(U1, E)[0]                                   # (8, 24, 3, 3)
            dsig = []
            for q in range(8):
                S = k.voigt_to_tensor(sigma_ip[q])
                L = dF[q] @ np.linalg.inv(F1[q])                         # (24, 3, 3), real
                d = 0.5 * (L + np.swapaxes(L, 1, 2))
                w = 0.5 * (L - np.swapaxes(L, 1, 2))
                dS = (w @ S - S @ w) - S[None] * np.trace(d, axis1=1, axis2=2)[:, None, None]
                dsig.append(voigt_eng(d) @ np.asarray(D_ip[q], dtype=number).T
                            + np.stack([dS[:, 0, 0], dS[:, 1, 1], dS[:, 2, 2], dS[:, 0, 1],
                                        dS[:, 0, 2], dS[:, 1, 2]], axis=1))
            dsig = np.array(dsig, dtype=number)                       # (8, 24, 6)
            K = (el.dforce_fixed_stress(U1, sigma_ip[None], E) + el.force(U1, dsig[None]))[0]
        return r, K, None, {"formulation": self.name, "n_ip": 8,
                            "integration": "mean_dilatation",
                            "tangent_measure": "exact_weak_linearization_through_fbar",
                            "history": "none", "rotation": "global isotropic total response"}
