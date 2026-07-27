"""2-node bar / truss backend — a small, non-solid proof of model agnosticism.

This exists to demonstrate that the framework is **not** secretly hard-coded to
C3D8 + crystal plasticity. It is a classic 3D small-strain, linear-elastic bar
(Euler truss) with axial stiffness only:

    F_internal,e = K_e u_e ,     K_e = (E A / L) [[ n n^T , -n n^T ],
                                                  [-n n^T ,  n n^T ]]

where ``n`` is the unit vector along the element and ``L`` its length. The tangent
is constant (linear element), so the Level-4 finite-difference check is exact.

DOFs: (UX, UY, UZ) per node — translational only, *no* rotations, *no*
integration-point stress, *no* material state. Section/material data (E, A) are
read from the element's ``properties`` slot (a MaterialBinding.section, a plain
dict, or the binding constants ``[E, A]``).

This is a proof backend, not a production truss element.
"""

from __future__ import annotations

import numpy as np

from .base import Formulation


def _prop(properties, key, aliases=(), default=None):
    """Pull a scalar section/material property from a flexible ``properties``
    object: a dict, a MaterialBinding with ``.section`` dict, or ``.constants``.

    The value is returned **unchanged** (may be a plain float OR a hypercomplex
    seed — Dual1 / OTI number — for sensitivity), so a seeded design parameter
    (e.g. E) flows through the residual and carries its derivative."""
    keys = (key,) + tuple(aliases)
    if properties is None:
        return default
    if isinstance(properties, dict):
        for k in keys:
            if k in properties:
                return properties[k]
    sec = getattr(properties, "section", None)
    if isinstance(sec, dict):
        for k in keys:
            if k in sec:
                return sec[k]
    return default


def _real_of(x):
    """Real value of a scalar that may be a plain float, a Dual1, or an OTI
    number — for diagnostics only (never forces a hypercomplex scalar through
    float())."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return float(getattr(x, "real", 0.0))


class Truss2(Formulation):
    name = "truss2"
    element_types = ("T3D2", "T2D2", "TRUSS2", "BAR2")
    dof_types = ("UX", "UY", "UZ")
    sign_convention = "residual"       # returns +F_internal contribution
    verification_levels = (0, 1, 4)
    supported_modes = ("formulation",)
    required_inputs = ("coords", "connectivity", "dofs", "section: E, A")
    optional_inputs = ()
    limitations = ("small-strain linear-elastic axial bar (Euler truss)",
                   "2 nodes, axial stiffness only, no bending")
    verification_tests = ("zero-field", "axial-force", "FD-tangent")
    notes = "Proof backend: shows the core is not C3D8/CP-specific."
    required_inputs_by_mode = {
        "formulation": ("coords", "connectivity", "dofs", "section: E, A")}
    optional_inputs_by_mode = {"formulation": ()}
    material_interface_needed = False
    state_requirements = "none"
    tangent_support = "analytic"
    verification_status = "implemented/simple"

    # design parameters this backend exposes for sensitivity analysis (section keys)
    sensitivity_parameters = ("E", "A")

    def n_nodes(self, element_type):
        return 2

    def stiffness(self, coords, E, A):
        coords = np.asarray(coords, dtype=float)
        d = coords[1] - coords[0]
        L = float(np.linalg.norm(d))
        if L <= 0.0:
            raise ValueError("truss2 element has zero length")
        n = d / L
        nn = np.outer(n, n)                       # (3,3) float geometry
        pattern = np.block([[nn, -nn], [-nn, nn]])   # (6,6) float
        k_ax = E * A / L                          # float OR hypercomplex (E seeded)
        K = k_ax * pattern                        # object array if k_ax is hyper
        return K, L, n

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        compute_tangent = options.get("compute_tangent", True)
        # keep u_e real (float) when possible; fall back to object dtype so a
        # hypercomplex u* (OTI order>=2 carries injected lower-order coeffs) works.
        try:
            u_e = np.asarray(dofs, dtype=float).reshape(6)
        except (TypeError, ValueError):
            u_e = np.asarray(list(dofs), dtype=object).reshape(6)

        E = _prop(properties, "E", ("youngs", "young", "E_modulus"), default=None)
        A = _prop(properties, "A", ("area", "section_area"), default=None)
        if E is None or A is None:
            # allow binding constants [E, A] as a last resort
            consts = list(getattr(properties, "constants", []) or [])
            if E is None and len(consts) >= 1:
                E = float(consts[0])
            if A is None and len(consts) >= 2:
                A = float(consts[1])
        if E is None or A is None:
            raise ValueError(
                "truss2 needs section properties E and A for element %r "
                "(provide via MaterialBinding.section={'E':..,'A':..} or a dict)."
                % (element_id,))

        K, L, n = self.stiffness(coords, E, A)
        r = K @ u_e                                    # F_internal (linear)
        tangent = K if compute_tangent else None
        # report the axial force for diagnostics/stress-driven cross-checks
        axial = _real_of(E) * _real_of(A) / L * float(
            n @ (np.asarray([_real_of(v) for v in u_e[3:]])
                 - np.asarray([_real_of(v) for v in u_e[:3]])))
        return r, tangent, material_state, {
            "formulation": self.name, "length": L, "axial_force": axial}
