"""1-DOF nonlinear spring — the first verified nonlinear residual backend.

The residual sensitivity method needs a backend whose derivatives are known in
closed form to prove the R^(1) → ``T U^(1) = -R^(1)`` pathway end to end. This is
that backend: a single grounded scalar DOF with a cubic force law

    R(u, k, f) = k u^3 - f

so at equilibrium ``u = (f/k)^(1/3)``, the tangent is ``T = dR/du = 3 k u^2`` and
the parameter sensitivity RHS is ``-dR/dk = -u^3``, giving analytically

    du/dk = (-u^3) / (3 k u^2) = -u / (3k).

The element residual is written with plain arithmetic so it is **dual-safe**: when
the parameter ``k`` (or ``f``) is a ``Dual1`` number, the returned residual carries
``dR/dk`` in its imaginary part — that is how ``core/rhs_provider.py`` extracts
R^(1) without any hand-coded derivative.

This is a proof backend (not production): one node, one DOF, parameters read from
the element's ``properties`` (a section dict or MaterialBinding.section with keys
``k`` and ``f``).
"""

from __future__ import annotations

from .base import Formulation


def _prop(properties, key, default=None):
    """Read a scalar parameter from a dict or a MaterialBinding.section dict.
    Returns the value unchanged (may be a float OR a Dual1 for AD)."""
    if properties is None:
        return default
    if isinstance(properties, dict):
        if key in properties:
            return properties[key]
    sec = getattr(properties, "section", None)
    if isinstance(sec, dict) and key in sec:
        return sec[key]
    return default


def _real_of(x):
    """Real value of a scalar that may be a plain float, a Dual1, or a
    hypercomplex OTI number — for diagnostics only (never forces a hypercomplex
    scalar through float())."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return float(getattr(x, "real", 0.0))


class NonlinearSpring1(Formulation):
    name = "nonlinear_spring1"
    element_types = ("SPRING1", "NLSPRING1", "NLS1")
    dof_types = ("U",)                         # one scalar DOF per node
    sign_convention = "residual"
    verification_levels = (0, 1, 4)
    supported_modes = ("formulation",)
    required_inputs = ("coords", "connectivity", "dofs", "section: k, f")
    optional_inputs = ()
    limitations = ("1 node, 1 scalar DOF (grounded cubic spring)",
                   "proof backend for the sensitivity pathway, not production")
    verification_tests = ("zero-field", "analytic-tangent", "sensitivity-RHS")
    notes = ("R = k u^3 - f. Dual-safe residual: with a Dual1 parameter the "
             "imaginary part is dR/dparam (used to build R^(1)).")
    required_inputs_by_mode = {"formulation": ("coords", "connectivity", "dofs",
                                               "section: k, f")}
    optional_inputs_by_mode = {"formulation": ()}
    material_interface_needed = False
    state_requirements = "none"
    tangent_support = "analytic"
    verification_status = "verified"

    # parameters this backend exposes for sensitivity analysis (section keys)
    sensitivity_parameters = ("k", "f")

    def n_nodes(self, element_type):
        return 1

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        compute_tangent = options.get("compute_tangent", True)
        u = dofs[0]                                 # scalar DOF (float or np.float)
        k = _prop(properties, "k", default=1.0)
        f = _prop(properties, "f", default=0.0)

        # residual: plain arithmetic -> dual-safe (k/f may be Dual1)
        r = k * u ** 3 - f
        tangent = None
        if compute_tangent:
            tangent = [[3.0 * k * u * u]]           # dR/du = 3 k u^2
        return [r], tangent, material_state, {
            "formulation": self.name, "u": _real_of(u)}
