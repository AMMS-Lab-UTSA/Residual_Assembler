"""2-node nonlinear axial bar — a small coupled FE backend for OTI verification.

Companion to ``nonlinear_spring1`` but with TWO coupled DOFs per element, so a
chain of these gives a genuine multi-element FE residual (not decoupled). Used to
prove the OTILib sensitivity path at the finite-element level, not just the
algebraic spring level.

Cubic axial law between the element's two nodal DOFs (u_a, u_b):

    d = u_b - u_a
    N = k * d^3                         (axial internal force)
    element residual = [-N, +N]         (dual/OTI-safe: plain arithmetic)
    element tangent  = 3 k d^2 [[ 1,-1],
                                [-1, 1]]

Parameter ``k`` is read from the element ``properties`` (section dict / binding).
Proof backend, not production.
"""

from __future__ import annotations

from .base import Formulation


def _prop(properties, key, default=None):
    if properties is None:
        return default
    if isinstance(properties, dict) and key in properties:
        return properties[key]
    sec = getattr(properties, "section", None)
    if isinstance(sec, dict) and key in sec:
        return sec[key]
    return default


class NonlinearBar1(Formulation):
    oti_differentiable = True   # pure generic arithmetic: OTI scalars survive
    name = "nonlinear_bar1"
    element_types = ("NLBAR2", "BAR1", "NLB2")
    dof_types = ("U",)                         # one axial DOF per node
    sign_convention = "residual"
    verification_levels = (0, 1, 4)
    supported_modes = ("formulation",)
    required_inputs = ("coords", "connectivity", "dofs", "section: k")
    optional_inputs = ()
    limitations = ("2 nodes, cubic axial law", "proof backend for OTI FE path")
    verification_tests = ("zero-field", "analytic-tangent", "oti-sensitivity")
    notes = "N = k (u_b - u_a)^3; dual/OTI-safe residual for multivariate sensitivity."
    required_inputs_by_mode = {"formulation": ("coords", "connectivity", "dofs",
                                               "section: k")}
    optional_inputs_by_mode = {"formulation": ()}
    material_interface_needed = False
    state_requirements = "none"
    tangent_support = "analytic"
    verification_status = "verified"
    sensitivity_parameters = ("k",)

    def n_nodes(self, element_type):
        return 2

    def eval_element(self, element_id, element_type, coords, dofs,
                     solution_state, material_state, properties,
                     time, dtime, fields, options):
        options = options or {}
        compute_tangent = options.get("compute_tangent", True)
        ua, ub = dofs[0], dofs[1]
        k = _prop(properties, "k", default=1.0)
        d = ub - ua
        N = k * d ** 3
        r = [-N, N]
        tangent = None
        if compute_tangent:
            # real branch only (parameters real when a tangent is requested)
            kd = 3.0 * float(k) * float(d) ** 2
            tangent = [[kd, -kd], [-kd, kd]]
        return r, tangent, material_state, {"formulation": self.name}
