"""PRIVATE partner residual provider (element level) — nonlinear bar chain.

Two nonlinear axial bar elements in series (nodes 0-1-2), node 0 fixed and node 2
given a prescribed displacement. Each element has its own stiffness parameter
(``k1``, ``k2``). Demonstrates the *element* provider level: the kit scatters
element residuals it never has to understand.

Element axial force (cubic law):  N = k * (u_b - u_a)^3
Element residual:                 [-N, +N]        (dual-safe arithmetic)
Element tangent:                  3 k d^2 * [[1,-1],[-1,1]]
"""

import numpy as np

from partner_kit.python.residual_provider import ElementResidualProvider


class PrivateBar(ElementResidualProvider):
    name = "private_bar"
    parameters = ("k1", "k2")
    ndof = 3                                   # u0, u1, u2

    def __init__(self, k1=1.0, k2=1.0, u2=2.0):
        self._k = {"k1": float(k1), "k2": float(k2)}
        self._u2 = float(u2)                   # prescribed tip displacement
        # element -> (nodes, parameter name)
        self._elems = {1: ([0, 1], "k1"), 2: ([1, 2], "k2")}

    def parameter_values(self):
        return dict(self._k)

    def free_mask(self):
        return np.array([False, True, False])  # node 0 fixed, node 2 prescribed

    def elements(self):
        return self._elems.keys()

    def element_dof_map(self, element_id):
        return self._elems[element_id][0]

    def eval_element_residual(self, element_id, element_dofs, element_state,
                              parameters, time, dtime):
        nodes, pname = self._elems[element_id]
        k = parameters[pname]
        ua, ub = element_dofs[0], element_dofs[1]
        d = ub - ua
        N = k * d ** 3
        r = [-N, N]
        # tangent (real branch only; parameters real here)
        try:
            kd = 3.0 * float(k) * float(d) ** 2
            tan = np.array([[kd, -kd], [-kd, kd]])
        except (TypeError, ValueError):
            tan = None                          # dual parameters -> skip tangent
        return r, tan, element_state
