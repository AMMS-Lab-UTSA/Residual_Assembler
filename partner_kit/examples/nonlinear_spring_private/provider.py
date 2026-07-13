"""PRIVATE partner residual provider — never leaves the partner machine.

Toy stand-in for a proprietary model: a 1-DOF cubic spring R(u,k,f) = k u^3 - f.
The residual is written scalar-generically, so the kit can pass a hypercomplex
``k`` and read dR/dk from the imaginary part — the partner shares no source, mesh,
or material model with us.
"""

import numpy as np

from partner_kit.python.residual_provider import GlobalResidualProvider


class PrivateSpring(GlobalResidualProvider):
    name = "private_spring"
    parameters = ("k", "f")
    ndof = 1

    def __init__(self, k=2.0, f=16.0):
        self._k = float(k)
        self._f = float(f)

    def parameter_values(self):
        return {"k": self._k, "f": self._f}

    def free_mask(self):
        return np.array([True])

    def eval_global_residual(self, global_solution, state, parameters, time, dtime):
        u = global_solution[0]
        k = parameters["k"]
        f = parameters["f"]
        return [k * u ** 3 - f]                 # scalar-generic (float or Dual1)

    def get_tangent(self, global_solution, state, parameters, time, dtime):
        u = global_solution[0]
        k = parameters["k"]
        return np.array([[3.0 * k * u * u]])    # dR/du = 3 k u^2
