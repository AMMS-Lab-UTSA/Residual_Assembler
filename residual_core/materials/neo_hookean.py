"""Compressible logarithmic neo-Hookean solid, real F and scalar-generic PROPS."""

import numpy as np

from .base import Material


class CompressibleNeoHookean(Material):
    """W = mu/2 (tr(F F.T)-3) - mu log(J) + lam/2 log(J)**2.

    DDSDDE:d = [mu (d b + b d) + lam tr(d) I]/J is the
    Jaumann rate of Kirchhoff stress divided by J, not a Cauchy rate.
    F0 is validated but unused: this is a total, isotropic, stateless law.
    All tensors are in fixed global axes; no DROT/history transport is needed.
    """

    name = "compressible_neo_hookean"
    kinematic_input = "deformation_gradient"
    stress_measure = "cauchy"
    tangent_measure = "kirchhoff_jaumann_over_j"
    response_kind = "isotropic_total_hyperelastic"
    input_variables = ("F0", "F1")
    parameters = ("mu", "lambda")
    supported_formulations = ("solid_c3d8_finite_strain",)
    limitations = ("real fixed geometry for parameter sensitivity",
                   "isotropic total hyperelasticity; no history or DROT update",
                   "mu > 0 and lambda >= 0; J > 0")

    def evaluate(self, kinematics, state_prev, binding, time, dtime, fields, options):
        if fields or (state_prev is not None and np.size(state_prev)):
            raise ValueError("compressible_neo_hookean supports no fields or history state")
        if len(binding.constants) != 2:
            raise ValueError("compressible_neo_hookean requires PROPS [mu, lambda]")
        mu, lam = binding.constants
        primal = [float(getattr(value, "real", value)) for value in (mu, lam)]
        if not np.all(np.isfinite(primal)) or primal[0] <= 0 or primal[1] < 0:
            raise ValueError("compressible_neo_hookean requires finite mu > 0, lambda >= 0")
        for name in ("F0", "F1"):
            deformation = np.asarray(kinematics[name])
            if (deformation.shape != (3, 3) or deformation.dtype.kind not in "fi"
                    or not np.all(np.isfinite(deformation))):
                raise ValueError("compressible_neo_hookean requires real finite 3x3 " + name)
            if np.linalg.det(deformation) <= 0:
                raise ValueError("compressible_neo_hookean requires det(%s) > 0" % name)
        deformation = np.asarray(kinematics["F1"])
        jacobian = np.linalg.det(deformation)
        left_cauchy_green = deformation @ deformation.T
        identity = np.eye(3)
        pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
        stress = np.array([(mu * (left_cauchy_green[row, column] - identity[row, column])
                            + lam * np.log(jacobian) * identity[row, column]) / jacobian
                           for row, column in pairs])
        columns = []
        for row, column in pairs:
            strain = np.zeros((3, 3))
            strain[row, column] = 1.0 if row == column else 0.5
            strain[column, row] = strain[row, column]
            rate = strain @ left_cauchy_green + left_cauchy_green @ strain
            columns.append([(mu * rate[first, second]
                             + lam * np.trace(strain) * identity[first, second]) / jacobian
                            for first, second in pairs])
        return stress, np.array(columns).T, np.zeros(0), {
            "material": self.name, "J": jacobian,
            "tangent_measure": self.tangent_measure, "history": "none"}