"""First-order OTILib parameter RHS at fixed real finite-strain geometry."""

from copy import copy

import numpy as np

from ..algebra.otilib_adapter import OtiContext
from ..core.assembler import Assembler
from ..core.sensitivity_package import AlgebraMetadata, SensitivityRHSResult


def parameter_slot(model, parameter):
    """Resolve an explicit material parameter; never silently return zero."""
    if "." not in parameter:
        raise ValueError("finite-strain parameter must be material.mu or material.lambda")
    material_name, key = parameter.rsplit(".", 1)
    binding = model.materials.get(material_name)
    names = getattr(getattr(binding, "material", None), "parameters", ())
    if key not in names:
        raise ValueError("unknown finite-strain material parameter %r" % parameter)
    return binding, names.index(key)


class FiniteStrainParameterRHS:
    name = "otilib-finite-strain-parameters"

    def evaluate_rhs(self, problem, real_solution, tangent, parameters, order=1):
        if order != 1:
            raise ValueError("finite-strain OTILib supports first-order material parameters only; "
                             "higher orders require displacement/geometry differentiation")
        if set(problem.model.element_formulation.values()) != {"solid_c3d8_finite_strain"}:
            raise ValueError("finite-strain OTILib requires a homogeneous C3D8 finite-strain model")
        params = list(parameters)
        if not params or len(set(params)) != len(params):
            raise ValueError("finite-strain OTILib requires distinct material parameters")
        context = OtiContext(len(params), 1)
        original = problem.model.materials
        seeded = {name: copy(binding) for name, binding in original.items()}
        for binding in seeded.values():
            binding.constants = list(binding.constants)
        for direction, parameter in enumerate(params, 1):
            binding, slot = parameter_slot(problem.model, parameter)
            name = parameter.rsplit(".", 1)[0]
            seeded[name].constants[slot] = context.seed(binding.constants[slot], direction)
        problem.model.materials = seeded
        try:
            residual, _, _ = Assembler(problem.model, problem.dof_manager, problem.forms).assemble(
                real_solution, compute_tangent=False)
        finally:
            problem.model.materials = original
        directions = context.order_directions(1)
        derivatives = np.array([[context.coeff(value, direction["exponents"])
                                 for direction in directions] for value in residual])
        result = SensitivityRHSResult(
            ndof=problem.dof_manager.ndof,
            algebra=AlgebraMetadata("otilib", n_bases=len(params), truncation_order=1),
            parameter_map={name: index + 1 for index, name in enumerate(params)}, max_order=1)
        result.set_residual_order(1, derivatives)
        result.diagnostics = {"provider": self.name, "hypercomplex_ready": True,
                              "scope": "first-order material parameters; fixed real geometry"}
        return result