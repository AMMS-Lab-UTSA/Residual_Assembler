"""Reproduce the bounded total-hyperelastic benchmark through the public API."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from residual_core import ResidualProblem
from residual_core.core.model import Model, Element
from residual_core.formulations import c3d8_kernel as kernel
from residual_core.materials.base import MaterialBinding
from residual_core.materials.neo_hookean import CompressibleNeoHookean


def reference_force(model, displacement, constants=(2.3, 4.1)):
    """Independent total-Lagrangian integral of P = dW/dF, not spatial B/D."""
    mu, lam = constants
    force = np.zeros(displacement.size)
    for element in model.elements.values():
        indices = np.array(element.connectivity) - 1
        coords = model.coords_of(element.connectivity)
        element_force = np.zeros((8, 3))
        for point in kernel.ABAQUS_C3D8_GAUSS.points:
            natural_gradient = kernel.shape_grad_natural(point)
            jacobian = coords.T @ natural_gradient
            gradient = np.linalg.solve(jacobian.T, natural_gradient.T).T
            deformation = np.eye(3) + displacement.reshape(-1, 3)[indices].T @ gradient
            inverse_transpose = np.linalg.inv(deformation).T
            piola = mu * (deformation - inverse_transpose) + lam * np.log(
                np.linalg.det(deformation)) * inverse_transpose
            element_force += gradient @ piola.T * np.linalg.det(jacobian)
        np.add.at(force.reshape(-1, 3), indices, element_force)
    return force


def make_problem():
    nodes = {}
    for layer in range(3):
        for corner, (height, depth) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
            nodes[4 * layer + corner + 1] = (
                layer + 0.08 * height * depth, height + 0.04 * layer * depth,
                depth + 0.03 * layer * height)
    model = Model(nodes=nodes)
    for index in range(2):
        offset = 4 * index
        connectivity = [offset + entry for entry in (1, 5, 6, 2, 4, 8, 7, 3)]
        model.elements[index + 1] = Element(index + 1, "C3D8", connectivity)
    model.element_formulation = {index: "solid_c3d8_finite_strain" for index in model.elements}
    model.element_material = {index: "solid" for index in model.elements}
    model.materials = {"solid": MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1], name="solid")}
    coordinates = np.array(list(nodes.values()))
    deformation = np.array([[1.3, 0.3, -0.1], [0.08, 0.85, 0.12], [0.02, -0.1, 1.12]])
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    displacement = (coordinates @ (rotation @ deformation - np.eye(3)).T).ravel()
    displacement[12:24] += np.array([0.02, -0.015, 0.01] * 4)
    for node in list(range(1, 5)) + list(range(9, 13)):
        for component in range(3):
            model.boundaries.append(SimpleNamespace(target=node, dof_start=component + 1,
                dof_end=component + 1, value=displacement[3 * (node - 1) + component], kind="fixed"))
    forces = reference_force(model, displacement)
    for node in range(5, 9):
        for component in range(3):
            model.cloads.append(SimpleNamespace(target=node, dof=component + 1,
                value=forces[3 * (node - 1) + component]))
    return ResidualProblem(model), displacement


def relative_error(actual, reference):
    return float(np.linalg.norm(actual - reference) / max(np.linalg.norm(reference), 1e-30))


def verify(out):
    from residual_core.algebra.otilib_adapter import otilib_status

    problem, solution = make_problem()
    output = Path(out)
    output.mkdir(parents=True, exist_ok=True)
    problem.save_neutral(str(output / "model.json"))
    config = {"mode": "material-replay", "options": {"solution": solution.tolist()}}
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    (output / "params.json").write_text(json.dumps({"mode": "material-replay",
        "parameters": ["solid.mu", "solid.lambda"], "backend": "otilib", "order": 1}, indent=2) + "\n")
    result = problem.result("material-replay", U=solution)
    expected = reference_force(problem.model, solution)
    actual_internal = result.R.copy()
    for load in problem.model.cloads:
        actual_internal[3 * (load.target - 1) + load.dof - 1] += load.value
    tangent_errors = []
    for step in (1e-4, 1e-5, 1e-6):
        numerical = np.column_stack([(problem.assemble("material-replay", U=solution + step * direction)
            - problem.assemble("material-replay", U=solution - step * direction)) / (2 * step)
            for direction in np.eye(solution.size)])
        tangent_errors.append({"step": step, "relative_error": relative_error(result.tangent.T, numerical)})
    initial = solution.copy()
    initial[result.free_mask] += 0.01
    converged = problem.solve_newton("material-replay", u0=initial)
    package = problem.sensitivity_package("material-replay", U=converged,
        parameters=["solid.mu", "solid.lambda"], generate_rhs=True, backend="otilib")
    sensitivity = package.solve(1)
    sensitivity_errors = []
    for step in (1e-3, 1e-4, 1e-5):
        numerical = problem.finite_difference_sensitivity("material-replay",
            ["solid.mu", "solid.lambda"], h=step, u0=converged)
        sensitivity_errors.append({"step": step, "relative_error": relative_error(sensitivity, numerical)})
    report = {"schema": "resasm-finite-strain-verification/1", "elements": 2, "ndof": solution.size,
        "law": "mu/2*(I1-3)-mu*log(J)+lambda/2*log(J)^2", "constants": [2.3, 4.1],
        "residual_reference": "independent reference-volume first-Piola quadrature",
        "residual_relative_error": relative_error(actual_internal, expected),
        "residual_norm": float(np.linalg.norm(result.R)),
        "free_residual_norm": result.free_residual_norm,
        "solution_relative_error": relative_error(converged, solution),
        "nodal_tangent_fd": tangent_errors, "oti_solution_fd": sensitivity_errors,
        "otilib": otilib_status(), "abaqus": "not run by this offline benchmark",
        "scope": "total isotropic hyperelasticity; first-order parameter sensitivity at real geometry"}
    report["passed"] = (report["residual_relative_error"] < 1e-12
        and report["free_residual_norm"] < 1e-11 and report["solution_relative_error"] < 1e-10
        and max(item["relative_error"] for item in tangent_errors) < 2e-8
        and max(item["relative_error"] for item in sensitivity_errors) < 2e-5)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(Path(__file__).parent / "verified"))
    args = parser.parse_args()
    report = verify(args.out)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["passed"] else 1)