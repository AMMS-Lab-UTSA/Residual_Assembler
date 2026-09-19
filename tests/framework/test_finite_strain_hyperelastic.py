"""Independent residual and energy checks for the bounded finite-strain path."""

import numpy as np
import pytest

from residual_core.formulations import c3d8_kernel as kernel
from residual_core.formulations.solid_c3d8_finite_strain import SolidC3D8FiniteStrain
from residual_core.materials.base import MaterialBinding
from residual_core.materials.neo_hookean import CompressibleNeoHookean


def evaluate(coords, displacement, constants=(2.3, 4.1), previous=None, material=None):
    binding = MaterialBinding(material or CompressibleNeoHookean(), list(constants))
    state = {} if previous is None else {"dofs_prev": previous}
    return SolidC3D8FiniteStrain().eval_element(
        1, "C3D8", coords, displacement, state, None, binding,
        (0, 1), 1, None, {"compute_tangent": True})[:2]


def distorted_state():
    coords = kernel.unit_cube_Xe().copy()
    coords[6] += [0.15, -0.08, 0.12]
    coords[1] += [0.04, 0.02, -0.05]
    deformation = np.array([[1.35, 0.32, -0.12], [0.1, 0.82, 0.2],
                            [0.05, -0.08, 1.16]])
    displacement = coords @ (deformation - np.eye(3)).T
    displacement[5] += [0.03, -0.04, 0.02]
    return coords, displacement.ravel()


def independent_stress(coords, displacement, mu=2.3, lam=4.1):
    stresses = []
    for point in kernel.ABAQUS_C3D8_GAUSS.points:
        _, _, gradient = kernel._b_from_coords(coords, point)
        deformation = np.eye(3) + displacement.reshape(8, 3).T @ gradient
        jacobian = np.linalg.det(deformation)
        stress = (mu * (deformation @ deformation.T - np.eye(3))
                  + lam * np.log(jacobian) * np.eye(3)) / jacobian
        stresses.append(stress[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]])
    return np.array(stresses)


def test_kirchhoff_jaumann_correction_against_residual_fd():
    coords, displacement = distorted_state()
    stresses = independent_stress(coords, displacement)
    moduli = []
    mu, lam = 2.3, 4.1
    for point, stress in zip(kernel.ABAQUS_C3D8_GAUSS.points, stresses):
        _, _, gradient = kernel._b_from_coords(coords, point)
        deformation = np.eye(3) + displacement.reshape(8, 3).T @ gradient
        jacobian = np.linalg.det(deformation)
        spatial = np.diag([2, 2, 2, 1, 1, 1]) * (
            mu - lam * np.log(jacobian)) / jacobian
        spatial[:3, :3] += lam / jacobian
        material = CompressibleNeoHookean()
        _, ddsdde, _, _ = material.evaluate({"F0": np.eye(3), "F1": deformation},
            np.zeros(0), MaterialBinding(material, [mu, lam]), (0, 1), 1, None, {})
        converted = kernel.kirchhoff_jaumann_to_spatial(ddsdde, stress)
        np.testing.assert_allclose(converted, spatial, atol=1e-14)
        moduli.append(converted)
    tangent = kernel.element_tangent(coords, displacement, moduli, stresses)
    errors = []
    for step in (1e-4, 1e-5, 1e-6):
        columns = []
        for direction in np.eye(24):
            plus = displacement + step * direction
            minus = displacement - step * direction
            columns.append((kernel.element_internal_force_finite_strain(
                coords, plus, independent_stress(coords, plus))
                - kernel.element_internal_force_finite_strain(
                    coords, minus, independent_stress(coords, minus))) / (2 * step))
        finite_difference = np.array(columns).T
        errors.append(np.linalg.norm(tangent - finite_difference) /
                      np.linalg.norm(finite_difference))
    assert max(errors) < 2e-8, errors


def test_material_formulation_fd_sweep():
    coords, displacement = distorted_state()
    force, tangent = evaluate(coords, displacement)
    np.testing.assert_allclose(force, kernel.element_internal_force_finite_strain(
        coords, displacement, independent_stress(coords, displacement)), atol=1e-13)
    np.testing.assert_allclose(tangent, tangent.T, atol=1e-13)
    errors = []
    for step in (1e-4, 1e-5, 1e-6):
        finite_difference = np.column_stack([
            (evaluate(coords, displacement + step * direction)[0]
             - evaluate(coords, displacement - step * direction)[0]) / (2 * step)
            for direction in np.eye(24)])
        errors.append(np.linalg.norm(tangent - finite_difference) / np.linalg.norm(tangent))
    assert max(errors) < 2e-8, errors


def test_rotation_covariance_and_history_independence():
    coords, displacement = distorted_state()
    angle = 1.1
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotated = (coords + displacement.reshape(8, 3)) @ rotation.T - coords
    force, tangent = evaluate(coords, displacement)
    rotated_force, rotated_tangent = evaluate(coords, rotated.ravel(), previous=displacement * 0.4)
    transform = np.kron(np.eye(8), rotation)
    np.testing.assert_allclose(rotated_force, transform @ force, atol=1e-13)
    np.testing.assert_allclose(np.linalg.norm(rotated_force), np.linalg.norm(force), atol=1e-13)
    np.testing.assert_allclose(rotated_tangent, transform @ tangent @ transform.T, atol=1e-13)
    rigid = coords @ rotation.T + [3, -2, 4] - coords
    np.testing.assert_allclose(evaluate(coords, rigid.ravel())[0], 0, atol=1e-13)


@pytest.mark.parametrize("attribute,value", [("kinematic_input", "small_strain"),
    ("stress_measure", "pk2"), ("tangent_measure", "ddsdde"),
    ("response_kind", "plastic"), ("n_state_vars", 1)])
def test_reject_incompatible_material(attribute, value):
    material = CompressibleNeoHookean()
    setattr(material, attribute, value)
    with pytest.raises(ValueError, match="requires|history"):
        evaluate(*distorted_state(), material=material)


def test_small_strain_limit():
    coords = kernel.unit_cube_Xe()
    displacement = coords @ np.array([[0.2, 0.3, 0], [0, -0.1, 0.2], [0, 0, 0.1]]).T
    mu, lam = 2.3, 4.1
    linear_modulus = np.diag([2, 2, 2, 1, 1, 1]) * mu
    linear_modulus[:3, :3] += lam
    linear_tangent = kernel.element_tangent(coords, None, linear_modulus, mode="small")
    zero_force, zero_tangent = evaluate(coords, np.zeros(24))
    np.testing.assert_allclose(zero_force, 0, atol=1e-14)
    np.testing.assert_allclose(zero_tangent, linear_tangent, atol=1e-14)
    force, _ = evaluate(coords, displacement.ravel() * 1e-7)
    np.testing.assert_allclose(force / 1e-7, linear_tangent @ displacement.ravel(),
                               rtol=2e-6, atol=1e-7)


def finite_problem():
    from residual_core import ResidualProblem
    from residual_core.core.model import Model, Element

    coords, displacement = distorted_state()
    model = Model(nodes={index + 1: tuple(point) for index, point in enumerate(coords)})
    model.elements = {1: Element(1, "C3D8", list(range(1, 9)))}
    model.element_formulation = {1: "solid_c3d8_finite_strain"}
    model.element_material = {1: "solid"}
    model.materials = {"solid": MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1], name="solid")}
    return ResidualProblem(model), displacement


def test_public_replay_and_neutral_roundtrip(tmp_path):
    from residual_core import ResidualProblem
    from residual_core.ui import cli

    problem, displacement = finite_problem()
    expected = evaluate(*distorted_state())
    actual = problem.assemble("material-replay", U=displacement, compute_tangent=True)
    for result, reference in zip(actual, expected):
        np.testing.assert_allclose(result, reference, atol=1e-13)
    path = tmp_path / "model.json"
    problem.save_neutral(str(path))
    loaded = ResidualProblem.from_neutral(str(path))
    for result, reference in zip(loaded.assemble("material-replay", U=displacement,
                                               compute_tangent=True), expected):
        np.testing.assert_allclose(result, reference, atol=1e-13)
    assert cli.main(["assemble", str(path), "--mode", "material-replay", "--tangent"]) == 0


def test_cli_incompatible_material_is_diagnostic(tmp_path, capsys):
    from residual_core.materials.elastic_adapter import IsotropicElastic
    from residual_core.ui import cli

    problem, _ = finite_problem()
    problem.model.materials["solid"].material = IsotropicElastic()
    path = tmp_path / "model.json"
    problem.save_neutral(str(path))
    assert cli.main(["assemble", str(path), "--mode", "material-replay"]) == 2
    assert "kinematic_input" in capsys.readouterr().err


def test_genuine_otilib_material_parameter_rhs_and_tangent():
    from residual_core.algebra import otilib_adapter as algebra

    if not algebra.otilib_available():
        import os
        assert os.environ.get("RUN_OTILIB_TESTS") != "1", algebra.otilib_status()
        pytest.skip("genuine OTILib unavailable")
    problem, displacement = finite_problem()
    package = problem.sensitivity_package(
        "material-replay", parameters=["solid.mu", "solid.lambda"], U=displacement,
        generate_rhs=True, backend="otilib")
    context = algebra.OtiContext(2, 1)
    _, tangent = evaluate(*distorted_state(), constants=(context.seed(2.3, 1),
                                                        context.seed(4.1, 2)))
    for parameter, exponents in enumerate(((1, 0), (0, 1))):
        errors = []
        for step in (1e-3, 1e-4, 1e-5):
            plus = np.array([2.3, 4.1])
            minus = plus.copy()
            plus[parameter] += step
            minus[parameter] -= step
            force_plus, tangent_plus = evaluate(*distorted_state(), constants=plus)
            force_minus, tangent_minus = evaluate(*distorted_state(), constants=minus)
            reference = (force_plus - force_minus) / (2 * step)
            errors.append(np.linalg.norm(package.sensitivity.R(1)[:, parameter] - reference)
                          / np.linalg.norm(reference))
            tangent_derivative = np.array([[context.deriv(value, exponents) for value in row]
                                           for row in tangent])
            np.testing.assert_allclose(tangent_derivative, (tangent_plus - tangent_minus) / (2 * step),
                                       rtol=1e-7, atol=1e-9)
        assert max(errors) < 1e-8, errors
    assert problem.model.materials["solid"].constants == [2.3, 4.1]
    with pytest.raises(ValueError, match="first-order"):
        problem.sensitivity_package("material-replay", parameters=["solid.mu"], U=displacement,
                                    generate_rhs=True, max_order=2)


def test_multi_element_benchmark_and_cli(tmp_path):
    import runpy
    from pathlib import Path
    from residual_core.algebra import otilib_adapter as algebra
    from residual_core.ui import cli

    if not algebra.otilib_available():
        import os
        assert os.environ.get("RUN_OTILIB_TESTS") != "1", algebra.otilib_status()
        pytest.skip("genuine OTILib unavailable")
    benchmark = Path(__file__).resolve().parents[2] / "examples/finite_strain_c3d8/benchmark.py"
    report = runpy.run_path(str(benchmark))["verify"](tmp_path)
    assert report["passed"], report
    assert cli.main(["--config", str(tmp_path / "config.json"), "assemble",
                     str(tmp_path / "model.json"), "--mode", "material-replay", "--tangent"]) == 0
    assert cli.main(["--config", str(tmp_path / "config.json"),
                     "sensitivity", str(tmp_path / "model.json"),
                     "--params", str(tmp_path / "params.json")]) == 0


def test_dual1_refuses_the_finite_strain_model_instead_of_printing_zeros(capsys):
    """Regression: `--backend dual1` never seeded mu and lambda (they are
    material constants, not section entries), printed d^1/e1 = d^1/e2 = 0 with
    'FD ..., rel 1.00e+00' beside them, and exited 0."""
    from pathlib import Path
    from residual_core.ui import cli

    verified = Path(__file__).resolve().parents[2] / "examples/finite_strain_c3d8/verified"
    code = cli.main(["--config", str(verified / "config.json"), "sensitivity",
                     str(verified / "model.json"), "--params", str(verified / "params.json"),
                     "--backend", "dual1"])
    captured = capsys.readouterr()
    assert code == 2
    assert "backend='dual1' cannot differentiate the finite-strain C3D8" in captured.err
    assert "d^1/" not in captured.out


def test_gui_uses_finite_cli_backend(tmp_path):
    import json
    pytest.importorskip("streamlit")
    from residual_core.app.streamlit_app import _run

    problem, solution = finite_problem()
    model = tmp_path / "model.json"
    config = tmp_path / "config.json"
    problem.save_neutral(str(model))
    config.write_text(json.dumps({"options": {"solution": solution.tolist()}}))
    output = tmp_path / "residual.npy"
    result = _run(["--config", str(config), "assemble", str(model),
                   "--mode", "material-replay", "--tangent", "--out", str(output)])
    assert result.code == 0, result.stderr
    np.testing.assert_allclose(np.load(output), problem.assemble("material-replay", U=solution), atol=1e-14)


@pytest.mark.parametrize("which", ["reference", "current", "previous", "nan"])
def test_invalid_geometry_is_refused(which):
    coords = kernel.unit_cube_Xe()
    displacement = np.zeros(24)
    previous = np.zeros(24)
    if which == "reference":
        coords[:, 0] *= -1
    elif which == "current":
        displacement.reshape(8, 3)[:, 0] = -2 * coords[:, 0]
    elif which == "previous":
        previous.reshape(8, 3)[:, 0] = -2 * coords[:, 0]
    else:
        displacement[0] = np.nan
    with pytest.raises(ValueError, match="Jacobian|real finite"):
        evaluate(coords, displacement, previous=previous)


@pytest.mark.parametrize("constants", [(0, 4.1), (2.3, -1), (2.3,), (np.nan, 4.1)])
def test_invalid_parameters_are_refused(constants):
    with pytest.raises(ValueError, match="requires"):
        evaluate(*distorted_state(), constants=constants)


def test_missing_material_tangent_is_not_silently_zero():
    class NoTangent(CompressibleNeoHookean):
        def evaluate(self, *args, **kwargs):
            stress, _, state, diagnostics = super().evaluate(*args, **kwargs)
            return stress, None, state, diagnostics

    with pytest.raises(ValueError, match="requires material DDSDDE"):
        evaluate(*distorted_state(), material=NoTangent())


def test_unsupported_parameter_and_history_state_are_refused():
    from residual_core.formulations.finite_strain_sensitivity import parameter_slot

    problem, solution = finite_problem()
    with pytest.raises(ValueError, match="unknown finite-strain"):
        parameter_slot(problem.model, "solid.typo")
    with pytest.raises(ValueError, match="history state"):
        SolidC3D8FiniteStrain().eval_element(1, "C3D8", distorted_state()[0], solution,
            {}, np.ones((8, 1)), problem.model.materials["solid"], (0, 1), 1, None, {})