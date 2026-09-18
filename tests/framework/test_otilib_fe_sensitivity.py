"""OTILib FE sensitivity — multi-element nonlinear bar chain (order 1).

Two coupled cubic bar elements (nodes 1-2-3), node 1 fixed, node 3 prescribed,
per-element stiffness parameters k1, k2. Proves the OTILib path at the
FINITE-ELEMENT level: element residuals evaluated with OTI scalars, coefficients
extracted and assembled into a global R^(1), sensitivity solved and matched to
finite differences from re-solving the real nonlinear FE problem.

The real-FE assembly/solve check runs ALWAYS (proves the coupled backend). The
OTI sensitivity check requires OTILib and SKIPS cleanly when absent.

Run: python tests/framework/test_otilib_fe_sensitivity.py
"""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.algebra import otilib_adapter as A
from residual_core import ResidualProblem
from residual_core.core.model import Model, Element
from residual_core.materials.base import MaterialBinding

U3 = 2.0                     # prescribed tip displacement


def _require_otilib():
    """Skip when OTILib is absent — unless RUN_OTILIB_TESTS=1 demands it."""
    if A.otilib_available():
        return
    assert os.environ.get("RUN_OTILIB_TESTS") != "1", (
        "RUN_OTILIB_TESTS=1 but OTILib is not installed; install: "
        "scripts/setup_otilib.sh (https://github.com/mauriaristi/otilib.git)")
    pytest.skip("OTILib not installed: %s"
                % A.otilib_status()["error"].splitlines()[0])


def _script_run(fn, *args):
    """Script-mode runner: a test passes unless it raises AssertionError."""
    try:
        fn(*args)
    except AssertionError as exc:
        print("  FAIL: %s" % exc)
        return False
    except pytest.skip.Exception as exc:
        print("  SKIP: %s" % exc)
    return True


def _report(checks):
    """Print every check, then assert that none failed (same condition as all())."""
    for label, v in checks:
        print("  %-42s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def _bar_chain(k1=1.0, k2=1.0):
    m = Model(nodes={1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)})
    m.elements[1] = Element(1, "NLBAR2", [1, 2])
    m.elements[2] = Element(2, "NLBAR2", [2, 3])
    m.element_material = {1: "barA", 2: "barB"}
    a = MaterialBinding(material=None, name="barA"); a.section = {"k": k1}
    b = MaterialBinding(material=None, name="barB"); b.section = {"k": k2}
    m.materials = {"barA": a, "barB": b}
    m.boundaries = [
        SimpleNamespace(target=1, dof_start=1, dof_end=1, value=0.0, kind="fixed"),
        SimpleNamespace(target=3, dof_start=1, dof_end=1, value=U3, kind="fixed"),
    ]
    return ResidualProblem(m)


def _real_fe_checks():
    p = _bar_chain()
    u = p.solve_newton("formulation", u0=[0.0, 0.5, U3])
    res = p.result("formulation", U=u)
    checks = [
        ("bar chain assembles (ndof=3)", res.ndof == 3),
        ("Newton converged R~0", res.free_residual_norm < 1e-9),
        ("u2 = 1 (symmetric k)", abs(u[1] - 1.0) < 1e-9),
        ("tangent backend-assembled", res.tangent.source == "backend-assembled"),
    ]
    return checks


def _oti_checks():
    p = _bar_chain()
    p.solve_newton("formulation", u0=[0.0, 0.5, U3])
    params = ["barA.k", "barB.k"]
    pkg = p.sensitivity_package(mode="formulation", parameters=params,
                                max_order=1, generate_rhs=True, backend="otilib")
    s = pkg.sensitivity
    U1 = pkg.solve(1)                    # (3, 2)
    # finite-difference reference by re-solving the real nonlinear FE problem
    fd = p.finite_difference_sensitivity("formulation", params, u0=[0.0, 0.5, U3])
    free = pkg.residual.free_mask
    rel = float(np.linalg.norm((U1 - fd)[free]) / max(np.linalg.norm(fd[free]), 1e-30))
    checks = [
        ("R^(1) shape (3,2)", np.asarray(s.R(1)).shape == (3, 2)),
        ("element-assembled OTI residual", s.diagnostics.get("hypercomplex_ready") is True),
        ("OTI sensitivity matches FD (rel<1e-4)", rel < 1e-4),
        ("prescribed rows zero", np.allclose(U1[~free], 0.0)),
    ]
    return checks


def test_real_fe_chain():
    """Coupled nonlinear bar chain assembles, solves and tangents — always runs."""
    _report(_real_fe_checks())


def test_oti_fe_sensitivity():
    """OTI-assembled R^(1) vs finite differences — requires OTILib."""
    _require_otilib()
    _report(_oti_checks())


@pytest.mark.parametrize("mode", ["small", "finite"])
def test_c3d8_stress_preserves_mixed_third_derivatives(mode):
    from residual_core.formulations.c3d8_kernel import (
        element_internal_force_finite_strain,
        element_internal_force_small_strain, unit_cube_Xe)

    _require_otilib()
    context = A.OtiContext(2, 3)
    amplitude = context.seed(2.0, 1) ** 2 * context.seed(3.0, 2)
    coordinates = unit_cube_Xe()
    stretches = np.array([1.2, 0.9, 1.1]) if mode == "finite" else np.ones(3)
    displacement = coordinates * (stretches - 1.0)
    stress = np.array([[amplitude * entry for entry in
                        (2.0, 3.0, 4.0, 0.5, 0.7, 0.9)]] * 8)
    if mode == "finite":
        force = element_internal_force_finite_strain(
            coordinates, displacement, stress)
    else:
        force = element_internal_force_small_strain(coordinates, stress)
    tensor = np.array([[2.0, 0.5, 0.7], [0.5, 3.0, 0.9], [0.7, 0.9, 4.0]])
    face_areas = np.prod(stretches) / stretches
    traction = (((2.0 * coordinates - 1.0) * face_areas) @ tensor.T / 4).ravel()
    for direction, factor in [((0, 0), 12.0), ((1, 0), 12.0),
                              ((0, 1), 4.0), ((2, 0), 6.0),
                              ((1, 1), 4.0), ((2, 1), 2.0)]:
        actual = np.array([context.deriv(value, direction) for value in force])
        np.testing.assert_allclose(actual, factor * traction, rtol=1e-13,
                                   atol=1e-13, err_msg=str(direction))


@pytest.mark.parametrize("mode", ["small", "finite"])
@pytest.mark.parametrize("assembly", ["force", "parameter_derivative"])
def test_c3d8_global_scatter_preserves_mixed_scalar_contributions(mode, assembly):
    from residual_core.formulations.c3d8_kernel import (
        assemble_global_internal_force, unit_cube_Xe)
    from residual_core.formulations.c3d8_sensitivity import assemble_dR_dp

    _require_otilib()
    context = A.OtiContext(2, 3)
    amplitude = context.seed(2.0, 1) ** 2 * context.seed(3.0, 2)
    coordinates = unit_cube_Xe()
    node_ids = [80, 20, 50, 10, 70, 40, 60, 30]
    connectivity = [(1, node_ids), (2, node_ids)]
    real_stress = np.tile([1.0, 0.0, 0.0, 0.0, 0.0, 0.0], (8, 1))
    oti_stress = np.array([[amplitude, 0.0, 0.0, 0.0, 0.0, 0.0]] * 8)
    stresses = {1: real_stress, 2: oti_stress}
    if assembly == "force":
        result, mapping = assemble_global_internal_force(
            node_ids, coordinates, connectivity, np.zeros(24), stresses, mode)
    else:
        result, mapping = assemble_dR_dp(
            node_ids, coordinates, connectivity, stresses, np.zeros(24), mode=mode)
    expected = np.zeros((8, 3))
    for local, node_id in enumerate(node_ids):
        expected[mapping[node_id], 0] = (2.0 * coordinates[local, 0] - 1.0) / 4
    assert mapping == {node_id: index for index, node_id in enumerate(sorted(node_ids))}
    for direction, factor in [((0, 0), 13.0), ((1, 0), 12.0), ((2, 1), 2.0)]:
        actual = np.array([context.deriv(value, direction) for value in result])
        np.testing.assert_allclose(actual, factor * expected.ravel(),
                                   rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("mode", ["small", "finite"])
def test_c3d8_material_stiffness_preserves_mixed_third_derivatives(mode):
    from residual_core.formulations.c3d8_kernel import element_tangent, unit_cube_Xe

    _require_otilib()
    context = A.OtiContext(2, 3)
    amplitude = context.seed(2.0, 1) ** 2 * context.seed(3.0, 2)
    coordinates = unit_cube_Xe()
    tangent = np.diag([amplitude * entry for entry in (2, 3, 4, 5, 6, 7)])
    stiffness = element_tangent(coordinates, np.zeros(24), tangent, mode=mode)
    strain_tensor = np.array([[0.01, 0.02, 0.03],
                              [0.02, 0.04, 0.05],
                              [0.03, 0.05, 0.06]])
    displacement = (coordinates @ strain_tensor).ravel()
    stress_tensor = np.array([[0.02, 0.2, 0.36],
                              [0.2, 0.12, 0.7],
                              [0.36, 0.7, 0.24]])
    traction = ((2 * coordinates - 1) @ stress_tensor.T / 4).ravel()
    force = stiffness @ displacement
    for direction, factor in [((0, 0), 12.0), ((1, 0), 12.0), ((2, 1), 2.0)]:
        actual = np.array([context.deriv(value, direction) for value in force])
        np.testing.assert_allclose(actual, factor * traction, rtol=1e-13, atol=1e-13)


def test_c3d8_geometric_stiffness_preserves_stress_derivatives():
    from residual_core.formulations.c3d8_kernel import element_tangent, unit_cube_Xe

    _require_otilib()
    context = A.OtiContext(2, 3)
    amplitude = context.seed(2.0, 1) ** 2 * context.seed(3.0, 2)
    coordinates = unit_cube_Xe()
    stress = np.array([[amplitude, 0.0, 0.0, 0.0, 0.0, 0.0]] * 8)
    stiffness = element_tangent(coordinates, np.zeros(24), np.zeros((6, 6)),
                                sigma_ip=stress, mode="finite")
    force = stiffness @ coordinates.ravel()
    traction = np.zeros((8, 3))
    traction[:, 0] = (2 * coordinates[:, 0] - 1) / 4
    for direction, factor in [((0, 0), 12.0), ((1, 0), 12.0), ((2, 1), 2.0)]:
        actual = np.array([context.deriv(value, direction) for value in force])
        np.testing.assert_allclose(actual, factor * traction.ravel(),
                                   rtol=1e-13, atol=1e-13)


def main():
    print("OTILib FE sensitivity (nonlinear bar chain)")
    results = [_script_run(test_real_fe_chain), _script_run(test_oti_fe_sensitivity)]
    ok = all(results)
    tail = "ALL PASS" if ok else "FAILURE"
    if not A.otilib_available() and ok:
        tail += " (real-FE only; OTI part skipped)"
    print("OVERALL: %s" % tail)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
