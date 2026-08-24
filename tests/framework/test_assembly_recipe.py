"""Path A — residual assembly from ingredients (the PRIMARY user story).

    The user should not provide R.
    The user should provide enough ingredients for us to build R.

These tests pin the recipe layer:

  * progressive inference — element type, backend, DOF map, BCs, integration rule
    and defaults are all derived from the mesh; the user types almost nothing;
  * the minimum-missing-input report;
  * the HONESTY GATE — a backend that can ASSEMBLE R but cannot carry an OTI
    number must be reported as `OTI-differentiate: NO` and must REFUSE to produce
    a sensitivity, rather than emit a plausible-looking wrong number.

The full assemble -> OTI -> solve run needs OTILib and is covered by
`test_assembly_recipe_runs` (skips cleanly without it; see
scripts/run_otilib_tests_wsl.sh).
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.algebra.otilib_adapter import otilib_available  # noqa: E402
from resasm_user.recipe import load_recipe, describe, sensitivity_capability  # noqa: E402
from resasm_user.config import ConfigError  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from external_sources import require_external_file  # noqa: E402

C3D8_INP = os.path.join(
    _ROOT, "sources", "permissive", "ngrilli_Oxford_Crystal_Plasticity",
    "ExampleInputFiles", "HCPnoTwin", "Compression111.inp")
SPRING_JSON = os.path.join(
    _ROOT, "residual_core", "examples", "minimal_nonlinear_spring_sensitivity",
    "model.json")


def _write(tmp_path, text, **files):
    d = str(tmp_path)
    with open(os.path.join(d, "resasm.yml"), "w") as fh:
        fh.write(text)
    for name, arr in files.items():
        np.save(os.path.join(d, name), np.asarray(arr, float))
    return os.path.join(d, "resasm.yml")


# --------------------------------------------------------------------------- #
# Progressive inference: the user types 3 lines; we infer the rest.
# --------------------------------------------------------------------------- #
def test_infers_everything_inferable_from_an_abaqus_mesh(tmp_path):
    require_external_file(
        C3D8_INP, "a real C3D8 Abaqus mesh to exercise recipe inference and the OTI honesty gate")
    cfg = _write(tmp_path, "problem:\n  name: c3d8_case\nmesh: %s\nparameters:\n"
                           "  - CPuranium.E\n" % C3D8_INP.replace("\\", "/"))
    r = load_recipe(cfg)

    # read straight off the mesh -- the user never typed any of this
    assert r.mesh.format == "abaqus_inp"
    assert r.mesh.element_types == ["C3D8"]
    assert r.mesh.n_nodes == 216 and r.mesh.n_elements == 125
    assert r.formulation.backends == {"C3D8": "solid_c3d8_finite_strain"}
    for key in ("mesh.format", "mesh.element_types", "formulation.backend",
                "formulation.integration", "dof_map", "constraints",
                "sensitivity.order", "sensitivity.backend", "output.dir"):
        assert key in r.inferred, key

    txt = describe(r)
    assert "C3D8 -> solid_c3d8_finite_strain" in txt
    assert "216 / 125" in txt


def test_reports_the_minimum_missing_ingredients(tmp_path):
    require_external_file(
        C3D8_INP, "a real C3D8 Abaqus mesh to exercise recipe inference and the OTI honesty gate")
    cfg = _write(tmp_path, "mesh: %s\nparameters:\n  - CPuranium.E\n"
                 % C3D8_INP.replace("\\", "/"))
    r = load_recipe(cfg)
    blob = " ".join(r.missing)
    assert "solution" in blob          # converged U
    assert "material" in blob          # material evaluator
    assert not r.blockers              # C3D8 itself IS supported


# --------------------------------------------------------------------------- #
# THE HONESTY GATE
# --------------------------------------------------------------------------- #
def test_c3d8_can_assemble_but_is_NOT_oti_differentiable(tmp_path):
    """Assembling R and differentiating R are different capabilities.

    The C3D8 kernels use numpy float arrays (core/voigt.py::isotropic_D raises on
    an OTI number), so parameters cannot be overloaded through them. The tool must
    say so, not pretend."""
    require_external_file(
        C3D8_INP, "a real C3D8 Abaqus mesh to exercise recipe inference and the OTI honesty gate")
    cfg = _write(tmp_path, "mesh: %s\nparameters:\n  - CPuranium.E\n"
                 % C3D8_INP.replace("\\", "/"))
    r = load_recipe(cfg)
    cap = sensitivity_capability(r)

    assert cap["can_assemble"] is True
    assert cap["oti_differentiable"] is False
    assert "C3D8 -> solid_c3d8_finite_strain" in cap["oti_blocked"]
    assert "numpy" in cap["reason"].lower()
    assert "OTI-differentiate R  : NO" in describe(r)


def test_run_REFUSES_a_sensitivity_it_cannot_compute(tmp_path):
    """Rather than emit a wrong number, `resasm run` must refuse and explain."""
    require_external_file(
        C3D8_INP, "a real C3D8 Abaqus mesh to exercise recipe inference and the OTI honesty gate")
    from resasm_user import run_from_config
    cfg = _write(tmp_path,
                 "mesh: %s\nsolution: U.npy\nmaterial:\n  type: elastic\n"
                 "  properties:\n    E: 210000.0\n    nu: 0.3\n"
                 "parameters:\n  - CPuranium.E\n" % C3D8_INP.replace("\\", "/"),
                 **{"U.npy": np.zeros(648)})
    with pytest.raises(ConfigError) as ei:
        run_from_config(cfg)
    msg = str(ei.value)
    assert "OTI-differentiate R : NO" in msg
    assert "assemble R          : yes" in msg      # it is honest about what DOES work
    assert "black-box" in msg                      # and names the routes that do work


def test_spring_model_IS_oti_differentiable(tmp_path):
    cfg = _write(tmp_path, "mesh: %s\nsolution: U.npy\nparameters:\n  - spring.k\n"
                 % SPRING_JSON.replace("\\", "/"), **{"U.npy": np.array([2.0])})
    r = load_recipe(cfg)
    cap = sensitivity_capability(r)
    assert r.formulation.backends == {"SPRING1": "nonlinear_spring1"}
    assert cap["can_assemble"] and cap["oti_differentiable"]
    assert not r.missing and not r.blockers
    assert "Ready to assemble." in describe(r)


# --------------------------------------------------------------------------- #
# The real thing: ingredients -> assemble R -> OTI -> solve  (needs OTILib)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not otilib_available(),
                    reason="OTILib not installed; see scripts/run_otilib_tests_wsl.sh")
def test_assembly_recipe_runs_end_to_end(tmp_path):
    """No residual function anywhere: only a mesh, a solution and a parameter."""
    from resasm_user import run_from_config
    cfg = _write(tmp_path, "problem:\n  name: spring_case\nmesh: %s\n"
                           "solution: U.npy\nparameters:\n  - spring.k\n"
                 % SPRING_JSON.replace("\\", "/"), **{"U.npy": np.array([2.0])})
    res = run_from_config(cfg)

    assert res.summary["mode"] == "assemble"
    assert res.summary["orders_solved"] == [1]

    z = np.load(os.path.join(res.private_dir,
                             "solution_sensitivities_order1.npz"))
    # exact: u = (f/k)^(1/3) with k=2, f=16 -> u=2, du/dk = -u/(3k) = -1/3
    assert abs(float(z["U_derivatives"].ravel()[0]) + 1.0 / 3.0) < 1e-8
    # the assembly path produces the same package as every other path
    for f in ("tangent.npz", "rhs_order1.npz", "direction_map_order1.json"):
        assert os.path.exists(os.path.join(res.private_dir, f)), f
    assert os.path.exists(os.path.join(res.public_dir, "summary.md"))
