"""The inspection commands, run as a user runs them, on models this repository ships.

``resasm backends``, ``template``, ``inspect-model``, ``requirements``,
``init-assembly`` and ``doctor`` only report; these tests hold each report to
what it describes. The backend listing must be the registry's own declaration,
and the claims it makes that can be executed (no modes, no tangent, an analytic
tangent) are executed. The model reports are checked against the parsed model
and, where they say something is ready, by running it.

Every command is called through ``residual_core.ui.cli.main`` in-process.
"""
import os
import re
from pathlib import Path

import numpy as np
import pytest

from residual_core.algebra.otilib_adapter import otilib_available
from residual_core.core import requirements as requirements_engine
from residual_core.formulations import c3d8_kernel as kernel
from residual_core.formulations.registry import build_formulation_registry
from residual_core.io.abaqus_inp_parser import parse_inp
from residual_core.materials.base import MaterialBinding
from residual_core.materials.registry import build_material_registry
from residual_core.ui import cli

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "residual_core" / "examples"
CUBE = EXAMPLES / "minimal_c3d8_stress_driven"
SPRING = EXAMPLES / "minimal_nonlinear_spring_sensitivity" / "model.json"
DECK = ROOT / "examples" / "presentation_request" / "Analysis.inp"


def run(capsys, *argv):
    code = cli.main([str(a) for a in argv])
    captured = capsys.readouterr()
    return code, captured.out, captured.err


# --------------------------------------------------------------------------- #
# resasm backends
# --------------------------------------------------------------------------- #
def _listing_blocks(text):
    """{backend name: [its indented lines]} from the ``resasm backends`` text."""
    body = text.split("===============================\n", 1)[1]
    blocks, materials = {}, None
    for chunk in body.strip().split("\n\n"):
        lines = chunk.splitlines()
        if lines[0].startswith("Materials: "):
            materials = lines[0][len("Materials: "):].split(", ")
            continue
        blocks[lines[0]] = lines[1:]
    return blocks, materials


def _declared_lines(spec):
    lines = ["  element types : %s" % ", ".join(spec.supported_element_types),
             "  dof types     : %s" % ", ".join(spec.dof_types),
             "  status        : %s" % spec.verification_status,
             "  available modes: %s" % (", ".join(spec.supported_modes) or "none"),
             "  material iface : %s" % spec.material_interface_needed,
             "  state         : %s" % spec.state_requirements,
             "  tangent       : %s" % spec.tangent_support]
    if spec.limitations:
        lines.append("  limitations   : %s" % "; ".join(spec.limitations))
    return lines


def test_backends_lists_every_registered_backend_as_it_is_declared(capsys):
    code, out, _ = run(capsys, "backends")
    assert code == 0
    assert out.startswith("Registered formulation backends\n")
    blocks, materials = _listing_blocks(out)
    registry = build_formulation_registry()
    assert list(blocks) == registry.names()
    for name in registry.names():
        assert blocks[name] == _declared_lines(registry.get(name).spec), name
    assert materials == build_material_registry().names()


def _status_of(blocks, name, field):
    prefix = "  %s" % field
    return next(line.split(":", 1)[1].strip() for line in blocks[name] if line.startswith(prefix))


def _c3d8_inputs():
    coords = kernel.unit_cube_Xe().copy()
    coords[6] += [0.12, -0.05, 0.08]
    displacement = 1e-3 * np.cos(np.arange(24.0))
    return coords, displacement


def _element_case(name):
    """A valid single-element evaluation for each backend whose listing claims
    a tangent: (coords, dofs, properties, fields, options)."""
    from residual_core.materials.elastic_adapter import IsotropicElastic
    from residual_core.materials.neo_hookean import CompressibleNeoHookean
    two_nodes = np.array([[0.0, 0.0, 0.0], [2.0, 1.0, 0.5]])
    if name == "truss2":
        return two_nodes, 1e-3 * np.array([0.3, -0.2, 0.1, 0.7, 0.4, -0.5]), \
            {"E": 210000.0, "A": 3.0}
    if name == "beam2":
        return two_nodes, 1e-3 * np.sin(np.arange(12.0)), \
            {"E": 210000.0, "A": 3.0, "Iz": 2.5, "Iy": 1.5}
    if name == "nonlinear_spring1":
        return np.zeros((1, 3)), np.array([1.3]), {"k": 2.0, "f": 16.0}
    if name == "nonlinear_bar1":
        return two_nodes, np.array([0.2, 0.9]), {"k": 3.0}
    if name == "solid_c3d8_small_strain":
        coords, displacement = _c3d8_inputs()
        return coords, displacement, MaterialBinding(IsotropicElastic(), [210000.0, 0.3])
    if name == "solid_c3d8_finite_strain":
        coords, displacement = _c3d8_inputs()
        return coords, 100 * displacement, MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1])
    raise KeyError(name)


def _evaluate(backend, coords, dofs, properties, compute_tangent=True, fields=None):
    etype = backend.element_types[0]
    r, K, _, _ = backend.eval_element(1, etype, coords, np.asarray(dofs, float), {}, None,
                                      properties, (0.0, 1.0), 1.0, fields,
                                      {"compute_tangent": compute_tangent})
    return np.asarray(r, float), (None if K is None else np.asarray(K, float))


def test_what_the_listing_claims_is_what_each_backend_does(capsys):
    """Run the claims: a backend listed with no modes cannot evaluate an
    element, the skeleton has no routine, a backend listed 'tangent: none'
    returns none, and every tangent listed as analytic matches central
    differences of that backend's own residual."""
    code, out, _ = run(capsys, "backends")
    assert code == 0
    blocks, _ = _listing_blocks(out)
    registry = build_formulation_registry()
    coords, displacement = _c3d8_inputs()

    no_modes = [n for n in blocks if _status_of(blocks, n, "available modes") == "none"]
    assert no_modes == ["shell_placeholder"]
    assert _status_of(blocks, "shell_placeholder", "status") == "contract-only"
    with pytest.raises(NotImplementedError):
        _evaluate(registry.get("shell_placeholder"), np.zeros((4, 3)), np.zeros(24), None)

    assert _status_of(blocks, "uel_direct", "status") == "skeleton"
    with pytest.raises(RuntimeError, match="no uel_fn"):
        _evaluate(registry.get("uel_direct"), np.zeros((8, 3)), np.zeros(24), None)

    no_tangent = [n for n in blocks if _status_of(blocks, n, "tangent") == "none"]
    assert set(no_tangent) == {"stress_driven_c3d8", "shell_placeholder"}
    stress = {"stress_ip": {1: np.tile([10.0, -4.0, 2.0, 1.5, 0.5, -0.7], (8, 1))}}
    force, tangent = _evaluate(registry.get("stress_driven_c3d8"), coords, displacement,
                               None, fields=stress)
    assert tangent is None and np.all(np.isfinite(force)) and np.abs(force).max() > 1

    analytic = [n for n in blocks if _status_of(blocks, n, "tangent").startswith("analytic")]
    assert analytic, "no backend claims an analytic tangent"
    for name in analytic:
        backend = registry.get(name)
        element_coords, dofs, properties = _element_case(name)
        _, tangent = _evaluate(backend, element_coords, dofs, properties)
        step = 1e-6 * max(1.0, float(np.abs(dofs).max()))
        columns = []
        for j in range(dofs.size):
            plus, minus = dofs.copy(), dofs.copy()
            plus[j] += step
            minus[j] -= step
            columns.append((_evaluate(backend, element_coords, plus, properties, False)[0]
                            - _evaluate(backend, element_coords, minus, properties, False)[0])
                           / (2 * step))
        differenced = np.array(columns).T
        scale = max(float(np.abs(tangent).max()), 1e-30)
        assert np.abs(tangent - differenced).max() / scale < 1e-6, name


# --------------------------------------------------------------------------- #
# resasm template
# --------------------------------------------------------------------------- #
def test_template_prints_the_declared_formulation_contract(capsys):
    registry = build_formulation_registry()
    for name in registry.names():
        code, out, _ = run(capsys, "template", "--formulation", name)
        assert code == 0
        declared = registry.get(name).spec.as_dict()
        lines = out.splitlines()
        assert lines[0] == "backend: %s" % name
        for key in ("kind", "verification_status", "supported_element_types",
                    "supported_modes", "tangent_support", "limitations"):
            if declared[key] or declared[key] is False:
                assert "  %-26s %s" % (key + ":", declared[key]) in lines, (name, key)
        for mode, items in declared["required_inputs_by_mode"].items():
            assert "      %-18s %s" % (mode + ":", items) in lines, (name, mode)


def test_template_prints_a_material_contract_with_its_measures(capsys):
    code, out, _ = run(capsys, "template", "--material", "compressible_neo_hookean")
    assert code == 0
    assert "  kind:                      material" in out
    assert "  supported_modes:           ['solid_c3d8_finite_strain']" in out
    assert "  required_inputs:           ['F0', 'F1']" in out
    assert "  tangent_support:           kirchhoff_jaumann_over_j" in out
    assert "input=deformation_gradient; stress=cauchy; tangent=kirchhoff_jaumann_over_j" in out
    code, out, _ = run(capsys, "template", "--material", "umat")
    assert code == 0
    assert "  required_inputs:           ['F0', 'F1']" in out
    assert "  tangent_support:           ddsdde" in out


def test_template_of_an_unknown_backend_exits_2_listing_the_known_ones(capsys):
    code, out, err = run(capsys, "template", "--formulation", "no_such_backend")
    assert code == 2 and out == ""
    assert err.startswith("unknown backend 'no_such_backend'; known formulations: ")
    for name in build_formulation_registry().names() + build_material_registry().names():
        assert name in err


# --------------------------------------------------------------------------- #
# resasm inspect-model
# --------------------------------------------------------------------------- #
def test_inspect_model_of_the_deck_names_what_is_missing_and_exits_1(capsys):
    code, out, _ = run(capsys, "inspect-model", DECK)
    assert code == 1
    parsed = parse_inp(str(DECK))
    assert "  nodes/elements  : %d / %d" % (len(parsed.nodes), len(parsed.elements)) in out
    assert "  element types   : C3D8" in out
    assert re.search(r"constraints +%d \*Boundary block\(s\) read from the mesh"
                     % len(parsed.boundaries), out)
    assert re.search(r"stimuli\.loads +%d \*Cload\(s\) read from the mesh"
                     % len(parsed.cloads), out)
    assert "  OTI-differentiate R  : NO\n     blocked: C3D8 -> solid_c3d8_finite_strain" in out
    missing = [line.split("] ", 1)[1].split(":", 1)[0]
               for line in out.splitlines() if line.startswith("  [missing] ")]
    assert missing == ["fields.solution", "material", "parameters"]
    assert "Ready to assemble." not in out


def test_inspect_model_with_every_ingredient_is_ready_and_exits_0(tmp_path, capsys):
    solution = tmp_path / "U.npy"
    np.save(solution, np.array([2.0]))
    code, out, _ = run(capsys, "inspect-model", SPRING, "--solution", solution,
                       "--param", "spring.k")
    assert code == 0
    assert "  solution field  : %s" % solution in out
    assert "  parameters      : spring.k" in out
    assert re.search(r"formulation\.backend +auto-selected per element type: "
                     r"SPRING1 -> nonlinear_spring1", out)
    assert "  assemble R           : yes\n  OTI-differentiate R  : yes" in out
    assert "[missing]" not in out
    assert out.rstrip().endswith("Ready to assemble.")


# --------------------------------------------------------------------------- #
# resasm requirements
# --------------------------------------------------------------------------- #
def test_requirements_names_the_one_missing_input_and_assemble_agrees(capsys):
    code, out, _ = run(capsys, "requirements", CUBE / "model.json", "--mode", "stress-driven")
    assert code == 0
    assert out.splitlines()[:5] == ["Cannot assemble in stress-driven mode.", "Available:",
                                    "  mesh: yes", "  solution field (U / U+rotation / T): yes",
                                    "  stress / resultant field: no"]
    assert out.count("provide ") == 1
    assert "Minimum missing input:\n  provide integration-point stress field S" in out
    # the report is what assemble acts on: it refuses with the same text
    code, blocked, err = run(capsys, "assemble", CUBE / "model.json", "--mode", "stress-driven")
    assert code == 2 and blocked == "" and err == out


def test_requirements_with_the_field_is_ready_and_assemble_runs(capsys):
    fields = CUBE / "fields.json"
    code, out, _ = run(capsys, "requirements", CUBE / "model.json", "--mode", "stress-driven",
                       "--fields", fields)
    assert code == 0
    assert out.startswith("Ready to assemble in stress-driven mode.\n")
    assert "  stress / resultant field: yes" in out
    assert out.rstrip().endswith("All minimum inputs are present.")
    code, assembled, _ = run(capsys, "assemble", CUBE / "model.json", "--mode", "stress-driven",
                             "--fields", fields)
    assert code == 0
    assert assembled.startswith("assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01")


def test_requirements_of_the_deck_in_replay_mode_asks_for_the_history(capsys):
    code, out, _ = run(capsys, "requirements", DECK, "--mode", "material-replay")
    assert code == 0
    assert out.startswith("Cannot assemble in material-replay mode.\n")
    assert "  solution history: no" in out
    assert "  material parameters (PROPS): yes" in out
    assert "Minimum missing input:\n  provide the solution history" in out


# --------------------------------------------------------------------------- #
# resasm init-assembly
# --------------------------------------------------------------------------- #
def test_init_assembly_writes_a_recipe_that_checks_and_runs(tmp_path, capsys):
    from resasm_user.recipe import load_recipe
    np.save(tmp_path / "U.npy", np.array([2.0]))
    recipe = tmp_path / "resasm.yml"
    code, out, _ = run(capsys, "init-assembly", "--model", SPRING, "--solution", "U.npy",
                       "--param", "spring.k", "--name", "spring_case", "--out", recipe)
    assert code == 0
    assert out.startswith("wrote %s\n" % recipe)
    assert recipe.read_text().splitlines() == [
        "problem:", "  name: spring_case", "", "mesh: %s" % SPRING, "",
        "solution: U.npy", "", "parameters:", "  - spring.k", "",
        "sensitivity:", "  order: 1", "  backend: otilib"]
    assert "Ready to assemble." in out
    assert "Next:\n  resasm check %s\n  resasm run   %s" % (recipe, recipe) in out

    loaded = load_recipe(str(recipe))
    assert loaded.name == "spring_case"
    assert loaded.parameters == ["spring.k"]
    assert loaded.formulation.backends == {"SPRING1": "nonlinear_spring1"}
    assert not loaded.missing and not loaded.blockers

    code, checked, _ = run(capsys, "check", recipe)
    assert code == 0, checked
    assert "Ready to assemble." in checked


@pytest.mark.skipif(not otilib_available(),
                    reason="OTILib not installed; the recipe's run needs it")
def test_the_recipe_init_assembly_writes_solves_the_spring(tmp_path, capsys):
    np.save(tmp_path / "U.npy", np.array([2.0]))
    recipe = tmp_path / "resasm.yml"
    assert run(capsys, "init-assembly", "--model", SPRING, "--solution", "U.npy",
               "--param", "spring.k", "--out", recipe)[0] == 0
    code, out, err = run(capsys, "run", recipe)
    assert code == 0, err
    private = tmp_path / "resasm_output" / "private"
    derivative = np.load(private / "solution_sensitivities_order1.npz")["U_derivatives"]
    # R = k u^3 - f with k = 2, f = 16: u = 2 and du/dk = -u / (3 k) = -1/3
    assert abs(float(derivative.ravel()[0]) + 1.0 / 3.0) < 1e-12


def test_init_assembly_without_a_model_exits_2(tmp_path, capsys):
    code, out, err = run(capsys, "init-assembly", "--out", tmp_path / "resasm.yml")
    assert code == 2 and out == ""
    assert "ERROR: --model <mesh file> is required" in err
    assert not (tmp_path / "resasm.yml").exists()


# --------------------------------------------------------------------------- #
# resasm doctor
# --------------------------------------------------------------------------- #
def _readiness(text):
    lines = text.split("per-mode readiness:\n", 1)[1].split("\n\n", 1)[0].splitlines()
    return {line.split()[0]: line.split(None, 1)[1] for line in lines}


def test_doctor_prints_one_readiness_line_per_mode_from_the_requirements(tmp_path, capsys):
    from residual_core import ResidualProblem
    template = tmp_path / "doctor.yml"
    code, out, _ = run(capsys, "doctor", CUBE / "model.json",
                       "--write-config-template", template)
    assert code == 0
    assert out.startswith("Model inspection summary\n")
    assert "  - C3D8 elements (x1): supported by" in out
    readiness = _readiness(out)
    assert list(readiness) == requirements_engine.known_modes()
    problem = ResidualProblem.from_neutral(str(CUBE / "model.json"))
    for mode, line in readiness.items():
        report = problem.requirements(mode)
        assert line == ("ready" if report.runnable
                        else "needs: %s" % report.minimum_next.display()), mode
    assert readiness["stress-driven"] == "needs: stress / resultant field"
    assert readiness["direct-residual"] == "needs: callable UEL adapter"
    assert out.rstrip().endswith("config template written to %s" % template)


def test_doctor_writes_a_config_template_that_reads_back(tmp_path, capsys):
    from residual_core.ui.config import Config
    written = tmp_path / "doctor.yml"
    assert run(capsys, "doctor", CUBE / "model.json", "--write-config-template", written)[0] == 0
    text = written.read_text()
    for key in ("mode", "odb", "subroutine", "formulation_policy", "material_backend",
                "material_parameters", "options"):
        assert re.search(r"^%s:" % key, text, re.MULTILINE), key
        assert key == "options" or re.search(r"^#   %s +: " % key, text, re.MULTILINE), key
    as_json = tmp_path / "doctor.json"
    assert run(capsys, "doctor", CUBE / "model.json", "--write-config-template", as_json)[0] == 0
    config = Config.load(str(as_json))
    assert config.mode is None and config.odb is None and config.subroutine is None
    assert run(capsys, "--config", as_json, "doctor", CUBE / "model.json")[0] == 0
    yaml = pytest.importorskip("yaml")
    assert yaml.safe_load(text)["formulation_policy"] == {}
    assert run(capsys, "--config", written, "doctor", CUBE / "model.json")[0] == 0
