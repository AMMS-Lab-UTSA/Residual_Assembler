"""Public-API test (audit item 5): the user touches ONLY ResidualProblem.

This module imports nothing from residual_core.formulations.* or
residual_core.materials.* — no backend classes, no DofManager, no Assembler. If a
new user can drive the framework with just the public facade, this passes.

    from residual_core import ResidualProblem
    problem = ResidualProblem.from_neutral("model.json")
    report  = problem.inspect()
    problem.attach_results("fields.json")
    R = problem.assemble(mode="stress-driven")

Run: python tests/framework/test_public_api.py
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# The ONLY framework import a user needs:
from residual_core import ResidualProblem

_EX = os.path.join(_ROOT, "residual_core", "examples")


def _ensure_examples():
    """Generate the minimal example assets if they are not present."""
    c3d8 = os.path.join(_EX, "minimal_c3d8_stress_driven", "model.json")
    if not os.path.exists(c3d8):
        from residual_core.examples import generate_minimal
        generate_minimal.main()


@pytest.fixture(scope="module", autouse=True)
def _examples():
    """Same asset bootstrap the script entry point does, for the pytest run."""
    _ensure_examples()


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


def _report(label, ok, detail=""):
    print("  %-38s -> %s" % (label, "PASS" if ok else "FAIL"))
    assert ok, "%s failed%s" % (label, (": " + detail) if detail else "")


def test_stress_driven_via_facade():
    model = os.path.join(_EX, "minimal_c3d8_stress_driven", "model.json")
    fields = os.path.join(_EX, "minimal_c3d8_stress_driven", "fields.json")
    problem = ResidualProblem.from_neutral(model)
    report = problem.inspect(print_it=False)
    ok = any(e.etype == "C3D8" for e in report.elements)
    problem.attach_results(fields)
    R = problem.assemble(mode="stress-driven")
    ok = ok and R.size == 24 and np.linalg.norm(R) > 1.0
    _report("C3D8 stress-driven via facade", ok,
            "R.size=%d (need 24), ||R||=%.3e (need > 1.0)" % (R.size, np.linalg.norm(R)))


def test_formulation_via_facade():
    model = os.path.join(_EX, "minimal_truss", "model.json")
    problem = ResidualProblem.from_neutral(model)
    R = problem.assemble(mode="formulation")
    ok = R.size == 9 and np.linalg.norm(R) < 1e-9   # zero at U=0
    _report("truss formulation via facade", ok,
            "R.size=%d (need 9), ||R||=%.3e (need < 1e-9)" % (R.size, np.linalg.norm(R)))


def test_requirements_via_facade():
    model = os.path.join(_EX, "minimal_c3d8_stress_driven", "model.json")
    problem = ResidualProblem.from_neutral(model)     # no field attached yet
    r = problem.requirements("stress-driven")
    ok = (not r.runnable) and r.minimum_next.key == "element_field"
    _report("requirements report via facade", ok,
            "runnable=%r minimum_next=%r (need not-runnable + 'element_field')"
            % (r.runnable, r.minimum_next.key))


def _uses_only_public_api():
    """Static guard: this test module's IMPORT statements reference no backend or
    internal engine classes — only the public ``residual_core.ResidualProblem``."""
    banned_mods = ("residual_core.formulations", "residual_core.materials",
                   "residual_core.core.assembler", "residual_core.core.dof_manager")
    hits = []
    with open(__file__, "r", encoding="utf-8") as fh:
        source = fh.read()
    for raw in source.splitlines():
        line = raw.strip()
        if not (line.startswith("from ") or line.startswith("import ")):
            continue
        for mod in banned_mods:
            if line.startswith("from " + mod) or line.startswith("import " + mod):
                hits.append(line)
    _report("imports only public API (no backends)", not hits,
            "banned imports: %s" % hits)


def main():
    print("Public-API test (facade only)")
    _ensure_examples()
    cases = [test_stress_driven_via_facade, test_formulation_via_facade,
             test_requirements_via_facade, _uses_only_public_api]
    results = [_script_run(fn) for fn in cases]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
