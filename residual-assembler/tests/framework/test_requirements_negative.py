"""Requirements-engine negative tests (audit item 2).

For each intentionally-incomplete model the engine must:
  - not crash,
  - not dump a giant generic checklist,
  - name exactly what is missing and what to provide next (the minimum item).

Cases:
  1. beam model with no section properties            -> section_properties
  2. solid stress-driven model with no stress field   -> element_field
  3. UMAT material-replay model with no UMAT source    -> material_model
  4. shell model with only shell_placeholder available -> recognised, no runnable mode
  5. UEL model with no callable UEL adapter            -> uel_routine

Run: python tests/framework/test_requirements_negative.py
"""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core import ResidualProblem
from residual_core.core.model import Model, Element
from residual_core.core.diagnostics import inspect_model
from residual_core.formulations.registry import build_formulation_registry
from residual_core.materials.registry import build_material_registry
from residual_core.materials.base import MaterialBinding


def _binding(name, section=None):
    b = MaterialBinding(material=None, name=name)
    b.section = section
    return b


def case_beam_no_section():
    m = Model(nodes={1: (0, 0, 0), 2: (1, 0, 0)})
    m.elements[1] = Element(1, "B31", [1, 2])
    m.element_material = {1: "sec"}
    m.materials = {"sec": _binding("sec", section=None)}   # NO section props
    p = ResidualProblem(m)
    r = p.requirements("formulation")
    ok = (not r.runnable) and r.minimum_next.key == "section_properties"
    txt = r.render()
    ok = ok and txt.count("provide") == 1 and "Minimum missing input" in txt
    return "beam / no section", ok, r.minimum_next.display()


def case_solid_no_field():
    m = Model(nodes={i: (0, 0, 0) for i in range(1, 9)})
    m.elements[1] = Element(1, "C3D8", list(range(1, 9)))
    m.element_material = {1: "solid"}
    m.materials = {"solid": _binding("solid", {"E": 210000.0, "nu": 0.3})}
    p = ResidualProblem(m)
    r = p.requirements("stress-driven")
    ok = (not r.runnable) and r.minimum_next.key == "element_field"
    return "solid stress-driven / no field", ok, r.minimum_next.display()


def case_umat_no_source():
    m = Model(nodes={i: (0, 0, 0) for i in range(1, 9)})
    m.elements[1] = Element(1, "C3D8", list(range(1, 9)))
    m.element_material = {1: "CPmat"}
    m.materials = {"CPmat": SimpleNamespace(name="CPmat", user_material=True,
                                            constants=[1.0], depvar=125)}
    p = ResidualProblem(m)
    # supply a solution history so the ONLY thing left missing is the material source
    p.set_solution(np.zeros(1))
    r = p.requirements("material-replay")
    ok = (not r.runnable) and r.minimum_next.key == "material_model"
    return "UMAT replay / no source", ok, r.minimum_next.display()


def case_shell_placeholder_only():
    m = Model(nodes={i: (0, 0, 0) for i in range(1, 5)})
    m.elements[1] = Element(1, "S4", [1, 2, 3, 4])
    forms = build_formulation_registry()
    rep = inspect_model(m, forms, build_material_registry())
    er = next(e for e in rep.elements if e.etype == "S4")
    # recognised by shell_placeholder, but not runnable, with a clear next step
    ok = (er.selected_backend == "shell_placeholder"
          and not er.implemented
          and er.status == "contract-only"
          and er.modes == []
          and bool(er.next_step))
    return "shell / placeholder only", ok, er.next_step


def case_uel_no_callable():
    m = Model(nodes={1: (0, 0, 0), 2: (1, 0, 0)})
    m.elements[1] = Element(1, "U1", [1, 2])
    p = ResidualProblem(m)
    r = p.requirements("direct-residual")
    ok = (not r.runnable) and r.minimum_next.key == "uel_routine"
    return "UEL / no callable", ok, r.minimum_next.display()


_CASES = (case_beam_no_section, case_solid_no_field, case_umat_no_source,
          case_shell_placeholder_only, case_uel_no_callable)


def test_incomplete_models_report_the_minimum_missing_item():
    """Every incomplete model: no crash, and the named minimum next input."""
    failures = []
    for fn in _CASES:
        try:
            label, ok, detail = fn()
        except Exception as exc:               # must never crash
            print("  %-34s -> CRASH: %s" % (fn.__name__, exc))
            failures.append("%s crashed: %s" % (fn.__name__, exc))
            continue
        print("  %-34s -> %s   (next: %s)"
              % (label, "PASS" if ok else "FAIL", detail))
        if not ok:
            failures.append(label)
    assert not failures, "failed cases: " + "; ".join(failures)


def main():
    print("Requirements-engine negative tests")
    try:
        test_incomplete_models_report_the_minimum_missing_item()
        ok = True
    except AssertionError as exc:
        print("  FAIL: %s" % exc)
        ok = False
    except pytest.skip.Exception as exc:
        print("  SKIP: %s" % exc)
        ok = True
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
