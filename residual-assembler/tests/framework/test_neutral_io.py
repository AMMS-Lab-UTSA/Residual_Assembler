"""Solver-neutral JSON round-trip tests (audit item 6).

Verifies the neutral format preserves the pieces a downstream user relies on:
nodes, elements + types, material/section tags, the ``user_material`` flag,
boundary/load metadata, and optional field-data references. Also checks that a
model loaded from an Abaqus .inp survives a save/load round-trip, and that
``resasm inspect model.json`` works from the CLI.

Run: python tests/framework/test_neutral_io.py
"""

import json
import os
import sys
import tempfile

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.core.model import Model, Element
from residual_core.io import neutral_model_io as nio
from residual_core.materials.base import MaterialBinding


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


def _hand_model():
    m = Model(nodes={1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0)})
    m.elements[1] = Element(1, "B31", [1, 2])
    m.element_formulation = {1: "beam2"}
    m.element_material = {1: "sec"}
    b = MaterialBinding(material=None, name="sec")
    b.section = {"E": 210000.0, "A": 100.0, "Iz": 833.0}
    m.materials = {"sec": b}
    m.boundaries = [{"target": 1, "dof_start": 1, "dof_end": 6, "value": 0.0,
                     "kind": "fixed"}]
    m.cloads = [{"target": 2, "dof": 2, "value": -1000.0}]
    m.meta = {"title": "beam", "field_refs": {"stress_ip": "fields.json"}}
    return m


def test_roundtrip_hand_model():
    m = _hand_model()
    d = tempfile.mkdtemp()
    path = os.path.join(d, "m.json")
    nio.save_model(m, path)
    q = nio.load_model(path)
    checks = [
        ("nodes preserved", q.nodes == m.nodes),
        ("element type preserved", q.elements[1].etype == "B31"),
        ("connectivity preserved", q.elements[1].connectivity == [1, 2]),
        ("formulation tag preserved", q.element_formulation[1] == "beam2"),
        ("material tag preserved", q.element_material[1] == "sec"),
        ("section preserved", q.materials["sec"].section["Iz"] == 833.0),
        ("boundary metadata preserved", q.boundaries[0].target == 1),
        ("cload metadata preserved", q.cloads[0].value == -1000.0),
        ("field ref preserved",
         q.meta.get("field_refs", {}).get("stress_ip") == "fields.json"),
    ]
    for label, v in checks:
        print("  %-38s -> %s" % (label, "PASS" if v else "FAIL"))
    failed = [label for label, v in checks if not v]
    assert not failed, "failed checks: " + "; ".join(failed)


def test_user_material_flag_roundtrip():
    """A UMAT (user material) must survive the round-trip as user_material."""
    from residual_core.io import abaqus_inp_parser
    from residual_core.core.model import from_abaqus
    inp = os.path.join(_ROOT, "sources", "permissive",
                       "ngrilli_Oxford_Crystal_Plasticity", "ExampleInputFiles",
                       "HCPnoTwin", "Compression111.inp")
    if not os.path.exists(inp):
        pytest.skip("user_material test — sample .inp not present")
    am = abaqus_inp_parser.parse_inp(inp)
    m = from_abaqus(am, formulation_policy=lambda et: None)
    m.materials = dict(am.materials)
    d = tempfile.mkdtemp()
    path = os.path.join(d, "cp.json")
    nio.save_model(m, path)
    q = nio.load_model(path)
    name = next(iter(q.materials))
    ok = bool(getattr(q.materials[name], "user_material", False))
    print("  %-38s -> %s" % ("user_material flag round-trips", "PASS" if ok else "FAIL"))
    assert ok, "user_material flag lost in round-trip for material %r" % name


def test_cli_inspect_json():
    from residual_core.ui import cli
    ex = os.path.join(_ROOT, "residual_core", "examples", "minimal_truss", "model.json")
    if not os.path.exists(ex):
        from residual_core.examples import generate_minimal
        generate_minimal.main()
    rc = cli.main(["inspect", ex])
    ok = rc == 0
    print("  %-38s -> %s" % ("resasm inspect model.json (rc=0)", "PASS" if ok else "FAIL"))
    assert ok, "resasm inspect model.json returned rc=%r (need 0)" % rc


def main():
    print("Neutral JSON round-trip tests")
    results = [_script_run(fn) for fn in (test_roundtrip_hand_model,
                                          test_user_material_flag_roundtrip,
                                          test_cli_inspect_json)]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
