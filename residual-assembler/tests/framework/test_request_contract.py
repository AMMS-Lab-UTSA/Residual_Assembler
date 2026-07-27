"""Residual request contract (the deck's OUTPUT-PARAMETER REQUESTS): up-front
validation, the guided template, and output-name mapping.

Pure Python -- no OTILib, no Abaqus, no material link. Run:
    python tests/framework/test_request_contract.py   (or: pytest ...)
"""
import json
import os
import sys
import tempfile

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, RA)
from residual_core.replay.request_contract import (           # noqa: E402
    validate_request, new_request_template, _output_kind, _state_names, SCHEMA)


def _materialize(rc, d):
    """Create the object (+ sibling completed contract) and record the request
    references, so existence checks pass. Contract declares E, nu, SIGY0, H and
    a single state variable EQPLAS."""
    obj = os.path.join(d, rc["material"])
    os.makedirs(os.path.dirname(obj) or d, exist_ok=True); open(obj, "w").close()
    contract = {"schema": "resasm_umat_oti_contract_v1",
                "dimensions": {"ntens": 6, "nprops": 4, "nstatev": 1, "nparam": 4},
                "parameters": [{"name": n} for n in ("E", "nu", "SIGY0", "H")],
                "history": {"path_dependent": True, "state": ["EQPLAS"]}}
    json.dump(contract, open(os.path.splitext(obj)[0] + ".json", "w"))
    rp = os.path.join(d, rc["record"]); os.makedirs(os.path.dirname(rp) or d, exist_ok=True); open(rp, "w").close()


def test_template_is_valid():
    with tempfile.TemporaryDirectory() as d:
        rc = new_request_template(); _materialize(rc, d)
        assert validate_request(rc, d) == [], "the built-in template must validate clean"


def test_missing_material_and_record():
    with tempfile.TemporaryDirectory() as d:
        probs = validate_request({"requests": [{"output": "U", "with_respect_to": ["E"]}]}, d)
        assert any("material" in p for p in probs)
        assert any("record" in p for p in probs)


def test_unknown_output_and_parameter():
    with tempfile.TemporaryDirectory() as d:
        rc = new_request_template(); _materialize(rc, d)
        rc["requests"][0]["output"] = "TEMPERATURE"          # not a known field
        rc["requests"][1]["with_respect_to"] = ["E", "BOGUS"]
        probs = validate_request(rc, d)
        assert any("TEMPERATURE" in p for p in probs)
        assert any("BOGUS" in p for p in probs)


def test_bad_increments_and_empty_requests():
    with tempfile.TemporaryDirectory() as d:
        rc = new_request_template(); _materialize(rc, d)
        rc["scope"]["increments"] = "SOMETIMES"
        probs = validate_request(rc, d)
        assert any("increments" in p for p in probs)
        rc2 = new_request_template(); _materialize(rc2, d); rc2["requests"] = []
        assert any("requests" in p for p in validate_request(rc2, d))


def test_output_kind_mapping():
    st = ["EQPLAS"]
    assert _output_kind("U", st) == "displacement"
    assert _output_kind("RF", st) == "reaction"
    assert _output_kind("S", st) == "stress"
    assert _output_kind("EQPLAS", st) == "state"
    assert _output_kind("NOPE", st) == ""


def test_state_names_from_contract():
    assert _state_names({"history": {"state": ["EQPLAS"]}}) == ["EQPLAS"]
    assert _state_names({"dimensions": {"nstatev": 2}}) == ["SDV1", "SDV2"]


def main():
    tests = [f for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    fails = 0
    for t in tests:
        try:
            t(); print("  ok   %s" % t.__name__)
        except AssertionError as e:
            fails += 1; print("  FAIL %s: %s" % (t.__name__, e))
    print("PASS: residual request contract" if not fails else "FAIL: %d" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
