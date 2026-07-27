"""Up-front sensitivity-job validation, the guided template, and the reduction
default fix (residual_core/replay/job.py).

Pure Python -- no OTILib, no Abaqus, no material link. Run:
    python tests/framework/test_job_validation.py     (or: pytest tests/framework/test_job_validation.py)
"""
import json
import os
import sys
import tempfile

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, RA)
from residual_core.replay.job import (              # noqa: E402
    validate_job, new_job_template, _param_names, _DEFAULT_REDUCTION, SCHEMA)


def _materialize(job, d):
    """Create the files a valid job references inside dir d (so existence
    checks pass), with a contract that declares parameters E and nu."""
    for rel in (job["material"]["oti_umat"], job["analysis"]["record"]):
        p = os.path.join(d, rel); os.makedirs(os.path.dirname(p), exist_ok=True); open(p, "w").close()
    cp = os.path.join(d, job["material"]["contract"]); os.makedirs(os.path.dirname(cp), exist_ok=True)
    json.dump({"parameters": [{"name": "E"}, {"name": "nu"}]}, open(cp, "w"))


def test_template_is_valid():
    with tempfile.TemporaryDirectory() as d:
        job = new_job_template(); _materialize(job, d)
        assert validate_job(job, d) == [], "the built-in template must validate clean"


def test_missing_blocks_reported():
    with tempfile.TemporaryDirectory() as d:
        probs = validate_job({}, d)
        assert any("schema" in p for p in probs)
        assert any("material" in p for p in probs)
        assert any("analysis" in p for p in probs)
        assert any("parameters" in p for p in probs)
        assert any("outputs" in p for p in probs)


def test_bad_component_and_reduction_and_type():
    with tempfile.TemporaryDirectory() as d:
        job = new_job_template(); _materialize(job, d)
        job["outputs"][0]["component"] = "U9"           # invalid for nodal_displacement
        job["outputs"][1]["reduction"] = "mean"          # not in the enum
        job["outputs"].append({"output_id": "x", "type": "bogus_type",
                               "region": {"node_set": "N"}, "component": "U1"})
        probs = validate_job(job, d)
        assert any("component" in p and "U9" in p for p in probs)
        assert any("reduction" in p and "mean" in p for p in probs)
        assert any("type" in p and "bogus_type" in p for p in probs)


def test_unknown_parameter_reported_early():
    with tempfile.TemporaryDirectory() as d:
        job = new_job_template(); _materialize(job, d)
        job["parameters"] = [{"name": "E"}, {"name": "BOGUS"}]
        probs = validate_job(job, d)
        assert any("BOGUS" in p for p in probs), "unknown parameter must be caught before the link"


def test_duplicate_output_id_and_missing_region():
    with tempfile.TemporaryDirectory() as d:
        job = new_job_template(); _materialize(job, d)
        job["outputs"][1]["output_id"] = job["outputs"][0]["output_id"]   # duplicate
        job["outputs"][0].pop("region")                                    # missing region
        probs = validate_job(job, d)
        assert any("duplicate output_id" in p for p in probs)
        assert any("region" in p for p in probs)


def test_missing_referenced_files():
    with tempfile.TemporaryDirectory() as d:
        job = new_job_template()                          # files NOT created
        probs = validate_job(job, d)
        assert any("not found" in p for p in probs)


def test_parameters_shorthand_accepted():
    assert _param_names({"parameters": ["E", "nu"]}) == ["E", "nu"]
    assert _param_names({"parameters": [{"name": "E"}, {"name": "nu"}]}) == ["E", "nu"]


def test_reduction_default_fix():
    # the bug: a reaction over several nodes with no explicit reduction used to
    # fall through to the first node's value; the default is now type-specific.
    assert _DEFAULT_REDUCTION["reaction_force"] == "sum"
    assert _DEFAULT_REDUCTION["integration_point_stress"] == "volume_average"
    assert _DEFAULT_REDUCTION["nodal_displacement"] == "average"
    assert (None or _DEFAULT_REDUCTION["reaction_force"]) == "sum"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    fails = 0
    for t in tests:
        try:
            t(); print("  ok   %s" % t.__name__)
        except AssertionError as e:
            fails += 1; print("  FAIL %s: %s" % (t.__name__, e))
    print("PASS: job validation" if not fails else "FAIL: %d test(s)" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
