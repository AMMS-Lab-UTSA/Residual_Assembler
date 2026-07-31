"""Interface-version drift guard.

residual_core/interface/versions.py is the single source of truth for every
schema tag. This test asserts that each module still AGREES with it, so a change
to any one literal (or to versions.py) is caught here instead of silently
diverging. Pure Python -- no OTILib, no Abaqus.

    python tests/framework/test_interface_versions.py     (or: pytest ...)
"""
import os
import sys

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, RA)
from residual_core.interface import versions as v          # noqa: E402


def test_declared_tags_match_versions():
    from residual_core.replay import job, record, objlink, package
    from residual_core.io import derivative_fields
    assert job.SCHEMA == v.SENSITIVITY_JOB
    assert record.SCHEMA == v.REPLAY_RECORD
    assert objlink.CONTRACT_SCHEMA == v.UMAT_OTI_CONTRACT
    assert package.SCHEMA == v.MATERIAL_PACKAGE
    assert derivative_fields.SCHEMA_ID == v.DERIVATIVE_FIELDS


def test_every_tag_is_versioned_and_unique():
    tags = list(v.SCHEMAS.values())
    assert all(t.startswith("resasm_") and "_v" in t for t in tags)
    assert len(set(tags)) == len(tags), "duplicate schema tags"


def test_contract_hash_present():
    h = v.contract_version()
    assert isinstance(h, str) and len(h) > 0, "combined ABI/package contract hash missing"


def test_check_raises_on_drift():
    v.check("resasm_sensitivity_job_v1", "resasm_sensitivity_job_v1")   # no raise
    try:
        v.check("resasm_sensitivity_job_v1", "something_else")
    except ValueError:
        return
    raise AssertionError("check() must raise on a mismatched tag")


def main():
    tests = [f for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    fails = 0
    for t in tests:
        try:
            t(); print("  ok   %s" % t.__name__)
        except AssertionError as e:
            fails += 1; print("  FAIL %s: %s" % (t.__name__, e))
    print("PASS: interface versions" if not fails else "FAIL: %d" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
