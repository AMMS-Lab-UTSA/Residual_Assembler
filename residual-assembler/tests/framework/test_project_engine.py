"""M5-A: the project/action engine -- resumable, idempotent, logged.

Offline by default (no Abaqus). Exercises the manifest + state machine, action
idempotency (skip-if-complete), failure recovery (a failed action keeps the last
good state), precondition enforcement, inspection, parameter configuration, and
the offline assemble -> validate -> report chain driven by a committed fixture as
the "exported" derivative fields. A gated test runs the full Abaqus pipeline.
"""

import json
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
_ASSET = os.path.join(_ROOT, "tests", "abaqus_derivative_export")
_FIX = os.path.join(_ASSET, "fixtures")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from resasm_user.project import (Project, State, run_action, run_pipeline,
                                 set_parameter_selection)

_INP = os.path.join(_ASSET, "fcc_cubic_c3d8.inp")
_UMAT = os.path.join(_ASSET, "cubic_oti_umat.for")
_SEL = [{"name": "C11", "index": 1}, {"name": "C12", "index": 2},
        {"name": "C44", "index": 3}]


def _new(tmp):
    return Project.create(os.path.join(tmp, "proj"), name="t", inp=_INP, umat=_UMAT)


# --------------------------------------------------------------------------- #
def test_manifest_and_inputs_copied():
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        assert os.path.exists(p.path("project.json"))
        assert os.path.exists(p.path("inputs", "fcc_cubic_c3d8.inp"))
        assert os.path.exists(p.path("inputs", "cubic_oti_umat.for"))
        assert p.state == State.NEW
        # originals untouched: the project holds copies
        assert p.input_path("inp") != _INP
        p2 = Project.load(p.root)              # round-trips
        assert p2.data["name"] == "t"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_precondition_and_failure_recovery():
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        # configure requires INSPECTED -> fails, state stays NEW (resume point kept)
        r = run_action(p, "configure_parameters")
        assert r.status == "failed"
        assert p.state == State.NEW
        run_pipeline(p, State.INSPECTED)
        assert p.state == State.INSPECTED
        # still no selection -> configure fails cleanly, keeps INSPECTED
        r = run_action(p, "configure_parameters")
        assert r.status == "failed" and p.state == State.INSPECTED
        assert os.path.exists(p.logfile("configure_parameters"))   # structured log
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_inspection_reports():
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        run_pipeline(p, State.INSPECTED)
        m = p.data["model"]
        assert m["element_types"] == {"C3D8": 1}
        assert m["integration"] == "selective_reduced"
        assert m["supported"] is True
        assert m["materials"]["FCCCUBIC"]["n_constants"] == 3
        u = p.data["umat"]
        assert u["already_oti"] is True and u["writes_ddsdde"] is True
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_configure_and_idempotency():
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        run_pipeline(p, State.INSPECTED)
        set_parameter_selection(p, _SEL)
        assert run_action(p, "configure_parameters").status == "complete"
        assert p.state == State.CONFIGURED
        params = p.data["parameters"]
        assert [x["oti_direction"] for x in params] == ["E1", "E2", "E3"]
        assert [x["sdv_range"] for x in params] == [[37, 42], [43, 48], [49, 54]]
        assert params[0]["value"] == 168400.0
        # idempotent: re-run is skipped
        assert run_action(p, "configure_parameters").status == "skipped"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_prepare_generates_artifacts():
    oti = os.path.exists(os.path.expanduser(
        os.environ.get("OTI_DIR", "~/MultiZ_f/oti")) + "/libotim4n1.a")
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        run_pipeline(p, State.INSPECTED)
        set_parameter_selection(p, _SEL)
        run_action(p, "configure_parameters")
        r = run_action(p, "prepare")
        if not oti:                             # env generation needs the OTI lib
            assert r.status == "failed"
            return
        assert r.status == "complete" and p.state == State.PREPARED
        layout = json.load(open(p.artifact("generated", "derivative_layout.json")))
        assert layout["sdv_layout"]["parameters"]["C44"] == [49, 54]
        assert layout["integration"] == "selective_reduced"
        assert os.path.exists(p.artifact("generated", "umat.for"))
        assert os.path.exists(p.artifact("generated", "abaqus_v6.env"))
        assert "*Depvar" in open(p.artifact("generated", "model.inp")).read() \
            or "*DEPVAR" in open(p.artifact("generated", "model.inp")).read().upper()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_offline_assemble_validate_report():
    """Stage a committed fixture as the exported fields and drive the downstream
    (offline) actions -- no Abaqus needed."""
    tmp = tempfile.mkdtemp(prefix="resasm_proj_")
    try:
        p = _new(tmp)
        run_pipeline(p, State.INSPECTED)
        set_parameter_selection(p, _SEL)
        run_action(p, "configure_parameters")
        # stage generated inputs the downstream actions need (bypassing Abaqus)
        shutil.copy(_INP, p.artifact("generated", "model.inp"))
        shutil.copy(os.path.join(_FIX, "fcc_fields.json"),
                    p.artifact("generated", "derivative_fields.json"))
        p.advance_to(State.EXTRACTED)
        p.save()
        assert run_action(p, "assemble_sensitivities").status == "complete"
        assert run_action(p, "run_validation").status == "complete"
        assert run_action(p, "generate_report").status == "complete"
        assert p.state == State.REPORTED
        man = json.load(open(p.artifact("results", "manifest.json")))
        assert set(man["du_da_norms"]) == {"C11", "C12", "C44"}
        val = json.load(open(p.artifact("results", "validation_report.json")))
        assert val["passed"] is True
        assert os.path.exists(p.artifact("results", "report.md"))
        assert os.path.exists(p.artifact("results", "project_manifest.json"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_full_pipeline_with_abaqus():
    """Gated: create -> inspect -> configure -> prepare -> run Abaqus -> export ->
    assemble -> validate -> report, all through the action engine."""
    import shutil as _sh
    oti = os.path.exists(os.path.expanduser(
        os.environ.get("OTI_DIR", "~/MultiZ_f/oti")) + "/libotim4n1.a")
    if not (_sh.which(os.environ.get("ABAQUS_CMD", "abaqus"))
            and os.environ.get("RESASM_RUN_ABAQUS") and oti):
        return _skip("needs RESASM_RUN_ABAQUS=1, Abaqus on PATH and the OTI lib")
    tmp = tempfile.mkdtemp(prefix="resasm_proj_abq_")   # paren-free for the UMAT compile
    try:
        p = _new(tmp)
        run_pipeline(p, State.INSPECTED)
        set_parameter_selection(p, _SEL)
        run_action(p, "configure_parameters")
        results = run_pipeline(p, State.REPORTED)
        assert all(r.ok for r in results), [(r.action, r.status, r.error) for r in results]
        assert p.state == State.REPORTED
        man = json.load(open(p.artifact("results", "manifest.json")))
        assert all(v == v for v in man["du_da_norms"].values())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


_IN_SCRIPT = False


def _skip(msg):
    if _IN_SCRIPT:
        print("  SKIP: %s" % msg)
        return
    import pytest
    pytest.skip(msg)


def main():
    global _IN_SCRIPT
    _IN_SCRIPT = True
    tests = [test_manifest_and_inputs_copied, test_precondition_and_failure_recovery,
             test_inspection_reports, test_configure_and_idempotency,
             test_prepare_generates_artifacts, test_offline_assemble_validate_report,
             test_full_pipeline_with_abaqus]
    ok = True
    for t in tests:
        try:
            t()
            print("  PASS: %s" % t.__name__)
        except Exception as exc:            # noqa: BLE001
            ok = False
            print("  FAIL: %s -> %s" % (t.__name__, exc))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
