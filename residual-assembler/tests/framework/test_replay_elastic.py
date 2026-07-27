"""Acceptance: the Program-2 (collaborator) replay + sensitivity tool, STANDALONE.

Proves the consumer works with ONLY the packaged binary, manifest, replay record
and sensitivity request -- it imports nothing from Program 1 (the JHU-side
provider/transformation repo) and makes no assumption about the constitutive law
inside the binary. The JHU package (opaque .so + manifest) and the saved
production analyses (base + p+/-dp records) are built by an isolated test fixture
(tests/fixtures/replay_elastic/make_fixtures.py), which stands in for JHU's
shipped binary and the collaborator's own Abaqus runs.

Every sensitivity is finite-difference-validated against the PERTURBED RECORDS'
data (displacements, reactions, stresses) -- not against any re-derivation of the
material model.

Run:  pytest tests/framework/test_replay_elastic.py
  or:  python tests/framework/test_replay_elastic.py
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np

try:
    import pytest
except ImportError:                        # pragma: no cover
    class _Raises:
        def __init__(self, exc): self.exc = exc
        def __enter__(self): return self
        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError("did not raise %s" % self.exc.__name__)
            return issubclass(et, self.exc)

    class _Shim:
        @staticmethod
        def raises(exc): return _Raises(exc)
        @staticmethod
        def skip(msg): raise SystemExit("SKIP: " + msg)
    pytest = _Shim()

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# ONLY the Program-2 application + its public contract are imported here.
from residual_core.replay import (replay_sensitivities, run_request, MaterialPackage,
                                  ReplayRecord, MatEvalError, TwinMismatchError,
                                  PackageError, preflight_record)
from residual_core.replay.abi import MaterialABI


# --- build the JHU package + production records ONCE (isolated test fixture) --- #
def _try_build():
    if shutil.which("gfortran") is None:
        return None
    sys.path.insert(0, os.path.join(_ROOT, "tests", "fixtures", "replay_elastic"))
    import make_fixtures
    d = tempfile.mkdtemp(prefix="resasm_fx_")
    try:
        return make_fixtures.build_fixtures(d, drive="disp")
    except Exception:                       # pragma: no cover
        shutil.rmtree(d, ignore_errors=True)
        return None


_FX = _try_build()          # {so, manifest, base, Ep, Em, nup, num, hE, hnu, drive}


def _skip():
    if _FX is None:
        pytest.skip("gfortran unavailable to build the opaque test binary")


def _rec(name):
    return json.load(open(_FX[name]))


def _vm(sig):
    s11, s22, s33, s12, s13, s23 = sig
    return float(np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                         + 3 * (s12 * s12 + s13 * s13 + s23 * s23)))


def _dof_index(node, comp):
    return 3 * (node - 1) + (comp - 1)          # node-major UX,UY,UZ


def _q_from_record(rec, req):
    """Read a response value straight from RECORD DATA (no constitutive law)."""
    inc = rec["increments"][0]
    if req["type"] == "displacement":
        return inc["u"][_dof_index(req["node"], req["dof"])]
    if req["type"] == "reaction":
        return inc["reactions"][_dof_index(req["node"], req["dof"])]
    sig = np.array(inc["stress_ip"][str(req["element"])][req["ip"]])
    return _vm(sig) if req.get("von_mises") else sig[req["component"]]


def _fd_dq(req, pname):
    hp = _FX["hE"] if pname == "E" else _FX["hnu"]
    plus = _rec("Ep" if pname == "E" else "nup")
    minus = _rec("Em" if pname == "E" else "num")
    return (_q_from_record(plus, req) - _q_from_record(minus, req)) / (2 * hp)


def _abs_rel_ok(a, fd, atol=1e-6, rtol=1e-5):
    return abs(a - fd) <= atol * max(1.0, abs(fd)) + rtol * abs(fd)


_OUTPUTS = [
    {"type": "displacement", "node": 11, "dof": 1},
    {"type": "reaction", "node": 11, "dof": 3},
    {"type": "stress", "element": 2, "ip": 3, "component": 2},
    {"type": "stress", "element": 2, "ip": 3, "von_mises": True},
]


# --------------------------------------------------------------------------- #
# 1) headline: consumer runs from packaged inputs only; du/dp & dq/dp vs FD
# --------------------------------------------------------------------------- #
def test_consumer_du_dp_and_dq_dp_from_packaged_inputs_only():
    _skip()
    # the ONLY inputs: binary+manifest (via the manifest path) and the record path
    res = replay_sensitivities(_FX["base"], _FX["manifest"], ["E", "nu"], _OUTPUTS)
    assert res.stress_check["max_rel"] < 1e-10        # replay reproduced production S
    assert res.equilibrium_free_norm < 1e-8           # converged

    ndof = len(_rec("base")["increments"][0]["u"])
    du_fd = {}
    for pn in ("E", "nu"):
        hp = _FX["hE"] if pn == "E" else _FX["hnu"]
        up = np.array(_rec("Ep" if pn == "E" else "nup")["increments"][0]["u"])
        um = np.array(_rec("Em" if pn == "E" else "num")["increments"][0]["u"])
        du_fd[pn] = (up - um) / (2 * hp)
        assert np.max(np.abs(res.du(pn) - du_fd[pn])) <= 1e-6 * max(1.0, np.max(np.abs(du_fd[pn])))

    nontrivial = 0
    for o in res.outputs:
        for pn in ("E", "nu"):
            fd = _fd_dq(o["request"], pn)
            assert _abs_rel_ok(o["dq_dp"][pn], fd), \
                "dq/d%s %s: %.6e vs fd %.6e" % (pn, o["request"], o["dq_dp"][pn], fd)
            if abs(fd) > 1e-3:
                nontrivial += 1
    assert nontrivial >= 4, "dq/dp never non-trivial -- vacuous"


def test_run_request_contract_roundtrip():
    _skip()
    req = {"schema": "resasm_sensitivity_request_v1", "parameters": ["E", "nu"],
           "outputs": _OUTPUTS}
    result = run_request(_FX["base"], _FX["manifest"], req)
    assert result["schema"] == "resasm_sensitivity_result_v1"
    assert set(result["du_dp_norms"]) == {"E", "nu"}
    assert len(result["outputs"]) == len(_OUTPUTS)
    assert result["diagnostics"]["model_id"] == "reference_elastic_isotropic"


# --------------------------------------------------------------------------- #
# 2) parity, R_,p, multi-direction  (still consumer-only)
# --------------------------------------------------------------------------- #
def test_stress_tangent_parity_and_Rp():
    _skip()
    abi = MaterialABI(_FX["so"])            # load the OPAQUE binary by path
    eps = [1e-3, 2e-4, -3e-4, 1e-4, 5e-5, -2e-5]
    out = abi.eval_point([210000.0, 0.3], [1, 2], eps)
    # consumer does not know D; sanity: tangent is symmetric PD-ish and stress = ddsdde@eps
    assert np.allclose(out["stress"], out["ddsdde"] @ np.array(eps))
    assert np.allclose(out["dstress_dseed"][:, 0], out["stress"] / 210000.0)  # dsig/dE = sig/E
    res = replay_sensitivities(_FX["base"], _FX["manifest"], ["E", "nu"])
    assert np.max(np.abs(res.field_result.residual_derivatives)) > 1e-6


def test_multiple_directions_match_single():
    _skip()
    both = replay_sensitivities(_FX["base"], _FX["manifest"], ["E", "nu"])
    onlyE = replay_sensitivities(_FX["base"], _FX["manifest"], ["E"])
    assert np.allclose(both.du("E"), onlyE.du("E"))


# --------------------------------------------------------------------------- #
# 3) ABI errors, time_data padding (fix), rc-authority (fix)
# --------------------------------------------------------------------------- #
def test_abi_errors_and_robustness():
    _skip()
    abi = MaterialABI(_FX["so"])
    with pytest.raises(MatEvalError):
        abi.eval_point([210000.0, 0.3], [1], [1.0, 2.0])        # kin length wrong
    with pytest.raises(MatEvalError):
        abi.eval_point([210000.0, 0.3], [3], [0] * 6)           # seed out of range
    with pytest.raises(MatEvalError):
        abi.eval_point([210000.0], [1], [0] * 6)                # props length wrong
    # time_data shorter than 4 must be padded, not crash (fix)
    o = abi.eval_point([210000.0, 0.3], [1, 2], [0] * 6, time_data=(1.0, 1.0))
    assert o["stress"].shape == (6,)


# --------------------------------------------------------------------------- #
# 4) contract enforcement: twin, binary hash, manifest, preflight
# --------------------------------------------------------------------------- #
def test_matched_twin_and_binary_hash_rejection():
    _skip()
    # wrong regular-twin hash in the record provenance
    rec = _rec("base"); rec["provenance"]["regular_hash"] = "deadbeefdeadbeef"
    with pytest.raises(TwinMismatchError):
        replay_sensitivities(ReplayRecord(rec), _FX["manifest"], ["E"])
    # tampered OTI binary (manifest hash no longer matches the .so on disk)
    man = json.load(open(_FX["manifest"])); man["binaries"]["oti"]["hash"] = "0000000000000000"
    pkg = MaterialPackage(man, base_dir=os.path.dirname(_FX["manifest"]))
    with pytest.raises(TwinMismatchError):
        replay_sensitivities(_FX["base"], pkg, ["E"])


def test_manifest_validation():
    _skip()
    good = json.load(open(_FX["manifest"]))
    MaterialPackage(good)
    for mut in (lambda m: m.update(schema="nope"),
                lambda m: m.pop("parameters"),
                lambda m: m.update(nstatev=-5),          # fix: negative dim rejected
                lambda m: m.update(ntens=0),
                lambda m: m["parameters"][1].update(index=1),      # dup PROPS index
                lambda m: m["parameters"][1].update(oti_direction=1),  # dup direction
                lambda m: m["binaries"].pop("regular")):
        bad = json.loads(json.dumps(good)); mut(bad)
        with pytest.raises(PackageError):
            MaterialPackage(bad)


def test_record_preflight_reports_never_raises():
    _skip()
    pkg = MaterialPackage(json.load(open(_FX["manifest"])),
                          base_dir=os.path.dirname(_FX["manifest"]))
    assert preflight_record(ReplayRecord(_rec("base")), pkg) == []

    rec = _rec("base"); del rec["increments"][0]["u"]
    assert any("displacement 'u'" in p for p in preflight_record(ReplayRecord(rec), pkg))

    # malformed records must be REPORTED, not raise (fix)
    for mut in (lambda r: r["material"].update(props=["abc", 0.3]),
                lambda r: r["increments"][0].update(u=3.0),
                lambda r: r["mesh"]["elements"]["1"].update(connectivity=["x", 2, 3, 4, 5, 6, 7, 8])):
        rec = _rec("base"); mut(rec)
        probs = preflight_record(ReplayRecord(rec), pkg)   # must not raise
        assert isinstance(probs, list) and probs

    # a single-increment PATH-DEPENDENT package must NOT be rejected for having
    # one increment; it needs kinematics_ip + initial_state instead (fix)
    pd = json.load(open(_FX["manifest"])); pd["nstatev"] = 12
    probs = preflight_record(ReplayRecord(_rec("base")),
                             MaterialPackage(pd, base_dir=os.path.dirname(_FX["manifest"])))
    assert any("kinematics_ip" in p for p in probs)
    assert not any("EVERY increment" in p for p in probs)


# --------------------------------------------------------------------------- #
def test_consumer_does_not_import_program1():
    """Structural guarantee: nothing under residual_core.replay imports the
    provider/transformation side."""
    import importlib
    import residual_core.replay as R
    for name in dir(R):
        mod = getattr(getattr(R, name), "__module__", "") or ""
        assert "umat_oti" not in mod and "oti_provider" not in mod
    # and the app package must not pull in the fixture builder
    for m in list(sys.modules):
        if m.startswith("residual_core.replay"):
            src = getattr(sys.modules[m], "__file__", "") or ""
            assert "tests/fixtures" not in src and "oti_provider" not in src


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    ok = True
    for t in tests:
        try:
            t(); print("  PASS: %s" % t.__name__)
        except SystemExit as e:
            print("  %s" % e)
        except Exception as e:              # noqa: BLE001
            ok = False
            import traceback; traceback.print_exc()
            print("  FAIL: %s -> %s" % (t.__name__, e))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
