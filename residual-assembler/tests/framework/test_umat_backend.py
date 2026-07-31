"""M5-B: the umat_oti transformation-backend adapter (slices 1 and 2).

Proves Residual Assembler can drive the separate UMAT_source_transformation
(umat_oti) project through the adapter to:

  * slice 1 -- turn an ordinary (non-OTI) UMAT into an OTI UMAT that produces the
    consistent tangent DDSDDE, and emit an HONEST resasm_umat_transform_v1 handoff
    manifest (no parameter derivatives claimed when none were generated); and
  * slice 2 (B2.2) -- generate the parameter-seeding contract so the transformed
    UMAT ALSO writes d(sigma)/da_i to SDVs and the real DDSDDE to STATEV[1..36],
    flip the manifest's parameter_derivatives to true, and VALIDATE it numerically
    at a material point against a central finite difference of the ORIGINAL UMAT.

The repositories stay separate (nothing GPL is copied in). The fixture UMAT here
is our own ordinary 3D elastic UMAT that reads E,nu from PROPS. The numerical
real-response check needs gfortran and is skipped (not failed) when it is absent.
"""

import json
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_FIX = os.path.join(_ROOT, "tests", "fixtures")
_UMAT = os.path.join(_FIX, "ordinary_elastic_3d.for")
_BASE = os.path.join(_FIX, "ordinary_elastic_3d_base.json")

_IN_SCRIPT = False


def _transformer_available():
    src = os.environ.get("UMAT_OTI_SRC") or os.path.expanduser(
        "~/Documents/UMAT_source_transformation/src")
    return (os.path.exists(_UMAT) and os.path.exists(_BASE)
            and os.path.isdir(os.path.join(src, "umat_oti")))


def _gfortran():
    return shutil.which("gfortran") is not None


def _skip(msg):
    if _IN_SCRIPT:
        print("  SKIP: %s" % msg)
        return True
    import pytest
    pytest.skip(msg)


def _minimal_c3d8_inp(path):
    """A trivial 1-element C3D8 with a 2-constant user material (E, nu)."""
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    lines = ["*Node"]
    for i, c in enumerate(coords, 1):
        lines.append("%d, %g, %g, %g" % (i, c[0], c[1], c[2]))
    lines += ["*Element, type=C3D8, elset=ALL", "1, 1,2,3,4,5,6,7,8",
              "*Solid Section, elset=ALL, material=ELAS",
              "*Material, name=ELAS", "*Depvar", "40,",
              "*User Material, constants=2", "210000., 0.3",
              "*Nset, nset=BASE", "1,2,3,4", "*Step, name=Step-1",
              "*Static", "*Boundary", "BASE, ENCASTRE", "*End Step"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
def test_backend_ddsdde_only_honest_manifest():
    """No parameters requested -> DDSDDE tangent only; the manifest must NOT claim
    any parameter derivatives, and validate_real_response has nothing to compare."""
    if not _transformer_available():
        return _skip("umat_oti transformer not available (set UMAT_OTI_SRC)")
    from resasm_user.umat_backend import UmatOtiBackend
    tmp = tempfile.mkdtemp(prefix="resasm_xform_")
    try:
        be = UmatOtiBackend()
        insp = be.inspect(_UMAT)
        assert insp.already_oti is False and insp.writes_ddsdde is True
        assert insp.props_refs == [1, 2]
        prop = be.propose_contract(_UMAT, [], base_config=_BASE)   # no parameters
        res = be.transform(prop, tmp)
        assert res.success is True, res.blockers
        assert all(res.semantic_checks.values())
        assert be.validate_transformation_semantics(res).passed is True
        man = res.manifest
        assert man["schema"] == "resasm_umat_transform_v1"
        assert man["capabilities"]["automatic_ddsdde"] is True
        assert man["capabilities"]["parameter_derivatives"] is False
        assert man["derivatives"]["parameters"] == {}
        assert man["abaqus"]["required_depvar"] == 0
        try:
            be.validate_real_response(res)
            assert False, "real-response must be NotImplemented with no derivatives"
        except NotImplementedError:
            pass
        print("  ddsdde-only: honest manifest (parameter_derivatives=False, depvar=0)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_backend_parameter_derivatives_generated_and_validated():
    """E + nu requested -> transformed UMAT writes d(sigma)/dE, d(sigma)/dnu and the
    real DDSDDE to SDVs; the manifest flips to parameter_derivatives=true and the
    numbers match a central finite difference of the ORIGINAL UMAT."""
    if not _transformer_available():
        return _skip("umat_oti transformer not available")
    from resasm_user.umat_backend import UmatOtiBackend, ParameterSelection
    tmp = tempfile.mkdtemp(prefix="resasm_xform_")
    try:
        be = UmatOtiBackend()
        prop = be.propose_contract(
            _UMAT, [ParameterSelection("E", 1), ParameterSelection("nu", 2)],
            base_config=_BASE)
        res = be.transform(prop, tmp)
        assert res.success is True, res.blockers
        assert all(res.semantic_checks.values())
        man = res.manifest
        assert man["capabilities"]["parameter_derivatives"] is True
        params = man["derivatives"]["parameters"]
        assert params["E"]["sdv_range"] == [37, 42]
        assert params["E"]["status"] == "generated"
        assert params["nu"]["sdv_range"] == [43, 48]
        assert man["derivatives"]["ddsdde"]["in_sdv"] is True
        assert man["derivatives"]["ddsdde"]["sdv_range"] == [1, 36]
        assert man["requested_parameter_derivatives"]["E"]["status"] == "generated"
        assert man["abaqus"]["required_depvar"] == 48         # 36 DDSDDE + 2*6 params
        if not _gfortran():
            return _skip("gfortran absent; skipping numerical real-response check")
        assert man["numerical_validation"]["status"] == "not_run"   # before the check
        rr = be.validate_real_response(res)
        assert rr.checks.get("available") is True, rr.message
        assert rr.passed is True, rr.message
        assert rr.checks["max_relerr"] < 1e-6
        assert rr.checks["ddsdde"]["max_relerr"] < 1e-6
        # the verdict is stamped into the handoff manifest (verified != generated)
        assert man["numerical_validation"]["status"] == "passed"
        assert man["numerical_validation"]["max_relerr"] < 1e-6
        print("  param derivs: dsig/dE, dsig/dnu, DDSDDE -> SDV; "
              "numerical relerr %.2e (%s)" % (rr.checks["max_relerr"], rr.message))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_backend_rejects_incoherent_parameter_requests():
    """Honesty guards: a PROPS index the UMAT never reads, or a base config missing
    the tangent variable, must be rejected up front rather than producing a silently
    zero derivative or a manifest that advertises an unwritten DDSDDE SDV block."""
    if not _transformer_available():
        return _skip("umat_oti transformer not available")
    from resasm_user.umat_backend import (UmatOtiBackend, ParameterSelection,
                                          UmatBackendError)
    be = UmatOtiBackend()
    # (a) index not referenced by the UMAT (fixture reads only PROPS(1),PROPS(2))
    try:
        be.propose_contract(_UMAT, [ParameterSelection("ghost", 3)], base_config=_BASE)
        assert False, "unread PROPS index must be rejected"
    except UmatBackendError:
        pass
    # (b) parameters requested but no tangent variable -> cannot pack DDSDDE to SDV
    base = json.load(open(_BASE))
    base["resasm"].pop("tangent_variable", None)
    tmp = tempfile.mkdtemp(prefix="resasm_cfg_")
    try:
        bt = os.path.join(tmp, "no_tangent.json")
        json.dump(base, open(bt, "w"))
        try:
            be.propose_contract(_UMAT, [ParameterSelection("E", 1)], base_config=bt)
            assert False, "missing tangent_variable must be rejected when params given"
        except UmatBackendError:
            pass
        print("  guards: unread-index and missing-tangent both rejected up front")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_project_transform_umat_action_generates_parameters():
    """The action engine drives the ordinary UMAT -> parameter-derivative OTI UMAT
    and records the honest manifest + numerical real-response verdict."""
    if not _transformer_available():
        return _skip("umat_oti transformer not available")
    from resasm_user.project import (Project, State, run_pipeline, run_action,
                                     set_parameter_selection, set_umat_oti_config)
    tmp = tempfile.mkdtemp(prefix="resasm_xproj_")
    try:
        inp = os.path.join(tmp, "model.inp")
        _minimal_c3d8_inp(inp)
        p = Project.create(os.path.join(tmp, "proj"), name="elastic_auto",
                           inp=inp, umat=_UMAT)
        run_pipeline(p, State.INSPECTED)
        assert p.data["umat"]["already_oti"] is False        # ordinary UMAT
        set_parameter_selection(p, [{"name": "E", "index": 1},
                                    {"name": "nu", "index": 2}])
        run_action(p, "configure_parameters")
        set_umat_oti_config(p, _BASE)
        r = run_action(p, "transform_umat")
        assert r.status == "complete", r.error
        assert os.path.exists(p.artifact("generated", "umat_transformed.for"))
        man = json.load(open(p.artifact("generated", "umat_transform_manifest.json")))
        assert man["transformation"]["success"] is True
        assert man["capabilities"]["parameter_derivatives"] is True
        assert man["abaqus"]["required_depvar"] == 48
        if _gfortran():
            rr = p.data.get("transform_real_response")
            assert rr and rr.get("available") is True and rr.get("passed") is True, rr
            print("  action: ordinary UMAT -> parameter derivatives, numerically "
                  "validated (max relerr %.2e)" % rr["max_relerr"])
        else:
            print("  action: ordinary UMAT -> parameter derivatives (gfortran absent)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    global _IN_SCRIPT
    _IN_SCRIPT = True
    ok = True
    for t in (test_backend_ddsdde_only_honest_manifest,
              test_backend_parameter_derivatives_generated_and_validated,
              test_backend_rejects_incoherent_parameter_requests,
              test_project_transform_umat_action_generates_parameters):
        try:
            t()
            print("  PASS: %s" % t.__name__)
        except Exception as exc:            # noqa: BLE001
            ok = False
            import traceback
            traceback.print_exc()
            print("  FAIL: %s -> %s" % (t.__name__, exc))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
