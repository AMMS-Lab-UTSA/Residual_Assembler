"""M2 acceptance: the field-driven residual sensitivity, end-to-end from files.

Builds the SAME two-element C3D8 E/nu case as M1, but drives it entirely through
the user-facing path:

    model.inp  +  derivative_fields.json  +  sensitivity.yaml
        -> run_field_recipe()  (what `resasm run` calls)
        -> results/{displacement_sensitivities.npz, residual_derivatives.npz,
                    tangent.npz, manifest.json}

The saved sensitivities must reproduce the direct M1 engine result to < 1e-12.
Also checks: string JSON element ids normalize, unordered IP keys are sorted,
a missing integration point errors, and a misspelled parameter errors.

Run:  pytest tests/framework/test_field_recipe_e2e.py
  or:  python tests/framework/test_field_recipe_e2e.py
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import pytest
except ImportError:                                # script mode without pytest
    class _Raises:
        def __init__(self, exc):
            self.exc = exc
        def __enter__(self):
            return self
        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError("did not raise %s" % self.exc.__name__)
            return issubclass(et, self.exc)

    class _PytestShim:
        @staticmethod
        def raises(exc):
            return _Raises(exc)
    pytest = _PytestShim()

# the M1 fixture (two stacked cubes, full-solve, synthetic fields) is reused
import test_field_sensitivity as m1
from residual_core.core.dof_manager import DofManager
from residual_core.core.field_sensitivity import solve_field_sensitivities
from residual_core.io.derivative_fields import load_derivative_fields
from resasm_user.field_recipe import run_field_recipe
from resasm_user.config import ConfigError

SDV_LAYOUT = {"ddsdde": [1, 36], "parameters": {"E": [37, 42], "nu": [43, 48]}}


# --------------------------------------------------------------------------- #
# serialize the in-memory M1 fixture to on-disk user files
# --------------------------------------------------------------------------- #
def _write_inp(path, model):
    lines = ["*Heading", "** two stacked C3D8 cubes (M2 e2e fixture)", "*Node"]
    for nid in sorted(model.nodes):
        x, y, z = model.nodes[nid]
        lines.append("%d, %.10g, %.10g, %.10g" % (nid, x, y, z))
    lines.append("*Element, type=C3D8")
    for eid in sorted(model.elements):
        conn = model.elements[eid].connectivity
        lines.append("%d, %s" % (eid, ", ".join(str(c) for c in conn)))
    lines.append("*Nset, nset=BASE")
    lines.append("1, 2, 3, 4")
    lines.append("*Boundary")
    lines.append("BASE, ENCASTRE")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def _statev_rows(tangent_fields, dsig, eid):
    """Per-IP SDV vectors: [ DDSDDE(36) | dE(6) | dnu(6) ]."""
    rows = {}
    for ip in range(8):
        rows[str(ip + 1)] = (list(np.asarray(tangent_fields[eid][ip]).ravel())
                             + list(dsig["E"][eid][ip]) + list(dsig["nu"][eid][ip]))
    return rows


def _write_fields_json(path, model, dm, u_star, tangent_fields, dsig,
                       scramble_ips=False, drop_ip=None):
    statev = {}
    for eid in sorted(model.elements):
        rows = _statev_rows(tangent_fields, dsig, eid)
        if scramble_ips:
            rows = {k: rows[k] for k in ["3", "1", "8", "2", "6", "4", "7", "5"]}
        if drop_ip is not None and eid == 1:
            rows = {k: v for k, v in rows.items() if k != str(drop_ip)}
        statev[str(eid)] = rows
    disp = {}
    for nid in sorted(model.nodes):
        gd = dm.node_dofs(nid)
        disp[str(nid)] = [float(u_star[gd[0]]), float(u_star[gd[1]]), float(u_star[gd[2]])]
    doc = {
        "metadata": {
            "schema": "resasm_derivative_fields_v1",
            "element_type": "C3D8",
            "kinematics": "small_strain",
            "voigt_order": ["11", "22", "33", "12", "13", "23"],
            "integration_point_order": "abaqus_label_ascending",
            "sdv_indexing": "abaqus_1_based",
            "sdv_layout": SDV_LAYOUT,
        },
        "displacements": disp,
        "statev": statev,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh)


def _write_yaml(path, params=("E", "nu"), out="results"):
    txt = (
        "analysis:\n"
        "  type: field_residual_sensitivity\n"
        "  kinematics: small_strain\n"
        "mesh:\n"
        "  format: abaqus\n"
        "  file: model.inp\n"
        "results:\n"
        "  format: resasm_derivative_fields_v1\n"
        "  file: derivative_fields.json\n"
        "parameters: [%s]\n"
        "assumptions:\n"
        "  parameter_independent_geometry: true\n"
        "  parameter_independent_loads: true\n"
        "  parameter_independent_boundaries: true\n"
        "output:\n"
        "  directory: %s\n" % (", ".join(params), out)
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(txt)


def _make_case(tmp, scramble_ips=False, drop_ip=None, params=("E", "nu")):
    """Write model.inp + derivative_fields.json + sensitivity.yaml into tmp and
    return (yaml_path, m1_result) where m1_result is the direct engine result."""
    model = m1._two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = m1._load_vector(model, dm)
    u_star = m1._full_solve(model, dm, m1.E0, m1.NU0, fext)
    tangent_fields, dsig = m1._synthetic_fields(model, dm, u_star, m1.E0, m1.NU0)
    res_m1 = solve_field_sensitivities(
        model=model, tangent_fields=tangent_fields,
        stress_derivative_fields=dsig, parameters=["E", "nu"])

    _write_inp(os.path.join(tmp, "model.inp"), model)
    _write_fields_json(os.path.join(tmp, "derivative_fields.json"), model, dm,
                       u_star, tangent_fields, dsig,
                       scramble_ips=scramble_ips, drop_ip=drop_ip)
    yaml_path = os.path.join(tmp, "sensitivity.yaml")
    _write_yaml(yaml_path, params=params)
    return yaml_path, res_m1


# --------------------------------------------------------------------------- #
# 1) headline end-to-end acceptance
# --------------------------------------------------------------------------- #
def test_recipe_matches_m1_engine():
    tmp = tempfile.mkdtemp(prefix="resasm_m2_")
    try:
        yaml_path, res_m1 = _make_case(tmp)
        result = run_field_recipe(yaml_path)
        assert result.ok

        out = os.path.join(tmp, "results")
        for f in ("displacement_sensitivities.npz", "residual_derivatives.npz",
                  "tangent.npz", "manifest.json"):
            assert os.path.exists(os.path.join(out, f)), "missing output %s" % f

        d = np.load(os.path.join(out, "displacement_sensitivities.npz"),
                    allow_pickle=True)
        for j, p in enumerate(["E", "nu"]):
            rec = d["du_d%s" % p]
            ref = res_m1.displacement_sensitivities[:, j]
            rel = np.max(np.abs(rec - ref)) / max(np.max(np.abs(ref)), 1e-30)
            print("  recipe vs M1  du/d%s rel(inf)=%.3e   (tol 1e-12)" % (p, rel))
            assert rel < 1e-12, "du/d%s rel=%.3e" % (p, rel)

        K = np.load(os.path.join(out, "tangent.npz"))["K"]
        assert np.max(np.abs(K - res_m1.K)) / max(np.max(np.abs(res_m1.K)), 1e-30) < 1e-12
        R = np.load(os.path.join(out, "residual_derivatives.npz"),
                    allow_pickle=True)["R"]
        assert np.allclose(R, res_m1.residual_derivatives, atol=0, rtol=1e-12)

        manifest = json.load(open(os.path.join(out, "manifest.json")))
        assert manifest["parameters"] == ["E", "nu"]
        assert manifest["prescribed_dof_sensitivity_is_zero"] is True
        assert manifest["n_prescribed"] == 12
        assert manifest["tangent_source"] == "exported_ddsdde"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# 2) string element ids + unordered IP keys normalize correctly
# --------------------------------------------------------------------------- #
def test_string_ids_and_unordered_ips():
    tmp = tempfile.mkdtemp(prefix="resasm_m2_")
    try:
        # JSON keys are strings already; force a scrambled IP order too.
        yaml_path, res_m1 = _make_case(tmp, scramble_ips=True)

        # the loader normalizes string eids -> int and sorts IP labels
        df = load_derivative_fields(os.path.join(tmp, "derivative_fields.json"))
        assert set(df.statev.keys()) == {1, 2}                 # ints, not "1"/"2"
        assert df.statev[1].shape == (8, 48)

        result = run_field_recipe(yaml_path)
        d = np.load(os.path.join(tmp, "results", "displacement_sensitivities.npz"),
                    allow_pickle=True)
        for j, p in enumerate(["E", "nu"]):
            ref = res_m1.displacement_sensitivities[:, j]
            rel = np.max(np.abs(d["du_d%s" % p] - ref)) / max(np.max(np.abs(ref)), 1e-30)
            assert rel < 1e-12, "scrambled-IP du/d%s rel=%.3e" % (p, rel)
        print("  scrambled IP order reproduces M1 to < 1e-12 (IP sort works)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# 3) honest failures through the recipe
# --------------------------------------------------------------------------- #
def test_missing_integration_point_errors():
    tmp = tempfile.mkdtemp(prefix="resasm_m2_")
    try:
        yaml_path, _ = _make_case(tmp, drop_ip=6)          # element 1 loses IP 6
        with pytest.raises(ConfigError):
            run_field_recipe(yaml_path)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_misspelled_parameter_errors():
    tmp = tempfile.mkdtemp(prefix="resasm_m2_")
    try:
        yaml_path, _ = _make_case(tmp, params=("E", "nuu"))   # 'nuu' not in layout
        with pytest.raises(ConfigError):
            run_field_recipe(yaml_path)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_parameter_dependent_assumption_rejected():
    tmp = tempfile.mkdtemp(prefix="resasm_m2_")
    try:
        yaml_path, _ = _make_case(tmp)
        # flip an assumption to false -> must be rejected
        with open(yaml_path, "r", encoding="utf-8") as fh:
            txt = fh.read()
        txt = txt.replace("parameter_independent_loads: true",
                          "parameter_independent_loads: false")
        with open(yaml_path, "w", encoding="utf-8") as fh:
            fh.write(txt)
        with pytest.raises(ConfigError):
            run_field_recipe(yaml_path)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
def main():
    tests = [test_recipe_matches_m1_engine,
             test_string_ids_and_unordered_ips,
             test_missing_integration_point_errors,
             test_misspelled_parameter_errors,
             test_parameter_dependent_assumption_rejected]
    ok = True
    for t in tests:
        try:
            t()
            print("  PASS: %s" % t.__name__)
        except Exception as exc:            # noqa: BLE001 - script-mode reporting
            ok = False
            print("  FAIL: %s -> %s" % (t.__name__, exc))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
