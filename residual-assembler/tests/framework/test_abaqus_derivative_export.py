"""M3 acceptance: real Abaqus ODB -> derivative fields -> residual sensitivities.

The heavy Abaqus run happened once; its outputs are committed as fixtures under
``tests/abaqus_derivative_export/fixtures/`` (from a real Abaqus 2021 analysis of
the nonuniform single-C3D8 job with the verification elastic UMAT). These tests
validate the whole pipeline OFFLINE against that real data, so the suite needs no
Abaqus.

What is proven here:

* Exporter/sidecar (M3.1/M3.2): the pure builder + ``resasm_sdv_layout_v1`` loader
  produce a valid ``resasm_derivative_fields_v1`` doc and reject malformed input.
* Tangent (M3.5): K assembled from the ODB-exported DDSDDE matches the analytic
  C3D8 tangent to < 1e-6 (the single-precision ODB stores SDV as float32, so
  ~1e-7 is the floor -- the spec's 1e-10 is unreachable via an ODB).
* Integration-point ordering (M3.3): per-IP stress recovered from the exported
  fields matches Abaqus S at the SAME IP label; an intentional IP permutation is
  detected. Independently, assembling the exported stress reproduces the applied
  load (equilibrium), which a permuted order breaks.
* Formulation (B-bar): Abaqus C3D8 uses selective reduced (mean-dilatation)
  volumetric integration; the B-bar recovery matches Abaqus S to ~1e-7 while full
  integration is off by O(1).
* End-to-end (M3.6): ``run_field_recipe`` on the exported file reproduces the
  Abaqus displacement sensitivities. du/dnu matches full central-difference Abaqus
  reruns to < 1e-5; du/dE to ~2e-5 (that column is intrinsically ~1e-8 in size, so
  it sits at the float32-ODB finite-difference floor -- confirmed, not a defect).

Run:  pytest tests/framework/test_abaqus_derivative_export.py
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
_ASSET = os.path.join(_ROOT, "tests", "abaqus_derivative_export")
_FIX = os.path.join(_ASSET, "fixtures")
for _p in (_ROOT,):
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

from residual_core.io import abaqus_inp_parser as P
from residual_core.core.model import from_abaqus
from residual_core.core.dof_manager import DofManager
from residual_core.core.constraints import partition
from residual_core.formulations import c3d8_kernel as kern
from residual_core.core.field_sensitivity import fields_from_statev
from residual_core.io.derivative_fields import (
    load_derivative_fields, DerivativeFieldError)
from residual_core.io.export_derivative_fields import (
    load_sdv_layout, build_derivative_fields, ExportError)
from resasm_user.field_recipe import run_field_recipe
from resasm_user.config import ConfigError

E0, NU0 = 210000.0, 0.3
H_E, H_NU = 840.0, 0.004           # the finite-difference steps used for the reruns
INP = os.path.join(_ASSET, "nonuniform_c3d8.inp")


def _uvec(fname, dm):
    j = json.load(open(os.path.join(_FIX, fname)))
    u = np.zeros(dm.ndof)
    for n, v in j.items():
        gd = dm.node_dofs(int(n))
        for i in range(3):
            u[gd[i]] = v[i]
    return u


def _model_dm():
    model = from_abaqus(P.parse_inp(INP))
    return model, DofManager(model.nodes.keys())


# --------------------------------------------------------------------------- #
# M3.1 / M3.2 -- exporter pure functions (no Abaqus)
# --------------------------------------------------------------------------- #
def test_load_sdv_layout_ok_and_bad():
    meta = load_sdv_layout(os.path.join(_ASSET, "derivative_layout.json"))
    assert meta["element_type"] == "C3D8"
    assert meta["integration"] == "selective_reduced"
    assert meta["sdv_layout"]["parameters"]["E"] == [37, 42]
    with pytest.raises(ExportError):
        load_sdv_layout({"schema": "wrong", "element_type": "C3D8"})


def test_build_derivative_fields_and_guards():
    layout = load_sdv_layout(os.path.join(_ASSET, "derivative_layout.json"))
    good = {str(ip): [0.0] * 48 for ip in range(1, 9)}
    disp = {str(n): [0.0, 0.0, 0.0] for n in range(1, 9)}
    doc = build_derivative_fields(disp, {"1": dict(good)}, layout,
                                  mesh_element_ids=[1])
    assert doc["metadata"]["schema"] == "resasm_derivative_fields_v1"
    assert doc["metadata"]["integration"] == "selective_reduced"
    # missing IP
    bad = {str(ip): [0.0] * 48 for ip in range(1, 8)}          # only 1..7
    with pytest.raises(ExportError):
        build_derivative_fields(disp, {"1": bad}, layout)
    # too-short SDV vector for the layout
    short = {str(ip): [0.0] * 40 for ip in range(1, 9)}
    with pytest.raises(ExportError):
        build_derivative_fields(disp, {"1": short}, layout)
    # element present in mesh but absent from statev
    with pytest.raises(ExportError):
        build_derivative_fields(disp, {"1": dict(good)}, layout,
                                mesh_element_ids=[1, 2])


def test_exported_fixture_is_valid_schema():
    df = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    assert df.element_type == "C3D8"
    assert df.integration == "selective_reduced"
    assert set(df.statev.keys()) == {1}
    assert df.statev[1].shape == (8, 48)
    assert df.metadata["provenance"]["odb"] == "nonuniform_c3d8.odb"


# --------------------------------------------------------------------------- #
# M3.5 -- tangent from exported DDSDDE vs analytic
# --------------------------------------------------------------------------- #
def test_tangent_from_odb_vs_analytic():
    model, dm = _model_dm()
    Xe = model.coords_of(model.elements[1].connectivity)
    df = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    tf, _sf = fields_from_statev(df.statev, df.sdv_layout, parameters=["E"])
    K_odb = kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=tf[1], mode="small")
    K_an = kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=kern.isotropic_D(E0, NU0),
                                mode="small")
    rel = np.max(np.abs(K_odb - K_an)) / np.max(np.abs(K_an))
    print("  M3.5 K_ODB vs analytic rel(inf)=%.3e  (float32 ODB floor ~1e-7)" % rel)
    assert rel < 1e-6                # not 1e-10: SDV is single precision in the ODB


# --------------------------------------------------------------------------- #
# M3.3 -- integration-point ordering (per-IP identifiable stress) + negative
# --------------------------------------------------------------------------- #
def _sigma_bbar_per_ip(Xe, ue):
    """Mean-dilatation (B-bar) stress at the 8 IPs from nodal displacement."""
    Bbar, dets = kern.bbar_b_matrices(Xe)
    D = kern.isotropic_D(E0, NU0)
    return np.array([D @ (Bbar[q] @ ue) for q in range(8)])


def test_ip_order_matches_and_permutation_detected():
    model, dm = _model_dm()
    el = model.elements[1]
    Xe = model.coords_of(el.connectivity)
    edofs = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
    ue = _uvec("U_base.json", dm)[edofs]
    sig_ker = _sigma_bbar_per_ip(Xe, ue)                       # kernel IP order
    S = json.load(open(os.path.join(_FIX, "S_ip.json")))
    sig_abq = np.array([S[str(i + 1)] for i in range(8)])      # Abaqus IP label 1..8

    scale = np.max(np.abs(sig_abq))
    rel_aligned = np.max(np.abs(sig_ker - sig_abq)) / scale
    perm = [0, 6, 2, 3, 4, 5, 1, 7]                            # swap IP2 <-> IP7
    rel_perm = np.max(np.abs(sig_ker - sig_abq[perm])) / scale
    print("  M3.3 aligned rel=%.3e   permuted(2<->7) rel=%.3e" % (rel_aligned, rel_perm))
    assert rel_aligned < 1e-4          # Abaqus IP k == kernel Gauss k (float32)
    assert rel_perm > 1e-2             # a wrong IP order is loudly detected


def test_equilibrium_confirms_ip_order():
    """Independent IP-order check: assembling the exported Abaqus stress with the
    kernel Gauss order reproduces the applied load; a permuted order does not."""
    model, dm = _model_dm()
    el = model.elements[1]
    Xe = model.coords_of(el.connectivity)
    S = json.load(open(os.path.join(_FIX, "S_ip.json")))
    sig = np.array([S[str(i + 1)] for i in range(8)])
    f_int = kern.element_internal_force_bbar(Xe, sig)          # B-bar, kernel order
    # applied Cload (element-local node-major dof indices)
    fext = np.zeros(24)
    for n, comps in {5: {0: 15., 2: 30.}, 6: {0: 25., 2: 10.},
                     7: {0: 5., 2: 40.}, 8: {0: 35., 2: 20.}}.items():
        for c, val in comps.items():
            fext[3 * (n - 1) + c] = val
    free = np.arange(12, 24)                                   # top-face dofs
    rel = np.linalg.norm(f_int[free] - fext[free]) / np.linalg.norm(fext[free])
    sig_perm = sig[[0, 6, 2, 3, 4, 5, 1, 7]]
    rel_perm = np.linalg.norm(kern.element_internal_force_bbar(Xe, sig_perm)[free]
                              - fext[free]) / np.linalg.norm(fext[free])
    print("  M3.3 equilibrium rel=%.3e   permuted rel=%.3e" % (rel, rel_perm))
    assert rel < 1e-5
    assert rel_perm > 1e-1


# --------------------------------------------------------------------------- #
# B-bar vs full integration (formulation)
# --------------------------------------------------------------------------- #
def test_bbar_matches_abaqus_full_int_does_not():
    model, dm = _model_dm()
    el = model.elements[1]
    Xe = model.coords_of(el.connectivity)
    edofs = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
    ue = _uvec("U_base.json", dm)[edofs]
    S = json.load(open(os.path.join(_FIX, "S_ip.json")))
    sig_abq = np.array([S[str(i + 1)] for i in range(8)])
    D = kern.isotropic_D(E0, NU0)
    pts = kern.ABAQUS_C3D8_GAUSS.points
    sig_full = np.array([D @ (kern.b_matrix_reference(Xe, pts[q])[0] @ ue)
                         for q in range(8)])
    sig_bbar = _sigma_bbar_per_ip(Xe, ue)
    rel_full = np.max(np.abs(sig_full - sig_abq)) / np.max(np.abs(sig_abq))
    rel_bbar = np.max(np.abs(sig_bbar - sig_abq)) / np.max(np.abs(sig_abq))
    print("  B-bar vs Abaqus rel=%.3e   full-int vs Abaqus rel=%.3e" % (rel_bbar, rel_full))
    assert rel_bbar < 1e-4             # selective reduced integration matches Abaqus C3D8
    assert rel_full > 0.1              # full integration does NOT


# --------------------------------------------------------------------------- #
# M3.6 -- full recipe vs Abaqus finite-difference reruns
# --------------------------------------------------------------------------- #
def test_recipe_vs_abaqus_finite_difference():
    model, dm = _model_dm()
    tmp = tempfile.mkdtemp(prefix="resasm_m3_")
    try:
        shutil.copy(INP, os.path.join(tmp, "model.inp"))
        shutil.copy(os.path.join(_FIX, "derivative_fields.json"),
                    os.path.join(tmp, "derivative_fields.json"))
        with open(os.path.join(tmp, "sensitivity.yaml"), "w") as fh:
            fh.write("analysis:\n  type: field_residual_sensitivity\n"
                     "  kinematics: small_strain\n"
                     "mesh:\n  format: abaqus\n  file: model.inp\n"
                     "results:\n  format: resasm_derivative_fields_v1\n"
                     "  file: derivative_fields.json\n"
                     "parameters: [E, nu]\n"
                     "assumptions:\n  parameter_independent_geometry: true\n"
                     "  parameter_independent_loads: true\n"
                     "  parameter_independent_boundaries: true\n"
                     "output:\n  directory: results\n")
        res = run_field_recipe(os.path.join(tmp, "sensitivity.yaml"))
        assert res.ok
        manifest = json.load(open(os.path.join(tmp, "results", "manifest.json")))
        assert manifest["integration"] == "selective_reduced"

        d = np.load(os.path.join(tmp, "results", "displacement_sensitivities.npz"),
                    allow_pickle=True)
        du_E_fd = (_uvec("U_Ep.json", dm) - _uvec("U_Em.json", dm)) / (2 * H_E)
        du_nu_fd = (_uvec("U_nup.json", dm) - _uvec("U_num.json", dm)) / (2 * H_NU)

        rel_E = np.max(np.abs(d["du_dE"] - du_E_fd)) / max(np.max(np.abs(du_E_fd)), 1e-30)
        rel_nu = np.max(np.abs(d["du_dnu"] - du_nu_fd)) / max(np.max(np.abs(du_nu_fd)), 1e-30)
        print("  M3.6 du/dnu vs Abaqus-FD rel(inf)=%.3e   ||fd||inf=%.3e" % (rel_nu, np.max(np.abs(du_nu_fd))))
        print("  M3.6 du/dE  vs Abaqus-FD rel(inf)=%.3e (h/E=%.3f) ||fd||inf=%.3e"
              % (rel_E, H_E / E0, np.max(np.abs(du_E_fd))))
        # du/dnu is well-resolved. du/dE is intrinsically ~1e-8, so the FD of two
        # single-precision ODB displacement vectors is precision-limited; this
        # particular step (h/E=0.004) sits just past the FD optimum. The
        # step-size sweep (below) shows du/dE FD reaches < 1e-5 at h/E~0.002, and
        # the accurate elastic oracle (below) pins du/dE to ~6e-8 -- so du/dE is
        # correct; only THIS finite-difference comparison is storage-limited.
        assert rel_nu < 1e-5
        assert rel_E < 5e-5              # storage-precision-limited FD at this step
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_dudE_accurate_elastic_oracle():
    """Precision-independent du/dE check: for linear elasticity u = c/E, so the
    exact sensitivity is du/dE = -u/E. Comparing the residual-method du/dE to
    -u_base/E avoids the finite-difference cancellation entirely and validates the
    E-column to the ODB float32 floor (unlike the FD comparison, which subtracts
    two nearly-equal single-precision displacement vectors)."""
    model, dm = _model_dm()
    df = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    tf, sf = fields_from_statev(df.statev, df.sdv_layout, parameters=["E", "nu"])
    from residual_core.core.field_sensitivity import solve_field_sensitivities
    res = solve_field_sensitivities(
        model=model, tangent_fields=tf, stress_derivative_fields=sf,
        parameters=["E", "nu"], dof_manager=dm, integration="selective_reduced")
    u_base = _uvec("U_base.json", dm)
    oracle = -u_base / E0
    rel = np.max(np.abs(res.du_da("E") - oracle)) / max(np.max(np.abs(oracle)), 1e-30)
    print("  du/dE(residual) vs -u/E (exact elastic oracle) rel(inf)=%.3e" % rel)
    assert rel < 1e-6              # ODB float32 floor -- NOT finite-difference-limited


def test_dudE_fd_step_size_plateau():
    """Demonstrate that the du/dE finite-difference error is a numerical plateau
    (single-precision ODB), not a formulation error: it is truncation-dominated at
    large steps (grows ~ (h/E)^2) and cancellation-dominated at tiny steps, with a
    minimum < 1e-5 near h/E~0.002. Data are real Abaqus E-perturbation reruns."""
    model, dm = _model_dm()
    df = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    tf, sf = fields_from_statev(df.statev, df.sdv_layout, parameters=["E", "nu"])
    from residual_core.core.field_sensitivity import solve_field_sensitivities
    res = solve_field_sensitivities(
        model=model, tangent_fields=tf, stress_derivative_fields=sf,
        parameters=["E", "nu"], dof_manager=dm, integration="selective_reduced")
    du_dE = res.du_da("E")

    sweep = json.load(open(os.path.join(_FIX, "fd_sweep_E.json")))["steps"]

    def _vec(d):
        w = np.zeros(dm.ndof)
        for n, v in d.items():
            gd = dm.node_dofs(int(n))
            for i in range(3):
                w[gd[i]] = v[i]
        return w

    curve = {}
    for h_str, pair in sweep.items():
        h = float(h_str)
        fd = (_vec(pair["Ep"]) - _vec(pair["Em"])) / (2 * h)
        curve[h] = np.max(np.abs(du_dE - fd)) / max(np.max(np.abs(fd)), 1e-30)
    print("  du/dE FD step-size sweep (h/E -> rel):")
    for h in sorted(curve):
        print("     h/E=%.4f  rel=%.3e" % (h / E0, curve[h]))
    best = min(curve.values())
    # the optimum step meets the original 1e-5 target ...
    assert best < 1e-5
    # ... and the error genuinely grows at large steps (truncation), confirming a
    # numerical plateau rather than a constant formulation offset.
    assert curve[max(curve)] > 5 * curve[420.0]


# --------------------------------------------------------------------------- #
# M4 -- OTI-parameter-seeded UMAT: automatic differentiation == analytic
# --------------------------------------------------------------------------- #
# fixtures/derivative_fields_oti.json is the export of the SAME nonuniform C3D8
# job run with elastic_oti_umat.for, which obtains dsigma/dE and dsigma/dnu by
# OTI seeding (OTIM4N1 directions E1, E2) instead of analytically. Because both
# UMATs implement identical constitutive equations, the OTI export must match the
# analytic M3 export to the ODB storage floor (here: bit-for-bit).
def test_oti_umat_matches_analytic_derivatives():
    oti = load_derivative_fields(os.path.join(_FIX, "derivative_fields_oti.json"))
    an = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    so, sa = oti.statev[1], an.statev[1]
    relD = np.max(np.abs(so[:, 0:36] - sa[:, 0:36])) / max(np.max(np.abs(sa[:, 0:36])), 1e-30)
    relE = np.max(np.abs(so[:, 36:42] - sa[:, 36:42])) / max(np.max(np.abs(sa[:, 36:42])), 1e-30)
    relN = np.max(np.abs(so[:, 42:48] - sa[:, 42:48])) / max(np.max(np.abs(sa[:, 42:48])), 1e-30)
    du = max(np.max(np.abs(np.array(oti.displacements[k]) - np.array(an.displacements[k])))
             for k in oti.displacements)
    print("  M4 OTI vs analytic: DDSDDE=%.3e  dsig/dE=%.3e  dsig/dnu=%.3e  U(abs)=%.3e"
          % (relD, relE, relN, du))
    assert relD < 1e-6 and relE < 1e-6 and relN < 1e-6   # OTI AD == analytic derivatives
    assert du < 1e-12                                     # real solution preserved


def test_oti_umat_reproduces_residual_sensitivities():
    from residual_core.core.field_sensitivity import solve_field_sensitivities
    model, dm = _model_dm()
    oti = load_derivative_fields(os.path.join(_FIX, "derivative_fields_oti.json"))
    an = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
    out = {}
    for name, df in (("oti", oti), ("analytic", an)):
        tf, sf = fields_from_statev(df.statev, df.sdv_layout, parameters=["E", "nu"])
        out[name] = solve_field_sensitivities(
            model=model, tangent_fields=tf, stress_derivative_fields=sf,
            parameters=["E", "nu"], dof_manager=dm, integration="selective_reduced")
    for p in ("E", "nu"):
        rel = np.max(np.abs(out["oti"].du_da(p) - out["analytic"].du_da(p))) / \
            max(np.max(np.abs(out["analytic"].du_da(p))), 1e-30)
        print("  M4 du/d%-2s OTI-pipeline vs analytic-pipeline rel(inf)=%.3e" % (p, rel))
        assert rel < 1e-6


# --------------------------------------------------------------------------- #
# gated: regenerate the fixtures from a live Abaqus run (opt-in, slow)
# --------------------------------------------------------------------------- #
_IN_SCRIPT = False        # set True by main() so the gated test can skip cleanly


def _abaqus_regen_enabled():
    return bool(shutil.which(os.environ.get("ABAQUS_CMD", "abaqus"))
                and os.environ.get("RESASM_RUN_ABAQUS"))


def test_live_abaqus_regenerates_fixture():
    """Run the real chain (Abaqus solve -> export) and confirm it reproduces the
    committed fixture to single-precision. Skipped unless RESASM_RUN_ABAQUS=1 and
    Abaqus is on PATH (the offline tests above already validate the pipeline)."""
    if not _abaqus_regen_enabled():
        msg = "set RESASM_RUN_ABAQUS=1 (and have Abaqus on PATH) to run the live chain"
        if _IN_SCRIPT:
            print("  SKIP: %s" % msg)
            return
        pytest.skip(msg)
    import subprocess
    work = tempfile.mkdtemp(prefix="resasm_m3_live_")
    try:
        rc = subprocess.call(["bash", os.path.join(_ASSET, "run_chain.sh"), work])
        assert rc == 0, "run_chain.sh failed (rc=%d)" % rc
        got = load_derivative_fields(os.path.join(work, "derivative_fields.json"))
        ref = load_derivative_fields(os.path.join(_FIX, "derivative_fields.json"))
        rel = np.max(np.abs(got.statev[1] - ref.statev[1])) / \
            max(np.max(np.abs(ref.statev[1])), 1e-30)
        print("  live regen vs fixture statev rel(inf)=%.3e" % rel)
        assert rel < 1e-5              # deterministic solve, single-precision ODB
    finally:
        shutil.rmtree(work, ignore_errors=True)


def test_live_abaqus_oti_umat_regenerates_fixture():
    """Run the live chain with the OTI-seeded UMAT (links the MultiZ_f OTIM4N1
    library) and confirm it reproduces the committed OTI fixture. Skipped unless
    RESASM_RUN_ABAQUS=1 and Abaqus is on PATH; also needs the OTI Fortran lib."""
    oti_lib = os.path.expanduser(
        os.environ.get("OTI_DIR", "~/MultiZ_f/oti")) + "/libotim4n1.a"
    if not (_abaqus_regen_enabled() and os.path.exists(oti_lib)):
        msg = "needs RESASM_RUN_ABAQUS=1, Abaqus on PATH, and the OTI lib (%s)" % oti_lib
        if _IN_SCRIPT:
            print("  SKIP: %s" % msg)
            return
        pytest.skip(msg)
    import subprocess
    work = tempfile.mkdtemp(prefix="resasm_m4_live_")
    try:
        rc = subprocess.call(["bash", os.path.join(_ASSET, "run_chain.sh"), work,
                              "elastic_oti_umat.for"])
        assert rc == 0, "run_chain.sh (OTI) failed (rc=%d)" % rc
        got = load_derivative_fields(os.path.join(work, "derivative_fields.json"))
        ref = load_derivative_fields(os.path.join(_FIX, "derivative_fields_oti.json"))
        rel = np.max(np.abs(got.statev[1] - ref.statev[1])) / \
            max(np.max(np.abs(ref.statev[1])), 1e-30)
        print("  live OTI regen vs OTI fixture statev rel(inf)=%.3e" % rel)
        assert rel < 1e-5
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    global _IN_SCRIPT
    _IN_SCRIPT = True
    tests = [test_load_sdv_layout_ok_and_bad, test_build_derivative_fields_and_guards,
             test_exported_fixture_is_valid_schema, test_tangent_from_odb_vs_analytic,
             test_ip_order_matches_and_permutation_detected,
             test_equilibrium_confirms_ip_order,
             test_bbar_matches_abaqus_full_int_does_not,
             test_recipe_vs_abaqus_finite_difference,
             test_dudE_accurate_elastic_oracle,
             test_dudE_fd_step_size_plateau,
             test_oti_umat_matches_analytic_derivatives,
             test_oti_umat_reproduces_residual_sensitivities,
             test_live_abaqus_regenerates_fixture,
             test_live_abaqus_oti_umat_regenerates_fixture]
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
