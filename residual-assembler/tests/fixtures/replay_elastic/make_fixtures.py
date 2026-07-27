"""FIXTURE BUILDER (test-only) — simulate the JHU package + a production analysis.

This is NOT part of the collaborator application. It stands in for two things a
real collaborator would receive/produce out-of-band:

  * the JHU-shipped package: an opaque OTI binary (.so) + material_manifest.json,
    built here from the isolated opaque-provider fixture;
  * the saved production analysis: converged replay records at the base and at
    p +/- dp for each parameter (simulated with an ordinary elastic solve, which
    stands in for running JHU's REGULAR .obj in Abaqus).

The Program-2 consumer under test then reads ONLY {binary, manifest, record,
request} through the public API — it never imports this module. The perturbed
records let the test finite-difference du/dp and dq/dp from analysis DATA, with
no assumption about the constitutive law inside the binary.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess

import numpy as np

from residual_core.core.model import Element, Model
from residual_core.core.dof_manager import DofManager
from residual_core.core.constraints import partition
from residual_core.formulations import c3d8_kernel as kern
from residual_core.replay.record import _BC

HERE = os.path.dirname(os.path.abspath(__file__))
PROVIDER_SRC = os.path.join(HERE, "opaque_provider", "elastic_reference.f90")
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # tests/fixtures/replay_elastic -> root
_CONTRACT = os.path.join(_REPO_ROOT, "residual_core", "replay", "contract", "CONTRACT_VERSION.json")
MODEL_ID = "reference_elastic_isotropic"
E0, NU0, TOTAL_FZ, UZ_TOP = 210000.0, 0.3, 1000.0, 0.01

_NODES = {
    1: (0., 0., 0.), 2: (1., 0., 0.), 3: (1., 1., 0.), 4: (0., 1., 0.),
    5: (0., 0., 1.), 6: (1., 0., 1.), 7: (1., 1., 1.), 8: (0., 1., 1.),
    9: (0., 0., 2.), 10: (1., 0., 2.), 11: (1., 1., 2.), 12: (0., 1., 2.),
}
_ELEMS = {1: [1, 2, 3, 4, 5, 6, 7, 8], 2: [5, 6, 7, 8, 9, 10, 11, 12]}
_BASE, _TOP = (1, 2, 3, 4), (9, 10, 11, 12)


def _model(drive):
    m = Model(nodes=dict(_NODES),
              elements={e: Element(e, "C3D8", list(c)) for e, c in _ELEMS.items()})
    m.element_formulation = {e: "solid_c3d8_small_strain" for e in _ELEMS}
    m.element_material = {e: "MAT" for e in _ELEMS}
    bcs = [_BC(target=n, kind="ENCASTRE") for n in _BASE]
    if drive == "disp":
        bcs += [_BC(target=n, dof_start=3, dof_end=3, value=UZ_TOP, kind="value") for n in _TOP]
    m.boundaries = bcs
    return m, DofManager(m.nodes.keys())


def _bc_json(drive):
    bcs = [{"target": n, "kind": "ENCASTRE"} for n in _BASE]
    if drive == "disp":
        bcs += [{"target": n, "dof": 3, "value": UZ_TOP, "kind": "value"} for n in _TOP]
    return bcs


def _cload(drive):
    return [{"node": n, "dof": 3, "value": TOTAL_FZ / len(_TOP)} for n in _TOP] if drive == "force" else []


def _fext(dm, drive):
    f = np.zeros(dm.ndof)
    for c in _cload(drive):
        f[dm.node_dofs(c["node"])[c["dof"] - 1]] += c["value"]
    return f


def _solve(m, dm, E, nu, fext):
    D = kern.isotropic_D(E, nu)
    K = np.zeros((dm.ndof, dm.ndof))
    for el in m.elements.values():
        Xe = m.coords_of(el.connectivity)
        ed = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
        K[np.ix_(ed, ed)] += kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=D, mode="small")
    free, pres_idx, pres = partition(m, dm)
    u = np.zeros(dm.ndof)
    for d, v in pres.items():
        u[d] = v
    rhs = fext[free] - (K[np.ix_(free, pres_idx)] @ u[pres_idx] if pres_idx.size else 0.0)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)
    return u, K


def _stress_and_reactions(m, dm, u, E, nu, fext):
    D = kern.isotropic_D(E, nu)
    pts = kern.ABAQUS_C3D8_GAUSS.points
    sig, fint = {}, np.zeros(dm.ndof)
    for el in m.elements.values():
        Xe = m.coords_of(el.connectivity)
        ed = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
        s = np.array([D @ (kern.b_matrix_reference(Xe, p)[0] @ u[ed]) for p in pts])
        sig[el.eid] = s
        fint[ed] += kern.element_internal_force_small_strain(Xe, s)
    return sig, (fint - fext)          # reactions = internal - external


def _record(drive, E, nu, regular_hash):
    m, dm = _model(drive)
    fext = _fext(dm, drive)
    u, _ = _solve(m, dm, E, nu, fext)
    sig, react = _stress_and_reactions(m, dm, u, E, nu, fext)
    return {
        "schema": "resasm_replay_record_v1",
        "provenance": {"model_id": MODEL_ID, "regular_hash": regular_hash,
                       "units": "MPa,mm",
                       "analysis_hash": hashlib.sha256(u.tobytes()).hexdigest()[:16]},
        "kinematics": "small_strain",
        "material": {"props": [E, nu]},
        "mesh": {"nodes": {str(k): list(v) for k, v in _NODES.items()},
                 "elements": {str(e): {"type": "C3D8", "connectivity": list(c)}
                              for e, c in _ELEMS.items()}, "ip_order": "abaqus_c3d8"},
        "dof": {"map": "node_major_uxuyuz"},
        "boundaries": _bc_json(drive), "loads": {"cload": _cload(drive)},
        "increments": [{"time": 1.0, "dtime": 1.0, "u": u.tolist(),
                        "stress_ip": {str(e): s.tolist() for e, s in sig.items()},
                        "reactions": react.tolist()}],
    }


def build_fixtures(outdir: str, drive: str = "disp", hE: float = None, hnu: float = 1e-6):
    """Build {opaque .so, material_manifest.json, base + perturbed records} into
    outdir. Returns a dict of paths + the FD steps used."""
    os.makedirs(outdir, exist_ok=True)
    hE = hE if hE is not None else 1e-6 * E0
    so = os.path.join(outdir, "libmat.so")
    subprocess.check_call(["gfortran", "-shared", "-fPIC", "-O2", PROVIDER_SRC, "-o", so])
    src_hash = hashlib.sha256(open(PROVIDER_SRC, "rb").read()).hexdigest()[:16]
    so_hash = hashlib.sha256(open(so, "rb").read()).hexdigest()[:16]
    contract_ver = json.load(open(_CONTRACT)).get("combined_hash", "")

    manifest = {
        "schema": "resasm_material_package_v1", "model_id": MODEL_ID,
        "kinematics": "small_strain", "ntens": 6, "nprops": 2, "nstatev": 0,
        "parameters": [{"name": "E", "index": 1, "oti_direction": 1, "units": "MPa"},
                       {"name": "nu", "index": 2, "oti_direction": 2, "units": "-"}],
        "state_layout": [], "outputs": ["stress", "consistent_tangent"],
        "sensitivity_capabilities": {"runtime_parameter_seeding": True,
                                     "maximum_directions": 2, "order": 1},
        "abi": {"symbol": "mat_eval_v1", "version": 1},
        "contract_version": contract_ver,
        "binaries": {"regular": {"hash": src_hash, "build_id": "fixture"},
                     "oti": {"path": "libmat.so", "hash": so_hash, "build_id": "fixture"}},
        "provenance": {"generated_by": "tests/fixtures/replay_elastic/make_fixtures.py"},
    }
    mpath = os.path.join(outdir, "material_manifest.json")
    json.dump(manifest, open(mpath, "w"), indent=2)

    paths = {"so": so, "manifest": mpath, "drive": drive, "hE": hE, "hnu": hnu}
    cases = {"base": (E0, NU0), "Ep": (E0 + hE, NU0), "Em": (E0 - hE, NU0),
             "nup": (E0, NU0 + hnu), "num": (E0, NU0 - hnu)}
    for name, (E, nu) in cases.items():
        p = os.path.join(outdir, "record_%s.json" % name)
        json.dump(_record(drive, E, nu, src_hash), open(p, "w"), indent=2)
        paths[name] = p
    return paths
