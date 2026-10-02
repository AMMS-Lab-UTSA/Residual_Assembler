"""The replay's finite-strain C3D8 operators are Abaqus's: mean dilatation.

Evidence: Abaqus 2021 job cc_noether_nlgeom_probe (corpus_campaign/batches/B2/
noether/abaqus_probe): one C3D8 with a DISTORTED reference shape, every node
driven through a 35-degree rotation about (1,2,3) and an inhomogeneous
stretch, 5 increments; a probe UMAT wrote DFGRD0/DFGRD1 at full precision,
reactions read from the ODB (float32). The frozen increment (3) is shared
with tests/formulations/test_finite_strain_history_abaqus_contract.py.

Pinned through the REPLAY API (residual_core.replay.kinematics):

* DFGRD0/1 = F_bar with the element MEAN volume change (to 1e-14);
* the reactions are the mean-dilatation B-bar force with weight w0 Jbar
  (to 2e-7, the ODB's resolution);
* the centroid form of commit 30a5ab4 misses both measurably;
* on a distorted 3x3x3 mesh under a linear field (uniform stress), interior
  nodes are in equilibrium to round-off with mean dilatation; the centroid
  form leaves a residual -- it fails the patch test.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from residual_core.replay import kinematics

pytestmark = pytest.mark.unit

_SOURCE = (Path(__file__).resolve().parents[1] / "formulations"
           / "test_finite_strain_history_abaqus_contract.py")
_spec = importlib.util.spec_from_file_location("_abaqus_probe_fixture", _SOURCE)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
PROBE = _module._PROBE


def _abaqus(name):
    a = np.array(PROBE[name])
    return a.reshape(8, 3, 3).transpose(0, 2, 1) if a.shape[1] == 9 else a


def _state():
    X = np.array(PROBE["X"])
    Uf = np.array(PROBE["U_final"])
    k, n = PROBE["k"], PROBE["increments"]
    return X, Uf * (k - 1) / n, Uf * k / n


def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(b)))


def test_dfgrd0_and_dfgrd1_are_the_mean_dilatation_fbar():
    X, U0, U1 = _state()
    assert _rel(kinematics.deformation_gradients(X[None], U0[None])[0], _abaqus("DFGRD0")) < 1e-14
    assert _rel(kinematics.deformation_gradients(X[None], U1[None])[0], _abaqus("DFGRD1")) < 1e-14


def test_the_reactions_are_the_mean_dilatation_bbar_force():
    X, _, U1 = _state()
    stress = np.array(PROBE["STRESSOUT"])
    force = kinematics.internal_force(X[None], U1[None], stress[None])[0].reshape(8, 3)
    assert _rel(force, np.array(PROBE["RF"])) < 2e-7
    B, w = kinematics.spatial_operators(X[None], U1[None])
    assert _rel(np.einsum("qai,qa,q->i", B[0], stress, w[0]).reshape(8, 3), force) < 1e-14


def test_the_centroid_form_is_measurably_not_abaqus():
    X, _, U1 = _state()
    bar = kinematics.deformation_gradients(X[None], U1[None], volume="centroid")[0]
    assert _rel(bar, _abaqus("DFGRD1")) > 1e-4
    force = kinematics.internal_force(X[None], U1[None], np.array(PROBE["STRESSOUT"])[None],
                                      integration="centroid")[0].reshape(8, 3)
    assert _rel(force, np.array(PROBE["RF"])) > 1e-3       # 6.8e-3 at this increment (2.1e-2 worst over the run)


def _patch_residual(integration):
    from residual_core.corpus.mesh import brick
    mesh = brick((3, 3, 3), distort=0.15)
    H = 0.05 * np.array([[1.0, 0.3, -0.2], [0.1, -0.4, 0.25], [-0.15, 0.2, 0.5]])
    U = mesh.coords @ H.T
    Xe = mesh.coords[mesh.conn]
    Ue = U[mesh.conn]
    stress = np.tile([120.0, -35.0, 60.0, 14.0, -22.0, 9.0], (len(mesh.conn), 8, 1))
    fe = kinematics.internal_force(Xe, Ue, stress, integration=integration)
    R = np.zeros(mesh.ndof)
    dofs = (3 * mesh.conn[:, :, None] + np.arange(3)).reshape(len(mesh.conn), 24)
    np.add.at(R, dofs, fe)
    scale = np.zeros(mesh.ndof)
    np.add.at(scale, dofs, np.abs(fe))
    interior = [a for a, x in enumerate(mesh.coords)
                if all(1e-9 < x[d] < mesh.size[d] - 1e-9 for d in range(3))]
    idx = np.array([3 * a + d for a in interior for d in range(3)])
    return float(np.max(np.abs(R[idx])) / np.max(scale))


def test_a_distorted_mesh_passes_the_patch_test_with_mean_dilatation():
    assert _patch_residual("selective_reduced") < 1e-14


def test_the_centroid_form_fails_the_patch_test_on_a_distorted_mesh():
    assert _patch_residual("centroid") > 1e-6


# ------------------------------------------------- R2 (B3 review): deck_replay
# core.deck_replay built NLGEOM STRAN and DFGRD1 from the plain F; Abaqus's
# C3D8 hands the UMAT the mean-dilatation ones. Pinned against the probe's own
# STRAN at the frozen increment k: the end-of-increment strain
# STRAN_k + DSTRAN_k (what deck_replay returns, accumulated from the undeformed
# state), and the start-of-increment STRAN_k = DROT_k e_(k-1) DROT_k^T.

def _replayed(k, integration="mean_dilatation"):
    from residual_core.core import deck_replay
    X, Uf = np.array(PROBE["X"]), np.array(PROBE["U_final"])
    history = [Uf * i / PROBE["increments"] for i in range(1, k + 1)]
    return deck_replay.hughes_winget_strain_at_points(X, history, integration=integration)


def test_deck_replay_nlgeom_stran_is_abaqus_mean_dilatation_stran():
    k = PROBE["k"]
    end = np.array(PROBE["STRAN"]) + np.array(PROBE["DSTRAN"])
    assert _rel(_replayed(k), end) < 1e-14
    from residual_core.formulations import c3d8_nlgeom as nl
    start = nl.strain_tensor_to_voigt(nl.rotate(_abaqus("DROT"),
                                                nl.strain_voigt_to_tensor(_replayed(k - 1))))
    assert _rel(start, np.array(PROBE["STRAN"])) < 1e-14


def test_the_plain_f_stran_reading_is_measurably_not_abaqus():
    k = PROBE["k"]
    end = np.array(PROBE["STRAN"]) + np.array(PROBE["DSTRAN"])
    assert _rel(_replayed(k, "full"), end) > 1e-2
    trace = _replayed(k)[:, :3].sum(axis=1)
    assert np.ptp(trace) < 1e-15                      # one volumetric strain per element


def test_deck_replay_dfgrd1_is_the_mean_dilatation_fbar():
    from residual_core.core import deck_replay
    X, _, U1 = _state()
    assert _rel(deck_replay.deformation_gradients_at_points(X, U1), _abaqus("DFGRD1")) < 1e-14
    assert _rel(deck_replay.deformation_gradients_at_points(X, U1, integration="full"),
                _abaqus("DFGRD1")) > 1e-3
