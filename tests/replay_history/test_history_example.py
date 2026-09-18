"""The committed Abaqus example, the command line, routing and a cross-check.

``examples/replay_history/j2_beam`` holds a real Abaqus 2021.HF5 run (12x4x2
C3D8 J2 beam, tip pushed 0.08 mm in 10 increments, tight convergence controls):
the deck, the single-precision ODB field export and central differences of
perturbed Abaqus reruns at h = 1e-2, 5e-3, 2e-3, 1e-3
(``regenerate_example.py`` rebuilds all of it). The replay runs offline.
"""
import json
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from residual_core.replay.history import HistoryEngine, run_history
from residual_core.replay.history_inputs import load_recorded_fields, read_history_model
from residual_core.ui.cli import main as resasm

from conftest import EXAMPLE, RA_ROOT

BEAM = EXAMPLE / "j2_beam"
EPS32 = float(np.finfo(np.float32).eps)


@pytest.fixture(scope="module")
def replayed(provider_factory, tmp_path_factory):
    obj, contract, _ = provider_factory("m3_j2")
    out = tmp_path_factory.mktemp("example") / "results"
    code = resasm(["history", "--model", str(BEAM / "Analysis.inp"), "--fields", str(BEAM / "fields.npz"),
                   "--material", str(obj), "--request", str(BEAM / "sensitivity_request.json"),
                   "--out", str(out)])
    assert code == 0
    return out, json.loads((out / "sensitivity_results.json").read_text())


def test_example_is_small():
    total = sum(path.stat().st_size for path in BEAM.iterdir())
    assert total < 400_000, total


def test_public_private_outputs_and_report(replayed):
    out, results = replayed
    assert {p.name for p in out.iterdir()} == {"sensitivity_results.json", "sensitivity_tables.csv",
                                               "run_report.txt", "sensitivity_shares.csv", "private"}
    report = (out / "run_report.txt").read_text()
    for line in ("Status: executed successfully", "Equilibrium passed: yes", "Tangent available: yes",
                 "Derivative calculated: yes", "Unsupported feature detected: none",
                 "Abaqus comparison available: yes (primal)", "Public and private outputs separated: yes"):
        assert line in report
    assert results["scope"]["increments_replayed"] == 10
    assert results["scope"]["integration_points"] == 768
    assert results["scope"]["parameters"] == ["E", "nu", "SIGY0", "H"]
    private = json.loads((out / "private" / "run_details.json").read_text())
    assert len(private["increments"]) == 10


def test_replayed_reaction_equals_the_odb(replayed):
    _, results = replayed
    data = np.load(BEAM / "fields.npz")
    tip = np.isclose(data["coords"][:, 0], data["coords"][:, 0].max())
    odb = data["RF"][:, tip, 1].astype(float).sum(axis=1)
    for row in results["results"]:
        if row["output"] == "tip_RF2":
            n = row["increment"]
            assert abs(row["value"] - odb[n]) <= 1e-5 * abs(odb[n])


def test_sensitivities_within_the_abaqus_fd_uncertainty(replayed):
    """Every output, parameter and increment: |OTI - FD| <= spread + float32 floor of the plateau pair."""
    _, results = replayed
    fd = json.loads((BEAM / "abaqus_fd.json").read_text())
    rows = {(r["output"], r["increment"]): r for r in results["results"]}
    steps = [str(h) for h in fd["steps"]]
    checked = zero = 0
    for output in ("tip_RF2", "midtop_U2", "tiptop_U1", "e1_ip1_S11"):
        scale = np.abs(np.array(fd["nominal"][output])).max()
        for j, name in enumerate(fd["parameters"]):
            p = abs(fd["props"][j])
            for n in range(1, 11):
                estimates = [fd["fd"][name][h][output][n] for h in steps]
                spread, k = min((abs(estimates[i] - estimates[i + 1]), i) for i in range(len(steps) - 1))
                reference = estimates[k + 1]
                ours = rows[(output, n)]["derivatives"][name]
                floor = EPS32 * abs(fd["nominal"][output][n]) / (2 * float(steps[k + 1]) * p)
                if abs(reference) * p / scale < 1e-4:          # zero reference (e.g. du/dE, elastic)
                    assert abs(ours - reference) * p / scale < 1e-4, (output, name, n)
                    zero += 1
                else:
                    assert abs(ours - reference) <= spread + floor, (output, name, n, ours, reference)
                    checked += 1
    assert checked + zero == 160 and checked >= 110          # 4 outputs x 4 parameters x 10 increments


def test_elastic_increments_obey_the_scaling_laws(replayed):
    """Increments 1-4 are elastic (no EQPLAS): dRF/dE = RF/E and d sigma/dE = sigma/E."""
    _, results = replayed
    rows = {(r["output"], r["increment"]): r for r in results["results"]}
    for n in range(1, 5):
        for output in ("tip_RF2", "e1_ip1_S11"):
            row = rows[(output, n)]
            assert row["derivatives"]["E"] * 2e5 == pytest.approx(row["value"], rel=2e-5)
            assert row["derivatives"]["SIGY0"] == 0.0 and row["derivatives"]["H"] == 0.0


def test_weighted_shares(replayed):
    _, results = replayed
    shares = results["weighted_shares"]
    for row in shares["rows"]:
        assert sum(row["field_share_percent"].values()) == pytest.approx(100.0)
        assert sum(row["scalar_share_percent"].values()) == pytest.approx(100.0)
    first = shares["rows"][0]["field_share_percent"]
    assert first["E"] > 90 and first["SIGY0"] == 0.0
    last = shares["rows"][-1]["field_share_percent"]
    assert last["SIGY0"] > last["E"]                       # plasticity takes over


def test_reequilibration_changes_little_on_a_tightly_converged_odb(provider_factory):
    _, _, material = provider_factory("m3_j2")
    model = read_history_model(BEAM / "Analysis.inp")
    fields = load_recorded_fields(BEAM / "fields.npz", model)
    engine = HistoryEngine(model, material)
    recorded = run_history(engine, fields=fields)
    polished = run_history(engine, fields=fields, reequilibrate=True)
    for a, b in zip(recorded.increments, polished.increments):
        assert b.residual_free_max <= 1e-10 * b.residual_scale
        assert np.abs(a.du - b.du).max() <= 1e-4 * np.abs(b.du).max()
        assert b.correction <= 1e-6


def test_routing_keeps_bounded_models_and_sends_the_rest_to_the_history_engine(provider_factory, tmp_path):
    from residual_core.ui.cmd_history import bounded_scope_reason
    from umat_oti.provider import build_provider
    from conftest import CONTRACTS, umat_repository
    built = build_provider(umat_repository() / CONTRACTS["m3_j2"][0], tmp_path / "provider")
    bounded = Namespace(model=RA_ROOT / "examples/presentation_request/Analysis.inp", material=built["object"],
                        mapping=None, request=RA_ROOT / "examples/presentation_request/sensitivity_request.json")
    assert bounded_scope_reason(bounded) is None
    beam = Namespace(model=BEAM / "Analysis.inp", material=built["object"], mapping=None,
                     request=BEAM / "sensitivity_request.json")
    assert bounded_scope_reason(beam)          # tight *Controls, nonzero BCs, sets, MISES ...
    for key, value in (("model", "/nonexistent.inp"), ):
        assert bounded_scope_reason(Namespace(**{**vars(bounded), key: value}))


def test_cross_check_with_the_bounded_j2_engine(provider_factory, tmp_path):
    """The bounded example solved by both engines (J2-specific chain vs UMAT_OTI_EVAL_TOTAL)."""
    from residual_core.replay.connected import solve_history
    from residual_core.replay.path_material import PathMaterial
    from residual_core.replay.presentation_inputs import read_model, target_nodes
    from residual_core.replay.record import ReplayRecord
    obj, contract, material = provider_factory("m3_j2")
    deck = RA_ROOT / "examples/presentation_request/Analysis.inp"
    model = read_model(deck)
    raw = {"schema": "resasm_replay_record_v1", "kinematics": "small_strain", "integration": "selective_reduced",
           "mesh": {"nodes": {str(node): list(coords) for node, coords in model.nodes.items()},
                    "elements": {"1": {"type": "C3D8", "connectivity": list(range(1, 9))}}},
           "material": {"props": [210000., .3, 250., 2000.]},
           "boundaries": [{"target": node, "dof": boundary.dof_start, "value": 0.}
                          for boundary in model.boundaries for node in target_nodes(model, boundary.target)],
           "loads": {"cload": [{"node": node, "dof": 1, "value": 75.} for node in (2, 3, 6, 7)]},
           "increments": [{"dt": .25, "load_factor": number / 4.} for number in range(1, 5)],
           "provenance": {"regular_source_hash": contract["regular_source_hash"]}}
    theirs = solve_history(ReplayRecord(raw), PathMaterial(str(obj), contract, str(tmp_path / "path")))
    engine = HistoryEngine(read_history_model(deck), material)
    ours = run_history(engine, times=np.linspace(0, 1, 5), rtol=1e-13)
    assert (ours.increments[-1].state > 0).any()
    for a, b in zip(theirs["increments"], ours.increments):
        np.testing.assert_allclose(b.u, np.asarray(a["u"]), rtol=0, atol=1e-12 * np.abs(b.u).max())
        du = np.asarray(a["du_dp"])
        np.testing.assert_allclose(b.du, du, rtol=0, atol=1e-9 * np.abs(du).max())
