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

from history_support import CONTRACTS, EXAMPLE, RA_ROOT, umat_repository

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
    assert "Command executed: yes: resasm history --model " in report
    for line in ("Status: executed successfully", "Equilibrium passed: yes", "Tangent available: yes",
                 "Derivative calculated: yes", "Unsupported feature detected: none",
                 "Abaqus comparison available: yes (primal)", "Public and private outputs separated: yes"):
        assert line in report
    assert results["scope"]["increments_replayed"] == 10
    assert results["scope"]["integration_points"] == 768
    assert results["scope"]["parameters"] == ["E", "nu", "SIGY0", "H"]
    private = json.loads((out / "private" / "run_details.json").read_text())
    assert len(private["increments"]) == 10
    assert _listed_public_files(report) == {p.name for p in out.iterdir() if p.is_file()}


def _listed_public_files(report):
    """The public files named by the report's 'Public and private outputs' line."""
    import ast
    line = next(line for line in report.splitlines()
                if line.startswith("Public and private outputs separated: yes: public "))
    listed = line.split("public ", 1)[1].rsplit("; private private/", 1)[0]
    return set(ast.literal_eval(listed))


def test_the_report_lists_every_public_file_written(provider_factory, tmp_path):
    """Regression: the line named only the three standard files, not
    sensitivity_shares.csv or fields.npz when a request made them public."""
    obj, _, _ = provider_factory("m3_j2")
    request = json.loads((BEAM / "sensitivity_request.json").read_text())
    request["full_field"] = True
    (tmp_path / "request.json").write_text(json.dumps(request))
    out = tmp_path / "results"
    assert resasm(["history", "--model", str(BEAM / "Analysis.inp"), "--fields", str(BEAM / "fields.npz"),
                   "--material", str(obj), "--request", str(tmp_path / "request.json"),
                   "--out", str(out)]) == 0
    written = {p.name for p in out.iterdir() if p.is_file()}
    assert {"sensitivity_shares.csv", "fields.npz"} <= written
    assert _listed_public_files((out / "run_report.txt").read_text()) == written


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


def test_every_increment_obeys_the_homogeneity_identity(provider_factory, tmp_path):
    """Euler's theorem for the J2 model, elastic AND plastic increments.

    With linear isotropic hardening the J2 stress update is homogeneous of
    degree one in (E, SIGY0, H) at fixed nu: scaling the three by lambda scales
    every stress for the same strain history, and the yield test f = q -
    (SIGY0 + H*eqplas) scales with it, so no branch changes. Under prescribed
    displacements the equilibrium displacements therefore do not change at
    all, and stresses and reactions scale by lambda. Differentiating at
    lambda = 1 gives, at EVERY increment,

        E dQ/dE + SIGY0 dQ/dSIGY0 + H dQ/dH = Q    (Q = RF, S, MISES)
        E dQ/dE + SIGY0 dQ/dSIGY0 + H dQ/dH = 0    (Q = U, eqplas)

    an exact identity the sensitivities must satisfy whatever the history, and
    one that no part of the engine uses. Re-equilibrated to double precision
    (free residual <= 1e-10 of its scale, pinned above), the identity holds to
    roundoff: measured 9.4e-15 of the largest term on 2026-09-18. The bound
    1e-11 leaves three orders for platform differences and still catches a
    single derivative wrong by 1e-6 (that alone leaves a residual of about
    1e-9 here). The recorded-state replay satisfies it only to the size of the
    ODB's single-precision residual (8e-6 measured), which is why this runs
    with --reequilibrate.
    """
    obj, _, _ = provider_factory("m3_j2")
    request = tmp_path / "request.json"
    request.write_text(json.dumps({
        "outputs": [
            {"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum", "domain": {"nset": "TIP"}},
            {"name": "midtop_U2", "field": "U", "component": 2, "reduction": "component", "domain": {"nodes": [59]}},
            {"name": "e1_ip1_S11", "field": "S", "component": 1, "reduction": "component",
             "domain": {"elements": [1], "points": [1]}},
            {"name": "mises_mean", "field": "MISES", "component": 1, "reduction": "volume_mean",
             "domain": {"elements": "ALL"}},
            {"name": "mises_root_max", "field": "MISES", "component": 1, "reduction": "max",
             "domain": {"elset": "ROOTEL"}},
            {"name": "eqplas_max", "field": "SDV", "component": 1, "reduction": "max", "domain": {"elements": "ALL"}}],
        "parameters": "ALL", "domain": {"nodes": "ALL", "elements": "ALL"}, "increments": "ALL"}))
    out = tmp_path / "results"
    assert resasm(["history", "--model", str(BEAM / "Analysis.inp"), "--fields", str(BEAM / "fields.npz"),
                   "--material", str(obj), "--request", str(request), "--out", str(out),
                   "--reequilibrate"]) == 0
    results = json.loads((out / "sensitivity_results.json").read_text())
    degree = {"RF": 1, "S": 1, "MISES": 1, "U": 0, "SDV": 0}
    plastic_increments = set()
    for row in results["results"]:
        weighted = [row["weighted"][name] for name in ("E", "SIGY0", "H")]
        largest = max(max(abs(w) for w in weighted), abs(row["value"]))
        residual = sum(weighted) - degree[row["field"]] * row["value"]
        assert abs(residual) <= 1e-11 * largest, (row["output"], row["increment"], residual, largest)
        if row["output"] == "eqplas_max" and row["value"] > 0:
            plastic_increments.add(row["increment"])
    # the identity must have been exercised where plasticity is active
    assert len(plastic_increments) >= 5, sorted(plastic_increments)


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
    built = build_provider(umat_repository() / CONTRACTS["m3_j2"][0], tmp_path / "provider")
    bounded = Namespace(model=RA_ROOT / "examples/presentation_request/Analysis.inp", material=built["object"],
                        mapping=None, request=RA_ROOT / "examples/presentation_request/sensitivity_request.json")
    assert bounded_scope_reason(bounded) is None
    beam = Namespace(model=BEAM / "Analysis.inp", material=built["object"], mapping=None,
                     request=BEAM / "sensitivity_request.json")
    assert bounded_scope_reason(beam)          # tight *Controls, nonzero BCs, sets, MISES ...
    # An input neither engine can read is the bounded engine's to report: it
    # owns the presentation interface's categorised, private diagnostics, and
    # test_presentation_request pins what it says. Only a READABLE model that
    # is out of the bounded scope is forwarded to the history engine.
    for key, value in (("model", "/nonexistent.inp"), ("material", "/nonexistent.obj")):
        assert bounded_scope_reason(Namespace(**{**vars(bounded), key: value})) is None


def test_another_provider_on_a_bounded_deck_goes_to_the_history_engine(provider_factory, tmp_path):
    """Regression: the one-element deck, which the bounded reader reads, with the
    real FCC provider and its own valid mapping stayed with the bounded engine,
    which then refused the material (Category: material_mapping, exit 2)
    instead of handing the model to the history engine."""
    from residual_core.ui.cmd_history import bounded_scope_reason
    fcc, fcc_contract, _ = provider_factory("m6_fcc")
    deck = RA_ROOT / "examples/presentation_request"
    args = Namespace(model=deck / "Analysis.inp", material=fcc, mapping=None,
                     request=deck / "sensitivity_request.json")
    reason = bounded_scope_reason(args)
    assert reason and "not the fingerprint-pinned m3_j2 J2 provider" in reason
    assert fcc_contract["model_id"] in reason
    # a damaged mapping of the pinned provider is still the bounded engine's to diagnose
    j2, j2_contract, _ = provider_factory("m3_j2")
    damaged = json.loads(json.dumps(j2_contract))
    damaged["layouts"]["voigt"].reverse()
    (tmp_path / "Mapping.json").write_text(json.dumps(damaged))
    assert bounded_scope_reason(Namespace(**{**vars(args), "material": j2,
                                             "mapping": tmp_path / "Mapping.json"})) is None


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
