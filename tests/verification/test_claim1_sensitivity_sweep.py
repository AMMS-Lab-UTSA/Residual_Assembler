"""Slide 13: both programs against centred FD of the ORIGINAL UMAT (claim 1)."""

import json

import pytest

from verification import claim1_sensitivity_sweep as claim1

pytestmark = [pytest.mark.integration, pytest.mark.fortran]


def _run(out, work, *arguments):
    assert claim1.main([*arguments, "--out", str(out), "--work", str(work)]) == 0
    return json.loads((out / "claim1_sensitivity_sweep.json").read_text())


def test_stateless_plastic_and_kinked_models_agree_with_fd(out_dirs):
    out, work = out_dirs
    payload = _run(out, work, "--model", "m1_elastic", "--model", "m3_j2",
                   "--model", "sweep_drucker_prager", "--increments", "60", "--skip-umat-sweep")
    rows = {row["model"]: row for row in payload["summary"]["rows"]}
    models = {record["model"]: record for record in payload["models"]}
    assert set(rows) == {"m1_elastic", "m3_j2", "sweep_drucker_prager"}
    for name, row in rows.items():
        # the slide's two metrics, and the stricter max-over-path ones
        assert row["p1_dsigma_dp"] < 1e-5, name
        assert row["p2_dsigmavm_dp"] < 1e-5, name
        assert row["p1_max_over_path"] < 1e-5, name
        assert row["p2_max_over_path"] < 1e-5, name
        assert models[name]["primal_parity_scaled"] < 1e-12, name
    # plasticity is actually exercised on this path
    assert models["m3_j2"]["inelastic_increments"] > 30
    # MARCH and the EVAL replay are two entry points of one object: identical
    assert models["m3_j2"]["crosscheck_march_vs_eval_replay"] < 1e-12
    # Drucker-Prager first yields exactly at an increment boundary; the reference
    # there is the one-sided stencil on the nominal branch, and it is listed
    per_parameter = models["sweep_drucker_prager"]["fd_reference"]["per_parameter"]
    assert all(flags for p in per_parameter.values() for flags in p["one_sided_increments"].values())


def test_umat_sweep_tool_and_transform_driver_cross_check(out_dirs):
    out, work = out_dirs
    payload = _run(out, work, "--model", "m5_cpflow", "--increments", "40")
    sweep = payload["umat_sweep_tool"]["m5_cpflow"]
    assert sweep["verdict"] == "succeeded"
    assert sorted(sweep["directions_verified"]) == sorted(sweep["parameters"])
    record = payload["models"][0]
    # Program 1's own transform driver and the compiled provider are separate
    # builds; on the identical path their DSIGMA_DP agree to round-off
    assert record["crosscheck_transform_driver_vs_replay"]["max_abs_difference_scaled"] < 1e-12
    assert payload["summary"]["models_both_below_1e-5"] == 1


@pytest.mark.slow
def test_all_twenty_models(out_dirs):
    out, work = out_dirs
    payload = _run(out, work)
    summary = payload["summary"]
    assert summary["models_attempted"] == 20
    assert summary["models_measured"] == 20
    assert summary["models_both_below_1e-5"] == 20
    assert summary["parameter_directions_in_passing_models"] == 84
    assert payload["umat_sweep_tool"]["_run"]["funnel"]["attempted"] == 20
