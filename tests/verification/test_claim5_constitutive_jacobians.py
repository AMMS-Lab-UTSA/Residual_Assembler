"""Slide 25: internal constitutive Jacobians of the ICP family (claim 5)."""

import json

import pytest

from verification import claim5_constitutive_jacobians as claim5

pytestmark = [pytest.mark.integration, pytest.mark.fortran]


def test_oti_agrees_with_fd_and_audits_the_hand_coded_jacobians(out_dirs):
    out, work = out_dirs
    assert claim5.main(["--model", "UMAT_VPDCL", "--model", "UMAT_NKH_1.02",
                        "--out", str(out), "--work", str(work)]) == 0
    payload = json.loads((out / "claim5_constitutive_jacobians.json").read_text())
    records = {r["model"]: r for r in payload["records"]}
    for record in records.values():
        assert record["status"] == "measured"
        for row in record["rows"]:
            assert row["oti_vs_fd_rel"] < 1e-8, (record["model"], row["symbol"])
    table = payload["summary"]["table"]
    assert table["UMAT_VPDCL"]["DETDG"] == "Exact"
    assert table["UMAT_NKH_1.02"]["DETDG"] == "Exact"
    # NKH's hand-coded ANP1P drops a (1-D) factor: it disagrees with FD, OTI does not
    anp1p = next(r for r in records["UMAT_NKH_1.02"]["rows"] if r["symbol"] == "ANP1P")
    assert anp1p["hand_vs_fd_rel"] > 1e-5
    assert anp1p["oti_vs_fd_rel"] < 1e-8
