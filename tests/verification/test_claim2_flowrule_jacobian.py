"""Slides 26-27: OTI vs the hand-coded analytical flow-rule Jacobian, FD, timing (claim 2).

The flow-rule sources are not part of this repository (see the claim script);
the test runs wherever they are available and says so where they are not.
"""

import json
import shutil

import pytest

from verification import claim2_flowrule_jacobian as claim2

SOURCE = claim2.default_source_dir()
pytestmark = [
    pytest.mark.integration, pytest.mark.fortran,
    pytest.mark.skipif(not (SOURCE / "computeFlowRule.f90").is_file(),
                       reason=f"external flow-rule sources not found in {SOURCE}; set CLAIM2_FLOWRULE_DIR"),
]


def test_oti_matches_hand_coded_jacobian_and_fd_is_orders_worse(out_dirs):
    out, work = out_dirs
    compiler = "ifort" if shutil.which("ifort") else "gfortran"
    assert claim2.main(["--calls", "2000", "--repetitions", "2", "--compiler", compiler,
                        "--out", str(out), "--work", str(work)]) == 0
    payload = json.loads((out / "claim2_flowrule_jacobian.json").read_text())
    result = payload["results"][compiler]
    summary = result["accuracy"]["summary"]
    assert summary["active_slip_systems"] == 12
    assert summary["oti_vs_analytical_max_rel"] < 1e-13
    assert summary["primal_oti_vs_analytical_max_rel"] < 1e-13
    # even at its best step, centred FD is several orders of magnitude worse
    assert summary["fd_centred_best_step_max_rel_smooth_variables"] > 1e4 * summary["oti_vs_analytical_max_rel"]
    assert summary["fd_forward_best_step_max_rel_all_variables"] > 1e6 * summary["oti_vs_analytical_max_rel"]
    # each slip system depends only on its own variables (FD Jacobian may perturb all at once)
    for entry in result["accuracy"]["per_variable"].values():
        assert entry["all_systems_vs_one_system_fd_max_rel"] < 1e-6
    ratios = result["timing"]["ratios"]
    assert ratios["fd_centred_known_step_over_oti"] > 1.0
    assert ratios["fd_centred_step_search_over_oti"] > ratios["fd_centred_known_step_over_oti"]
