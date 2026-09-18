"""Slides 28-32: residual-method sensitivities of the m5_cpflow C3D8 shear test (claim 3)."""

import json
import shutil

import numpy as np
import pytest

from presentation import claim3_cp_residual_c3d8 as claim3

pytestmark = [pytest.mark.integration, pytest.mark.fortran]


def test_single_element_matches_analytic_chain_rule_mesh_and_fd(out_dirs):
    out, work = out_dirs
    assert claim3.main(["--refine", "2", "--mesh", "2", "--timing-reps", "1", "--kernel-reps", "5",
                        "--out", str(out), "--work", str(work)]) == 0
    record = json.loads((out / "claim3_cp_residual_c3d8.json").read_text())
    single = record["single_element"]
    assert single["sensitivities_in_one_run"] == 6
    assert single["max_du_dp"] < 1e-12
    assert single["enriched_vs_plain_primal_max_rel"] < 1e-13
    analytic = record["analytic_chain_rule"]
    assert analytic["branches"] == ["flow"]
    assert analytic["primal_sigma_vm_max_rel"] < 1e-13
    assert analytic["nrmse_worst"] < 1e-12
    # the analytic reference itself agrees with FD of its own primal
    assert max(analytic["self_check_vs_fd_of_analytic_primal_best_nrmse"].values()) < 1e-7
    # best-case central FD of the ORIGINAL (h chosen with the answer) is still worse than OTI
    for entry in record["central_fd"]["per_parameter"].values():
        assert entry["best_nrmse"] > entry["oti_nrmse"]
        assert entry["best_nrmse"] < 1e-6
    assert record["larger_mesh"]["dvm_mesh_vs_single_max_rel"] < 1e-12
    assert record["larger_mesh"]["max_du_dp"] < 1e-12
    kernel = record["cost"]["compiled_kernel"]["normalised_to_plain"]
    assert kernel["hypad_march"] > 1.0
    assert len(record["series"]["sigma_vm"]) == 50
    assert np.isclose(record["single_element"]["sigma_vm_final_MPa"], record["series"]["sigma_vm"][-1])


@pytest.mark.abaqus
def test_abaqus_single_element_primal(tmp_path):
    if not shutil.which("abaqus"):
        pytest.fail("abaqus-marked test run without Abaqus on PATH")
    from presentation import claim3_abaqus
    from presentation import c3d8_residual as C
    from residual_core.replay.path_material import PathMaterial
    from presentation.common import build_provider
    v2, props, _ = claim3.material()
    build = build_provider(v2, tmp_path / "provider")
    mat = PathMaterial(build["object"], build["contract"], workdir=str(tmp_path / "link"))
    problem = C.simple_shear_problem(1, claim3.GAMMA, 50, 1.0)
    vm = claim3.element_series(problem, C.solve(mat, props, problem, mode="regular"))["vm"]
    names = list(mat.params)
    report = claim3_abaqus.run(props=props, names=names, pidx=[i - 1 for i in mat.props_index],
                               n_increments=50, dt=1.0, gamma=claim3.GAMMA, work=tmp_path / "abaqus",
                               python_vm=vm, oti_dvm=np.zeros((50, 6)), analytic_dvm=np.ones((50, 6)),
                               steps=(), mesh=1)
    assert report["jobs"][0]["completed"]
    assert report["nominal"]["increments"] == 50
    assert report["nominal"]["vs_python_c3d8_max_rel"] < 1e-12
