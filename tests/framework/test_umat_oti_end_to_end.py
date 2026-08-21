"""End-to-end integration test: UMAT-OTI J2 driver → JSONL stream → ResAsm.

This test exercises the complete SoftwareX bridge on the reference J2 case:

1. Run the UMAT-OTI J2 parameter-sensitivity driver
   (``umat_oti.validation.parameter_sensitivity``) to obtain
   ``DSIGMA_DP`` and ``DSTATEV_DP`` for the SoftwareX loading path.
2. Emit an ``umat-oti-driver-contract/1.1`` JSON contract plus a JSONL
   increment stream via ``umat_oti.reports.driver_contract``.
3. Load both on the Residual Assembler side through
   ``residual_core.materials.umat_oti_driver``.
4. Assemble ``dR/dp`` on a single 1D truss "element" whose stress equals the
   J2 axial stress at each increment. Verify that the two sides are
   arithmetically consistent (i.e. the assembled ``dR/dp`` equals a
   hand-derived analytical form built from the same ``DSIGMA_DP``).

The test does not need Abaqus, OTILib, or a Fortran compiler.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("umat_oti")

from umat_oti.reports.driver_contract import (
    build_softwarex_j2_contract,
    write_j2_stream,
)
from umat_oti.validation.j2_reference import J2Parameters, build_softwarex_j2_path
from umat_oti.validation.parameter_sensitivity import (
    ParameterMap,
    StateMap,
    compute_j2_parameter_sensitivities,
)

from residual_core.core.umat_oti_sensitivity import assemble_dRdp
from residual_core.materials.umat_oti_driver import (
    iter_stream,
    load_driver_contract,
)


PARAMS = J2Parameters()
PATH = build_softwarex_j2_path()


def _emit_bridge_artifacts(tmp_path: Path) -> tuple[Path, list[dict]]:
    """Produce the driver contract + JSONL stream from UMAT-OTI."""
    run = compute_j2_parameter_sensitivities(
        params=PARAMS,
        path=PATH,
        parameter_map=ParameterMap.softwarex_default(),
        state_map=StateMap.softwarex_default(),
        fd_step_relative=1.0e-6,
    )
    records = []
    for inc in run.increments:
        records.append(
            {
                "increment": inc.increment,
                "stress": list(inc.stress),
                "statev": list(inc.statev),
                "dsigma_dp": [list(row) for row in inc.dsigma_dp],
                "dstatev_dp": [list(row) for row in inc.dstatev_dp],
            }
        )
    stream_path = tmp_path / "j2_stream.jsonl"
    write_j2_stream(records, stream_path)
    contract = build_softwarex_j2_contract(stream_path=stream_path.name)
    contract_path = tmp_path / "j2_contract.json"
    contract.write(contract_path)
    return contract_path, records


def test_bridge_end_to_end_dRdp_matches_analytical_truss(tmp_path: Path):
    contract_path, expected_records = _emit_bridge_artifacts(tmp_path)
    contract = load_driver_contract(contract_path)
    increments = list(iter_stream(contract))
    assert [inc.increment for inc in increments] == [
        rec["increment"] for rec in expected_records
    ], "stream ordering must match the driver output"

    # Bridge to a 1D truss element with a single IP:
    #   B = [-1, 1] / L, weight * detJ = A * L (one-point rule).
    #   dR_e/dp = A * L * B^T * (dsigma11/dp)  =  A * [-dsigma11/dp, dsigma11/dp]
    L = 2.0
    A = 3.0
    B = np.array([[-1.0 / L, 1.0 / L]])   # ntens=1, ndof=2
    weight_detJ = A * L

    # Focus on the last increment (post-yield: history-dependent).
    inc = increments[-1]
    dsigma_axial_row = np.array(inc.dsigma_dp[0])  # first stress component
    element_contributions = [
        {
            "B_at_ip": [B],
            "weight_x_detJ_at_ip": [weight_detJ],
            # For a 1-D bridge we only use the axial (11) component.
            "dsigma_dp_at_ip": [dsigma_axial_row.reshape(1, -1)],
            "element_dof_map": [0, 1],
        }
    ]
    dRdp = assemble_dRdp(n_global_dofs=2, element_contributions=element_contributions)
    expected = np.stack([-A * dsigma_axial_row, A * dsigma_axial_row])
    assert dRdp.shape == (2, len(dsigma_axial_row))
    assert np.allclose(dRdp, expected)


def test_bridge_end_to_end_stream_shapes_are_valid(tmp_path: Path):
    """Every record in the stream must satisfy the contract's shape checks."""
    contract_path, _ = _emit_bridge_artifacts(tmp_path)
    contract = load_driver_contract(contract_path)
    increments = list(iter_stream(contract))
    for inc in increments:
        assert len(inc.stress) == contract.ntens
        assert len(inc.statev) == contract.nstatv
        assert len(inc.dsigma_dp) == contract.ntens
        for row in inc.dsigma_dp:
            assert len(row) == len(contract.parameters)
        assert len(inc.dstatev_dp) == contract.nstatv
        for row in inc.dstatev_dp:
            assert len(row) == len(contract.parameters)
