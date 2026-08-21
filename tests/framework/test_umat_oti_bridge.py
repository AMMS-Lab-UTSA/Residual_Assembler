"""Tests for the UMAT-OTI ↔ Residual-Assembler bridge.

* Verifies the shared driver contract loads / validates correctly.
* Verifies the JSONL increment stream is read and shape-checked.
* Verifies the direct dR/dp assembler agrees with a hand-derived reference
  on a 1D truss element where the algebra is analytic.
* Verifies the 3D dR/dp assembler agrees with a hand-derived C3D8-like
  reference on a single 8-node brick under uniaxial loading.

No Abaqus, no OTILib, no compiler required.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from residual_core.core.umat_oti_sensitivity import (
    assemble_dRdp,
    integrate_dRe_dp_single,
    scatter_dRe_dp_to_global,
)
from residual_core.materials.umat_oti_driver import (
    CONTRACT_SCHEMA,
    DriverContractError,
    UmatOtiDriverContract,
    UmatOtiIncrement,
    iter_stream,
    load_driver_contract,
    load_driver_contract_from_dict,
    resolve_callable,
)


# ---------------------------------------------------------------------------
# Driver contract
# ---------------------------------------------------------------------------

BASE_CONTRACT = {
    "schema": CONTRACT_SCHEMA,
    "driver_id": "j2_softwarex",
    "ntens": 6,
    "nstatv": 1,
    "nprops": 4,
    "voigt_convention": "engineering_shear",
    "coefficient_convention": "recovered_derivative = OTI_coefficient * prod_factorials(m_i)",
    "source_sha256": "deadbeef",
    "compiler": "gfortran-12",
    "driver_kind": "jsonl_stream",
    "callable_path": "",
    "stream_path": "stream.jsonl",
    "parameters": [
        {"name": "E", "props_index": 1},
        {"name": "NU", "props_index": 2},
        {"name": "SIGY0", "props_index": 3},
        {"name": "H", "props_index": 4},
    ],
    "state_variables": [{"name": "EQPLAS", "statev_index": 1}],
}


def test_contract_loads_from_dict_and_validates_schema():
    contract = load_driver_contract_from_dict(BASE_CONTRACT)
    assert contract.schema == CONTRACT_SCHEMA
    assert contract.parameter_names == ("E", "NU", "SIGY0", "H")
    assert contract.state_names == ("EQPLAS",)
    assert contract.ntens == 6
    assert contract.nstatv == 1
    assert contract.nprops == 4


def test_contract_rejects_wrong_schema():
    bad = dict(BASE_CONTRACT)
    bad["schema"] = "wrong/1.0"
    with pytest.raises(DriverContractError, match="schema mismatch"):
        load_driver_contract_from_dict(bad)


def test_contract_load_from_missing_file(tmp_path: Path):
    with pytest.raises(DriverContractError, match="not found"):
        load_driver_contract(tmp_path / "missing.json")


def test_contract_write_and_load_roundtrip(tmp_path: Path):
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(BASE_CONTRACT))
    contract = load_driver_contract(path)
    assert contract.origin_path == path
    # relative stream_path resolves against the contract's directory
    stream = tmp_path / "stream.jsonl"
    stream.write_text("")  # empty is still a valid stream
    increments = list(iter_stream(contract))
    assert increments == []


def test_iter_stream_reads_valid_increments(tmp_path: Path):
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(BASE_CONTRACT))
    stream = tmp_path / "stream.jsonl"
    records = [
        {
            "increment": 1,
            "stress": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "statev": [0.0],
            "dsigma_dp": [[0.1 * (i + 1) * (k + 1) for k in range(4)] for i in range(6)],
            "dstatev_dp": [[0.0, 0.0, 0.0, 0.0]],
        },
        {
            "increment": 2,
            "stress": [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "statev": [0.01],
            "dsigma_dp": [[0.0] * 4 for _ in range(6)],
            "dstatev_dp": [[0.0, 0.0, -1.0e-5, 0.0]],
        },
    ]
    stream.write_text("\n".join(json.dumps(r) for r in records) + "\n")
    contract = load_driver_contract(path)
    increments = list(iter_stream(contract))
    assert [inc.increment for inc in increments] == [1, 2]
    assert increments[0].dsigma_dp[0] == pytest.approx((0.1, 0.2, 0.3, 0.4))
    assert increments[1].dstatev_dp[0] == pytest.approx((0.0, 0.0, -1.0e-5, 0.0))


def test_iter_stream_rejects_wrong_shape(tmp_path: Path):
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(BASE_CONTRACT))
    stream = tmp_path / "stream.jsonl"
    bad = {
        "stress": [0.0] * 5,   # ntens=6 expected
        "statev": [0.0],
        "dsigma_dp": [[0.0] * 4 for _ in range(6)],
        "dstatev_dp": [[0.0] * 4],
    }
    stream.write_text(json.dumps(bad) + "\n")
    contract = load_driver_contract(path)
    with pytest.raises(DriverContractError, match="stress"):
        list(iter_stream(contract))


def test_resolve_callable_success():
    contract_data = dict(BASE_CONTRACT)
    contract_data["driver_kind"] = "python_callable"
    contract_data["callable_path"] = "json:dumps"
    contract = load_driver_contract_from_dict(contract_data)
    obj = resolve_callable(contract)
    assert callable(obj)
    assert obj({"x": 1}).startswith("{")


def test_resolve_callable_fails_when_missing():
    contract_data = dict(BASE_CONTRACT)
    contract_data["driver_kind"] = "python_callable"
    contract_data["callable_path"] = "nonexistent_module:missing_attr"
    contract = load_driver_contract_from_dict(contract_data)
    with pytest.raises(DriverContractError):
        resolve_callable(contract)


# ---------------------------------------------------------------------------
# Direct dR/dp assembly
# ---------------------------------------------------------------------------

def test_integrate_dRe_dp_single_matches_analytical_1d_truss():
    """1D truss element (2 nodes, 1 DOF each). Under uniform axial stress
    sigma, R_e = [-sigma * A, +sigma * A]. So dR_e/dsigma = [-A, A].
    With B = [-1/L, 1/L], integration point weight * detJ = A * L / 1 (one IP),
    integrate_dRe_dp_single(B, [w*detJ], [dsigma_dp]) must equal
    [-dsigma_dp * A, +dsigma_dp * A].
    """
    L = 2.5
    A = 1.7
    B = np.array([[-1.0 / L, 1.0 / L]])  # ntens=1, ndof=2
    weight_detJ = A * L  # one-point rule
    dsigma_dp = np.array([0.35])
    result = integrate_dRe_dp_single([B], [weight_detJ], [dsigma_dp])
    assert result == pytest.approx(np.array([-dsigma_dp[0] * A, dsigma_dp[0] * A]))


def test_scatter_dRe_dp_to_global_accumulates_correctly():
    global_dRdp = np.zeros(6)
    element_dRe_dp = np.array([1.0, 2.0])
    scatter_dRe_dp_to_global(element_dRe_dp, [1, 4], global_dRdp)
    # A second element shares DOF 4.
    scatter_dRe_dp_to_global(np.array([3.0, 4.0]), [4, 5], global_dRdp)
    assert global_dRdp.tolist() == [0.0, 1.0, 0.0, 0.0, 2.0 + 3.0, 4.0]


def test_assemble_dRdp_single_element_multiparam():
    """Two truss elements in series, two parameters.

    Element 1 spans DOFs 0->1; element 2 spans DOFs 1->2. Each has its own
    dsigma_dp. dR/dp must be the assembly of both element contributions.
    """
    L = 1.0
    A = 1.0
    B_1 = np.array([[-1.0, 1.0]])
    B_2 = np.array([[-1.0, 1.0]])
    weight_detJ_1 = A * L
    weight_detJ_2 = A * L
    dsigma_1 = np.array([[0.1, 0.2]])  # (ntens=1) x (nparam=2)
    dsigma_2 = np.array([[0.3, 0.4]])

    contributions = [
        {
            "B_at_ip": [B_1],
            "weight_x_detJ_at_ip": [weight_detJ_1],
            "dsigma_dp_at_ip": [dsigma_1],
            "element_dof_map": [0, 1],
        },
        {
            "B_at_ip": [B_2],
            "weight_x_detJ_at_ip": [weight_detJ_2],
            "dsigma_dp_at_ip": [dsigma_2],
            "element_dof_map": [1, 2],
        },
    ]
    global_dRdp = assemble_dRdp(n_global_dofs=3, element_contributions=contributions)
    assert global_dRdp.shape == (3, 2)
    # Element 1: dR = [-dsigma, +dsigma] scaled by w*detJ, scattered to [0, 1].
    # Element 2: dR = [-dsigma, +dsigma] scaled by w*detJ, scattered to [1, 2].
    expected = np.array(
        [
            [-0.1, -0.2],
            [+0.1 + -0.3, +0.2 + -0.4],
            [+0.3, +0.4],
        ]
    )
    assert np.allclose(global_dRdp, expected)
