"""Real compiled-OTI-Fortran-driver -> Residual Assembler test (Priority 7).

This is not a JSONL-only roundtrip test. It:

1. Emits and compiles the SoftwareX PROPS-seeded OTI Fortran J2 driver
   from :mod:`umat_oti.fortran_emit.parameter_sensitivity_j2`, which
   requires a real ``gfortran`` on ``PATH``.
2. Runs the compiled binary to produce ``DSIGMA_DP_OTI.csv`` /
   ``DSTATEV_DP_OTI.csv`` from actual OTI Fortran GETIM extraction.
3. Translates the CSV output into the versioned
   ``umat-oti-driver-contract/1.1`` JSONL increment stream expected by the
   Residual Assembler bridge (:mod:`residual_core.materials.umat_oti_driver`).
4. Loads the stream through the ResAsm bridge and assembles ``dR/dp``
    for a real eight-point C3D8 element using the production kernel.
5. Compares the assembled ``dR/dp`` against a **full-residual finite
   difference** obtained by replaying the full loading history with each
   parameter perturbed positive/negative and re-computing the same
   element's residual from the Python J2 reference.

This test would silently pass on a Python-FD-only bridge (Priority 0
regression) because the OTI CSV would then be a copy of the FD one. It
does not: it insists on going through the compiled Fortran driver
end-to-end.

The test is skipped only when gfortran is missing (environmental blocker).
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest


REQUIRES_GFORTRAN = pytest.mark.skipif(
    shutil.which("gfortran") is None,
    reason="gfortran not on PATH (environmental blocker).",
)


def _require_umat_oti():
    try:
        import umat_oti  # noqa: F401
    except ImportError:
        pytest.skip("umat_oti not importable in this test environment")


def _unit_c3d8_coordinates() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    )


def _element_dof_map() -> list[int]:
    return list(range(24))


def _write_stream_from_oti_csv(
    oti_dsigma_csv: Path,
    oti_dstatev_csv: Path,
    oti_primal_csv: Path,
    stream_path: Path,
) -> None:
    """Convert the compiled driver's CSV output into the driver-contract JSONL."""
    # Read primal first to get stress/statev.
    with oti_primal_csv.open("r", encoding="utf-8") as fh:
        primal = list(csv.DictReader(fh))
    primal_by_inc = {int(row["increment"]): row for row in primal}
    # Group DSIGMA_DP rows by increment.
    with oti_dsigma_csv.open("r", encoding="utf-8") as fh:
        dsigma = list(csv.DictReader(fh))
    dsigma_by_inc: dict[int, list[list[float]]] = {}
    for row in dsigma:
        inc = int(row["increment"])
        vec = [
            float(row["E"]),
            float(row["NU"]),
            float(row["SIGY0"]),
            float(row["H"]),
        ]
        dsigma_by_inc.setdefault(inc, []).append(vec)
    with oti_dstatev_csv.open("r", encoding="utf-8") as fh:
        dstatev = list(csv.DictReader(fh))
    dstatev_by_inc: dict[int, list[list[float]]] = {}
    for row in dstatev:
        inc = int(row["increment"])
        vec = [
            float(row["E"]),
            float(row["NU"]),
            float(row["SIGY0"]),
            float(row["H"]),
        ]
        dstatev_by_inc.setdefault(inc, []).append(vec)
    with stream_path.open("w", encoding="utf-8") as fh:
        for inc in sorted(primal_by_inc):
            record = {
                "increment": inc,
                "stress": [
                    float(primal_by_inc[inc][f"stress_{i}"]) for i in range(1, 7)
                ],
                "statev": [float(primal_by_inc[inc]["EQPLAS"])],
                "dsigma_dp": dsigma_by_inc[inc],
                "dstatev_dp": dstatev_by_inc[inc],
                "source": "oti_fortran_driver",
            }
            fh.write(json.dumps(record))
            fh.write("\n")


def _full_residual_fd(Xe: np.ndarray) -> np.ndarray:
    """Compute dR_e/dp via full loading-history FD of the Python J2 model.

    This uses :mod:`umat_oti.validation.j2_reference` -- the SAME algorithm
    the OTI Fortran driver implements -- so the comparison is between:

    * the assembled dR_e/dp obtained by pushing the *Fortran* driver's
      DSIGMA_DP through the ResAsm bridge, and
    * a finite-difference reference computed from the *Python* reference
      model at the same material point.

    Primal parity between the Fortran driver and the Python reference is
    proven independently in ``tests/test_j2_oti_fortran_driver.py``, so
    any disagreement here indicates a bug in the bridge or the assembler,
    not in the primal stress update.
    """
    from umat_oti.validation.j2_reference import J2Parameters, build_softwarex_j2_path, run_path
    from umat_oti.validation.parameter_sensitivity import ParameterMap
    from residual_core.formulations.c3d8_kernel import element_internal_force_small_strain

    params = J2Parameters()
    path = build_softwarex_j2_path()
    param_map = ParameterMap.softwarex_default()

    baseline = run_path(params, path)
    dp = 1.0e-6
    fd = np.zeros((24, len(param_map.entries)))
    for k, (name, _) in enumerate(param_map.entries):
        # Perturb positive and negative around the operating point.
        p0 = _current(params, name)
        step = dp * abs(p0) if abs(p0) > 1.0 else dp
        params_plus = params.with_replaced(name, p0 + step)
        params_minus = params.with_replaced(name, p0 - step)
        plus = run_path(params_plus, path)
        minus = run_path(params_minus, path)
        stress_plus = np.array(plus[-1].stress)
        stress_minus = np.array(minus[-1].stress)
        residual_plus = element_internal_force_small_strain(Xe, np.tile(stress_plus, (8, 1)))
        residual_minus = element_internal_force_small_strain(Xe, np.tile(stress_minus, (8, 1)))
        fd[:, k] = (residual_plus - residual_minus) / (2.0 * step)
    return fd


def _current(params, name: str) -> float:
    upper = name.upper()
    if upper == "E":
        return params.E
    if upper == "NU":
        return params.nu
    if upper == "SIGY0":
        return params.SIGY0
    if upper == "H":
        return params.H
    raise ValueError(name)


@REQUIRES_GFORTRAN
def test_compiled_oti_j2_drives_bridge_and_dRdp_matches_full_residual_fd(tmp_path: Path):
    _require_umat_oti()
    from umat_oti.fortran_emit.parameter_sensitivity_j2 import (
        compile_j2_oti_build,
        generate_j2_oti_build,
        run_j2_oti_driver,
    )
    from residual_core.core.umat_oti_sensitivity import assemble_dRdp
    from residual_core.materials.umat_oti_driver import (
        iter_stream,
        load_driver_contract,
    )
    from umat_oti.reports.driver_contract import build_softwarex_j2_contract

    # 1) Build and run the compiled Fortran driver.
    build_dir = tmp_path / "oti_build"
    layout = generate_j2_oti_build(build_dir)
    exe = compile_j2_oti_build(layout)
    result = run_j2_oti_driver(exe)
    assert result.returncode == 0, result.stderr

    # 2) Translate the CSVs into the JSONL stream + driver contract.
    stream_path = tmp_path / "j2_stream.jsonl"
    _write_stream_from_oti_csv(
        result.dsigma_csv, result.dstatev_csv, result.primal_csv, stream_path
    )
    contract = build_softwarex_j2_contract(stream_path=stream_path.name)
    contract_path = tmp_path / "driver_contract.json"
    contract.write(contract_path)

    # 3) Load through the ResAsm bridge; keep the LAST increment (fully
    # plastic, history-dependent) as the design point.
    loaded = load_driver_contract(contract_path)
    increments = list(iter_stream(loaded))
    assert increments, "empty stream"
    last = increments[-1]
    dsigma_dp = np.array(last.dsigma_dp)   # shape (6, 4)
    assert dsigma_dp.shape == (6, 4)

    # 4) Assemble dR_e/dp through the production eight-point C3D8 kernel.
    from residual_core.formulations.c3d8_kernel import ABAQUS_C3D8_GAUSS, b_matrix_reference

    Xe = _unit_c3d8_coordinates()
    B_at_ip = []
    weight_detJ_at_ip = []
    for point, weight in zip(ABAQUS_C3D8_GAUSS.points, ABAQUS_C3D8_GAUSS.weights):
        B, detJ = b_matrix_reference(Xe, point)
        B_at_ip.append(B)
        weight_detJ_at_ip.append(weight * detJ)
    contributions = [
        {
            "B_at_ip": B_at_ip,
            "weight_x_detJ_at_ip": weight_detJ_at_ip,
            "dsigma_dp_at_ip": [dsigma_dp] * 8,
            "element_dof_map": _element_dof_map(),
        }
    ]
    assembled = assemble_dRdp(n_global_dofs=24, element_contributions=contributions)
    assert assembled.shape == (24, 4)

    # 5) Full-residual FD reference computed from the Python J2 model at
    # the same last increment.
    fd_reference = _full_residual_fd(Xe)
    assert fd_reference.shape == (24, 4)

    # 6) Compare. The scale for each column varies wildly (E and NU rows
    # differ by ~1e5), so use a scale-aware relative tolerance per entry.
    max_rel = 0.0
    for dof in range(24):
        for k in range(4):
            a = assembled[dof, k]
            b = fd_reference[dof, k]
            scale = max(abs(a), abs(b), 1.0)
            rel = abs(a - b) / scale
            if rel > max_rel:
                max_rel = rel
    assert max_rel < 1.0e-4, max_rel
