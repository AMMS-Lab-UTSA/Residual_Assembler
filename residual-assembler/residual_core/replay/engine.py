"""Stage-2 orchestration: replay record + OTI package + request -> sensitivities.

Ties the pieces together and REUSES the existing, tested solver
(:func:`residual_core.core.field_sensitivity.solve_field_sensitivities`): the
replay only fills that solver's per-IP field inputs by calling the binary, then
this module assembles the equilibrium residual, solves for du/dp, and propagates
to the requested outputs dq/dp.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Union

import numpy as np

from ..core.dof_manager import DofManager
from ..core.constraints import partition
from ..core.field_sensitivity import solve_field_sensitivities
from .abi import MaterialABI
from .package import MaterialPackage
from .record import ReplayRecord, preflight_record, RecordError
from .driver import (replay_elastic_fields, verify_replay_stress,
                     assemble_internal_force)
from .outputs import OutputContext, compute_outputs


@dataclass
class ReplayResult:
    model_id: str
    parameters: List[str]
    du_dp: Dict[str, np.ndarray]              # {param -> (ndof,)}
    outputs: List[Dict[str, Any]]             # dq/dp per requested response
    equilibrium_free_norm: float              # ||R_free|| at the converged state
    equilibrium_norm: float
    stress_check: Dict[str, float]            # replay vs production stress
    diagnostics: Dict[str, Any] = field(default_factory=dict)
    field_result: Any = None                  # the underlying FieldSensitivityResult

    def du(self, param: str) -> np.ndarray:
        return self.du_dp[param]

    def to_dict(self) -> Dict[str, Any]:
        """Serialise to the resasm_sensitivity_result_v1 contract (du/dp is
        summarised by norm; full arrays stay on the object / .du(param))."""
        return {
            "schema": "resasm_sensitivity_result_v1",
            "model_id": self.model_id,
            "parameters": list(self.parameters),
            "du_dp_norms": {p: float(np.linalg.norm(v)) for p, v in self.du_dp.items()},
            "outputs": self.outputs,
            "equilibrium": {"free_norm": self.equilibrium_free_norm,
                            "norm": self.equilibrium_norm},
            "stress_check": self.stress_check,
            "diagnostics": self.diagnostics,
        }


def replay_sensitivities(
    record: Union[str, ReplayRecord],
    package: Union[str, MaterialPackage],
    parameters: List[str],
    outputs: Optional[List[Mapping[str, Any]]] = None,
    *,
    abi: Optional[MaterialABI] = None,
) -> ReplayResult:
    t0 = time.perf_counter()
    if isinstance(record, str):
        record = ReplayRecord.load(record)
    if isinstance(package, str):
        package = MaterialPackage.load(package)

    # 1. the record must be mathematically sufficient for THIS package
    problems = preflight_record(record, package)
    if problems:
        raise RecordError("replay record is insufficient:\n  - " + "\n  - ".join(problems))

    # 2. load the OTI binary through the ABI; confirm it matches its manifest and
    #    that it is the matched twin of the model that produced the record
    if abi is None:
        abi = MaterialABI(package.binary_path())
    package.check_contract_version()      # shared ABI/manifest contract not drifted
    package.check_binary_hash()           # the .so on disk IS the one the manifest names
    package.check_against_abi(abi)        # the binary's self-reported dims match
    package.check_twin(record.provenance)  # OTI is the matched twin of the record's producer

    for p in parameters:                      # fail early on an unknown parameter
        package.param_index(p)

    # 3. build the model, replay every IP through the binary
    model = record.build_model()
    dm = DofManager(model.nodes.keys())
    t_replay = time.perf_counter()
    rep = replay_elastic_fields(model, record, package, abi, list(parameters), dm)
    t_replay = time.perf_counter() - t_replay
    stress_check = verify_replay_stress(rep["stress_replay"], record)

    # 4. reuse the tested solver: K du/dp = -R_,p
    t_solve = time.perf_counter()
    fs = solve_field_sensitivities(
        model=model,
        tangent_fields=rep["tangent_fields"],
        stress_derivative_fields=rep["stress_derivative_fields"],
        solution=record.converged_u(),
        parameters=list(parameters))
    t_solve = time.perf_counter() - t_solve

    # 5. equilibrium residual R = F_int - F_ext at the converged state
    f_int = assemble_internal_force(model, rep["stress_replay"], dm)
    f_ext = record.external_load(dm)
    R = f_int - f_ext
    free_mask, _pres, _vals = partition(model, dm)
    eqn_free = float(np.linalg.norm(R[free_mask]))
    eqn = float(np.linalg.norm(R))

    du_dp = {p: fs.du_da(p) for p in parameters}

    # 6. propagate to requested outputs dq/dp
    out_results: List[Dict[str, Any]] = []
    if outputs:
        ctx = OutputContext(
            model=model, dof_manager=dm, params=list(parameters),
            u=record.converged_u(), K=fs.K, R_da=fs.residual_derivatives,
            du_dp=fs.displacement_sensitivities, f_ext=f_ext, f_int=f_int,
            tangent_fields=rep["tangent_fields"],
            stress_derivative_fields=rep["stress_derivative_fields"],
            stress_replay=rep["stress_replay"])
        out_results = compute_outputs(list(outputs), ctx)

    diagnostics = {
        "model_id": package.model_id,
        "regular_hash": package.regular_hash(),
        "abi_model_id": abi.model_id,
        "n_elements": len(model.elements),
        "n_dof": dm.ndof,
        "n_parameters": len(parameters),
        "timings_s": {"replay": t_replay, "solve": t_solve,
                      "total": time.perf_counter() - t0},
    }
    return ReplayResult(
        model_id=package.model_id, parameters=list(parameters), du_dp=du_dp,
        outputs=out_results, equilibrium_free_norm=eqn_free, equilibrium_norm=eqn,
        stress_check=stress_check, diagnostics=diagnostics, field_result=fs)


def run_request(
    record: Union[str, ReplayRecord],
    package: Union[str, MaterialPackage],
    request: Mapping[str, Any],
) -> Dict[str, Any]:
    """Run a resasm_sensitivity_request_v1 and return a
    resasm_sensitivity_result_v1 dict. ``record``/``package`` may be paths or
    objects; the request may also carry their paths (request wins if given)."""
    if request.get("schema") not in (None, "resasm_sensitivity_request_v1"):
        raise ValueError("request schema must be resasm_sensitivity_request_v1, got %r"
                         % request.get("schema"))
    if request.get("record"):
        record = request["record"]
    if request.get("material_package"):
        package = request["material_package"]
    params = list(request.get("parameters") or [])
    if not params:
        raise ValueError("request.parameters is required (at least one parameter)")
    outputs = list(request.get("outputs") or [])
    res = replay_sensitivities(record, package, params, outputs)
    return res.to_dict()
