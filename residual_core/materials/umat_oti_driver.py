"""UMAT-OTI driver bridge for the Residual Assembler.

Consumes an ``umat-oti-driver-contract/1.1`` payload emitted by the
UMAT-OTI project (see ``umat_oti.reports.driver_contract``) and exposes it
to the assembler as a :class:`~residual_core.materials.base.Material` that
carries the *point-wise parameter sensitivities* ``DSIGMA_DP`` / ``DSTATEV_DP``
alongside the usual (STRESS, DDSDDE, STATEV) triple.

Two ingestion modes are supported by the shared contract:

* ``driver_kind == "jsonl_stream"`` — the material driver is a pre-computed
  stream of increments (one JSON object per line) with the fields
  ``stress``, ``statev``, ``dsigma_dp``, ``dstatev_dp``, and optionally
  ``ddsdde``. This is the mode used by the SoftwareX offline verification
  path: the UMAT-OTI J2 reference driver emits the stream, ResAsm replays
  it, and dR/dp is assembled from ``B^T · DSIGMA_DP`` at each integration
  point.

* ``driver_kind == "python_callable"`` — the material driver is an
  importable Python object at ``callable_path``. Not exercised in this
  offline suite (it needs a compiled OTI-seeded UMAT); the loader validates
  the contract and hands the dotted path back for the caller to resolve.

The bridge is intentionally minimal: it does not decode Fortran ABIs and it
does not try to drive a compiled UMAT-OTI ``.so``. Those live outside this
module.
"""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple


CONTRACT_SCHEMA = "umat-oti-driver-contract/1.1"


class DriverContractError(RuntimeError):
    """Raised when a driver contract is malformed or the schema mismatches."""


@dataclass(frozen=True)
class ParameterEntry:
    name: str
    props_index: int


@dataclass(frozen=True)
class StateEntry:
    name: str
    statev_index: int


@dataclass
class UmatOtiDriverContract:
    """In-memory view of an ``umat-oti-driver-contract/1.1`` payload."""

    schema: str
    driver_id: str
    ntens: int
    nstatv: int
    nprops: int
    parameters: Tuple[ParameterEntry, ...]
    state_variables: Tuple[StateEntry, ...]
    voigt_convention: str
    compiler: str
    coefficient_convention: str
    source_sha256: str
    driver_kind: str
    callable_path: str
    stream_path: str
    origin_path: Optional[Path] = field(default=None)

    @property
    def parameter_names(self) -> Tuple[str, ...]:
        return tuple(p.name for p in self.parameters)

    @property
    def state_names(self) -> Tuple[str, ...]:
        return tuple(s.name for s in self.state_variables)


def load_driver_contract(path: str | Path) -> UmatOtiDriverContract:
    """Load and validate a driver-contract JSON file."""
    p = Path(path)
    if not p.is_file():
        raise DriverContractError(f"contract file not found: {p}")
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise DriverContractError(f"invalid JSON in {p}: {exc}") from exc
    return _build_contract(data, origin_path=p)


def load_driver_contract_from_dict(data: Dict[str, Any]) -> UmatOtiDriverContract:
    """Load a contract from an in-memory dict (used in tests)."""
    return _build_contract(data, origin_path=None)


def _build_contract(data: Dict[str, Any], *, origin_path: Optional[Path]) -> UmatOtiDriverContract:
    schema = data.get("schema")
    if schema != CONTRACT_SCHEMA:
        raise DriverContractError(
            f"driver contract schema mismatch: expected {CONTRACT_SCHEMA!r}, "
            f"got {schema!r}"
        )
    try:
        return UmatOtiDriverContract(
            schema=schema,
            driver_id=str(data.get("driver_id", "")),
            ntens=int(data["ntens"]),
            nstatv=int(data["nstatv"]),
            nprops=int(data["nprops"]),
            parameters=tuple(
                ParameterEntry(name=str(p["name"]), props_index=int(p["props_index"]))
                for p in data.get("parameters", [])
            ),
            state_variables=tuple(
                StateEntry(name=str(s["name"]), statev_index=int(s["statev_index"]))
                for s in data.get("state_variables", [])
            ),
            voigt_convention=str(data.get("voigt_convention", "engineering_shear")),
            compiler=str(data.get("compiler", "")),
            coefficient_convention=str(data.get("coefficient_convention", "")),
            source_sha256=str(data.get("source_sha256", "")),
            driver_kind=str(data.get("driver_kind", "jsonl_stream")),
            callable_path=str(data.get("callable_path", "")),
            stream_path=str(data.get("stream_path", "")),
            origin_path=origin_path,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise DriverContractError(f"malformed contract payload: {exc}") from exc


# ---------------------------------------------------------------------------
# Increment stream
# ---------------------------------------------------------------------------

@dataclass
class UmatOtiIncrement:
    """One increment as delivered by the driver stream."""

    increment: Optional[int]
    stress: Tuple[float, ...]           # shape (ntens,)
    statev: Tuple[float, ...]           # shape (nstatv,)
    dsigma_dp: Tuple[Tuple[float, ...], ...]     # shape (ntens, nparam)
    dstatev_dp: Tuple[Tuple[float, ...], ...]    # shape (nstatv, nparam)
    ddsdde: Optional[Tuple[Tuple[float, ...], ...]] = None


def iter_stream(contract: UmatOtiDriverContract) -> Iterator[UmatOtiIncrement]:
    """Yield validated increments from the contract's JSONL stream."""
    if contract.driver_kind != "jsonl_stream":
        raise DriverContractError(
            f"iter_stream() requires driver_kind='jsonl_stream'; got "
            f"{contract.driver_kind!r}"
        )
    stream_path = _resolve_stream_path(contract)
    if not stream_path.is_file():
        raise DriverContractError(f"stream file not found: {stream_path}")
    nparam = len(contract.parameters)
    with stream_path.open("r", encoding="utf-8") as fh:
        for line_index, raw in enumerate(fh, start=1):
            raw = raw.strip()
            if not raw:
                continue
            try:
                data = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise DriverContractError(
                    f"{stream_path}:{line_index} invalid JSON: {exc}"
                ) from exc
            yield _build_increment(data, contract, nparam, line_index, stream_path)


def _resolve_stream_path(contract: UmatOtiDriverContract) -> Path:
    p = Path(contract.stream_path)
    if p.is_absolute() or contract.origin_path is None:
        return p
    return (contract.origin_path.parent / p).resolve()


def _build_increment(
    data: Dict[str, Any],
    contract: UmatOtiDriverContract,
    nparam: int,
    line_index: int,
    stream_path: Path,
) -> UmatOtiIncrement:
    def _expect_len(label: str, vector: Any, expected: int) -> Tuple[float, ...]:
        if not isinstance(vector, list):
            raise DriverContractError(f"{stream_path}:{line_index} {label} must be a list")
        if len(vector) != expected:
            raise DriverContractError(
                f"{stream_path}:{line_index} {label} has length {len(vector)}, "
                f"expected {expected}"
            )
        return tuple(float(v) for v in vector)

    def _expect_matrix(label: str, mat: Any, rows: int, cols: int) -> Tuple[Tuple[float, ...], ...]:
        if not isinstance(mat, list) or len(mat) != rows:
            raise DriverContractError(
                f"{stream_path}:{line_index} {label} must be a list of length {rows}"
            )
        result: List[Tuple[float, ...]] = []
        for row in mat:
            if not isinstance(row, list) or len(row) != cols:
                raise DriverContractError(
                    f"{stream_path}:{line_index} {label} row must be a list of length {cols}"
                )
            result.append(tuple(float(v) for v in row))
        return tuple(result)

    stress = _expect_len("stress", data.get("stress"), contract.ntens)
    statev = _expect_len("statev", data.get("statev"), contract.nstatv)
    dsigma = _expect_matrix("dsigma_dp", data.get("dsigma_dp"), contract.ntens, nparam)
    dstatev = _expect_matrix("dstatev_dp", data.get("dstatev_dp"), contract.nstatv, nparam)
    ddsdde: Optional[Tuple[Tuple[float, ...], ...]] = None
    if "ddsdde" in data:
        ddsdde = _expect_matrix("ddsdde", data["ddsdde"], contract.ntens, contract.ntens)
    increment = data.get("increment")
    return UmatOtiIncrement(
        increment=int(increment) if isinstance(increment, (int, float)) else None,
        stress=stress,
        statev=statev,
        dsigma_dp=dsigma,
        dstatev_dp=dstatev,
        ddsdde=ddsdde,
    )


# ---------------------------------------------------------------------------
# Callable-based drivers (used when the OTI-seeded UMAT is available)
# ---------------------------------------------------------------------------

def resolve_callable(contract: UmatOtiDriverContract):
    """Resolve ``contract.callable_path`` to a Python callable.

    ``callable_path`` uses the standard ``package.module:attribute`` notation.
    """
    if contract.driver_kind != "python_callable":
        raise DriverContractError(
            f"resolve_callable() requires driver_kind='python_callable'; got "
            f"{contract.driver_kind!r}"
        )
    if not contract.callable_path:
        raise DriverContractError("contract.callable_path is empty")
    if ":" in contract.callable_path:
        module_name, attr_name = contract.callable_path.split(":", 1)
    elif "." in contract.callable_path:
        module_name, attr_name = contract.callable_path.rsplit(".", 1)
    else:
        raise DriverContractError(
            f"callable_path must be 'module.path:attr' or 'module.path.attr'; "
            f"got {contract.callable_path!r}"
        )
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise DriverContractError(
            f"could not import driver module {module_name!r}: {exc}"
        ) from exc
    try:
        obj = getattr(module, attr_name)
    except AttributeError as exc:
        raise DriverContractError(
            f"module {module_name!r} has no attribute {attr_name!r}"
        ) from exc
    if not callable(obj):
        raise DriverContractError(
            f"callable_path {contract.callable_path!r} did not resolve to a callable"
        )
    return obj


__all__ = [
    "CONTRACT_SCHEMA",
    "DriverContractError",
    "ParameterEntry",
    "StateEntry",
    "UmatOtiDriverContract",
    "UmatOtiIncrement",
    "iter_stream",
    "load_driver_contract",
    "load_driver_contract_from_dict",
    "resolve_callable",
]
