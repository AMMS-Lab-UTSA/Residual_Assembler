"""UMAT transformation backend adapter (M5-B, slice 1).

Residual Assembler does NOT transform Fortran itself. It calls the separate
``UMAT_source_transformation`` project (Python package ``umat_oti``) through a
small adapter interface, and consumes a single authoritative *handoff manifest*
(``resasm_umat_transform_v1``) from which everything downstream is derived
(*Depvar, SDV layout, patched .inp, recipe, compiler env). The repositories stay
separate; only this adapter knows the transformer's API.

This slice wires the transformer's existing, verified capability -- transform an
ordinary UMAT into an OTI UMAT that produces the consistent tangent ``DDSDDE`` by
automatic differentiation. Parameter seeding (writing ``dsigma/da_i`` per selected
PROPS to SDVs) is a codegen extension in the transformer repo and is the next
slice; the handoff manifest already carries the ``parameters`` shape so the
consumer does not change when it lands.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

HANDOFF_SCHEMA = "resasm_umat_transform_v1"
_DDSDDE_SLOTS = 36
_PER_PARAM = 6


class UmatBackendError(Exception):
    """A user-facing problem locating or invoking the UMAT transformer."""


# --------------------------------------------------------------------------- #
# data model (transformer-agnostic)
# --------------------------------------------------------------------------- #
@dataclass
class ParameterSelection:
    name: str
    index: int                        # 1-based PROPS index
    source: str = "PROPS"


@dataclass
class UmatInspection:
    props_refs: List[int] = field(default_factory=list)
    statev_max: int = 0
    writes_ddsdde: bool = False
    already_oti: bool = False
    subroutines: List[str] = field(default_factory=list)


@dataclass
class TransformationProposal:
    config_path: str
    ntens: int
    config: Dict[str, Any] = field(default_factory=dict)
    parameters: List[ParameterSelection] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    original_umat: Optional[str] = None            # absolute path to the ordinary UMAT
    props_values: List[float] = field(default_factory=list)
    tangent_variable: str = ""


@dataclass
class TransformationResult:
    success: bool
    transformed_umat: Optional[str] = None      # the *combined* Abaqus-ready source
    transformed_plain: Optional[str] = None      # the UMAT-only transformed source
    generated_files: List[str] = field(default_factory=list)
    blockers: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    semantic_checks: Dict[str, Any] = field(default_factory=dict)
    report: Dict[str, Any] = field(default_factory=dict)
    manifest: Optional[Dict[str, Any]] = None
    # context needed for the numerical real-response check (set in transform())
    output_dir: Optional[str] = None
    original_umat: Optional[str] = None
    ntens: int = 0
    parameters: List[ParameterSelection] = field(default_factory=list)
    props_values: List[float] = field(default_factory=list)


def _row_major(i: int, j: int, ntens: int) -> int:
    """1-based row-major flat index of tangent entry (i,j) in a ntens x ntens
    block: STATEV(ntens*(i-1)+j), matching the M4 DDSDDE SDV convention."""
    return ntens * (i - 1) + j


def build_parameter_contracts(ntens: int, params: List[ParameterSelection],
                              tangent_var: str, stress_line: int
                              ) -> List[Dict[str, Any]]:
    """Turn (parameter, PROPS index) selections into umat_oti extra-Jacobian
    contracts that (a) seed each PROPS(index) and extract dsigma/da_i into its
    STATEV block, and (b) copy the real tangent into STATEV[1..36] (row-major),
    reproducing the hand-written M4 SDV layout automatically.

    DDSDDE -> SDV uses a ``real_copy_map`` from the promoted tangent variable, so
    it needs no OTI direction (unlike the GETIM parameter extractions)."""
    ddsdde_components = [
        {"target_indices": [_row_major(i, j, ntens)], "output_indices": [i, j],
         "seed_direction_offset": 0}
        for i in range(1, ntens + 1) for j in range(1, ntens + 1)
    ] if tangent_var else []
    contracts: List[Dict[str, Any]] = []
    for k, p in enumerate(params):
        lo = _DDSDDE_SLOTS + 1 + k * _PER_PARAM
        param_components = [
            {"target_indices": [lo + i - 1], "output_indices": [i],
             "seed_direction_offset": 0}
            for i in range(1, ntens + 1)
        ]
        additional = [{
            "target_variable": "STATEV", "from_output_variable": "STRESS",
            "after_line": stress_line, "extract_kind": "component_map",
            "components": param_components,
        }]
        if k == 0 and tangent_var:
            additional.append({
                "target_variable": "STATEV", "from_output_variable": tangent_var,
                "after_line": stress_line, "extract_kind": "real_copy_map",
                "components": ddsdde_components,
            })
        contracts.append({
            "id": "dsigma_d%s" % p.name,
            "seed": {"variable": p.source, "shape": "vector", "directions": 1,
                     "components": [[int(p.index)]],
                     "operating_point_expression": "%s(%d)" % (p.source, p.index)},
            "output": {"variable": "STRESS", "shape": "vector"},
            "internal_use": {"replace_variable": p.source},
            "additional_extractions": additional,
        })
    return contracts


@dataclass
class TransformationValidation:
    passed: bool
    checks: Dict[str, Any] = field(default_factory=dict)
    message: str = ""


# --------------------------------------------------------------------------- #
# interface
# --------------------------------------------------------------------------- #
class UmatTransformationBackend:
    """The adapter Residual Assembler codes against; ``UmatOtiBackend`` is the
    concrete implementation. A future backend (or a different transformer) only
    has to satisfy this interface."""

    def inspect(self, umat_path: str) -> UmatInspection:
        raise NotImplementedError

    def propose_contract(self, umat_path: str,
                         parameters: List[ParameterSelection],
                         base_config: Optional[Dict[str, Any]] = None,
                         ) -> TransformationProposal:
        raise NotImplementedError

    def transform(self, proposal: TransformationProposal,
                  output_dir: str) -> TransformationResult:
        raise NotImplementedError

    def validate_transformation_semantics(self, result: TransformationResult
                                          ) -> TransformationValidation:
        """Structural check: the generated source has the expected shape (real
        stress extracted before/independent of the derivative extraction, etc.).
        This does NOT prove numerically that the transformed UMAT preserves the
        original response -- for that see :meth:`validate_real_response`."""
        raise NotImplementedError

    def validate_real_response(self, result: TransformationResult
                               ) -> TransformationValidation:
        """NUMERICAL original-vs-transformed comparison (material-point or Abaqus):
        sigma_original ~ sigma_OTI.R and u_original ~ u_OTI. Not the semantic
        check."""
        raise NotImplementedError


# --------------------------------------------------------------------------- #
# umat_oti implementation
# --------------------------------------------------------------------------- #
def _default_umat_oti_src() -> Optional[str]:
    cand = os.environ.get("UMAT_OTI_SRC") or os.path.expanduser(
        "~/Documents/UMAT_source_transformation/src")
    return cand if os.path.isdir(os.path.join(cand, "umat_oti")) else None


class UmatOtiBackend(UmatTransformationBackend):
    def __init__(self, src_dir: Optional[str] = None):
        self.src_dir = src_dir or _default_umat_oti_src()
        if not self.src_dir:
            raise UmatBackendError(
                "the umat_oti transformer package was not found. Set UMAT_OTI_SRC "
                "to the 'src' folder of UMAT_source_transformation (a checkout of "
                "github santiagarcia/UMAT_source_transformation).")

    # -- lazy import of the transformer's Python API --------------------- #
    def _cli(self):
        if self.src_dir not in sys.path:
            sys.path.insert(0, self.src_dir)
        try:
            from umat_oti.cli_json import main as cli_main
            from umat_oti.core.config_loader import load_project_config_json
        except Exception as exc:                # noqa: BLE001
            raise UmatBackendError("cannot import umat_oti from %s: %s"
                                   % (self.src_dir, exc))
        return cli_main, load_project_config_json

    # -- interface -------------------------------------------------------- #
    def inspect(self, umat_path: str) -> UmatInspection:
        import re
        with open(umat_path, "r", encoding="utf-8", errors="replace") as fh:
            src = fh.read()
        up = src.upper()
        return UmatInspection(
            props_refs=sorted(set(int(m) for m in re.findall(r"PROPS\(\s*(\d+)", up))),
            statev_max=max([int(m) for m in re.findall(r"STATEV\(\s*(\d+)", up)] or [0]),
            writes_ddsdde="DDSDDE" in up,
            already_oti=("USE OTIM" in up or "ONUMM" in up),
            subroutines=re.findall(r"^\s*SUBROUTINE\s+(\w+)", src, re.I | re.M))

    def propose_contract(self, umat_path: str,
                         parameters: List[ParameterSelection],
                         base_config: Optional[str] = None,
                         ) -> TransformationProposal:
        if base_config is None:
            raise UmatBackendError(
                "propose_contract needs a base umat_oti config for this UMAT "
                "(auto-config generation from the scanner is a later slice). "
                "Provide the transformer's contract JSON path for the material.")
        config_path = str(base_config)
        # Read the RAW user-format config (NOT the elaborated form): we augment it
        # with parameter-seeding contracts and hand the merged config back to the
        # transformer's own loader, which elaborates it exactly once.
        with open(config_path, "r", encoding="utf-8") as fh:
            base = json.load(fh)
        ntens = int(base.get("ntens")
                    or (base.get("transformation_settings", {}) or {}).get("ntens")
                    or 6)
        resasm = dict(base.get("resasm") or {})
        tangent_var = str(resasm.get("tangent_variable") or "").upper()
        stress_line = int(resasm.get("stress_update_line") or 0)
        props_values = [float(v) for v in (resasm.get("props_values") or [])]
        params = list(parameters)

        config = dict(base)
        config.pop("resasm", None)                 # adapter metadata, not a transformer key
        umat_abs = os.path.abspath(umat_path)
        config["source"] = {"file": umat_abs}      # absolute -> resolves regardless of cwd
        notes = []
        if params:
            # The SDV layout constants (_DDSDDE_SLOTS=36, _PER_PARAM=6) and the whole
            # C3D8/selective-reduced pipeline assume ntens=6; refuse anything else
            # rather than emit a manifest whose advertised ranges do not match the
            # ntens-parameterised writes.
            if ntens != 6:
                raise UmatBackendError(
                    "parameter seeding currently supports ntens=6 (C3D8) only; got "
                    "ntens=%d. The SDV layout (DDSDDE[1-36] + 6/param) is fixed to it."
                    % ntens)
            if not stress_line:
                raise UmatBackendError(
                    "base config needs resasm.stress_update_line (the source line of "
                    "the dependent-variable update) so the derivative extractions can "
                    "be placed after the real write-back")
            # DDSDDE must be packed into STATEV[1..36] for the residual method to
            # reassemble the global tangent K from the ODB; that copy is emitted only
            # when the tangent variable is known.
            if not tangent_var:
                raise UmatBackendError(
                    "base config needs resasm.tangent_variable (the UMAT's tangent "
                    "matrix, e.g. DDS) so the real DDSDDE can be packed into "
                    "STATEV[1..36]; parameter sensitivities alone cannot assemble K")
            # A requested PROPS index the UMAT never reads would extract an identically
            # zero (and un-catchable) derivative; reject it up front with the actual
            # PROPS references found in the source.
            refs = set(self.inspect(umat_abs).props_refs)
            bad = [p for p in params
                   if str(p.source).upper() == "PROPS" and refs and p.index not in refs]
            if bad:
                raise UmatBackendError(
                    "requested PROPS index %s not read by the UMAT (it references "
                    "PROPS%s); a derivative there would be a silent zero"
                    % (", ".join(str(p.index) for p in bad), sorted(refs)))
            if props_values:
                over = [p for p in params if p.index > len(props_values)]
                if over:
                    raise UmatBackendError(
                        "requested PROPS index %s exceeds resasm.props_values length %d"
                        % (", ".join(str(p.index) for p in over), len(props_values)))
            config["extra_jacobian_contracts"] = build_parameter_contracts(
                ntens, params, tangent_var, stress_line)
            notes.append(
                "parameter derivatives d(sigma)/d{%s} -> STATEV[%d..], "
                "real DDSDDE -> STATEV[1..%d]"
                % (", ".join(p.name for p in params), _DDSDDE_SLOTS + 1, _DDSDDE_SLOTS))
        return TransformationProposal(
            config_path=config_path, ntens=ntens, config=config,
            parameters=params, notes=notes, original_umat=umat_abs,
            props_values=props_values, tangent_variable=tangent_var)

    def transform(self, proposal: TransformationProposal,
                  output_dir: str) -> TransformationResult:
        cli_main, _ = self._cli()
        os.makedirs(output_dir, exist_ok=True)
        # Write the augmented (parameter-seeded) config and drive the transformer
        # through its own Python entry (load + merge anchors + transform + combined
        # source + reports), then read the structured report.
        contract_path = os.path.join(output_dir, "contract.json")
        with open(contract_path, "w", encoding="utf-8") as fh:
            json.dump(proposal.config, fh, indent=2)
        import io
        import contextlib
        with contextlib.redirect_stdout(io.StringIO()):
            cli_main(["--config", contract_path, "--out", output_dir])
        rep_path = os.path.join(output_dir, "transform_report.json")
        if not os.path.exists(rep_path):
            raise UmatBackendError("transformer produced no report in %s" % output_dir)
        with open(rep_path, "r", encoding="utf-8") as fh:
            report = json.load(fh)
        gen = [os.path.join(output_dir, f) for f in os.listdir(output_dir)]
        combined = _first_match(gen, "_combined")
        plain = report.get("transformed_source") or _first_match(gen, "_oti.f")
        src_file = (report.get("source_file") or report.get("source")
                    or _config_source_file(proposal.config) or "umat")
        result = TransformationResult(
            success=bool(report.get("success", report.get("transform_success"))),
            transformed_umat=str(combined) if combined else None,
            transformed_plain=str(plain) if plain else None,
            generated_files=sorted(gen),
            blockers=list(report.get("blockers", []) or []),
            warnings=list(report.get("warnings", []) or []),
            semantic_checks=dict(report.get("semantic_checks", {})),
            report=report,
            output_dir=output_dir,
            original_umat=proposal.original_umat,
            ntens=proposal.ntens,
            parameters=list(proposal.parameters),
            props_values=list(proposal.props_values))
        result.manifest = self._manifest(proposal, result, src_file)
        return result

    def validate_transformation_semantics(self, result: TransformationResult
                                          ) -> TransformationValidation:
        """Structural (not numerical): the transformer's semantic checks confirm
        the generated source has the expected shape -- real stress extracted
        before/independent of the DDSDDE extraction, seed before stress update,
        old assignments disabled, etc. This does NOT prove the transformed UMAT
        reproduces the original response; that is validate_real_response (a real
        material-point / Abaqus comparison), reserved for the parameter-seeding
        slice where derivative numbers exist to compare."""
        checks = dict(result.semantic_checks)
        must = ("ddsdde_output_present",
                "real_stress_extraction_before_ddsdde_extraction",
                "stress_oti_update_before_real_stress_extraction")
        passed = bool(result.success) and all(checks.get(k, False) for k in must)
        return TransformationValidation(
            passed=passed, checks=checks,
            message="transformation semantics%s consistent" % ("" if passed else " NOT"))

    def validate_real_response(self, result: TransformationResult
                               ) -> TransformationValidation:
        """NUMERICAL original-vs-transformed comparison at a material point: compile
        BOTH UMATs with a shared standalone driver and check each SDV-extracted
        d(sigma)/da_i against a central finite difference of the ORIGINAL UMAT's
        real stress response (real gfortran, real numbers). Returns passed=False
        with the detail dict when the numbers disagree; a NotImplementedError only
        when there is nothing to compare (no generated parameter derivatives)."""
        gen = (result.manifest or {}).get("derivatives", {}).get("parameters", {})
        if not gen:
            raise NotImplementedError(
                "no generated parameter derivatives to validate; run a "
                "parameter-seeded transform first (validate_transformation_semantics "
                "covers the structural check)")
        if not (result.original_umat and result.transformed_plain and result.output_dir):
            raise UmatBackendError(
                "real-response validation needs original_umat + transformed_plain + "
                "output_dir on the result (set by transform())")
        from resasm_user.umat_matpoint import material_point_derivative_check
        params = [{"name": name, "index": int(info["index"]),
                   "sdv_range": list(info["sdv_range"])}
                  for name, info in gen.items()]
        ddsdde = (result.manifest or {}).get("derivatives", {}).get("ddsdde", {})
        dd_range = ddsdde.get("sdv_range") if ddsdde.get("in_sdv") else None
        verdict = material_point_derivative_check(
            original_for=result.original_umat,
            transformed_plain=result.transformed_plain,
            output_dir=result.output_dir,
            ntens=result.ntens,
            props=list(result.props_values),
            parameters=params,
            ddsdde_sdv_range=dd_range)
        if not verdict.get("available"):
            _record_numerical_validation(result.manifest, "unavailable", verdict)
            return TransformationValidation(
                passed=False, checks=verdict,
                message="real-response check unavailable: %s"
                        % verdict.get("reason", "unknown"))
        passed = bool(verdict.get("passed"))
        _record_numerical_validation(
            result.manifest, "passed" if passed else "failed", verdict)
        return TransformationValidation(
            passed=passed, checks=verdict,
            message="d(sigma)/da vs finite-difference max relerr = %.3e (tol %.1e)"
                    % (verdict.get("max_relerr", float("nan")),
                       verdict.get("rel_tol", 0.0)))

    # -- handoff manifest ------------------------------------------------- #
    def _manifest(self, proposal, result, src_file) -> Dict[str, Any]:
        # Honesty: report only what was actually GENERATED. When the transform
        # succeeds and parameters were seeded, the transformed UMAT writes the real
        # DDSDDE tangent to STATEV[1..36] and each d(sigma)/da_i to its STATEV block;
        # the manifest advertises exactly those SDV ranges so downstream never
        # prepares against SDV blocks that were not written.
        ok = bool(result.success)
        generated_params: Dict[str, Any] = {}
        if ok and proposal.parameters:
            for k, p in enumerate(proposal.parameters):
                lo = _DDSDDE_SLOTS + 1 + k * _PER_PARAM
                generated_params[p.name] = {
                    "source": p.source, "index": int(p.index), "status": "generated",
                    "target": "STATEV", "sdv_range": [lo, lo + _PER_PARAM - 1],
                    "oti_direction": "E%d" % (proposal.ntens + 1 + k),
                }
        requested = {p.name: {"source": p.source, "index": int(p.index),
                              "status": ("generated" if p.name in generated_params
                                         else "not_generated")}
                     for p in proposal.parameters}
        # The DDSDDE->STATEV[1..36] copy is emitted by build_parameter_contracts ONLY
        # when a tangent variable is known (real_copy_map from it); gate the manifest
        # on the SAME condition so it never advertises an SDV block that was never
        # written. propose_contract already requires tangent_variable with params, so
        # this is defence-in-depth.
        ddsdde_in_sdv = bool(generated_params and proposal.tangent_variable)
        required_depvar = (_DDSDDE_SLOTS if ddsdde_in_sdv else 0) \
            + _PER_PARAM * len(generated_params)
        return {
            "schema": HANDOFF_SCHEMA,
            "transformer": {"name": "umat-oti", "src": self.src_dir,
                            "contract_hash": _hash(proposal.config)},
            "source": {"original": os.path.basename(src_file),
                       "transformed": (os.path.basename(result.transformed_plain)
                                       if result.transformed_plain else None),
                       "combined": (os.path.basename(result.transformed_umat)
                                    if result.transformed_umat else None)},
            "capabilities": {"automatic_ddsdde": ok,
                             "parameter_derivatives": bool(generated_params)},
            "abaqus": {"required_depvar": required_depvar,     # generated SDVs only
                       "required_field_outputs": ["U", "S", "SDV"],
                       "compiler_form": "free"},
            "derivatives": {
                "ddsdde": {"status": "generated" if ok else "failed",
                           "target": "STATEV" if ddsdde_in_sdv else "DDSDDE",
                           "in_sdv": ddsdde_in_sdv,
                           "sdv_range": [1, _DDSDDE_SLOTS] if ddsdde_in_sdv else None},
                "parameters": generated_params,      # ONLY generated derivatives
            },
            "requested_parameter_derivatives": requested,   # requested == generated on success
            # STRUCTURAL codegen ("generated") is distinct from NUMERICAL agreement.
            # This block starts "not_run"; the numerical real-response check
            # (validate_real_response) records its verdict here so a downstream
            # consumer can tell a numerically-verified derivative from an unverified
            # one. See update_numerical_validation().
            "numerical_validation": {"status": "not_run"},
            "transformation": {"success": result.success,
                               "warnings": result.warnings,
                               "semantic_checks": result.semantic_checks},
        }


# --------------------------------------------------------------------------- #
def _config_source_file(config: Dict[str, Any]) -> Optional[str]:
    src = config.get("source")
    if isinstance(src, dict):
        # the loader resolves the compact 'source.file' to an absolute path here
        return (src.get("selected_umat_file") or src.get("file") or src.get("path"))
    return str(src) if src else None


def _record_numerical_validation(manifest: Optional[Dict[str, Any]],
                                 status: str, verdict: Dict[str, Any]) -> None:
    """Stamp the numerical real-response verdict into the handoff manifest so a
    consumer reading only the manifest can distinguish a numerically-verified
    derivative from an unverified/skipped one (status: passed|failed|unavailable)."""
    if not manifest:
        return
    entry: Dict[str, Any] = {"status": status}
    if verdict.get("available"):
        entry.update({"max_relerr": verdict.get("max_relerr"),
                      "rel_tol": verdict.get("rel_tol"),
                      "method": verdict.get("note")})
    else:
        entry["reason"] = verdict.get("reason")
    manifest["numerical_validation"] = entry


def _first_match(paths, needle):
    for p in paths:
        if needle in str(p):
            return p
    return None


def _hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str)
                          .encode("utf-8")).hexdigest()[:16]


def write_handoff_manifest(path: str, manifest: Dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
