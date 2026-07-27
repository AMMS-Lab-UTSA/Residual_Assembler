#!/usr/bin/env python3
"""PROGRAM 1 CLI -- UMAT-to-OTI transformation to a distributable .obj.

    python -m oti_provider.umat_transform build   <compact_contract.json>
    python -m oti_provider.umat_transform validate-all

Public product (per material):
    umat_<name>_oti.obj     -- relocatable object with THREE Fortran symbols:
                                 UMAT               (the original, unchanged),
                                 UMAT_OTI_INTERNAL  (the transformed OTI UMAT),
                                 UMAT_OTI_EVAL      (offline entry returning the
                                                     derivative arrays SEPARATELY),
                               plus the bundled OTI runtime.
    umat_<name>_oti.json    -- completed interface contract (symbols + layouts).

DSIGMA_DP / DSTATEV_DP are returned as dedicated arguments of UMAT_OTI_EVAL and
are NEVER stored in STATEV; STATEV_OUT is the physical material state only.

The compact user contract carries only: schema, source, kinematics, dimensions,
parameters (name + props_index), derivative, history, output. Everything else
(stress-update line, DDSDDE block, promoted variables, OTI directions) is
inferred automatically; parameter order fixes the derivative-column order.

Reuses the existing umat_oti transformer (src/umat_oti). No mesh/assembly/solve.
"""
import ctypes
import glob
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import contextlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
if os.path.join(REPO, "src") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "src"))
CONTRACT_DIR = os.path.join(HERE, "contract")
ABA_PARAM = "      implicit real*8(a-h,o-z)\n      parameter (nprecd=2)\n"
TRANSFORM_SCHEMA = "resasm_umat_transform_v2"

from umat_oti.runtime import (                                          # noqa: E402
    LibraryLoadError,
    binary_metadata,
    load_shared_library,
    prefer_static_runtime,
    provenance_block,
    shared_library_suffix,
    static_fortran_link_flags,
)
from umat_oti.validation import fd_reference as fdref                   # noqa: E402

#: methodology version for the validation record -- bump when the FD reference
#: method or metric definitions change, so archived results stay identifiable.
VALIDATION_METHOD_VERSION = "adaptive-fd-v1"

#: Precise definition of every metric the self-check reports. Embedded verbatim
#: in each validation JSON so a number is never separated from what it means.
METRIC_DEFINITIONS = {
    "stress_parity_max_rel": {
        "name": "primal stress equivalence (max-norm relative)",
        "formula": "max_over(path,components) |sigma_OTI - sigma_orig| / max(|sigma_orig|)",
        "aggregation": "max component, max over the full validation path",
        "normalization": "max |sigma_orig| over the path (floor 1e-30)",
        "kind": "equivalence to the original UMAT (same code path); analytical, not FD",
        "pass_if": "<= 1e-8",
    },
    "statev_parity_max_rel": {
        "name": "history/state equivalence (max-norm relative)",
        "formula": "max_over(path,components) |STATEV_OTI - STATEV_orig| / max(|STATEV_orig|)",
        "aggregation": "max component, max over the full validation path",
        "normalization": "max |STATEV_orig| over the path (floor 1e-30)",
        "kind": "equivalence to the original UMAT; analytical, not FD",
        "pass_if": "<= 1e-8",
    },
    "ddsdde_parity_max_rel": {
        "name": "consistent tangent vs converged FD of DSTRAN",
        "formula": "max |DDSDDE_OTI - dSigma/dDSTRAN_FD| / max(|dSigma/dDSTRAN_FD|)",
        "aggregation": "max component at the last (plastic) increment",
        "normalization": "max |FD tangent| (floor 1e-30)",
        "kind": "OTI-FD relative discrepancy (FD wrt strain increment)",
        "pass_if": "<= 1e-4",
    },
    "dsigma_max_rel": {
        "name": "DSIGMA_DP vs converged FD (OTI-FD relative discrepancy)",
        "formula": "max_c max_over(path,components) |DSIGMA_DP_OTI[:,c] - F*(c)| / max(|F*(c)|)",
        "aggregation": "max component, max over the full path, worst over parameters",
        "normalization": "per parameter: max |converged FD column| over the path (floor 1e-30)",
        "near_zero_reference": "denominator floored at 1e-30; a near-zero column inflates the "
                               "ratio and is reported, never silently zeroed",
        "fd_reference": "parameter-scaled centered difference, step selected per parameter from "
                        "the FD sequence (two-sided plateau), Richardson-checked; ladder 1e-3..1e-8",
        "kind": "OTI-FD relative discrepancy (no analytical derivative available)",
        "pass_if": "<= 1e-4  AND FD reference converged",
    },
    "dstatev_max_rel": {
        "name": "DSTATEV_DP vs converged FD (OTI-FD relative discrepancy)",
        "formula": "max_c max_over(path,components) |DSTATEV_DP_OTI[:,c] - F*(c)| / max(|F*(c)|)",
        "aggregation": "max component, max over the full path, worst over parameters",
        "normalization": "per parameter: max |converged FD column| over the path (floor 1e-30)",
        "fd_reference": "same canonical converged FD as dsigma_max_rel",
        "kind": "OTI-FD relative discrepancy",
        "pass_if": "<= 1e-4  AND FD reference converged",
    },
}


class ContractError(ValueError):
    """A transform contract is malformed, or lacks what validation needs."""


# --------------------------------------------------------------------------- #
# contract discovery and normalisation
#
# The material developer writes a COMPACT contract carrying only genuine
# decisions: the UMAT source, the interface sizes, which PROPS to differentiate,
# and which derivatives are wanted.  Everything else -- OTI directions and order,
# derivative array shapes, output filenames, ABI/binary metadata, tolerances,
# the FD ladder, and the completed detailed contract -- is inferred or generated.
#
#   {
#     "material": "j2_linear_hardening",          # optional label
#     "source": "umat_j2.for",                    # string, or {main, additional_files}
#     "kinematics": "small_strain",               # small_strain | finite_strain
#     "dimensions": {"stress_components": 6,       # (aliases: ntens / nprops / nstatev)
#                    "material_properties": 4,
#                    "state_variables": 1},
#     "parameters": {"E": 1, "nu": 2, ...},       # name -> PROPS index (nparam INFERRED)
#     "derivatives": {"stress": true, "statev": "auto"}
#   }
#
# Property VALUES do not belong here -- they are a validation operating point, not
# a transformation decision.  They live in a separate verification case
# ({"props": [...], "loading_path": "..."}); see resolve_props().
#
# Three input shapes are accepted and normalised to ONE internal representation:
#   * new compact  -- `parameters` is a dict (the documented format above);
#   * compact v2   -- `parameters` is a list of {name, props_index} (older);
#   * verbose      -- `interface` + `derivative_requests` (m2_elastic3d).
# --------------------------------------------------------------------------- #

#: dimension aliases: the clear names, plus the raw Fortran names for back-compat.
_DIM_ALIAS = {
    "stress_components": "ntens", "ntens": "ntens",
    "material_properties": "nprops", "nprops": "nprops",
    "state_variables": "nstatev", "nstatev": "nstatev",
}
#: keys that must never be written by the user -- nparam is inferred.
_FORBIDDEN_DIM = ("nparam", "n_param", "n_parameters", "parameters_count", "num_parameters")


def _norm_dimensions(dims, path):
    if not isinstance(dims, dict) or not dims:
        raise ContractError("%s: 'dimensions' is required (stress_components, "
                            "material_properties, state_variables)" % path)
    for bad in _FORBIDDEN_DIM:
        if bad in dims:
            raise ContractError("%s: dimensions.%s must not be set -- the number of "
                                "differentiated parameters (nparam) is inferred as "
                                "len(parameters)" % (path, bad))
    out = {}
    for k, v in dims.items():
        key = _DIM_ALIAS.get(k)
        if key is None:
            raise ContractError("%s: unknown dimension %r (use stress_components, "
                                "material_properties, state_variables)" % (path, k))
        out[key] = int(v)
    if "ntens" not in out or "nprops" not in out:
        raise ContractError("%s: dimensions must give stress_components and "
                            "material_properties" % path)
    out.setdefault("nstatev", 0)
    return out


def _norm_parameters(params, nprops, path):
    """Accept dict {name: index} or list [{name, props_index}] -> ordered list.

    The listed order fixes each parameter's OTI direction (1..nparam).  Python
    dicts preserve insertion order, so the developer's ordering is honoured.
    """
    if isinstance(params, dict):
        items = [{"name": str(n), "props_index": int(i)} for n, i in params.items()]
    elif isinstance(params, list):
        items = []
        for p in params:
            if "name" not in p or ("props_index" not in p and "index" not in p):
                raise ContractError("%s: each parameter needs a name and a PROPS index" % path)
            items.append({"name": str(p["name"]),
                          "props_index": int(p.get("props_index", p.get("index")))})
    else:
        raise ContractError("%s: 'parameters' must be {name: index} or a list of "
                            "{name, props_index}" % path)
    if not items:
        raise ContractError("%s: at least one parameter to differentiate is required" % path)
    idxs = [p["props_index"] for p in items]
    out_of_range = [i for i in idxs if not (1 <= i <= nprops)]
    if out_of_range:
        raise ContractError("%s: PROPS index %s outside 1..%d (material_properties)"
                            % (path, out_of_range, nprops))
    if len(set(idxs)) != len(idxs):
        raise ContractError("%s: duplicate PROPS index in parameters" % path)
    names = [p["name"] for p in items]
    if len(set(names)) != len(names):
        raise ContractError("%s: duplicate parameter name" % path)
    for k, p in enumerate(items):
        p["oti_direction"] = k + 1          # inferred, never user-supplied
    return items


def _norm_source(src, path):
    """Accept 'umat.for' or {main[/main_file], additional_files} -> internal form."""
    if isinstance(src, str):
        return {"main_file": src, "additional_files": []}
    if isinstance(src, dict):
        main = src.get("main") or src.get("main_file")
        if not main:
            raise ContractError("%s: source needs a 'main' (or 'main_file') UMAT file" % path)
        return {"main_file": str(main),
                "additional_files": list(src.get("additional_files") or [])}
    raise ContractError("%s: 'source' must be a filename or {main, additional_files}" % path)


def _norm_derivatives(deriv, nstatev, path):
    """Requested derivative outputs -> {stress, statev, path_dependent}.

    ``statev: "auto"`` (the default) needs DSTATEV_DP exactly when the material
    carries state.  ``true`` forces it (requires state); ``false`` forbids it
    (only consistent when the material is stateless -- a material with state is
    path-dependent and its stress derivative needs the state derivative chained).
    """
    deriv = deriv or {}
    if deriv.get("stress", True) is not True:
        raise ContractError("%s: derivatives.stress must be true -- the stress "
                            "parameter derivative is the framework's purpose" % path)
    statev = deriv.get("statev", "auto")
    if statev == "auto":
        path_dependent = nstatev > 0
    elif statev is True:
        if nstatev <= 0:
            raise ContractError("%s: derivatives.statev is true but state_variables=0" % path)
        path_dependent = True
    elif statev is False:
        if nstatev > 0:
            raise ContractError("%s: derivatives.statev is false but state_variables=%d>0; "
                                "a material with state is path-dependent -- use \"auto\""
                                % (path, nstatev))
        path_dependent = False
    else:
        raise ContractError("%s: derivatives.statev must be true, false, or \"auto\"" % path)
    return {"stress": True, "statev": statev, "path_dependent": path_dependent}


def _reject_manual_nparam(raw, nparam, path):
    for key in _FORBIDDEN_DIM:
        if key in raw:
            declared = int(raw[key])
            raise ContractError("%s: %r must not be set -- nparam is inferred as "
                                "len(parameters)=%d (a manual value can contradict the "
                                "parameter list; declared %d)" % (path, key, nparam, declared))


def _finalize(*, material, source_raw, kinematics, dims_raw, params_raw,
              derivatives_raw, validation_props, options, path):
    """Build the one internal representation every input shape normalises to."""
    dims = _norm_dimensions(dims_raw, path)
    params = _norm_parameters(params_raw, dims["nprops"], path)
    src = _norm_source(source_raw, path)
    deriv = _norm_derivatives(derivatives_raw, dims["nstatev"], path)
    options = options or {}
    kin = options.get("kinematics") or kinematics or "small_strain"
    if kin not in ("small_strain", "finite_strain"):
        raise ContractError("%s: kinematics must be small_strain or finite_strain" % path)
    out = {
        "schema": TRANSFORM_SCHEMA,
        "_material": material,
        "source": src,
        "kinematics": kin,
        "dimensions": {k: dims[k] for k in ("ntens", "nprops", "nstatev")},
        "parameters": params,
        "_derivatives": deriv,
    }
    if validation_props is not None:                 # embedded props (back-compat only)
        out["_verification"] = {"props_values": validation_props}
    override = {}
    if options.get("stress_update_line"):
        override["stress_update_line"] = int(options["stress_update_line"])
    if options.get("ddsdde_block"):
        override["ddsdde_block"] = str(options["ddsdde_block"])
    if options.get("entry_point"):
        override["entry_point"] = str(options["entry_point"])
    if override:
        out["_auto_override"] = override
    return out


def normalize_contract(raw, path):
    """Return (internal_contract, input_format) for any accepted input shape."""
    if not isinstance(raw, dict):
        raise ContractError("%s: contract must be a JSON object" % path)
    if raw.get("schema") == "resasm_umat_oti_contract_v1":
        raise ContractError("%s: this is a COMPLETED (generated) contract, not an input "
                            "transformation contract" % path)
    # verbose legacy: interface + derivative_requests
    if "interface" in raw and "derivative_requests" in raw:
        return _from_verbose(raw, path), "verbose"
    if "parameters" not in raw:
        raise ContractError("%s: not a transformation contract (no 'parameters')" % path)
    _reject_manual_nparam(raw, len(raw["parameters"]), path)
    validation_props = (raw.get("validation") or {}).get("props_values")
    fmt = "compact" if isinstance(raw["parameters"], dict) else "compact-v2"
    internal = _finalize(
        material=raw.get("material"),
        source_raw=raw.get("source"),
        kinematics=raw.get("kinematics"),
        dims_raw=raw.get("dimensions"),
        params_raw=raw["parameters"],
        derivatives_raw=raw.get("derivatives") or raw.get("derivative"),
        validation_props=validation_props,
        options=raw.get("options"),
        path=path,
    )
    return internal, fmt


def _from_verbose(raw, path):
    """Normalise the legacy verbose layout (interface + derivative_requests)."""
    iface = raw["interface"]
    reqs = raw["derivative_requests"]
    if not reqs:
        raise ContractError("%s: derivative_requests is empty" % path)
    comps = (reqs[0].get("seed") or {}).get("components") or []
    if not comps:
        raise ContractError("%s: derivative_requests[0].seed.components is empty" % path)
    params = [{"name": c["name"], "props_index": int(c["index"])}
              for c in sorted(comps, key=lambda c: int(c.get("oti_direction", 0)))]
    prov = raw.get("resasm_provider") or {}
    options = {}
    if prov.get("stress_update_line"):
        options["stress_update_line"] = int(prov["stress_update_line"])
    if prov.get("ddsdde_block"):
        options["ddsdde_block"] = str(prov["ddsdde_block"])
    return _finalize(
        material=raw.get("material"),
        source_raw=raw["source"],
        kinematics=iface.get("kinematics"),
        dims_raw={"ntens": iface["ntens"], "nprops": iface["nprops"],
                  "nstatev": iface.get("nstatev", 0)},
        params_raw=params,
        derivatives_raw={"stress": True, "statev": "auto"},
        validation_props=prov.get("props_values"),
        options=options,
        path=path,
    )


#: verification-case filenames (they carry props, never parameters).
_VERIFY_NAMES = ("verification.json", "verify.json")


def _looks_like_transform_contract(raw):
    """True if *raw* is an input transformation contract (any accepted shape)."""
    if not isinstance(raw, dict):
        return False
    if raw.get("schema") == "resasm_umat_oti_contract_v1":   # completed output
        return False
    if "props" in raw and "parameters" not in raw:           # a verification case
        return False
    if "interface" in raw and "derivative_requests" in raw:  # verbose
        return True
    # compact (new or v2): parameters + something that pins the interface
    return "parameters" in raw and ("source" in raw or "dimensions" in raw)


def discover_contracts(materials_dir):
    """Every transform contract under *materials_dir*, found by structure.

    Recognises all accepted input shapes (compact / compact-v2 / verbose) by
    content, not by a fixed filename or schema tag, so a new-format contract that
    omits the schema is still found.  Completed output contracts (``umat_*``) and
    verification cases (``props`` without ``parameters``) are excluded.

    Returns ``(entries, skipped)``: entries are ``{"material", "path", "format",
    "contract"}``; skipped are ``{"material", "reason"}`` -- nothing vanishes
    without a stated reason.
    """
    entries, skipped = [], []
    for name in sorted(os.listdir(materials_dir)):
        mdir = os.path.join(materials_dir, name)
        if not os.path.isdir(mdir):
            continue
        matches, notes = [], []
        for cand in sorted(glob.glob(os.path.join(mdir, "*.json"))):
            base = os.path.basename(cand)
            if base.startswith("umat_") or base in _VERIFY_NAMES or base.endswith(".verify.json"):
                continue
            try:
                raw = json.load(open(cand))
            except Exception as exc:
                notes.append("%s unreadable (%s)" % (base, exc))
                continue
            if _looks_like_transform_contract(raw):
                matches.append((cand, raw))
        if not matches:
            skipped.append({"material": name,
                            "reason": "; ".join(notes) or "no transformation contract found"})
            continue
        preferred = {"contract.json": 0, "transform_contract_v2.json": 1}
        matches.sort(key=lambda m: preferred.get(os.path.basename(m[0]), 2))
        path, raw = matches[0]
        try:
            contract, fmt = normalize_contract(raw, path)
        except ContractError as exc:
            skipped.append({"material": name, "reason": str(exc)})
            continue
        entries.append({"material": name, "path": path, "format": fmt, "contract": contract})
    return entries, skipped


# --------------------------------------------------------------------------- #
# the VERIFICATION CASE: property values + loading history for validation.
#
# Property values describe WHERE to test the transform, not HOW to transform it,
# so they live outside the transformation contract, in a small verification case:
#
#     {"props": [210000.0, 0.30, 250.0, 1000.0],
#      "loading_path": "uniaxial_tension.json"}   # optional
#
# Resolution precedence:  --props  ->  --verify <file>  ->  a sibling
# verification.json / *.verify.json  ->  props embedded in an older contract
# (back-compat)  ->  an explicit error.  There is deliberately no [1.0]*nprops
# fallback: validating at placeholder parameters certifies the transform at an
# operating point the model never sees (for crystal plasticity it diverges to
# NaN).  A missing operating point is a defect, raised as one.
# --------------------------------------------------------------------------- #
def find_verification_case(cdir):
    """The verification-case file beside a contract, or None."""
    if not cdir:
        return None
    for base in _VERIFY_NAMES:
        p = os.path.join(cdir, base)
        if os.path.exists(p):
            return p
    extra = sorted(glob.glob(os.path.join(cdir, "*.verify.json")))
    return extra[0] if extra else None


def _read_verification_case(path):
    raw = json.load(open(path))
    props = raw.get("props", raw.get("props_values"))
    return props, raw.get("loading_path")


def resolve_props(contract, nprops, cli_props=None, cdir=None, verify_path=None):
    """Return (props_values, origin, loading_path) for validation."""
    loading_path = None
    if cli_props is not None:
        values, origin = list(cli_props), "--props"
    elif verify_path:
        values, loading_path = _read_verification_case(verify_path)
        origin = "verification case %s (--verify)" % os.path.basename(verify_path)
    else:
        sibling = find_verification_case(cdir)
        if sibling:
            values, loading_path = _read_verification_case(sibling)
            origin = "verification case %s" % os.path.basename(sibling)
        else:
            values = (contract.get("_verification") or {}).get("props_values")
            origin = "contract (embedded validation props, legacy)"
            if values is None:
                raise ContractError(
                    "no property values to validate at.\n"
                    "         add a verification case beside the contract:\n"
                    "             verification.json = {\"props\": [ ... %d entries ... ]}\n"
                    "         or pass  --props v1,v2,...  (or --verify <file>).\n"
                    "         Property VALUES are a validation operating point, not a\n"
                    "         transformation decision, so they live outside the contract."
                    % nprops)

    if values is None:
        raise ContractError("%s: no 'props' found" % origin)
    try:
        values = [float(v) for v in values]
    except (TypeError, ValueError) as exc:
        raise ContractError("%s: property values are not all numbers (%s)" % (origin, exc))
    if len(values) != nprops:
        raise ContractError("%s: %d property values but the contract declares "
                            "material_properties=%d" % (origin, len(values), nprops))
    bad = [i + 1 for i, v in enumerate(values) if not math.isfinite(v)]
    if bad:
        raise ContractError("%s: non-finite value at PROPS index %s"
                            % (origin, ", ".join(str(b) for b in bad)))
    return values, origin, loading_path


# --------------------------------------------------------------------------- #
# auto-inference from the fixed-form UMAT source (no user line numbers)
# --------------------------------------------------------------------------- #
def _auto_infer(umat_path, dependent="STRESS", tangent="DDSDDE"):
    lines = open(umat_path).read().splitlines()
    stress_line = ddsdde_line = None
    for i, ln in enumerate(lines, 1):
        code = ln[6:] if len(ln) > 6 and ln[:1] in " 0123456789" else ln
        if re.match(r"\s*%s\s*\(" % dependent, code) and "=" in code:
            stress_line = i                     # last STRESS(...) = ...
        if re.match(r"\s*%s\s*\(" % tangent, code) and "=" in code:
            ddsdde_line = i
    if stress_line is None:
        raise ValueError("could not locate a %s(...) = assignment in %s" % (dependent, umat_path))
    if ddsdde_line is None:
        raise ValueError("could not locate a %s(...) = assignment in %s" % (tangent, umat_path))
    return {"stress_update_line": stress_line, "ddsdde_block": "%d-%d" % (ddsdde_line, ddsdde_line)}


# --------------------------------------------------------------------------- #
# compact v2 contract -> the transformer's internal config
# --------------------------------------------------------------------------- #
def _transform_config(contract, umat_abs, auto):
    dims = contract["dimensions"]
    ntens = int(dims["ntens"])
    params = contract["parameters"]              # [{name, props_index}], order == direction
    # each parameter's d(sigma)/d(p) -> a scratch STATEV block (demuxed later)
    contracts = []
    for k, p in enumerate(params):
        lo = ntens * k + 1
        comps = [{"target_indices": [lo + i - 1], "output_indices": [i],
                  "seed_direction_offset": 0} for i in range(1, ntens + 1)]
        contracts.append({
            "id": "dsigma_d%s" % p["name"],
            "seed": {"variable": "PROPS", "shape": "vector", "directions": 1,
                     "components": [[int(p["props_index"])]],
                     "operating_point_expression": "PROPS(%d)" % int(p["props_index"])},
            "output": {"variable": "STRESS", "shape": "vector"},
            "internal_use": {"replace_variable": "PROPS"},
            "additional_extractions": [{
                "target_variable": "STATEV", "from_output_variable": "STRESS",
                "after_line": auto["stress_update_line"], "extract_kind": "component_map",
                "components": comps}]})
    return {
        "name": contract["_name"],
        "source": {"file": umat_abs},
        "ntens": ntens, "order": 1,
        "jacobian": {"independent": "DSTRAN", "dependent": "STRESS", "target": "DDSDDE"},
        "promote": ["PROPS", "STRESS", "DDSDDE"],   # closure auto-adds the rest
        "replace": {"ddsdde_block": [auto["ddsdde_block"]]},
        "extra_jacobian_contracts": contracts,
    }


# --------------------------------------------------------------------------- #
# path-dependent (NSTATEV>0) transform config: seed PROPS + STATEV_in + STRESS_in
# so the eval can chain total dsigma/dp and dstatev/dp across increments.
#
# Scratch STATEV layout (1-based):
#   physical            1 .. NS
#   then one (NT+NS) block per seed contract, contracts in this order:
#     PROPS(param k)    k = 1..nparam
#     STATEV(j)         j = 1..NS        (state-transition wrt incoming state)
#     STRESS(i)         i = 1..NT        (state-transition wrt incoming stress)
#   within a block: d(sigma)/d(seed)  -> base+1 .. base+NT   (output 1..NT)
#                   d(statev)/d(seed)  -> base+NT+1 .. base+NT+NS (output 1..NS)
# --------------------------------------------------------------------------- #
def _transform_config_path(contract, umat_abs, auto):
    dims = contract["dimensions"]
    ntens = int(dims["ntens"]); nstatev = int(dims["nstatev"])
    params = contract["parameters"]
    sl = auto["stress_update_line"]
    blk = ntens + nstatev

    def _sig_ext(base):
        return [{"target_indices": [base + 1 + i], "output_indices": [i + 1],
                 "seed_direction_offset": 0} for i in range(ntens)]

    def _sv_ext(base):
        return [{"target_indices": [base + ntens + 1 + m], "output_indices": [m + 1],
                 "seed_direction_offset": 0} for m in range(nstatev)]

    seeds = ([("PROPS", int(p["props_index"]), "dsig_dp_%s" % p["name"]) for p in params]
             + [("STATEV", j + 1, "dsig_dsv%d" % (j + 1)) for j in range(nstatev)]
             + [("STRESS", i + 1, "dsig_dst%d" % (i + 1)) for i in range(ntens)])
    contracts = []
    for c, (var, comp, cid) in enumerate(seeds):
        base = nstatev + c * blk
        contracts.append({
            "id": cid,
            "seed": {"variable": var, "shape": "vector", "directions": 1,
                     "components": [[comp]],
                     "operating_point_expression": "%s(%d)" % (var, comp)},
            "output": {"variable": "STRESS", "shape": "vector"},
            "internal_use": {"replace_variable": var},
            "additional_extractions": [
                {"target_variable": "STATEV", "from_output_variable": "STRESS",
                 "after_line": sl, "extract_kind": "component_map", "components": _sig_ext(base)},
                {"target_variable": "STATEV", "from_output_variable": "STATEV",
                 "after_line": sl, "extract_kind": "component_map", "components": _sv_ext(base)}]})
    return {
        "name": contract["_name"],
        "source": {"file": umat_abs},
        "ntens": ntens, "order": 1,
        "jacobian": {"independent": "DSTRAN", "dependent": "STRESS", "target": "DDSDDE"},
        "promote": ["PROPS", "STRESS", "STATEV", "DDSDDE"],
        "replace": {"ddsdde_block": [auto["ddsdde_block"]]},
        "extra_jacobian_contracts": contracts,
    }


# --------------------------------------------------------------------------- #
# generated Fortran: UMAT_OTI_EVAL, demuxing the scratch STATEV into separate
# DSIGMA_DP / DSTATEV_DP (physical STATEV stays clean)
# --------------------------------------------------------------------------- #
def _gen_eval(ntens, nprops, nstatev, nparam):
    nscratch = nstatev + ntens * nparam
    return r"""! AUTO-GENERATED offline OTI evaluation entry point. Do not edit.
      SUBROUTINE UMAT_OTI_EVAL(STRESS, STATEV, DDSDDE, STRAN, DSTRAN,
     1 TIME, DTIME, TEMP, DTEMP, PROPS, NPROPS, NTENS, NSTATV,
     2 NPARAM, DSIGMA_DP, DSTATEV_DP)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
      DIMENSION STRAN(NTENS), DSTRAN(NTENS), TIME(2), PROPS(NPROPS)
      DIMENSION DSIGMA_DP(NTENS,NPARAM), DSTATEV_DP(NSTATV,NPARAM)
      DIMENSION STATEV_SCR(%(NSCR)d), PREDEF(1), DPRED(1)
      DIMENSION DDSDDT(NTENS), DRPLDE(NTENS), COORDS(3), DROT(3,3)
      DIMENSION DFGRD0(3,3), DFGRD1(3,3)
      CHARACTER*80 CMNAME
C     zero ALL scratch passed to the UMAT (uninitialised scratch -> NaN).
      SSE=0.D0; SPD=0.D0; SCD=0.D0; RPL=0.D0; DRPLDT=0.D0
      DO I = 1, %(NSCR)d
        STATEV_SCR(I) = 0.D0
      END DO
      DO I = 1, NSTATV
        STATEV_SCR(I) = STATEV(I)
      END DO
      DO I = 1, NTENS
        DDSDDT(I)=0.D0
        DRPLDE(I)=0.D0
      END DO
      DO I = 1, 3
        COORDS(I)=0.D0
      END DO
      PREDEF(1)=0.D0
      DPRED(1)=0.D0
      DO I = 1, 3
        DO J = 1, 3
          DROT(I,J)=0.D0
          DFGRD0(I,J)=0.D0
          DFGRD1(I,J)=0.D0
        END DO
        DROT(I,I)=1.D0
        DFGRD0(I,I)=1.D0
        DFGRD1(I,I)=1.D0
      END DO
      PNEWDT = 1.D0
      CMNAME = 'MAT'
      CALL UMAT_OTI_INTERNAL(STRESS, STATEV_SCR, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,
     3 3, NTENS-3, NTENS, %(NSCR)d, PROPS, NPROPS, COORDS, DROT, PNEWDT,
     4 1.D0, DFGRD0, DFGRD1, 1, 1, 1, 1, 1, 1)
C     physical state back (clean); derivative blocks OUT separately.
      DO I = 1, NSTATV
        STATEV(I) = STATEV_SCR(I)
      END DO
      DO J = 1, NPARAM
        DO I = 1, NTENS
          DSIGMA_DP(I,J) = STATEV_SCR(NSTATV + NTENS*(J-1) + I)
        END DO
        DO I = 1, NSTATV
          DSTATEV_DP(I,J) = 0.D0
        END DO
      END DO
      RETURN
      END
""" % {"NSCR": nscratch}


# --------------------------------------------------------------------------- #
# path-dependent eval: demux local partials from scratch and CHAIN with the
# incoming carried sensitivities (DSIGMA_DP_IN, DSTATEV_DP_IN) to return the
# TOTAL dsigma/dp and dstatev/dp for this increment. STATEV_OUT is physical.
# --------------------------------------------------------------------------- #
def _gen_eval_path(ntens, nprops, nstatev, nparam):
    nseed = nparam + nstatev + ntens
    blk = ntens + nstatev
    nscratch = nstatev + nseed * blk
    return r"""! AUTO-GENERATED offline OTI evaluation entry point (path-dependent). Do not edit.
      SUBROUTINE UMAT_OTI_EVAL(STRESS, STATEV, DDSDDE, STRAN, DSTRAN,
     1 TIME, DTIME, TEMP, DTEMP, PROPS, NPROPS, NTENS, NSTATV,
     2 NPARAM, DSIGMA_DP, DSTATEV_DP, DSIGMA_DP_IN, DSTATEV_DP_IN)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
      DIMENSION STRAN(NTENS), DSTRAN(NTENS), TIME(2), PROPS(NPROPS)
      DIMENSION DSIGMA_DP(NTENS,NPARAM), DSTATEV_DP(NSTATV,NPARAM)
      DIMENSION DSIGMA_DP_IN(NTENS,NPARAM), DSTATEV_DP_IN(NSTATV,NPARAM)
      DIMENSION STATEV_SCR(%(NSCR)d), PREDEF(1), DPRED(1)
      DIMENSION DDSDDT(NTENS), DRPLDE(NTENS), COORDS(3), DROT(3,3)
      DIMENSION DFGRD0(3,3), DFGRD1(3,3)
      CHARACTER*80 CMNAME
      INTEGER CP, CV, CS, IB, JB
      SSE=0.D0; SPD=0.D0; SCD=0.D0; RPL=0.D0; DRPLDT=0.D0
      DO I = 1, %(NSCR)d
        STATEV_SCR(I) = 0.D0
      END DO
      DO I = 1, NSTATV
        STATEV_SCR(I) = STATEV(I)
      END DO
      DO I = 1, NTENS
        DDSDDT(I)=0.D0
        DRPLDE(I)=0.D0
      END DO
      DO I = 1, 3
        COORDS(I)=0.D0
      END DO
      PREDEF(1)=0.D0
      DPRED(1)=0.D0
      DO I = 1, 3
        DO J = 1, 3
          DROT(I,J)=0.D0
          DFGRD0(I,J)=0.D0
          DFGRD1(I,J)=0.D0
        END DO
        DROT(I,I)=1.D0
        DFGRD0(I,I)=1.D0
        DFGRD1(I,I)=1.D0
      END DO
      PNEWDT = 1.D0
      CMNAME = 'MAT'
      CALL UMAT_OTI_INTERNAL(STRESS, STATEV_SCR, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,
     3 3, NTENS-3, NTENS, %(NSCR)d, PROPS, NPROPS, COORDS, DROT, PNEWDT,
     4 1.D0, DFGRD0, DFGRD1, 1, 1, 1, 1, 1, 1)
C     physical state back (clean).
      DO I = 1, NSTATV
        STATEV(I) = STATEV_SCR(I)
      END DO
C     Scratch block base for seed-contract C0 (0-based): NSTATV + C0*(NTENS+NSTATV)
C     sigma-partial slot for comp A:  base + A ; statev-partial slot for comp M: base + NTENS + M
C     Contract order: props k=1..NPARAM, then statev j=1..NSTATV, then stress i=1..NTENS.
      DO K = 1, NPARAM
        CP = NSTATV + (K-1)*(NTENS+NSTATV)
        DO IB = 1, NTENS
          DSIGMA_DP(IB,K) = STATEV_SCR(CP + IB)
        END DO
        DO IB = 1, NSTATV
          DSTATEV_DP(IB,K) = STATEV_SCR(CP + NTENS + IB)
        END DO
C       + incoming-state contributions (state-transition Jacobians)
        DO JB = 1, NSTATV
          CV = NSTATV + (NPARAM + JB - 1)*(NTENS+NSTATV)
          DO IB = 1, NTENS
            DSIGMA_DP(IB,K) = DSIGMA_DP(IB,K)
     1        + STATEV_SCR(CV + IB) * DSTATEV_DP_IN(JB,K)
          END DO
          DO IB = 1, NSTATV
            DSTATEV_DP(IB,K) = DSTATEV_DP(IB,K)
     1        + STATEV_SCR(CV + NTENS + IB) * DSTATEV_DP_IN(JB,K)
          END DO
        END DO
C       + incoming-stress contributions
        DO JB = 1, NTENS
          CS = NSTATV + (NPARAM + NSTATV + JB - 1)*(NTENS+NSTATV)
          DO IB = 1, NTENS
            DSIGMA_DP(IB,K) = DSIGMA_DP(IB,K)
     1        + STATEV_SCR(CS + IB) * DSIGMA_DP_IN(JB,K)
          END DO
          DO IB = 1, NSTATV
            DSTATEV_DP(IB,K) = DSTATEV_DP(IB,K)
     1        + STATEV_SCR(CS + NTENS + IB) * DSIGMA_DP_IN(JB,K)
          END DO
        END DO
      END DO
      RETURN
      END
""" % {"NSCR": nscratch}


# --------------------------------------------------------------------------- #
# EFFICIENT whole-path march: seed ONLY the parameters (+ DSTRAN for the tangent)
# and keep the OTI STRESS/STATEV alive across increments, so the OTI state carries
# the accumulated derivative automatically -- no per-increment re-seeding and no
# state-transition-Jacobian chaining. ntens+nparam directions instead of
# ntens+nparam+nstatev+ntens, and the whole path runs in one Fortran call.
# --------------------------------------------------------------------------- #
def _transform_config_march(contract, umat_abs, auto):
    dims = contract["dimensions"]; ntens = int(dims["ntens"])
    params = contract["parameters"]
    contracts = [{
        "id": "seed_%s" % p["name"],
        "seed": {"variable": "PROPS", "shape": "vector", "directions": 1,
                 "components": [[int(p["props_index"])]],
                 "operating_point_expression": "PROPS(%d)" % int(p["props_index"])},
        "output": {"variable": "STRESS", "shape": "vector"},
        "internal_use": {"replace_variable": "PROPS"},
        "additional_extractions": []} for p in params]
    return {
        "name": contract["_name"], "source": {"file": umat_abs},
        "ntens": ntens, "order": 1,
        "jacobian": {"independent": "DSTRAN", "dependent": "STRESS", "target": "DDSDDE"},
        "promote": ["PROPS", "STRESS", "STATEV", "DDSDDE"],
        "replace": {"ddsdde_block": [auto["ddsdde_block"]]},
        "extra_jacobian_contracts": contracts,
    }


def _gen_march(internal_text, ntens, nstatev, nparam):
    """Rewrite a param-only transformed internal into UMAT_OTI_MARCH: a routine
    that marches the whole recorded strain path in one call with the OTI state
    carried across increments, returning DSIGMA_DP and DDSDDE at each increment."""
    lines = internal_text.split("\n")

    def is_cont(l):
        if not l or l[:1] in "Cc*!":     # comment lines are never continuations
            return False
        return len(l) > 5 and l[5] not in (" ", "0", "\t")

    ui = next(i for i, l in enumerate(lines) if l.strip().startswith("USE otim"))
    uend = ui
    while uend + 1 < len(lines) and is_cont(lines[uend + 1]):
        uend += 1
    use_block = "\n".join(lines[ui:uend + 1])
    DECL_KW = ("DIMENSION", "TYPE", "INTEGER", "REAL", "DOUBLE", "CHARACTER", "PARAMETER",
               "LOGICAL", "COMPLEX", "IMPLICIT", "DATA", "COMMON", "SAVE", "EXTERNAL",
               "INTRINSIC", "USE", "INCLUDE", "EQUIVALENCE")

    def code_of(l):
        if not l.strip() or l[:1] in "Cc*!":
            return ""
        return (l[6:] if (len(l) > 6 and l[:1] in " 0123456789") else l.lstrip()).strip()

    def statements(lo, hi):
        i = lo
        while i < hi:
            stmt = [lines[i]]; j = i + 1
            while j < hi and is_cont(lines[j]):
                stmt.append(lines[j]); j += 1
            if code_of(lines[i]).upper().startswith("DO "):     # group DO..END DO (nested)
                depth = 1
                while j < hi and depth > 0:
                    stmt.append(lines[j]); cj = code_of(lines[j]).upper()
                    if cj.startswith("DO "):
                        depth += 1
                    elif cj.startswith("END DO") or cj.startswith("ENDDO"):
                        depth -= 1
                    j += 1
            yield i, stmt
            i = j

    entry_re = re.compile(r"(STRESS|STATEV|PROPS|DSTRAN|DDSDDE)_OTI\(OTI_I\)\s*=\s*"
                          r"(STRESS|STATEV|PROPS|DSTRAN|DDSDDE)\(")
    entryA = next(i for i, l in enumerate(lines) if re.match(r"\s*DO\s+OTI_I\b", l))
    seed_lines = [l for l in lines if "+ OTI_E" in l]
    dstran_seeds = [l for l in seed_lines if "DSTRAN_OTI" in l]
    param_seeds = [l for l in seed_lines if "PROPS_OTI" in l]
    last_seed = max(i for i, l in enumerate(lines) if "+ OTI_E" in l)
    copyback = next(i for i, l in enumerate(lines) if "Copy real-valued OTIS outputs" in l)
    body = lines[last_seed + 1:copyback]

    # classify the pre-entry statements: pure declarations vs zero-init (reset each
    # increment) vs pre-loop executables (e.g. a constant RMAT initialisation).
    decl, zero_loop, preseed, seen_exec = [], [], [], False
    for i, stmt in statements(0, entryA):
        up = " ".join(x.strip() for x in stmt).upper().replace(" ", "")
        if ui <= i <= uend or up.startswith("SUBROUTINE") or up.startswith("INCLUDE") or "STRESS(NTENS)" in up:
            continue
        code0 = code_of(stmt[0]).upper()
        if stmt[0][:1] in "Cc*!" or not stmt[0].strip():        # comment / blank
            if not seen_exec:
                decl += stmt
            continue
        if any(code0.startswith(k) for k in DECL_KW):           # declaration (incl. PARAMETER)
            decl += stmt
            continue
        # executable statement
        seen_exec = True
        if ("=0.0D0" in up or "=0.D0" in up) and "_OTI" in up:  # OTI zero-init
            if not any(p in up for p in ("PROPS_OTI(", "STATEV_OTI(", "STRESS_OTI(")):
                zero_loop += stmt                               # reset non-persistent temps
        else:
            preseed += stmt                                     # constant pre-loop init (RMAT)
    nl = chr(10)
    strip_s = nl.join("            STRESS_OTI(OTI_I)=STRESS_OTI(OTI_I)-"
                      "GETIM(STRESS_OTI(OTI_I),%d)*OTI_E%d" % (d, d) for d in range(1, ntens + 1))
    strip_v = nl.join("            STATEV_OTI(OTI_I)=STATEV_OTI(OTI_I)-"
                      "GETIM(STATEV_OTI(OTI_I),%d)*OTI_E%d" % (d, d) for d in range(1, ntens + 1))
    return f"""      SUBROUTINE UMAT_OTI_MARCH(PROPS,NPROPS,PATH,NPATH,DTARR,NTENS,
     1 NSTATV,NPARAM,DSIG,STROUT,DDOUT)
{use_block}
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS),PATH(NTENS,NPATH),DTARR(NPATH),
     1 DSIG(NTENS,NPARAM,NPATH),STROUT(NTENS),DDOUT(NTENS,NTENS,NPATH)
{nl.join(decl)}
      INTEGER KSTEP, IDIR, IPD, NDI, NSHR
      NDI=3
      NSHR=NTENS-3
      DO OTI_I=1,NPROPS
         PROPS_OTI(OTI_I)=PROPS(OTI_I)
      END DO
{nl.join(param_seeds)}
      DO OTI_I=1,NTENS
         STRESS_OTI(OTI_I)=0.0D0
      END DO
      DO OTI_I=1,NSTATV
         STATEV_OTI(OTI_I)=0.0D0
      END DO
{nl.join(preseed)}
      DO KSTEP=1,NPATH
         DTIME=DTARR(KSTEP)
{nl.join(zero_loop)}
         DO OTI_I=1,NTENS
            DSTRAN_OTI(OTI_I)=PATH(OTI_I,KSTEP)
         END DO
{nl.join(dstran_seeds)}
{nl.join(body)}
         DO OTI_I=1,NTENS
            DO IPD=1,NPARAM
               DSIG(OTI_I,IPD,KSTEP)=GETIM(STRESS_OTI(OTI_I),NTENS+IPD)
            END DO
            DO OTI_J=1,NTENS
               DDOUT(OTI_I,OTI_J,KSTEP)=GETIM(STRESS_OTI(OTI_I),OTI_J)
            END DO
         END DO
         DO OTI_I=1,NTENS
{strip_s}
         END DO
         DO OTI_I=1,NSTATV
{strip_v}
         END DO
      END DO
      DO OTI_I=1,NTENS
         STROUT(OTI_I)=REAL(STRESS_OTI(OTI_I))
      END DO
      RETURN
      END
"""


# --------------------------------------------------------------------------- #
def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def _gfc(args):
    subprocess.check_call([os.environ.get("FC") or "gfortran"] + args)


def _validation_lib(out, stem):
    """Path of an internal validation library, with the platform's suffix."""
    return os.path.join(out, stem + shared_library_suffix())


def _link_and_load(objects, lib_path):
    """Link a validation library and load it, cross-platform.

    The dynamic link is the primary path on both platforms -- it is what the
    collaborator's own toolchain will do with the shipped object.  Folding the
    gfortran runtime in statically is kept as a *fallback* (and as an explicit
    opt-in via ``UMAT_OTI_STATIC_FORTRAN=1``) so a machine with an unusual
    runtime layout still validates instead of failing outright.
    """
    def _link(static):
        flags = static_fortran_link_flags() if static else []
        _gfc(["-shared"] + flags + list(objects) + ["-o", lib_path])

    static_first = prefer_static_runtime()
    _link(static_first)
    try:
        return load_shared_library(lib_path)
    except LibraryLoadError as first:
        if static_first:
            raise
        try:
            _link(True)
            handle = load_shared_library(lib_path)
        except Exception:
            raise first from None
        sys.stderr.write("  note: dynamic gfortran runtime unresolved; "
                         "relinked %s with a static runtime\n" % os.path.basename(lib_path))
        return handle


def build(contract_path, props_values=None, contract=None, layout=None, verify_path=None):
    if contract is None:
        contract, layout = normalize_contract(json.load(open(contract_path)), contract_path)
    else:
        contract = dict(contract)
    cdir = os.path.dirname(os.path.abspath(contract_path))
    # output filenames come from the directory; `material` is an optional label.
    contract["_name"] = os.path.basename(cdir)
    umat_abs = os.path.join(cdir, contract["source"]["main_file"])
    if not os.path.exists(umat_abs):
        raise ContractError("%s: source.main_file %r does not exist"
                            % (contract_path, contract["source"]["main_file"]))
    dims = contract["dimensions"]
    ntens, nprops, nstatev = int(dims["ntens"]), int(dims["nprops"]), int(dims["nstatev"])
    params = contract["parameters"]
    # parameters/dimensions were already validated in normalize_contract; nparam
    # is inferred here (never read from the contract).
    nparam = len(params)

    # the operating point validation runs at, from the verification case (or
    # --props / --verify), decided BEFORE anything is built
    props_values, props_origin, loading_path = resolve_props(
        contract, nprops, props_values, cdir=cdir, verify_path=verify_path)

    out = os.path.join(cdir, "build"); os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "ABA_PARAM.INC"), "w").write(ABA_PARAM)

    # 1) transform (the real transformer)
    auto = _auto_infer(umat_abs)
    auto.update(contract.get("_auto_override") or {})
    # path-dependence is the derivatives decision (auto -> nstatev>0), resolved
    # in normalize_contract; DSTATEV_DP is chained exactly when it is set.
    path_dependent = (contract.get("_derivatives") or {}).get("path_dependent", nstatev > 0)
    cfg = (_transform_config_path(contract, umat_abs, auto) if path_dependent
           else _transform_config(contract, umat_abs, auto))
    json.dump(cfg, open(os.path.join(out, "transform_config.json"), "w"), indent=2)
    from umat_oti.cli_json import run_config_transform
    from pathlib import Path
    with contextlib.redirect_stdout(io.StringIO()):
        summary, _ = run_config_transform(Path(os.path.join(out, "transform_config.json")), Path(out))
    assert summary.get("transform_success"), "P1-05/10 transform: %s" % summary.get("status_category")
    transformed_for = summary["transformed_source"]

    # 2) rename SUBROUTINE UMAT -> UMAT_OTI_INTERNAL so the original UMAT can co-exist
    txt = open(transformed_for).read()
    txt = re.sub(r"(SUBROUTINE\s+)UMAT(\b)", r"\1UMAT_OTI_INTERNAL\2", txt, count=1)
    internal_for = os.path.join(out, "umat_oti_internal.for")
    open(internal_for, "w").write(txt)

    # 3) generate UMAT_OTI_EVAL (path-dependent eval chains dstatev/dp across increments)
    eval_for = os.path.join(out, "umat_oti_eval.for")
    open(eval_for, "w").write(_gen_eval_path(ntens, nprops, nstatev, nparam) if path_dependent
                              else _gen_eval(ntens, nprops, nstatev, nparam))

    # 3b) EFFICIENT whole-path march (path-dependent only): a second, param-only
    #     transform -> UMAT_OTI_MARCH that carries the OTI state across increments.
    march_for = march_module = None
    if path_dependent:
        mout = os.path.join(out, "march"); os.makedirs(mout, exist_ok=True)
        open(os.path.join(mout, "ABA_PARAM.INC"), "w").write(ABA_PARAM)
        mcfg = _transform_config_march(contract, umat_abs, auto)
        json.dump(mcfg, open(os.path.join(mout, "transform_config.json"), "w"), indent=2)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                msum, _ = run_config_transform(Path(os.path.join(mout, "transform_config.json")), Path(mout))
            if msum.get("transform_success"):
                march_src = _gen_march(open(msum["transformed_source"]).read(), ntens, nstatev, nparam)
                march_for = os.path.join(mout, "umat_oti_march.for")
                open(march_for, "w").write(march_src)
                march_module = [f for f in os.listdir(mout) if re.match(r"otim\d+n1\.f90$", f)][0]
        except Exception:
            march_for = None

    # 4) compile everything and bundle into ONE relocatable .obj (ld -r)
    module = [f for f in os.listdir(out) if re.match(r"otim\d+n1\.f90$", f)][0]
    objs = []
    for src in ("master_parameters.f90", "real_utils.f90", module):
        o = os.path.join(out, src + ".o")
        _gfc(["-c", "-fPIC", "-ffree-form", "-ffree-line-length-none",
              os.path.join(out, src), "-J" + out, "-o", o]); objs.append(o)
    for src in (internal_for, eval_for, umat_abs):
        o = os.path.join(out, os.path.basename(src) + ".o")
        _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + out,
              src, "-J" + out, "-o", o]); objs.append(o)
    if march_for and march_module != module:
        # the march's (smaller) OTI module -- master_parameters/real_utils .mod are
        # reused from -I out so their objects are not duplicated in the bundle.
        mo = os.path.join(mout, march_module + ".o")
        _gfc(["-c", "-fPIC", "-ffree-form", "-ffree-line-length-none", "-I" + out,
              os.path.join(mout, march_module), "-J" + mout, "-o", mo]); objs.append(mo)
    if march_for:
        mo = os.path.join(mout, "umat_oti_march.for.o")
        _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + mout, "-I" + out,
              march_for, "-J" + mout, "-o", mo]); objs.append(mo)
    obj = os.path.join(cdir, "umat_%s_oti.obj" % contract["_name"])
    subprocess.check_call(["ld", "-r", "-o", obj] + objs)

    # 5) completed contract JSON
    contract_ver = json.load(open(os.path.join(CONTRACT_DIR, "CONTRACT_VERSION.json")))["combined_hash"]
    completed = {
        "schema": "resasm_umat_oti_contract_v1",
        "model_id": "umat_%s_oti" % contract["_name"],
        "material": contract.get("_material") or contract["_name"],
        "kinematics": contract["kinematics"],
        # nparam and the derivative array shapes are INFERRED/generated here from
        # len(parameters) and the interface sizes -- they are never user input.
        "dimensions": {"ntens": ntens, "nprops": nprops, "nstatev": nstatev, "nparam": nparam},
        "derived": {
            "nparam": nparam,
            "DSIGMA_DP_shape": [ntens, nparam],
            "DSTATEV_DP_shape": [nstatev, nparam] if path_dependent else None,
            "oti_directions": {p["name"]: k + 1 for k, p in enumerate(params)},
            "derivatives_requested": contract.get("_derivatives", {"stress": True, "statev": "auto"}),
        },
        "symbols": {
            "regular_umat": "umat",
            "oti_internal": "umat_oti_internal",
            "oti_eval": "umat_oti_eval_",
            "oti_eval_signature": (["STRESS(NTENS)", "STATEV(NSTATV)", "DDSDDE(NTENS,NTENS)",
                                   "STRAN(NTENS)", "DSTRAN(NTENS)", "TIME(2)", "DTIME",
                                   "TEMP", "DTEMP", "PROPS(NPROPS)", "NPROPS", "NTENS",
                                   "NSTATV", "NPARAM", "DSIGMA_DP(NTENS,NPARAM)",
                                   "DSTATEV_DP(NSTATV,NPARAM)"]
                                   + (["DSIGMA_DP_IN(NTENS,NPARAM)", "DSTATEV_DP_IN(NSTATV,NPARAM)"]
                                      if nstatev > 0 else [])),
        },
        "replay": {"mode": "path_marching" if nstatev > 0 else "single_evaluation",
                   "carry": (["DSIGMA_DP", "DSTATEV_DP"] if nstatev > 0 else [])},
        "march": ({"symbol": "umat_oti_march_", "directions": ntens + nparam,
                   "signature": ["PROPS(NPROPS)", "NPROPS", "PATH(NTENS,NPATH)", "NPATH",
                                 "DTARR(NPATH)", "NTENS", "NSTATV", "NPARAM",
                                 "DSIG(NTENS,NPARAM,NPATH)", "STROUT(NTENS)",
                                 "DDOUT(NTENS,NTENS,NPATH)"]}
                  if march_for else None),
        "layouts": {"DSIGMA_DP": "fortran(NTENS,NPARAM)", "DSTATEV_DP": "fortran(NSTATV,NPARAM)",
                    "DDSDDE": "fortran(NTENS,NTENS)", "voigt": ["11", "22", "33", "12", "13", "23"]},
        "parameters": [{"name": p["name"], "props_index": int(p["props_index"]),
                        "oti_direction": k + 1} for k, p in enumerate(params)],
        "history": {"path_dependent": nstatev > 0,
                    "dstatev_dp": "returned" if nstatev > 0 else "not_applicable"},
        "contract_version": contract_ver,
        "object": {"file": os.path.basename(obj), "sha256": _sha(obj)},
        "regular_source_hash": _sha(umat_abs),
    }
    # binary metadata: what this compiled object IS and how it was built, so a
    # collaborator on a different platform is rejected before linking (item 2).
    # The contract is platform-independent; the object it names is not.
    completed["binary"] = binary_metadata(
        obj, abi_version=contract_ver, source_hash=_sha(umat_abs),
        transform_hash=summary.get("transform_hash", ""),
        extra={"object_file": os.path.basename(obj), "object_sha256": _sha(obj)})
    cjson = os.path.join(cdir, "umat_%s_oti.json" % contract["_name"])
    json.dump(completed, open(cjson, "w"), indent=2)

    # 6) INDEPENDENT validation (non-circular): compare against the separately
    #    compiled ORIGINAL regular UMAT. Path-dependent materials are validated
    #    over a MULTI-INCREMENT strain path (DSIGMA_DP and DSTATEV_DP).
    val = (_validate_path(out, obj, umat_abs, ntens, nprops, nstatev, params, props_values)
           if path_dependent
           else _validate(out, obj, umat_abs, ntens, nprops, nstatev, params, props_values))

    # provenance: everything needed to reproduce and audit this record (item 4)
    val["method_version"] = VALIDATION_METHOD_VERSION
    val["metric_definitions"] = METRIC_DEFINITIONS
    val["props_origin"] = props_origin
    val["loading_path"] = loading_path
    val["provenance"] = provenance_block(
        repos={"umat_oti": REPO},
        props_values=props_values,
        ladder=list(fdref.DEFAULT_LADDER),
        step_selection="two-sided plateau (max of adjacent successive changes), "
                       "Richardson-checked; finest rung self-corroborates only above the "
                       "roundoff floor 1e-12",
        tolerances=(val.get("fd") or {}).get("tolerances", {}),
        metric_defs={k: v["formula"] for k, v in METRIC_DEFINITIONS.items()},
        status="PASS" if val["passed"] else "FAIL",
        reason="; ".join(val.get("failures", [])) or "all checks within tolerance",
        extra={"props_origin": props_origin, "loading_path": loading_path,
               "input_format": layout},
    )
    completed["validation"] = val
    json.dump(completed, open(cjson, "w"), indent=2)

    def _num(x):
        return "%.2e" % x if isinstance(x, float) and math.isfinite(x) else "NON-FINITE"

    print("=" * 74)
    print(" Program 1  umat_%s_oti  (object: %s)" % (contract["_name"], os.path.basename(obj)))
    print("=" * 74)
    print(" contract / format           : %s  (%s)" % (os.path.basename(contract_path), layout))
    print(" material / nparam           : %s / %d (inferred from parameters)"
          % (contract.get("_material") or contract["_name"], nparam))
    print(" validated at PROPS          : [%s]" % ", ".join("%g" % v for v in props_values))
    print("   source of those values    : %s" % props_origin)
    print(" symbols in .obj             : UMAT, UMAT_OTI_INTERNAL, UMAT_OTI_EVAL")
    print(" transformer semantic checks : %s" % ("all pass" if all(summary.get("semantic_checks", {}).values()) else "SEE REPORT"))
    if path_dependent:
        print(" path-dependent validation   : %d increments (elastic->plastic->unload->reload)" % val["n_increments"])
        print(" stress parity   vs orig     : %s" % _num(val["stress_parity_max_rel"]))
        print(" STATEV parity   vs orig     : %s" % _num(val["statev_parity_max_rel"]))
        print(" DDSDDE vs FD(orig dstran)   : %s" % _num(val["ddsdde_parity_max_rel"]))
        print(" DDSDDE vs regular UMAT      : %s  (elastic-tangent UMAT -> expected)" % _num(val["ddsdde_vs_regular_rel"]))
        print(" DSIGMA_DP  vs converged FD  : %s" % _num(val["dsigma_max_rel"]))
        print(" DSTATEV_DP vs converged FD  : %s" % _num(val["dstatev_max_rel"]))
    else:
        print(" stress real parity  vs orig : %s" % _num(val["stress_parity_max_rel"]))
        print(" DDSDDE parity       vs orig : %s" % _num(val["ddsdde_parity_max_rel"]))
        print(" DSIGMA_DP vs converged FD   : %s" % _num(val["dsigma_max_rel"]))
    print(" per-parameter FD evidence   :")
    print(fdref.format_report(val["fd"], indent="   "))
    for reason in val["failures"]:
        print("   FAIL: %s" % reason)
    print(" -> %s" % ("PASS" if val["passed"] else "FAIL"))
    return 0 if val["passed"] else 1


# --------------------------------------------------------------------------- #
# validation drivers: build one .so around UMAT_OTI_EVAL and one around the
# ORIGINAL UMAT, evaluate both at material points, FD the ORIGINAL.
# --------------------------------------------------------------------------- #
_ORIG_DRV = r"""
      subroutine orig_eval(props, np, dstrain, nt, stress, ddsdde) bind(C, name="orig_eval")
      use iso_c_binding
      real(c_double), intent(in) :: props(*), dstrain(*)
      integer(c_int), value :: np, nt
      real(c_double), intent(out) :: stress(*), ddsdde(*)
      real*8 :: STRESS_(nt), STATEV_(1), DDSDDE_(nt,nt), STRAN(nt), DSTRAN(nt)
      real*8 :: TIME(2), PREDEF(1), DPRED(1), PROPS_(np), COORDS(3), DROT(3,3)
      real*8 :: DFGRD0(3,3), DFGRD1(3,3), DDSDDT(nt), DRPLDE(nt)
      character*80 CMNAME
      integer i,j
      do i=1,nt; STRESS_(i)=0.d0; STRAN(i)=0.d0; DSTRAN(i)=dstrain(i); end do
      do i=1,np; PROPS_(i)=props(i); end do
      do i=1,3; do j=1,3; DROT(i,j)=0.d0; DFGRD0(i,j)=0.d0; DFGRD1(i,j)=0.d0; end do
        DROT(i,i)=1.d0; DFGRD0(i,i)=1.d0; DFGRD1(i,i)=1.d0; end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT(STRESS_,STATEV_,DDSDDE_,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     1 STRAN,DSTRAN,TIME,1.d0,0.d0,0.d0,PREDEF,DPRED,CMNAME,
     2 3,nt-3,nt,1,PROPS_,np,COORDS,DROT,1.d0,1.d0,DFGRD0,DFGRD1,1,1,1,1,1,1)
      do i=1,nt; stress(i)=STRESS_(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      end subroutine
"""

_OTI_DRV = r"""
      subroutine oti_eval(props, np, dstrain, nt, ns, nprm, stress, ddsdde, dsig, dstv) bind(C, name="oti_eval")
      use iso_c_binding
      real(c_double), intent(in) :: props(*), dstrain(*)
      integer(c_int), value :: np, nt, ns, nprm
      real(c_double), intent(out) :: stress(*), ddsdde(*), dsig(*), dstv(*)
      real*8 :: STRESS_(nt), STATEV_(max(ns,1)), DDSDDE_(nt,nt), STRAN(nt), DSTRAN(nt)
      real*8 :: TIME(2), PROPS_(np), DSIGMA_DP(nt,nprm), DSTATEV_DP(max(ns,1),nprm)
      integer i,j
      do i=1,nt; STRESS_(i)=0.d0; STRAN(i)=0.d0; DSTRAN(i)=dstrain(i); end do
      do i=1,max(ns,1); STATEV_(i)=0.d0; end do
      do i=1,np; PROPS_(i)=props(i); end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT_OTI_EVAL(STRESS_,STATEV_,DDSDDE_,STRAN,DSTRAN,TIME,1.d0,0.d0,0.d0,
     1 PROPS_,np,nt,ns,nprm,DSIGMA_DP,DSTATEV_DP)
      do i=1,nt; stress(i)=STRESS_(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      do j=1,nprm; do i=1,nt; dsig((i-1)*nprm+j)=DSIGMA_DP(i,j); end do; end do
      do j=1,nprm; do i=1,ns; dstv((i-1)*nprm+j)=DSTATEV_DP(i,j); end do; end do
      end subroutine
"""


# marching drivers for path-dependent validation: both accept incoming STRESS +
# STATEV and return the updated ones, so a strain PATH can be marched.
_ORIG_DRV_PATH = r"""
      subroutine orig_path(props,np,strin,statin,ns,dstran,nt,
     1 strout,ddsdde,statout) bind(C,name="orig_path")
      use iso_c_binding
      real(c_double),intent(in)::props(*),strin(*),statin(*),dstran(*)
      integer(c_int),value::np,ns,nt
      real(c_double),intent(out)::strout(*),ddsdde(*),statout(*)
      real*8 STRESS(nt),STATEV(max(ns,1)),DDSDDE_(nt,nt),STRAN(nt),DSTRN(nt)
      real*8 TIME(2),PREDEF(1),DPRED(1),PROPS_(np),COORDS(3),DROT(3,3)
      real*8 DFGRD0(3,3),DFGRD1(3,3),DDSDDT(nt),DRPLDE(nt)
      real*8 SSE,SPD,SCD,RPL,DRPLDT
      character*80 CMNAME
      integer i,j
      SSE=0.d0; SPD=0.d0; SCD=0.d0; RPL=0.d0; DRPLDT=0.d0
      do i=1,nt; STRESS(i)=strin(i); STRAN(i)=0.d0; DSTRN(i)=dstran(i)
        DDSDDT(i)=0.d0; DRPLDE(i)=0.d0; end do
      do i=1,max(ns,1); STATEV(i)=0.d0; end do
      do i=1,ns; STATEV(i)=statin(i); end do
      do i=1,np; PROPS_(i)=props(i); end do
      do i=1,3; COORDS(i)=0.d0; do j=1,3; DROT(i,j)=0.d0
        DFGRD0(i,j)=0.d0; DFGRD1(i,j)=0.d0; end do
        DROT(i,i)=1.d0; DFGRD0(i,i)=1.d0; DFGRD1(i,i)=1.d0; end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT(STRESS,STATEV,DDSDDE_,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     1 STRAN,DSTRN,TIME,1.d0,0.d0,0.d0,PREDEF,DPRED,CMNAME,
     2 3,nt-3,nt,max(ns,1),PROPS_,np,COORDS,DROT,1.d0,1.d0,DFGRD0,DFGRD1,
     3 1,1,1,1,1,1)
      do i=1,nt; strout(i)=STRESS(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      do i=1,ns; statout(i)=STATEV(i); end do
      end subroutine
"""

_OTI_DRV_PATH = r"""
      subroutine oti_path(props,np,strin,statin,ns,dsigin,dstvin,
     1 dstran,nt,nprm,strout,ddsdde,statout,dsigout,dstvout)
     2 bind(C,name="oti_path")
      use iso_c_binding
      real(c_double),intent(in)::props(*),strin(*),statin(*),dstran(*)
      real(c_double),intent(in)::dsigin(*),dstvin(*)
      integer(c_int),value::np,ns,nt,nprm
      real(c_double),intent(out)::strout(*),ddsdde(*),statout(*),dsigout(*),dstvout(*)
      real*8 STRESS(nt),STATEV(max(ns,1)),DDSDDE_(nt,nt),STRAN(nt),DSTRN(nt)
      real*8 TIME(2),PROPS_(np),DSIGMA_DP(nt,nprm),DSTATEV_DP(max(ns,1),nprm)
      real*8 DSIGMA_DP_IN(nt,nprm),DSTATEV_DP_IN(max(ns,1),nprm)
      integer i,j
      do i=1,nt; STRESS(i)=strin(i); STRAN(i)=0.d0; DSTRN(i)=dstran(i); end do
      do i=1,max(ns,1); STATEV(i)=0.d0; end do
      do i=1,ns; STATEV(i)=statin(i); end do
      do i=1,np; PROPS_(i)=props(i); end do
      do j=1,nprm; do i=1,nt; DSIGMA_DP_IN(i,j)=dsigin((i-1)*nprm+j); end do; end do
      do j=1,nprm; do i=1,ns; DSTATEV_DP_IN(i,j)=dstvin((i-1)*nprm+j); end do; end do
      TIME(1)=0.d0; TIME(2)=0.d0
      call UMAT_OTI_EVAL(STRESS,STATEV,DDSDDE_,STRAN,DSTRN,TIME,1.d0,0.d0,0.d0,
     1 PROPS_,np,nt,ns,nprm,DSIGMA_DP,DSTATEV_DP,DSIGMA_DP_IN,DSTATEV_DP_IN)
      do i=1,nt; strout(i)=STRESS(i); end do
      do i=1,nt; do j=1,nt; ddsdde((i-1)*nt+j)=DDSDDE_(i,j); end do; end do
      do i=1,ns; statout(i)=STATEV(i); end do
      do j=1,nprm; do i=1,nt; dsigout((i-1)*nprm+j)=DSIGMA_DP(i,j); end do; end do
      do j=1,nprm; do i=1,ns; dstvout((i-1)*nprm+j)=DSTATEV_DP(i,j); end do; end do
      end subroutine
"""


def _validate_path(out, obj, umat_abs, ntens, nprops, nstatev, params, props_values):
    """Multi-increment material-point validation for a path-dependent material.

    March the OTI-chained eval over a strain PATH (elastic -> plastic -> partial
    unload -> reload) and compare TOTAL DSIGMA_DP and DSTATEV_DP at every
    increment against the canonical converged finite-difference reference
    (:mod:`umat_oti.validation.fd_reference`) built from the separately compiled
    ORIGINAL UMAT marched over the SAME path.  Non-circular: FD perturbs PROPS
    of the ORIGINAL, and the FD step is chosen from the FD sequence alone.
    """
    od = os.path.join(out, "orig_path.for"); open(od, "w").write(_ORIG_DRV_PATH)
    _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + out, od, "-o", od + ".o"])
    orig_umat_o = os.path.join(out, os.path.basename(umat_abs) + ".o")
    lo = _link_and_load([od + ".o", orig_umat_o], _validation_lib(out, "orig"))
    td = os.path.join(out, "oti_path.for"); open(td, "w").write(_OTI_DRV_PATH)
    _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + out, td, "-o", td + ".o"])
    lt = _link_and_load([td + ".o", obj], _validation_lib(out, "oti"))

    D = ctypes.c_double; IP = ctypes.POINTER(D)
    lo.lib.orig_path.argtypes = [IP, ctypes.c_int, IP, IP, ctypes.c_int, IP, ctypes.c_int, IP, IP, IP]
    lt.lib.oti_path.argtypes = [IP, ctypes.c_int, IP, IP, ctypes.c_int, IP, IP, IP,
                                ctypes.c_int, ctypes.c_int, IP, IP, IP, IP, IP]
    nprm = len(params)
    p0 = list(props_values)
    ns = nstatev

    def orig_step(props, strin, statin, dstran):
        pr = (D * nprops)(*props); si = (D * ntens)(*strin)
        st = (D * max(ns, 1))(*(list(statin) + [0.0] * (max(ns, 1) - ns))); de = (D * ntens)(*dstran)
        so = (D * ntens)(); dd = (D * (ntens * ntens))(); sto = (D * max(ns, 1))()
        lo.lib.orig_path(pr, nprops, si, st, ns, de, ntens, so, dd, sto)
        return np.array(so), np.array(sto)[:ns], np.array(dd).reshape(ntens, ntens)

    def oti_step(props, strin, statin, dsigin, dstvin, dstran):
        pr = (D * nprops)(*props); si = (D * ntens)(*strin)
        st = (D * max(ns, 1))(*(list(statin) + [0.0] * (max(ns, 1) - ns))); de = (D * ntens)(*dstran)
        dsi = (D * (ntens * nprm))(*dsigin.flatten()); dsv = (D * (ns * nprm))(*dstvin.flatten())
        so = (D * ntens)(); dd = (D * (ntens * ntens))(); sto = (D * max(ns, 1))()
        dso = (D * (ntens * nprm))(); dsvo = (D * (ns * nprm))()
        lt.lib.oti_path(pr, nprops, si, st, ns, dsi, dsv, de, ntens, nprm, so, dd, sto, dso, dsvo)
        return (np.array(so), np.array(sto)[:ns], np.array(dd).reshape(ntens, ntens),
                np.array(dso).reshape(ntens, nprm), np.array(dsvo).reshape(ns, nprm))

    full = [1.0e-3, -3.0e-4, -3.0e-4, 2.0e-4, 1.0e-4, -5.0e-5]
    unld = [-6.0e-4, 2.0e-4, 2.0e-4, -1.0e-4, 0.0, 0.0]
    dpath = [full, [x * 0.7 for x in full], unld, full, [x * 0.5 for x in full]]
    dpath = [row[:ntens] for row in dpath]

    # OTI-chained march (once) and ORIGINAL base march (once)
    oti_dsig, oti_dstv, oti_str, base_str, base_stv, base_dd, oti_dd = [], [], [], [], [], [], []
    stro = np.zeros(ntens); stato = np.zeros(ns); dsig = np.zeros((ntens, nprm)); dstv = np.zeros((ns, nprm))
    sb = np.zeros(ntens); stb = np.zeros(ns)
    for d in dpath:
        stro, stato, ddo, dsig, dstv = oti_step(p0, stro, stato, dsig, dstv, d)
        sb, stb, ddb = orig_step(p0, sb, stb, d)
        oti_dsig.append(dsig.copy()); oti_dstv.append(dstv.copy()); oti_str.append(stro.copy())
        oti_dd.append(ddo.copy()); base_str.append(sb.copy()); base_stv.append(stb.copy()); base_dd.append(ddb.copy())
    # every primal quantity must be finite before anything is compared
    failures = []
    try:
        fdref.require_finite("original STRESS along the path", base_str)
        fdref.require_finite("original STATEV along the path", base_stv)
        fdref.require_finite("OTI STRESS along the path", oti_str)
        fdref.require_finite("OTI DSIGMA_DP along the path", oti_dsig)
        fdref.require_finite("OTI DSTATEV_DP along the path", oti_dstv)
    except fdref.NonFiniteResult as exc:
        failures.append(str(exc))

    sp_rel = max(fdref.relative_error(oti_str[n], base_str[n]) for n in range(len(dpath)))
    # physical STATEV parity: re-march OTI and compare updated state to the original
    sv_rel = 0.0
    strc = np.zeros(ntens); stc = np.zeros(ns); dsg = np.zeros((ntens, nprm)); dtv = np.zeros((ns, nprm))
    for n, d in enumerate(dpath):
        strc, stc, _, dsg, dtv = oti_step(p0, strc, stc, dsg, dtv, d)
        sv_rel = max(sv_rel, fdref.relative_error(stc, base_stv[n]))
    # informational: OTI consistent tangent vs the regular UMAT's returned DDSDDE
    # (may be elastic-only for viscoplastic models -> a documented difference).
    dd_rel_reg = max(fdref.relative_error(oti_dd[n], base_dd[n]) for n in range(len(dpath)))

    # tangent check: OTI DDSDDE vs the ORIGINAL stress differenced wrt DSTRAN at
    # the last (plastic) increment, incoming state held fixed.  Same canonical
    # ladder as the parameter derivatives, scaled by the strain increment itself
    # so components that are legitimately zero are not perturbed by a full unit.
    nlast = len(dpath) - 1
    s_in = base_str[nlast - 1] if nlast >= 1 else np.zeros(ntens)
    st_in = base_stv[nlast - 1] if nlast >= 1 else np.zeros(ns)
    dlast = np.array(dpath[nlast], dtype=float)
    dscale = float(np.max(np.abs(dlast))) or 1.0

    def _tangent_response(dvec):
        s, _, _ = orig_step(p0, s_in, st_in, list(dvec))
        return {"stress": fdref.require_finite("original STRESS (tangent probe)", s)}

    fd_T = np.zeros((ntens, ntens))
    tangent_rows = []
    for j in range(ntens):
        st_j = fdref.parameter_studies(_tangent_response, dlast, j + 1,
                                       parameter="DSTRAN%d" % (j + 1), scale=dscale)["stress"]
        fdref.apply_convergence_gate(st_j, 1e-5)
        if st_j.reference is None:
            failures.append("tangent FD for DSTRAN(%d): %s" % (j + 1, st_j.note))
            fd_T[:, j] = np.nan
        else:
            fd_T[:, j] = st_j.reference
        tangent_rows.append(st_j.to_json())
    dd_rel = fdref.relative_error(oti_dd[nlast], fd_T)

    # DSIGMA_DP / DSTATEV_DP vs the canonical converged FD reference
    def _path_response(props):
        s = np.zeros(ntens); st = np.zeros(ns)
        S, V = [], []
        for d in dpath:
            s, st, _ = orig_step(props, s, st, d)
            S.append(s.copy()); V.append(st.copy())
        return {"stress": fdref.require_finite("original STRESS (FD probe)", S),
                "state": fdref.require_finite("original STATEV (FD probe)", V)}

    fd_report = fdref.verify_derivatives(
        _path_response, p0, params,
        {"stress": np.array(oti_dsig), "state": np.array(oti_dstv)},
        tolerances={"stress": 1e-4, "state": 1e-4},
    )
    fd_report["tangent"] = tangent_rows
    best_dsig = fd_report["worst_rel"]["stress"]
    best_dstv = fd_report["worst_rel"]["state"]

    scalars_ok = all(math.isfinite(v) for v in (sp_rel, sv_rel, dd_rel, dd_rel_reg))
    if not scalars_ok:
        failures.append("a parity/tangent measure evaluated to a non-finite number")
    if sp_rel > 1e-8:
        failures.append("stress parity %.2e exceeds 1e-8" % sp_rel)
    if sv_rel > 1e-8:
        failures.append("STATEV parity %.2e exceeds 1e-8" % sv_rel)
    if dd_rel > 1e-4:
        failures.append("DDSDDE vs FD %.2e exceeds 1e-4" % dd_rel)
    failures += [r["reason"] for r in fd_report["parameters"] if r["status"] != "PASS"]

    return {"stress_parity_max_rel": sp_rel, "statev_parity_max_rel": sv_rel,
            "ddsdde_parity_max_rel": dd_rel, "ddsdde_vs_regular_rel": dd_rel_reg,
            "dsigma_max_rel": best_dsig, "dstatev_max_rel": best_dstv,
            "alldir_vs_single": 0.0, "n_increments": len(dpath), "props_values": p0,
            "path_dependent": True, "fd": fd_report, "failures": failures,
            "passed": bool(not failures and scalars_ok and fd_report["passed"])}


def _validate(out, obj, umat_abs, ntens, nprops, nstatev, params, props_values):
    # build original-UMAT .so
    od = os.path.join(out, "orig_drv.for"); open(od, "w").write(_ORIG_DRV)
    _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + out, od, "-o", od + ".o"])
    orig_umat_o = os.path.join(out, os.path.basename(umat_abs) + ".o")   # already compiled
    lo = _link_and_load([od + ".o", orig_umat_o], _validation_lib(out, "orig"))
    # build the OTI-eval library (links the bundled .obj)
    td = os.path.join(out, "oti_drv.for"); open(td, "w").write(_OTI_DRV)
    _gfc(["-c", "-fPIC", "-ffixed-form", "-ffixed-line-length-none", "-I" + out, td, "-o", td + ".o"])
    lt = _link_and_load([td + ".o", obj], _validation_lib(out, "oti"))

    D = ctypes.c_double; IP = ctypes.POINTER(D)
    lo.lib.orig_eval.argtypes = [IP, ctypes.c_int, IP, ctypes.c_int, IP, IP]
    lt.lib.oti_eval.argtypes = [IP, ctypes.c_int, IP, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                IP, IP, IP, IP]
    nprm = len(params)
    p0 = list(props_values)

    def orig(props, eps):
        pr = (D * nprops)(*props); ep = (D * ntens)(*eps)
        sg = (D * ntens)(); dd = (D * (ntens * ntens))()
        lo.lib.orig_eval(pr, nprops, ep, ntens, sg, dd)
        return np.array(sg), np.array(dd).reshape(ntens, ntens)

    def oti(props, eps):
        pr = (D * nprops)(*props); ep = (D * ntens)(*eps)
        sg = (D * ntens)(); dd = (D * (ntens * ntens))()
        dsi = (D * (ntens * nprm))(); dsv = (D * (max(nstatev, 1) * nprm))()
        lt.lib.oti_eval(pr, nprops, ep, ntens, nstatev, nprm, sg, dd, dsi, dsv)
        return np.array(sg), np.array(dd).reshape(ntens, ntens), np.array(dsi).reshape(ntens, nprm)

    eps = [1e-3, 2e-4, -3e-4, 1e-4, 5e-5, -2e-5][:ntens]
    so, do = orig(p0, eps)
    st, dt, dsig = oti(p0, eps)

    failures = []
    try:
        fdref.require_finite("original STRESS", so)
        fdref.require_finite("original DDSDDE", do)
        fdref.require_finite("OTI STRESS", st)
        fdref.require_finite("OTI DDSDDE", dt)
        fdref.require_finite("OTI DSIGMA_DP", dsig)
    except fdref.NonFiniteResult as exc:
        failures.append(str(exc))

    sp_rel = fdref.relative_error(st, so)
    dd_rel = fdref.relative_error(dt, do)

    def _response(props):
        s, _ = orig(props, eps)
        return {"stress": fdref.require_finite("original STRESS (FD probe)", s)}

    fd_report = fdref.verify_derivatives(_response, p0, params, {"stress": dsig},
                                         tolerances={"stress": 1e-5})
    best = fd_report["worst_rel"]["stress"]

    scalars_ok = math.isfinite(sp_rel) and math.isfinite(dd_rel)
    if not scalars_ok:
        failures.append("a parity measure evaluated to a non-finite number")
    if sp_rel > 1e-8:
        failures.append("stress parity %.2e exceeds 1e-8" % sp_rel)
    if dd_rel > 1e-8:
        failures.append("DDSDDE parity %.2e exceeds 1e-8" % dd_rel)
    failures += [r["reason"] for r in fd_report["parameters"] if r["status"] != "PASS"]

    return {"stress_parity_max_rel": sp_rel, "ddsdde_parity_max_rel": dd_rel,
            "dsigma_max_rel": best, "alldir_vs_single": 0.0,
            "props_values": p0, "path_dependent": False,
            "fd": fd_report, "failures": failures,
            "passed": bool(not failures and scalars_ok and fd_report["passed"])}


def _main(argv):
    if not argv or argv[0] not in ("build", "validate-all"):
        print(__doc__); return 2

    if argv[0] == "build":
        props = None
        if "--props" in argv:
            props = [float(x) for x in argv[argv.index("--props") + 1].split(",")]
        verify = argv[argv.index("--verify") + 1] if "--verify" in argv else None
        try:
            return build(argv[1], props, verify_path=verify)
        except (ContractError, LibraryLoadError) as exc:
            print("ERROR: %s" % exc)
            return 1

    # validate-all: every material with a transformation contract (any format)
    mats = os.path.join(HERE, "materials")
    entries, skipped = discover_contracts(mats)
    print("=" * 74)
    print(" contract discovery under %s" % os.path.relpath(mats, REPO))
    print("=" * 74)
    print(" discovered : %d" % len(entries))
    for e in entries:
        print("     %-24s %-28s (%s)" % (e["material"], os.path.basename(e["path"]), e["format"]))
    if skipped:
        print(" skipped    : %d" % len(skipped))
        for s in skipped:
            print("     %-24s %s" % (s["material"], s["reason"]))
    print()

    executed, failed, errored = [], [], []
    for e in entries:
        try:
            rc = build(e["path"], None, contract=e["contract"], layout=e["format"])
        except (ContractError, LibraryLoadError) as exc:
            print("=" * 74)
            print(" Program 1  %s" % e["material"])
            print("=" * 74)
            print(" ERROR: %s" % exc)
            print(" -> FAIL")
            errored.append((e["material"], str(exc).splitlines()[0]))
            continue
        executed.append(e["material"])
        if rc:
            failed.append(e["material"])

    print("=" * 74)
    print(" validate-all summary")
    print("=" * 74)
    print(" discovered : %d" % (len(entries) + len(skipped)))
    print(" executed   : %d  (%s)" % (len(executed), ", ".join(executed) or "-"))
    print(" passed     : %d" % (len(executed) - len(failed)))
    print(" failed     : %d  (%s)" % (len(failed), ", ".join(failed) or "-"))
    print(" errored    : %d  (%s)" % (len(errored), ", ".join(m for m, _ in errored) or "-"))
    for material, reason in errored:
        print("     %-24s %s" % (material, reason))
    print(" skipped    : %d  (%s)" % (len(skipped), ", ".join(s["material"] for s in skipped) or "-"))
    return 1 if (failed or errored) else 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
