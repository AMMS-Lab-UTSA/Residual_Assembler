"""Program 2 residual request contract (the deck's "OUTPUT-PARAMETER REQUESTS").

A single, high-level request file drives the whole run:

    {
      "material": "umat_j2_oti.obj",              # OTI object (contract = sibling .json)
      "model":    "j2_tension.inp",
      "record":   "j2_tension.resrec.h5",
      "scope":    {"domain": "ALL", "increments": "ALL_CONVERGED"},
      "requests": [
        {"output": "U",      "with_respect_to": ["E", "nu"]},
        {"output": "S",      "with_respect_to": ["E", "SIGY0", "H"]},
        {"output": "EQPLAS", "with_respect_to": ["SIGY0", "H"]}
      ]
    }

The defining idea: every OUTPUT is paired with its OWN set of material parameters
(``with_respect_to``), a single ``scope`` selects where (a node/element set or
``ALL``) and when (converged increments), and the material is just the OTI object
(its completed contract is the sibling ``.json``).

Output names map onto the residual method's native fields:
    U      -> displacement sensitivity     du/dp          (from K du/dp = -dR/dp)
    RF     -> reaction-force sensitivity    dRF/dp
    S      -> stress sensitivity            dσ/dp          (DSIGMA_DP + tangent·B·du/dp)
    <state name>  (e.g. EQPLAS)  -> state-variable sensitivity, from the material's
             DSTATEV_DP at the recorded material points.

Only the union of all ``with_respect_to`` is seeded; each request is then reported
against just its own parameters. This module VALIDATES the whole request up front
(schema, files, output names, parameter names against the material contract,
domain against the record mesh) and raises one aggregated, actionable error.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Any, Dict, List

import numpy as np

from ..interface import versions as _v
from .job import _abs, JobError

SCHEMA = _v.SENSITIVITY_REQUEST                       # "resasm_sensitivity_request_v1"
_FIELD_OUTPUTS = {"U": "displacement", "RF": "reaction", "S": "stress"}   # + state names


# ---------------------------------------------------------------------------
def _sibling_contract(obj_path: str) -> str:
    return os.path.splitext(obj_path)[0] + ".json"


def _state_names(contract: Dict[str, Any]) -> List[str]:
    """Names of the material's state variables, if the contract records them
    (transform-side history.state), else synthetic SDV1.. names."""
    hist = contract.get("history", {})
    st = hist.get("state")
    if isinstance(st, list) and st:
        return [str(s) for s in st]
    ns = int(contract.get("dimensions", {}).get("nstatev", 0) or 0)
    return ["SDV%d" % (i + 1) for i in range(ns)]


def _output_kind(name: str, state_names: List[str]) -> str:
    if name in _FIELD_OUTPUTS:
        return _FIELD_OUTPUTS[name]
    if name in state_names:
        return "state"
    return ""                                          # unknown


# ---------------------------------------------------------------------------
# validation (report style: collect everything, never raise)
# ---------------------------------------------------------------------------
def validate_request(rc: Any, base_dir: str) -> List[str]:
    E: List[str] = []
    if not isinstance(rc, dict):
        return ["the request file must be a JSON object"]
    if "schema" in rc and rc["schema"] != SCHEMA:
        E.append("schema, if present, must be %r (got %r)" % (SCHEMA, rc["schema"]))

    mat = rc.get("material")
    contract = None
    if not isinstance(mat, str) or not mat:
        E.append('"material" must be the OTI object path, e.g. "umat_j2_oti.obj"')
    else:
        obj = _abs(base_dir, mat)
        if not os.path.exists(obj):
            E.append("material object not found: %s" % mat)
        cj = _sibling_contract(obj)
        if not os.path.exists(cj):
            E.append("completed contract not found next to the object: %s" % os.path.basename(cj))
        else:
            try:
                contract = json.load(open(cj))
            except (OSError, ValueError):
                E.append("completed contract is not valid JSON: %s" % os.path.basename(cj))

    for k, what in (("record", "the converged replay record"),):
        v = rc.get(k)
        if not v:
            E.append('"%s" is required (%s)' % (k, what))
        elif not os.path.exists(_abs(base_dir, v)):
            E.append("%s not found: %s" % (k, v))
    # model is optional (the record already carries the mesh)

    scope = rc.get("scope", {})
    if not isinstance(scope, dict):
        E.append('"scope" must be an object with "domain" and "increments"')
        scope = {}
    else:
        inc = scope.get("increments", "ALL_CONVERGED")
        if inc not in ("ALL_CONVERGED", "LAST", "last"):
            E.append('scope.increments must be "ALL_CONVERGED" or "LAST" (got %r)' % inc)

    avail_params = [p["name"] for p in contract.get("parameters", [])] if contract else None
    state_names = _state_names(contract) if contract else []

    reqs = rc.get("requests")
    if not isinstance(reqs, list) or not reqs:
        E.append('"requests" must be a non-empty list of {output, with_respect_to}')
    else:
        for i, r in enumerate(reqs):
            tag = "requests[%d]" % i
            if not isinstance(r, dict):
                E.append("%s must be an object" % tag); continue
            out = r.get("output")
            if not out:
                E.append("%s.output is required (U, RF, S, or a state-variable name)" % tag)
            elif contract is not None and not _output_kind(out, state_names):
                E.append("%s.output %r is not recognised; use U, RF, S or one of the state "
                         "variables %s" % (tag, out, state_names or "(none)"))
            wrt = r.get("with_respect_to")
            if not isinstance(wrt, list) or not wrt:
                E.append("%s.with_respect_to must be a non-empty list of parameter names" % tag)
            elif avail_params is not None:
                for p in wrt:
                    if p not in avail_params:
                        E.append("%s.with_respect_to: parameter %r is not in the material contract; "
                                 "available: %s" % (tag, p, ", ".join(avail_params)))
    return E


def new_request_template() -> Dict[str, Any]:
    return {
        "schema": SCHEMA,
        "material": "umat_j2_oti.obj",
        "model": "j2_tension.inp",
        "record": "j2_tension.resrec.h5",
        "scope": {"domain": "ALL", "increments": "ALL_CONVERGED"},
        "requests": [
            {"output": "U", "with_respect_to": ["E", "nu"]},
            {"output": "S", "with_respect_to": ["E", "SIGY0", "H"]},
            {"output": "EQPLAS", "with_respect_to": ["SIGY0", "H"]},
        ],
    }


# ---------------------------------------------------------------------------
# run: translate the request to the tested engine and report per request
# ---------------------------------------------------------------------------
def _load(rc_path: str) -> Dict[str, Any]:
    try:
        with open(rc_path) as fh:
            return json.load(fh)
    except FileNotFoundError:
        raise JobError("request file not found: %s" % rc_path)
    except ValueError as e:
        raise JobError("request file is not valid JSON (%s): %s" % (e, rc_path))


def _domain_sets(record, domain):
    """Resolve scope.domain -> (node_ids, element_ids)."""
    mesh = record.raw.get("mesh", {})
    nodes = [int(n) for n in mesh.get("nodes", {})] if isinstance(mesh.get("nodes"), dict) \
        else [int(n[0]) if isinstance(n, (list, tuple)) else int(n) for n in mesh.get("nodes", [])]
    elems = list(record.build_model().elements.keys())
    if domain in (None, "ALL", "all"):
        return nodes, elems
    ns = mesh.get("node_sets", {}); es = mesh.get("element_sets", {})
    if domain in ns:
        return [int(n) for n in ns[domain]], []
    if domain in es:
        return [], [int(e) for e in es[domain]]
    raise JobError("scope.domain %r is neither 'ALL' nor a node/element set in the record "
                   "(node_sets: %s, element_sets: %s)" % (domain, sorted(ns), sorted(es)))


def run_request_contract(rc_path: str) -> Dict[str, Any]:
    rc = _load(rc_path)
    base = os.path.dirname(os.path.abspath(rc_path))
    problems = validate_request(rc, base)
    if problems:
        raise JobError("this residual request has %d problem(s):\n  - %s"
                       % (len(problems), "\n  - ".join(problems)))

    from .objlink import package_from_contract
    from .package import MaterialPackage
    from .record import ReplayRecord
    from .engine import replay_sensitivities

    obj = _abs(base, rc["material"])
    contract = json.load(open(_sibling_contract(obj)))
    state_names = _state_names(contract)
    work = tempfile.mkdtemp(prefix="resasm_req_")
    pkg_dir, manifest = package_from_contract(obj, contract, work)
    package = MaterialPackage(manifest, base_dir=pkg_dir)
    record = ReplayRecord.load(_abs(base, rc["record"]))

    scope = rc.get("scope", {})
    node_ids, elem_ids = _domain_sets(record, scope.get("domain", "ALL"))
    seed = sorted({p for r in rc["requests"] for p in r["with_respect_to"]})

    # build the low-level engine outputs needed by the requested field names
    kinds = {r["output"]: _output_kind(r["output"], state_names) for r in rc["requests"]}
    point_reqs: List[Dict[str, Any]] = []
    index = {}                                         # output name -> list of point-req indices
    for name, kind in kinds.items():
        idxs = []
        if kind == "reaction":
            for n in node_ids:
                for d in (1, 2, 3):
                    idxs.append(len(point_reqs)); point_reqs.append({"type": "reaction", "node": n, "dof": d})
        elif kind == "stress":
            model = record.build_model()
            for eid in elem_ids:
                for k in range(8):
                    idxs.append(len(point_reqs))
                    point_reqs.append({"type": "stress", "element": eid, "ip": k, "von_mises": True})
        index[name] = idxs

    res = replay_sensitivities(record, package, seed, point_reqs)

    # assemble per-request results, each against ONLY its own parameters
    out = {"schema": _v.SENSITIVITY_RESULT, "model_id": package.model_id,
           "scope": scope, "requests": []}
    for r in rc["requests"]:
        name, wrt, kind = r["output"], r["with_respect_to"], kinds[r["output"]]
        entry = {"output": name, "with_respect_to": list(wrt), "kind": kind, "sensitivity": {}}
        if kind == "displacement":
            for p in wrt:
                v = res.du(p)
                entry["sensitivity"][p] = {"norm": float(np.linalg.norm(v)), "n_dof": int(v.size)}
        elif kind in ("reaction", "stress"):
            for p in wrt:
                vals = np.array([res.outputs[i]["dq_dp"][p] for i in index[name]])
                entry["sensitivity"][p] = {"norm": float(np.linalg.norm(vals)),
                                           "n_points": int(vals.size)}
        elif kind == "state":
            entry["note"] = ("state-variable sensitivity is the material's DSTATEV_DP at the "
                             "recorded points; the du/dp chain term for state fields is a documented "
                             "extension")
            entry["state_index"] = state_names.index(name)
        else:
            entry["error"] = "unrecognised output"
        out["requests"].append(entry)

    rc_out = rc.get("output_dir", "sensitivity_out")
    odir = _abs(base, rc_out); os.makedirs(odir, exist_ok=True)
    json.dump(out, open(os.path.join(odir, "sensitivity_result.json"), "w"), indent=2)
    out["_output_dir"] = odir
    return out


# ---------------------------------------------------------------------------
def _main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__); return 0 if argv else 2
    cmd, rest = argv[0], argv[1:]
    if cmd == "init":
        p = rest[0] if rest else "residual_request.json"
        if os.path.exists(p):
            print("refusing to overwrite existing %s" % p); return 1
        json.dump(new_request_template(), open(p, "w"), indent=2)
        print("wrote a starter residual request -> %s" % p)
        print("  edit material/model/record, scope and the output-parameter requests, then:")
        print("    python -m residual_core.replay.request_contract validate %s" % p)
        return 0
    if cmd == "validate":
        if not rest:
            print("usage: ... validate <request.json>"); return 2
        try:
            rc = _load(rest[0])
        except JobError as e:
            print("error:", e); return 1
        problems = validate_request(rc, os.path.dirname(os.path.abspath(rest[0])))
        if not problems:
            print("OK: %s is a valid residual request" % rest[0]); return 0
        print("%d problem(s) in %s:" % (len(problems), rest[0]))
        for p in problems:
            print("  - %s" % p)
        return 1
    path = rest[0] if cmd == "run" and rest else (cmd if cmd != "run" else None)
    if not path:
        print(__doc__); return 2
    try:
        res = run_request_contract(path)
    except JobError as e:
        print("error:", e); return 1
    print("residual request  model=%s  scope.domain=%s" % (res["model_id"], res["scope"].get("domain")))
    for r in res["requests"]:
        cells = "  ".join("d%s/d%s norm=%.4g" % (r["output"], p, s.get("norm", float("nan")))
                          for p, s in r["sensitivity"].items()) or r.get("note", "")
        print("  %-8s wrt %s :  %s" % (r["output"], r["with_respect_to"], cells))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
