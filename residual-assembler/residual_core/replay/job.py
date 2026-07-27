"""Program 2 job runner: one sensitivity_job.json drives the whole run.

Consumes the direct OTI .obj + completed contract, a converged analysis record,
requested parameters and physical outputs (with region + reduction), replays the
material, assembles R / DR_DP / K, solves K DU_DP = -DR_DP, evaluates the
requested DQ_DP, and writes a result package. The FE analysis is REPLAYED, never
re-solved.

    python -m residual_core.replay.job init [job.json]     # write a starter job
    python -m residual_core.replay.job validate <job.json> # check it, no run
    python -m residual_core.replay.job run <job.json>       # (or just <job.json>)

The job is validated up front: schema, required blocks, referenced-file existence,
parameter names against the material contract, and per-output type / component /
region / reduction. All problems are reported at once with actionable messages,
before the (expensive) object link -- never a raw KeyError or Fortran dump.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
from typing import Any, Dict, List

import numpy as np

from ..formulations import c3d8_kernel as kern
from ..interface import versions as _v
from .objlink import package_from_contract
from .record import ReplayRecord
from .engine import replay_sensitivities

SCHEMA = _v.SENSITIVITY_JOB
_UCOMP = {"U1": 1, "U2": 2, "U3": 3}
_RCOMP = {"RF1": 1, "RF2": 2, "RF3": 3}
_SCOMP = {"S11": 0, "S22": 1, "S33": 2, "S12": 3, "S13": 4, "S23": 5}

# per output type: the allowed component set and the DEFAULT reduction
_ALLOWED_COMPONENTS = {"nodal_displacement": _UCOMP, "reaction_force": _RCOMP,
                       "integration_point_stress": _SCOMP}
_DEFAULT_REDUCTION = {"nodal_displacement": "average", "reaction_force": "sum",
                      "integration_point_stress": "volume_average"}
_REDUCTIONS = {"average", "sum", "volume_average", "max_abs", "none"}
_OUTPUT_TYPES = set(_ALLOWED_COMPONENTS)


class JobError(Exception):
    pass


def _abs(base: str, rel: str) -> str:
    return rel if os.path.isabs(rel) else os.path.join(base, rel)


def _param_names(job: Dict[str, Any]) -> List[str]:
    """Accept parameters as [{"name": "E"}, ...] or the shorthand ["E", ...]."""
    out = []
    for p in (job.get("parameters") or []):
        out.append(p if isinstance(p, str) else (p.get("name") if isinstance(p, dict) else None))
    return out


# ---------------------------------------------------------------------------
# up-front validation (report style: collect everything, never raise)
# ---------------------------------------------------------------------------
def validate_job(job: Any, base_dir: str) -> List[str]:
    """Return a list of human-readable problems with a sensitivity job (empty if
    valid). Checks structure, referenced-file existence, parameter names against
    the material contract, and each output's type/component/region/reduction."""
    E: List[str] = []
    if not isinstance(job, dict):
        return ["the job file must be a JSON object"]
    if job.get("schema") != SCHEMA:
        E.append("schema must be %r, got %r" % (SCHEMA, job.get("schema")))

    mat = job.get("material")
    if not isinstance(mat, dict):
        E.append("'material' is required: an object with 'oti_umat' (.obj) and 'contract' (.json)")
        mat = {}
    else:
        for k, what in (("oti_umat", "the OTI .obj"), ("contract", "the completed contract .json")):
            v = mat.get(k)
            if not v:
                E.append("material.%s is required (%s)" % (k, what))
            elif not os.path.exists(_abs(base_dir, v)):
                E.append("material.%s not found: %s" % (k, v))

    an = job.get("analysis")
    if not isinstance(an, dict):
        E.append("'analysis' is required: an object with a 'record' (the converged run)")
    else:
        v = an.get("record")
        if not v:
            E.append("analysis.record is required (the converged analysis record)")
        elif not os.path.exists(_abs(base_dir, v)):
            E.append("analysis.record not found: %s" % v)

    pn = _param_names(job)
    if not pn:
        E.append('parameters must be a non-empty list, e.g. [{"name": "E"}, {"name": "nu"}]')
    elif any(not x for x in pn):
        E.append("every entry in 'parameters' needs a name")
    else:
        cpath = _abs(base_dir, mat["contract"]) if mat.get("contract") else None
        if cpath and os.path.exists(cpath):
            try:
                avail = [p["name"] for p in json.load(open(cpath)).get("parameters", [])]
                for name in pn:
                    if name not in avail:
                        E.append("parameter %r is not in the material contract; available: %s"
                                 % (name, ", ".join(avail)))
            except (OSError, ValueError, KeyError, TypeError):
                pass  # a malformed contract is reported by its own loader downstream

    outs = job.get("outputs")
    if not isinstance(outs, list) or not outs:
        E.append("'outputs' must be a non-empty list of requested responses")
    else:
        seen = set()
        for i, o in enumerate(outs):
            tag = "outputs[%d]" % i
            if not isinstance(o, dict):
                E.append("%s must be an object" % tag); continue
            oid = o.get("output_id")
            if not oid:
                E.append("%s.output_id is required" % tag)
            elif oid in seen:
                E.append("duplicate output_id %r" % oid)
            else:
                seen.add(oid)
            typ = o.get("type")
            if typ not in _OUTPUT_TYPES:
                E.append("%s.type %r invalid; use one of %s" % (tag, typ, sorted(_OUTPUT_TYPES)))
                continue
            comp = o.get("component")
            allowed = _ALLOWED_COMPONENTS[typ]
            if comp not in allowed:
                E.append("%s.component %r invalid for %s; use one of %s"
                         % (tag, comp, typ, sorted(allowed)))
            red = o.get("reduction")
            if red is not None and red not in _REDUCTIONS:
                E.append("%s.reduction %r invalid; use one of %s" % (tag, red, sorted(_REDUCTIONS)))
            region = o.get("region", {})
            if not isinstance(region, dict):
                E.append("%s.region must be an object" % tag)
            elif typ in ("nodal_displacement", "reaction_force"):
                if "node" not in region and "node_set" not in region:
                    E.append("%s.region needs a 'node' id or a 'node_set' name" % tag)
            else:
                if "element" not in region and "element_set" not in region:
                    E.append("%s.region needs an 'element' id or an 'element_set' name" % tag)
    return E


def new_job_template() -> Dict[str, Any]:
    """A valid, ready-to-edit resasm_sensitivity_job_v1 skeleton."""
    return {
        "schema": SCHEMA,
        "material": {"oti_umat": "material/umat_oti.obj", "contract": "material/umat_oti.json"},
        "analysis": {"model": "analysis/model.inp", "record": "analysis/record.odb"},
        "location": {"step": "LOAD", "increment": "last"},
        "parameters": [{"name": "E"}, {"name": "nu"}],
        "outputs": [
            {"output_id": "tip_displacement", "type": "nodal_displacement",
             "region": {"node_set": "LOADED_FACE"}, "component": "U1", "reduction": "average"},
            {"output_id": "support_reaction", "type": "reaction_force",
             "region": {"node_set": "FIXED_FACE"}, "component": "RF1", "reduction": "sum"},
        ],
        "result": {"output_dir": "sensitivity_out", "summary_json": True, "summary_csv": True},
    }


def _nodes_of(record, region) -> List[int]:
    if "node" in region:
        return [int(region["node"])]
    ns = record.raw.get("mesh", {}).get("node_sets", {})
    name = region.get("node_set")
    if name not in ns:
        raise JobError("node_set %r not in record mesh.node_sets %r" % (name, sorted(ns)))
    return [int(n) for n in ns[name]]


def _elems_of(record, region) -> List[int]:
    if "element" in region:
        return [int(region["element"])]
    es = record.raw.get("mesh", {}).get("element_sets", {})
    name = region.get("element_set")
    if name not in es:
        raise JobError("element_set %r not in record mesh.element_sets %r" % (name, sorted(es)))
    return [int(e) for e in es[name]]


def _load_job(job_path: str) -> Dict[str, Any]:
    try:
        with open(job_path) as fh:
            return json.load(fh)
    except FileNotFoundError:
        raise JobError("job file not found: %s" % job_path)
    except ValueError as e:
        raise JobError("job file is not valid JSON (%s): %s" % (e, job_path))


def run_job(job_path: str) -> Dict[str, Any]:
    job = _load_job(job_path)
    base = os.path.dirname(os.path.abspath(job_path))
    problems = validate_job(job, base)
    if problems:
        raise JobError("this sensitivity job has %d problem(s):\n  - %s"
                       % (len(problems), "\n  - ".join(problems)))

    def _p(rel):
        return _abs(base, rel)

    contract = json.load(open(_p(job["material"]["contract"])))
    work = tempfile.mkdtemp(prefix="resasm_job_")
    pkg_dir, manifest = package_from_contract(_p(job["material"]["oti_umat"]), contract, work)
    from .package import MaterialPackage
    package = MaterialPackage(manifest, base_dir=pkg_dir)
    record = ReplayRecord.load(_p(job["analysis"]["record"]))
    params = _param_names(job)

    # expand every requested output into element/point member requests
    point_reqs: List[Dict[str, Any]] = []
    plan: List[Dict[str, Any]] = []
    for o in job["outputs"]:
        typ, region = o["type"], o.get("region", {})
        # type-specific default reduction (fixes the "none" no-op default);
        # an explicit reduction -- including "none" -- is always honoured.
        red = o.get("reduction") or _DEFAULT_REDUCTION[typ]
        members: List[int] = []
        if typ == "nodal_displacement":
            comp = _UCOMP[o["component"]]
            for n in _nodes_of(record, region):
                idx = len(point_reqs)
                point_reqs.append({"type": "displacement", "node": n, "dof": comp})
                members.append(idx)
            plan.append({"output_id": o["output_id"], "reduction": red, "members": members, "weights": None})
        elif typ == "reaction_force":
            comp = _RCOMP[o["component"]]
            for n in _nodes_of(record, region):
                idx = len(point_reqs)
                point_reqs.append({"type": "reaction", "node": n, "dof": comp})
                members.append(idx)
            plan.append({"output_id": o["output_id"], "reduction": red, "members": members, "weights": None})
        elif typ == "integration_point_stress":
            comp = _SCOMP[o["component"]]
            weights = []
            model = record.build_model()
            for eid in _elems_of(record, region):
                el = model.elements[eid]
                Xe = model.coords_of(el.connectivity)
                for k in range(kern.ABAQUS_C3D8_GAUSS.points.shape[0]):
                    _, detJ = kern.b_matrix_reference(Xe, kern.ABAQUS_C3D8_GAUSS.points[k])
                    w = kern.ABAQUS_C3D8_GAUSS.weights[k]
                    idx = len(point_reqs)
                    point_reqs.append({"type": "stress", "element": eid, "ip": k, "component": comp})
                    members.append(idx)
                    weights.append(detJ * w)
            plan.append({"output_id": o["output_id"], "reduction": red,
                         "members": members, "weights": np.array(weights)})

    res = replay_sensitivities(record, package, params, point_reqs)

    # reduce
    outputs = []
    for pl in plan:
        vals = [res.outputs[i]["value"] for i in pl["members"]]
        grads = {p: np.array([res.outputs[i]["dq_dp"][p] for i in pl["members"]]) for p in params}
        red = pl["reduction"]
        if red in ("average", "none") and len(vals) == 1:
            value = vals[0]; dq = {p: float(grads[p][0]) for p in params}
        elif red == "average":
            value = float(np.mean(vals)); dq = {p: float(np.mean(grads[p])) for p in params}
        elif red == "sum":
            value = float(np.sum(vals)); dq = {p: float(np.sum(grads[p])) for p in params}
        elif red == "max_abs":
            j = int(np.argmax(np.abs(vals))); value = float(vals[j]); dq = {p: float(grads[p][j]) for p in params}
        elif red == "volume_average":
            w = pl["weights"]; W = float(np.sum(w))
            value = float(np.dot(w, vals) / W); dq = {p: float(np.dot(w, grads[p]) / W) for p in params}
        else:  # "none" over several members: no single reduction is defined
            value = float(vals[0]); dq = {p: float(grads[p][0]) for p in params}
        outputs.append({"output_id": pl["output_id"], "value": value, "dq_dp": dq})

    result = {
        "schema": _v.SENSITIVITY_RESULT,
        "job_id": job.get("job", {}).get("job_id", os.path.basename(job_path)),
        "model_id": package.model_id, "parameters": params,
        "du_dp_norms": {p: float(np.linalg.norm(res.du(p))) for p in params},
        "outputs": outputs,
        "equilibrium": {"free_norm": res.equilibrium_free_norm, "norm": res.equilibrium_norm},
        "stress_check": res.stress_check, "diagnostics": res.diagnostics,
    }

    rc = job.get("result", {})
    odir = _p(rc.get("output_dir", "sensitivity_out"))
    os.makedirs(odir, exist_ok=True)
    if rc.get("summary_json", True):
        json.dump(result, open(os.path.join(odir, "sensitivity_result.json"), "w"), indent=2)
    if rc.get("summary_csv"):
        with open(os.path.join(odir, "summary.csv"), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["output_id", "value"] + ["dq_d%s" % p for p in params])
            for o in outputs:
                w.writerow([o["output_id"], o["value"]] + [o["dq_dp"][p] for p in params])
    if rc.get("full_displacement_sensitivity_fields"):
        np.savez(os.path.join(odir, "du_dp_fields.npz"), **{p: res.du(p) for p in params})
    return result


def _main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__); return 0 if argv else 2
    cmd, rest = argv[0], argv[1:]

    if cmd == "init":
        out = rest[0] if rest else "sensitivity_job.json"
        if os.path.exists(out):
            print("refusing to overwrite existing %s" % out); return 1
        json.dump(new_job_template(), open(out, "w"), indent=2)
        print("wrote a starter sensitivity job -> %s" % out)
        print("  edit the material/analysis paths, parameters and outputs, then:")
        print("    python -m residual_core.replay.job validate %s" % out)
        print("    python -m residual_core.replay.job run %s" % out)
        return 0

    if cmd == "validate":
        if not rest:
            print("usage: ... validate <job.json>"); return 2
        try:
            job = _load_job(rest[0])
        except JobError as e:
            print("error:", e); return 1
        problems = validate_job(job, os.path.dirname(os.path.abspath(rest[0])))
        if not problems:
            print("OK: %s is a valid sensitivity job" % rest[0]); return 0
        print("%d problem(s) in %s:" % (len(problems), rest[0]))
        for p in problems:
            print("  - %s" % p)
        return 1

    # "run <job>", or the backward-compatible bare "<job>"
    path = rest[0] if cmd == "run" and rest else (cmd if cmd != "run" else None)
    if not path:
        print(__doc__); return 2
    try:
        res = run_job(path)
    except JobError as e:
        print("error:", e); return 1
    print("job %s  model=%s" % (res["job_id"], res["model_id"]))
    for o in res["outputs"]:
        print("  %-22s value=%+.6e  %s" % (o["output_id"], o["value"],
              {k: round(v, 6) for k, v in o["dq_dp"].items()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
