"""Requested outputs of a history replay: fields, reductions, weighted shares.

Request (``sensitivity_request.json``): the four keys of the bounded
presentation request -- ``outputs``, ``parameters``, ``domain``,
``increments`` -- with these extensions:

- ``field``: ``U``, ``RF`` (nodal), ``S``, ``SDV``, ``MISES`` (integration
  points; MISES is the von Mises stress, one component);
- a domain may name sets: ``{"nset": "TIP"}``, ``{"elset": "ROOTEL"}``, and
  restrict integration points with ``"points": [1, ..., 8]``;
- ``reduction``: ``component`` (exactly one location), ``sum``, ``mean``,
  ``volume_mean`` (integration-point volume weights det(J) w_q), ``max``,
  ``min``, ``L2``. A max/min is differentiable when the extreme value is
  attained at one location, or at several (e.g. mirror points of a symmetric
  mesh) whose derivatives agree; otherwise the output is refused;
- ``parameters`` may be ``"ALL"`` (the provider's order);
- optional ``weighted_shares``: ``{"field": "MISES", "domain": {...}}`` adds
  the per-increment weighted sensitivity shares (see :func:`weighted_shares`);
- optional ``full_field``: ``true`` writes every field and its parameter
  derivatives at every increment to ``fields.npz``.

Weighted sensitivity (slide 29/33: "a weighted sensitivity, so parameters
with different units can be compared"), recovered from
``results/cp_residual_sensitivities.py`` and ``results/fcc_crystal_results.py``
of the ``cross-platform-hardening`` branch: for a scalar q and parameter p_j,
``W_j = |p_j dq/dp_j|`` and the share is ``100 W_j / sum_k W_k`` at every
increment. For a field the same weight is aggregated over the domain with
the integration-point volumes:

    W_j(n) = sum_i V_i |p_j d q_i(n)/d p_j| / sum_i V_i,   share_j = W_j / sum_k W_k.

The scalar form applied to the volume-averaged q is reported next to it.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List

import numpy as np

from .history import HistoryEngine, HistoryResult, von_mises, von_mises_gradient

NODAL = {"U": 3, "RF": 3}
POINT = {"S": 6, "SDV": None, "MISES": 1}
REDUCTIONS = ("component", "sum", "mean", "volume_mean", "max", "min", "L2")
PUBLIC_FILES = ("sensitivity_results.json", "sensitivity_tables.csv", "run_report.txt")


class RequestError(ValueError):
    pass


def load_request(path):
    try:
        request = json.loads(Path(path).read_text())
    except json.JSONDecodeError as error:
        raise RequestError("sensitivity_request.json is not valid JSON: %s" % error) from error
    return request


def validate_request(request, parameters, nstatev):
    if not isinstance(request, dict):
        raise RequestError("the request must be a JSON object")
    required = {"outputs", "parameters", "domain", "increments"}
    missing = required - set(request)
    if missing:
        raise RequestError("request lacks %s" % sorted(missing))
    extra = set(request) - required - {"weighted_shares", "full_field"}
    if extra:
        raise RequestError("unsupported request keys %s" % sorted(extra))
    selected = parameters if request["parameters"] == "ALL" else request["parameters"]
    if (not isinstance(selected, list) or not selected or len(set(selected)) != len(selected)
            or any(name not in parameters for name in selected)):
        raise RequestError("parameters must be ALL or a unique nonempty list drawn from %s" % parameters)
    _check_domain(request["domain"], "request domain")
    outputs = request["outputs"]
    if not isinstance(outputs, list):
        raise RequestError("outputs must be a list")
    names = []
    for output in outputs:
        if not isinstance(output, dict) or not {"name", "field", "component", "reduction"} <= set(output):
            raise RequestError("each output needs name, field, component, reduction")
        if set(output) - {"name", "field", "component", "reduction", "domain"}:
            raise RequestError("unsupported output keys in %s" % output.get("name"))
        names.append(output["name"])
        field = output["field"]
        width = NODAL.get(field) or ({"S": 6, "SDV": nstatev, "MISES": 1}.get(field))
        if width is None:
            raise RequestError("unsupported field %r; use U, RF, S, SDV, MISES" % field)
        component = output["component"]
        if component != "ALL" and (type(component) is not int or not 1 <= component <= width):
            raise RequestError("%s: component must be ALL or 1..%d" % (output["name"], width))
        if output["reduction"] not in REDUCTIONS:
            raise RequestError("%s: reduction must be one of %s" % (output["name"], REDUCTIONS))
        if output["reduction"] == "volume_mean" and field in NODAL:
            raise RequestError("%s: volume_mean applies to integration-point fields" % output["name"])
        if "domain" in output:
            _check_domain(output["domain"], output["name"])
    if len(set(names)) != len(names) or not all(isinstance(n, str) and n for n in names):
        raise RequestError("output names must be unique nonempty strings")
    increments = request["increments"]
    if increments not in ("ALL", "LAST") and (
            not isinstance(increments, list) or not increments or len(set(increments)) != len(increments)
            or any(type(n) is not int or n < 1 for n in increments)):
        raise RequestError("increments must be ALL, LAST or a unique list of one-based numbers")
    shares = request.get("weighted_shares")
    if shares is not None:
        if not isinstance(shares, dict) or set(shares) - {"field", "domain"} or shares.get("field", "MISES") not in POINT:
            raise RequestError("weighted_shares takes field (S/SDV/MISES, default MISES) and domain")
        if "domain" in shares:
            _check_domain(shares["domain"], "weighted_shares")
    if "full_field" in request and not isinstance(request["full_field"], bool):
        raise RequestError("full_field must be true or false")
    return selected


def _check_domain(domain, where):
    allowed = {"nodes", "elements", "nset", "elset", "points"}
    if not isinstance(domain, dict) or set(domain) - allowed:
        raise RequestError("%s: domain keys must be among %s" % (where, sorted(allowed)))
    for key in ("nodes", "elements"):
        if key in domain and domain[key] != "ALL" and (
                not isinstance(domain[key], list) or not domain[key]
                or any(type(v) is not int for v in domain[key])):
            raise RequestError("%s: %s must be ALL or a list of ids" % (where, key))
    if "points" in domain and (not isinstance(domain["points"], list) or not domain["points"]
                               or any(type(p) is not int or not 1 <= p <= 8 for p in domain["points"])):
        raise RequestError("%s: points must be a list of integration points 1..8" % where)


class Fields:
    """Stacked history fields (increment axis first, increment 0 = virgin)."""

    def __init__(self, engine: HistoryEngine, result: HistoryResult):
        self.engine, self.result = engine, result
        model = engine.model
        nn, ne = len(model.node_ids), engine.ne
        npar = len(result.parameters)
        count = len(result.increments) + 1
        stack = lambda name, shape: np.concatenate(
            [np.zeros((1,) + shape)] + [getattr(r, name)[None] for r in result.increments])
        self.times = result.times
        self.U = stack("u", (engine.ndof,)).reshape(count, nn, 3)
        self.RF = np.zeros((count, engine.ndof))
        for n, record in enumerate(result.increments, 1):
            self.RF[n, engine.constrained] = record.reaction[engine.constrained]
        self.RF = self.RF.reshape(count, nn, 3)
        self.S = stack("stress", (ne, 8, 6))
        self.SDV = stack("state", (ne, 8, engine.material.nstatev))
        self.MISES = von_mises(self.S)[..., None]
        if result.sensitivities:
            self.dU = stack("du", (engine.ndof, npar)).reshape(count, nn, 3, npar)
            self.dRF = np.zeros((count, engine.ndof, npar))
            for n, record in enumerate(result.increments, 1):
                self.dRF[n, engine.constrained] = record.dreaction[engine.constrained]
            self.dRF = self.dRF.reshape(count, nn, 3, npar)
            self.dS = stack("dstress", (ne, 8, 6, npar))
            self.dSDV = stack("dstate", (ne, 8, engine.material.nstatev, npar))
            self.dMISES = np.einsum("neqa,neqam->neqm", von_mises_gradient(self.S), self.dS)[..., None, :]

    def value(self, field):
        return getattr(self, field)

    def derivative(self, field):
        return getattr(self, "d" + field)


def _nodes_in(domain, model):
    if "nset" in domain:
        ids = model.node_set(domain["nset"])
    elif domain.get("nodes", "ALL") == "ALL":
        ids = list(model.node_ids)
    else:
        ids = domain["nodes"]
    position = {int(node): i for i, node in enumerate(model.node_ids)}
    missing = [n for n in ids if int(n) not in position]
    if missing:
        raise RequestError("unknown node ids %s" % missing[:5])
    return [int(n) for n in ids], np.array([position[int(n)] for n in ids])


def _elements_in(domain, model):
    if "elset" in domain:
        ids = model.element_set(domain["elset"])
    elif domain.get("elements", "ALL") == "ALL":
        ids = list(model.element_ids)
    else:
        ids = domain["elements"]
    position = {int(e): i for i, e in enumerate(model.element_ids)}
    missing = [e for e in ids if int(e) not in position]
    if missing:
        raise RequestError("unknown element ids %s" % missing[:5])
    points = np.array(domain.get("points", list(range(1, 9)))) - 1
    return [int(e) for e in ids], np.array([position[int(e)] for e in ids]), points


def select(fields: Fields, field, component, domain, increment):
    """Values (m,), derivatives (m, npar), volume weights (m,) or None, labels."""
    model = fields.engine.model
    if field in NODAL:
        ids, positions = _nodes_in(domain, model)
        value = fields.value(field)[increment][positions]                        # (k, 3)
        derivative = fields.derivative(field)[increment][positions]             # (k, 3, npar)
        labels = [(n, c + 1) for n in ids for c in range(3)]
        weights = None
    else:
        ids, positions, points = _elements_in(domain, model)
        value = fields.value(field)[increment][positions][:, points]             # (k, q, c)
        derivative = fields.derivative(field)[increment][positions][:, points]  # (k, q, c, npar)
        weights = np.repeat(fields.engine.w[positions][:, points][..., None], value.shape[-1], axis=-1)
        labels = [(e, int(q) + 1, c + 1) for e in ids for q in points for c in range(value.shape[-1])]
    width = value.shape[-1]
    if component != "ALL":
        value = value[..., component - 1]
        derivative = derivative[..., component - 1, :]
        if weights is not None:
            weights = weights[..., component - 1]
        labels = [label for label in labels if label[-1] == component]
    npar = derivative.shape[-1]
    return (value.reshape(-1), derivative.reshape(-1, npar),
            None if weights is None else weights.reshape(-1), labels, width)


def reduce(values, derivatives, weights, reduction):
    """Scalar and its parameter gradient; refuses non-differentiable extremes."""
    if not values.size:
        raise RequestError("empty output domain")
    if reduction == "component":
        if values.size != 1:
            raise RequestError("component reduction needs exactly one location (got %d)" % values.size)
        return float(values[0]), derivatives[0], {}
    if reduction == "sum":
        return float(values.sum()), derivatives.sum(axis=0), {}
    if reduction == "mean":
        return float(values.mean()), derivatives.mean(axis=0), {}
    if reduction == "volume_mean":
        total = weights.sum()
        return float(weights @ values / total), weights @ derivatives / total, {}
    if reduction == "L2":
        norm = float(np.linalg.norm(values))
        if norm == 0:
            raise RequestError("L2 of a zero field is not differentiable")
        return norm, values @ derivatives / norm, {}
    if reduction in ("max", "min"):
        extreme = values.max() if reduction == "max" else values.min()
        scale = max(abs(extreme), np.abs(values).max(), 1e-300)
        tied = np.flatnonzero(np.abs(values - extreme) <= 1e-9 * scale)
        gradient = derivatives[tied[0]]
        spread = np.abs(derivatives[tied] - gradient).max() if len(tied) > 1 else 0.0
        if spread > 1e-6 * max(np.abs(gradient).max(), 1e-300):
            raise RequestError("%s is attained at %d locations whose derivatives differ; the %s is not "
                               "differentiable there" % (reduction, len(tied), reduction))
        return float(extreme), gradient, {"arg_locations": int(len(tied)), "arg_index": int(tied[0])}
    raise RequestError("unsupported reduction %r" % reduction)


def weighted_shares(fields: Fields, parameter_values, field="MISES", domain=None):
    """Per-increment weighted shares (percent), field-aggregated and scalar forms."""
    domain = domain or {"elements": "ALL"}
    rows = []
    p = np.abs(np.asarray(parameter_values, dtype=float))
    for n in range(1, len(fields.times)):
        values, derivatives, weights, _, width = select(fields, field, 1 if field == "MISES" else "ALL",
                                                        domain, n)
        weighted = np.abs(derivatives) * p                                 # (m, npar)
        aggregate = weights @ weighted / weights.sum()
        scalar = np.abs(weights @ derivatives / weights.sum()) * p
        rows.append({"increment": n, "time": float(fields.times[n]),
                     "field_weight": aggregate, "field_share": 100 * aggregate / max(aggregate.sum(), 1e-300),
                     "scalar_weight": scalar, "scalar_share": 100 * scalar / max(scalar.sum(), 1e-300),
                     "volume_mean": float(weights @ values / weights.sum())})
    return rows


def increments_of(request, count):
    selected = request["increments"]
    numbers = list(range(1, count + 1)) if selected == "ALL" else [count] if selected == "LAST" else selected
    if any(n > count for n in numbers):
        raise RequestError("requested increment exceeds the %d replayed increments" % count)
    return numbers


def scalar_rows(fields: Fields, request, parameters) -> List[Dict]:
    result = fields.result
    columns = [result.parameters.index(name) for name in parameters]
    values_p = result.parameter_values[columns]
    rows = []
    for number in increments_of(request, len(result.increments)):
        for output in request["outputs"]:
            domain = output.get("domain", request["domain"])
            values, derivatives, weights, _, _ = select(fields, output["field"], output["component"],
                                                         domain, number)
            value, gradient, info = reduce(values, derivatives[:, columns], weights, output["reduction"])
            rows.append({"output": output["name"], "field": output["field"],
                         "component": output["component"], "reduction": output["reduction"],
                         "domain": domain, "increment": number, "time": float(fields.times[number]),
                         "value": value, "derivatives": dict(zip(parameters, gradient.tolist())),
                         "weighted": dict(zip(parameters, (gradient * values_p).tolist())), **info})
    return rows


def write_outputs(out, *, fields: Fields, request, parameters, report: Dict, private: Dict):
    """Public: the three slide-12 files (+ shares/full field when requested). Private: the rest."""
    out = Path(out)
    public_rows = scalar_rows(fields, request, parameters)
    result = fields.result
    summary = {
        "schema": "resasm_history_results_v1", "status": "completed", "request": request,
        "scope": report["scope"], "metadata": report["metadata"], "results": public_rows,
    }
    shares = None
    if request.get("weighted_shares") is not None:
        spec = request["weighted_shares"]
        shares = weighted_shares(fields, result.parameter_values, spec.get("field", "MISES"),
                                 spec.get("domain"))
        summary["weighted_shares"] = {
            "definition": ("field_share_j(n) = W_j/sum_k W_k, W_j = sum_i V_i |p_j dq_i/dp_j| / sum_i V_i "
                           "(q = %s over the domain, V_i integration-point volumes); scalar_share uses "
                           "|p_j d(mean_V q)/dp_j|" % spec.get("field", "MISES")),
            "parameters": result.parameters,
            "rows": [{"increment": r["increment"], "time": r["time"], "volume_mean": r["volume_mean"],
                      "field_share_percent": dict(zip(result.parameters, r["field_share"].tolist())),
                      "scalar_share_percent": dict(zip(result.parameters, r["scalar_share"].tolist()))}
                     for r in shares]}
        with (out / "sensitivity_shares.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["increment", "time", "parameter", "field_share_percent", "scalar_share_percent",
                             "field_weight", "volume_mean"])
            for r in shares:
                for j, name in enumerate(result.parameters):
                    writer.writerow([r["increment"], r["time"], name, r["field_share"][j],
                                     r["scalar_share"][j], r["field_weight"][j], r["volume_mean"]])
    (out / PUBLIC_FILES[0]).write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    with (out / PUBLIC_FILES[1]).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["output", "field", "component", "reduction", "domain", "increment", "time",
                         "value", "parameter", "derivative", "weighted_derivative"])
        for row in public_rows:
            for name in parameters:
                writer.writerow([row["output"], row["field"], row["component"], row["reduction"],
                                 json.dumps(row["domain"], sort_keys=True), row["increment"], row["time"],
                                 row["value"], name, row["derivatives"][name], row["weighted"][name]])
    if request.get("full_field"):
        model = fields.engine.model
        np.savez_compressed(out / "fields.npz", parameters=np.array(result.parameters),
                            parameter_values=result.parameter_values, time=fields.times,
                            node_labels=model.node_ids, elem_labels=model.element_ids,
                            U=fields.U, dU=fields.dU, RF=fields.RF, dRF=fields.dRF, S=fields.S,
                            dS=fields.dS, SDV=fields.SDV, dSDV=fields.dSDV, MISES=fields.MISES[..., 0],
                            dMISES=fields.dMISES[..., 0, :], ip_volume=fields.engine.w)
    (out / PUBLIC_FILES[2]).write_text(render_report(report))
    private_dir = out / "private"
    private_dir.mkdir(exist_ok=True)
    (private_dir / "run_details.json").write_text(json.dumps(private, indent=1, allow_nan=False) + "\n")
    return summary, shares


def render_report(report):
    k = report["report_fields"]
    lines = ["Status: %s" % report["status"]]
    for key in ("command executed", "residual assembled", "equilibrium checked", "equilibrium passed",
                "tangent available", "tangent verified", "derivative calculated", "derivative verified",
                "reference resolved", "Abaqus comparison available", "unsupported feature detected",
                "public and private outputs separated"):
        lines.append("%s: %s" % (key[0].upper() + key[1:], k[key]))
    lines.append("")
    for key, value in report["details"].items():
        lines.append("%s: %s" % (key, value))
    return "\n".join(lines) + "\n"
