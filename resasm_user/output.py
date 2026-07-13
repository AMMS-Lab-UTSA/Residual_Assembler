"""Write the sensitivity results as a private/ (full arrays) and public/ (safe to
share) package.

Public NEVER contains: source, mesh, the full residual vector, the full tangent
matrix, state variables, or private parameter values (only names + norms).
"""

from __future__ import annotations

import csv
import json
import os
import time as _time
from typing import Any, Dict, List

import numpy as np


def write_outputs(out_dir: str, cfg, result: Dict[str, Any]) -> Dict[str, str]:
    priv = os.path.join(out_dir, "private")
    pub = os.path.join(out_dir, "public")
    os.makedirs(priv, exist_ok=True)
    os.makedirs(pub, exist_ok=True)

    names: List[str] = result["parameter_names"]
    param_map = {n: i + 1 for i, n in enumerate(names)}
    R_orders = result["R_orders"]
    U_orders = result["U_orders"]
    labels = result["direction_labels"]        # {p: [ {label, exponents, recovery_factor} ]}
    ndof = int(result["ndof"])
    free = np.asarray(result["free_mask"], bool)

    # ---- private ---------------------------------------------------------
    _dump_json(os.path.join(priv, "metadata.json"), {
        "name": cfg.name, "unknowns": ndof, "order": cfg.order,
        "backend": cfg.backend, "residual_type": cfg.residual_type,
        "tangent_source": result["tangent_source"],
        "parameters": names, "created": _now(),
        "engine": result.get("diagnostics", {}).get("engine")})
    _dump_json(os.path.join(priv, "parameter_map.json"), param_map)
    _dump_json(os.path.join(priv, "dof_map.json"),
               {"ndof": ndof, "free_dofs": [int(i) for i in np.where(free)[0]],
                "prescribed_dofs": [int(i) for i in np.where(~free)[0]]})
    np.savez(os.path.join(priv, "residual_real.npz"), residual=result["R_real"])
    if result.get("tangent") is not None:
        np.savez(os.path.join(priv, "tangent.npz"), tangent=result["tangent"])
    for p, arr in R_orders.items():
        np.savez(os.path.join(priv, "rhs_order%d.npz" % p), rhs=-arr, residual=arr)
    for p, arr in U_orders.items():
        np.savez(os.path.join(priv, "solution_sensitivities_order%d.npz" % p), U=arr)
    _dump_json(os.path.join(priv, "validation_full.json"), result["validation"])

    # ---- public (safe to share) -----------------------------------------
    ranking = _parameter_ranking(names, U_orders)
    _write_csv(os.path.join(pub, "parameter_ranking.csv"),
               ["rank", "parameter", "order1_sensitivity_norm"],
               [[r["rank"], r["parameter"], "%.6e" % r["norm"]] for r in ranking])

    norm_rows = []
    for p in sorted(U_orders):
        for col, lab in enumerate(labels[p]):
            norm_rows.append([p, lab["label"],
                              "%.6e" % float(np.linalg.norm(U_orders[p][:, col])),
                              "%.6e" % float(np.linalg.norm(R_orders[p][:, col]))])
    _write_csv(os.path.join(pub, "sensitivity_norms.csv"),
               ["order", "direction", "solution_sensitivity_norm", "rhs_norm"],
               norm_rows)

    _dump_json(os.path.join(pub, "timing.json"), result.get("timing", {}))
    v = result["validation"]
    summary = {
        # null (not 0.0) when the residual was never available -- e.g. black-box.
        "residual_free_norm": v.get("residual_free_norm"),
        "residual_free_norm_available": v.get("residual_free_norm") is not None,
        "tangent_source": result["tangent_source"],
        "orders_solved": sorted(U_orders.keys()),
        "parameters": names,
        "status": v.get("status", "ok"),
        # accurately named: this is a residual-derivative (RHS) FD check, not a
        # re-solve of the nonlinear problem.
        "rhs_finite_difference_check": v.get("rhs_finite_difference_check"),
        "solution_finite_difference_check": None,   # not implemented
    }
    if v.get("residual_free_norm_reason"):
        summary["residual_free_norm_reason"] = v["residual_free_norm_reason"]
    if v.get("notes"):
        summary["notes"] = v["notes"]
    _dump_json(os.path.join(pub, "validation_summary.json"), summary)
    _write_summary_md(os.path.join(pub, "summary.md"), cfg, result, ranking)

    return {"private": priv, "public": pub}


def _parameter_ranking(names, U_orders):
    rows = []
    U1 = U_orders.get(1)
    for i, n in enumerate(names):
        norm = float(np.linalg.norm(U1[:, i])) if U1 is not None and i < U1.shape[1] else 0.0
        rows.append({"parameter": n, "norm": norm})
    rows.sort(key=lambda r: r["norm"], reverse=True)
    for k, r in enumerate(rows, start=1):
        r["rank"] = k
    return rows


def _write_summary_md(path, cfg, result, ranking):
    v = result["validation"]

    # Truthfulness: when the residual was never exposed (black-box), print an
    # explicit n/a. Printing 0.000e+00 would imply we verified equilibrium.
    rfn = v.get("residual_free_norm")
    if rfn is None:
        reason = v.get("residual_free_norm_reason",
                       "black-box did not expose real residual")
        rfn_cell = "n/a — %s" % reason
    else:
        rfn_cell = "%.3e" % float(rfn)

    lines = [
        "# Sensitivity summary — %s" % cfg.name,
        "",
        "Public report (safe to share). No source, mesh, full residual/tangent,",
        "or state variables are included.",
        "",
        "## Run",
        "",
        "| item | value |",
        "|---|---|",
        "| unknowns (DOFs) | %d |" % result["ndof"],
        "| parameters | %s |" % ", ".join(result["parameter_names"]),
        "| derivative order | %d |" % cfg.order,
        "| backend | %s |" % cfg.backend,
        "| tangent source | %s |" % result["tangent_source"],
        "| residual free-DOF norm | %s |" % rfn_cell,
        "| orders solved | %s |" % ", ".join(str(p) for p in sorted(result["U_orders"])),
        "",
    ]
    if rfn is None:
        lines += [
            "> **Equilibrium was NOT verified by this run.** The residual vector was "
            "never available (%s), so no residual norm could be computed. Nothing "
            "here should be read as confirming the supplied solution is converged."
            % v.get("residual_free_norm_reason",
                    "black-box did not expose real residual"),
            "",
        ]
    lines += [
        "## Parameter ranking (order-1 solution sensitivity magnitude)",
        "",
        "| rank | parameter | ||du/da|| |",
        "|---|---|---|",
    ]
    for r in ranking:
        lines.append("| %d | %s | %.6e |" % (r["rank"], r["parameter"], r["norm"]))

    fd = v.get("rhs_finite_difference_check")
    if fd:
        lines += ["", "## RHS finite-difference cross-check", "",
                  "Checks **d(residual)/d(parameter) at the fixed solution u**, solved "
                  "with the *same* tangent used by the hypercomplex path. It validates "
                  "the generated RHS and that solve. It does **not** re-solve the "
                  "nonlinear problem at perturbed parameters, so it is *not* a "
                  "solution-level finite-difference validation and cannot detect an "
                  "error in the tangent itself.",
                  "",
                  "| quantity | value |", "|---|---|",
                  "| max relative error | %.3e |" % fd.get("max_rel_error", float("nan")),
                  "| status | %s |" % fd.get("status", "n/a")]
    for n in v.get("notes", []) or []:
        lines += ["", "> %s" % n]

    lines += ["", "## Files", "",
              "- `sensitivity_norms.csv` — per-direction solution/RHS norms",
              "- `parameter_ranking.csv` — parameters ranked by influence",
              "- `validation_summary.json` — machine-readable status",
              "- private arrays are under `../private/` (kept local).", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def _dump_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(_jsonable(obj), fh, indent=2)


def _write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def _now():
    return _time.strftime("%Y-%m-%dT%H:%M:%S")
