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


def _derivative_label(exps, names) -> str:
    """Human-readable derivative label, e.g. (2,0)+[k,f] -> 'd2/dk2';
    (1,1) -> 'd2/dk_df'; (1,) -> 'd/dk'."""
    p = int(sum(exps))
    num = "d" if p == 1 else "d%d" % p
    den = "_".join("d%s%s" % (names[i], "" if k == 1 else str(int(k)))
                   for i, k in enumerate(exps) if k)
    return "%s/%s" % (num, den)


def direction_columns(labels_p, names) -> List[Dict[str, Any]]:
    """Column metadata for one order: index, exponents, label, recovery_factor.

    ``recovery_factor`` = prod_i (kappa_i!). It converts a raw OTI **Taylor
    coefficient** into the true partial **derivative**:

        derivative = recovery_factor * coefficient

    It is 1 for every order-1 direction and for mixed directions whose exponents
    are all 1, and > 1 whenever a parameter is repeated (e.g. d2/dk2 -> 2! = 2).
    """
    cols = []
    for i, d in enumerate(labels_p):
        exps = [int(e) for e in d["exponents"]]
        cols.append({
            "index": i,
            "exponents": exps,
            "label": _derivative_label(exps, names),
            "oti_label": d["label"],
            "recovery_factor": float(d["recovery_factor"]),
            "parameter_names": names,
        })
    return cols


def _factors(labels_p) -> np.ndarray:
    """(ncols,) array of recovery factors, aligned with the array columns."""
    return np.array([float(d["recovery_factor"]) for d in labels_p], dtype=float)


def _exponents(labels_p) -> np.ndarray:
    """(ncols, m) integer exponent matrix, aligned with the array columns."""
    return np.array([[int(e) for e in d["exponents"]] for d in labels_p], dtype=int)


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
    # Only write the real residual if it actually exists. In black-box mode the
    # provider never returns it, so writing a zero vector would fake equilibrium
    # in the private package (the same lie we refuse to print in public/).
    if result.get("R_real") is not None:
        np.savez(os.path.join(priv, "residual_real.npz"),
                 residual=np.asarray(result["R_real"], float))
    else:
        _dump_json(os.path.join(priv, "residual_real_UNAVAILABLE.json"), {
            "residual_available": False,
            "reason": result["validation"].get(
                "residual_free_norm_reason",
                "black-box did not expose real residual"),
            "note": "No residual_real.npz is written: the provider never returned "
                    "the real residual vector, so there is nothing to store. A zero "
                    "vector would falsely imply equilibrium was verified."})
    if result.get("tangent") is not None:
        np.savez(os.path.join(priv, "tangent.npz"), tangent=result["tangent"])
    # ---- per-order arrays -------------------------------------------------
    # The OTI evaluation yields raw TAYLOR COEFFICIENTS. The true partial
    # derivative is  derivative = recovery_factor * coefficient, with
    # recovery_factor = prod_i (kappa_i!). This is 1 at order 1 (so order-1
    # numbers are unchanged), but 2! = 2 for d2/dk2, 3! = 6 for d3/dk3, ...
    # We export BOTH, plus the direction map, so nothing has to be inferred.
    for p, arr in R_orders.items():
        fac = _factors(labels[p])                    # (ncols,)
        exps = _exponents(labels[p])                 # (ncols, m)
        np.savez(
            os.path.join(priv, "rhs_order%d.npz" % p),
            # legacy keys (RAW COEFFICIENTS -- kept for backward compatibility)
            residual=arr, rhs=-arr,
            # explicit coefficient arrays
            residual_coefficients=arr, rhs_coefficients=-arr,
            # recovered partial derivatives
            residual_derivatives=arr * fac, rhs_derivatives=(-arr) * fac,
            direction_exponents=exps, recovery_factors=fac)
    for p, arr in U_orders.items():
        fac = _factors(labels[p])
        exps = _exponents(labels[p])
        np.savez(
            os.path.join(priv, "solution_sensitivities_order%d.npz" % p),
            # legacy key (RAW COEFFICIENTS -- kept for backward compatibility)
            U=arr,
            U_coefficients=arr,
            U_derivatives=arr * fac,
            direction_exponents=exps, recovery_factors=fac)
    for p in sorted(labels):
        _dump_json(os.path.join(priv, "direction_map_order%d.json" % p), {
            "order": p,
            "parameter_names": names,
            "note": "derivative = recovery_factor * coefficient "
                    "(recovery_factor = prod_i kappa_i!)",
            "columns": direction_columns(labels[p], names)})
    _dump_json(os.path.join(priv, "validation_full.json"), result["validation"])

    # ---- public (safe to share) -----------------------------------------
    # PUBLIC REPORTS QUOTE RECOVERED DERIVATIVES, not raw coefficients. At order 1
    # every recovery factor is 1, so order-1 numbers are numerically unchanged; at
    # order >= 2 this is what stops repeated directions (d2/dk2 etc.) from being
    # under-reported by a factorial.
    U_deriv = {p: U_orders[p] * _factors(labels[p]) for p in U_orders}
    R_deriv = {p: R_orders[p] * _factors(labels[p]) for p in R_orders}

    ranking = _parameter_ranking(names, U_deriv)      # order-1: factors are 1
    _write_csv(os.path.join(pub, "parameter_ranking.csv"),
               ["rank", "parameter", "order1_sensitivity_norm"],
               [[r["rank"], r["parameter"], "%.6e" % r["norm"]] for r in ranking])

    norm_rows = []
    for p in sorted(U_deriv):
        for col, cmeta in enumerate(direction_columns(labels[p], names)):
            norm_rows.append([
                p, cmeta["label"], cmeta["oti_label"],
                "%g" % cmeta["recovery_factor"],
                "%.6e" % float(np.linalg.norm(U_deriv[p][:, col])),
                "%.6e" % float(np.linalg.norm(R_deriv[p][:, col]))])
    _write_csv(os.path.join(pub, "sensitivity_norms.csv"),
               ["order", "direction", "oti_direction", "recovery_factor",
                "solution_sensitivity_norm", "rhs_norm"],
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

    lines += [
        "",
        "> **Convention.** All numbers in this public report are **recovered partial"
        " derivatives**, not raw OTI Taylor coefficients:"
        " `derivative = recovery_factor * coefficient`, with"
        " `recovery_factor = prod_i (kappa_i!)`. It is 1 for every order-1 direction"
        " (so order-1 values are identical either way) and greater than 1 for"
        " repeated directions (e.g. `d2/dk2` -> 2! = 2). `sensitivity_norms.csv`"
        " lists the factor used for each direction. The private package ships both"
        " conventions (`*_coefficients` and `*_derivatives`) plus"
        " `direction_map_order<p>.json`.",
    ]

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
