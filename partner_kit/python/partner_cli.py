"""partner_cli — run the residual-sensitivity kit locally.

    python -m partner_kit.python.partner_cli --config config.json

Loads the partner's residual provider (a local Python class or a black-box
command), obtains a converged solution, generates R^(1) = dR/dp via the
hypercomplex seam (or the black-box executable), acquires the tangent T if
available, solves ``T U^(1) = -R^(1)``, validates locally, and writes:

  <out>/sensitivity_package/   PRIVATE — full arrays, may stay on the partner box
  <out>/public_report/         SHAREABLE — norms, rankings, validation summary only

Nothing proprietary is required to leave the partner machine. See
``docs/privacy_model.md`` and ``docs/output_contract.md``.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import importlib.util
import json
import os
import sys
import time as _time
from typing import Any, Dict, Optional

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_KIT_PARENT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
if _KIT_PARENT not in sys.path:
    sys.path.insert(0, _KIT_PARENT)

from partner_kit.python import sensitivity as _sens
from partner_kit.python import validators as _val
from partner_kit.python.blackbox_runner import BlackBoxRunner
from partner_kit.python.residual_provider import BlackBoxResidualProvider


# --------------------------------------------------------------------------- #
def _load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    if path.lower().endswith((".yml", ".yaml")):
        try:
            import yaml
            return yaml.safe_load(text) or {}
        except Exception:
            pass
    return json.loads(text)


def _load_provider(cfg: Dict[str, Any], base_dir: str):
    if "blackbox" in cfg:
        bb = cfg["blackbox"]
        command = bb["command"]
        # resolve any script path in the command relative to the config dir so
        # the executable is found from the runner's temp working directory
        if isinstance(command, list):
            command = [(os.path.abspath(os.path.join(base_dir, tok))
                        if (isinstance(tok, str) and os.path.exists(os.path.join(base_dir, tok)))
                        else tok) for tok in command]
        runner = BlackBoxRunner(command, workdir=bb.get("workdir"),
                                timeout=bb.get("timeout"))
        free = np.asarray(bb["free_mask"], bool) if bb.get("free_mask") else None
        return BlackBoxResidualProvider(runner, ndof=int(bb["ndof"]),
                                        parameters=list(bb["parameters"]),
                                        free_mask=free, name=bb.get("name", "blackbox"),
                                        parameter_values=bb.get("parameter_values"))
    args = cfg.get("provider_args", {}) or {}
    if cfg.get("provider_file"):
        pf = cfg["provider_file"]
        pf = pf if os.path.isabs(pf) else os.path.join(base_dir, pf)
        spec = importlib.util.spec_from_file_location("_partner_provider", pf)
        mod = importlib.util.module_from_spec(spec)
        if base_dir not in sys.path:
            sys.path.insert(0, base_dir)
        spec.loader.exec_module(mod)
        cls = getattr(mod, cfg["provider_class"])
        return cls(**args)
    dotted = cfg["provider"]
    mod_name, cls_name = dotted.split(":")
    cls = getattr(importlib.import_module(mod_name), cls_name)
    return cls(**args)


# --------------------------------------------------------------------------- #
def run_kit(cfg: Dict[str, Any], base_dir: str = ".") -> Dict[str, Any]:
    t0 = _time.perf_counter()
    provider = _load_provider(cfg, base_dir)
    params = list(cfg.get("parameters") or provider.parameters)
    provider.parameters = tuple(params)
    order = int(cfg.get("order", 1))
    time_t = tuple(cfg.get("time", (0.0, 0.0)))
    dtime = float(cfg.get("dtime", 0.0))
    algebra = cfg.get("algebra", "dual1")
    out_dir = cfg.get("output_dir", "out")
    out_dir = out_dir if os.path.isabs(out_dir) else os.path.join(base_dir, out_dir)
    state = None
    timings: Dict[str, float] = {}
    errors: list = []

    # ---- solution (supplied or solved locally) --------------------------
    sol = cfg.get("solution")
    if isinstance(sol, str):
        U = np.load(sol if os.path.isabs(sol) else os.path.join(base_dir, sol))
    elif sol is not None:
        U = np.asarray(sol, float)
    elif cfg.get("solve", False) and not isinstance(provider, BlackBoxResidualProvider):
        ts = _time.perf_counter()
        U = _sens.newton_solve(provider, u0=cfg.get("u0"))
        timings["newton_solve_s"] = _time.perf_counter() - ts
    else:
        U = np.zeros(provider.ndof)
        errors.append("no solution supplied and solve=false; using zeros (R may be nonzero)")

    # ---- R^(1) ----------------------------------------------------------
    ts = _time.perf_counter()
    R1, rhs_diag, raw = _sens.generate_rhs(provider, U, order=order, state=state,
                                           time=time_t, dtime=dtime, algebra=algebra)
    timings["generate_rhs_s"] = _time.perf_counter() - ts

    # ---- real residual + tangent ----------------------------------------
    R0 = None
    if not isinstance(provider, BlackBoxResidualProvider):
        R0, _ = _sens.assemble_real(provider, U, state, provider.parameter_values(),
                                    time_t, dtime)
    elif raw is not None and raw.get("residual_real") is not None:
        R0 = np.asarray(raw["residual_real"], float)

    ext_T = None
    if cfg.get("tangent_file"):
        tf = cfg["tangent_file"]
        tf = tf if os.path.isabs(tf) else os.path.join(base_dir, tf)
        data = np.load(tf)
        ext_T = data[data.files[0]]
    ts = _time.perf_counter()
    if isinstance(provider, BlackBoxResidualProvider) and raw is not None and raw.get("tangent") is not None:
        T, T_source = np.asarray(raw["tangent"], float), "blackbox-response"
    else:
        try:
            T, T_source = _sens.acquire_tangent(
                provider, U, state, provider.parameter_values(), time_t, dtime,
                external=ext_T)
        except Exception as exc:
            T, T_source = None, "unavailable"
            errors.append("tangent acquisition failed: %s" % exc)
    timings["acquire_tangent_s"] = _time.perf_counter() - ts

    free = provider.free_mask()

    # ---- solve T U^(1) = -R^(1) -----------------------------------------
    U1 = None
    if T is not None:
        ts = _time.perf_counter()
        try:
            U1 = _sens.solve_sensitivity(T, R1, free)
            timings["solve_s"] = _time.perf_counter() - ts
        except Exception as exc:
            errors.append("sensitivity solve failed: %s" % exc)

    # ---- validation -----------------------------------------------------
    checks = _val.run_all(provider, U, R1, T, U1, state=state, time=time_t,
                          dtime=dtime, free_mask=free)
    timings["total_s"] = _time.perf_counter() - t0

    ctx = dict(provider=provider, params=params, order=order, algebra=algebra,
               U=U, R0=R0, R1=R1, T=T, T_source=T_source, U1=U1, free=free,
               checks=checks, rhs_diag=rhs_diag, timings=timings, errors=errors)
    _write_private(out_dir, ctx)
    _write_public(out_dir, ctx, share_arrays=cfg.get("share_full_arrays", False))
    return ctx


# --------------------------------------------------------------------------- #
def _write_private(out_dir, ctx):
    d = os.path.join(out_dir, "sensitivity_package")
    os.makedirs(d, exist_ok=True)
    provider = ctx["provider"]; order = ctx["order"]
    pmap = {p: i + 1 for i, p in enumerate(ctx["params"])}

    metadata = {
        "schema": "resasm-partner-sensitivity/1",
        "target_system": "T U^(p) = -R^(p)",
        "provider": provider.name, "provider_kind": provider.kind(),
        "ndof": int(provider.ndof), "order": order, "algebra": ctx["algebra"],
        "parameters": list(ctx["params"]),
        "tangent_source": ctx["T_source"],
        "tangent_available": ctx["T"] is not None,
        "system_ready": ctx["U1"] is not None,
        "hypercomplex_ready": ctx["rhs_diag"].get("hypercomplex_ready", None),
        "rhs_diagnostics": ctx["rhs_diag"],
        "timings_s": ctx["timings"],
    }
    _json(os.path.join(d, "metadata.json"), metadata)
    free = ctx["free"]
    _json(os.path.join(d, "dof_map.json"),
          {"ndof": int(provider.ndof),
           "free_mask": (None if free is None else [bool(x) for x in np.asarray(free)])})
    _json(os.path.join(d, "parameter_map.json"), pmap)
    if ctx["R0"] is not None:
        np.savez(os.path.join(d, "residual_real.npz"), R=np.asarray(ctx["R0"], float))
    if ctx["T"] is not None:
        np.savez(os.path.join(d, "tangent.npz"), T=np.asarray(ctx["T"], float))
    else:
        _json(os.path.join(d, "tangent_operator_info.json"),
              {"available": False, "source": ctx["T_source"],
               "note": "no dense tangent; supply get_tangent / apply_tangent / "
                       "an external tangent.npz, or export the RHS for us to solve"})
    np.savez(os.path.join(d, "rhs_order_%d.npz" % order),
             R=np.asarray(ctx["R1"], float), rhs=-np.asarray(ctx["R1"], float))
    if ctx["U1"] is not None:
        np.savez(os.path.join(d, "solution_sensitivities_order_%d.npz" % order),
                 U=np.asarray(ctx["U1"], float))
    with open(os.path.join(d, "validation_report.md"), "w", encoding="utf-8") as fh:
        fh.write(_validation_md(ctx["checks"], private=True))
    _json(os.path.join(d, "diagnostics.json"),
          {"rhs_diagnostics": ctx["rhs_diag"], "errors": ctx["errors"],
           "timings_s": ctx["timings"]})


def _write_public(out_dir, ctx, share_arrays=False):
    d = os.path.join(out_dir, "public_report")
    os.makedirs(d, exist_ok=True)
    params = ctx["params"]; U1 = ctx["U1"]; R1 = ctx["R1"]

    with open(os.path.join(d, "validation_summary.md"), "w", encoding="utf-8") as fh:
        fh.write(_validation_md(ctx["checks"], private=False))
    _json(os.path.join(d, "timing_summary.json"), ctx["timings"])

    # sensitivity norms (per parameter): prefer ||dU/dp||; fall back to ||R^(1)_col||
    norms = []
    for i, p in enumerate(params):
        if U1 is not None:
            val = float(np.linalg.norm(U1[:, i]))
            metric = "||dU/dp||"
        else:
            val = float(np.linalg.norm(R1[:, i]))
            metric = "||R^(1)_col|| (RHS proxy; no tangent)"
        norms.append((p, val, metric))
    with open(os.path.join(d, "sensitivity_norms.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["parameter", "norm", "metric"])
        for p, v, mtr in norms:
            w.writerow([p, "%.8e" % v, mtr])

    ranked = sorted(norms, key=lambda x: x[1], reverse=True)
    with open(os.path.join(d, "parameter_ranking.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["rank", "parameter", "norm", "metric"])
        for r, (p, v, mtr) in enumerate(ranked, 1):
            w.writerow([r, p, "%.8e" % v, mtr])

    _json(os.path.join(d, "errors.json"), {"errors": ctx["errors"]})

    if share_arrays:            # opt-in only
        np.savez(os.path.join(d, "shared_sensitivities.npz"),
                 U=None if U1 is None else np.asarray(U1, float))


# ---- small helpers --------------------------------------------------------- #
def _json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


def _validation_md(checks, private):
    L = ["# Validation " + ("Report" if private else "Summary"), "",
         "overall: %s" % _val.overall(checks), "",
         "| check | value | tol | result | detail |",
         "|---|---|---|---|---|"]
    for c in checks:
        val = "%.3e" % c["value"] if isinstance(c["value"], (int, float)) else "—"
        tol = "%g" % c["tolerance"] if isinstance(c["tolerance"], (int, float)) else "—"
        res = {True: "PASS", False: "FAIL", None: "pending"}[c["passed"]]
        L.append("| %s | %s | %s | %s | %s |" % (c["name"], val, tol, res, c["detail"]))
    if not private:
        L += ["", "_This summary contains no mesh, source, or full residual/tangent._"]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="partner-sensitivity",
                                description="Local residual-sensitivity kit (partner-side)")
    p.add_argument("--config", required=True)
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    base_dir = os.path.dirname(os.path.abspath(args.config))
    cfg = _load_config(args.config)
    ctx = run_kit(cfg, base_dir=base_dir)
    if not args.quiet:
        out = os.path.join(cfg.get("output_dir", "out"))
        print("provider        : %s (%s)" % (ctx["provider"].name, ctx["provider"].kind()))
        print("parameters      : %s" % ctx["params"])
        print("tangent source  : %s" % ctx["T_source"])
        print("R^(1) shape     : %s" % (np.asarray(ctx["R1"]).shape,))
        print("hypercomplex ok : %s" % ctx["rhs_diag"].get("hypercomplex_ready"))
        if ctx["U1"] is not None:
            for i, pn in enumerate(ctx["params"]):
                col = ctx["U1"][:, i]
                scal = col[0] if col.size == 1 else float(np.linalg.norm(col))
                print("  dU/d(%-10s) = %+.6e" % (pn, scal))
        print("validation      : %s" % _val.overall(ctx["checks"]))
        print("outputs         : %s/{sensitivity_package,public_report}" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
