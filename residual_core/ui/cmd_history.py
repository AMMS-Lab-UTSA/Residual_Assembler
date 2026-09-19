"""``resasm history``: the collaborator workflow of slide 12 for any provider.

    resasm history --model Analysis.inp --odb Analysis.odb \\
        --material OTI_UMAT.obj --request sensitivity_request.json --out results

(``--fields fields.npz`` replaces ``--odb`` when the ODB was exported
already.) Writes the public ``sensitivity_results.json``,
``sensitivity_tables.csv`` and ``run_report.txt`` (plus
``sensitivity_shares.csv`` and ``fields.npz`` when the request asks for them)
and keeps the link directory, the ODB export and the run details under
``results/private``.

``route_request`` is the dispatcher for ``resasm request``: models inside the
bounded presentation scope (pinned m3_j2 provider, zero-valued boundaries,
concentrated loads, the four-key request) stay on that path; everything else
goes to this engine. See docs/REPLAY_HISTORY.md for the one-line wiring.
"""
from __future__ import annotations

import hashlib
import json
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from ..replay.history import HistoryEngine, ReplayMismatch, run_history
from ..replay.history_inputs import UnsupportedFeature, load_recorded_fields, read_history_model
from ..replay.history_material import HistoryMaterial, ProviderError
from ..replay.history_outputs import (PUBLIC_FILES, Fields, RequestError, load_request,
                                      render_report, validate_request, write_outputs)
from ..replay.history_verify import summarize_fd, tangent_check, whole_model_fd

EXPORTER = Path(__file__).resolve().parents[1] / "replay" / "odb_export_npz.py"

#: Adjacent FD steps must agree to this (relative) for the reference to count
#: as resolved: "Reference resolved: yes". Only a resolved reference can verify.
PLATEAU_SPREAD_LIMIT = 1e-4


def register(subparsers):
    parser = subparsers.add_parser(
        "history", help="history replay for any provider: sensitivities from Analysis.inp + "
                        "Analysis.odb + OTI_UMAT.obj + sensitivity_request.json")
    parser.add_argument("--model", required=True, type=Path, help="Analysis.inp")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--odb", type=Path, help="Analysis.odb (exported with Abaqus Python)")
    source.add_argument("--fields", type=Path, help="an existing odb_export_npz.py export")
    parser.add_argument("--material", required=True, type=Path, help="OTI_UMAT.obj")
    parser.add_argument("--mapping", type=Path, help="completed contract (default: beside the object)")
    parser.add_argument("--request", required=True, type=Path, help="sensitivity_request.json")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--reequilibrate", action="store_true",
                        help="Newton-polish every recorded increment to double-precision equilibrium")
    parser.add_argument("--verify", choices=("none", "tangent", "fd"), default="none",
                        help="independent checks with the ORIGINAL UMAT: tangent spot check, or "
                             "tangent + whole-model central finite differences (re-equilibrated in Python)")
    parser.add_argument("--fd-steps", default="1e-3,3e-4,1e-4,3e-5,1e-5", help="relative FD step ladder")
    parser.add_argument("--abaqus", default="abaqus", help="Abaqus launcher for the ODB export")
    parser.set_defaults(func=run)


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def find_mapping(material: Path, mapping=None):
    candidates = [Path(mapping)] if mapping else [
        path for path in (material.with_suffix(".json"), material.parent / "Mapping.json") if path.is_file()]
    if not candidates:
        raise ProviderError("missing completed contract: put <object>.json or Mapping.json beside the "
                            "object, or pass --mapping")
    contracts = [(path, json.loads(path.read_text())) for path in candidates]
    if any(contract != contracts[0][1] for _, contract in contracts[1:]):
        raise ProviderError("the sidecar contracts beside the object differ; choose one with --mapping")
    return contracts[0]


def export_odb(odb: Path, destination: Path, abaqus="abaqus"):
    executable = shutil.which(str(abaqus))
    if executable is None:
        raise ValueError("Abaqus launcher %r not found: exporting Analysis.odb needs Abaqus Python "
                         "(or pass --fields with an existing export)" % str(abaqus))
    command = [executable, "python", str(EXPORTER), "--", str(odb.resolve()), str(destination.resolve())]
    completed = subprocess.run(command, capture_output=True, text=True, cwd=str(destination.parent))
    (destination.parent / "odb_export.log").write_text(completed.stdout + completed.stderr)
    if completed.returncode or not destination.is_file():
        raise ValueError("ODB export failed (exit %s): %s\n%s" % (
            completed.returncode, " ".join(command), completed.stdout + completed.stderr))
    return command


def _report_skeleton(command):
    return {"status": "failed", "report_fields": {
        "command executed": "yes: %s" % command, "residual assembled": "no",
        "equilibrium checked": "no", "equilibrium passed": "no", "tangent available": "no",
        "tangent verified": "not run", "derivative calculated": "no", "derivative verified": "not run",
        "reference resolved": "not applicable (no reference was run)",
        "Abaqus comparison available": "no", "unsupported feature detected": "none",
        # until write_outputs names what it wrote, only the report is public
        "public and private outputs separated": "yes: public %s; private private/" % [PUBLIC_FILES[2]]},
        "details": {}, "scope": {}, "metadata": {}}


def run_history_request(*, model, material, request, out, odb=None, fields=None, mapping=None,
                        reequilibrate=False, verify="none", fd_steps=(1e-3, 3e-4, 1e-4, 3e-5, 1e-5),
                        abaqus="abaqus", command=None):
    model, material, request, out = map(Path, (model, material, request, out))
    command = command or "resasm history"
    if out.exists() and any(out.iterdir()):
        raise ValueError("output directory must be empty or new; refusing to mix results")
    private = out / "private"
    private.mkdir(parents=True, exist_ok=True)
    report = _report_skeleton(command)
    started = time.perf_counter()
    try:
        for path, suffix in ((model, ".inp"), (material, ".obj"), (request, ".json")) + (
                ((Path(odb), ".odb"),) if odb else ((Path(fields), ".npz"),)):
            if path.suffix.lower() != suffix or not path.is_file():
                raise ValueError("required %s input does not exist or has the wrong suffix: %s" % (suffix, path))
        mapping_path, contract = find_mapping(material, mapping)
        deck = read_history_model(model)
        if odb:
            export_command = export_odb(Path(odb), private / "fields.npz", abaqus)
            fields_path = private / "fields.npz"
        else:
            export_command, fields_path = None, Path(fields)
        recorded = load_recorded_fields(fields_path, deck)
        provider = HistoryMaterial(material, contract, str(private / "link"))
        specification = load_request(request)
        parameters = validate_request(specification, provider.params, provider.nstatev)
        engine = HistoryEngine(deck, provider)
        report["report_fields"]["residual assembled"] = (
            "yes: C3D8 selective-reduced (B-bar), %d elements, %d integration points, %d DOF, "
            "%d increments, sparse assembly" % (engine.ne, engine.ne * 8, engine.ndof, recorded.increments))
        result = run_history(engine, fields=recorded, reequilibrate=reequilibrate)
        worst = max(result.increments, key=lambda r: r.residual_free_max / r.residual_limit)
        parity = {key: max(r.parity[key] for r in result.increments)
                  for key in result.increments[0].parity}

        def ratio(name):
            return max(r.parity["%s_max_abs" % name] / r.parity["%s_limit" % name] for r in result.increments)
        report["report_fields"].update({
            "equilibrium checked": "yes: free-DOF residual of the %s state at all %d increments"
                                   % ("re-equilibrated" if reequilibrate else "recorded", len(result.increments)),
            "equilibrium passed": "yes: max|R_free| = %.3e N at increment %d (limit %.3e N there)"
                                  % (worst.residual_free_max, worst.number, worst.residual_limit),
            "tangent available": "yes: DDSDDE = dSTRESS/dDSTRAN from the provider's OTI strain directions",
            "derivative calculated": "yes: total-history du/dp, dRF/dp, dS/dp, dSDV/dp, dMISES/dp for %d "
                                     "parameters %s at %d increments" % (len(result.parameters),
                                                                         result.parameters, len(result.increments)),
            "Abaqus comparison available": (
                "yes (primal): replayed S, SDV and RF reproduce the ODB at every integration point and "
                "increment; largest error/limit ratios: stress %.3f (max |dS| %.3e MPa), reaction %.3f "
                "(max |dRF| %.3e N)%s; no Abaqus derivative reference is part of this request" % (
                    ratio("stress"), parity["stress_max_abs"], ratio("reaction"), parity["reaction_max_abs"],
                    ", state %.3f (max |dSDV| %.3e)" % (ratio("state"), parity["state_max_abs"])
                    if "state_max_abs" in parity else ""))})
        verification = {}
        derivative_verified = False
        if verify in ("tangent", "fd"):
            check = tangent_check(engine, result, sorted({1, len(result.increments) // 2 or 1,
                                                          len(result.increments)}))
            verification["tangent"] = check
            report["report_fields"]["tangent verified"] = (
                "%s: max relative error %.2e vs central FD of the ORIGINAL UMAT at %d points "
                "(FD plateau spread %.2e)" % ("yes" if check["max_relative_error"] < 1e-5 else "NO",
                                                check["max_relative_error"], check["points_checked"],
                                                check["fd_plateau_spread"]))
        if verify == "fd":
            fd = whole_model_fd(engine, recorded.time, parameters=parameters, steps=tuple(fd_steps))
            summary = summarize_fd(fd)
            live = [(name, key, row) for name, rows in summary.items() for key, row in rows.items()
                    if row["nonzero_increments"]]
            worst_error = max((row["max_error"] for _, _, row in live), default=0.0)
            worst_spread = max((row["max_spread"] for _, _, row in live), default=0.0)
            zero_oti = max((row["zero_reference_max_weighted_error"] for rows in summary.values()
                            for row in rows.values()), default=0.0)
            reference = fd.pop("reference_result")
            replay_vs_solve = max(
                float(np.abs(a.du - b.du).max() / max(np.abs(b.du).max(), 1e-300))
                for a, b in zip(result.increments, reference.increments))
            verification["whole_model_fd"] = {"summary": summary, "steps": fd["steps"],
                                              "definition": fd["definition"], "seconds": fd["seconds"],
                                              "replay_vs_python_equilibrium_du_rel": replay_vs_solve}
            # an unresolved reference (no plateau) verifies nothing, however
            # close OTI and FD happen to be
            resolved = worst_spread < PLATEAU_SPREAD_LIMIT
            derivative_verified = (resolved and worst_error <= max(1e-6, 2 * worst_spread)
                                   and zero_oti <= 1e-6)
            details = (
                "whole-model central FD of the ORIGINAL UMAT re-equilibrated in Python; worst "
                "nonzero-derivative error %.2e (plateau spread %.2e); zero references: |OTI - FD| <= %.1e "
                "on the field scale; the ODB-driven du/dp differs from the Python-equilibrium du/dp by "
                "%.2e (relative)" % (worst_error, worst_spread, zero_oti, replay_vs_solve))
            report["report_fields"]["derivative verified"] = (
                "not verified: the reference did not resolve (largest plateau spread %.2e >= %.0e); %s"
                % (worst_spread, PLATEAU_SPREAD_LIMIT, details) if not resolved else
                "%s: %s" % ("yes" if derivative_verified else "NO", details))
            report["report_fields"]["reference resolved"] = (
                "yes: a plateau (adjacent steps of %s agreeing to %.1e) for every nonzero derivative"
                % (fd["steps"], worst_spread) if resolved else
                "partially: largest plateau spread %.2e" % worst_spread)
        results_fields = Fields(engine, result)
        timings = {key: round(value, 3) for key, value in result.timings.items()}
        timings["wall_total"] = round(time.perf_counter() - started, 3)
        report["status"] = "executed successfully"
        report["scope"] = {"element_type": "C3D8", "integration": result.integration,
                           "elements": engine.ne, "integration_points": engine.ne * 8,
                           "dof": engine.ndof, "increments_replayed": len(result.increments),
                           "parameters": result.parameters, "requested_parameters": parameters,
                           "mode": result.mode, "provider": contract.get("model_id"),
                           "kinematics": "small_strain", "step": deck.step_name}
        report["metadata"] = {
            "sensitivity_semantics": "total equilibrated history (dSTRESS, dSTATEV chained through "
                                     "UMAT_OTI_EVAL_TOTAL; du_c/dp = 0 on prescribed DOFs)",
            # True only when the derivatives passed the whole-model FD check
            "verified": derivative_verified, "verification": verification,
            "max_scaled_free_residual": result.max_scaled_residual,
            "parity": parity, "tolerances": result.tolerances, "timings_s": timings,
            "parameter_values": dict(zip(result.parameters, result.parameter_values.tolist())),
            "input_sha256": {"Analysis.inp": _digest(model), "OTI_UMAT.obj": provider.object_sha256,
                             "contract": _digest(mapping_path), "sensitivity_request.json": _digest(request),
                             ("Analysis.odb" if odb else "fields.npz"): _digest(odb or fields)},
            "regular_source_hash": contract.get("regular_source_hash")}
        report["details"] = {
            "Mode": result.mode, "Parameters (provider order)": result.parameters,
            "Timings (s)": timings, "Tolerances": result.tolerances,
            "Production analysis rerun": "no", "Material source read or transformed": "no"}
        write_outputs(out, fields=results_fields, request=specification, parameters=parameters,
                      report=report, private={"export_command": export_command, "fields": str(fields_path),
                                              "increments": [{
                                                  "number": r.number, "time": r.time,
                                                  "iterations": r.iterations, "correction": r.correction,
                                                  "residual_free_max": r.residual_free_max,
                                                  "residual_limit": r.residual_limit,
                                                  "sensitivity_residual": r.sensitivity_residual,
                                                  "parity": r.parity} for r in result.increments]})
        return report
    except UnsupportedFeature as error:
        report["report_fields"]["unsupported feature detected"] = "yes: %s" % error
        report["details"]["Error"] = str(error)
        (out / PUBLIC_FILES[2]).write_text(render_report(report))
        raise
    except (ReplayMismatch, RequestError, ProviderError, ValueError, OSError, KeyError) as error:
        report["details"]["Error"] = "%s: %s" % (type(error).__name__, error)
        (out / PUBLIC_FILES[2]).write_text(render_report(report))
        raise ValueError(str(error)) from error


def _command_line(name, args, keys):
    """The command as parsed (not sys.argv, which is the host process when called in-process)."""
    words = ["resasm", name]
    for key in keys:
        value = getattr(args, key, None)
        if value is None or value is False:
            continue
        flag = "--" + key.replace("_", "-")
        words += [flag] if value is True else [flag, str(value)]
    return " ".join(shlex.quote(word) for word in words)


def run(args):
    command = _command_line("history", args, ("model", "odb", "fields", "material", "mapping", "request",
                                              "out", "reequilibrate", "verify", "fd_steps", "abaqus"))
    try:
        report = run_history_request(
            model=args.model, material=args.material, request=args.request, out=args.out,
            odb=args.odb, fields=args.fields, mapping=args.mapping, reequilibrate=args.reequilibrate,
            verify=args.verify, fd_steps=tuple(float(x) for x in args.fd_steps.split(",")),
            abaqus=args.abaqus, command=command)
    except (ValueError, OSError) as error:
        print("history replay failed: %s" % error, file=sys.stderr)
        return 2
    scope = report["scope"]
    print("history replay: %d increments, %d integration points, %d parameters; %s"
          % (scope["increments_replayed"], scope["integration_points"], len(scope["parameters"]),
             report["report_fields"]["equilibrium passed"]))
    print("results: %s" % (Path(args.out) / PUBLIC_FILES[0]))
    return 0


# ----------------------------------------------------------------- routing
def bounded_scope_reason(args):
    """None when the bounded presentation engine can run the request, else why not."""
    from ..replay.presentation import NotThePinnedProvider, mapping_for
    from ..replay.presentation_inputs import read_model
    try:
        model = read_model(Path(args.model))
        mapping_for(Path(args.material), getattr(args, "mapping", None))
    except NotThePinnedProvider as error:
        # a readable deck and a valid mapping of another provider: the bounded
        # engine would only refuse the material, the history engine takes any
        return str(error)
    except (ValueError, OSError, KeyError):
        # Inputs the bounded engine cannot read are its to diagnose: it owns
        # the presentation interface's categorised, private-by-default error
        # reports, so an unreadable deck or a missing mapping is reported by it
        # rather than being forwarded to the history engine as a scope reason.
        try:
            read_model(Path(args.model))
        except (ValueError, OSError, KeyError) as error:
            if "unsupported" in str(error).lower() and Path(args.model).is_file():
                return str(error)
        return None
    if any(boundary.value != 0 for boundary in model.boundaries):
        return "nonzero prescribed displacements"
    if not model.cloads:
        return "no concentrated loads (displacement-driven history)"
    try:
        request = json.loads(Path(args.request).read_text())
    except (ValueError, OSError) as error:
        return "request: %s" % error
    if set(request) != {"outputs", "parameters", "domain", "increments"}:
        return "request uses keys beyond outputs/parameters/domain/increments"
    for output in request.get("outputs", []):
        if output.get("field") not in ("U", "RF", "S", "SDV") or output.get("reduction") not in (
                "component", "sum", "mean", "L2", "max"):
            return "output %r needs the history engine (field/reduction)" % output.get("name")
        if set(output.get("domain", {})) - {"nodes", "elements"}:
            return "output %r uses set or point domains" % output.get("name")
    if set(request.get("domain", {})) - {"nodes", "elements"}:
        return "request domain uses sets or points"
    return None


def route_request(args):
    """Dispatch ``resasm request`` to the bounded engine or to the history engine."""
    reason = bounded_scope_reason(args)
    if reason is None:
        from .cmd_request import run as bounded
        return bounded(args)
    print("resasm request: outside the bounded presentation scope (%s); using the history "
          "replay engine" % reason)
    command = _command_line("request", args, ("model", "odb", "material", "mapping", "request", "out",
                                              "abaqus", "validate"))
    try:
        report = run_history_request(
            model=args.model, material=args.material, request=args.request, out=args.out,
            odb=args.odb, mapping=getattr(args, "mapping", None),
            verify="fd" if getattr(args, "validate", False) else "none",
            abaqus=getattr(args, "abaqus", "abaqus"), command=command)
    except (ValueError, OSError) as error:
        print("request failed: %s" % error, file=sys.stderr)
        return 2
    print("request executed: history engine, %d increments; verified=%s"
          % (report["scope"]["increments_replayed"], report["metadata"]["verified"]))
    return 0
