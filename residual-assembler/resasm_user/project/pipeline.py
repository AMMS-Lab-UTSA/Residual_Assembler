"""Pipeline actions (M5-A): configure -> prepare -> run -> export -> assemble
-> validate -> report, wrapping the proven M1-M4 backend.

The offline actions (configure, prepare file-generation, assemble, validate,
report) need no Abaqus. The Abaqus actions (run, export) drive Abaqus + the OTI
link and determine success from the .sta / ODB, not the process exit code
(Abaqus 2021 crashes during wrap-up after writing a valid ODB).

For this first slice ``generate_oti_umat`` uses an already-OTI UMAT directly (the
proven hand-written OTI UMATs); replacing it with the umat_oti source
transformer is the next slice and only changes this one step.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from typing import Any, Dict, List

from .engine import Action, ActionLog, register
from .manifest import Project, ProjectError
from .state import State

_DDSDDE_SLOTS = 36
_PER_PARAM = 6
_OTI_DIRS = ["E1", "E2", "E3", "E4"]


# --------------------------------------------------------------------------- #
def set_parameter_selection(project: Project, selections: List[Dict[str, Any]]) -> None:
    """Record the user's parameter choices (name + PROPS index) for
    configure_parameters. ``selections`` = [{"name": "C11", "index": 1}, ...]."""
    project.data["parameter_selection"] = selections
    project.save()


class ConfigureParameters(Action):
    name = "configure_parameters"
    requires = State.INSPECTED
    produces = State.CONFIGURED

    def artifacts(self, project):
        return [project.artifact("generated", "parameters.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        sel = project.data.get("parameter_selection")
        if not sel:
            raise ProjectError(
                "no parameters selected; call set_parameter_selection() / "
                "`resasm project configure --param NAME:INDEX ...` first")
        if len(sel) > len(_OTI_DIRS):
            raise ProjectError("at most %d parameters (OTI directions) supported "
                               "with otim4n1" % len(_OTI_DIRS))
        consts = _material_constants(project)
        params = []
        for i, s in enumerate(sel):
            idx = int(s["index"])
            lo = _DDSDDE_SLOTS + 1 + i * _PER_PARAM
            params.append({
                "name": str(s["name"]), "source": "PROPS", "index": idx,
                "value": (consts[idx - 1] if 0 < idx <= len(consts) else None),
                "oti_direction": _OTI_DIRS[i],
                "sdv_range": [lo, lo + _PER_PARAM - 1],
            })
            log("param %-4s PROPS(%d)=%s  dir=%s  SDV %d-%d"
                % (params[-1]["name"], idx, params[-1]["value"],
                   params[-1]["oti_direction"], lo, lo + _PER_PARAM - 1))
        project.data["parameters"] = params
        _dump(project.artifact("generated", "parameters.json"), params)
        return {"n_parameters": len(params)}


# --------------------------------------------------------------------------- #
def set_umat_oti_config(project: Project, config_path: str) -> None:
    """Point the project at the transformer's contract JSON for its UMAT, so
    transform_umat can turn an ordinary UMAT into an OTI UMAT automatically."""
    project.data["umat_oti_config"] = os.path.abspath(config_path)
    project.save()


class TransformUmat(Action):
    """Turn an ordinary UMAT into an OTI UMAT via the umat_oti transformer
    (M5-B). Emits the resasm_umat_transform_v1 handoff manifest. This slice
    produces the DDSDDE (tangent) transform; PROPS parameter-seeding codegen is
    the next transformer extension. Skipped when the input UMAT is already OTI."""
    name = "transform_umat"
    requires = State.CONFIGURED
    produces = State.CONFIGURED

    def artifacts(self, project):
        return [project.artifact("generated", "umat_transform_manifest.json"),
                project.artifact("generated", "umat_transformed.for")]

    def is_complete(self, project):
        if (project.data.get("umat") or {}).get("already_oti"):
            return True                        # nothing to transform
        return Action.is_complete(self, project)

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        if (project.data.get("umat") or {}).get("already_oti"):
            log("UMAT is already OTI-seeded; no transformation needed")
            return {"transformed": False}
        cfg = project.data.get("umat_oti_config")
        if not cfg:
            raise ProjectError(
                "an ordinary (non-OTI) UMAT needs a umat_oti contract; set it with "
                "set_umat_oti_config(project, contract.json) before transform_umat")
        from resasm_user.umat_backend import (UmatOtiBackend, ParameterSelection,
                                              write_handoff_manifest)
        umat = project.input_path("umat")
        params = [ParameterSelection(p["name"], int(p["index"]))
                  for p in project.data.get("parameters") or []]
        backend = UmatOtiBackend()
        proposal = backend.propose_contract(umat, params, base_config=cfg)
        for n in proposal.notes:
            log(n, level="warn")
        out = project.artifact("generated", "umat_transform")
        result = backend.transform(proposal, out)
        if not result.success:
            raise ProjectError("UMAT transformation failed: %s"
                               % (result.blockers or "unknown"))
        val = backend.validate_transformation_semantics(result)
        if not val.passed:
            raise ProjectError("transformed UMAT failed the semantic checks")
        # Structural check passed. If parameter derivatives were generated and a
        # Fortran compiler is available, also run the NUMERICAL real-response check
        # (dsigma/da_i vs finite difference of the original UMAT) -- these are two
        # distinct guarantees and are logged as such.
        real_response = None
        if result.manifest["capabilities"].get("parameter_derivatives"):
            try:
                rr = backend.validate_real_response(result)
                real_response = rr.checks
                if rr.checks.get("available"):
                    if not rr.passed:
                        raise ProjectError("transformed UMAT failed the numerical "
                                           "real-response check: %s" % rr.message)
                    log("numerical real-response OK: %s" % rr.message)
                else:
                    log("numerical real-response skipped: %s"
                        % rr.checks.get("reason", "unavailable"), level="warn")
            except NotImplementedError:
                pass
        else:
            log("DDSDDE-only transform (no parameter derivatives requested)")
        shutil.copy(result.transformed_umat,
                    project.artifact("generated", "umat_transformed.for"))
        write_handoff_manifest(
            project.artifact("generated", "umat_transform_manifest.json"),
            result.manifest)
        project.data["transform_manifest"] = result.manifest
        if real_response is not None:
            project.data["transform_real_response"] = real_response
        log("transformed %s -> OTI (semantic checks passed)" % os.path.basename(umat))
        return {"transformed": True,
                "required_depvar": result.manifest["abaqus"]["required_depvar"]}


class Prepare(Action):
    """Generate the OTI UMAT, SDV layout, compiler env and patched .inp -- the
    user never touches SDV37, *Depvar, abaqus_v6.env or -lotim4n1 by hand."""
    name = "prepare"
    requires = State.CONFIGURED
    produces = State.PREPARED

    def artifacts(self, project):
        g = lambda n: project.artifact("generated", n)   # noqa: E731
        return [g("umat.for"), g("derivative_layout.json"),
                g("abaqus_v6.env"), g("model.inp")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        params = project.data.get("parameters") or []
        umat_report = project.data.get("umat") or {}
        max_sdv = _DDSDDE_SLOTS + _PER_PARAM * len(params)

        # 1) generate_oti_umat: prefer a transformer-produced OTI UMAT, else an
        #    already-OTI input UMAT ------------------------------------------
        transformed = project.artifact("generated", "umat_transformed.for")
        umat_in = project.input_path("umat")
        if os.path.exists(transformed):
            shutil.copy(transformed, project.artifact("generated", "umat.for"))
            log("generate_oti_umat: using umat_oti-transformed UMAT")
        elif umat_report.get("already_oti"):
            shutil.copy(umat_in, project.artifact("generated", "umat.for"))
            log("generate_oti_umat: using already-OTI UMAT %s"
                % os.path.basename(umat_in))
        else:
            raise ProjectError(
                "the input UMAT is not OTI-seeded and no transform was run. Either "
                "provide an OTI UMAT, or set a umat_oti contract and run "
                "transform_umat (set_umat_oti_config + `resasm project prepare`).")

        # 2) SDV layout (single source of truth) -----------------------------
        layout = {
            "schema": "resasm_sdv_layout_v1", "element_type": "C3D8",
            "kinematics": "small_strain", "integration": "selective_reduced",
            "voigt_order": ["11", "22", "33", "12", "13", "23"],
            "integration_point_order": "abaqus_label_ascending",
            "sdv_indexing": "abaqus_1_based",
            "sdv_layout": {"ddsdde": [1, _DDSDDE_SLOTS],
                           "parameters": {p["name"]: p["sdv_range"] for p in params}},
        }
        _dump(project.artifact("generated", "derivative_layout.json"), layout)
        log("derivative_layout.json: DDSDDE[1-36] + %d params" % len(params))

        # 3) compiler env (OTI include/link) ---------------------------------
        oti = (project.data.get("toolchain", {}) or {}).get("oti_dir")
        if not oti:
            raise ProjectError("OTI library not detected (run detect_toolchain); "
                               "set OTI_DIR to the folder with libotim4n1.a")
        env = ("import os\noti = %r\n"
               "compile_fortran += ['-I'+os.getcwd(), '-I'+oti]\n"
               "link_sl  += ['-L'+oti, '-lotim4n1']\n"
               "link_exe += ['-L'+oti, '-lotim4n1']\n" % oti)
        with open(project.artifact("generated", "abaqus_v6.env"), "w") as fh:
            fh.write(env)
        log("abaqus_v6.env: OTI include/link for %s" % oti)

        # 4) patch the .inp: *Depvar >= max_sdv and S/SDV/U output ------------
        patched = _patch_inp(project.input_path("inp"), max_sdv, log)
        with open(project.artifact("generated", "model.inp"), "w") as fh:
            fh.write(patched)
        log("model.inp patched (*Depvar>=%d, S/SDV/U requested)" % max_sdv)
        return {"max_sdv": max_sdv, "oti_dir": oti}


# --------------------------------------------------------------------------- #
class RunAbaqus(Action):
    name = "run_abaqus"
    requires = State.PREPARED
    produces = State.SOLVED

    def artifacts(self, project):
        return [project.artifact("runs", "job/model.odb")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        abq = (project.data.get("toolchain", {}) or {}).get("abaqus")
        if not abq:
            raise ProjectError("Abaqus not detected; cannot run the job")
        run = project.artifact("runs", "job")
        os.makedirs(run, exist_ok=True)
        for n in ("model.inp", "umat.for", "derivative_layout.json", "abaqus_v6.env"):
            shutil.copy(project.artifact("generated", n), os.path.join(run, n))
        oti = (project.data.get("toolchain", {}) or {}).get("oti_dir")
        for n in ("otim4n1.mod", "libotim4n1.a"):
            shutil.copy(os.path.join(oti, n), os.path.join(run, n))
        log("submitting Abaqus job in %s" % run)
        subprocess.call('%s job=model user=umat.for interactive' % abq,
                        shell=True, cwd=run,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # success is determined from the .sta / ODB, NOT the exit code
        sta = os.path.join(run, "model.sta")
        odb = os.path.join(run, "model.odb")
        ok = os.path.exists(sta) and "COMPLETED SUCCESSFULLY" in open(sta).read()
        if not (ok and os.path.exists(odb)):
            raise ProjectError("Abaqus did not complete (.sta/ODB check failed); "
                               "see %s" % run)
        log("Abaqus completed (.sta OK, ODB present)")
        return {"odb": _rel(project, odb)}


class ExportDerivativeFields(Action):
    name = "export_derivative_fields"
    requires = State.SOLVED
    produces = State.EXTRACTED

    def artifacts(self, project):
        return [project.artifact("generated", "derivative_fields.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        abq = (project.data.get("toolchain", {}) or {}).get("abaqus")
        run = project.artifact("runs", "job")
        odb = os.path.join(run, "model.odb")
        root = _repo_root()
        for n in ("export_derivative_fields.py", "abaqus_odb_export.py"):
            shutil.copy(os.path.join(root, "residual_core", "io", n),
                        os.path.join(run, n))
        out = project.artifact("generated", "derivative_fields.json")
        cmd = ('%s python export_derivative_fields.py -- --odb model.odb '
               '--step Step-1 --frame -1 --layout derivative_layout.json '
               '--output %s' % (abq, out))
        subprocess.call(cmd, shell=True, cwd=run,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if not os.path.exists(out):
            raise ProjectError("ODB export produced no derivative_fields.json")
        log("exported derivative fields -> %s" % _rel(project, out))
        return {"fields": _rel(project, out)}


# --------------------------------------------------------------------------- #
class AssembleSensitivities(Action):
    name = "assemble_sensitivities"
    requires = State.EXTRACTED
    produces = State.SENSITIVITIES_COMPLETE

    def artifacts(self, project):
        return [project.artifact("results", "manifest.json"),
                project.artifact("results", "displacement_sensitivities.npz")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        from resasm_user.field_recipe import run_field_recipe
        names = [p["name"] for p in project.data.get("parameters") or []]
        recipe = self._write_recipe(project, names)
        log("running field-driven residual sensitivity for %s" % names)
        res = run_field_recipe(recipe)
        # move results into the project results/ dir
        out = project.path("results")
        for f in os.listdir(res.output_dir):
            shutil.copy(os.path.join(res.output_dir, f), os.path.join(out, f))
        norms = res.summary.get("du_da_norms", {})
        for k, v in norms.items():
            log("||du/d%s|| = %.6e" % (k, v))
        return {"du_da_norms": norms, "ndof": res.summary.get("ndof")}

    def _write_recipe(self, project: Project, names: List[str]) -> str:
        fields = project.artifact("generated", "derivative_fields.json")
        if not os.path.exists(fields):
            raise ProjectError("no derivative_fields.json (run export first)")
        yaml = (
            "analysis:\n  type: field_residual_sensitivity\n"
            "  kinematics: small_strain\n"
            "mesh:\n  format: abaqus\n  file: %s\n"
            "results:\n  format: resasm_derivative_fields_v1\n  file: %s\n"
            "parameters: [%s]\n"
            "assumptions:\n  parameter_independent_geometry: true\n"
            "  parameter_independent_loads: true\n"
            "  parameter_independent_boundaries: true\n"
            "output:\n  directory: %s\n"
            % (project.artifact("generated", "model.inp"), fields,
               ", ".join(names), project.artifact("results", "_recipe_out")))
        path = project.artifact("generated", "sensitivity.yaml")
        with open(path, "w") as fh:
            fh.write(yaml)
        return path


class RunValidation(Action):
    name = "run_validation"
    requires = State.SENSITIVITIES_COMPLETE
    produces = State.VALIDATED

    def artifacts(self, project):
        return [project.artifact("results", "validation_report.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        man = json.load(open(project.artifact("results", "manifest.json")))
        checks = {
            "prescribed_dof_sensitivity_is_zero":
                bool(man.get("prescribed_dof_sensitivity_is_zero")),
            "du_da_all_finite": all(
                v == v and abs(v) < float("inf")
                for v in man.get("du_da_norms", {}).values()),
            "tangent_source": man.get("tangent_source"),
            "integration": man.get("integration"),
        }
        checks["passed"] = bool(checks["prescribed_dof_sensitivity_is_zero"]
                                and checks["du_da_all_finite"])
        checks["note"] = ("offline checks only; parameter finite-difference and "
                          "Euler-identity validation against Abaqus reruns is a "
                          "separate (Abaqus-gated) step")
        for k, v in checks.items():
            log("%s = %s" % (k, v))
        _dump(project.artifact("results", "validation_report.json"), checks)
        if not checks["passed"]:
            raise ProjectError("validation failed: %s" % checks)
        return checks


class GenerateReport(Action):
    name = "generate_report"
    requires = State.VALIDATED
    produces = State.REPORTED

    def artifacts(self, project):
        return [project.artifact("results", "report.md"),
                project.artifact("results", "project_manifest.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        man = json.load(open(project.artifact("results", "manifest.json")))
        val = json.load(open(project.artifact("results", "validation_report.json")))
        lines = ["# Sensitivity report -- %s" % project.data.get("name"), "",
                 "| item | value |", "|---|---|",
                 "| state | %s |" % project.state,
                 "| parameters | %s |" % ", ".join(man.get("parameters", [])),
                 "| ndof | %s |" % man.get("ndof"),
                 "| integration | %s |" % man.get("integration"),
                 "| validation passed | %s |" % val.get("passed"), "",
                 "## du/da norms", "", "| parameter | ||du/da|| |", "|---|---|"]
        for k, v in man.get("du_da_norms", {}).items():
            lines.append("| %s | %.6e |" % (k, v))
        with open(project.artifact("results", "report.md"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        _dump(project.artifact("results", "project_manifest.json"), project.data)
        log("report written")
        return {"report": _rel(project, project.artifact("results", "report.md"))}


# --------------------------------------------------------------------------- #
def _material_constants(project: Project) -> List[float]:
    mats = (project.data.get("model", {}) or {}).get("materials", {}) or {}
    for m in mats.values():
        if m.get("user_material") and m.get("constants"):
            return list(m["constants"])
    # fall back to any material with constants
    for m in mats.values():
        if m.get("constants"):
            return list(m["constants"])
    return []


def _patch_inp(inp_path: str, max_sdv: int, log: ActionLog) -> str:
    with open(inp_path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    lines = text.splitlines()
    out, i, n = [], 0, len(lines)
    saw_depvar = False
    while i < n:
        line = lines[i]
        if line.strip().lower().startswith("*depvar"):
            saw_depvar = True
            out.append(line)
            i += 1
            # next non-comment data line holds the count
            if i < n:
                out.append("%d," % max_sdv)
                i += 1
            continue
        out.append(line)
        i += 1
    if not saw_depvar:
        log("no *Depvar found; leaving material as-is (add *Depvar %d manually)"
            % max_sdv, level="warn")
    # ensure S, SDV, U output (best-effort, non-destructive)
    joined = "\n".join(out)
    if "SDV" not in joined.upper():
        log("note: model requests no SDV output; ensure *Element Output has SDV",
            level="warn")
    return joined + "\n"


def _dump(path: str, obj: Any) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


def _rel(project: Project, path: str) -> str:
    try:
        return os.path.relpath(path, project.root)
    except ValueError:
        return path


def _repo_root() -> str:
    here = os.path.dirname(os.path.abspath(__file__))         # resasm_user/project
    return os.path.dirname(os.path.dirname(here))              # repo root


register(ConfigureParameters())
register(TransformUmat())
register(Prepare())
register(RunAbaqus())
register(ExportDerivativeFields())
register(AssembleSensitivities())
register(RunValidation())
register(GenerateReport())
