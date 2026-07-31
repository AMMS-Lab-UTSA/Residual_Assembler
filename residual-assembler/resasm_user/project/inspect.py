"""Inspection actions (M5-A): scan the model, the UMAT and the toolchain.

These produce the readiness report -- element types, materials, PROPS values and
count, *Depvar, existing SDV use, supported/unsupported model features, whether a
C3D8 model needs the Abaqus B-bar path, and whether Abaqus / Intel Fortran / the
OTI library are available. Everything here is offline (no Abaqus required).
"""

from __future__ import annotations

import json
import os
import re
import shutil
from typing import Any, Dict

from .engine import Action, ActionLog, register
from .manifest import Project, ProjectError
from .state import State

_C3D8 = "C3D8"


# --------------------------------------------------------------------------- #
class DetectToolchain(Action):
    name = "detect_toolchain"
    requires = State.NEW
    produces = State.NEW

    def artifacts(self, project):
        return [project.artifact("generated", "toolchain.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        abq = os.environ.get("ABAQUS_CMD", "abaqus")
        oti = os.path.expanduser(os.environ.get("OTI_DIR", "~/MultiZ_f/oti"))
        tc = {
            "abaqus": shutil.which(abq),
            "ifort": shutil.which("ifort"),
            "oti_dir": oti if os.path.exists(os.path.join(oti, "libotim4n1.a")) else None,
            "oti_lib": os.path.join(oti, "libotim4n1.a")
            if os.path.exists(os.path.join(oti, "libotim4n1.a")) else None,
        }
        for k, v in tc.items():
            log("%-8s : %s" % (k, v if v else "NOT FOUND"))
        project.data["toolchain"] = tc
        _dump(project.artifact("generated", "toolchain.json"), tc)
        return tc


# --------------------------------------------------------------------------- #
class InspectModel(Action):
    name = "inspect_model"
    requires = State.NEW
    produces = State.NEW

    def artifacts(self, project):
        return [project.artifact("generated", "model_inspection.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        from residual_core.io import abaqus_inp_parser as P
        inp = project.input_path("inp")
        if not inp:
            raise ProjectError("project has no .inp input")
        am = P.parse_inp(inp)

        etypes: Dict[str, int] = {}
        for e in am.elements.values():
            etypes[e.etype.upper()] = etypes.get(e.etype.upper(), 0) + 1

        mats = {}
        for name, m in am.materials.items():
            mats[name] = {
                "user_material": bool(getattr(m, "user_material", False)),
                "depvar": getattr(m, "depvar", None),
                "n_constants": len(list(getattr(m, "constants", []) or [])),
                "constants": [float(c) for c in (getattr(m, "constants", []) or [])],
            }

        only_c3d8 = set(etypes) == {_C3D8}
        unsupported_elems = sorted(set(etypes) - {_C3D8})
        blockers = []
        if not only_c3d8:
            blockers.append("unsupported element type(s): %s (C3D8 only)"
                            % ", ".join(unsupported_elems))
        if getattr(am, "equations", None):
            blockers.append("*Equation constraints present (not supported)")
        if getattr(am, "dsloads", None):
            blockers.append("*Dsload distributed loads present (not supported)")

        report = {
            "element_types": etypes,
            "n_nodes": len(am.nodes), "n_elements": len(am.elements),
            "materials": mats,
            "n_boundaries": len(getattr(am, "boundaries", []) or []),
            "n_cloads": len(getattr(am, "cloads", []) or []),
            "integration": "selective_reduced" if only_c3d8 else "unknown",
            "supported": not blockers,
            "blockers": blockers,
        }
        for line in ("elements: %s" % etypes,
                     "materials: %s" % list(mats),
                     "integration: %s" % report["integration"],
                     "supported: %s" % report["supported"]):
            log(line)
        for b in blockers:
            log("blocker: %s" % b, level="warn")
        project.data["model"] = report
        _dump(project.artifact("generated", "model_inspection.json"), report)
        return report


# --------------------------------------------------------------------------- #
class InspectUmat(Action):
    name = "inspect_umat"
    requires = State.NEW
    produces = State.INSPECTED

    def artifacts(self, project):
        return [project.artifact("generated", "umat_inspection.json")]

    def perform(self, project: Project, log: ActionLog) -> Dict[str, Any]:
        umat = project.input_path("umat")
        if not umat:
            report = {"present": False,
                      "note": "no UMAT provided; a user material needs one"}
            log("no UMAT input", level="warn")
            project.data["umat"] = report
            _dump(project.artifact("generated", "umat_inspection.json"), report)
            return report

        with open(umat, "r", encoding="utf-8", errors="replace") as fh:
            src = fh.read()
        up = src.upper()
        props = sorted(set(int(m) for m in re.findall(r"PROPS\(\s*(\d+)\s*\)", up)))
        statev = sorted(set(int(m) for m in re.findall(r"STATEV\(\s*(\d+)\s*\)", up)))
        subs = re.findall(r"^\s*SUBROUTINE\s+(\w+)", src, re.IGNORECASE | re.MULTILINE)
        report = {
            "present": True,
            "props_refs": props,
            "n_props_referenced": (max(props) if props else 0),
            "statev_refs_max": (max(statev) if statev else 0),
            "writes_ddsdde": "DDSDDE" in up,
            "already_oti": "USE OTIM" in up or "ONUMM" in up,
            "subroutines": subs,
        }
        log("PROPS referenced up to %d; STATEV up to %d; DDSDDE=%s; OTI=%s"
            % (report["n_props_referenced"], report["statev_refs_max"],
               report["writes_ddsdde"], report["already_oti"]))
        project.data["umat"] = report
        _dump(project.artifact("generated", "umat_inspection.json"), report)
        return report


def _dump(path: str, obj: Any) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


register(DetectToolchain())
register(InspectModel())
register(InspectUmat())
