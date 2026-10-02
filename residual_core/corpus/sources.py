"""Where a verified corpus case's source, constants and family come from.

Nothing here is invented. The registry says which cases reached
``fully_verified``; the pass that verified them (``store_verification.jsonl``)
carries the manifest Abaqus ran -- PROPS, NSTATV, kinematics -- read from the
author's own deck; the reviewed family classification names the material
family. Paths default to the workspace layout and can be redirected with
``CORPUS_WORKSPACE``.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

__all__ = ["CorpusCase", "CorpusPaths", "load_case", "verified_keys", "paths"]


@dataclass(frozen=True)
class CorpusPaths:
    workspace: Path

    @property
    def registry(self) -> Path:
        return self.workspace / "final-umat" / "paper_results" / "corpus" / "corpus_registry.json"

    @property
    def verification(self) -> Path:
        return self.workspace / "corpus_run" / "pass16" / "results" / "store_verification.jsonl"

    @property
    def families(self) -> Path:
        return self.workspace / "corpus_run" / "material_families_checked_E.json"

    @property
    def discovery(self) -> Path:
        return self.workspace / "discovery_cache"

    @property
    def umat_repo(self) -> Path:
        return Path(os.environ.get("UMAT_OTI_REPO", str(self.workspace / "final-umat")))

    def available(self) -> bool:
        return self.registry.is_file() and self.verification.is_file() and self.discovery.is_dir()


def paths() -> CorpusPaths:
    return CorpusPaths(Path(os.environ.get("CORPUS_WORKSPACE", "/home/ammslab3/softwarex_work")))


@dataclass
class CorpusCase:
    key: str
    source_id: str
    source_path: Path
    source_form: str
    family: str
    terminal_state: str
    verification_fingerprint: str
    kinematics: str                 # 'small strain' | 'finite'
    ntens: int
    nstatv: int
    props: List[float]
    material_provenance: str
    element_type: str
    manifest_notes: str = ""
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def finite(self) -> bool:
        return self.kinematics.strip().lower().startswith("finite")


@lru_cache(maxsize=4)
def _registry(path: str) -> Dict[str, Dict[str, Any]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {record["key"]: record for record in data["records"]}


@lru_cache(maxsize=4)
def _verification(path: str) -> Dict[str, Dict[str, Any]]:
    out = {}
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("key"):
                out[record["key"]] = {"manifest": record.get("manifest", {}),
                                      "kinematics": record.get("kinematics"),
                                      "fingerprint": record.get("fingerprint")}
    return out


@lru_cache(maxsize=4)
def _families(path: str) -> Dict[str, str]:
    if not Path(path).is_file():
        return {}
    rows = json.loads(Path(path).read_text(encoding="utf-8")).get("rows", [])
    return {row["source_id"]: row.get("family", "") for row in rows}


def verified_keys(where: Optional[CorpusPaths] = None) -> List[str]:
    where = where or paths()
    return sorted(key for key, record in _registry(str(where.registry)).items()
                  if record.get("terminal_state") == "fully_verified")


def load_case(key: str, where: Optional[CorpusPaths] = None) -> CorpusCase:
    where = where or paths()
    registry = _registry(str(where.registry))
    if key not in registry:
        raise KeyError("no registry record with key %s" % key)
    record = registry[key]
    run = _verification(str(where.verification)).get(key)
    if run is None:
        raise KeyError("%s has no record in %s" % (key, where.verification))
    manifest = run["manifest"]
    source = where.discovery / record["cache_path"]
    if not source.is_file():
        raise FileNotFoundError("source not in the acquisition cache: %s" % source)
    return CorpusCase(
        key=key, source_id=record["source_id"], source_path=source,
        source_form=manifest.get("source_form") or record.get("source_form") or "fixed",
        family=_families(str(where.families)).get(record["source_id"], ""),
        terminal_state=record.get("terminal_state", ""),
        verification_fingerprint=record.get("verification_fingerprint", ""),
        kinematics=manifest.get("kinematics") or record.get("kinematics") or "",
        ntens=int(manifest.get("ntens") or record.get("ntens") or 0),
        nstatv=int(manifest.get("nstatv") or record.get("nstatv") or 0),
        props=[float(v) for v in manifest.get("props", [])],
        material_provenance=manifest.get("material_provenance", ""),
        element_type=manifest.get("element_type", ""),
        manifest_notes=manifest.get("notes", ""),
        extra={"initial_statev": manifest.get("initial_statev", []),
               "initial_state_from_user_subroutine": bool(
                   manifest.get("initial_state_from_user_subroutine")),
               "isothermal_temperature": manifest.get("isothermal_temperature"),
               "bundle": manifest.get("bundle", []),
               "deck": record.get("deck", "")})
