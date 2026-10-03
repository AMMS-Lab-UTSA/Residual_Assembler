"""Where a verified corpus case's source, constants and family come from.

Nothing here is invented. The registry says which cases reached
``fully_verified``; the pass that verified them (``store_verification.jsonl``)
carries the manifest Abaqus ran -- PROPS, NSTATV, kinematics -- read from the
author's own deck; the reviewed family classification names the material
family. All of it lives outside this repository, in the workspace named by the
environment variable ``CORPUS_WORKSPACE`` -- the folder holding ``final-umat/``
(or set ``UMAT_OTI_REPO``), ``corpus_run/`` and ``discovery_cache/``. There is
no default: reading corpus data with the variable unset raises
:class:`CorpusWorkspaceUnset`.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

__all__ = ["CORPUS_WORKSPACE_ENV", "CorpusCase", "CorpusPaths", "CorpusWorkspaceUnset",
           "key_for_source", "load_case", "verified_keys", "paths"]

CORPUS_WORKSPACE_ENV = "CORPUS_WORKSPACE"


class CorpusWorkspaceUnset(RuntimeError):
    """Corpus data was needed and ``CORPUS_WORKSPACE`` does not name a folder."""


def _unset(what: str) -> CorpusWorkspaceUnset:
    return CorpusWorkspaceUnset(
        "%s is outside this repository; set %s to the workspace folder that holds "
        "final-umat/, corpus_run/ and discovery_cache/" % (what, CORPUS_WORKSPACE_ENV))


@dataclass(frozen=True)
class CorpusPaths:
    workspace: Optional[Path]

    def _root(self, what: str) -> Path:
        if self.workspace is None:
            raise _unset(what)
        return self.workspace

    @property
    def registry(self) -> Path:
        return self._root("the corpus registry") / "final-umat" / "paper_results" / "corpus" / "corpus_registry.json"

    @property
    def verification(self) -> Path:
        return self._root("the corpus verification records") / "corpus_run" / "pass16" / "results" / "store_verification.jsonl"

    @property
    def families(self) -> Path:
        return self._root("the corpus family classification") / "corpus_run" / "material_families_checked_E.json"

    @property
    def discovery(self) -> Path:
        return self._root("the corpus acquisition cache") / "discovery_cache"

    @property
    def umat_repo(self) -> Path:
        explicit = os.environ.get("UMAT_OTI_REPO", "").strip()
        if explicit:
            return Path(explicit).expanduser()
        if self.workspace is None:
            raise CorpusWorkspaceUnset(
                "the UMAT-OTI checkout is unknown; set UMAT_OTI_REPO, or %s to the workspace "
                "folder that holds final-umat/" % CORPUS_WORKSPACE_ENV)
        return self.workspace / "final-umat"

    def available(self) -> bool:
        if self.workspace is None:
            return False
        return self.registry.is_file() and self.verification.is_file() and self.discovery.is_dir()


def paths() -> CorpusPaths:
    """The corpus workspace from ``CORPUS_WORKSPACE`` (``workspace`` is None when unset)."""
    value = os.environ.get(CORPUS_WORKSPACE_ENV, "").strip()
    return CorpusPaths(Path(value).expanduser() if value else None)


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


def key_for_source(source_id: str, sha256: str,
                   where: Optional[CorpusPaths] = None) -> str:
    """The registry key of the source with this id and content.

    A store key is derived from the transform fingerprint, so every re-freeze
    renames it. The source id and the sha256 of the acquired file do not
    change, so a caller that means one particular material names it by those.
    """
    where = where or paths()
    keys = [key for key, record in _registry(str(where.registry)).items()
            if record.get("source_id") == source_id and record.get("sha256") == sha256]
    if len(keys) != 1:
        raise KeyError("expected one registry record for %s at sha256 %s, found %d"
                       % (source_id, sha256, len(keys)))
    return keys[0]


def load_case(key: str, where: Optional[CorpusPaths] = None) -> CorpusCase:
    where = where or paths()
    registry = _registry(str(where.registry))
    if key not in registry:
        raise KeyError("no registry record with key %s" % key)
    record = registry[key]
    # The registry names the run its verdict came from; keys are only valid in
    # that run, so read the manifest there, not in a fixed older pass.
    verification = where.verification
    if record.get("verification_source"):
        verification = (where._root("the corpus verification records") / "corpus_run"
                        / record["verification_source"])
    run = _verification(str(verification)).get(key)
    if run is None:
        raise KeyError("%s has no record in %s" % (key, verification))
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
