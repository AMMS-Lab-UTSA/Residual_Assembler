"""Material package manifest: the boundary contract for a JHU material.

A material package is ``{material_package.json + oti binary}``. The manifest is
the load-bearing half: it is the ONLY place parameter names, PROPS slots, OTI
directions, state layout, kinematics and units are pinned, so the residual tool
knows exactly what it is differentiating. It also carries the identity of the
matched REGULAR twin, so the tool can refuse to replay a record produced by a
different model.

This module validates the manifest, resolves the binary, checks it against the
loaded ABI descriptor, and enforces the matched-twin rule.
"""

from __future__ import annotations

import hashlib
import json
import os
from typing import Any, Dict, List, Mapping, Optional

from residual_core.runtime import IncompatibleBinaryError, check_binary_compatibility
from .abi import KIN_SMALL_STRAIN, KIN_FINITE_STRAIN, MaterialABI

SCHEMA = "resasm_material_package_v1"
_KIN_CODE = {"small_strain": KIN_SMALL_STRAIN, "finite_strain": KIN_FINITE_STRAIN}
_CONTRACT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "contract")


def expected_contract_version() -> str:
    """The combined hash of THIS tool's copy of the shared contract (ABI header +
    manifest schema), or '' if unavailable. Used to detect contract drift."""
    try:
        with open(os.path.join(_CONTRACT_DIR, "CONTRACT_VERSION.json")) as fh:
            return str(json.load(fh).get("combined_hash", ""))
    except (OSError, ValueError):
        return ""


class PackageError(Exception):
    """The material package manifest is missing/invalid/inconsistent."""


class TwinMismatchError(PackageError):
    """The OTI binary is not the matched twin of the model that produced the
    replay record (different model_id or regular-binary hash). Replaying would
    silently differentiate the wrong constitutive law."""


class MaterialPackage:
    def __init__(self, manifest: Mapping[str, Any], base_dir: str = "."):
        self.raw = dict(manifest)
        self.base_dir = base_dir
        self._validate()

    # -- construction ----------------------------------------------------- #
    @classmethod
    def load(cls, path: str) -> "MaterialPackage":
        try:
            with open(path, "r", encoding="utf-8") as fh:
                manifest = json.load(fh)
        except (OSError, ValueError) as exc:
            raise PackageError("cannot read material package %r: %s" % (path, exc))
        return cls(manifest, base_dir=os.path.dirname(os.path.abspath(path)))

    # -- validation ------------------------------------------------------- #
    _REQUIRED = ("model_id", "kinematics", "ntens", "nprops", "nstatev",
                 "parameters", "abi", "binaries")

    def _validate(self) -> None:
        m = self.raw
        if m.get("schema") != SCHEMA:
            raise PackageError("manifest schema must be %r, got %r"
                               % (SCHEMA, m.get("schema")))
        for k in self._REQUIRED:
            if k not in m:
                raise PackageError("manifest is missing required key %r" % k)
        if m["kinematics"] not in _KIN_CODE:
            raise PackageError("kinematics must be one of %r, got %r"
                               % (sorted(_KIN_CODE), m["kinematics"]))
        # dimensions must be coherent -- a negative nstatev would silently defeat
        # the record preflight's path-dependent checks, etc.
        for key, lo in (("ntens", 1), ("nprops", 1), ("nstatev", 0)):
            try:
                val = int(m[key])
            except (TypeError, ValueError):
                raise PackageError("%s must be an integer, got %r" % (key, m[key]))
            if val < lo:
                raise PackageError("%s must be >= %d, got %d" % (key, lo, val))
        params = m["parameters"]
        if not isinstance(params, list) or not params:
            raise PackageError("parameters must be a non-empty list")
        seen_name, seen_index, seen_dir = set(), set(), set()
        for p in params:
            for k in ("name", "index", "oti_direction"):
                if k not in p:
                    raise PackageError("parameter %r missing %r" % (p, k))
            name, idx, dr = p["name"], int(p["index"]), int(p["oti_direction"])
            if not (1 <= idx <= int(m["nprops"])):
                raise PackageError("parameter %r index %d out of [1, nprops=%d]"
                                   % (name, idx, m["nprops"]))
            if dr < 1:
                raise PackageError("parameter %r oti_direction must be >= 1, got %d"
                                   % (name, dr))
            if name in seen_name:
                raise PackageError("duplicate parameter name %r" % name)
            if idx in seen_index:
                raise PackageError("two parameters share PROPS index %d" % idx)
            if dr in seen_dir:
                raise PackageError("two parameters share OTI direction %d" % dr)
            seen_name.add(name); seen_index.add(idx); seen_dir.add(dr)
        bins = m["binaries"]
        for side in ("regular", "oti"):
            if side not in bins:
                raise PackageError("binaries.%s is required (matched-twin rule)" % side)
            if "hash" not in bins[side]:
                raise PackageError("binaries.%s.hash is required" % side)
        if "path" not in bins["oti"]:
            raise PackageError("binaries.oti.path (the loadable .so) is required")

    # -- accessors -------------------------------------------------------- #
    @property
    def model_id(self) -> str:
        return str(self.raw["model_id"])

    @property
    def ntens(self) -> int:
        return int(self.raw["ntens"])

    @property
    def nprops(self) -> int:
        return int(self.raw["nprops"])

    @property
    def nstatev(self) -> int:
        return int(self.raw["nstatev"])

    @property
    def kinematics_code(self) -> int:
        return _KIN_CODE[self.raw["kinematics"]]

    @property
    def parameters(self) -> Dict[str, Dict[str, Any]]:
        return {p["name"]: dict(p) for p in self.raw["parameters"]}

    @property
    def param_names(self) -> List[str]:
        return [p["name"] for p in self.raw["parameters"]]

    def param_index(self, name: str) -> int:
        """1-based PROPS index of a named parameter (raises if unknown)."""
        for p in self.raw["parameters"]:
            if p["name"] == name:
                return int(p["index"])
        raise PackageError(
            "parameter %r is not in the material package (has: %s)"
            % (name, ", ".join(self.param_names)))

    def binary_path(self) -> str:
        p = self.raw["binaries"]["oti"]["path"]
        return p if os.path.isabs(p) else os.path.join(self.base_dir, p)

    def regular_hash(self) -> str:
        return str(self.raw["binaries"]["regular"]["hash"])

    # -- consistency checks ---------------------------------------------- #
    def check_against_abi(self, abi: MaterialABI) -> None:
        """The loaded binary's self-reported dimensions/id must match the
        manifest. A mismatch means the manifest describes a different binary."""
        d = abi.desc
        problems = []
        if abi.model_id != self.model_id:
            problems.append("model_id: binary=%r manifest=%r" % (abi.model_id, self.model_id))
        if d.ntens != self.ntens:
            problems.append("ntens: binary=%d manifest=%d" % (d.ntens, self.ntens))
        if d.nprops != self.nprops:
            problems.append("nprops: binary=%d manifest=%d" % (d.nprops, self.nprops))
        if d.nstatev != self.nstatev:
            problems.append("nstatev: binary=%d manifest=%d" % (d.nstatev, self.nstatev))
        if d.kinematics != self.kinematics_code:
            problems.append("kinematics: binary=%d manifest=%d"
                            % (d.kinematics, self.kinematics_code))
        if problems:
            raise PackageError(
                "loaded binary does not match its manifest: " + "; ".join(problems))

    def check_contract_version(self) -> None:
        """Refuse a manifest built against a different version of the shared
        ABI/manifest contract (drift between the provider and this loader)."""
        declared = self.raw.get("contract_version")
        expected = expected_contract_version()
        if declared and expected and declared != expected:
            raise PackageError(
                "manifest contract_version %r != this tool's contract %r; the "
                "shared ABI/manifest contract has drifted -- rebuild against a "
                "matching contract version" % (declared, expected))

    def check_binary_hash(self) -> None:
        """Confirm the loaded binary IS the one the manifest describes, by hashing
        the file and comparing to binaries.oti.hash, and confirm it is
        loadable on THIS platform before any linking/loading is attempted.
        Without these, a substituted, stale, or foreign-platform binary would
        fail late and opaquely."""
        declared = self.raw["binaries"]["oti"].get("hash")
        if not declared:
            raise PackageError("binaries.oti.hash is required to verify binary identity")
        path = self.binary_path()
        try:
            actual = hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
        except OSError as exc:
            raise PackageError("cannot hash OTI binary %r: %s" % (path, exc))
        if actual != declared:
            raise TwinMismatchError(
                "loaded OTI binary hash %r != manifest binaries.oti.hash %r -- the "
                "binary on disk is not the one this manifest describes" % (actual, declared))
        # platform compatibility: reject a foreign build with a clear diagnostic
        # (the manifest 'binary' block is optional; the format probe always runs).
        try:
            check_binary_compatibility(path, self.raw.get("binary"))
        except IncompatibleBinaryError as exc:
            raise PackageError(str(exc)) from exc

    def check_twin(self, record_provenance: Mapping[str, Any]) -> None:
        """Refuse to replay a record whose producing REGULAR model is not the
        matched twin of this OTI binary."""
        rec_model = record_provenance.get("model_id")
        rec_hash = record_provenance.get("regular_hash")
        if rec_model != self.model_id:
            raise TwinMismatchError(
                "record was produced by model_id=%r but the OTI package is %r; "
                "these are not a matched pair" % (rec_model, self.model_id))
        if rec_hash is None:
            raise TwinMismatchError(
                "record has no regular_hash provenance; cannot confirm the OTI "
                "binary is the matched twin of the model that produced it")
        if str(rec_hash) != self.regular_hash():
            raise TwinMismatchError(
                "record's regular binary hash %r != this package's regular twin "
                "hash %r; the OTI binary is not built from the same source as the "
                "model that ran the analysis" % (rec_hash, self.regular_hash()))
