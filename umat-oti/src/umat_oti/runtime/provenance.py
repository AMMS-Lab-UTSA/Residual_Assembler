"""Binary metadata and build provenance for the distributable material object.

A compiled material object is **platform-specific**: an object built by the
Windows toolchain is a PE/COFF ``.obj`` that links into a ``.dll``; one built by
the Linux toolchain is an ELF ``.o`` that links into a ``.so``. One binary
cannot serve both. The *contract* (JSON) is platform-independent; the *binary*
is not.

This module records, next to the completed contract, exactly what a binary is
and how it was built, so the collaborator's tool (Program 2) can reject an
incompatible package **before** it tries to link it -- with a clear diagnostic
rather than a linker error.
"""
from __future__ import annotations

import os
import platform
import re
import struct
import subprocess
import sys
import time

from umat_oti.runtime.libload import fortran_toolchain, shared_library_suffix

#: bump when the metadata *shape* changes (not on every build)
BINARY_METADATA_SCHEMA = "resasm_binary_metadata_v1"


def _elf_or_coff(path: str) -> str:
    """Classify a compiled object/library by its magic bytes."""
    try:
        with open(path, "rb") as fh:
            magic = fh.read(4)
    except OSError:
        return "unknown"
    if magic[:4] == b"\x7fELF":
        return "ELF"
    if magic[:2] == b"MZ":
        return "PE"                       # a linked .dll / .exe
    # COFF object (relocatable .obj / .o): 16-bit machine id, no universal magic.
    # gnu-ld -r output on Windows is COFF; on Linux it is ELF (caught above).
    return "COFF"


def _normalized_os() -> str:
    s = sys.platform
    if s == "win32":
        return "windows"
    if s == "darwin":
        return "macos"
    return "linux"


def _binary_format_for_os() -> str:
    """The object/loadable format this platform's toolchain produces."""
    return "coff/pe" if _normalized_os() == "windows" else "elf"


def _compiler_id() -> tuple[str, str]:
    tc = fortran_toolchain()
    version = tc.version or "unknown"
    m = re.search(r"(\d+\.\d+\.\d+)", version)
    return ("gfortran", m.group(1) if m else version)


def build_id(*hashes: str) -> str:
    """A short, deterministic build identifier from the inputs' hashes.

    Deterministic (no timestamp / randomness) so an unchanged source + contract
    on the same platform reproduces the same id.
    """
    import hashlib
    h = hashlib.sha256()
    h.update(("%s|%s|%s|" % (_normalized_os(), platform.machine().lower(), "|".join(hashes))).encode())
    return h.hexdigest()[:16]


def binary_metadata(object_path: str, *, abi_version, source_hash: str,
                    transform_hash: str = "", extra=None) -> dict:
    """The full binary-metadata block for a completed contract / manifest.

    Records OS, architecture, compiler + version, binary format (probed from the
    file's magic bytes), ABI/interface version, a deterministic build id, and the
    source/transform hashes that produced it.
    """
    compiler, compiler_version = _compiler_id()
    fmt_probe = _elf_or_coff(object_path) if os.path.exists(object_path) else "unknown"
    meta = {
        "schema": BINARY_METADATA_SCHEMA,
        "os": _normalized_os(),
        "arch": platform.machine().lower() or "unknown",
        "pointer_bits": struct.calcsize("P") * 8,
        "compiler": compiler,
        "compiler_version": compiler_version,
        "binary_format": _binary_format_for_os(),        # expected format for this OS
        "object_format_probed": fmt_probe,               # what the file actually is
        "loadable_suffix": shared_library_suffix(),      # .dll / .so
        "abi_version": abi_version,
        "build_id": build_id(source_hash, transform_hash),
        "source_hash": source_hash,
        "transform_hash": transform_hash,
    }
    if extra:
        meta.update(extra)
    return meta


def build_environment() -> dict:
    """Platform / interpreter / compiler snapshot for a provenance block."""
    compiler, compiler_version = _compiler_id()
    tc = fortran_toolchain()
    return {
        "os": _normalized_os(),
        "platform": platform.platform(),
        "arch": platform.machine().lower() or "unknown",
        "python_version": sys.version.split()[0],
        "numpy_version": _numpy_version(),
        "compiler": compiler,
        "compiler_version": compiler_version,
        "compiler_path": tc.executable or "",
    }


def _numpy_version() -> str:
    try:
        import numpy
        return numpy.__version__
    except Exception:                     # pragma: no cover
        return "unknown"


def git_commit(repo_dir: str) -> str:
    """The HEAD commit SHA of *repo_dir*, or a clear sentinel if unavailable."""
    try:
        out = subprocess.run(["git", "-C", repo_dir, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=15)
        sha = (out.stdout or "").strip()
        if out.returncode == 0 and sha:
            dirty = subprocess.run(["git", "-C", repo_dir, "status", "--porcelain"],
                                   capture_output=True, text=True, timeout=15)
            return sha + ("-dirty" if (dirty.stdout or "").strip() else "")
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown (no committed git state)"


def timestamp_utc(epoch: float | None = None) -> str:
    """ISO-8601 UTC timestamp.

    ``epoch`` may be supplied by the caller (workflow scripts forbid an argless
    clock); when omitted this falls back to ``time.time()`` for interactive CLI
    use, which is not on the resumable path.
    """
    t = time.time() if epoch is None else float(epoch)
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))


def provenance_block(*, repos: dict, props_values, ladder, step_selection: str,
                     tolerances: dict, metric_defs: dict, status: str, reason: str = "",
                     epoch: float | None = None, extra=None) -> dict:
    """A complete provenance record to embed in a validation JSON.

    ``repos`` maps a label -> repo directory; each resolves to a commit SHA.
    """
    block = {
        "schema": "resasm_validation_provenance_v1",
        "git_commits": {name: git_commit(path) for name, path in repos.items()},
        "environment": build_environment(),
        "props_values": list(props_values),
        "fd_ladder": list(ladder),
        "step_selection": step_selection,
        "comparison_tolerances": dict(tolerances),
        "metric_definitions": dict(metric_defs),
        "timestamp_utc": timestamp_utc(epoch),
        "status": status,
        "reason": reason,
    }
    if extra:
        block.update(extra)
    return block
