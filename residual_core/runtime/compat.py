"""Detect an incompatible material package BEFORE linking it.

A distributed material object is platform-specific: a Windows build is a COFF
``.obj`` (links into a ``.dll``); a Linux build is an ELF ``.o`` (links into a
``.so``). Handing the wrong one to the local toolchain fails at link time with an
opaque ``ld`` error. This module catches the mismatch first -- from the object's
own magic bytes and the ``binary`` metadata Program 1 records in the completed
contract -- and raises a clear, actionable diagnostic instead.

The check is authoritative on the *actual* object format (magic bytes) and
cross-checks the *declared* metadata (os / arch / pointer width) when present, so
it works even for an older package that predates the metadata block.
"""
from __future__ import annotations

import os
import platform
import struct
import sys


class IncompatibleBinaryError(RuntimeError):
    """The material object cannot be linked on this platform."""


def _current_os() -> str:
    if sys.platform == "win32":
        return "windows"
    if sys.platform == "darwin":
        return "macos"
    return "linux"


def object_kind(path: str) -> str:
    """Classify a compiled object/library by its leading bytes.

    Returns ``elf`` | ``pe`` | ``archive`` | ``coff`` | ``unknown``.  ELF, PE and
    Unix archives have unambiguous magic; a relocatable Windows COFF object has
    no universal magic, so "not ELF / not PE / not archive, but a real object" is
    reported as ``coff``.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        return "unknown"
    if head[:4] == b"\x7fELF":
        return "elf"
    if head[:2] == b"MZ":
        return "pe"
    if head[:8] == b"!<arch>\n":
        return "archive"
    if len(head) >= 2:
        # COFF relocatable object: first 2 bytes are the machine id. Common ids:
        # 0x8664 (amd64) -> b'\x64\x86', 0x014c (x86) -> b'\x4c\x01', 0xAA64 (arm64).
        machine = struct.unpack("<H", head[:2])[0]
        if machine in (0x8664, 0x014C, 0xAA64, 0x01C4, 0x0200):
            return "coff"
    return "unknown"


# A platform accepts a whole *family* of object formats, not one kind: a
# relocatable object to be linked (COFF/.obj on Windows, ELF/.o on Linux) and an
# already-linked loadable library (PE/.dll on Windows, ELF/.so on Linux) are both
# native. The cross-family cases -- ELF on Windows, COFF/PE on Linux -- are the
# real incompatibilities.
_OS_FORMAT_FAMILY = {
    "windows": {"coff", "pe"},
    "linux": {"elf"},
    "macos": {"macho"},
}


def _accepted_kinds() -> set:
    return _OS_FORMAT_FAMILY.get(_current_os(), {"elf"})


def _family_label() -> str:
    return "/".join(sorted(k.upper() for k in _accepted_kinds()))


def _current_pointer_bits() -> int:
    return struct.calcsize("P") * 8


def check_binary_compatibility(object_path: str, binary_meta=None) -> dict:
    """Verify *object_path* can be linked on this platform, or raise.

    ``binary_meta`` is the completed contract's ``binary`` block (or None for an
    older package). Returns a small dict describing what was checked, for logging.

    Raises :class:`IncompatibleBinaryError` with a diagnostic naming the produced
    platform, the current platform, and the format mismatch.
    """
    if not os.path.exists(object_path):
        raise IncompatibleBinaryError("material object does not exist: %s" % object_path)

    cur_os = _current_os()
    kind = object_kind(object_path)
    accepted = _accepted_kinds()

    # authoritative check: the object's actual format family vs this OS's family.
    # Both a relocatable object (COFF/ELF) and a linked library (PE/ELF-shared)
    # of the native family are fine; only cross-family is rejected. An 'unknown'
    # or 'archive' kind is not blocked -- let the toolchain try.
    linkable = kind in accepted or kind in ("unknown", "archive")

    produced_os = None
    problems = []
    if binary_meta:
        produced_os = binary_meta.get("os")
        if produced_os and produced_os != cur_os:
            problems.append("built for OS %r but running on %r" % (produced_os, cur_os))
        bits = binary_meta.get("pointer_bits")
        if bits and int(bits) != _current_pointer_bits():
            problems.append("built %d-bit but running %d-bit" % (int(bits), _current_pointer_bits()))
        arch = (binary_meta.get("arch") or "").lower()
        cur_arch = (platform.machine() or "").lower()
        # amd64 and x86_64 are the same target under different spellings
        norm = {"amd64": "x86_64", "x64": "x86_64", "arm64": "aarch64"}
        if arch and cur_arch and norm.get(arch, arch) != norm.get(cur_arch, cur_arch):
            problems.append("built for arch %r but running on %r" % (arch, cur_arch))

    if not linkable:
        problems.insert(0, "object format is %s but this platform (%s) uses %s"
                        % (kind.upper(), cur_os, _family_label()))

    if problems:
        produced = (" produced by: os=%s arch=%s compiler=%s %s, format=%s\n"
                    % (binary_meta.get("os"), binary_meta.get("arch"),
                       binary_meta.get("compiler"), binary_meta.get("compiler_version"),
                       binary_meta.get("binary_format")) if binary_meta else "")
        raise IncompatibleBinaryError(
            "material object is not compatible with this platform.\n"
            "  object            : %s\n"
            "  object format     : %s (from magic bytes)\n"
            "%s"
            "  this platform     : os=%s arch=%s pointer=%d-bit, uses %s\n"
            "  problem(s)        : %s\n"
            "  -> obtain the material object built by the %s toolchain, or rebuild it there.\n"
            "     The JSON contract is platform-independent; the compiled object is not."
            % (object_path, kind.upper(), produced, cur_os, (platform.machine() or "?").lower(),
               _current_pointer_bits(), _family_label(), "; ".join(problems),
               (produced_os or cur_os)))

    return {"object_kind": kind, "os": cur_os, "accepted_kinds": sorted(accepted),
            "declared_os": produced_os, "compatible": True}
