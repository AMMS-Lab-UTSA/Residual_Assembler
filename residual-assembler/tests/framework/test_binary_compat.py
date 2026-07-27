"""Program 2 rejects a platform-incompatible material object BEFORE linking."""
import os
import struct
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.runtime import (
    IncompatibleBinaryError,
    check_binary_compatibility,
    object_kind,
)

IS_WINDOWS = sys.platform == "win32"


def _write(path, data):
    with open(path, "wb") as fh:
        fh.write(data)
    return str(path)


def test_object_kind_detects_elf_and_coff(tmp_path):
    elf = _write(tmp_path / "a.o", b"\x7fELF" + b"\x00" * 60)
    coff = _write(tmp_path / "b.obj", struct.pack("<H", 0x8664) + b"\x00" * 60)
    pe = _write(tmp_path / "c.dll", b"MZ" + b"\x00" * 60)
    assert object_kind(elf) == "elf"
    assert object_kind(coff) == "coff"
    assert object_kind(pe) == "pe"


def test_rejects_wrong_format_for_this_platform(tmp_path):
    # build the object of the OTHER platform's format and confirm rejection
    if IS_WINDOWS:
        wrong = _write(tmp_path / "linux.o", b"\x7fELF" + b"\x00" * 60)   # ELF on Windows
        meta = {"os": "linux", "arch": "x86_64", "pointer_bits": 64, "binary_format": "elf"}
    else:
        wrong = _write(tmp_path / "win.obj", struct.pack("<H", 0x8664) + b"\x00" * 60)  # COFF on Linux
        meta = {"os": "windows", "arch": "amd64", "pointer_bits": 64, "binary_format": "coff/pe"}
    with pytest.raises(IncompatibleBinaryError) as exc:
        check_binary_compatibility(wrong, meta)
    msg = str(exc.value)
    assert "not compatible" in msg
    assert ("ELF" in msg or "COFF" in msg)
    assert "contract is platform-independent" in msg


def test_accepts_native_relocatable_object(tmp_path):
    # the .obj hand-off: a relocatable object to be linked (COFF / ELF-object)
    if IS_WINDOWS:
        native = _write(tmp_path / "win.obj", struct.pack("<H", 0x8664) + b"\x00" * 60)
        meta = {"os": "windows", "arch": "amd64", "pointer_bits": 64}
    else:
        native = _write(tmp_path / "linux.o", b"\x7fELF" + b"\x00" * 60)
        meta = {"os": "linux", "arch": "x86_64", "pointer_bits": 64}
    assert check_binary_compatibility(native, meta)["compatible"] is True


def test_accepts_native_linked_library(tmp_path):
    # the material.resmat package: an already-linked loadable library
    # (PE/.dll on Windows, ELF-shared/.so on Linux). Both belong to the native
    # format family and must be accepted, not just relocatable objects.
    if IS_WINDOWS:
        native = _write(tmp_path / "libmat.so", b"MZ" + b"\x00" * 60)     # PE named .so
        meta = {"os": "windows", "arch": "amd64", "pointer_bits": 64}
    else:
        native = _write(tmp_path / "libmat.so", b"\x7fELF" + b"\x00" * 60)
        meta = {"os": "linux", "arch": "x86_64", "pointer_bits": 64}
    assert check_binary_compatibility(native, meta)["compatible"] is True


def test_pointer_width_mismatch_rejected(tmp_path):
    if IS_WINDOWS:
        native = _write(tmp_path / "win.obj", struct.pack("<H", 0x8664) + b"\x00" * 60)
        meta = {"os": "windows", "arch": "amd64", "pointer_bits": 32}
    else:
        native = _write(tmp_path / "linux.o", b"\x7fELF" + b"\x00" * 60)
        meta = {"os": "linux", "arch": "x86_64", "pointer_bits": 32}
    with pytest.raises(IncompatibleBinaryError) as exc:
        check_binary_compatibility(native, meta)
    assert "32-bit" in str(exc.value)


def test_missing_metadata_still_probes_format(tmp_path):
    # no metadata block (older package) -> still catch a wrong-format object
    if IS_WINDOWS:
        wrong = _write(tmp_path / "linux.o", b"\x7fELF" + b"\x00" * 60)
    else:
        wrong = _write(tmp_path / "win.obj", struct.pack("<H", 0x8664) + b"\x00" * 60)
    with pytest.raises(IncompatibleBinaryError):
        check_binary_compatibility(wrong, None)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
