"""Unit tests for the cross-platform library loader.

These do not require a build: they check the platform-aware behaviour and that a
failed load produces the structured diagnostic (platform, compiler, library,
missing dependency) rather than a bare OSError.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from umat_oti.runtime import (
    LibraryLoadError,
    fortran_toolchain,
    load_shared_library,
    shared_library_suffix,
)


def test_suffix_matches_platform():
    assert shared_library_suffix() == (".dll" if sys.platform == "win32" else ".so")


def test_toolchain_discovers_gfortran_when_present():
    tc = fortran_toolchain()
    # gfortran may or may not be installed in a given CI job; when it is, its
    # runtime directory must be discovered (that is the whole point on Windows).
    if tc.executable:
        assert tc.runtime_dirs, "gfortran found but no runtime directory discovered"
        assert os.path.isfile(tc.executable) or os.path.exists(tc.executable)


def test_missing_library_raises_structured_error(tmp_path):
    with pytest.raises(LibraryLoadError) as exc:
        load_shared_library(str(tmp_path / ("does_not_exist" + shared_library_suffix())))
    assert "never produced" in str(exc.value) or "does_not_exist" in str(exc.value)


def test_diagnostic_mentions_platform_and_library(tmp_path):
    # a file that exists but is not a valid library -> load fails with the
    # structured diagnostic, not a bare OSError
    bogus = tmp_path / ("bogus" + shared_library_suffix())
    bogus.write_bytes(b"not a real shared object")
    with pytest.raises(LibraryLoadError) as exc:
        load_shared_library(str(bogus))
    msg = str(exc.value)
    assert sys.platform in msg
    assert str(bogus) in msg


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
