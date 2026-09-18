"""Runtime helpers: cross-platform loading of the libraries Program 2 links."""

from residual_core.runtime.compat import (
    IncompatibleBinaryError,
    check_binary_compatibility,
    object_kind,
)
from residual_core.runtime.libload import (
    LibraryLoadError,
    LoadedLibrary,
    fortran_runtime_dirs,
    fortran_toolchain,
    library_dependencies,
    load_shared_library,
    missing_dependencies,
    prefer_static_runtime,
    shared_library_suffix,
    static_fortran_link_flags,
)

__all__ = [
    "IncompatibleBinaryError",
    "check_binary_compatibility",
    "object_kind",
    "LibraryLoadError",
    "LoadedLibrary",
    "fortran_runtime_dirs",
    "fortran_toolchain",
    "library_dependencies",
    "load_shared_library",
    "missing_dependencies",
    "prefer_static_runtime",
    "shared_library_suffix",
    "static_fortran_link_flags",
]
