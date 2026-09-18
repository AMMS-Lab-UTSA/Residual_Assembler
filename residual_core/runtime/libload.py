"""Cross-platform loading of the gfortran-built libraries Program 2 links.

Program 2 links the collaborator's distributed ``.obj`` into a small loadable
library and calls it through :mod:`ctypes`.  Windows and Linux differ in two
ways that matter:

* **File name.**  A shared library is ``.so`` on Linux and ``.dll`` on Windows.

* **Dependency resolution.**  A gfortran-built library imports ``libgfortran``.
  On Linux ``ld.so`` finds it through the usual search path.  On Windows,
  Python 3.8+ deliberately stopped honouring ``PATH`` for *dependent* DLLs, so a
  ``libgfortran-5.dll`` next to ``gfortran.exe`` -- on ``PATH``, findable by
  every other tool -- is invisible to ``ctypes.CDLL``, and the load fails with a
  ``FileNotFoundError`` naming *our* library rather than the missing dependency.
  The directory must be registered with :func:`os.add_dll_directory`, and the
  handle that returns must stay alive while the library is loaded.

Nothing here hard-codes a compiler location: the runtime directories are
discovered by asking the ``gfortran`` that is actually on ``PATH``.

This mirrors ``umat_oti.runtime.libload`` in Program 1.  The two programs ship
independently -- a collaborator has Program 2 and a compiled object, never
Program 1 -- so the module is duplicated rather than imported across the
boundary.  Both are small, and the shared C-ABI version check already guards the
interface that actually has to agree.
"""
from __future__ import annotations

import ctypes
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache

IS_WINDOWS = sys.platform == "win32"

#: Runtime libraries a gfortran-built object may import.  Used both to ask the
#: compiler where its runtime lives and to explain a failed load.
_RUNTIME_PROBES = (
    ("libgfortran-5.dll", "libgfortran.so.5"),
    ("libquadmath-0.dll", "libquadmath.so.0"),
    ("libgcc_s_seh-1.dll", "libgcc_s.so.1"),
    ("libwinpthread-1.dll", None),
)

#: Windows DLLs that always come from the OS; never reported as missing.
_SYSTEM_DLL_RE = re.compile(
    r"^(api-ms-win-|ext-ms-win-|kernel32|kernelbase|ntdll|msvcrt|ucrtbase|"
    r"user32|advapi32|ole32|oleaut32|shell32|rpcrt4|sechost|combase|bcrypt|"
    r"vcruntime|python\d)",
    re.IGNORECASE,
)


class LibraryLoadError(RuntimeError):
    """Raised when a freshly built validation library cannot be loaded.

    The message is the diagnostic: platform, compiler, library, dependency
    list and -- when it can be determined -- the dependency that is missing.
    """


def shared_library_suffix() -> str:
    """``.dll`` on Windows, ``.so`` elsewhere."""
    return ".dll" if IS_WINDOWS else ".so"


@dataclass(frozen=True)
class FortranToolchain:
    """The gfortran actually in use, as discovered at runtime."""

    executable: str | None
    version: str
    runtime_dirs: tuple[str, ...] = field(default_factory=tuple)

    def describe(self) -> str:
        if not self.executable:
            return "gfortran: NOT FOUND on PATH (set FC to point at one)"
        return "gfortran: %s\n  version: %s" % (self.executable, self.version)


def _run(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return ""
    return (out.stdout or "") + (out.stderr or "")


def _compiler_executable() -> str | None:
    return shutil.which(os.environ.get("FC") or "gfortran")


def _print_file_name(exe: str, name: str) -> str | None:
    """``gfortran -print-file-name=X`` -> absolute path, or None if unknown.

    gcc echoes the bare name back when it cannot place the file, so an answer
    only counts when it is absolute and exists.
    """
    out = _run([exe, "-print-file-name=" + name]).strip().splitlines()
    if not out:
        return None
    cand = out[0].strip()
    if cand and cand != name and os.path.isabs(cand) and os.path.exists(cand):
        return cand
    return None


def _search_dirs(exe: str) -> list[str]:
    """Directories from ``gfortran -print-search-dirs``."""
    dirs: list[str] = []
    for line in _run([exe, "-print-search-dirs"]).splitlines():
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        if key.strip() not in ("libraries", "programs"):
            continue
        value = value.strip().lstrip("=")
        for part in value.split(os.pathsep):
            part = os.path.normpath(part.strip())
            if part and os.path.isdir(part):
                dirs.append(part)
    return dirs


@lru_cache(maxsize=1)
def fortran_toolchain() -> FortranToolchain:
    """Discover the active gfortran and the directories holding its runtime.

    Discovery order, all dynamic:

    1. the directory of the ``gfortran`` resolved from ``FC`` or ``PATH``
       (where MinGW keeps ``libgfortran-5.dll``);
    2. whatever ``gfortran -print-file-name=<runtime>`` reports for each known
       runtime library;
    3. ``gfortran -print-search-dirs``.
    """
    exe = _compiler_executable()
    if not exe:
        return FortranToolchain(None, "unknown", ())

    version = ""
    for line in _run([exe, "--version"]).splitlines():
        if line.strip():
            version = line.strip()
            break

    dirs: list[str] = []

    def _add(path: str | None) -> None:
        if not path:
            return
        path = os.path.normpath(path)
        if os.path.isdir(path) and path not in dirs:
            dirs.append(path)

    _add(os.path.dirname(os.path.abspath(exe)))
    for win_name, posix_name in _RUNTIME_PROBES:
        probe = win_name if IS_WINDOWS else posix_name
        if probe:
            found = _print_file_name(exe, probe)
            if found:
                _add(os.path.dirname(found))
    for d in _search_dirs(exe):
        _add(d)

    return FortranToolchain(exe, version or "unknown", tuple(dirs))


def fortran_runtime_dirs() -> tuple[str, ...]:
    """Directories that must be searchable for the gfortran runtime."""
    return fortran_toolchain().runtime_dirs


def static_fortran_link_flags() -> list[str]:
    """Flags that fold the gfortran runtime into the library itself.

    Offered as a *fallback* (and as an explicit opt-in through
    ``UMAT_OTI_STATIC_FORTRAN=1``), not as the primary mechanism: the shipped
    product is a relocatable object that the collaborator links themselves, so
    the dynamic path is the one that has to work.
    """
    return ["-static-libgfortran", "-static-libgcc"]


def prefer_static_runtime() -> bool:
    return os.environ.get("UMAT_OTI_STATIC_FORTRAN", "").strip().lower() in ("1", "true", "yes", "on")


def library_dependencies(path: str) -> list[str]:
    """Names of the shared libraries *path* imports (best effort)."""
    objdump = shutil.which("objdump")
    names: list[str] = []
    if objdump:
        text = _run([objdump, "-p", path])
        pattern = r"DLL Name:\s*(\S+)" if IS_WINDOWS else r"NEEDED\s+(\S+)"
        names = re.findall(pattern, text)
    if not names and not IS_WINDOWS:
        readelf = shutil.which("readelf")
        if readelf:
            names = re.findall(r"NEEDED.*\[(.+?)\]", _run([readelf, "-d", path]))
    seen: list[str] = []
    for n in names:
        if n not in seen:
            seen.append(n)
    return seen


def _visible_dirs(lib_path: str, registered: list[str]) -> list[str]:
    """Directories ctypes can actually resolve dependencies from."""
    dirs = [os.path.dirname(os.path.abspath(lib_path))] + list(registered)
    if IS_WINDOWS:
        root = os.environ.get("SystemRoot", r"C:\Windows")
        dirs += [os.path.join(root, "System32"), root]
    else:
        dirs += [d for d in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep) if d]
        dirs += ["/lib", "/lib64", "/usr/lib", "/usr/lib64", "/usr/local/lib"]
        dirs += ["/usr/lib/x86_64-linux-gnu", "/lib/x86_64-linux-gnu"]
    return [d for d in dirs if d and os.path.isdir(d)]


def missing_dependencies(lib_path: str, registered: list[str]) -> list[str]:
    """Imported libraries that are not resolvable from the visible directories."""
    dirs = _visible_dirs(lib_path, registered)
    missing = []
    for dep in library_dependencies(lib_path):
        if IS_WINDOWS and _SYSTEM_DLL_RE.match(dep):
            continue
        if not IS_WINDOWS and dep.startswith(("libc.so", "libm.so", "libdl.so", "libpthread.so", "ld-linux")):
            continue
        if not any(os.path.exists(os.path.join(d, dep)) for d in dirs):
            missing.append(dep)
    return missing


@dataclass
class LoadedLibrary:
    """A loaded library plus the DLL-directory handles that keep it loadable.

    On Windows the cookies returned by :func:`os.add_dll_directory` must stay
    alive: closing them removes the directory from the dependent-DLL search
    path, which breaks any *later* lazy resolution.  Holding them on the object
    ties their lifetime to the library's.
    """

    lib: ctypes.CDLL
    path: str
    _cookies: list = field(default_factory=list, repr=False)

    def __getattr__(self, item):          # delegate lo.orig_path(...) etc.
        return getattr(self.lib, item)

    def close(self) -> None:
        for cookie in self._cookies:
            try:
                cookie.close()
            except Exception:
                pass
        self._cookies.clear()


#: Cookies are also parked here so a library outlives the LoadedLibrary handle
#: if a caller drops it -- the process keeps the module mapped either way.
_RETAINED: list = []


def _diagnostic(lib_path: str, exc: BaseException, registered: list[str]) -> str:
    tc = fortran_toolchain()
    deps = library_dependencies(lib_path)
    missing = missing_dependencies(lib_path, registered)
    size = os.path.getsize(lib_path) if os.path.exists(lib_path) else -1

    lines = [
        "could not load the generated validation library",
        "",
        "  platform          : %s (%s), Python %s" % (sys.platform, os.name, sys.version.split()[0]),
        "  %s" % tc.describe().replace("\n", "\n  "),
        "  library           : %s" % lib_path,
        "  library exists    : %s (%d bytes)" % (os.path.exists(lib_path), size),
        "  imports           : %s" % (", ".join(deps) if deps else "(objdump unavailable)"),
    ]
    if missing:
        lines += [
            "  MISSING           : %s" % ", ".join(missing),
            "",
            "  -> the library built fine; its dependency could not be found.",
        ]
    else:
        lines += ["", "  -> every dependency resolved; the loader still refused the library."]
    lines += [
        "  searched          :",
        *["      %s" % d for d in _visible_dirs(lib_path, registered)],
    ]
    if IS_WINDOWS:
        lines += [
            "",
            "  Note: Python 3.8+ ignores PATH when resolving dependent DLLs, so a",
            "  gfortran runtime on PATH is not enough -- its directory must be passed",
            "  to os.add_dll_directory().  That is what this loader does; if the",
            "  directory above is wrong, point FC at the gfortran you mean to use.",
        ]
    lines += [
        "",
        "  Retry with the runtime linked statically:  UMAT_OTI_STATIC_FORTRAN=1",
        "",
        "  underlying error   : %s: %s" % (type(exc).__name__, exc),
    ]
    return "\n".join(lines)


def load_shared_library(path: str, *, extra_dirs=()) -> LoadedLibrary:
    """Load *path* with the gfortran runtime made resolvable.

    Raises :class:`LibraryLoadError` -- carrying platform, compiler, library and
    missing-dependency detail -- instead of a bare ``OSError``.
    """
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise LibraryLoadError(
            "validation library was never produced: %s\n  %s" % (path, fortran_toolchain().describe())
        )

    registered: list[str] = []
    cookies: list = []
    if IS_WINDOWS:
        for d in list(extra_dirs) + list(fortran_runtime_dirs()):
            if not os.path.isdir(d) or d in registered:
                continue
            try:
                cookies.append(os.add_dll_directory(d))
            except OSError:
                continue
            registered.append(d)
    else:
        registered = [d for d in extra_dirs if os.path.isdir(d)]

    try:
        lib = ctypes.CDLL(path)
    except OSError as exc:
        for cookie in cookies:
            try:
                cookie.close()
            except Exception:
                pass
        raise LibraryLoadError(_diagnostic(path, exc, registered)) from exc

    _RETAINED.extend(cookies)
    return LoadedLibrary(lib, path, cookies)
