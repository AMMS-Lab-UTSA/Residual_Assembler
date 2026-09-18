"""Shared plumbing for the presentation-claim reproductions.

Nothing in here computes a derivative. It locates the two repositories, records
provenance, builds compiled OTI providers with the Program-1 command, and
provides the centred finite-difference ladder used as the independent reference
by several claims. Every claim module owns its own mathematics and verdicts.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Sequence

import numpy as np

RA_ROOT = Path(__file__).resolve().parents[1]
PRESENTATION = RA_ROOT / "verification"
EXPECTED = PRESENTATION / "expected"

#: Relative steps of the centred-difference ladder, h = step * max(|p|, 1).
FD_LADDER = (1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7)


# --------------------------------------------------------------------------- #
# repositories and environment
# --------------------------------------------------------------------------- #
def umat_repo() -> Path:
    """The UMAT_source_transformation checkout that provides Program 1."""
    env = os.environ.get("UMAT_OTI_REPO")
    if env:
        return Path(env).resolve()
    import umat_oti  # noqa: WPS433 - resolved lazily on purpose
    return Path(umat_oti.__file__).resolve().parents[2]


def import_origins() -> Dict[str, str]:
    """Where the three packages were actually imported from (for provenance)."""
    import residual_core
    import umat_oti
    origins = {"residual_core": residual_core.__file__, "umat_oti": umat_oti.__file__}
    try:
        import pyoti
        origins["pyoti"] = pyoti.__file__
    except Exception as exc:  # noqa: BLE001 - reported, never hidden
        origins["pyoti"] = f"unavailable: {exc}"
    return origins


def _same_package(imported_init: str, source_package: Path) -> bool:
    """True when an imported package is a byte-identical copy of the checkout's.

    A package installed from this very checkout (``pip install <checkout>``)
    lives in site-packages, not under the checkout, yet runs the same code.
    Every Python file on either side must exist on the other with the same
    bytes; anything else is a different version and is refused.
    """
    installed = Path(imported_init).resolve().parent
    if not source_package.is_dir():
        return False
    ours = {p.relative_to(source_package): p for p in source_package.rglob("*.py")
            if "__pycache__" not in p.parts}
    theirs = {p.relative_to(installed): p for p in installed.rglob("*.py")
              if "__pycache__" not in p.parts}
    return bool(ours) and ours.keys() == theirs.keys() and all(
        ours[name].read_bytes() == theirs[name].read_bytes() for name in ours)


def require_worktree_imports() -> Dict[str, str]:
    """Refuse to run against packages that are not this checkout's code."""
    origins = import_origins()
    ra = str(RA_ROOT)
    um = str(umat_repo())
    if not (origins["residual_core"].startswith(ra)
            or _same_package(origins["residual_core"], RA_ROOT / "residual_core")):
        raise SystemExit(f"residual_core imported from {origins['residual_core']}, not {ra}; "
                         "set PYTHONPATH=<RA>:<UMAT>/src:<otilib build>")
    if not (origins["umat_oti"].startswith(um)
            or _same_package(origins["umat_oti"], Path(um) / "src" / "umat_oti")):
        raise SystemExit(f"umat_oti imported from {origins['umat_oti']}, not {um}, and it is not "
                         "the same code; set UMAT_OTI_REPO and PYTHONPATH consistently")
    return origins


def _first_line(command: Sequence[str]) -> str:
    try:
        out = subprocess.run(list(command), capture_output=True, text=True, timeout=60)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0] if text else ""
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def git_state(path: Path) -> Dict[str, object]:
    def git(*args):
        try:
            return subprocess.check_output(["git", "-C", str(path), *args], text=True,
                                           stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return ""
    status = git("status", "--short")
    return {"path": str(path), "commit": git("rev-parse", "HEAD"),
            "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty_files": len([line for line in status.splitlines() if line.strip()])}


def provenance() -> Dict[str, object]:
    cpu = ""
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "timestamp_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "cpu": cpu,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "gfortran": _first_line(["gfortran", "--version"]),
        "ifort": _first_line(["ifort", "--version"]) if shutil.which("ifort") else "not on PATH",
        "abaqus": shutil.which("abaqus") or "not on PATH",
        "repositories": {"residual_assembler": git_state(RA_ROOT),
                         "umat_source_transformation": git_state(umat_repo())},
        "imports": import_origins(),
        "platform": platform.platform(),
    }


# --------------------------------------------------------------------------- #
# output locations
# --------------------------------------------------------------------------- #
def default_out() -> Path:
    """Small machine-readable results (JSON/CSV). Git-ignored."""
    return Path(os.environ.get("PRESENTATION_OUT", PRESENTATION / "results"))


def default_work() -> Path:
    """Large intermediate files (builds, FD runs, Abaqus jobs). Never in the repo."""
    env = os.environ.get("PRESENTATION_WORK")
    if env:
        return Path(env)
    return Path(tempfile.gettempdir()) / f"presentation_work_{os.getuid()}"


def fresh_dir(path: Path) -> Path:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def write_json(path: Path, payload) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(payload), indent=2, sort_keys=False) + "\n",
                    encoding="utf-8")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# Program 1: compiled OTI provider
# --------------------------------------------------------------------------- #
def build_provider(contract_path: Path, out_dir: Path) -> Dict[str, object]:
    """Build one relocatable OTI provider with the Program-1 command.

    Runs ``python -m umat_oti.provider build`` (the module behind the
    ``umat-oti-provider`` console script) in a subprocess with the current
    interpreter and import path, so the build is exactly what a developer runs.
    """
    out_dir = fresh_dir(Path(out_dir))
    command = [sys.executable, "-m", "umat_oti.provider", "build", str(contract_path),
               "--out", str(out_dir)]
    started = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, env=os.environ.copy())
    seconds = time.perf_counter() - started
    record: Dict[str, object] = {"command": " ".join(command[1:]).replace(str(RA_ROOT), "<RA>"),
                                 "returncode": completed.returncode,
                                 "seconds": round(seconds, 3)}
    if completed.returncode != 0:
        record["error"] = (completed.stdout + completed.stderr).strip()[-2000:]
        return record
    built = json.loads(completed.stdout)
    contract = json.loads(Path(built["contract"]).read_text(encoding="utf-8"))
    record.update(object=built["object"], contract_path=built["contract"],
                  build_dir=built["build_dir"], contract=contract,
                  object_sha256=sha256(Path(built["object"])))
    return record


# --------------------------------------------------------------------------- #
# independent reference: centred finite differences with a plateau
# --------------------------------------------------------------------------- #
def fd_ladder(evaluate: Callable[[List[float]], np.ndarray], props: Sequence[float],
              index: int, steps: Iterable[float] = FD_LADDER) -> Dict[float, np.ndarray]:
    """Centred differences of ``evaluate`` in PROPS(index+1) for every step.

    ``evaluate`` must run the ORIGINAL (untransformed) model; nothing here is
    allowed to see the quantity under test.
    """
    base = float(props[index])
    scale = abs(base) if base != 0.0 else 1.0
    out: Dict[float, np.ndarray] = {}
    for step in steps:
        h = step * scale
        plus = list(props); plus[index] = base + h
        minus = list(props); minus[index] = base - h
        out[float(step)] = (np.asarray(evaluate(plus), float)
                            - np.asarray(evaluate(minus), float)) / (2.0 * h)
    return out


def plateau(ladder: Mapping[float, np.ndarray]) -> Dict[str, object]:
    """Pick the step where two consecutive centred differences agree best.

    The choice uses only the finite differences themselves (never the value
    being verified). Returns the selected estimate, its step, and the
    discrepancy to the neighbouring step as the reference's own uncertainty.
    """
    steps = sorted(ladder, reverse=True)
    scale = max(float(np.max(np.abs(ladder[s]))) for s in steps) or 1.0
    diffs = []
    for coarse, fine in zip(steps, steps[1:]):
        diffs.append((float(np.max(np.abs(ladder[coarse] - ladder[fine]))) / scale, coarse, fine))
    best = min(diffs)
    return {"step": best[1], "value": ladder[best[1]],
            "uncertainty_scaled": best[0], "neighbour_step": best[2],
            "pair_discrepancies_scaled": {f"{c:g}->{f:g}": d for d, c, f in diffs}}


def format_e(value, digits: int = 1) -> str:
    if value is None:
        return "n/a"
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value == 0.0:
        return "0"
    return f"{value:.{digits}e}"
