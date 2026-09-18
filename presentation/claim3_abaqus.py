"""Abaqus part of claim 3: the single-C3D8 m5_cpflow shear test solved by Abaqus.

The deck is the configuration the residual method is run on (same nodes, same
boundary conditions, same increments), with the repository's unchanged
``parameter_sensitivity/models/m5_cpflow/umat.for``. Stresses are read from the
ASCII results file (``*FILE FORMAT, ASCII``), which carries double precision;
the ODB (single precision) is not used for derivatives.

Jobs run one at a time, named ``claudeP_c3_*``. A job counts as completed only
when its .sta says "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" (Abaqus 2021.HF5
aborts in teardown after writing complete results, so the exit code is not a
verdict).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from presentation import c3d8_residual as C
from presentation.common import umat_repo

JOB_PREFIX = os.environ.get("PRESENTATION_ABAQUS_PREFIX", "claudeP_c3")
COMPLETED = "THE ANALYSIS HAS COMPLETED SUCCESSFULLY"


def deck(props: Sequence[float], n_increments: int, dt: float, gamma: float, mesh: int = 1) -> str:
    """Abaqus input for the shear test. mesh=1: the claim-3 single element; mesh>1:
    the homogeneous (every boundary node prescribed) larger-mesh test."""
    problem = C.simple_shear_problem(mesh, gamma, 1, 1.0, all_boundary=mesh > 1)
    lines = ["*HEADING", f"claudeP claim 3: m5_cpflow simple shear, {mesh}x{mesh}x{mesh} C3D8",
             "*NODE, NSET=ALLN"]
    for node, (x, y, z) in enumerate(problem.coords, 1):
        lines.append(f"{node}, {float(x)!r}, {float(y)!r}, {float(z)!r}")
    lines.append("*ELEMENT, TYPE=C3D8, ELSET=EALL")
    for element, conn in enumerate(problem.conn, 1):
        lines.append(f"{element}, " + ", ".join(str(n + 1) for n in conn))
    lines += ["*SOLID SECTION, ELSET=EALL, MATERIAL=CPFLOW", "*MATERIAL, NAME=CPFLOW",
              "*DEPVAR", "1", f"*USER MATERIAL, CONSTANTS={len(props)}"]
    values = [repr(float(p)) for p in props]
    for start in range(0, len(values), 8):
        lines.append(", ".join(values[start:start + 8]))
    lines += ["*STEP, INC=100000", "*STATIC, DIRECT", f"{dt!r}, {dt * n_increments!r}", "*BOUNDARY"]
    # every prescribed dof with its final value; *STATIC ramps it linearly (RAMP)
    for dof, value in zip(problem.prescribed, problem.values[0]):
        lines.append(f"{dof // 3 + 1}, {dof % 3 + 1}, {dof % 3 + 1}, {float(value)!r}")
    lines += ["*OUTPUT, FIELD, FREQUENCY=1", "*ELEMENT OUTPUT", "S, SDV", "*NODE OUTPUT", "U, RF",
              "*FILE FORMAT, ASCII", "*EL FILE, POSITION=INTEGRATION POINTS, FREQUENCY=1", "S, SDV",
              "*END STEP"]
    return "\n".join(lines) + "\n"


def run_job(job: str, inp_text: str, work: Path, umat: Path, timeout: int = 3600) -> Dict[str, object]:
    job_dir = work / job
    if job_dir.exists():
        shutil.rmtree(job_dir)
    job_dir.mkdir(parents=True)
    (job_dir / f"{job}.inp").write_text(inp_text)
    command = ["abaqus", f"job={job}", f"input={job}.inp", f"user={umat}", "double=both",
               "cpus=1", "interactive"]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=job_dir, capture_output=True, text=True, timeout=timeout)
    seconds = time.perf_counter() - started
    (job_dir / "abaqus_stdout.log").write_text(completed.stdout + completed.stderr)
    sta = job_dir / f"{job}.sta"
    ok = sta.is_file() and COMPLETED in sta.read_text(errors="replace")
    return {"job": job, "command": " ".join(command), "returncode": completed.returncode,
            "completed": ok, "seconds": round(seconds, 1), "directory": str(job_dir)}


# --------------------------------------------------------------------------- #
# ASCII results file (.fil)
# --------------------------------------------------------------------------- #
_TOKEN = re.compile(r"I(\d\d)|D(.{22})|A(.{8})", re.S)


def read_fil(path: Path) -> List[List]:
    """Records of an ASCII Abaqus results file as lists of int/float/str."""
    text = "".join(line.rstrip("\n").ljust(80)[:80] for line in path.read_text().splitlines())
    records = []
    for chunk in text.split("*")[1:]:
        values, pos = [], 0
        while pos < len(chunk):
            char = chunk[pos]
            if char == "I":
                width = int(chunk[pos + 1:pos + 3])
                values.append(int(chunk[pos + 3:pos + 3 + width]))
                pos += 3 + width
            elif char == "D":
                values.append(float(chunk[pos + 1:pos + 23].replace("D", "E")))
                pos += 23
            elif char == "A":
                values.append(chunk[pos + 1:pos + 9])
                pos += 9
            elif char == " ":
                pos += 1
            else:
                raise ValueError(f"unexpected character {char!r} in {path}")
        if values:
            records.append(values)
    return records


def stresses_from_fil(path: Path):
    """Per increment: {(element, ip): stress(6)} and {(element, ip): sdv list}."""
    increments, current, key = [], None, None
    for record in read_fil(path):
        rtype = record[1]
        if rtype == 2000:                       # increment start
            current = {"S": {}, "SDV": {}}
            increments.append(current)
        elif rtype == 1 and current is not None:   # element header: element, ip, ...
            key = (record[2], record[3])
        elif rtype == 11 and current is not None:  # stress components
            current["S"][key] = np.array(record[2:8], float)
        elif rtype == 5 and current is not None:   # SDV
            current["SDV"][key] = list(record[2:])
    return [inc for inc in increments if inc["S"]]


def mises_series(increments) -> np.ndarray:
    """Volume-average-free check: every IP must agree (homogeneous field); returns
    per increment the mean sigma_vM over IPs and the IP spread."""
    out = []
    for inc in increments:
        stress = np.array(list(inc["S"].values()))
        vm = C.mises_field(stress)
        out.append((vm.mean(), vm.max() - vm.min()))
    return np.array(out)


def run(*, props, names, pidx, n_increments, dt, gamma, work: Path, python_vm, oti_dvm, analytic_dvm,
        steps=(1e-4, 1e-5), mesh: int = 4) -> Dict[str, object]:
    if not shutil.which("abaqus"):
        return {"status": "not_run: abaqus not on PATH"}
    work.mkdir(parents=True, exist_ok=True)
    umat = umat_repo() / "parameter_sensitivity" / "models" / "m5_cpflow" / "umat.for"
    report: Dict[str, object] = {"umat": str(umat), "work": str(work), "jobs": [],
                                 "judged_by": f".sta contains '{COMPLETED}'"}

    def solve(tag, prop, m=1):
        job = f"{JOB_PREFIX}_{tag}"
        result = run_job(job, deck(prop, n_increments, dt, gamma, m), work, umat)
        report["jobs"].append(result)
        if not result["completed"]:
            return None
        return mises_series(stresses_from_fil(Path(result["directory"]) / f"{job}.fil"))

    nominal = solve("nominal", props)
    if nominal is None:
        report["status"] = "nominal job did not complete"
        return report
    vm = nominal[:, 0]
    report["nominal"] = {
        "increments": int(len(vm)), "sigma_vm_final": float(vm[-1]),
        "ip_spread_max": float(nominal[:, 1].max()),
        "vs_python_c3d8_max_rel": float(np.max(np.abs(vm - python_vm)) / np.max(python_vm)),
    }
    fd = {}
    for step in steps:
        for k, name in enumerate(names):
            base = props[pidx[k]]
            h = step * abs(base)
            vals = []
            for sign, label in ((+1, "p"), (-1, "m")):
                prop = list(props); prop[pidx[k]] = base + sign * h
                series = solve(f"{name}_{label}{step:g}".replace("-", "m").replace(".", "p"), prop)
                vals.append(None if series is None else series[:, 0])
            if any(v is None for v in vals):
                continue
            fd.setdefault(name, {})[step] = (vals[0] - vals[1]) / (2.0 * h)
    comparison = {}
    for k, name in enumerate(names):
        if name not in fd:
            continue
        comparison[name] = {}
        for step, values in fd[name].items():
            scale = np.max(np.abs(analytic_dvm[:, k]))
            comparison[name][f"{step:g}"] = {
                "abaqus_fd_vs_oti_nrmse": float(np.sqrt(np.mean((values - oti_dvm[:, k]) ** 2)) / scale),
                "abaqus_fd_vs_analytic_nrmse": float(np.sqrt(np.mean((values - analytic_dvm[:, k]) ** 2)) / scale),
            }
    report["fd_of_abaqus"] = comparison
    if comparison:
        best = {name: min(v.values(), key=lambda e: e["abaqus_fd_vs_oti_nrmse"]) for name, v in comparison.items()}
        report["fd_of_abaqus_best_step_vs_oti_worst"] = max(b["abaqus_fd_vs_oti_nrmse"] for b in best.values())
    if mesh > 1:
        big = solve(f"mesh{mesh}", props, mesh)
        if big is not None:
            report["larger_mesh"] = {
                "mesh": f"{mesh}x{mesh}x{mesh}", "ip_spread_max": float(big[:, 1].max()),
                "vs_single_element_max_rel": float(np.max(np.abs(big[:, 0] - vm)) / np.max(vm))}
    report["status"] = "completed" if all(j["completed"] for j in report["jobs"]) else "some jobs failed"
    return report
