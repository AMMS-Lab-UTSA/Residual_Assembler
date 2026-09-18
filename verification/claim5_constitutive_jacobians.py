#!/usr/bin/env python3
"""Slide 25: the internal constitutive Jacobians of the ICP UMAT family.

What the slide says
    For HIN, NKH_1.02, PCL, PCLI, PCLI_R, PCLK, PCO, VPDCL, VPDCL_R and VPDCO, the
    hand-coded constitutive Jacobians FJAC, DETDG, GDIA, ANP1P, BNP1P, CEVPI were
    "compared ... from the original results and the transformed ones": 14
    entries "Exact" (zero measured difference) and 5 "Pass" (within 1e-5), 19 in
    all (manuscript claim TABLE-3).

What existed
    No executed record of that 19-entry table exists in either repository (the
    compliance matrix calls the old Table-3 artefact "header-only").
    ``tools/run_internal_jacobian_round.py`` verifies only the Newton-update
    Jacobian FJAC, and only where a property vector is declared.

What this script measures
    For every (UMAT, symbol) pair that occurs in the source, at a converged state
    of the model's own local Newton solve (seeded at the iterate it converged to,
    one increment only, as the repository's internal-Jacobian probe does):

      hand-coded  the value the ORIGINAL source computes for the symbol,
      OTI         d(response)/d(iterate) extracted by the OTI build of the same
                  injected source (``transform_umat_for_parameter_sensitivity``),
      FD          centred differences of the ORIGINAL build (step ladder,
                  plateau-selected), the independent reference.

    The derivative each symbol stands for is declared, as a derivative request
    (seed, response, target), after reading the sources:
      FJAC  = dFGAM/dGAM_PAR   (discovered: the Newton update GAM_PAR -= FGAM/FJAC)
      DETDG = dPHIINV/dGAM_PAR (Perzyna overstress inverse)
      ANP1P = dANP1/dGAM_PAR,  BNP1P = dBNP1/dGAM_PAR  (NKH hardening moduli)
      GDIA  = dDIAG/dGAM_PAR   (diagonal, component by component)
    CEVPI (HIN) is not a derivative with respect to a local iterate: it is the
    inverse of HIN's local viscoplastic operator and the UMAT returns
    DDSDDE = (1-D) CEVPI, so it is reported from the slide-8 DDSDDE comparison
    (claim 4) and not re-measured here.

    Classification as on the slide: "Exact" = OTI and hand-coded identical;
    "Pass" = relative difference <= 1e-5; otherwise the measured difference is
    reported. OTI and hand-coded are both also scored against FD.

    Material: the probe vector the paired Abaqus validation uses for these
    sources (unit constants, 0.3 for Poisson-like names; not a physical
    material), with the uniaxial strain path of the repository's VPDCL entry.

    python verification/claim5_constitutive_jacobians.py
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from verification.common import (  # noqa: E402
    default_out, default_work, format_e, fresh_dir, provenance, require_worktree_imports, umat_repo,
    write_json,
)

UMATS = ("UMAT_HIN", "UMAT_NKH_1.02", "UMAT_PCL", "UMAT_PCLI", "UMAT_PCLI_R", "UMAT_PCLK",
         "UMAT_PCO", "UMAT_VPDCL", "UMAT_VPDCL_R", "UMAT_VPDCO")
SYMBOLS = ("FJAC", "DETDG", "GDIA", "ANP1P", "BNP1P", "CEVPI")
#: response each hand-coded symbol is the derivative of, w.r.t. the local iterate
RESPONSE = {"DETDG": "PHIINV", "ANP1P": "ANP1", "BNP1P": "BNP1", "GDIA": "DIAG"}
PASS_TOLERANCE = 1.0e-5
STATE_HEADROOM = 32
#: The ICP sources read their hardening table through KUHARD(..., NVALUE, PROPS(3)),
#: i.e. PROPS(3..2+2*NVALUE), past the largest literal PROPS index the dimension
#: inference sees; the probe vector is extended with the same unit constants so
#: those reads stay in bounds (Abaqus would hand the routine whatever follows).
PROPS_HEADROOM = 8
#: Declared probe loading, tried in order until the local solve is entered. Path A
#: is the uniaxial strain increment the repository declares for UMAT_VPDCL
#: (parameter_sensitivity/internal_jacobian_sources.json). With the unit probe
#: constants it can stay elastic, so path B drives the same component to twice
#: the unit yield strain. Path C is pure shear, for sources whose probe Poisson
#: ratio is capped near 0.5 (PCLI, PCLI_R read ENU=MIN(PROPS(2),ENUMAX), which the
#: probe's Poisson-name rule does not recognise, so PROPS(2)=1.0 is capped and
#: uniaxial strain produces an almost purely hydrostatic stress).
PROBE_PATHS = (("A: 4 x 1e-3 uniaxial strain", 0, 1.0e-3), ("B: 4 x 0.5 uniaxial strain", 0, 0.5),
               ("C: 4 x 1.0 engineering shear strain (component 4)", 3, 1.0))
#: NTENS override with its reason. The benchmark contract declares 6, but the
#: source reads its back stress with DO K1=10,2*NTENS+5 into XBACK(NTENS): in
#: bounds only for NTENS=4 (a -fcheck=bounds build stops at XBACK(7) with NTENS=6,
#: and both archived Abaqus runs of it at NTENS=6 failed).
#: PROPS overrides with their reason (1-based index -> value).
PROPS_OVERRIDE = {
    # IF(PROPS(1).NE.0) THTA=PROPS(1) ELSE THTA=TEMP; DTHTA=DTEMP: with the probe's
    # PROPS(1)=1.0, DTHTA is used at line 111 (thermal strain) without ever being
    # set (the snan gate traps there), so the temperature is taken from TEMP instead.
    "UMAT_NKH_1.02": ({1: 0.0}, "PROPS(1)=0: otherwise DTHTA is read uninitialised (line 111)"),
}
NTENS_OVERRIDE = {
    "UMAT_VPDCL_R": (4, "back-stress read DO K1=10,2*NTENS+5 is in bounds only for NTENS=4"),
    # PCLI/PCLI_R at the contract's NTENS=6: the uniaxial probe returns NaN in
    # STATEV(14) and the shear probe never yields; both sources assign the
    # shear stiffness for component 4 only (DO K1=4,4), i.e. plane strain.
    "UMAT_PCLI": (4, "at NTENS=6 the probe returns NaN state and never yields; shear stiffness set for component 4 only"),
    "UMAT_PCLI_R": (4, "at NTENS=6 the probe returns NaN state and never yields; shear stiffness set for component 4 only"),
}
FD_STEPS = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
SLIDE = {
    "UMAT_HIN": {"CEVPI": "Exact"},
    "UMAT_NKH_1.02": {"FJAC": "Pass", "DETDG": "Exact", "ANP1P": "Exact", "BNP1P": "Exact"},
    "UMAT_PCL": {"FJAC": "Exact"}, "UMAT_PCLI": {"FJAC": "Exact", "GDIA": "Exact"},
    "UMAT_PCLI_R": {"FJAC": "Exact"}, "UMAT_PCLK": {"FJAC": "Exact"},
    "UMAT_PCO": {"FJAC": "Exact", "GDIA": "Exact"}, "UMAT_VPDCL": {"FJAC": "Pass", "DETDG": "Exact"},
    "UMAT_VPDCL_R": {"FJAC": "Pass", "DETDG": "Exact"},
    "UMAT_VPDCO": {"FJAC": "Pass", "DETDG": "Exact", "GDIA": "Pass"},
}


def _code_lines(text: str) -> List[str]:
    return [line for line in text.splitlines() if line[:1] not in ("C", "c", "*", "!")]


def symbols_present(text: str) -> List[str]:
    body = "\n".join(_code_lines(text)).upper()
    return [s for s in SYMBOLS if re.search(rf"\b{s}\b", body)]


def probe_material(source: Path, ntens_hint: int):
    """The paired-validation probe: dimensions inferred from the source, unit
    constants with 0.3 for Poisson-like names (umat_oti.validation.job_builder)."""
    from umat_oti.validation import job_builder as jb
    text = source.read_text(errors="replace")
    ntens, _ = jb.infer_validation_ntens_from_source(text, fallback_ntens=ntens_hint)
    nstatv, nprops = jb.infer_validation_dimensions_from_source(text, ntens=ntens)
    block = jb._props_lines(nprops + PROPS_HEADROOM, text)
    props = [float(v) for v in block.replace("\n", "").split(",") if v.strip()]
    return ntens, nstatv, props


def resolved_source(source: Path, work: Path) -> Tuple[Path, Dict[str, object]]:
    """The source, or its resolved routine closure when it calls helpers it does not define."""
    from umat_oti.transform.dependency_resolution import (
        DependencyResolutionError, combined_source, resolve_closure,
    )
    try:
        graph = resolve_closure(source, entry="UMAT", roots=[source.parent])
    except DependencyResolutionError as exc:
        return source, {"closure": "unresolvable", "detail": str(exc)[:300]}
    if graph.missing:
        return source, {"closure": "missing", "missing": [m.symbol for m in graph.missing]}
    if not graph.is_multi_file:
        return source, {"closure": "single_file"}
    target = work / f"{source.stem}_resolved.for"
    target.write_text(combined_source(graph), encoding="utf-8")
    return target, {"closure": "multi_file", "resolved": sorted(graph.resolved)}


def uninitialised_read_gate(source: Path, props: List[float], path: List[List[float]], ntens: int,
                            nstatv: int, work: Path) -> Dict[str, object]:
    """Run the ORIGINAL once with every local REAL initialised to a signalling NaN
    and invalid operations trapped. A trap means the source computes with memory
    it never wrote (or out of bounds), so its answer depends on the build."""
    from umat_oti.validation.parameter_sensitivity_validation import ABA_PARAM, driver_source
    work.mkdir(parents=True, exist_ok=True)
    for name in ("aba_param.inc", "ABA_PARAM.INC", "ABA_PARAM.inc", "aba_param.INC"):
        (work / name).write_text(ABA_PARAM)
    (work / "gate_driver.f90").write_text(driver_source(ntens=ntens, nstatv=nstatv, nprops=len(props)))
    flags = ["-O0", "-g", "-finit-real=snan", "-ffpe-trap=invalid,zero", "-fcheck=bounds"]
    compile_umat = subprocess.run(["gfortran", *flags, "-ffixed-form", "-ffixed-line-length-none",
                                   "-std=legacy", "-I", str(work), "-c", str(source), "-o", "gate_umat.o"],
                                  cwd=work, capture_output=True, text=True)
    link = subprocess.run(["gfortran", *flags, "gate_driver.f90", "gate_umat.o", "-o", "gate"],
                          cwd=work, capture_output=True, text=True) if compile_umat.returncode == 0 else compile_umat
    if link.returncode != 0:
        return {"status": "gate_not_built", "detail": (link.stderr or "")[-300:]}
    stdin = " ".join(repr(float(v)) for v in props) + f"\n{len(path)}\n" + "\n".join(
        " ".join(repr(float(v)) for v in row) for row in path) + "\n"
    run = subprocess.run([str(work / "gate")], input=stdin, capture_output=True, text=True, cwd=work)
    if run.returncode == 0:
        return {"status": "clean", "flags": " ".join(flags)}
    where = re.findall(r"at (\S+\.for):(\d+)", run.stderr) or re.findall(r"At line (\d+) of file (\S+)", run.stderr)
    return {"status": "trapped", "flags": " ".join(flags), "where": where[:3],
            "message": (run.stderr.strip().splitlines() or [""])[0][:200]}


def _insert_records(text: str, solve, slots, pairs: List[Tuple[str, str]], first_slot: int) -> str:
    """Add STATEV records of every (response, symbol) pair beside the probe's own."""
    lines = text.splitlines()
    anchor = f"STATEV({slots.residual})={solve.residual}"
    out, inserted = [], 0
    for line in lines:
        out.append(line)
        if line.strip().upper() == anchor.upper():
            indent = line[:len(line) - len(line.lstrip())]
            slot = first_slot
            for response, symbol in pairs:
                out.append(f"{indent}STATEV({slot})={response}")
                out.append(f"{indent}STATEV({slot + 1})={symbol}")
                slot += 2
            inserted += 1
    if not inserted:
        raise RuntimeError("probe record anchor not found in the injected source")
    return "\n".join(out) + "\n"


def measure(model: str, source: Path, symbols: List[str], props: List[float], ntens: int,
            nstatv: int, path: List[List[float]], work: Path) -> Dict[str, object]:
    from umat_oti.transform.internal_jacobian import discover_local_solves
    from umat_oti.transform.local_jacobian_probe import (
        PROBE_SENTINEL, inject_local_solve_probe, plan_probe_slots,
    )
    from umat_oti.transform.parameter_sensitivity_transform import (
        GenericPSContract, transform_umat_for_parameter_sensitivity,
    )
    from umat_oti.validation.parameter_sensitivity_validation import (
        build_original_driver, centered_fd, read_oti_csv, replay,
    )

    text = source.read_text(errors="replace")
    record: Dict[str, object] = {"model": model, "source": str(source)}
    solves = discover_local_solves(text)
    if not solves:
        record["status"] = "no_local_solve"
        return record
    solve = solves[0]
    record["solve"] = solve.as_dict()
    pairs: List[Tuple[str, str, str, int]] = []   # (symbol, response expr, symbol expr, component)
    for symbol in symbols:
        if symbol == "FJAC":
            pairs.append(("FJAC", solve.residual, solve.jacobian, 0))
        elif symbol == "GDIA":
            for i in range(1, ntens + 1):
                pairs.append(("GDIA", f"DIAG({i},{i})", f"GDIA({i},{i})", i))
        elif symbol in RESPONSE:
            pairs.append((symbol, RESPONSE[symbol], symbol, 0))
    nprops = len(props)
    ndi, nshr = 3, ntens - 3
    # headroom: a source may address state slots past the count inference finds
    # (UMAT_VPDCL_R corrupts the heap otherwise); unused slots stay zero
    nstatv = nstatv + STATE_HEADROOM
    work = fresh_dir(work)
    gate = uninitialised_read_gate(source, props, path, ntens, nstatv, work / "gate")
    record["uninitialised_read_gate"] = gate
    if gate["status"] != "clean":
        record["status"] = "original_reads_unwritten_memory" if gate["status"] == "trapped" else "gate_not_built"
        return record
    plain_exe = build_original_driver(source, work / "plain", ntens=ntens, nstatv=nstatv, nprops=nprops)
    plain = replay(plain_exe, props, path, ntens=ntens, nstatv=nstatv)

    chosen = None
    for offset in (0, 8, 32):
        base = plan_probe_slots(nstatv=nstatv, nprops=nprops, offset=offset)
        slots = dataclasses.replace(base, nstatv=base.nstatv + 2 * len(pairs))
        observe = inject_local_solve_probe(text, solve, slots, target_increment=1, override_iterate=False)
        observed_text = _insert_records(observe.source, solve, slots,
                                        [(r, s) for _, r, s, _ in pairs], base.nstatv + 1)
        observed_src = work / f"observe{offset}.for"
        observed_src.write_text(observed_text)
        exe = build_original_driver(observed_src, work / f"observe{offset}", ntens=ntens,
                                    nstatv=slots.nstatv, nprops=slots.nprops)
        run = replay(exe, list(props) + [0.0], path, ntens=ntens, nstatv=slots.nstatv)
        intact = sum(1 for row in run.statev if row[slots.sentinel - 1] == PROBE_SENTINEL)
        counted = sum(1 for row in run.statev if row[slots.counter - 1] >= 1.0)
        if intact and counted:
            chosen = (base, slots, run)
            break
    if chosen is None:
        record["status"] = "probe_not_placed"
        return record
    base, slots, observed = chosen
    drift = max((abs(a - b) for ra, rb in zip(plain.stress, observed.stress) for a, b in zip(ra, rb)), default=0.0)
    record["recording_stress_drift"] = drift
    if drift != 0.0:
        record["status"] = "recording_perturbs_primal"
        return record
    iterates = [row[slots.iterate - 1] for row in observed.statev]
    evaluations = [row[slots.counter - 1] for row in observed.statev]
    entered = [(abs(iterates[i]), i + 1) for i in range(len(iterates)) if evaluations[i] >= 2.0 and iterates[i] != 0.0]
    if not entered:
        record["status"] = "local_solve_not_entered_on_this_path"
        record["residual_evaluations"] = evaluations
        return record
    target = max(entered)[1]
    gamma = iterates[target - 1]
    record.update(target_increment=target, converged_iterate=gamma)

    seeded = inject_local_solve_probe(text, solve, slots, target_increment=target, override_iterate=True)
    seeded_text = _insert_records(seeded.source, solve, slots, [(r, s) for _, r, s, _ in pairs], base.nstatv + 1)
    seeded_src = work / "seeded.for"
    seeded_src.write_text(seeded_text)
    probe_props = list(props) + [gamma]
    seeded_exe = build_original_driver(seeded_src, work / "seeded", ntens=ntens, nstatv=slots.nstatv,
                                       nprops=slots.nprops)
    seeded_run = replay(seeded_exe, probe_props, path, ntens=ntens, nstatv=slots.nstatv)
    took = (seeded_run.statev[target - 1][slots.counter - 1] >= 2.0
            and abs(seeded_run.statev[target - 1][slots.iterate - 1] - gamma) <= 1e-12 * max(abs(gamma), 1.0))
    if not took:
        record["status"] = "seed_did_not_take_effect"
        return record

    ps_dir = work / "oti_probe"
    contract = GenericPSContract(
        name=f"{_safe(model)}_constitutive_jacobians", umat_source_path=seeded_src,
        parameters=(("GSEED", slots.seed_props),), parameter_values=(gamma,),
        state_variables=tuple((f"SDV{i + 1}", i + 1) for i in range(slots.nstatv)),
        ntens=ntens, nstatv=slots.nstatv, ndi=ndi, nshr=nshr,
        dstran_per_increment=tuple(path[0]), n_increments=len(path), static_props=tuple(probe_props))
    try:
        transform_umat_for_parameter_sensitivity(contract=contract, output_dir=ps_dir)
    except Exception as exc:  # noqa: BLE001 - reported
        record["status"] = "oti_transform_failed"
        record["reason"] = f"{type(exc).__name__}: {exc}"[:400]
        return record
    build = subprocess.run(["make"], cwd=ps_dir, capture_output=True, text=True)
    run = subprocess.run([str(ps_dir / "ps_driver")], cwd=ps_dir, capture_output=True, text=True) \
        if build.returncode == 0 else None
    if run is None or run.returncode != 0:
        record["status"] = "oti_build_or_run_failed"
        record["reason"] = ((build.stderr if build.returncode else run.stderr) or "")[:400]
        return record
    oti = read_oti_csv(ps_dir / "DSTATEV_DP_OTI.csv")
    ladder = {step: centered_fd(seeded_exe, probe_props, path, ntens=ntens, nstatv=slots.nstatv,
                                props_indices=[slots.seed_props], rel_step=step)[slots.seed_props]["dstatev"]
              for step in FD_STEPS}
    rows = []
    slot = base.nstatv + 1
    for symbol, response, expression, component in pairs:
        y_slot, j_slot = slot, slot + 1
        slot += 2
        hand = seeded_run.statev[target - 1][j_slot - 1]
        entry = oti.get((target, y_slot)) or {}
        oti_value = entry.get("GSEED", entry.get("gseed"))
        fd_values = [ladder[s][target - 1][y_slot - 1] for s in FD_STEPS]
        response = seeded_run.statev[target - 1][y_slot - 1]
        # consecutive-step agreement, but never below the round-off floor of the
        # finer step, eps*|Y|/(2h): two steps that quantise the response to the
        # same bits agree exactly and are not thereby converged
        scores = []
        for i, (a, b) in enumerate(zip(fd_values, fd_values[1:])):
            gap = abs(a - b) / max(abs(a), abs(b), 1e-300)
            h_fine = FD_STEPS[i + 1] * max(abs(gamma), 1e-300)
            floor = 2.220446049250313e-16 * max(abs(response), 1e-300) / (2.0 * h_fine) / max(abs(b), 1e-300)
            scores.append((max(gap, floor), gap, i))
        _score, best_gap, best_i = min(scores)
        fd = fd_values[best_i]
        rows.append({
            "symbol": symbol, "component": component, "response": response,
            "hand_coded": hand, "oti": oti_value, "fd": fd, "fd_step": FD_STEPS[best_i],
            "fd_plateau_gap": best_gap, "fd_plateau_score": _score,
            "oti_vs_hand_rel": _rel(oti_value, hand), "oti_vs_fd_rel": _rel(oti_value, fd),
            "hand_vs_fd_rel": _rel(hand, fd), "response_value": seeded_run.statev[target - 1][y_slot - 1],
        })
    record["rows"] = rows
    record["status"] = "measured"
    return record


def configurations(model: str, probe_props: List[float], ntens: int):
    """Material/path candidates, most authoritative first."""
    out = []
    try:
        from umat_oti.validation.actual_umat_higher_order_generic import MODELS
        if model in MODELS:
            spec = MODELS[model]
            out.append((f"repository material (actual_umat_higher_order_generic.MODELS['{model}'])",
                        [float(v) for v in spec.props], [list(spec.increments[0])] * len(spec.increments)))
    except ImportError:
        pass
    for label, component, magnitude in PROBE_PATHS:
        increment = [0.0] * ntens
        increment[component] = magnitude
        out.append((f"paired-validation probe constants, path {label}", list(probe_props),
                    [increment] * 4))
    return out


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


def _rel(a, b):
    if a is None or b is None:
        return None
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def verdict(rows: List[dict]) -> Dict[str, object]:
    """Slide classification for one (UMAT, symbol): worst over its components."""
    if not rows or any(r["oti"] is None for r in rows):
        return {"cell": "not measured"}
    worst = max(r["oti_vs_hand_rel"] for r in rows)
    cell = "Exact" if worst == 0.0 else ("Pass" if worst <= PASS_TOLERANCE else f"Differs {worst:.1e}")
    return {"cell": cell, "oti_vs_hand_rel_max": worst,
            "oti_vs_fd_rel_max": max(r["oti_vs_fd_rel"] for r in rows),
            "hand_vs_fd_rel_max": max(r["hand_vs_fd_rel"] for r in rows),
            "fd_plateau_gap_max": max(r["fd_plateau_gap"] for r in rows)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", action="append", dest="models")
    parser.add_argument("--claim4", type=Path, default=None,
                        help="claim4_benchmark_ddsdde.json for the CEVPI (HIN DDSDDE) cell")
    parser.add_argument("--out", type=Path, default=default_out())
    parser.add_argument("--work", type=Path, default=default_work() / "claim5")
    args = parser.parse_args(argv)
    origins = require_worktree_imports()
    work = args.work.resolve(); work.mkdir(parents=True, exist_ok=True)
    claim4_path = args.claim4 or (args.out / "claim4_benchmark_ddsdde.json")
    claim4 = json.loads(claim4_path.read_text()) if claim4_path.is_file() else None
    table, records = {}, []
    for model in (args.models or UMATS):
        source = umat_repo() / "UMATs" / "UMATs" / "ICP" / f"{model}.for"
        text = source.read_text(errors="replace")
        present = symbols_present(text)
        config = json.loads((umat_repo() / "benchmarks" / f"{model}.json").read_text())
        ntens_hint = NTENS_OVERRIDE.get(model, (int(config.get("ntens", 6)), ""))[0]
        ntens, nstatv, props = probe_material(source, ntens_hint)
        for index, value in PROPS_OVERRIDE.get(model, ({}, ""))[0].items():
            props[index - 1] = value
        print(f"[claim5] {model}: symbols {present}; ntens={ntens} nstatv={nstatv} nprops={len(props)}", flush=True)
        model_work = fresh_dir(work / _safe(model))
        executed, closure = resolved_source(source, model_work)
        record = {"model": model, "symbols_present": present, "ntens": ntens, "nstatv": nstatv,
                  "closure": closure,
                  "ntens_override": NTENS_OVERRIDE.get(model, (None, None))[1],
                  "props_override": PROPS_OVERRIDE.get(model, (None, None))[1]}
        local = [s for s in present if s != "CEVPI"]
        record["attempts"] = []
        if local:
            for label, cfg_props, cfg_path in configurations(model, props, ntens):
                attempt = {"material_and_path": label}
                try:
                    attempt.update(measure(model, executed, local, cfg_props, ntens, nstatv, cfg_path,
                                           model_work / f"attempt{len(record['attempts'])}"))
                except Exception as exc:  # noqa: BLE001 - recorded, never skipped
                    attempt["status"] = "error"
                    attempt["reason"] = f"{type(exc).__name__}: {exc}"[:600]
                record["attempts"].append({k: v for k, v in attempt.items() if k != "rows"})
                if attempt.get("status") == "measured":
                    record.update(attempt)
                    record["props"], record["path"] = cfg_props, {
                        "dstran_per_increment": cfg_path[0], "n_increments": len(cfg_path)}
                    break
            else:
                record["status"] = record["attempts"][-1].get("status")
        cells = {}
        for symbol in SYMBOLS:
            if symbol not in present:
                cells[symbol] = {"cell": "-"}
            elif symbol == "CEVPI":
                cells[symbol] = _cevpi_cell(model, claim4)
            else:
                rows = [r for r in record.get("rows", []) if r["symbol"] == symbol]
                cells[symbol] = verdict(rows) if rows else {"cell": f"not measured ({record.get('status')})"}
        record["cells"] = cells
        table[model] = {s: cells[s]["cell"] for s in SYMBOLS}
        print("          " + "  ".join(f"{s}={cells[s]['cell']}" for s in SYMBOLS), flush=True)
        records.append(record)
    measured = [(m, s) for m, row in table.items() for s, c in row.items() if c not in ("-",) and not c.startswith("not")]
    summary = {
        "symbol_model_pairs_in_sources": sum(1 for row in table.values() for c in row.values() if c != "-"),
        "pairs_measured": len(measured),
        "exact": sum(1 for m, s in measured if table[m][s].split(" ")[0] == "Exact"),
        "pass": sum(1 for m, s in measured if table[m][s].split(" ")[0] == "Pass"),
        "differs": [f"{m}/{s}: {table[m][s]}" for m, s in measured if table[m][s].startswith("Differs")],
        "not_measured": [f"{m}/{s}: {c}" for m, row in table.items() for s, c in row.items() if c.startswith("not")],
        "slide": {"entries": 19, "exact": 14, "pass": 5},
        "slide_cells_not_on_slide": [f"{m}/{s}" for m, row in table.items() for s, c in row.items()
                                     if c != "-" and s not in SLIDE.get(m, {})],
        "table": table, "slide_table": SLIDE,
    }
    payload = {"claim": "slide 25", "summary": summary, "records": records,
               "provenance": provenance(), "imports": origins}
    write_json(args.out.resolve() / "claim5_constitutive_jacobians.json", payload)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("table", "slide_table")}, indent=1))
    return 0


def _cevpi_cell(model: str, claim4) -> Dict[str, object]:
    if not claim4:
        return {"cell": "not measured (claim 4 result absent)"}
    rows = [r for r in claim4.get("rows", []) if r.get("case") == model]
    main = next((r for r in rows if not r.get("variant")), None)
    variant = next((r for r in rows if r.get("variant") and r.get("compared")), None)
    chosen, note = main, "claim 4 Abaqus DDSDDE = (1-D) CEVPI, committed contract"
    if (not main or not main.get("compared")) and variant:
        chosen, note = variant, ("claim 4 Abaqus DDSDDE = (1-D) CEVPI, documented contract variant: "
                                 + variant["variant"])
    if not chosen or not chosen.get("compared"):
        return {"cell": f"not measured (claim 4: {(main or {}).get('status')})",
                "reason": "; ".join((main or {}).get("blockers", []))[:400]}
    rel = chosen["compared"].get("ddsdde_max_rel")
    cell = "Exact" if rel == 0.0 else ("Pass" if rel is not None and rel <= PASS_TOLERANCE else f"Differs {rel}")
    return {"cell": cell + (" (variant)" if chosen is variant else ""), "source": note, "ddsdde_max_rel": rel}


if __name__ == "__main__":
    sys.exit(main())
