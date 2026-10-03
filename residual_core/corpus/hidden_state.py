"""Does the ORIGINAL routine read memory it never wrote? (Vera B1 review, A)

Finite differences of the original are an independent reference only if its
outputs are a function of its arguments. Restoring every argument does not
establish that: a local read before it is assigned (or a SAVEd value) makes
the output depend on whatever the memory held. Two checks, both on ALL outputs
and both bit-exact:

* replays (in ``runner``): h = 0 calls of the original after the perturbed
  ones, and a nominal re-run of the whole history after every perturbed
  re-solve;
* this probe: the source (exactly as the provider compiled it) is rebuilt
  twice with gfortran -O0 -fcheck=bounds and opposite initialisation of every
  local (``-finit-real=snan -finit-integer=-77777 -finit-logical=true`` vs
  ``-finit-real=zero -finit-integer=0 -finit-logical=false``, both
  ``-finit-local-zero`` siblings), then driven through the SAME call sequence.
  A routine that never reads an unassigned local returns identical bits from
  both; any difference is a trip, and the source is ``not_attempted``.

The probe libraries contain only the original routine (and the provider's own
Abaqus utility stubs); they never see OTI arithmetic.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from residual_core.runtime import load_shared_library

from .provider import _SHIM, CorpusProvider, MaterialCallError

__all__ = ["FLAG_SETS", "build_probe", "probe"]

FLAG_SETS = {
    "snan": ["-finit-real=snan", "-finit-integer=-77777", "-finit-logical=true",
             "-finit-character=35"],
    "zero": ["-finit-real=zero", "-finit-integer=0", "-finit-logical=false",
             "-finit-character=32"],
}
_BASE = ["gfortran", "-O0", "-fPIC", "-std=legacy", "-ffree-line-length-none", "-fcheck=bounds"]
_REGULAR_SHIM = _SHIM[:_SHIM.index("subroutine corpus_total")] + r"""
subroutine corpus_flush() bind(C, name="corpus_flush")
  flush(6)
end subroutine corpus_flush
"""


def _support_dir(record: dict) -> Optional[Path]:
    build = Path(record["object"]).parent
    for candidate in sorted(build.glob("build-*")):
        if (candidate / "abaqus_stubs.f90").is_file() and (candidate / "aba_param.inc").is_file():
            return candidate
    return None


class _ProbeProvider(CorpusProvider):
    """CorpusProvider over a probe library (``regular`` only)."""

    def __init__(self, record, case, workdir, library):
        self.record, self.case, self.workdir = record, case, Path(workdir)
        self.slots = [int(p["props_index"]) for p in record["parameters"]]
        self.nparam = len(self.slots)
        self.ntens, self.nstatv = int(case.ntens), int(case.nstatv)
        self.nprops = len(case.props)
        self.finite = bool(case.finite)
        self.has_sdvini = False
        self.reads_drot = False
        self._handle = load_shared_library(str(library))
        self.lib = self._handle.lib
        self.lib.corpus_regular.restype = None
        self.library_path = str(library)


def build_probe(record: dict, case, out_dir: Path, which: str):
    """Compile the provider's (adapted) source with FLAG_SETS[which]; returns a
    provider-like object with ``regular``, or raises RuntimeError(reason)."""
    out_dir = Path(out_dir) / ("probe_" + which)
    support = _support_dir(record)
    if support is None:
        raise RuntimeError("provider build directory has no abaqus_stubs.f90/aba_param.inc")
    source_dir = Path(record["object"]).parents[1] / "src"
    sources = [p for p in source_dir.iterdir() if p.suffix.lower() in (".f", ".for", ".f90", ".f77")]
    if len(sources) != 1:
        raise RuntimeError("expected one adapted source in %s, found %d" % (source_dir, len(sources)))
    source = sources[0]
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    for name in ("aba_param.inc", "ABA_PARAM.INC", "ABA_PARAM.inc", "aba_param.INC"):
        shutil.copy2(support / "aba_param.inc", out_dir / name)
    shutil.copy2(support / "abaqus_stubs.f90", out_dir / "abaqus_stubs.f90")
    shutil.copy2(source, out_dir / source.name)
    (out_dir / "probe_shim.f90").write_text(_REGULAR_SHIM, encoding="utf-8")
    form = ["-ffixed-form", "-ffixed-line-length-none"] if case.source_form == "fixed" \
        else ["-ffree-form"]
    flags = _BASE + FLAG_SETS[which]
    commands = [
        flags + ["-c", "abaqus_stubs.f90", "-o", "stubs.o"],
        flags + form + ["-I", ".", "-c", source.name, "-o", "original.o"],
        flags + ["-ffree-form", "-c", "probe_shim.f90", "-o", "shim.o"],
        ["gfortran", "-shared", "-Wl,--no-undefined", "shim.o", "original.o", "stubs.o",
         "-o", "probe_%s.so" % which],
    ]
    for command in commands:
        done = subprocess.run(command, cwd=out_dir, capture_output=True, text=True)
        if done.returncode:
            raise RuntimeError("probe build (%s) failed: %s" % (which, (done.stderr or done.stdout)[-1500:]))
    return _ProbeProvider(record, case, out_dir, out_dir / ("probe_%s.so" % which)), commands


def _masked(out: Dict, undefined_statev: Sequence[int]) -> Dict:
    """``out`` with the declared-undefined STATEV entries set to zero."""
    if not undefined_statev or "state" not in out:
        return out
    out = dict(out)
    state = np.array(out["state"], copy=True)
    state[..., [i - 1 for i in undefined_statev]] = 0.0
    out["state"] = state
    return out


def _call(lib, kw):
    try:
        return lib.regular(**kw)
    except MaterialCallError as error:
        return {"error": np.frombuffer(str(error).encode(), dtype=np.uint8)}


def _differing(a: Dict, b: Dict) -> List[str]:
    out = []
    for name in sorted(set(a) | set(b)):
        x, y = a.get(name), b.get(name)
        if x is None or y is None or np.ascontiguousarray(x).tobytes() != np.ascontiguousarray(y).tobytes():
            out.append(name)
    return out


#: two incoming values written into a declared-undefined STATEV entry to show
#: the routine never reads it (the outputs must not move by a single bit)
UNREAD_PROBE_VALUES = (0.0, -1.2345e30, 1.2345e30, float("nan"))


def probe(record: dict, case, out_dir: Path, calls: List[Dict],
          undefined_statev: Sequence[int] = ()) -> dict:
    """Run ``calls`` (keyword dicts for ``CorpusProvider.regular``) through the
    snan- and zero-initialised builds; compare every output bit for bit.

    ``undefined_statev``: STATEV entries (1-based) the upstream verification
    run established as undefined in the original (its zero/snan/inf init builds
    differ there and in no STRESS or DDSDDE entry). Those entries are left out
    of the comparison ONLY IF a second check shows the routine never reads them:
    every call is repeated on the zero build with each value of
    ``UNREAD_PROBE_VALUES`` written into those incoming entries, and every other
    output must stay bit-identical. Otherwise the source trips as before.
    """
    undefined_statev = sorted(int(i) for i in undefined_statev
                              if 1 <= int(i) <= int(getattr(case, "nstatv", 0) or 0))
    try:
        snan, commands = build_probe(record, case, out_dir, "snan")
        zero, _ = build_probe(record, case, out_dir, "zero")
    except RuntimeError as error:
        return {"status": "not_run", "reason": str(error)}
    differing, masked_only, read = [], 0, []
    for index, kw in enumerate(calls):
        a, b = _call(snan, kw), _call(zero, kw)
        names = _differing(a, b)
        if names and undefined_statev:
            confined = _differing(_masked(a, undefined_statev), _masked(b, undefined_statev))
            if not confined:
                masked_only += 1
            names = confined
        differing.extend({"call": index, "output": name} for name in names)
        if undefined_statev:
            ref = _masked(b, undefined_statev)
            for value in UNREAD_PROBE_VALUES:
                state = np.array(kw["state"], dtype=float, copy=True)
                state[..., [i - 1 for i in undefined_statev]] = value
                moved = _differing(ref, _masked(_call(zero, dict(kw, state=state)), undefined_statev))
                read.extend({"call": index, "output": name, "incoming_value": value}
                            for name in moved)
        if len(differing) + len(read) > 20:
            break
    trip = bool(differing or read)
    statement = ("the original routine rebuilt with opposite initialisation of every local "
                 "returns %s outputs on the same %d calls"
                 % ("DIFFERENT" if differing else "bit-identical", len(calls)))
    if undefined_statev:
        statement += ("%s STATEV%s (declared undefined in the original by the verification run) "
                      "%s" % (" except" if not differing else "; ", undefined_statev,
                              "is READ by the routine: changing its incoming value moves "
                              "other outputs" if read else
                              "is never read: writing %s into it moves no other output by a bit"
                              % (list(UNREAD_PROBE_VALUES),)))
    return {"status": "trip" if trip else "clean", "calls": len(calls),
            "differences": differing[:20],
            "declared_undefined_statev": undefined_statev,
            "calls_differing_only_in_declared_undefined": masked_only,
            "declared_undefined_read": read[:20],
            "flags": {k: _BASE + v for k, v in FLAG_SETS.items()},
            "statement": statement}
