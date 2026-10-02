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
from typing import Dict, List, Optional

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


def probe(record: dict, case, out_dir: Path, calls: List[Dict]) -> dict:
    """Run ``calls`` (keyword dicts for ``CorpusProvider.regular``) through the
    snan- and zero-initialised builds; compare every output bit for bit."""
    try:
        snan, commands = build_probe(record, case, out_dir, "snan")
        zero, _ = build_probe(record, case, out_dir, "zero")
    except RuntimeError as error:
        return {"status": "not_run", "reason": str(error)}
    differing = []
    for index, kw in enumerate(calls):
        outs = []
        for lib in (snan, zero):
            try:
                outs.append(lib.regular(**kw))
            except MaterialCallError as error:
                outs.append({"error": np.frombuffer(str(error).encode(), dtype=np.uint8)})
        for name in sorted(set(outs[0]) | set(outs[1])):
            a, b = outs[0].get(name), outs[1].get(name)
            if a is None or b is None or np.ascontiguousarray(a).tobytes() != np.ascontiguousarray(b).tobytes():
                differing.append({"call": index, "output": name})
        if len(differing) > 20:
            break
    return {"status": "trip" if differing else "clean", "calls": len(calls),
            "differences": differing[:20],
            "flags": {k: _BASE + v for k, v in FLAG_SETS.items()},
            "statement": "the original routine rebuilt with opposite initialisation of every "
                         "local returns %s outputs on the same %d calls"
                         % ("DIFFERENT" if differing else "bit-identical", len(calls))}
