"""Locate Program 1 (`umat-oti`) and import its canonical validation code.

The figure and table scripts under `results/` compare Program 2's assembled
sensitivities against Program 1's material derivatives, so they need both the
built material objects *and* the finite-difference methodology Program 1 uses.
Re-implementing that methodology here is how the two ended up disagreeing about
the FCC crystal model: the transformer's own check called it a failure while the
deck's table called it 2.5e-09, because each rolled its own FD step.  There is
now one implementation, and it lives in Program 1.

Resolution order for the Program 1 checkout:

1. ``UMAT_OTI_ROOT`` if set;
2. a sibling ``umat-oti/`` next to this repository (the two-repo layout);
3. the historical ``~/Documents/UMAT_source_transformation``.
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_RA = os.path.dirname(_HERE)

_CANDIDATES = [
    os.environ.get("UMAT_OTI_ROOT"),
    os.path.join(os.path.dirname(_RA), "umat-oti"),
    os.path.expanduser("~/Documents/UMAT_source_transformation"),
]


def program1_root() -> str:
    """Path to the umat-oti checkout, or raise with the places tried."""
    tried = []
    for cand in _CANDIDATES:
        if not cand:
            continue
        cand = os.path.abspath(os.path.expanduser(cand))
        tried.append(cand)
        if os.path.isdir(os.path.join(cand, "oti_provider", "materials")):
            return cand
    raise RuntimeError(
        "could not locate the umat-oti (Program 1) checkout.\n"
        "  set UMAT_OTI_ROOT, or place it as a sibling of this repository.\n"
        "  tried:\n%s" % "\n".join("    " + t for t in tried))


def fd_reference():
    """Program 1's canonical FD-reference module, imported from its checkout."""
    src = os.path.join(program1_root(), "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from umat_oti.validation import fd_reference as module
    return module


P1 = program1_root()
