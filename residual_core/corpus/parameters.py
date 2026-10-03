"""Which PROPS slot is which material constant, and the domain it is declared on.

A finite-difference ladder step p +/- h is evidence only where the ORIGINAL is
defined: Poisson's ratio on (-1, 1/2), a modulus above 0. A step that leaves
that domain (nu = 0.499 + 0.05 crosses 1/2, where the bulk modulus changes
sign) evaluates a different material, or none. Such steps are dropped from the
ladder and tagged ``outside_parameter_domain``, for every source alike (Vera,
B8 review of f7be16bc).

The props map is read from the source itself: an assignment ``NAME = PROPS(i)``
whose NAME is in ``DECLARED`` gives slot i that name's domain. A slot no
declared name reaches has no declared domain and keeps its whole ladder. The
map, with the source line each entry came from, goes into every record.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

__all__ = ["DECLARED", "props_map", "domain_of", "inside"]

#: name (upper case) -> (domain label, lower bound, upper bound); bounds are open
_POISSON = ("poisson_ratio", -1.0, 0.5)
_MODULUS = ("positive_modulus", 0.0, float("inf"))
DECLARED: Dict[str, Tuple[str, float, float]] = {
    **{n: _POISSON for n in ("NU", "ENU", "XNU", "ANU", "PNU", "POISSON", "POISSON_RATIO",
                             "POISSION_RATIO", "PR")},
    **{n: _MODULUS for n in ("E", "EMOD", "YOUNG", "YOUNGS", "YM", "ELASTIC_MODULUS", "C0",
                             "C10", "MU", "G", "GMOD", "SHEAR_MODULUS", "K", "BULK",
                             "BULK_MODULUS")},
}

_ASSIGN = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*PROPS\s*\(\s*(\d+)\s*\)\s*(?:!.*)?$",
                     re.IGNORECASE)


def props_map(source_text: str, source_form: str = "fixed") -> List[dict]:
    """Every ``NAME = PROPS(i)`` assignment, with its line and, when NAME is
    declared, its domain."""
    out = []
    for number, line in enumerate(source_text.splitlines(), start=1):
        if source_form == "fixed" and line[:1] in ("c", "C", "*", "!"):
            continue
        if source_form != "fixed" and line.lstrip().startswith("!"):
            continue
        match = _ASSIGN.match(line)
        if not match:
            continue
        name, index = match.group(1).upper(), int(match.group(2))
        entry = {"index": index, "name": name, "line": number, "code": line.strip()}
        if name in DECLARED:
            label, lo, hi = DECLARED[name]
            entry["domain"] = {"kind": label, "open_interval": [lo, hi]}
        out.append(entry)
    return out


def domain_of(mapping: List[dict], slot: int) -> Optional[Tuple[float, float, str]]:
    """(lo, hi, basis) of PROPS(slot), or None. Two declared names with
    different domains on one slot: the intersection."""
    hits = [m for m in mapping if m["index"] == slot and m.get("domain")]
    if not hits:
        return None
    lo = max(m["domain"]["open_interval"][0] for m in hits)
    hi = min(m["domain"]["open_interval"][1] for m in hits)
    basis = "; ".join("%s (line %d: %s)" % (m["domain"]["kind"], m["line"], m["code"]) for m in hits)
    return lo, hi, basis


def inside(value: float, domain: Optional[Tuple[float, float, str]]) -> bool:
    if domain is None:
        return True
    return domain[0] < value < domain[1]
