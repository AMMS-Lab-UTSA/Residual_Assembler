"""Abaqus ``.inp`` parser for the Residual Assembler.

Implements ``parse_inp(path) -> AbaqusModel`` per ``residual_core/CONTRACT.md``
section 3.  The parser is deliberately narrow: it understands exactly the
keywords needed to assemble a residual for a single-instance, no-transform
C3D8 + UMAT model (the Step-1 milestone), and *tolerates* everything else
without crashing (unrecognised keywords are recorded in
``model.unsupported_keywords``; unparseable data lines are recorded in
``model.warnings``).

Design notes / conventions honoured:

* Keyword lines start with a single ``*``; comment lines start with ``**``;
  keyword and parameter names are case-insensitive; parameters are
  comma-separated ``key`` or ``key=value``.
* Part / Instance / Assembly are flattened into ONE global model.  Node and
  element ids are taken as global exactly as written (valid for these
  single-instance, no-transform files).  A non-identity ``*Instance``
  transform raises ``NotImplementedError`` rather than being silently dropped.
* ``element_material`` is resolved from a *snapshot* of the elset membership
  taken at the moment the ``*Solid Section`` is parsed.  This is required for
  correctness: an elset name (e.g. ``Set-1``) may be defined once inside the
  part (all elements -> used by the section) and again at assembly level with
  a strided subset.  Resolving lazily against the final (last-wins) dict would
  give the wrong element set for the section.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union


# --------------------------------------------------------------------------- #
# Data model (fields per CONTRACT.md section 3)
# --------------------------------------------------------------------------- #
@dataclass
class Element:
    eid: int
    etype: str
    connectivity: List[int]


@dataclass
class Material:
    name: str
    user_material: bool = False
    constants: List[float] = field(default_factory=list)
    depvar: Optional[int] = None


@dataclass
class Section:
    elset: str
    material: str
    kind: str = "solid"


@dataclass
class Boundary:
    target: Union[str, int]
    dof_start: int
    dof_end: int
    value: float
    kind: str
    amplitude: Optional[str] = None


@dataclass
class Cload:
    target: Union[str, int]
    dof: int
    value: float
    amplitude: Optional[str] = None


@dataclass
class AbaqusModel:
    nodes: Dict[int, Tuple[float, float, float]] = field(default_factory=dict)
    elements: Dict[int, Element] = field(default_factory=dict)
    node_sets: Dict[str, List[int]] = field(default_factory=dict)
    element_sets: Dict[str, List[int]] = field(default_factory=dict)
    materials: Dict[str, Material] = field(default_factory=dict)
    sections: List[Section] = field(default_factory=list)
    element_material: Dict[int, str] = field(default_factory=dict)
    boundaries: List[Boundary] = field(default_factory=list)
    cloads: List[Cload] = field(default_factory=list)
    dsloads: List[dict] = field(default_factory=list)
    equations: List[dict] = field(default_factory=list)
    includes: List[str] = field(default_factory=list)
    instances: List[str] = field(default_factory=list)
    unsupported_keywords: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------ #
    @property
    def node_id_to_index(self) -> Dict[int, int]:
        """Map node id -> 0-based index, nodes sorted ascending by id."""
        return {nid: i for i, nid in enumerate(sorted(self.nodes))}

    def material_by_name(self, name: str) -> Optional[Material]:
        """Case-insensitive material lookup (Abaqus names are case-folded)."""
        key = name.strip().upper()
        for mname, mat in self.materials.items():
            if mname.upper() == key:
                return mat
        return None

    def to_dict(self) -> dict:
        """Return a plain, JSON-serialisable representation."""
        return {
            "nodes": {str(nid): list(xyz) for nid, xyz in self.nodes.items()},
            "elements": {
                str(eid): {
                    "eid": el.eid,
                    "etype": el.etype,
                    "connectivity": list(el.connectivity),
                }
                for eid, el in self.elements.items()
            },
            "node_sets": {k: list(v) for k, v in self.node_sets.items()},
            "element_sets": {k: list(v) for k, v in self.element_sets.items()},
            "materials": {
                k: {
                    "name": m.name,
                    "user_material": m.user_material,
                    "constants": list(m.constants),
                    "depvar": m.depvar,
                }
                for k, m in self.materials.items()
            },
            "sections": [
                {"elset": s.elset, "material": s.material, "kind": s.kind}
                for s in self.sections
            ],
            "element_material": {str(k): v for k, v in self.element_material.items()},
            "boundaries": [
                {
                    "target": b.target,
                    "dof_start": b.dof_start,
                    "dof_end": b.dof_end,
                    "value": b.value,
                    "kind": b.kind,
                    "amplitude": b.amplitude,
                }
                for b in self.boundaries
            ],
            "cloads": [
                {
                    "target": c.target,
                    "dof": c.dof,
                    "value": c.value,
                    "amplitude": c.amplitude,
                }
                for c in self.cloads
            ],
            "dsloads": list(self.dsloads),
            "equations": list(self.equations),
            "includes": list(self.includes),
            "instances": list(self.instances),
            "unsupported_keywords": list(self.unsupported_keywords),
            "warnings": list(self.warnings),
            "node_id_to_index": {str(k): v for k, v in self.node_id_to_index.items()},
        }


# --------------------------------------------------------------------------- #
# Symmetry boundary keyword -> fixed DOF range
# --------------------------------------------------------------------------- #
_SYMM_DOFS = {
    "XSYMM": (1, 1),
    "YSYMM": (2, 2),
    "ZSYMM": (3, 3),
    "XASYMM": (2, 6),   # tolerated; rows other than 1 fixed (not the milestone)
    "YASYMM": (1, 6),
    "ZASYMM": (1, 6),
    "ENCASTRE": (1, 6),
    "PINNED": (1, 3),
}

_MAX_INCLUDE_DEPTH = 64


# --------------------------------------------------------------------------- #
# Low-level line helpers
# --------------------------------------------------------------------------- #
def _is_comment(stripped: str) -> bool:
    return stripped[:2] == "**"


def _is_keyword(stripped: str) -> bool:
    return stripped[:1] == "*" and stripped[:2] != "**"


def _parse_keyword_line(stripped: str) -> Tuple[str, str, Dict[str, Any]]:
    """Split a ``*Keyword, key=value, flag`` line.

    Returns ``(keyword_norm, keyword_raw, params)`` where ``keyword_norm`` is
    lower-cased with collapsed internal whitespace, ``keyword_raw`` preserves
    the original spelling, and ``params`` keys are lower-cased (values keep
    their original case; bare flags map to ``True``).
    """
    body = stripped[1:]  # drop leading '*'
    parts = body.split(",")
    keyword_raw = parts[0].strip()
    keyword_norm = " ".join(keyword_raw.lower().split())
    params: Dict[str, Any] = {}
    for p in parts[1:]:
        p = p.strip()
        if not p:
            continue
        if "=" in p:
            k, v = p.split("=", 1)
            params[k.strip().lower()] = v.strip()
        else:
            params[p.lower()] = True
    return keyword_norm, keyword_raw, params


def _tokens(line: str) -> List[str]:
    """Comma-split a data line, stripping whitespace and dropping empties
    (handles trailing commas)."""
    return [t.strip() for t in line.split(",") if t.strip() != ""]


def _to_float(tok: str) -> float:
    return float(tok)


def _to_int(tok: str) -> int:
    # tolerate ids written as floats ("12.0")
    if tok.endswith(".") or "." in tok or "e" in tok.lower():
        return int(round(float(tok)))
    return int(tok)


def _target(tok: str) -> Union[str, int]:
    """A BC/load target is an int node id if it parses as one, else a set name."""
    try:
        return int(tok)
    except ValueError:
        return tok


# --------------------------------------------------------------------------- #
# Stage 1 -- read raw lines and splice *Include files
# --------------------------------------------------------------------------- #
def _read_with_includes(path: str, model: AbaqusModel, stack: List[str]) -> List[str]:
    apath = os.path.abspath(path)
    if apath in stack:
        raise RecursionError(
            "Circular *Include detected: %s already open in %s" % (apath, stack)
        )
    if len(stack) >= _MAX_INCLUDE_DEPTH:
        raise RecursionError("*Include nesting exceeds %d levels" % _MAX_INCLUDE_DEPTH)

    stack.append(apath)
    base = os.path.dirname(apath)
    out: List[str] = []
    try:
        with open(apath, "r", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                s = raw.rstrip("\r\n")
                st = s.lstrip()
                if _is_keyword(st):
                    kw_norm, _kw_raw, params = _parse_keyword_line(st)
                    if kw_norm == "include":
                        inc = params.get("input") or params.get("inp")
                        if not inc:
                            model.warnings.append(
                                "*Include without input= parameter: %r" % s
                            )
                            continue
                        inc = inc.strip().strip('"').strip("'")
                        inc_path = inc if os.path.isabs(inc) else os.path.join(base, inc)
                        model.includes.append(os.path.abspath(inc_path))
                        # A missing / unreadable *Include target must be
                        # non-fatal (CONTRACT.md sec.3 "never crash"): record a
                        # warning and splice in nothing, continuing the parse.
                        try:
                            out.extend(_read_with_includes(inc_path, model, stack))
                        except OSError as exc:
                            model.warnings.append(
                                "*Include target not found: %s (%s)" % (inc_path, exc)
                            )
                        continue
                out.append(s)
    finally:
        stack.pop()
    return out


# --------------------------------------------------------------------------- #
# Stage 2 -- group lines into (keyword, params, data-lines) blocks
# --------------------------------------------------------------------------- #
def _iter_blocks(lines: List[str]):
    """Yield ``(keyword_norm, keyword_raw, params, data_lines)`` blocks.

    Comment lines and blank lines are skipped everywhere.  Data lines are the
    non-keyword, non-comment lines following a keyword up to the next keyword.
    """
    i = 0
    n = len(lines)
    # skip anything before the first keyword
    while i < n:
        st = lines[i].lstrip()
        if _is_keyword(st):
            break
        i += 1

    while i < n:
        st = lines[i].lstrip()
        if not _is_keyword(st):
            i += 1
            continue
        kw_norm, kw_raw, params = _parse_keyword_line(st)
        i += 1
        data: List[str] = []
        while i < n:
            s = lines[i]
            sst = s.lstrip()
            if sst == "":
                i += 1
                continue
            if _is_comment(sst):
                i += 1
                continue
            if _is_keyword(sst):
                break
            data.append(s)
            i += 1
        yield kw_norm, kw_raw, params, data


def _merge_continuations(data: List[str]) -> List[str]:
    """Merge Abaqus data-line continuations (a line ending in ',' continues
    onto the next data line)."""
    merged: List[str] = []
    buf = ""
    for line in data:
        piece = line.strip()
        if buf:
            buf = buf + "," + piece
        else:
            buf = piece
        if buf.endswith(","):
            buf = buf[:-1]  # keep accumulating; drop the trailing comma marker
            continue
        merged.append(buf)
        buf = ""
    if buf:
        merged.append(buf)
    return merged


# --------------------------------------------------------------------------- #
# parse_inp
# --------------------------------------------------------------------------- #
def parse_inp(path: str) -> AbaqusModel:
    model = AbaqusModel()
    lines = _read_with_includes(path, model, [])

    current_material: Optional[Material] = None
    pending_sections: List[Section] = []  # sections whose elset was not yet known
    seen_setnames_node: set = set()
    seen_setnames_elem: set = set()

    def _record_unsupported(kw_raw: str) -> None:
        if kw_raw not in model.unsupported_keywords:
            model.unsupported_keywords.append(kw_raw)

    def _resolve_section(sec: Section) -> bool:
        """Resolve element_material for a section against the CURRENT elset
        membership snapshot.  Returns True if the elset was found."""
        members = model.element_sets.get(sec.elset)
        if members is None:
            # try case-insensitive match
            for name, mem in model.element_sets.items():
                if name.upper() == sec.elset.upper():
                    members = mem
                    break
        if members is None:
            return False
        for eid in members:
            model.element_material[eid] = sec.material
        return True

    for kw_norm, kw_raw, params, data in _iter_blocks(lines):
        try:
            if kw_norm == "node":
                _handle_node(model, data)

            elif kw_norm == "element":
                _handle_element(model, params, data, seen=seen_setnames_elem)

            elif kw_norm == "nset":
                _handle_set(
                    model, params, data, node=True, seen=seen_setnames_node
                )

            elif kw_norm == "elset":
                _handle_set(
                    model, params, data, node=False, seen=seen_setnames_elem
                )

            elif kw_norm == "solid section":
                sec = Section(
                    elset=str(params.get("elset", "")),
                    material=str(params.get("material", "")),
                    kind="solid",
                )
                model.sections.append(sec)
                if not _resolve_section(sec):
                    pending_sections.append(sec)

            elif kw_norm == "material":
                name = str(params.get("name", "material_%d" % (len(model.materials) + 1)))
                mat = Material(name=name)
                model.materials[name] = mat
                current_material = mat

            elif kw_norm == "depvar":
                _handle_depvar(model, current_material, data)

            elif kw_norm == "user material":
                _handle_user_material(model, current_material, params, data)

            elif kw_norm == "boundary":
                _handle_boundary(model, params, data)

            elif kw_norm == "cload":
                _handle_cload(model, params, data)

            elif kw_norm == "dsload":
                _handle_dsload(model, params, data)

            elif kw_norm == "equation":
                _handle_equation(model, data)

            elif kw_norm == "instance":
                _handle_instance(model, params, data)

            elif kw_norm in (
                "part",
                "end part",
                "assembly",
                "end assembly",
                "end instance",
            ):
                # structural markers -- flatten into one global model
                pass

            else:
                _record_unsupported(kw_raw)
        except NotImplementedError:
            raise
        except Exception as exc:  # never crash on a real file
            model.warnings.append(
                "Failed to parse *%s block: %s" % (kw_raw, exc)
            )

    # Final pass: resolve any section whose elset became known only later.
    for sec in pending_sections:
        if not _resolve_section(sec):
            model.warnings.append(
                "Section elset %r not found; element_material unresolved" % sec.elset
            )

    return model


# --------------------------------------------------------------------------- #
# Keyword handlers
# --------------------------------------------------------------------------- #
def _handle_node(model: AbaqusModel, data: List[str]) -> None:
    for line in data:
        tok = _tokens(line)
        if len(tok) < 2:
            model.warnings.append("Skipping malformed *Node line: %r" % line)
            continue
        nid = _to_int(tok[0])
        coords = [0.0, 0.0, 0.0]
        for k in range(min(3, len(tok) - 1)):
            coords[k] = _to_float(tok[1 + k])
        if nid in model.nodes:
            model.warnings.append("Duplicate node id %d (last wins)" % nid)
        model.nodes[nid] = (coords[0], coords[1], coords[2])


def _handle_element(
    model: AbaqusModel, params: dict, data: List[str], seen: Optional[set] = None
) -> None:
    etype = str(params.get("type", "")).upper()
    elset_name = params.get("elset")
    eids: List[int] = []
    for line in _merge_continuations(data):
        tok = _tokens(line)
        if len(tok) < 2:
            model.warnings.append("Skipping malformed *Element line: %r" % line)
            continue
        eid = _to_int(tok[0])
        conn = [_to_int(t) for t in tok[1:]]
        if eid in model.elements:
            model.warnings.append("Duplicate element id %d (last wins)" % eid)
        model.elements[eid] = Element(eid=eid, etype=etype, connectivity=conn)
        eids.append(eid)
    # An `elset=` on *Element implicitly defines that element set (Abaqus
    # semantics): a later *Solid Section, elset=NAME must be able to resolve it.
    # Register with the SAME last-wins / collision semantics as _handle_set so
    # that a later explicit *Elset, elset=NAME stays consistent.
    if elset_name is not None:
        name = str(elset_name)
        if seen is not None:
            if name in seen:
                model.warnings.append(
                    "Elset %r redefined (part/assembly scope collision); "
                    "last definition wins" % name
                )
            seen.add(name)
        model.element_sets[name] = eids


def _handle_set(
    model: AbaqusModel, params: dict, data: List[str], node: bool, seen: set
) -> None:
    target = model.node_sets if node else model.element_sets
    name = params.get("nset") if node else params.get("elset")
    if name is None:
        model.warnings.append("*%s without name" % ("Nset" if node else "Elset"))
        return
    name = str(name)
    ids: List[int] = []
    if params.get("generate"):
        for line in data:
            tok = _tokens(line)
            if len(tok) < 2:
                continue
            start = _to_int(tok[0])
            end = _to_int(tok[1])
            step = _to_int(tok[2]) if len(tok) >= 3 else 1
            if step == 0:
                step = 1
            ids.extend(range(start, end + 1, step))  # inclusive
    else:
        for line in data:
            for t in _tokens(line):
                try:
                    ids.append(_to_int(t))
                except ValueError:
                    # non-integer token in a set (e.g. a nested set name) -- keep raw
                    model.warnings.append(
                        "Non-integer token %r in set %r ignored" % (t, name)
                    )
    if name in seen:
        model.warnings.append(
            "%s %r redefined (part/assembly scope collision); last definition wins"
            % ("Nset" if node else "Elset", name)
        )
    seen.add(name)
    target[name] = ids


def _handle_depvar(
    model: AbaqusModel, current: Optional[Material], data: List[str]
) -> None:
    if current is None:
        model.warnings.append("*Depvar outside a *Material block")
        return
    for line in data:
        tok = _tokens(line)
        if tok:
            current.depvar = _to_int(tok[0])
            return


def _handle_user_material(
    model: AbaqusModel, current: Optional[Material], params: dict, data: List[str]
) -> None:
    if current is None:
        model.warnings.append("*User Material outside a *Material block")
        return
    current.user_material = True
    consts: List[float] = []
    for line in data:
        for t in _tokens(line):
            try:
                consts.append(_to_float(t))
            except ValueError:
                model.warnings.append(
                    "Non-numeric *User Material constant %r ignored" % t
                )
    current.constants = consts
    # 'constants=N' is advisory; warn on mismatch but keep what we parsed.
    n = params.get("constants")
    if n is not None:
        try:
            if int(n) != len(consts):
                model.warnings.append(
                    "*User Material constants=%s but %d values parsed for %r"
                    % (n, len(consts), current.name)
                )
        except ValueError:
            pass


def _handle_boundary(model: AbaqusModel, params: dict, data: List[str]) -> None:
    amp = params.get("amplitude")
    if amp is not None:
        amp = str(amp)
    for line in data:
        tok = _tokens(line)
        if len(tok) < 2:
            model.warnings.append("Skipping malformed *Boundary line: %r" % line)
            continue
        target = _target(tok[0])
        second = tok[1]
        sym = second.strip().upper()
        if sym in _SYMM_DOFS:
            d0, d1 = _SYMM_DOFS[sym]
            model.boundaries.append(
                Boundary(
                    target=target,
                    dof_start=d0,
                    dof_end=d1,
                    value=0.0,
                    kind=sym,
                    amplitude=amp,
                )
            )
        else:
            dof_start = _to_int(tok[1])
            dof_end = _to_int(tok[2]) if len(tok) >= 3 else dof_start
            value = _to_float(tok[3]) if len(tok) >= 4 else 0.0
            model.boundaries.append(
                Boundary(
                    target=target,
                    dof_start=dof_start,
                    dof_end=dof_end,
                    value=value,
                    kind="value",
                    amplitude=amp,
                )
            )


def _handle_cload(model: AbaqusModel, params: dict, data: List[str]) -> None:
    amp = params.get("amplitude")
    if amp is not None:
        amp = str(amp)
    for line in data:
        tok = _tokens(line)
        if len(tok) < 3:
            model.warnings.append("Skipping malformed *Cload line: %r" % line)
            continue
        model.cloads.append(
            Cload(
                target=_target(tok[0]),
                dof=_to_int(tok[1]),
                value=_to_float(tok[2]),
                amplitude=amp,
            )
        )


def _handle_dsload(model: AbaqusModel, params: dict, data: List[str]) -> None:
    for line in data:
        tok = _tokens(line)
        if not tok:
            continue
        entry = {
            "target": _target(tok[0]),
            "type": tok[1] if len(tok) >= 2 else None,
            "values": [],
        }
        for t in tok[2:]:
            try:
                entry["values"].append(_to_float(t))
            except ValueError:
                entry["values"].append(t)
        model.dsloads.append(entry)


def _handle_equation(model: AbaqusModel, data: List[str]) -> None:
    # *Equation: first data line = number of terms; following lines list the
    # terms as (node/set, dof, coefficient) triples.  Parsed minimally (rows).
    lines = _merge_continuations(data)
    if not lines:
        return
    try:
        nterms = _to_int(_tokens(lines[0])[0])
    except (ValueError, IndexError):
        nterms = None
    terms = []
    for line in lines[1:]:
        tok = _tokens(line)
        i = 0
        while i + 2 < len(tok) + 1 and i + 1 < len(tok):
            try:
                terms.append(
                    {
                        "target": _target(tok[i]),
                        "dof": _to_int(tok[i + 1]),
                        "coeff": _to_float(tok[i + 2]) if i + 2 < len(tok) else 0.0,
                    }
                )
            except (ValueError, IndexError):
                break
            i += 3
    model.equations.append({"nterms": nterms, "terms": terms})


def _handle_instance(model: AbaqusModel, params: dict, data: List[str]) -> None:
    name = str(params.get("name", "instance_%d" % (len(model.instances) + 1)))
    model.instances.append(name)
    # Data lines between *Instance and *End Instance encode a translation
    # (line 1) and optional rotation (line 2).  A non-identity transform is
    # NOT supported -- we must not silently ignore it.
    vals: List[float] = []
    for line in data:
        for t in _tokens(line):
            try:
                vals.append(_to_float(t))
            except ValueError:
                pass
    if vals and any(abs(v) > 0.0 for v in vals):
        raise NotImplementedError(
            "*Instance %r has a non-zero translation/rotation transform "
            "(%s); flattening with a real transform is not implemented. "
            "This parser only supports single-instance, no-transform models."
            % (name, vals)
        )


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def _example_path(*parts: str) -> str:
    here = os.path.dirname(os.path.abspath(__file__))     # residual_core/io
    root = os.path.dirname(os.path.dirname(here))          # repo root
    return os.path.join(root, "sources", "permissive",
                        "ngrilli_Oxford_Crystal_Plasticity", *parts)


def _selftest() -> None:
    inp = _example_path("ExampleInputFiles", "HCPnoTwin", "Compression111.inp")
    model = parse_inp(inp)

    # --- structural assertions ------------------------------------------- #
    assert len(model.nodes) == 216, "expected 216 nodes, got %d" % len(model.nodes)
    assert len(model.elements) == 125, (
        "expected 125 elements, got %d" % len(model.elements)
    )
    etypes = {e.etype.upper() for e in model.elements.values()}
    assert etypes == {"C3D8"}, "expected all C3D8, got %s" % etypes

    # --- material -------------------------------------------------------- #
    mat = model.material_by_name("CPURANIUM")
    assert mat is not None, "material CPURANIUM (case-insensitive) not found"
    assert mat.user_material is True, "expected user_material=True"
    assert mat.depvar == 125, "expected depvar==125, got %r" % mat.depvar
    assert len(mat.constants) == 11, (
        "expected 11 constants, got %d" % len(mat.constants)
    )
    assert mat.constants[0] == 0.0, "expected constants[0]==0.0, got %r" % mat.constants[0]

    # --- element_material resolution (section snapshot) ------------------ #
    assert model.element_material[1].upper() == "CPURANIUM", (
        "element_material[1]=%r" % model.element_material.get(1)
    )
    assert all(
        v.upper() == "CPURANIUM" for v in model.element_material.values()
    ), "not all elements resolved to CPURANIUM"
    assert len(model.element_material) == 125, (
        "expected 125 elements resolved, got %d" % len(model.element_material)
    )

    # --- node_id_to_index ------------------------------------------------ #
    n2i = model.node_id_to_index
    assert n2i[1] == 0 and n2i[216] == 215, "node_id_to_index mapping wrong"

    # --- boundary conditions --------------------------------------------- #
    def _find_bc(target, kind=None):
        for b in model.boundaries:
            if str(b.target) == target and (kind is None or b.kind == kind):
                return b
        return None

    b2 = _find_bc("Set-2", "XSYMM")
    assert b2 is not None and b2.dof_start == 1 and b2.dof_end == 1, "Set-2 XSYMM"
    b1 = _find_bc("Set-1", "YSYMM")
    assert b1 is not None and b1.dof_start == 2 and b1.dof_end == 2, "Set-1 YSYMM"
    b3 = _find_bc("Set-3", "ZSYMM")
    assert b3 is not None and b3.dof_start == 3 and b3.dof_end == 3, "Set-3 ZSYMM"
    b4 = _find_bc("Set-4", "value")
    assert (
        b4 is not None
        and b4.dof_start == 3
        and b4.dof_end == 3
        and b4.value == 1.0
        and b4.amplitude == "Amp-1"
    ), "Set-4 prescribed dof3=1.0 amp Amp-1: got %r" % (b4,)

    # --- JSON round-trip -------------------------------------------------- #
    json.dumps(model.to_dict())

    print("PASS  Compression111.inp")
    print("  nodes=%d  elements=%d (types=%s)"
          % (len(model.nodes), len(model.elements), ",".join(sorted(etypes))))
    print("  material=%s  user_material=%s  depvar=%s  constants=%d  constants[0]=%s"
          % (mat.name, mat.user_material, mat.depvar, len(mat.constants),
             mat.constants[0]))
    print("  element_material[1]=%s  resolved=%d"
          % (model.element_material[1], len(model.element_material)))
    print("  boundaries=%d  node_sets=%d  element_sets=%d  unsupported=%s"
          % (len(model.boundaries), len(model.node_sets),
             len(model.element_sets), model.unsupported_keywords))
    if model.warnings:
        print("  warnings(%d): %s%s"
              % (len(model.warnings), model.warnings[:3],
                 " ..." if len(model.warnings) > 3 else ""))


def _smoke(rel_parts: Tuple[str, ...]) -> None:
    inp = _example_path(*rel_parts)
    if not os.path.exists(inp):
        print("SKIP  %s (not found)" % os.path.join(*rel_parts))
        return
    model = parse_inp(inp)
    etypes = {e.etype.upper() for e in model.elements.values()}
    print("PASS  %s" % os.path.join(*rel_parts))
    print("  nodes=%d  elements=%d  element-types=%s  materials=%d"
          % (len(model.nodes), len(model.elements),
             ",".join(sorted(etypes)) or "-", len(model.materials)))
    print("  element_material resolved=%d  sections=%d  boundaries=%d  includes=%d"
          % (len(model.element_material), len(model.sections),
             len(model.boundaries), len(model.includes)))
    if model.unsupported_keywords:
        print("  unsupported=%s" % model.unsupported_keywords)
    if model.warnings:
        print("  warnings=%d (e.g. %s)" % (len(model.warnings), model.warnings[0]))


if __name__ == "__main__":
    _selftest()
    print()
    _smoke(("InterfaceElementGeneration", "PyCiGen", "Job-1.inp"))
    print()
    _smoke(("ExampleInputFiles", "DiscreteTwin", "Job-1.inp"))
