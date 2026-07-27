"""Neutral model exchange format (JSON).

A solver-independent serialization of the core ``Model`` so a model can be built
once (e.g. from an Abaqus ``.inp``) and re-loaded without re-parsing, or authored
by hand / another front-end. The format is deliberately physics-free and stores:

    - nodes (coordinates)
    - elements + element types + connectivity
    - element -> formulation / material bindings, node & element sets
    - materials as neutral descriptors (name, user_material flag, constants /
      PROPS, state size, section dict) — material *instances* are not serialized;
      the material registry re-binds concrete behavior on load
    - boundary / load / equation metadata (as attribute-accessible records)
    - optional ``field_refs``: references (paths) to external field-data exports
      (e.g. {"stress_ip": "fields.json"}) — the field data itself is not inlined

    save_model(model, path)
    model = load_model(path)   # -> core.model.Model with bindings unresolved
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from typing import Any, Dict

from ..core.model import Model, Element
from ..materials.base import MaterialBinding


SCHEMA = "resasm-neutral-model/1"


def _obj_to_dict(o: Any) -> Any:
    if o is None or isinstance(o, (str, int, float, bool)):
        return o
    if isinstance(o, dict):
        return {k: _obj_to_dict(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_obj_to_dict(v) for v in o]
    if hasattr(o, "__dict__"):
        return {k: _obj_to_dict(v) for k, v in vars(o).items()
                if not k.startswith("_")}
    return str(o)


def _material_to_dict(m: Any) -> Dict[str, Any]:
    if isinstance(m, MaterialBinding):
        return {"kind": "binding",
                "name": m.name,
                "material": getattr(m.material, "name", None),
                "constants": list(m.constants or []),
                "n_state_vars": int(m.n_state_vars or 0),
                "section": _obj_to_dict(m.section)}
    # parser Material or raw dict
    return {"kind": "raw",
            "name": getattr(m, "name", None),
            "user_material": bool(getattr(m, "user_material", False)),
            "constants": list(getattr(m, "constants", []) or []),
            "depvar": getattr(m, "depvar", None),
            "data": _obj_to_dict(m) if not hasattr(m, "name") else None}


def model_to_dict(model: Model) -> Dict[str, Any]:
    return {
        "schema": SCHEMA,
        "ndim": model.ndim,
        "nodes": {str(nid): list(xyz) for nid, xyz in model.nodes.items()},
        "elements": {str(eid): {"eid": el.eid, "etype": el.etype,
                                "connectivity": list(el.connectivity)}
                     for eid, el in model.elements.items()},
        "node_sets": {k: list(v) for k, v in model.node_sets.items()},
        "element_sets": {k: list(v) for k, v in model.element_sets.items()},
        "element_formulation": {str(k): v for k, v in model.element_formulation.items()},
        "element_material": {str(k): v for k, v in model.element_material.items()},
        "materials": {name: _material_to_dict(m) for name, m in model.materials.items()},
        "boundaries": [_obj_to_dict(b) for b in model.boundaries],
        "cloads": [_obj_to_dict(c) for c in model.cloads],
        "equations": [_obj_to_dict(e) for e in model.equations],
        # optional references to external field-data exports (not the data itself)
        "field_refs": dict(model.meta.get("field_refs", {})),
        "meta": {k: v for k, v in _obj_to_dict(model.meta).items()
                 if k != "field_refs"},
    }


def save_model(model: Model, path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(model_to_dict(model), fh, indent=2)


def _material_from_dict(md: Dict[str, Any], name: str):
    """Rebuild a material descriptor that still exposes ``user_material`` /
    ``constants`` as attributes (so the inspector and wizard work unchanged)."""
    if md.get("kind") == "binding":
        b = MaterialBinding(material=None, constants=list(md.get("constants", [])),
                            n_state_vars=int(md.get("n_state_vars", 0)),
                            name=md.get("name") or name)
        b.section = md.get("section")
        return b
    # raw parser-style material -> namespace with the expected attributes
    return SimpleNamespace(
        name=md.get("name") or name,
        user_material=bool(md.get("user_material", False)),
        constants=list(md.get("constants", []) or []),
        depvar=md.get("depvar"),
        section=md.get("section"))


def _record_from_dict(d: Any):
    """Rebuild a boundary/load/equation record as an attribute-accessible object."""
    if isinstance(d, dict):
        return SimpleNamespace(**d)
    return d


def dict_to_model(d: Dict[str, Any]) -> Model:
    if d.get("schema") != SCHEMA:
        raise ValueError("unexpected schema %r (want %r)" % (d.get("schema"), SCHEMA))
    m = Model()
    m.ndim = int(d.get("ndim", 3))
    m.nodes = {int(nid): tuple(xyz) for nid, xyz in d["nodes"].items()}
    for eid_s, el in d["elements"].items():
        eid = int(eid_s)
        m.elements[eid] = Element(el["eid"], el["etype"], list(el["connectivity"]))
    m.node_sets = {k: list(v) for k, v in d.get("node_sets", {}).items()}
    m.element_sets = {k: list(v) for k, v in d.get("element_sets", {}).items()}
    m.element_formulation = {int(k): v for k, v in d.get("element_formulation", {}).items()}
    m.element_material = {int(k): v for k, v in d.get("element_material", {}).items()}
    m.materials = {name: _material_from_dict(md, name)
                   for name, md in d.get("materials", {}).items()}
    m.boundaries = [_record_from_dict(b) for b in d.get("boundaries", [])]
    m.cloads = [_record_from_dict(c) for c in d.get("cloads", [])]
    m.equations = [_record_from_dict(e) for e in d.get("equations", [])]
    m.meta = dict(d.get("meta", {}))
    field_refs = d.get("field_refs") or {}
    if field_refs:
        m.meta["field_refs"] = dict(field_refs)
    m.meta.setdefault("source", "neutral_model_io")
    return m


def load_model(path: str) -> Model:
    with open(path, "r", encoding="utf-8") as fh:
        return dict_to_model(json.load(fh))
