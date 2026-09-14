"""Constraints: Dirichlet boundary conditions and (parsed) linear equations.

F_constraints in R = F_int - F_ext + F_constraints. For displacement-driven
problems the prescribed DOFs are removed from the free set and the assembled
internal force at those DOFs is the reaction. This module resolves *Boundary
specs (numeric ranges + symmetry keywords) into (global_dof -> value) and the
free/prescribed partition. Linear *Equation constraints are parsed but not yet
applied (documented limitation).

Consumes duck-typed boundary objects with attributes: target, dof_start,
dof_end, value, kind, amplitude -- so it is not tied to the Abaqus parser.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

# Abaqus symmetry/antisymmetry -> constrained TRANSLATIONAL dofs (1=x,2=y,3=z)
_SYMM_DOF = {"XSYMM": [1], "YSYMM": [2], "ZSYMM": [3],
             "XASYMM": [2, 3], "YASYMM": [1, 3], "ZASYMM": [1, 2],
             "ENCASTRE": [1, 2, 3], "PINNED": [1, 2, 3]}


def resolve_target_nodes(model, target) -> List[int]:
    if isinstance(target, int):
        return [target]
    if target in model.node_sets:
        return list(model.node_sets[target])
    for name, ids in model.node_sets.items():
        if name.upper() == str(target).upper():
            return list(ids)
    try:
        return [int(target)]
    except (TypeError, ValueError):
        return []


def active_boundaries(model, step=None) -> List:
    """The boundary conditions in force during one step, in Abaqus's order.

    Abaqus does not apply a deck's boundary blocks all at once. Conditions
    written before the first ``*STEP`` are initial; each step then either
    carries the previous set forward and modifies it (``OP=MOD``, the
    default) or replaces it outright (``OP=NEW``). A four-step verification
    deck read as one flat list therefore applies the LAST step's
    displacements to every increment of the first -- silently, because a
    dictionary of prescribed values has no way to say two steps disagreed.

    ``step=None`` keeps the historical behaviour -- every boundary in the
    file, last one wins -- so callers that never had a step to pass are not
    quietly given a different answer than before.
    """
    if step is None:
        return list(model.boundaries)
    step = int(step)
    active: List = []
    for index in range(0, step + 1):
        of_this_step = [b for b in model.boundaries
                        if int(getattr(b, "step", 0)) == index]
        if not of_this_step:
            continue
        if any(str(getattr(b, "op", "MOD")).upper() == "NEW"
               for b in of_this_step):
            active = list(of_this_step)
        else:
            active.extend(of_this_step)
    return active


def dirichlet_dofs(model, dof_manager, step=None) -> Dict[int, float]:
    """{global_dof_index -> prescribed value}. Symmetry BCs contribute 0.0.

    ``step`` selects which step's conditions are in force; see
    :func:`active_boundaries` for why that is not the same as all of them.
    """
    out: Dict[int, float] = {}
    for b in active_boundaries(model, step):
        nids = resolve_target_nodes(model, b.target)
        if getattr(b, "kind", "value") in _SYMM_DOF:
            dofs = _SYMM_DOF[b.kind]
            value = 0.0
        else:
            dofs = list(range(int(b.dof_start), int(b.dof_end) + 1))
            value = float(b.value)
        for nid in nids:
            if nid not in dof_manager.node_index:
                continue
            for d in dofs:
                if dof_manager.has_dof(nid, d):
                    out[dof_manager.dof_index(nid, d)] = value
    return out


def partition(model, dof_manager, step=None) -> Tuple[np.ndarray, np.ndarray, Dict[int, float]]:
    """Return (free_mask, prescribed_idx, prescribed_values_dict)."""
    pres = dirichlet_dofs(model, dof_manager, step)
    ndof = dof_manager.ndof
    free_mask = np.ones(ndof, dtype=bool)
    idx = np.array(sorted(pres.keys()), dtype=int)
    if idx.size:
        free_mask[idx] = False
    return free_mask, idx, pres
