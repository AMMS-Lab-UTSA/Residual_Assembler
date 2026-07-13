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


def dirichlet_dofs(model, dof_manager) -> Dict[int, float]:
    """{global_dof_index -> prescribed value}. Symmetry BCs contribute 0.0."""
    out: Dict[int, float] = {}
    for b in model.boundaries:
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


def partition(model, dof_manager) -> Tuple[np.ndarray, np.ndarray, Dict[int, float]]:
    """Return (free_mask, prescribed_idx, prescribed_values_dict)."""
    pres = dirichlet_dofs(model, dof_manager)
    ndof = dof_manager.ndof
    free_mask = np.ones(ndof, dtype=bool)
    idx = np.array(sorted(pres.keys()), dtype=int)
    if idx.size:
        free_mask[idx] = False
    return free_mask, idx, pres
