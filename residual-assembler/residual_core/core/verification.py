"""Generic verification harness (levels defined in docs/verification_strategy.md).

These operate on the formulation-agnostic Assembler, so any backend inherits
them. The finite-difference tangent check (Level 4) in particular is fully
backend-independent: it perturbs the global DOFs, finite-differences the
assembled residual, and compares to the assembled tangent -- catching
DDSDDE->element-tangent mapping errors for ANY formulation/material.

Level ladder (each backend lists which apply, base.py::verification_levels):
  0  parser / model inventory
  1  zero-field / zero-load residual
  2  rigid-body / null-mode
  3  manufactured solution / patch test
  4  finite-difference tangent
  5  solver comparison vs Abaqus (exported fields -> residual -> reactions)
  6  material replay vs Abaqus ODB (stress/state history)
  7  full residual replay (material replay + assembly -> reactions)
"""

from __future__ import annotations

from typing import Optional

import numpy as np


def zero_field_residual(assembler, fields=None, options=None):
    """Level 1: at U=0 with a zero stress field (or a material that returns zero
    at zero strain), the assembled internal force must vanish; R = -F_external."""
    ndof = assembler.dm.ndof
    R, _, diag = assembler.assemble(np.zeros(ndof), fields=fields, options=options)
    return float(np.max(np.abs(R))), diag


def rigid_body_translation(assembler, translation=(1.0, -2.0, 0.5),
                           time=(0.0, 0.0), dtime=0.0, fields=None, options=None):
    """Level 2: a uniform rigid-body TRANSLATION produces zero strain, hence (for
    material-update backends with no external load) zero internal force. Returns
    max|R|. Applicable to material-driven formulations (Mode 2); not to
    stress-driven (Mode 1), whose stress is externally imposed.

    Note: a full rigid-ROTATION objectivity check (force rotates with R) is
    verified at the element-kernel level for the finite-strain path; see
    tests/cp_c3d8_umat/tangent_fd_check and c3d8_kernel.
    """
    t = np.asarray(translation, float)
    ndpn = assembler.dm.ndof_per_node
    if t.size != ndpn:
        t = np.resize(t, ndpn)
    U = np.tile(t, assembler.dm.ndof // ndpn)
    R, _, diag = assembler.assemble(U, time=time, dtime=dtime, fields=fields,
                                    options=options)
    return float(np.max(np.abs(R))), diag


def finite_difference_tangent(assembler, U, h_rel=1e-6, time=(0.0, 0.0),
                              dtime=0.0, fields=None, options=None,
                              solution_state=None, dof_subset=None):
    """Level 4: central-difference the assembled residual wrt U and compare to
    the assembled tangent K. Returns (abs_err, rel_fro, max_entry, K, Kfd).

    dof_subset limits the columns perturbed (for large meshes); default = all.
    NOTE: material state is read from the StateManager committed state and is
    NOT advanced here (the perturbations are about a fixed committed state), so
    this checks d R / d U at frozen history -- the consistent-tangent question.
    """
    U = np.asarray(U, float)
    ndof = U.size
    R0, K, _ = assembler.assemble(U, time=time, dtime=dtime, fields=fields,
                                  options=options, compute_tangent=True,
                                  solution_state=solution_state)
    if K is None:
        raise ValueError("backend returned no tangent; Level-4 not applicable")
    cols = range(ndof) if dof_subset is None else list(dof_subset)
    scale = max(np.linalg.norm(U) / max(np.sqrt(ndof), 1.0), 1.0)
    h = h_rel * scale
    Kfd = np.zeros((ndof, len(list(cols)))) if dof_subset is not None else np.zeros((ndof, ndof))
    for c, j in enumerate(cols):
        Up = U.copy(); Up[j] += h
        Um = U.copy(); Um[j] -= h
        Rp, _, _ = assembler.assemble(Up, time=time, dtime=dtime, fields=fields,
                                      options=options, solution_state=solution_state)
        Rm, _, _ = assembler.assemble(Um, time=time, dtime=dtime, fields=fields,
                                      options=options, solution_state=solution_state)
        Kfd[:, c] = (Rp - Rm) / (2 * h)
    Kref = K if dof_subset is None else K[:, list(cols)]
    abs_err = float(np.max(np.abs(Kfd - Kref)))
    rel_fro = float(np.linalg.norm(Kfd - Kref) / max(np.linalg.norm(Kref), 1e-30))
    return abs_err, rel_fro, abs_err, K, Kfd


def reaction_check(assembler, R, reactions_ref=None):
    """Level 5 helper: split R into free residual and reactions at prescribed
    DOFs; if a reference reaction vector is given, compare (both signs)."""
    free_res, pres_idx, reac = assembler.split(R)
    out = {"free_residual_Linf": float(np.max(np.abs(free_res)) if free_res.size else 0.0),
           "n_free": int(free_res.size), "n_prescribed": int(pres_idx.size)}
    if reactions_ref is not None and pres_idx.size:
        ref = np.asarray(reactions_ref, float)[pres_idx]
        ep, em = np.linalg.norm(reac - ref), np.linalg.norm(reac + ref)
        out["reaction_sign"] = "+" if ep <= em else "-"
        out["reaction_abs_err"] = float(min(ep, em))
        out["reaction_rel_err"] = float(min(ep, em) / (np.linalg.norm(ref) + 1e-30))
    return out
