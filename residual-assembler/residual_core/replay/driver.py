"""Replay driver: turn a converged primal record + the OTI binary into the exact
per-IP fields the existing ``field_sensitivity`` solver consumes.

For every element and every integration point it derives the material-point
kinematics from the mesh and the converged displacement, calls the binary
through the C ABI, and collects the real tangent (DDSDDE) and the parameter
derivative ``d(sigma)/d(p)`` -- so the OTI arithmetic stays entirely inside the
binary and the tool only ever sees plain numbers.

Small-strain elastic v1: one evaluation per IP at the converged strain. The
per-IP loop and returned field shapes are the same ones a path-dependent,
increment-marching replay will fill, so the solver downstream never changes.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

from ..formulations import c3d8_kernel as kern
from .abi import MaterialABI

N_IP = 8
NTENS = 6


def replay_elastic_fields(
    model,
    record,
    package,
    abi: MaterialABI,
    params: List[str],
    dof_manager,
) -> Dict[str, Any]:
    """Replay every IP at the converged state.

    Returns a dict with:
      tangent_fields            {eid -> (8,6,6)}   DDSDDE  (for K)
      stress_derivative_fields  {pname -> {eid -> (8,6)}}  d sigma/d p (for R_,p)
      stress_replay             {eid -> (8,6)}     reconstructed primal stress
      seed_indices              the 1-based PROPS indices, in ``params`` order
    """
    u = record.converged_u()
    props = record.props
    seed_indices = [package.param_index(name) for name in params]
    pts = kern.ABAQUS_C3D8_GAUSS.points

    tangent_fields: Dict[int, np.ndarray] = {}
    stress_replay: Dict[int, np.ndarray] = {}
    dsig: Dict[str, Dict[int, np.ndarray]] = {p: {} for p in params}

    for el in model.elements.values():
        Xe = model.coords_of(el.connectivity)
        edofs = np.asarray(dof_manager.element_dofs(el.connectivity,
                                                    ("UX", "UY", "UZ")), int)
        ue = u[edofs]
        D_ip = np.zeros((N_IP, NTENS, NTENS))
        sig_ip = np.zeros((N_IP, NTENS))
        dsig_ip = {p: np.zeros((N_IP, NTENS)) for p in params}
        for k in range(N_IP):
            B = kern.b_matrix_reference(Xe, pts[k])[0]      # (6,24)
            eps = B @ ue
            out = abi.eval_point(props, seed_indices, eps)
            sig_ip[k] = out["stress"]
            D_ip[k] = out["ddsdde"]
            for c, p in enumerate(params):
                dsig_ip[p][k] = out["dstress_dseed"][:, c]
        tangent_fields[el.eid] = D_ip
        stress_replay[el.eid] = sig_ip
        for p in params:
            dsig[p][el.eid] = dsig_ip[p]

    return {
        "tangent_fields": tangent_fields,
        "stress_derivative_fields": dsig,
        "stress_replay": stress_replay,
        "seed_indices": seed_indices,
    }


def verify_replay_stress(stress_replay: Mapping[int, np.ndarray],
                         record) -> Dict[str, float]:
    """Compare the replayed primal stress against the production stress carried
    in the record (if any). Returns {max_abs, max_rel, checked} -- the core
    'reconstructed stress agrees with the saved production result' check."""
    prod = record.stress_ip()
    if prod is None:
        return {"checked": 0.0, "max_abs": float("nan"), "max_rel": float("nan")}
    max_abs = 0.0
    scale = 0.0
    for eid, sig in stress_replay.items():
        if eid not in prod:
            continue
        d = np.abs(sig - prod[eid])
        max_abs = max(max_abs, float(d.max()))
        scale = max(scale, float(np.abs(prod[eid]).max()))
    return {"checked": 1.0, "max_abs": max_abs,
            "max_rel": max_abs / scale if scale > 1e-30 else max_abs}


def assemble_internal_force(model, stress_ip: Mapping[int, np.ndarray],
                            dof_manager) -> np.ndarray:
    """Global internal force F_int = A_e int B^T sigma dV from per-IP stress."""
    f = np.zeros(dof_manager.ndof, dtype=float)
    for el in model.elements.values():
        Xe = model.coords_of(el.connectivity)
        edofs = np.asarray(dof_manager.element_dofs(el.connectivity,
                                                    ("UX", "UY", "UZ")), int)
        f[edofs] += kern.element_internal_force_small_strain(Xe, stress_ip[el.eid])
    return f
