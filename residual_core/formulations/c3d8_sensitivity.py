"""Structural parameter sensitivities on the C3D8 element path.

The element kernel assembles a residual. This turns it into a sensitivity: the
derivative of the nodal displacement field with respect to a material parameter,
obtained by differentiating equilibrium rather than by re-solving the model.

Equilibrium at the converged solution is

    R(u, p) = f_int(u, p) - f_ext = 0

so differentiating with respect to a material parameter p, with external load
independent of p, gives

    dR/dp = (partial R / partial u) du/dp + (partial R / partial p) = 0
          =>  K du/dp = - partial R / partial p

The explicit term is the one the material supplies:

    partial R / partial p = integral over the element of B^T (d sigma / d p) dV

which is the internal-force assembly again with ``d sigma / d p`` in place of
``sigma``. That is the whole point of the arrangement: the same quadrature and
the same B matrices are reused, and the only new input is the stress sensitivity
the transformed UMAT already produces at every integration point.

Sign convention, stated once because it is the easiest thing to get wrong:
``dR_dp`` below is ``+partial R / partial p``, and the solve applies the minus.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np

from residual_core.formulations.c3d8_kernel import (
    ABAQUS_C3D8_GAUSS, b_matrix_reference, b_matrix_spatial,
)

__all__ = [
    "element_dR_dp",
    "assemble_dR_dp",
    "apply_dirichlet",
    "solve_du_dp",
    "SensitivityResult",
]


def element_dR_dp(Xe, dsigma_dp_ip, Ue=None, *, mode: str = "small",
                  gauss=ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """Explicit parameter derivative of one element's internal force.

        dR_e/dp = sum_k B_k^T (d sigma_k / d p) detJ_k w_k

    ``dsigma_dp_ip`` is (8, 6): the Voigt stress sensitivity at each of the
    eight integration points, which is exactly what the transformed UMAT
    reports. The B matrices and quadrature are the element's own, so a
    sensitivity cannot silently use a different geometry from the residual it
    belongs to.
    """
    Xe = np.asarray(Xe, dtype=float)
    dsigma = np.asarray(dsigma_dp_ip, dtype=float)
    if dsigma.shape != (len(gauss.weights), 6):
        raise ValueError(
            f"dsigma_dp_ip must be ({len(gauss.weights)}, 6), got {dsigma.shape}")
    if mode == "finite":
        if Ue is None:
            raise ValueError("finite-strain mode needs the element displacement")
        xe = Xe + np.asarray(Ue, dtype=float).reshape(8, 3)
    out = np.zeros(24, dtype=float)
    for k, weight in enumerate(gauss.weights):
        if mode == "finite":
            B, detJ = b_matrix_spatial(xe, gauss.points[k])
        else:
            B, detJ = b_matrix_reference(Xe, gauss.points[k])
        out += (B.T @ dsigma[k]) * detJ * weight
    return out


def assemble_dR_dp(node_ids: Sequence[int], coords, connectivity,
                   dsigma_dp_all: Mapping[int, np.ndarray], U=None, *,
                   mode: str = "small") -> tuple[np.ndarray, dict[int, int]]:
    """Scatter element parameter derivatives into a global vector.

    Mirrors ``assemble_global_internal_force`` exactly, including its global DOF
    ordering, so the assembled sensitivity and the assembled residual are
    indexed identically and the linear solve is meaningful.
    """
    sorted_ids = sorted(node_ids)
    node_id_to_index = {nid: i for i, nid in enumerate(sorted_ids)}
    coords = np.asarray(coords, dtype=float)
    id_to_row = {nid: i for i, nid in enumerate(list(node_ids))}
    out = np.zeros(3 * len(sorted_ids), dtype=float)

    is_dict_U = isinstance(U, dict)
    U_arr = None if (is_dict_U or U is None) else np.asarray(U, dtype=float)

    for eid, conn in connectivity:
        Xe = np.array([coords[id_to_row[nid]] for nid in conn], dtype=float)
        Ue = None
        if mode == "finite":
            if is_dict_U:
                Ue = np.array([U[nid] for nid in conn], dtype=float)
            elif U_arr is not None:
                Ue = np.array([U_arr[3 * node_id_to_index[nid]:
                                     3 * node_id_to_index[nid] + 3]
                               for nid in conn], dtype=float)
        element = element_dR_dp(Xe, dsigma_dp_all[eid], Ue, mode=mode)
        for local, nid in enumerate(conn):
            gi = node_id_to_index[nid]
            out[3 * gi:3 * gi + 3] += element[3 * local:3 * local + 3]
    return out, node_id_to_index


def apply_dirichlet(K: np.ndarray, rhs: np.ndarray,
                    constrained_dofs: Iterable[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Eliminate constrained degrees of freedom by removing their rows/columns.

    The prescribed displacement does not depend on the material parameter, so
    its sensitivity is zero and elimination is exact rather than a penalty
    approximation. Returns the reduced system and the free-DOF index array.
    """
    ndof = rhs.shape[0]
    constrained = np.asarray(sorted(set(int(d) for d in constrained_dofs)), dtype=int)
    if constrained.size and (constrained.min() < 0 or constrained.max() >= ndof):
        raise ValueError("a constrained DOF index is outside the system")
    free = np.setdiff1d(np.arange(ndof), constrained)
    return K[np.ix_(free, free)], rhs[free], free


@dataclass
class SensitivityResult:
    """du/dp on the full DOF vector, with the evidence behind the solve."""

    du_dp: np.ndarray
    dR_dp: np.ndarray
    free_dofs: np.ndarray
    constrained_dofs: np.ndarray
    residual_norm: float
    condition_estimate: Optional[float] = None
    notes: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "ndof": int(self.du_dp.shape[0]),
            "free_dofs": int(self.free_dofs.size),
            "constrained_dofs": int(self.constrained_dofs.size),
            "linear_solve_residual_norm": float(self.residual_norm),
            "condition_estimate": (None if self.condition_estimate is None
                                   else float(self.condition_estimate)),
            "max_abs_du_dp": float(np.max(np.abs(self.du_dp))) if self.du_dp.size else 0.0,
            **self.notes,
        }


def solve_du_dp(K: np.ndarray, dR_dp: np.ndarray,
                constrained_dofs: Iterable[int], *,
                estimate_condition: bool = True) -> SensitivityResult:
    """Solve ``K du/dp = -dR/dp`` on the free degrees of freedom.

    The minus is applied here, once, so no caller has to remember it. Constrained
    DOFs keep a sensitivity of exactly zero: their displacement is prescribed and
    does not move when a material parameter does.

    The residual of the linear solve is returned rather than assumed small; a
    near-singular reduced stiffness produces a large one instead of a plausible
    answer.
    """
    K = np.asarray(K, dtype=float)
    dR_dp = np.asarray(dR_dp, dtype=float)
    if K.shape[0] != K.shape[1] or K.shape[0] != dR_dp.shape[0]:
        raise ValueError(
            f"K {K.shape} and dR/dp {dR_dp.shape} describe different systems")

    Kff, rhs_f, free = apply_dirichlet(K, -dR_dp, constrained_dofs)
    constrained = np.setdiff1d(np.arange(K.shape[0]), free)
    if free.size == 0:
        return SensitivityResult(
            du_dp=np.zeros_like(dR_dp), dR_dp=dR_dp, free_dofs=free,
            constrained_dofs=constrained, residual_norm=0.0,
            notes={"note": "every degree of freedom is constrained"})

    solution = np.linalg.solve(Kff, rhs_f)
    du_dp = np.zeros_like(dR_dp)
    du_dp[free] = solution
    residual = float(np.linalg.norm(Kff @ solution - rhs_f))
    condition = float(np.linalg.cond(Kff)) if estimate_condition else None
    return SensitivityResult(
        du_dp=du_dp, dR_dp=dR_dp, free_dofs=free, constrained_dofs=constrained,
        residual_norm=residual, condition_estimate=condition,
        notes={"sign_convention": "dR_dp is +dR/dp; the solve applies the minus"})
