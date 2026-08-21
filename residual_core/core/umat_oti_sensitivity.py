"""Direct dR/dp assembly from UMAT-OTI point-wise sensitivities.

Given per-integration-point sensitivities of the Cauchy stress with respect
to material parameters (``DSIGMA_DP``), the internal-force residual
sensitivity at each parameter ``p_k`` is::

    dR_e/dp_k = integral_V ( B^T · dsigma/dp_k ) dV

Discretising with quadrature reduces this to a sum over integration points
of ``w * detJ * B^T · dsigma/dp_k``. This module implements exactly that
sum for a *single* parameter column at a time, and a small wrapper that
loops over the columns and scatters into the global dR/dp array using a
supplied element-DOF map.

The routine is deliberately formulation-agnostic: the caller supplies the
element ``B`` matrices, weights, and DOF maps. This mirrors the way the
existing residual assembler is structured (see ``core/assembler.py``).
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

import numpy as np


def integrate_dRe_dp_single(
    B_at_ip: Sequence[np.ndarray],
    weight_x_detJ_at_ip: Sequence[float],
    dsigma_dp_at_ip: Sequence[np.ndarray],
) -> np.ndarray:
    """Assemble the *element* residual sensitivity dR_e / dp_k.

    Parameters
    ----------
    B_at_ip
        Sequence of ``ntens x ndof_e`` strain-displacement matrices, one per
        integration point in element order.
    weight_x_detJ_at_ip
        Sequence of scalar weights ``w * detJ`` for each integration point.
    dsigma_dp_at_ip
        Sequence of ``ntens`` vectors giving ``dsigma/dp_k`` at each
        integration point.

    Returns
    -------
    numpy.ndarray of shape ``(ndof_e,)`` holding ``sum_ip w*detJ * B^T dsigma/dp_k``.
    """
    if not B_at_ip:
        raise ValueError("B_at_ip must not be empty")
    if not (len(B_at_ip) == len(weight_x_detJ_at_ip) == len(dsigma_dp_at_ip)):
        raise ValueError(
            "B_at_ip, weight_x_detJ_at_ip and dsigma_dp_at_ip must have the "
            "same length"
        )
    ndof_e = B_at_ip[0].shape[1]
    accum = np.zeros(ndof_e, dtype=float)
    for B, w, ds in zip(B_at_ip, weight_x_detJ_at_ip, dsigma_dp_at_ip):
        if B.shape[1] != ndof_e:
            raise ValueError("B matrices must all have the same column count")
        if ds.shape[0] != B.shape[0]:
            raise ValueError(
                "dsigma_dp row count must match B row count (ntens)"
            )
        accum += float(w) * (B.T @ np.asarray(ds, dtype=float))
    return accum


def scatter_dRe_dp_to_global(
    dRe_dp: np.ndarray,
    element_dof_map: Sequence[int],
    global_dRdp: np.ndarray,
) -> None:
    """Scatter an element residual sensitivity into the global ``dR/dp`` array.

    ``element_dof_map`` gives, for each element DOF (in the same order as the
    columns of ``B``), the global DOF index. ``global_dRdp`` is modified
    in place.
    """
    if dRe_dp.shape[0] != len(element_dof_map):
        raise ValueError("element_dof_map length must match dRe_dp length")
    for local_dof, global_dof in enumerate(element_dof_map):
        global_dRdp[global_dof] += dRe_dp[local_dof]


def assemble_dRdp(
    *,
    n_global_dofs: int,
    element_contributions: Iterable[dict],
) -> np.ndarray:
    """Assemble the full ``dR/dp`` array from a list of element contributions.

    ``element_contributions`` yields one dict per element with keys:

    ``"B_at_ip"``
        Sequence of ``(ntens, ndof_e)`` numpy arrays.
    ``"weight_x_detJ_at_ip"``
        Sequence of quadrature weights times detJ.
    ``"dsigma_dp_at_ip"``
        Sequence of ``(ntens, nparam)`` matrices (per integration point).
    ``"element_dof_map"``
        Sequence of ``ndof_e`` global-DOF indices.

    Returns a ``(n_global_dofs, nparam)`` numpy array.
    """
    elements = list(element_contributions)
    if not elements:
        raise ValueError("no element contributions supplied")
    first_dsigma = elements[0]["dsigma_dp_at_ip"][0]
    nparam = np.asarray(first_dsigma).shape[1]
    global_dRdp = np.zeros((n_global_dofs, nparam), dtype=float)
    for element in elements:
        B_list = element["B_at_ip"]
        weights = element["weight_x_detJ_at_ip"]
        dsigma_at_ips = [np.asarray(m, dtype=float) for m in element["dsigma_dp_at_ip"]]
        dof_map = element["element_dof_map"]
        ndof_e = B_list[0].shape[1]
        element_dRdp = np.zeros((ndof_e, nparam), dtype=float)
        for k in range(nparam):
            dsigma_k = [ds[:, k] for ds in dsigma_at_ips]
            element_dRdp[:, k] = integrate_dRe_dp_single(B_list, weights, dsigma_k)
        for k in range(nparam):
            scatter_dRe_dp_to_global(element_dRdp[:, k], dof_map, global_dRdp[:, k])
    return global_dRdp


__all__ = [
    "assemble_dRdp",
    "integrate_dRe_dp_single",
    "scatter_dRe_dp_to_global",
]
