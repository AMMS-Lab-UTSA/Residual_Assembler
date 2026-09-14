"""dR/dq: the residual's derivative with respect to the internal state.

The residual of a history-dependent material is a function of three things,
not two::

    R(u, p, q) = integral_V B^T sigma(eps(u), p, q_n) dV - f_ext

``q_n`` is the state the material carried INTO the increment -- accumulated
plastic strain, a damage variable, a viscous internal strain. It is an
independent argument of the residual: freeze ``u`` and ``p``, move ``q_n``,
and the residual moves. Three things need that derivative and none of them can
get it from ``dR/du``: a time-continuous adjoint, which differentiates the
whole path and therefore needs the coupling between one increment's state and
the next; a sensitivity that treats an initial state as an unknown; and any
attribution of a wrong residual to the state rather than to the strain.

The assembly is the same integral as everything else in this file's
neighbourhood, with a different thing in the middle::

    dR_e/dq_(k,j) = w_k detJ_k B_k^T (d sigma_k / d q_(k,j))

with one structural fact that makes it different from ``dR/dp``: a parameter
is shared by every integration point and a STATE IS NOT. ``q_(k,j)`` is the
``j``-th state variable at integration point ``k``, and it moves the stress at
that point only, so ``dR/dq`` is block-sparse by construction -- one block of
``ndof_e`` rows per integration point. Assembling it as though the state were
global would produce a dense matrix whose extra entries are all wrong and
whose diagonal blocks are right, which is exactly the kind of error a
single-column check does not see.

Column ordering is stated once and never inferred: integration-point major,
state-variable minor, ``column = k * n_state + j``, matching the way the
StateManager lays out a ``(n_ip, n_state)`` slab in memory.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np

from residual_core.formulations.c3d8_kernel import (
    ABAQUS_C3D8_GAUSS, b_matrix_reference, b_matrix_spatial)

__all__ = ["element_dR_dq", "assemble_dR_dq", "StateColumns",
           "state_column", "StateSensitivityError"]


class StateSensitivityError(ValueError):
    """A state derivative that cannot be assembled, said in terms of why."""


def state_column(ip: int, state: int, n_state: int) -> int:
    """Which column of ``dR/dq`` holds state ``state`` at point ``ip``."""
    return int(ip) * int(n_state) + int(state)


@dataclass
class StateColumns:
    """What each column of an assembled ``dR/dq`` is, so nobody has to guess."""

    n_ip: int
    n_state: int
    elements: tuple = ()
    order: str = "integration-point major, state-variable minor"

    def __len__(self) -> int:
        return self.n_ip * self.n_state * max(len(self.elements), 1)

    def column_of(self, element_id, ip: int, state: int) -> int:
        if not self.elements:
            return state_column(ip, state, self.n_state)
        try:
            offset = list(self.elements).index(element_id)
        except ValueError:
            raise StateSensitivityError(
                f"element {element_id!r} has no state columns in this "
                f"assembly; it holds {list(self.elements)}") from None
        return (offset * self.n_ip * self.n_state
                + state_column(ip, state, self.n_state))

    def describe(self, column: int) -> dict:
        block = self.n_ip * self.n_state
        element = (list(self.elements)[column // block] if self.elements
                   else None)
        within = column % block
        return {"element": element, "ip": within // self.n_state,
                "state": within % self.n_state}


def element_dR_dq(Xe, dsigma_dq_at_ip, Ue=None, *, mode: str = "small",
                  gauss=ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """One element's ``dR/dq``: ``(ndof_e, n_ip * n_state)``.

    ``dsigma_dq_at_ip`` is ``(n_ip, ntens, n_state)`` -- at each integration
    point, how that point's stress moves when that point's own state moves.
    The block structure is enforced here rather than left to the caller: a
    state at point ``k`` writes only into the columns of point ``k``.
    """
    Xe = np.asarray(Xe, dtype=float)
    dsigma = np.asarray(dsigma_dq_at_ip, dtype=float)
    n_ip = len(gauss.weights)
    if dsigma.ndim != 3 or dsigma.shape[0] != n_ip:
        raise StateSensitivityError(
            f"dsigma_dq_at_ip must be (n_ip={n_ip}, ntens, n_state); got "
            f"{dsigma.shape}. A state derivative shared across integration "
            f"points is not a state derivative.")
    n_state = dsigma.shape[2]
    if mode == "finite":
        if Ue is None:
            raise StateSensitivityError(
                "finite-strain mode needs the element displacement")
        xe = Xe + np.asarray(Ue, dtype=float).reshape(-1, 3)
    ndof_e = 3 * Xe.shape[0]
    out = np.zeros((ndof_e, n_ip * n_state), dtype=float)
    for k, weight in enumerate(gauss.weights):
        if mode == "finite":
            B, detJ = b_matrix_spatial(xe, gauss.points[k])
        else:
            B, detJ = b_matrix_reference(Xe, gauss.points[k])
        block = (B.T @ dsigma[k]) * detJ * weight        # (ndof_e, n_state)
        for j in range(n_state):
            out[:, state_column(k, j, n_state)] = block[:, j]
    return out


def assemble_dR_dq(node_ids: Sequence[int], coords, connectivity,
                   dsigma_dq_all: Mapping[int, np.ndarray], U=None, *,
                   mode: str = "small", gauss=ABAQUS_C3D8_GAUSS):
    """Scatter every element's ``dR/dq`` into a global ``(ndof, ncols)``.

    The global DOF ordering is the same one ``assemble_global_internal_force``
    and ``assemble_dR_dp`` use -- three displacements per node, nodes ascending
    by id -- so the assembled residual, its displacement derivative and its
    state derivative are indexed identically and can be used in one equation.

    Returns ``(dR_dq, node_index, columns)``; ``columns`` is the
    :class:`StateColumns` that says what every column is.
    """
    sorted_ids = sorted(node_ids)
    node_index = {nid: i for i, nid in enumerate(sorted_ids)}
    coords = np.asarray(coords, dtype=float)
    id_to_row = {nid: i for i, nid in enumerate(list(node_ids))}
    elements = [eid for eid, _conn in connectivity]
    first = np.asarray(dsigma_dq_all[elements[0]], dtype=float)
    n_ip, n_state = first.shape[0], first.shape[2]
    columns = StateColumns(n_ip=n_ip, n_state=n_state, elements=tuple(elements))
    out = np.zeros((3 * len(sorted_ids), n_ip * n_state * len(elements)),
                   dtype=float)

    is_dict_U = isinstance(U, dict)
    U_arr = None if (is_dict_U or U is None) else np.asarray(U, dtype=float)

    for offset, (eid, conn) in enumerate(connectivity):
        Xe = np.array([coords[id_to_row[nid]] for nid in conn], dtype=float)
        Ue = None
        if mode == "finite":
            if is_dict_U:
                Ue = np.array([U[nid] for nid in conn], dtype=float)
            elif U_arr is not None:
                Ue = np.array([U_arr[3 * node_index[nid]:3 * node_index[nid] + 3]
                               for nid in conn], dtype=float)
        block = element_dR_dq(Xe, dsigma_dq_all[eid], Ue, mode=mode, gauss=gauss)
        if block.shape[1] != n_ip * n_state:
            raise StateSensitivityError(
                f"element {eid} has {block.shape[1]} state columns and the "
                f"first element has {n_ip * n_state}; a global state vector "
                f"cannot be indexed if elements disagree about its length")
        base = offset * n_ip * n_state
        for local, nid in enumerate(conn):
            gi = node_index[nid]
            out[3 * gi:3 * gi + 3, base:base + n_ip * n_state] += \
                block[3 * local:3 * local + 3, :]
    return out, node_index, columns
