"""Global DOF numbering — formulation-agnostic, heterogeneous-DOF capable.

A node may carry different degrees of freedom depending on which formulation(s)
touch it (a beam node has rotations, a solid node does not, a thermal node has
temperature). This manager therefore supports **per-node DOF-type sets**, while
staying fully backward compatible with the original uniform displacement-only
numbering.

Two construction paths
----------------------
1. Uniform (original behaviour, unchanged)::

        DofManager(node_ids)                        # UX,UY,UZ on every node
        DofManager(node_ids, ("T",))                # single-field example

2. Heterogeneous, inferred from the formulations bound to each element::

        DofManager.for_model(model, formulations)

   Each node's DOF set is the ordered union of the ``dof_types`` declared by the
   formulations of the elements that reference it. Nodes referenced by no
   assembled element receive zero DOFs.

Global DOF layout: nodes are ordered ascending by id and each is assigned a
contiguous block equal to the number of DOF types it carries. The mapping from a
(node, dof-type) pair to a global index is dense and gap-free.

The core assembler stays agnostic: it asks the formulation for its ``dof_types``
and calls ``element_dofs(connectivity, dof_types)`` to gather exactly those DOFs
in the formulation's own node-major order.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple


class DofManager:
    def __init__(self, node_ids: Iterable[int],
                 dof_types: Tuple[str, ...] = ("UX", "UY", "UZ"),
                 node_dof_types: Optional[Dict[int, Sequence[str]]] = None):
        self.node_ids = sorted(node_ids)

        if node_dof_types is None:
            # ---- uniform mode (original behaviour) ----------------------
            self.uniform = True
            self.dof_types = tuple(dof_types)
            self.ndof_per_node = len(self.dof_types)
            self._node_types: Dict[int, Tuple[str, ...]] = {
                nid: self.dof_types for nid in self.node_ids}
        else:
            # ---- heterogeneous mode -------------------------------------
            self.uniform = False
            self.dof_types = tuple(dof_types)       # a "reference" default only
            self.ndof_per_node = None               # varies by node
            self._node_types = {
                nid: tuple(node_dof_types.get(nid, ()))
                for nid in self.node_ids}

        # dense global numbering: contiguous block per node
        self.node_index: Dict[int, int] = {}
        self._node_base: Dict[int, int] = {}
        self._type_offset: Dict[int, Dict[str, int]] = {}
        base = 0
        for i, nid in enumerate(self.node_ids):
            self.node_index[nid] = i
            self._node_base[nid] = base
            types = self._node_types[nid]
            self._type_offset[nid] = {t: j for j, t in enumerate(types)}
            base += len(types)
        self.ndof = base

    # ------------------------------------------------------------------ #
    @classmethod
    def for_model(cls, model, formulations: Dict[str, object]) -> "DofManager":
        """Build a heterogeneous manager from a model + formulation registry.

        Each node's DOF-type set is the ordered union of the ``dof_types`` of the
        formulations bound to the elements that reference it.
        """
        node_types: Dict[int, List[str]] = {}
        for eid, el in model.elements.items():
            fk = model.element_formulation.get(eid)
            if fk is None:
                continue
            form = formulations.get(fk)
            if form is None:
                continue
            for nid in el.connectivity:
                bucket = node_types.setdefault(nid, [])
                for dt in form.dof_types:
                    if dt not in bucket:
                        bucket.append(dt)
        return cls(model.nodes.keys(), node_dof_types=node_types)

    # ------------------------------------------------------------------ #
    def dof_types_of(self, nid: int) -> Tuple[str, ...]:
        return self._node_types[nid]

    def n_node_dofs(self, nid: int) -> int:
        return len(self._node_types[nid])

    def node_dofs(self, nid: int) -> List[int]:
        base = self._node_base[nid]
        return [base + j for j in range(len(self._node_types[nid]))]

    def element_dofs(self, connectivity: Iterable[int],
                     dof_types: Optional[Sequence[str]] = None) -> List[int]:
        """Global DOF indices for an element's connectivity, node-major.

        If ``dof_types`` is given, only those DOFs are gathered per node, in the
        requested order (a formulation asking for exactly its own DOFs). If None,
        every DOF carried by each node is returned (original uniform behaviour).
        """
        out: List[int] = []
        if dof_types is None:
            for nid in connectivity:
                out.extend(self.node_dofs(nid))
            return out

        for nid in connectivity:
            base = self._node_base[nid]
            offs = self._type_offset[nid]
            for t in dof_types:
                if t not in offs:
                    raise KeyError(
                        "node %r does not carry DOF type %r (it has %r); the "
                        "formulation's DOFs are not provided at this node"
                        % (nid, t, self._node_types[nid]))
                out.append(base + offs[t])
        return out

    def dof_index(self, nid: int, local_dof: int) -> int:
        """Global index of a 1-based local DOF (Abaqus convention: 1=x,2=y,3=z,
        4=rx,5=ry,6=rz for mechanics). Resolves against the node's own DOF set."""
        base = self._node_base[nid]
        n = len(self._node_types[nid])
        j = local_dof - 1
        if not (0 <= j < n):
            raise IndexError("node %r has %d dofs; local_dof %d out of range"
                             % (nid, n, local_dof))
        return base + j

    def has_dof(self, nid: int, local_dof: int) -> bool:
        n = len(self._node_types.get(nid, ()))
        return 1 <= local_dof <= n
