"""Formulation-agnostic global assembler.

The assembler knows only how to (1) look up the formulation bound to each
element, (2) gather that element's coordinates, DOFs, material state and
properties, (3) call the formulation's ``eval_element`` and (4) scatter the
returned residual/tangent into the global system, then add external loads. It
contains NO element mathematics and NO material knowledge -- swap in any
Formulation backend and it just works.

    R(U) = F_internal(U, state) - F_external          (constraints handled by
           |__ sum of element eval_element residuals   core/constraints.py split)
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from . import loads
from . import constraints as _constraints


class Assembler:
    def __init__(self, model, dof_manager, formulations: Dict[str, object],
                 state_manager=None):
        self.model = model
        self.dm = dof_manager
        self.formulations = formulations          # name -> Formulation instance
        self.state = state_manager

    def assemble(self, U, time=(0.0, 0.0), dtime=0.0, fields=None, options=None,
                 compute_tangent=False, load_factor=1.0, solution_state=None,
                 collect_elements=False):
        """Return (R (ndof,), K (ndof,ndof) or None, diagnostics).

        If ``collect_elements`` is True, ``diagnostics['element_residuals']`` maps
        each assembled element id to ``(edofs, r_e)`` — the scatter map and the
        element's residual contribution (used by the output/sensitivity package).
        """
        ndof = self.dm.ndof
        U = np.asarray(U, dtype=float)
        R = np.zeros(ndof)
        K = np.zeros((ndof, ndof)) if compute_tangent else None
        opts = dict(options or {})
        opts.setdefault("compute_tangent", compute_tangent)
        base_ss = dict(solution_state or {})
        U_prev = np.asarray(base_ss["U_prev"], float) if "U_prev" in base_ss else None

        diag = {"elements": 0, "skipped_no_formulation": 0, "formulations": {},
                "tangent_contributions": 0}
        if collect_elements:
            diag["element_residuals"] = {}
        for eid, el in self.model.elements.items():
            fk = self.model.element_formulation.get(eid)
            if fk is None:
                diag["skipped_no_formulation"] += 1
                continue
            form = self.formulations.get(fk)
            if form is None:
                raise KeyError("element %r bound to unregistered formulation %r"
                               % (eid, fk))
            # Ask the DOF manager for exactly this formulation's DOFs, in the
            # formulation's own node-major order. In a uniform displacement model
            # this reduces to the original all-node-DOFs gather; in a mixed model
            # it selects the right subset from each node's heterogeneous DOF set.
            edofs = np.array(
                self.dm.element_dofs(el.connectivity,
                                     getattr(form, "dof_types", None)),
                dtype=int)
            coords = self.model.coords_of(el.connectivity)
            u_e = U[edofs]
            ss = dict(base_ss)
            if U_prev is not None:
                ss["dofs_prev"] = U_prev[edofs]
            mat_state = self.state.get(eid) if self.state is not None else None
            mname = self.model.element_material.get(eid)
            props = self.model.materials.get(mname)

            r_e, k_e, s_new, edic = form.eval_element(
                eid, el.etype, coords, u_e, ss, mat_state, props,
                time, dtime, fields, opts)

            np.add.at(R, edofs, np.asarray(r_e, float))
            if compute_tangent and k_e is not None:
                K[np.ix_(edofs, edofs)] += np.asarray(k_e, float)
                diag["tangent_contributions"] += 1
            if self.state is not None and s_new is not None:
                self.state.set_trial(eid, s_new)
            if collect_elements:
                diag["element_residuals"][eid] = (edofs.copy(),
                                                  np.asarray(r_e, float).copy())
            diag["elements"] += 1
            diag["formulations"][fk] = diag["formulations"].get(fk, 0) + 1

        F_ext = loads.external_force(self.model, self.dm, load_factor)
        R = R - F_ext
        return R, K, diag

    # ------------------------------------------------------------------ #
    def split(self, R):
        """Split a global vector into (free_residual, prescribed_idx,
        reaction_at_prescribed) using the model's Dirichlet BCs."""
        free_mask, pres_idx, _ = _constraints.partition(self.model, self.dm)
        return R[free_mask], pres_idx, (R[pres_idx] if pres_idx.size else np.zeros(0))
