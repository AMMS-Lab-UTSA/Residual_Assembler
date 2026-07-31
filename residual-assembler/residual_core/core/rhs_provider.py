"""Sensitivity RHS providers — generate R^(p) from a residual evaluation.

The ``SensitivityPackage`` can *store* R^(p) and solve ``T U^(p) = -R^(p)``; a
provider is what actually *produces* R^(p). This module defines the provider
interface and one **legacy** implementation that uses first-order dual numbers
(``algebra/dual1.py``) to evaluate the residual with a perturbed design parameter
and read the derivative from the imaginary part — no hand-coded derivatives.

    R^(1)[:, i] = Im_i[ r(u, a + eps_i) ] = dR/da_i

which is exactly the right-hand-side of the first-order sensitivity system
(Aristizabal et al.): ``T (du/da_i) = -dR/da_i``.

Scope of this Dual1 provider: **order 1 only**, and only for **dual-safe**
formulation backends (those whose ``eval_element`` residual is written in plain
arithmetic, e.g. ``nonlinear_spring1``). It is a preliminary smoke test, not the
production path. The production backend is arbitrary-order OTILib —
``core/oti_rhs_provider.py`` with ``algebra/otilib_adapter.py`` — which seeds all
parameters simultaneously and solves order-by-order, injecting lower-order U into
U* before advancing. Backends whose residual uses NumPy kernels that reject
``Dual1`` are reported as not hypercomplex-ready — the contract degrades cleanly
rather than fabricating a derivative.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

import numpy as np

from ..algebra.dual1 import seed, imag_part
from .sensitivity_package import SensitivityRHSResult, AlgebraMetadata


class SensitivityRHSProvider(ABC):
    """Produce a populated :class:`SensitivityRHSResult` for a problem."""

    name = "abstract"
    max_supported_order = 0

    @abstractmethod
    def evaluate_rhs(self, problem, real_solution, tangent, parameters,
                     order: int = 1) -> SensitivityRHSResult:
        raise NotImplementedError


def _split_parameter(name: str):
    """'material.key' -> ('material', 'key'). A bare key maps to (None, key)."""
    if "." in name:
        mat, key = name.rsplit(".", 1)
        return mat, key
    return None, name


class DualNumberRHSProvider(SensitivityRHSProvider):
    """First-order RHS via dual numbers (proves the hypercomplex pathway).

    Perturbs one design parameter at a time along a dual imaginary direction,
    evaluates every element residual at the real converged solution, and scatters
    the imaginary parts into the corresponding column of R^(1).
    """

    name = "dual1"
    max_supported_order = 1

    def evaluate_rhs(self, problem, real_solution, tangent, parameters,
                     order: int = 1) -> SensitivityRHSResult:
        if order != 1:
            raise NotImplementedError(
                "DualNumberRHSProvider supports order 1 only (got %r). A full OTI "
                "backend is required for higher orders." % order)

        dm = problem.dof_manager
        model = problem.model
        forms = problem.forms
        params = list(parameters)
        m = len(params)
        ndof = dm.ndof
        U = np.asarray(real_solution, float)

        R1 = np.zeros((ndof, m))
        diagnostics: Dict[str, Any] = {"provider": self.name,
                                       "hypercomplex_ready": True, "notes": []}

        for col, pname in enumerate(params):
            matname, key = _split_parameter(pname)
            self._fill_column(model, forms, dm, U, R1, col, matname, key,
                              diagnostics)

        algebra = AlgebraMetadata(algebra="dual1", n_bases=m, truncation_order=1)
        parameter_map = {name: i + 1 for i, name in enumerate(params)}
        rhs = SensitivityRHSResult(ndof=ndof, algebra=algebra,
                                   parameter_map=parameter_map, max_order=1)
        rhs.set_residual_order(1, R1)
        rhs.diagnostics = diagnostics
        return rhs

    # ------------------------------------------------------------------ #
    def _fill_column(self, model, forms, dm, U, R1, col, matname, key, diag):
        for eid, el in model.elements.items():
            fk = model.element_formulation.get(eid)
            if fk is None:
                continue
            form = forms.get(fk)
            if form is None:
                continue
            edofs = np.array(dm.element_dofs(el.connectivity, form.dof_types), int)
            coords = model.coords_of(el.connectivity)
            u_e = U[edofs]
            mname = model.element_material.get(eid)
            props = self._perturbed_properties(model, mname, matname, key)
            try:
                r_e, _k, _s, _d = form.eval_element(
                    eid, el.etype, coords, u_e, {}, None, props,
                    (0.0, 0.0), 0.0, None, {"compute_tangent": False})
            except (TypeError, ValueError) as exc:
                # backend not dual-safe: leave column zero, record honestly
                diag["hypercomplex_ready"] = False
                diag["notes"].append(
                    "element %r (%s) not hypercomplex-ready: %s"
                    % (eid, fk, exc))
                continue
            for j, val in enumerate(r_e):
                R1[edofs[j], col] += imag_part(val)

    @staticmethod
    def _perturbed_properties(model, mname, matname, key):
        """Return the element's properties with ``key`` dual-seeded iff this
        element uses the targeted material and carries that parameter."""
        binding = model.materials.get(mname)
        sec = binding if isinstance(binding, dict) else getattr(binding, "section", None)
        if mname == matname and isinstance(sec, dict) and key in sec:
            dual = dict(sec)
            dual[key] = seed(sec[key])
            return dual
        return binding


def default_rhs_provider() -> SensitivityRHSProvider:
    return DualNumberRHSProvider()
