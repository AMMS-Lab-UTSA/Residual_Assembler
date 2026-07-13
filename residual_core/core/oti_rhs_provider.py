"""OTILib RHS provider — the production, arbitrary-order sensitivity backend.

This is the **main** provider (Dual1 is legacy, order-1 only). It runs the
residual-method order loop of Aristizabal et al. using genuine OTILib scalars:

    initialize u* = real converged u   (all imaginary coefficients zero)
    seed parameters a* = a_i + e_i     (all selected parameters, simultaneously)
    for p = 1..q:
        evaluate R* = R(u*, a*)         (through the SAME formulation backend)
        extract R^(p)                   (all order-p imaginary coefficients)
        solve  T U^(p) = -R^(p)
        inject U^(p) into u*            (before evaluating the next order)
    return all U^(p) via a SensitivityRHSResult

The order-by-order injection is essential: for p > 1 the residual's p-th order
coefficients depend on the already-solved lower-order coefficients carried in u*.
This is NOT hand-coded derivative extraction — the residual is genuinely evaluated
with OTI numbers and the coefficients are read from the result.

Works through the standard formulation-backend interface (element or global), so
it is not special-cased to any one backend; any *dual/OTI-safe* residual (plain
scalar arithmetic) is supported. NumPy-kernel backends (e.g. C3D8) are not
OTI-safe and are reported as such rather than silently mis-evaluated.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from .rhs_provider import SensitivityRHSProvider, _split_parameter
from .sensitivity_package import SensitivityRHSResult, AlgebraMetadata
from ..algebra.otilib_adapter import (OtiContext, otilib_available, otilib_status,
                                      OtiUnavailableError)


class OtiLibRHSProvider(SensitivityRHSProvider):
    name = "otilib"
    max_supported_order = 64          # limited only by OTILib truncation order

    def evaluate_rhs(self, problem, real_solution, tangent, parameters,
                     order: int = 1) -> SensitivityRHSResult:
        if not otilib_available():
            raise OtiUnavailableError(otilib_status()["error"])
        if tangent is None or getattr(tangent, "T", None) is None:
            raise RuntimeError("OtiLibRHSProvider needs a dense tangent T to run "
                               "the order loop (T U^(p) = -R^(p)).")

        params = list(parameters)
        m = len(params)
        q = int(order)
        ndof = problem.dof_manager.ndof
        U = np.asarray(real_solution, float)
        T = np.asarray(tangent.T, float)
        free = tangent.free_mask
        free = np.ones(ndof, bool) if free is None else np.asarray(free, bool)

        ctx = OtiContext(num_bases=m, order=q)
        parameter_map = {p: i + 1 for i, p in enumerate(params)}

        # a* : each selected parameter seeded along its own basis, simultaneously
        real_vals = self._real_parameter_values(problem, params)
        a_star = {p: ctx.seed(real_vals[p], parameter_map[p]) for p in params}

        # u* : real part = converged solution, imaginary coefficients zero
        u_star = [ctx.scalar(float(U[j])) for j in range(ndof)]

        result = SensitivityRHSResult(
            ndof=ndof, algebra=AlgebraMetadata("otilib", n_bases=m, truncation_order=q),
            parameter_map=parameter_map, max_order=q)
        diagnostics: Dict[str, Any] = {"provider": self.name, "orders": [],
                                       "hypercomplex_ready": True,
                                       "backend": otilib_status()}
        U_orders: Dict[int, np.ndarray] = {}

        Tff = T[np.ix_(free, free)]
        for p in range(1, q + 1):
            R_oti = self._assemble_oti_residual(problem, ctx, u_star, a_star, params,
                                                parameter_map)
            dirs = ctx.order_directions(p)
            Np = len(dirs)
            R_p = np.zeros((ndof, Np))
            for col, d in enumerate(dirs):
                for j in range(ndof):
                    R_p[j, col] = ctx.coeff(R_oti[j], d["exponents"])
            # solve T U^(p) = -R^(p) on the free partition
            U_p = np.zeros((ndof, Np))
            U_p[free, :] = np.linalg.solve(Tff, -R_p[free, :])
            # inject U^(p) into u* along the matching order-p directions
            for col, d in enumerate(dirs):
                for j in range(ndof):
                    u_star[j] = ctx.set_coeff(u_star[j], d["exponents"], U_p[j, col])
            result.set_residual_order(p, R_p)
            U_orders[p] = U_p
            diagnostics["orders"].append({"p": p, "n_directions": Np})

        result.diagnostics = diagnostics
        # stash solved U^(p) for convenience / cross-check (package can re-solve too)
        result._solved_U = U_orders          # type: ignore[attr-defined]
        return result

    # ------------------------------------------------------------------ #
    def _real_parameter_values(self, problem, params) -> Dict[str, float]:
        vals = {}
        for pname in params:
            matname, key = _split_parameter(pname)
            binding = problem.model.materials.get(matname)
            sec = binding if isinstance(binding, dict) else getattr(binding, "section", None)
            vals[pname] = float(sec[key]) if isinstance(sec, dict) and key in sec else 0.0
        return vals

    def _assemble_oti_residual(self, problem, ctx, u_star, a_star, params, pmap):
        """Evaluate the global residual with OTI scalars: R*[j] as OTI numbers.

        One evaluation carries ALL seeded parameters simultaneously (the OTI
        multivariate advantage), so every order-p direction is extractable from a
        single residual per order."""
        model = problem.model
        forms = problem.forms
        dm = problem.dof_manager
        R = [ctx.scalar(0.0) for _ in range(dm.ndof)]
        for eid, el in model.elements.items():
            fk = model.element_formulation.get(eid)
            if fk is None:
                continue
            form = forms.get(fk)
            if form is None:
                continue
            edofs = np.array(dm.element_dofs(el.connectivity, form.dof_types), int)
            coords = model.coords_of(el.connectivity)
            u_e = [u_star[d] for d in edofs]
            mname = model.element_material.get(eid)
            props = self._oti_properties(model, mname, params, pmap, a_star, ctx)
            r_e, _k, _s, _d = form.eval_element(
                eid, el.etype, coords, u_e, {}, None, props,
                (0.0, 0.0), 0.0, None, {"compute_tangent": False})
            for j, val in enumerate(r_e):
                R[edofs[j]] = R[edofs[j]] + val
        return R

    @staticmethod
    def _oti_properties(model, mname, params, pmap, a_star, ctx):
        """Element properties dict with the targeted parameters replaced by their
        seeded OTI scalars (others left real)."""
        binding = model.materials.get(mname)
        sec = binding if isinstance(binding, dict) else getattr(binding, "section", None)
        if not isinstance(sec, dict):
            return binding
        out = dict(sec)
        for pname in params:
            mat, key = _split_parameter(pname)
            if mat == mname and key in out:
                out[key] = a_star[pname]
        return out


def default_oti_provider() -> OtiLibRHSProvider:
    return OtiLibRHSProvider()
