"""Partner-side residual provider contract.

A partner implements ONE of these locally. The kit never sees the mesh, the
material model, or the source — it only calls the declared methods. The residual
functions must be **scalar-generic**: a ``parameters[name]`` value may be a plain
float OR a hypercomplex seed; write the residual with ordinary arithmetic and it
carries the derivative(s). The production backend is OTILib (arbitrary order,
multi-parameter); ``hypercomplex.Dual1`` is only a legacy first-order smoke test.

Three provider levels (implement one):

  ElementResidualProvider  — element-by-element residual (+ optional tangent/state)
  GlobalResidualProvider   — whole-model residual vector
  BlackBoxResidualProvider — an executable behind a JSON contract (no linking)

Tangent is optional and can be supplied three ways:
  - ``get_tangent(...)``            dense T
  - ``apply_tangent(..., x)``       matrix-free T·x
  - an external ``tangent.npz`` file passed to the CLI

Nothing here requires exposing implementation details. Metadata is limited to a
provider name, the ordered parameter names, and the DOF count.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np


class ResidualProvider(ABC):
    """Base metadata every provider declares (no proprietary details)."""

    name: str = "provider"
    parameters: Sequence[str] = ()      # ordered design-parameter names
    ndof: int = 0                       # number of global DOFs

    def parameter_values(self) -> Dict[str, float]:
        """Current (nominal) parameter values. Override to supply real values."""
        return {p: 0.0 for p in self.parameters}

    def free_mask(self) -> Optional[np.ndarray]:
        """Boolean mask of free DOFs (True=free). None => all DOFs free.
        Prescribed (Dirichlet) DOFs carry zero design sensitivity."""
        return None

    # optional tangent hooks (return None if not provided) -----------------
    def get_tangent(self, U, state, parameters, time, dtime) -> Optional[np.ndarray]:
        return None

    def apply_tangent(self, U, state, parameters, x, time, dtime) -> Optional[np.ndarray]:
        return None

    def kind(self) -> str:
        return type(self).__name__


class GlobalResidualProvider(ResidualProvider):
    """Partner returns the whole global residual vector."""

    @abstractmethod
    def eval_global_residual(self, global_solution, state, parameters, time, dtime):
        """Return the residual as a length-``ndof`` sequence. Scalar-generic:
        entries may be floats or ``Dual1`` depending on ``parameters``."""
        raise NotImplementedError


class ElementResidualProvider(ResidualProvider):
    """Partner returns element residuals; the kit scatters them."""

    @abstractmethod
    def elements(self) -> Iterable[Any]:
        """Iterable of element ids."""
        raise NotImplementedError

    @abstractmethod
    def element_dof_map(self, element_id) -> List[int]:
        """Global DOF indices touched by this element (scatter map)."""
        raise NotImplementedError

    def element_state(self, element_id):
        """Previous state for this element (or None)."""
        return None

    @abstractmethod
    def eval_element_residual(self, element_id, element_dofs, element_state,
                              parameters, time, dtime):
        """Return ``(r_e, tangent_optional, state_optional)``. Scalar-generic:
        ``r_e`` entries may be floats or ``Dual1``."""
        raise NotImplementedError


class BlackBoxResidualProvider(ResidualProvider):
    """Partner runs an executable behind a JSON request/response contract.

    The kit writes a request (solution, parameter values, seeded directions,
    order) and the executable returns the residual coefficients R^(p) (and
    optionally the tangent). See ``blackbox_runner.py`` and
    ``docs/residual_provider_contract.md``. The partner's code stays entirely
    behind the executable boundary; the kit performs no scalar overloading.
    """

    def __init__(self, runner, ndof: int, parameters: Sequence[str],
                 free_mask: Optional[np.ndarray] = None, name: str = "blackbox",
                 parameter_values: Optional[Dict[str, float]] = None):
        self._runner = runner
        self.ndof = int(ndof)
        self.parameters = tuple(parameters)
        self._free = free_mask
        self.name = name
        self._values = dict(parameter_values or {})

    def parameter_values(self):
        vals = {p: 0.0 for p in self.parameters}
        vals.update(self._values)
        return vals

    def free_mask(self):
        return self._free

    def eval_rhs(self, order, global_solution, parameter_values,
                 seed_directions, time, dtime) -> Dict[str, Any]:
        """Ask the executable for R^(p) (and optionally T). Returns the parsed
        response dict with at least ``residual_coefficients`` (ndof x m).

        The request also carries the OTI configuration (basis count, truncation
        order, direction map) so a partner whose executable handles OTI internally
        can return higher-order coefficients."""
        m = len(seed_directions)
        direction_map = {str(p): _order_exponents(m, p) for p in range(1, int(order) + 1)}
        request = {
            "schema": "resasm-partner-request/1",
            "order": int(order),
            "solution": list(map(float, global_solution)),
            "parameters": {k: float(v) for k, v in parameter_values.items()},
            "seed_directions": seed_directions,     # {param: basis_index}
            "basis_count": m,                        # OTI m
            "truncation_order": int(order),          # OTI nt
            "direction_map": direction_map,          # order -> list of exponent vectors
            "time": list(time), "dtime": float(dtime),
        }
        return self._runner.run(request)


def _order_exponents(m, p):
    """Canonical order-p exponent multi-indices for m bases (list of lists)."""
    from itertools import combinations_with_replacement
    out = []
    if m <= 0 or p <= 0:
        return out
    for combo in combinations_with_replacement(range(1, m + 1), p):
        exps = [0] * m
        for b in combo:
            exps[b - 1] += 1
        out.append(exps)
    return out
