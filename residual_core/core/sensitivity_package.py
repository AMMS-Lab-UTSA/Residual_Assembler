"""Sensitivity right-hand-sides for the HYPAD residual method.

The endgame: emit the ingredients of

    T U^(p) = -R^(p)                                        (Aristizabal et al.)

for arbitrary derivative order ``p``. ``T`` (the tangent at the converged real
solution) and the order-0 residual come from ``core/results.py``; this module
defines the **p-th order right-hand-sides** ``R^(p)`` and ``rhs^(p) = -R^(p)``,
together with the algebra metadata (OTI basis count, truncation order), the
parameter→imaginary-direction map, and the per-order derivative-direction map.

No OTI/HYPAD algebra is implemented here. The container is a *contract*: it fixes
the shapes, orderings and maps so a hypercomplex residual evaluation can drop its
``R^(p)`` columns in later via ``set_residual_order(p, array)`` and everything
downstream (rhs, solve, export) already works. Column ``j`` of ``R^(p)`` is the
p-th order imaginary coefficient of the assembled residual along the j-th
order-p imaginary direction (see ``direction_map``).
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from itertools import combinations_with_replacement
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .results import _Serializable, _write_bundle, ResidualAssemblyResult, \
    StateResult, ValidationReport, TangentSource


# --------------------------------------------------------------------------- #
# algebra metadata + direction enumeration
# --------------------------------------------------------------------------- #
@dataclass
class AlgebraMetadata:
    """Hypercomplex algebra descriptor (OTI by default). Physics-free."""
    algebra: str = "OTI"               # OTI | dual | multidual | complex | ...
    n_bases: int = 0                   # m: number of imaginary bases (design params perturbed)
    truncation_order: int = 0          # n_t: maximum derivative order representable

    def n_coefficients(self) -> int:
        """Total OTI coefficients N = C(m + n_t, m) (Eq. 19)."""
        m, nt = self.n_bases, self.truncation_order
        if m == 0:
            return 1
        return math.comb(m + nt, m)

    def n_directions(self, p: int) -> int:
        """Number of order-p imaginary directions N^(p) = C(p + m - 1, p) (Eq. 20)."""
        m = self.n_bases
        if m == 0 or p == 0:
            return 1 if p == 0 else 0
        return math.comb(p + m - 1, p)

    def to_dict(self) -> Dict[str, Any]:
        return {"algebra": self.algebra, "n_bases": self.n_bases,
                "truncation_order": self.truncation_order,
                "n_coefficients": self.n_coefficients()}


def order_p_directions(m: int, p: int) -> List[Tuple[Tuple[int, ...], str, int]]:
    """Enumerate the order-p imaginary directions for ``m`` bases.

    Returns a list of ``(exponents, label, recovery_factor)`` where:
      - ``exponents`` is a length-m tuple of basis exponents (sum == p),
      - ``label`` is a human tag like ``e1^2`` or ``e1*e2``,
      - ``recovery_factor`` = prod(kappa_i!) — multiply the extracted imaginary
        coefficient by this to recover the true partial derivative (Eq. 23-28).
    """
    out: List[Tuple[Tuple[int, ...], str, int]] = []
    if m <= 0 or p <= 0:
        return out
    for combo in combinations_with_replacement(range(1, m + 1), p):
        exps = [0] * m
        for b in combo:
            exps[b - 1] += 1
        parts = []
        for i, k in enumerate(exps, start=1):
            if k == 1:
                parts.append("e%d" % i)
            elif k > 1:
                parts.append("e%d^%d" % (i, k))
        factor = 1
        for k in exps:
            factor *= math.factorial(k)
        out.append((tuple(exps), "*".join(parts), factor))
    return out


# --------------------------------------------------------------------------- #
# 4. Sensitivity RHS output
# --------------------------------------------------------------------------- #
@dataclass
class SensitivityRHSResult(_Serializable):
    ndof: int
    algebra: AlgebraMetadata
    parameter_map: Dict[str, int] = field(default_factory=dict)   # param name -> basis index (1-based)
    max_order: int = 1
    _R: Dict[int, Optional[np.ndarray]] = field(default_factory=dict)  # p -> (ndof, N^(p)) or None
    diagnostics: Dict[str, Any] = field(default_factory=dict)      # provider notes (source, readiness)

    def __post_init__(self):
        for p in range(1, self.max_order + 1):
            self._R.setdefault(p, None)

    # ---- direction / parameter maps -------------------------------------
    def direction_map(self, p: int) -> List[Dict[str, Any]]:
        dirs = order_p_directions(self.algebra.n_bases, p)
        inv = {v: k for k, v in self.parameter_map.items()}
        out = []
        for col, (exps, label, factor) in enumerate(dirs):
            params = []
            for i, k in enumerate(exps, start=1):
                if k:
                    params.extend([inv.get(i, "e%d" % i)] * k)
            out.append({"column": col, "label": label, "exponents": list(exps),
                        "parameters": params, "recovery_factor": factor})
        return out

    def expected_shape(self, p: int) -> Tuple[int, int]:
        return (self.ndof, self.algebra.n_directions(p))

    # ---- HYPAD plug-in point --------------------------------------------
    def set_residual_order(self, p: int, R_p: np.ndarray) -> "SensitivityRHSResult":
        """Attach the assembled order-p residual matrix R^(p) (from a HYPAD
        residual evaluation). Shape must be (ndof, N^(p))."""
        R_p = np.asarray(R_p, float)
        exp = self.expected_shape(p)
        if R_p.shape != exp:
            raise ValueError("R^(%d) has shape %s, expected %s"
                             % (p, R_p.shape, exp))
        self._R[p] = R_p
        self.max_order = max(self.max_order, p)
        return self

    def R(self, p: int) -> Optional[np.ndarray]:
        return self._R.get(p)

    def rhs(self, p: int) -> Optional[np.ndarray]:
        """rhs^(p) = -R^(p) (Eq. 75). None until R^(p) is supplied."""
        R_p = self._R.get(p)
        return None if R_p is None else -R_p

    def available_orders(self) -> List[int]:
        return sorted(p for p, v in self._R.items() if v is not None)

    # ---- serialization --------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        orders = {}
        for p in sorted(self._R):
            orders[str(p)] = {
                "expected_shape": list(self.expected_shape(p)),
                "n_directions": self.algebra.n_directions(p),
                "provided": self._R[p] is not None,
                "direction_map": self.direction_map(p),
            }
        return {
            "kind": "SensitivityRHSResult",
            "target_system": "T U^(p) = -R^(p)",
            "ndof": self.ndof,
            "algebra": self.algebra.to_dict(),
            "parameter_map": dict(self.parameter_map),
            "max_order": self.max_order,
            "available_orders": self.available_orders(),
            "orders": orders,
            "provider_diagnostics": self.diagnostics,
            "note": ("R^(p) columns are the p-th order imaginary coefficients of "
                     "the hypercomplex residual; supply them via "
                     "set_residual_order(p, array). rhs^(p) = -R^(p)."),
        }

    def arrays(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        for p, R_p in self._R.items():
            if R_p is not None:
                out["R_order_%d" % p] = R_p
                out["rhs_order_%d" % p] = -R_p
        return out

    def to_markdown(self) -> str:
        d = self.to_dict()
        L = ["# Sensitivity RHS  (T U^(p) = -R^(p))", "",
             "- ndof: %d" % self.ndof,
             "- algebra: %s  (bases m=%d, truncation n_t=%d, N=%d coeffs)"
             % (self.algebra.algebra, self.algebra.n_bases,
                self.algebra.truncation_order, self.algebra.n_coefficients()),
             "- parameter → basis: %s" % json.dumps(self.parameter_map),
             "- available orders (R^(p) supplied): %s" % (self.available_orders() or "none yet"),
             "", "## Per-order structure"]
        for p in sorted(self._R):
            od = d["orders"][str(p)]
            L.append("- **order %d**: R^(%d) shape %s, %d directions, provided=%s"
                     % (p, p, od["expected_shape"], od["n_directions"], od["provided"]))
            labels = ", ".join(x["label"] for x in od["direction_map"])
            if labels:
                L.append("    - directions: %s" % labels)
        L += ["", "> Supply R^(p) from a hypercomplex (OTI) residual evaluation to "
              "populate the right-hand-sides. The tangent T comes from the "
              "ResidualAssemblyResult; then U^(p) solves T U^(p) = -R^(p)."]
        return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# top-level package — the framework's deliverable
# --------------------------------------------------------------------------- #
@dataclass
class SensitivityPackage(_Serializable):
    """Bundles the answer to: *can I assemble this residual, what is missing, and
    (if ready) what are R, T and the RHS matrices for sensitivities?*"""
    mode: str
    runnable: bool
    residual: Optional[ResidualAssemblyResult] = None
    sensitivity: Optional[SensitivityRHSResult] = None
    state: Optional[StateResult] = None
    validation: Optional[ValidationReport] = None
    inspection: Optional[Dict[str, Any]] = None
    minimum_missing: Optional[str] = None

    # ---- system readiness -----------------------------------------------
    def tangent_available(self) -> bool:
        return bool(self.residual and self.residual.tangent
                    and self.residual.tangent.available)

    def system_ready(self, p: int = 1) -> bool:
        return bool(self.runnable and self.tangent_available()
                    and self.sensitivity is not None
                    and self.sensitivity.R(p) is not None)

    def solve(self, p: int = 1) -> np.ndarray:
        """Solve T U^(p) = -R^(p) on the free partition (prescribed DOFs have
        zero design-sensitivity). Requires T and R^(p) to be available."""
        if not self.tangent_available():
            raise RuntimeError("tangent T is unavailable (%s); cannot solve "
                               "T U^(p) = -R^(p)."
                               % (self.residual.tangent.source if self.residual
                                  and self.residual.tangent else "no tangent"))
        R_p = self.sensitivity.R(p) if self.sensitivity else None
        if R_p is None:
            raise RuntimeError("R^(%d) not supplied — attach it via "
                               "sensitivity.set_residual_order(%d, array) from a "
                               "HYPAD residual evaluation." % (p, p))
        T = self.residual.tangent.T
        free = self.residual.free_mask
        rhs = -R_p
        U = np.zeros_like(rhs)
        try:
            if free is None:
                U[:, :] = np.linalg.solve(T, rhs)
            else:
                Tff = T[np.ix_(free, free)]
                U[free, :] = np.linalg.solve(Tff, rhs[free, :])
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "tangent is singular on the free partition (%s) — the model has "
                "rigid-body / zero-stiffness modes (e.g. missing boundary "
                "conditions, or a pin-jointed truss with unconstrained transverse "
                "DOFs). A converged, properly-constrained FE problem yields a "
                "non-singular T." % exc)
        return U

    # ---- serialization --------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return {
            "kind": "SensitivityPackage",
            "mode": self.mode,
            "runnable": self.runnable,
            "minimum_missing": self.minimum_missing,
            "tangent_available": self.tangent_available(),
            "system_ready_order_1": self.system_ready(1),
            "target_system": "T U^(p) = -R^(p)",
            "residual": self.residual.to_dict() if self.residual else None,
            "tangent": (self.residual.tangent.to_dict()
                        if self.residual and self.residual.tangent else None),
            "sensitivity": self.sensitivity.to_dict() if self.sensitivity else None,
            "state": self.state.to_dict() if self.state else None,
            "validation": self.validation.to_dict() if self.validation else None,
            "inspection": self.inspection,
        }

    def arrays(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self.residual:
            for k, v in self.residual.arrays().items():
                out["residual_" + k] = v
        if self.sensitivity:
            for k, v in self.sensitivity.arrays().items():
                out["sensitivity_" + k] = v
        if self.state:
            for k, v in self.state.arrays().items():
                out["state_" + k] = v
        return out

    def to_markdown(self) -> str:
        L = ["# Sensitivity Package", "",
             "Answers: *can I assemble this residual, what is missing, and what "
             "are R, T and the RHS matrices for sensitivities?*", "",
             "- mode: `%s`" % self.mode,
             "- **runnable**: %s" % self.runnable]
        if not self.runnable:
            L.append("- minimum missing input: %s" % self.minimum_missing)
        L += ["- tangent available: %s" % self.tangent_available(),
              "- system `T U^(p) = -R^(p)` ready (order 1): %s" % self.system_ready(1),
              ""]
        if self.residual:
            L += ["## R (real residual)",
                  "- ||R_free|| = %.6e, ||reactions|| = %.6e"
                  % (self.residual.free_residual_norm, self.residual.reaction_norm),
                  ""]
        if self.residual and self.residual.tangent:
            L += ["## T (tangent)",
                  "- source: %s, available: %s"
                  % (self.residual.tangent.source, self.residual.tangent.available),
                  ""]
        if self.sensitivity:
            L += ["## R^(p) / rhs^(p)",
                  "- algebra: %s (m=%d, n_t=%d)"
                  % (self.sensitivity.algebra.algebra,
                     self.sensitivity.algebra.n_bases,
                     self.sensitivity.algebra.truncation_order),
                  "- orders present: %s"
                  % (self.sensitivity.available_orders() or "none (contract only)"),
                  ""]
        if self.validation:
            L += ["## Validation", "overall: %s" % self.validation.overall_passed, ""]
        return "\n".join(L) + "\n"

    def save(self, prefix: str) -> Dict[str, str]:
        """Write the whole package: a top-level bundle plus per-part files."""
        written = _write_bundle(prefix, self.to_dict(), self.arrays(), self.to_markdown())
        if self.residual:
            written["residual"] = self.residual.save(prefix + "_residual")
        if self.sensitivity:
            written["sensitivity"] = self.sensitivity.save(prefix + "_sensitivity")
        if self.validation:
            written["validation"] = self.validation.save(prefix + "_validation")
        if self.state:
            written["state"] = self.state.save(prefix + "_state")
        return written
