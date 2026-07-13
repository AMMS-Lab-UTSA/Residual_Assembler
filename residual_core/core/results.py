"""Serializable output objects — the framework's *product*.

The verification zoo is a proving ground; the deliverable is the set of numerical
ingredients a HYPAD (hypercomplex) residual-sensitivity solve needs:

    T U^(p) = -R^(p)

This module defines the result containers for the **real (order-0) side** of that
system — the assembled residual, the tangent, state/history, and validation — plus
a common serialization contract (``.json`` metadata + ``.npz`` arrays + ``.md``
report). The p-th order right-hand-sides live in ``sensitivity_package.py`` and
are designed so an OTI/HYPAD backend can plug in later without changing this API.

See ``docs/output_contract.md`` for the full contract.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Framework sign convention (documented once, referenced everywhere).
SIGN_CONVENTION = ("residual: R = F_internal - F_external (+F_constraints); "
                   "each element returns +F_internal. Abaqus reports RHS = -R "
                   "and AMATRX = dR/du.")


class TangentSource:
    """Where the tangent T came from (contract section 3)."""
    SOLVER_EXPORTED = "solver-exported"
    BACKEND_ASSEMBLED = "backend-assembled"
    FINITE_DIFFERENCE = "finite-difference"
    HYPAD = "hypad"
    UNAVAILABLE = "unavailable"


# --------------------------------------------------------------------------- #
# serialization helpers
# --------------------------------------------------------------------------- #
def dof_ordering(dof_manager) -> List[str]:
    """Global DOF labels in assembly order, e.g. ['1:UX','1:UY',...,'12:RZ']."""
    labels: List[str] = []
    for nid in dof_manager.node_ids:
        for t in dof_manager._node_types[nid]:
            labels.append("%s:%s" % (nid, t))
    return labels


def _write_bundle(prefix: str, meta: Dict[str, Any],
                  arrays: Dict[str, np.ndarray], markdown: str) -> Dict[str, str]:
    """Write ``prefix.json`` (metadata), ``prefix.npz`` (arrays, if any) and
    ``prefix.md`` (human report). Returns the written paths."""
    written = {}
    d = os.path.dirname(prefix)
    if d and not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)
    with open(prefix + ".json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
    written["json"] = prefix + ".json"
    real_arrays = {k: v for k, v in arrays.items() if v is not None}
    if real_arrays:
        np.savez(prefix + ".npz", **real_arrays)
        written["npz"] = prefix + ".npz"
    with open(prefix + ".md", "w", encoding="utf-8") as fh:
        fh.write(markdown)
    written["md"] = prefix + ".md"
    return written


class _Serializable:
    """Mixin: subclasses implement to_dict(), arrays(), to_markdown()."""

    def to_dict(self) -> Dict[str, Any]:
        raise NotImplementedError

    def arrays(self) -> Dict[str, np.ndarray]:
        return {}

    def to_markdown(self) -> str:
        return "# %s\n\n```json\n%s\n```\n" % (
            type(self).__name__, json.dumps(self.to_dict(), indent=2))

    def to_json(self, path: str) -> str:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)
        return path

    def to_npz(self, path: str) -> Optional[str]:
        a = {k: v for k, v in self.arrays().items() if v is not None}
        if not a:
            return None
        np.savez(path, **a)
        return path

    def save(self, prefix: str) -> Dict[str, str]:
        """Write prefix.json + prefix.npz + prefix.md. Returns written paths."""
        return _write_bundle(prefix, self.to_dict(), self.arrays(), self.to_markdown())


# --------------------------------------------------------------------------- #
# 3. Tangent output
# --------------------------------------------------------------------------- #
@dataclass
class TangentResult(_Serializable):
    T: Optional[np.ndarray] = None
    source: str = TangentSource.UNAVAILABLE
    dof_labels: List[str] = field(default_factory=list)
    free_mask: Optional[np.ndarray] = None
    prescribed_idx: Optional[np.ndarray] = None
    note: str = ""

    @property
    def available(self) -> bool:
        return self.T is not None

    @property
    def ndof(self) -> int:
        if self.T is not None:
            return int(self.T.shape[0])
        return len(self.dof_labels)

    @property
    def symmetric(self) -> Optional[bool]:
        if self.T is None:
            return None
        return bool(np.allclose(self.T, self.T.T, atol=1e-8))

    def to_dict(self) -> Dict[str, Any]:
        n_pres = int(self.prescribed_idx.size) if self.prescribed_idx is not None else 0
        return {
            "kind": "TangentResult",
            "available": self.available,
            "source": self.source,
            "ndof": self.ndof,
            "shape": list(self.T.shape) if self.T is not None else None,
            "symmetric": self.symmetric,
            "n_free": self.ndof - n_pres,
            "n_prescribed": n_pres,
            "dof_ordering": self.dof_labels,
            "note": self.note,
        }

    def arrays(self) -> Dict[str, np.ndarray]:
        return {"T": self.T, "free_mask": self.free_mask,
                "prescribed_idx": self.prescribed_idx}

    def to_markdown(self) -> str:
        d = self.to_dict()
        L = ["# Tangent (T)", "",
             "- available: %s" % d["available"],
             "- source: %s" % d["source"],
             "- shape: %s" % d["shape"],
             "- symmetric: %s" % d["symmetric"],
             "- free / prescribed DOFs: %d / %d" % (d["n_free"], d["n_prescribed"])]
        if self.note:
            L.append("- note: %s" % self.note)
        if not self.available:
            L.append("")
            L.append("> T is unavailable for this mode. Provide a solver-exported "
                     "tangent, use a backend that assembles one, or enable "
                     "finite-difference to complete `T U^(p) = -R^(p)`.")
        return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# 2. Real residual output
# --------------------------------------------------------------------------- #
@dataclass
class ResidualAssemblyResult(_Serializable):
    mode: str
    R: np.ndarray
    dof_labels: List[str] = field(default_factory=list)
    free_mask: Optional[np.ndarray] = None
    prescribed_idx: Optional[np.ndarray] = None
    sign_convention: str = SIGN_CONVENTION
    exported_reactions: Optional[np.ndarray] = None     # RF mapped to global dofs
    tangent: Optional[TangentResult] = None
    element_ids: List[int] = field(default_factory=list)
    _elem_dofs: List[np.ndarray] = field(default_factory=list)
    _elem_res: List[np.ndarray] = field(default_factory=list)
    diagnostics: Dict[str, Any] = field(default_factory=dict)

    # ---- derived quantities -----
    @property
    def ndof(self) -> int:
        return int(self.R.size)

    @property
    def free_residual(self) -> np.ndarray:
        return self.R[self.free_mask] if self.free_mask is not None else self.R

    @property
    def reactions(self) -> np.ndarray:
        if self.prescribed_idx is not None and self.prescribed_idx.size:
            return self.R[self.prescribed_idx]
        return np.zeros(0)

    @property
    def free_residual_norm(self) -> float:
        return float(np.linalg.norm(self.free_residual))

    @property
    def reaction_norm(self) -> float:
        return float(np.linalg.norm(self.reactions))

    def reaction_comparison(self, tol: float = 1e-6) -> Optional[Dict[str, Any]]:
        """Compare computed reactions to an exported RF (sign-aware)."""
        if self.exported_reactions is None or self.prescribed_idx is None \
                or not self.prescribed_idx.size:
            return None
        reac = self.reactions
        rf = self.exported_reactions[self.prescribed_idx]
        denom = max(float(np.linalg.norm(rf)), 1e-30)
        rel_same = float(np.linalg.norm(reac - rf) / denom)
        rel_opp = float(np.linalg.norm(reac + rf) / denom)
        rel = min(rel_same, rel_opp)
        return {"rel_same_sign": rel_same, "rel_opposite_sign": rel_opp,
                "best_rel_error": rel, "sign": "R = +RF" if rel_same <= rel_opp
                else "R = -RF", "tolerance": tol, "passed": rel < tol}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "kind": "ResidualAssemblyResult",
            "mode": self.mode,
            "ndof": self.ndof,
            "sign_convention": self.sign_convention,
            "free_residual_norm": self.free_residual_norm,
            "max_abs_residual": float(np.max(np.abs(self.R))) if self.R.size else 0.0,
            "reaction_norm": self.reaction_norm,
            "n_free_dofs": int(self.free_mask.sum()) if self.free_mask is not None else self.ndof,
            "n_prescribed_dofs": int(self.prescribed_idx.size) if self.prescribed_idx is not None else 0,
            "n_element_contributions": len(self.element_ids),
            "reaction_comparison": self.reaction_comparison(),
            "tangent": self.tangent.to_dict() if self.tangent else None,
            "dof_ordering_len": len(self.dof_labels),
            "diagnostics": {k: v for k, v in self.diagnostics.items()
                            if k not in ("element_residuals",)},
        }

    def arrays(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {"R": self.R}
        if self.free_mask is not None:
            out["free_mask"] = self.free_mask
        if self.prescribed_idx is not None:
            out["prescribed_idx"] = self.prescribed_idx
        if self.exported_reactions is not None:
            out["exported_reactions"] = self.exported_reactions
        if self.dof_labels:
            out["dof_ordering"] = np.array(self.dof_labels, dtype=object)
        # element residual contributions (CSR-like: values + dofs + pointers)
        if self.element_ids:
            ptr = [0]
            vals, dofs = [], []
            for ed, r in zip(self._elem_dofs, self._elem_res):
                dofs.append(np.asarray(ed, int))
                vals.append(np.asarray(r, float))
                ptr.append(ptr[-1] + len(r))
            out["element_ids"] = np.array(self.element_ids, dtype=int)
            out["element_res_ptr"] = np.array(ptr, dtype=int)
            out["element_res_dofs"] = (np.concatenate(dofs) if dofs
                                       else np.zeros(0, int))
            out["element_res_values"] = (np.concatenate(vals) if vals
                                         else np.zeros(0))
        if self.tangent is not None:
            for k, v in self.tangent.arrays().items():
                if v is not None:
                    out["tangent_" + k] = v
        return out

    def to_markdown(self) -> str:
        d = self.to_dict()
        L = ["# Residual Assembly Result", "",
             "- mode: `%s`" % self.mode,
             "- ndof: %d  (free %d / prescribed %d)"
             % (d["ndof"], d["n_free_dofs"], d["n_prescribed_dofs"]),
             "- sign convention: %s" % self.sign_convention,
             "- **||R_free||** = %.6e  (≈ 0 at a converged solution)"
             % d["free_residual_norm"],
             "- max|R| = %.6e" % d["max_abs_residual"],
             "- ||reactions|| = %.6e" % d["reaction_norm"],
             "- element contributions captured: %d" % d["n_element_contributions"]]
        rc = d["reaction_comparison"]
        if rc:
            L += ["", "## Reaction comparison (vs exported RF)",
                  "- best relative error: %.6e (%s)" % (rc["best_rel_error"], rc["sign"]),
                  "- passed (tol %g): %s" % (rc["tolerance"], rc["passed"])]
        L.append("")
        L.append("## Tangent")
        L.append(self.tangent.to_markdown() if self.tangent
                 else "T unavailable for this mode.\n")
        return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# 5. State / history output
# --------------------------------------------------------------------------- #
@dataclass
class StateResult(_Serializable):
    n_state_vars: int = 0
    increment_index: Optional[int] = None
    time: Tuple[float, float] = (0.0, 0.0)
    dtime: float = 0.0
    history_replay_status: str = "not-applicable"
    state_variables: Dict[int, np.ndarray] = field(default_factory=dict)   # eid -> (n_ip,n_state)
    state_sensitivities: Dict[int, np.ndarray] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "kind": "StateResult",
            "n_state_vars": self.n_state_vars,
            "increment_index": self.increment_index,
            "time": list(self.time),
            "dtime": self.dtime,
            "history_replay_status": self.history_replay_status,
            "n_elements_with_state": len(self.state_variables),
            "has_state_sensitivities": bool(self.state_sensitivities),
        }

    def arrays(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self.state_variables:
            eids = sorted(self.state_variables)
            out["state_element_ids"] = np.array(eids, dtype=int)
            out["state_values"] = np.stack([np.asarray(self.state_variables[e], float)
                                            for e in eids])
        if self.state_sensitivities:
            eids = sorted(self.state_sensitivities)
            out["state_sens_element_ids"] = np.array(eids, dtype=int)
            out["state_sens_values"] = np.stack(
                [np.asarray(self.state_sensitivities[e], float) for e in eids])
        return out

    def to_markdown(self) -> str:
        d = self.to_dict()
        return ("# State / History\n\n"
                "- n_state_vars: %d\n- increment index: %s\n- time: %s  dtime: %s\n"
                "- history replay: %s\n- elements with state: %d\n"
                "- state sensitivities: %s\n"
                % (d["n_state_vars"], d["increment_index"], d["time"], d["dtime"],
                   d["history_replay_status"], d["n_elements_with_state"],
                   d["has_state_sensitivities"]))


# --------------------------------------------------------------------------- #
# 6. Validation output
# --------------------------------------------------------------------------- #
@dataclass
class ValidationCheck:
    name: str
    quantity: str
    value: Optional[float]
    tolerance: Optional[float]
    passed: Optional[bool]
    detail: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"name": self.name, "quantity": self.quantity, "value": self.value,
                "tolerance": self.tolerance, "passed": self.passed,
                "detail": self.detail}


@dataclass
class ValidationReport(_Serializable):
    checks: List[ValidationCheck] = field(default_factory=list)

    def add(self, name, quantity, value=None, tolerance=None, passed=None, detail=""):
        self.checks.append(ValidationCheck(name, quantity, value, tolerance, passed, detail))
        return self

    def add_residual_vs_solver(self, free_norm, tol=1e-6):
        self.add("residual-vs-solver", "||R_free||", float(free_norm), tol,
                 float(free_norm) < tol, "free-DOF residual at the solver solution")
        return self

    def add_reaction_comparison(self, comparison: Optional[Dict[str, Any]]):
        if comparison is None:
            self.add("reaction-force", "rel error", None, None, None,
                     "no exported RF provided")
        else:
            self.add("reaction-force", "rel error", comparison["best_rel_error"],
                     comparison["tolerance"], comparison["passed"], comparison["sign"])
        return self

    def add_umat_replay(self, stress_rel=None, statev_rel=None, tol=1e-4):
        self.add("umat-replay-STRESS", "rel error", stress_rel, tol,
                 (stress_rel is not None and stress_rel < tol),
                 "material-point STRESS replay vs ODB (history in order)")
        self.add("umat-replay-STATEV", "rel error", statev_rel, tol,
                 (statev_rel is not None and statev_rel < tol),
                 "state-variable replay vs ODB")
        return self

    def add_uel_rhs_amatrx(self, rhs_rel=None, amatrx_rel=None, tol=1e-6):
        self.add("uel-RHS", "rel error", rhs_rel, tol,
                 (rhs_rel is not None and rhs_rel < tol),
                 "element_residual = -RHS")
        self.add("uel-AMATRX", "rel error", amatrx_rel, tol,
                 (amatrx_rel is not None and amatrx_rel < tol),
                 "AMATRX = dR/du")
        return self

    @property
    def overall_passed(self) -> Optional[bool]:
        decided = [c.passed for c in self.checks if c.passed is not None]
        if not decided:
            return None
        return all(decided)

    def to_dict(self) -> Dict[str, Any]:
        return {"kind": "ValidationReport",
                "overall_passed": self.overall_passed,
                "checks": [c.to_dict() for c in self.checks]}

    def to_markdown(self) -> str:
        L = ["# Validation Report", "",
             "overall: %s" % self.overall_passed, "",
             "| check | quantity | value | tol | passed | detail |",
             "|---|---|---|---|---|---|"]
        for c in self.checks:
            val = "%.3e" % c.value if isinstance(c.value, (int, float)) else "—"
            tol = "%g" % c.tolerance if isinstance(c.tolerance, (int, float)) else "—"
            pas = {True: "PASS", False: "FAIL", None: "pending"}[c.passed]
            L.append("| %s | %s | %s | %s | %s | %s |"
                     % (c.name, c.quantity, val, tol, pas, c.detail))
        return "\n".join(L) + "\n"
