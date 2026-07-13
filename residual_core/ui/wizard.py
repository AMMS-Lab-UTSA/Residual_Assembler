"""High-level workflow + the ``ResidualProblem`` facade.

This is the "ask for the least data possible, auto-detect the rest" layer. A user
does::

    from residual_core import ResidualProblem
    p = ResidualProblem.from_abaqus("model.inp")
    p.inspect()                       # what can be assembled, what's missing
    p.attach_results("results.odb")   # or a dict / CSV / JSON of fields
    R = p.assemble(mode="stress-driven")

or for a material replay::

    p = ResidualProblem.from_abaqus("model.inp")
    p.attach_subroutine("umat.for")
    R = p.assemble(mode="material-replay")

The facade instantiates the low-level pieces (Model, DofManager, Assembler,
registries) itself, chooses a formulation policy appropriate to the requested
mode, and reports the *minimum* missing input when it cannot proceed. No physics
lives here — it orchestrates registries and the agnostic core.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, Optional

import numpy as np

from ..core.model import from_abaqus as _from_abaqus, Model
from ..core.dof_manager import DofManager
from ..core.assembler import Assembler
from ..core.state_manager import StateManager
from ..core import requirements as _req
from ..core.diagnostics import inspect_model, InspectionReport
from ..formulations.registry import build_formulation_registry
from ..materials.registry import build_material_registry
from ..materials.base import MaterialBinding
from ..io import abaqus_inp_parser
from ..io import neutral_model_io


class ResidualProblem:
    def __init__(self, model: Model, formulation_registry=None,
                 material_registry=None, config=None):
        self.model = model
        self.forms = formulation_registry or build_formulation_registry()
        self.mats = material_registry or build_material_registry()
        self.config = config
        self._fields: Dict[str, Any] = {}       # e.g. {'stress_ip': {...}}
        self._subroutine: Optional[str] = None
        self._solution: Optional[np.ndarray] = None
        self._raw_materials = dict(model.materials)   # keep detection info

    # ---- construction ---------------------------------------------------
    @classmethod
    def from_abaqus(cls, path: str, config=None) -> "ResidualProblem":
        am = abaqus_inp_parser.parse_inp(path)
        # bind nothing yet — the policy is chosen per requested mode at assemble
        model = _from_abaqus(am, formulation_policy=lambda et: None)
        # keep the parser materials (carry the user_material flag for detection)
        model.materials = dict(am.materials)
        prob = cls(model, config=config)
        prob._abaqus = am
        return prob

    @classmethod
    def from_neutral(cls, path: str, config=None) -> "ResidualProblem":
        model = neutral_model_io.load_model(path)
        return cls(model, config=config)

    def save_neutral(self, path: str) -> None:
        neutral_model_io.save_model(self.model, path)

    # ---- attachments ----------------------------------------------------
    def attach_results(self, source: Any) -> "ResidualProblem":
        """Attach exported solver fields for stress-driven mode.

        ``source`` may be a dict already shaped as ``{'stress_ip': {eid:(nip,6)}}``
        (or just ``{eid:(nip,6)}``), or a path to a JSON export. Binary ODB
        reading is delegated to ``io/abaqus_odb_export`` and requires Abaqus.
        """
        if isinstance(source, dict):
            self._fields = source if "stress_ip" in source else {"stress_ip": source}
        elif isinstance(source, str) and source.lower().endswith(".json"):
            with open(source, "r", encoding="utf-8") as fh:
                d = json.load(fh)
            stress = d.get("stress_ip", d)
            self._fields = {"stress_ip": {int(k): np.asarray(v, float)
                                          for k, v in stress.items()}}
        else:
            raise NotImplementedError(
                "attach_results: binary ODB reading needs Abaqus (see "
                "io/abaqus_odb_export.py). Provide a JSON/dict field export offline.")
        return self

    def attach_subroutine(self, path: str) -> "ResidualProblem":
        self._subroutine = path
        return self

    def set_solution(self, U) -> "ResidualProblem":
        self._solution = np.asarray(U, float)
        return self

    # ---- inspection & requirements -------------------------------------
    def inspect(self, print_it: bool = True) -> InspectionReport:
        rep = inspect_model(self.model, self.forms, self.mats,
                            attached_subroutine=self._subroutine is not None)
        if print_it:
            print(rep.render())
        return rep

    def _availability(self, mode: str) -> Dict[str, bool]:
        m = _req.canonical_mode(mode)
        has_mesh = bool(self.model.nodes) and bool(self.model.elements)
        has_stress = "stress_ip" in self._fields
        has_sub = self._subroutine is not None
        # a material is "built-in resolvable" if it is not a user material
        any_builtin = any(not getattr(o, "user_material", False)
                          for o in self._raw_materials.values()) or not self._raw_materials
        any_user = any(getattr(o, "user_material", False)
                       for o in self._raw_materials.values())
        has_params = any(list(getattr(o, "constants", []) or [])
                         for o in self._raw_materials.values())
        if m == "stress-driven":
            return {"mesh": has_mesh,
                    "dof_field": True,          # zeros acceptable for small config
                    "element_field": has_stress}
        if m == "material-replay":
            material_ok = has_sub or (any_builtin and not any_user)
            return {"mesh": has_mesh,
                    "solution_history": self._solution is not None or (any_builtin and not any_user),
                    "material_model": material_ok,
                    "material_parameters": has_params or (any_builtin and not any_user),
                    "state_prev": (not any_user) or has_sub,
                    "time_increments": True}
        if m == "direct-residual":
            return {"mesh": has_mesh, "element_dofs": True,
                    "uel_routine": has_sub}
        if m == "formulation":
            has_sections = any(getattr(o, "section", None) for o in self.model.materials.values())
            return {"formulation_backend": True, "section_properties": has_sections}
        return {}

    def requirements(self, mode: str) -> "_req.RequirementsReport":
        return _req.evaluate_requirements(mode, self._availability(mode))

    # ---- formulation policy per mode -----------------------------------
    def _policy_for_mode(self, mode: str):
        m = _req.canonical_mode(mode)
        want = {"stress-driven": "stress-driven",
                "material-replay": "material-replay",
                "formulation": "formulation",
                "direct-residual": "direct-residual"}.get(m)

        def policy(etype: str):
            for name in self.forms.find_for_element(etype):
                spec = getattr(self.forms.get(name), "spec", None)
                if spec and want in spec.supported_modes:
                    return name
            return None
        return policy

    def _bind_materials_for_replay(self):
        """Attach a built-in elastic binding for non-user materials (offline)."""
        from ..materials.elastic_adapter import IsotropicElastic
        bindings = {}
        for name, obj in self._raw_materials.items():
            if getattr(obj, "user_material", False):
                # CP/UMAT: needs compiled Fortran; leave as-is (assemble will error)
                bindings[name] = obj
            else:
                consts = list(getattr(obj, "constants", []) or []) or [200000.0, 0.3]
                bindings[name] = MaterialBinding(IsotropicElastic(), consts, 0, name)
        return bindings

    # ---- assemble -------------------------------------------------------
    def _run(self, mode: str, U=None, compute_tangent: bool = False,
             collect_elements: bool = False):
        """Bind formulations/materials for ``mode``, build the DOF manager and
        assembler, and assemble. Returns (R, K, diag, dm). Raises a clear,
        minimal error if a required input is missing."""
        report = self.requirements(mode)
        if not report.runnable:
            nxt = report.minimum_next
            raise RuntimeError(
                "cannot assemble in mode '%s': minimum missing input -> %s. %s"
                % (report.mode, nxt.label,
                   ("Recommendation: provide %s." % (nxt.recommendation or nxt.label))))

        m = _req.canonical_mode(mode)
        policy = self._policy_for_mode(mode)
        # (re)bind element formulations under this mode's policy
        self.model.element_formulation = {}
        for eid, el in self.model.elements.items():
            fk = policy(el.etype)
            if fk is not None:
                self.model.element_formulation[eid] = fk

        if m == "material-replay":
            self.model.materials = self._bind_materials_for_replay()
            # honest guard: user-material replay needs compiled Fortran here
            for name, obj in self.model.materials.items():
                if not isinstance(obj, MaterialBinding):
                    raise RuntimeError(
                        "material-replay for user material '%s' needs the compiled "
                        "UMAT (Intel ifort + Abaqus). Use mode='stress-driven' with "
                        "an exported field for an offline residual." % name)

        dm = DofManager.for_model(self.model, self.forms)
        nstate = 0
        asm = Assembler(self.model, dm, self.forms, StateManager(8, nstate))
        if U is None:
            U = self._solution if self._solution is not None else np.zeros(dm.ndof)
        U = np.asarray(U, float)

        opts = {"config": "small"} if m == "stress-driven" else {}
        R, K, diag = asm.assemble(U, fields=self._fields or None,
                                  compute_tangent=compute_tangent, options=opts,
                                  collect_elements=collect_elements)
        self.last_diagnostics = diag
        self.dof_manager = dm
        return R, K, diag, dm

    def assemble(self, mode: str = "stress-driven", U=None,
                 compute_tangent: bool = False):
        """Assemble the global residual for ``mode``. Raises a clear, minimal
        error if a required input is missing."""
        R, K, _diag, _dm = self._run(mode, U, compute_tangent=compute_tangent)
        return (R, K) if compute_tangent else R

    # ---- structured output / sensitivity package ------------------------
    def result(self, mode: str = "stress-driven", U=None,
               compute_tangent: bool = True):
        """Assemble and return a :class:`ResidualAssemblyResult` (R, reactions,
        element contributions, and the tangent T when the backend can assemble
        one). This is the serializable real-side output (contract §2, §3)."""
        from ..core.results import (ResidualAssemblyResult, TangentResult,
                                    TangentSource, dof_ordering)
        from ..core import constraints as _constraints

        R, K, diag, dm = self._run(mode, U, compute_tangent=compute_tangent,
                                   collect_elements=True)
        free_mask, pres_idx, _pres = _constraints.partition(self.model, dm)
        labels = dof_ordering(dm)

        # a tangent is only real if at least one element actually contributed one
        tangent_meaningful = K is not None and diag.get("tangent_contributions", 0) > 0
        if tangent_meaningful:
            tangent = TangentResult(
                T=K, source=TangentSource.BACKEND_ASSEMBLED, dof_labels=labels,
                free_mask=free_mask, prescribed_idx=pres_idx,
                note="assembled from the bound formulation backends")
        else:
            tangent = TangentResult(
                source=TangentSource.UNAVAILABLE, dof_labels=labels,
                free_mask=free_mask, prescribed_idx=pres_idx,
                note=("this mode (%s) provides no material/element tangent; "
                      "supply a solver-exported T or use finite-difference" % mode))

        exported = self._exported_reactions_vector(dm)
        er = diag.get("element_residuals", {})
        eids = sorted(er)
        res = ResidualAssemblyResult(
            mode=_req.canonical_mode(mode), R=R, dof_labels=labels,
            free_mask=free_mask, prescribed_idx=pres_idx, tangent=tangent,
            exported_reactions=exported,
            element_ids=eids,
            _elem_dofs=[er[e][0] for e in eids],
            _elem_res=[er[e][1] for e in eids],
            diagnostics={k: v for k, v in diag.items() if k != "element_residuals"})
        return res

    def _exported_reactions_vector(self, dm):
        """Map attached exported reactions {node: [rf...]} onto a global vector."""
        reactions = self._fields.get("reactions") if self._fields else None
        if not reactions:
            return None
        rf = np.zeros(dm.ndof)
        for node_s, comps in reactions.items():
            nid = int(node_s)
            try:
                gdofs = dm.node_dofs(nid)
            except Exception:
                continue
            for k, val in enumerate(comps):
                if k < len(gdofs):
                    rf[gdofs[k]] = float(val)
        return rf

    def solve_newton(self, mode: str = "formulation", u0=None,
                     max_iter: int = 100, tol: float = 1e-11):
        """Solve the real nonlinear problem R(u)=0 by Newton-Raphson on the free
        partition (step 1 of the residual method). Returns the converged u and
        stores it as the problem's solution. A nonzero initial guess is used by
        default (the cubic spring has zero tangent at u=0)."""
        from ..core import constraints as _constraints
        _R0, _K0, _d0, dm = self._run(mode, None, compute_tangent=False)
        U = np.ones(dm.ndof) if u0 is None else np.asarray(u0, float).copy()
        for _ in range(max_iter):
            R, K, _diag, dm = self._run(mode, U, compute_tangent=True)
            free, _pres, _ = _constraints.partition(self.model, dm)
            free_R = R[free] if free.any() else R
            if float(np.linalg.norm(free_R)) < tol:
                break
            Kff = K[np.ix_(free, free)] if free.any() else K
            dU = np.linalg.solve(Kff, -free_R)
            U = U.copy()
            if free.any():
                U[free] += dU
            else:
                U += dU
        self._solution = U
        self.dof_manager = dm
        return U

    def finite_difference_sensitivity(self, mode: str, parameters,
                                      h: float = 1e-6, u0=None):
        """du/da_i by central-ish finite differences: re-solve the real problem
        with each parameter perturbed and difference the converged solutions.
        Returns an (ndof, m) array. Restores the model and solution afterwards."""
        params = list(parameters)
        base = self.solve_newton(mode, u0=u0)
        saved_solution = self._solution
        ndof = base.size
        fd = np.zeros((ndof, len(params)))
        for col, pname in enumerate(params):
            matname, key = (pname.rsplit(".", 1) if "." in pname else (None, pname))
            binding = self.model.materials.get(matname)
            sec = binding if isinstance(binding, dict) else getattr(binding, "section", None)
            if not (isinstance(sec, dict) and key in sec):
                continue
            v0 = float(sec[key])
            dv = h * max(1.0, abs(v0))
            sec[key] = v0 + dv
            try:
                u_plus = self.solve_newton(mode, u0=base)
            finally:
                sec[key] = v0
            fd[:, col] = (u_plus - base) / dv
        self._solution = saved_solution
        return fd

    def sensitivity_package(self, mode: str = "stress-driven",
                            parameters=None, max_order: int = 1,
                            algebra: str = "OTI", U=None,
                            generate_rhs: bool = False, rhs_provider=None,
                            fd_check: bool = False, fd_step: float = 1e-6,
                            backend: str = "otilib"):
        """Build the HYPAD deliverable: the ingredients of ``T U^(p) = -R^(p)``.

        With ``generate_rhs=False`` (default) this emits the assembled residual +
        tangent (real side) and an *empty* R^(p) contract (shapes, parameter→
        direction maps, algebra metadata) ready for a hypercomplex evaluation to
        fill later.

        With ``generate_rhs=True`` a provider actually computes R^(p) from a
        residual evaluation. ``backend`` selects it:
          - ``"otilib"`` (default, production): genuine OTILib scalars, arbitrary
            order, all parameters seeded simultaneously, order-by-order u* update.
            Requires OTILib installed — otherwise a clear error is raised (there is
            **no** silent Dual1 fallback).
          - ``"dual1"``: the legacy first-order smoke test (order 1 only).

        ``fd_check`` adds a finite-difference cross-check to the validation report.
        If the model cannot be assembled, returns a package flagged
        ``runnable=False`` with the single minimum missing input.
        """
        from ..core.sensitivity_package import (SensitivityPackage,
                                                SensitivityRHSResult, AlgebraMetadata)
        from ..core.results import ValidationReport, StateResult

        insp = inspect_model(self.model, self.forms, self.mats,
                             attached_subroutine=self._subroutine is not None).as_dict()
        report = self.requirements(mode)
        if not report.runnable:
            nxt = report.minimum_next
            return SensitivityPackage(
                mode=_req.canonical_mode(mode), runnable=False, inspection=insp,
                minimum_missing=(nxt.recommendation or nxt.label) if nxt else None)

        params = list(parameters) if parameters else self._default_parameters()

        # Step 1: obtain a converged real solution when we will generate R^(p).
        if generate_rhs and U is None and self._solution is None:
            U = self.solve_newton(mode)
        u_real = U if U is not None else self._solution

        res = self.result(mode, U=u_real, compute_tangent=True)

        if generate_rhs:
            provider = rhs_provider or self._select_rhs_provider(backend, max_order)
            rhs = provider.evaluate_rhs(self, u_real, res.tangent, params,
                                        order=max_order)
        else:
            algebra_meta = AlgebraMetadata(algebra=algebra, n_bases=len(params),
                                           truncation_order=max_order)
            rhs = SensitivityRHSResult(
                ndof=res.ndof, algebra=algebra_meta,
                parameter_map={n: i + 1 for i, n in enumerate(params)},
                max_order=max_order)

        validation = ValidationReport()
        validation.add_residual_vs_solver(res.free_residual_norm)
        validation.add_reaction_comparison(res.reaction_comparison())

        state = StateResult(history_replay_status=(
            "not-applicable (stateless/offline)" if _req.canonical_mode(mode)
            != "material-replay" else "single-increment (history pending)"))

        pkg = SensitivityPackage(
            mode=_req.canonical_mode(mode), runnable=True, residual=res,
            sensitivity=rhs, state=state, validation=validation, inspection=insp)

        if generate_rhs and fd_check and pkg.system_ready(1):
            U1 = pkg.solve(1)
            fd = self.finite_difference_sensitivity(mode, params, h=fd_step,
                                                    u0=u_real)
            denom = max(float(np.linalg.norm(fd)), 1e-30)
            rel = float(np.linalg.norm(U1 - fd) / denom)
            validation.add("sensitivity-vs-finite-difference", "rel error",
                           rel, 1e-4, rel < 1e-4,
                           "solved U^(1) vs re-solved finite differences")

        return pkg

    def _select_rhs_provider(self, backend: str, max_order: int):
        """Pick the RHS provider. OTILib is the production backend; Dual1 is the
        legacy first-order smoke test. Missing OTILib is reported cleanly — there
        is NO silent fallback to Dual1."""
        b = (backend or "otilib").lower()
        if b in ("otilib", "oti", "hypad"):
            from ..core.oti_rhs_provider import OtiLibRHSProvider
            from ..algebra.otilib_adapter import otilib_available, otilib_status
            if not otilib_available():
                raise RuntimeError(otilib_status()["error"])
            return OtiLibRHSProvider()
        if b in ("dual1", "dual"):
            if max_order > 1:
                raise RuntimeError(
                    "backend='dual1' supports order 1 only (legacy smoke test). "
                    "Use backend='otilib' for order >= 2.")
            from ..core.rhs_provider import default_rhs_provider
            return default_rhs_provider()
        raise ValueError("unknown sensitivity backend %r (use 'otilib' or 'dual1')"
                         % backend)

    def _default_parameters(self):
        """Design parameters detected from the model's materials (names of
        section/constant entries). Falls back to a generic single parameter."""
        names = []
        for mname, obj in self.model.materials.items():
            sec = getattr(obj, "section", None)
            if isinstance(sec, dict):
                names.extend("%s.%s" % (mname, k) for k in sec)
            consts = list(getattr(obj, "constants", []) or [])
            names.extend("%s.PROPS[%d]" % (mname, i) for i in range(len(consts)))
        return names or ["parameter_1"]
