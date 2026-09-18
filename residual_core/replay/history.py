"""General path-dependent residual-sensitivity replay for C3D8 meshes.

For ANY provider object (NPARAM parameters, NSTATV state variables, the
provider's parameter order) this module reconstructs a converged small-strain
static history and differentiates it with respect to the material parameters,
increment by increment, with the total-history recursion:

At increment n (displacements u_n, prescribed values u_c(t_n) independent of
p, loads independent of p) the residual is

    R_n(u_n, p) = sum_e sum_q w_q B_q^T sigma_n,q - F(t_n)

where sigma_n,q = UMAT(sigma_{n-1,q}, xi_{n-1,q}, STRAN = B_q u_{n-1},
DSTRAN = B_q (u_n - u_{n-1}), p). Differentiating R_n = 0 on the free DOFs,

    K_ff du_f/dp = - dR_f/dp |_(u_n fixed),     du_c/dp = 0,

with K = sum w B^T DDSDDE B (DDSDDE = dSTRESS/dDSTRAN from the provider's OTI
strain directions) and

    dR/dp |_(u_n fixed) = sum w B^T dsigma^A/dp,

where dsigma^A/dp is ONE call of UMAT_OTI_EVAL_TOTAL whose parameter
directions carry the incoming history derivatives dsigma_{n-1}/dp and
dxi_{n-1}/dp, dSTRAN/dp = B du_{n-1}/dp and dDSTRAN/dp = -B du_{n-1}/dp (the
"-B du_prev/dp" term: u_n is held fixed). First-order OTI is linear in the
seeds, so after the solve the total update derivatives follow exactly as

    dsigma_n/dp = dsigma^A/dp + DDSDDE   B du_n/dp,
    dxi_n/dp    = dxi^A/dp    + dxi/dDSTRAN B du_n/dp,

and the reactions as dRF_c/dp = dR_c/dp|_(u fixed) + K_cf du_f/dp.

Assembly is sparse (scipy.sparse); K_ff is factorised once per increment and
the factorisation is reused for every parameter column. No per-increment
dense K is formed or stored.

Two ways to obtain u_n:

- ``replay``: the recorded (ODB) displacements drive the kinematics. The
  replayed stress/state are checked against the recording at every point and
  increment and the free-DOF residual against its admissible size; any
  mismatch raises :class:`ReplayMismatch`. ``reequilibrate=True`` additionally
  Newton-polishes each increment from the recorded state to double-precision
  equilibrium (the correction is reported and bounded by the recording's
  single-precision resolution and Abaqus's force tolerance).
- ``solve``: Newton equilibrium in Python at given step times (used for the
  independent whole-model finite-difference reference with the ORIGINAL UMAT
  and for models without a recording).
"""
from __future__ import annotations

import time as _time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..formulations import c3d8_kernel as _kernel
from .history_inputs import PRECISION_FLOOR, HistoryModel, RecordedFields, UnsupportedFeature
from .history_material import HistoryMaterial
from .kinematics import batch_operators

#: Abaqus/Standard's default force-residual criterion R_n^alpha (5e-3 of the
#: time-averaged force); the largest residual an ODB increment may carry.
ABAQUS_RESIDUAL_RATIO = 5e-3


class ReplayMismatch(ValueError):
    """The recorded history is not reproduced by the replay within tolerance."""


def von_mises(stress):
    s = np.asarray(stress)
    return np.sqrt(0.5 * ((s[..., 0] - s[..., 1]) ** 2 + (s[..., 1] - s[..., 2]) ** 2
                          + (s[..., 2] - s[..., 0]) ** 2)
                   + 3.0 * (s[..., 3] ** 2 + s[..., 4] ** 2 + s[..., 5] ** 2))


def von_mises_gradient(stress):
    """d(sigma_vM)/d(sigma) for Voigt stress (tensor shear components); 0 where sigma_vM=0."""
    s = np.asarray(stress, dtype=float)
    vm = von_mises(s)
    mean = s[..., :3].sum(axis=-1, keepdims=True) / 3.0
    gradient = np.concatenate([1.5 * (s[..., :3] - mean), 3.0 * s[..., 3:]], axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        gradient = np.where(vm[..., None] > 0, gradient / vm[..., None], 0.0)
    return gradient


class HistoryEngine:
    """Mesh operators, sparse assembly and the provider for one model."""

    def __init__(self, model: HistoryModel, material: HistoryMaterial,
                 integration: str = "selective_reduced"):
        if material.nprops != len(model.props):
            raise UnsupportedFeature("the deck has %d USER MATERIAL constants; the provider contract "
                                     "declares NPROPS=%d" % (len(model.props), material.nprops))
        if material.nstatev != model.depvar:
            raise UnsupportedFeature("the deck declares *Depvar %d; the provider contract declares "
                                     "NSTATV=%d" % (model.depvar, material.nstatev))
        self.model, self.material, self.integration = model, material, integration
        self.ne = len(model.element_ids)
        self.nq = 8
        self.ndof = model.ndof
        X = model.coords[model.connectivity]                       # (ne, 8, 3)
        self.B, self.w = batch_operators(X, integration)             # (ne,8,6,24), (ne,8)
        self.edofs = (3 * model.connectivity[:, :, None] + np.arange(3)).reshape(self.ne, 24)
        self.rows = np.repeat(self.edofs, 24, axis=1).ravel()
        self.cols = np.tile(self.edofs, (1, 24)).ravel()
        self.scatter = sp.csr_matrix((np.ones(self.ne * 24), (self.edofs.ravel(), np.arange(self.ne * 24))),
                                     shape=(self.ndof, self.ne * 24))
        shape = np.array([_kernel.shape_functions(point) for point in _kernel.ABAQUS_C3D8_GAUSS.points])
        self.ip_coords = np.einsum("qa,eai->eqi", shape, X)
        self.volume = self.w.sum(axis=1)
        self.celent = np.repeat(self.volume ** (1.0 / 3.0), 8)
        self.noel = np.repeat(model.element_ids, 8).astype(np.int32)
        self.npt = np.tile(np.arange(1, 9), self.ne).astype(np.int32)
        constrained = np.zeros(self.ndof, dtype=bool)
        constrained[model.prescribed_dofs] = True
        self.constrained = model.prescribed_dofs
        self.free = np.flatnonzero(~constrained)
        if not len(self.free):
            raise ValueError("no free degrees of freedom")
        self.nparam = material.nparam

    # ---------------------------------------------------------------- kinematics
    def strain(self, u):
        """B u at every point -> (ne, 8, 6); u has shape (ndof,)."""
        return np.einsum("eqai,ei->eqa", self.B, u[self.edofs], optimize=True)

    def strain_param(self, du):
        """B du/dp -> (ne, 8, 6, nparam); du has shape (ndof, nparam)."""
        return np.einsum("eqai,eim->eqam", self.B, du[self.edofs], optimize=True)

    def force(self, stress):
        local = np.einsum("eqai,eqa,eq->ei", self.B, stress, self.w, optimize=True)
        return self.scatter @ local.ravel()

    def force_param(self, dstress):
        local = np.einsum("eqai,eqam,eq->eim", self.B, dstress, self.w, optimize=True)
        return self.scatter @ local.reshape(self.ne * 24, -1)

    def stiffness(self, ddsdde):
        local = np.einsum("eqai,eqab,eqbj,eq->eij", self.B, ddsdde, self.B, self.w, optimize=True)
        return sp.csr_matrix((local.ravel(), (self.rows, self.cols)), shape=(self.ndof, self.ndof))

    def flat(self, array):
        return array.reshape(self.ne * self.nq, *array.shape[2:])

    def unflat(self, array):
        return array.reshape(self.ne, self.nq, *array.shape[1:])


@dataclass
class IncrementRecord:
    number: int
    time: float
    u: np.ndarray
    stress: np.ndarray
    state: np.ndarray
    reaction: np.ndarray              # internal force minus external load, all DOFs
    residual_free_max: float
    residual_scale: float
    residual_limit: float
    iterations: int = 0
    du: Optional[np.ndarray] = None
    dstress: Optional[np.ndarray] = None
    dstate: Optional[np.ndarray] = None
    dreaction: Optional[np.ndarray] = None
    sensitivity_residual: float = float("nan")
    parity: Dict[str, float] = field(default_factory=dict)
    correction: float = 0.0


@dataclass
class HistoryResult:
    mode: str
    parameters: List[str]
    parameter_values: np.ndarray
    props: np.ndarray
    times: np.ndarray
    increments: List[IncrementRecord]
    integration: str
    material_path: str
    timings: Dict[str, float]
    tolerances: Dict[str, str]
    sensitivities: bool

    @property
    def max_scaled_residual(self) -> float:
        return max(r.residual_free_max / r.residual_scale for r in self.increments)

    def stacked(self, name):
        return np.stack([getattr(r, name) for r in self.increments])


def _tolerances_text():
    return {
        "stress_parity": ("|S_replay - S_odb| <= 8 eps32 max|S_n| + 8 eps32 ||D||_inf ||B||_inf max|U_n| "
                          "(storage rounding of S plus the effect of the single-precision U that drives "
                          "the kinematics; eps32 = 2^-23)"),
        "state_parity": "|SDV_replay - SDV_odb| <= 8 eps32 (max|SDV_n| + |d SDV| bound from U rounding)",
        "reaction_parity": "|RF_replay - RF_odb| <= 8 eps32 max|RF_n| + ||K||_inf 8 eps32 max|U_n|",
        "free_residual": ("max|R_free| <= 5e-3 * mean|nodal force| (Abaqus R_n^alpha default) "
                          "+ ||K_ff||_inf 8 eps32 max|U_n| (single-precision U)"),
        "solve_equilibrium": "max|R_free| <= rtol * max(max|R_c|, max|F|) with rtol given to solve",
    }


def run_history(engine: HistoryEngine, *, fields: Optional[RecordedFields] = None,
                times: Optional[np.ndarray] = None, props: Optional[np.ndarray] = None,
                sensitivities: bool = True, material_path: str = "oti",
                reequilibrate: bool = False, rtol: float = 1e-11, max_iterations: int = 60,
                newton_tangent: str = "oti", strict: bool = True, progress=None) -> HistoryResult:
    """Replay (``fields`` given) or solve (``times`` given) the history.

    ``material_path="regular"`` evaluates only the ORIGINAL UMAT (no
    sensitivities); ``newton_tangent`` chooses the Newton matrix in solve mode
    ("oti": the provider's consistent tangent, "material": whatever DDSDDE the
    ORIGINAL UMAT returns). The converged state depends only on the residual.
    ``strict=False`` records replay parity without raising (diagnostics only;
    the command-line and request paths always run strict).
    """
    if (fields is None) == (times is None):
        raise ValueError("give either recorded fields (replay) or step times (solve)")
    if material_path not in ("oti", "regular"):
        raise ValueError("material_path must be 'oti' or 'regular'")
    if sensitivities and material_path != "oti":
        raise ValueError("sensitivities need the OTI provider path")
    model, material = engine.model, engine.material
    props = np.array(model.props if props is None else props, dtype=float)
    replay = fields is not None
    step_times = fields.time if replay else np.asarray(times, dtype=float)
    if step_times[0] != 0 or np.any(np.diff(step_times) <= 0):
        raise ValueError("step times must start at 0 and increase strictly")
    ne, nq, npar, ns = engine.ne, engine.nq, engine.nparam, material.nstatev
    npts = ne * nq
    stress = np.zeros((npts, 6))
    state = np.zeros((npts, ns))
    dstress = np.zeros((npts, 6, npar))
    dstate = np.zeros((npts, ns, npar))
    u_prev = np.zeros(engine.ndof)
    du_prev = np.zeros((engine.ndof, npar))
    coords = engine.ip_coords.reshape(npts, 3)
    records: List[IncrementRecord] = []
    timings = {"material": 0.0, "assembly": 0.0, "factorisation_and_solve": 0.0, "checks": 0.0}
    started = _time.perf_counter()
    b_norm = np.abs(engine.B).sum(axis=-1).max()
    start_stiffness = None          # tangent at the start of the next increment (solve mode)

    for number in range(1, len(step_times)):
        t_prev, t_now = step_times[number - 1], step_times[number]
        dt = t_now - t_prev
        prescribed = model.prescribed_at(t_now)
        external = model.load_at(t_now)
        strain_prev = engine.flat(engine.strain(u_prev))
        if sensitivities:
            stran_dp = engine.flat(engine.strain_param(du_prev))
        common = dict(time=np.array([t_prev, t_prev]), dtime=dt, coords=coords,
                      celent=engine.celent, noel=engine.noel, npt=engine.npt, kstep=1, kinc=number)

        def evaluate(u, need_tangent=True):
            tic = _time.perf_counter()
            dstran = engine.flat(engine.strain(u)) - strain_prev
            if material_path == "oti" and (sensitivities or need_tangent):
                if sensitivities:
                    out = material.total(props, stress, state, strain_prev, dstran,
                                         dstress_in=dstress, dstate_in=dstate, stran_dp=stran_dp,
                                         dstran_dp=-stran_dp, **common)
                else:
                    zero = np.zeros((npts, 6, npar))
                    out = material.total(props, stress, state, strain_prev, dstran,
                                         dstress_in=zero, dstate_in=np.zeros((npts, ns, npar)),
                                         stran_dp=zero, dstran_dp=zero, **common)
            else:
                out = material.regular(props, stress, state, strain_prev, dstran, **common)
                if need_tangent and newton_tangent == "oti" and material_path == "regular":
                    zero = np.zeros((npts, 6, npar))
                    out["ddsdde"] = material.total(props, stress, state, strain_prev, dstran,
                                                   dstress_in=zero, dstate_in=np.zeros((npts, ns, npar)),
                                                   stran_dp=zero, dstran_dp=zero, **common)["ddsdde"]
            timings["material"] += _time.perf_counter() - tic
            tic = _time.perf_counter()
            internal = engine.force(engine.unflat(out["stress"]))
            stiffness = engine.stiffness(engine.unflat(out["ddsdde"])) if need_tangent else None
            timings["assembly"] += _time.perf_counter() - tic
            return out, internal - external, stiffness

        iterations = 0
        correction = 0.0
        if replay:
            u = fields.U[number].reshape(-1).copy()
            u[engine.constrained] = prescribed
            u_recorded = u.copy()
            out, residual, stiffness = evaluate(u)
        else:
            # Linear predictor (as Abaqus's first iteration of an increment): move the
            # free DOFs with the start-of-increment tangent so the trial state the
            # UMAT sees is near equilibrium, not a jump at the constrained DOFs only.
            if start_stiffness is None:
                start_stiffness = evaluate(u_prev)[2]
            tic = _time.perf_counter()
            jump = np.zeros(engine.ndof)
            jump[engine.constrained] = prescribed - u_prev[engine.constrained]
            rhs = (external - model.load_at(t_prev) - start_stiffness @ jump)[engine.free]
            u = u_prev + jump
            u[engine.free] += spla.splu(start_stiffness[engine.free][:, engine.free].tocsc()).solve(rhs)
            timings["factorisation_and_solve"] += _time.perf_counter() - tic
            out, residual, stiffness = evaluate(u)
        scale = max(np.abs(residual[engine.constrained]).max(initial=0.0),
                    np.abs(external).max(initial=0.0), 1e-300)
        if not np.all(np.isfinite(residual)):
            raise ReplayMismatch("increment %d: the %s UMAT returned non-finite stresses at the %s state"
                                 % (number, material_path.upper(), "recorded" if replay else "predicted"))
        if not replay or reequilibrate:
            for iterations in range(1, max_iterations + 1):
                if np.abs(residual[engine.free]).max() <= rtol * scale:
                    break
                tic = _time.perf_counter()
                matrix = stiffness[engine.free][:, engine.free].tocsc()
                step = spla.splu(matrix).solve(-residual[engine.free])
                timings["factorisation_and_solve"] += _time.perf_counter() - tic
                baseline = np.linalg.norm(residual[engine.free])
                for halving in range(12):
                    trial = u.copy()
                    trial[engine.free] += 0.5 ** halving * step
                    trial_out, trial_residual, trial_stiffness = evaluate(trial)
                    if np.linalg.norm(trial_residual[engine.free]) < baseline:
                        break
                else:
                    # no step along the Newton direction lowers the residual: the model
                    # has reached its roundoff floor (or the tangent is wrong); stop
                    # instead of iterating on noise
                    raise ReplayMismatch(
                        "increment %d: Newton stagnated at max|R_free| = %.3e (requested %.1e x %.3e = "
                        "%.3e); a tolerance below this model's roundoff floor cannot be met"
                        % (number, np.abs(residual[engine.free]).max(), rtol, scale, rtol * scale))
                if not np.all(np.isfinite(trial_residual)):
                    raise ReplayMismatch("increment %d: Newton reached non-finite stresses" % number)
                u, out, residual, stiffness = trial, trial_out, trial_residual, trial_stiffness
                scale = max(np.abs(residual[engine.constrained]).max(initial=0.0),
                            np.abs(external).max(initial=0.0), 1e-300)
            else:
                raise ReplayMismatch("increment %d: Newton did not reach max|R_free| <= %.1e x %.3e "
                                     "in %d iterations (now %.3e)" % (number, rtol, scale, max_iterations,
                                                                      np.abs(residual[engine.free]).max()))
            if replay:
                correction = float(np.abs(u - u_recorded).max())

        tic = _time.perf_counter()
        residual_free = float(np.abs(residual[engine.free]).max())
        u_max = float(np.abs(u).max())
        stiffness_norm = float(abs(stiffness).sum(axis=1).max()) if stiffness is not None else 0.0
        record = IncrementRecord(number=number, time=float(t_now), u=u.copy(),
                                 stress=engine.unflat(out["stress"]).copy(),
                                 state=engine.unflat(out["state"]).copy(), reaction=residual.copy(),
                                 residual_free_max=residual_free, residual_scale=scale,
                                 residual_limit=rtol * scale, iterations=iterations,
                                 correction=correction)
        if replay:
            failures = _check_replay(record, fields, number, engine, out, stiffness_norm, b_norm,
                                     u_max, reequilibrate, rtol, max(r.correction for r in records + [record]))
            if failures and strict:
                raise ReplayMismatch("increment %d does not reproduce the recording: %s"
                                     % (number, "; ".join(failures)))
        timings["checks"] += _time.perf_counter() - tic

        if sensitivities:
            tic = _time.perf_counter()
            rhs = engine.force_param(engine.unflat(out["dstress_dp"]))          # (ndof, npar)
            k_ff = stiffness[engine.free][:, engine.free].tocsc()
            factor = spla.splu(k_ff)
            du = np.zeros((engine.ndof, npar))
            du[engine.free] = factor.solve(np.ascontiguousarray(-rhs[engine.free]))
            dreaction = rhs + stiffness @ du
            record.sensitivity_residual = float(np.abs(dreaction[engine.free]).max()
                                                / max(np.abs(rhs).max(), 1e-300))
            timings["factorisation_and_solve"] += _time.perf_counter() - tic
            tic = _time.perf_counter()
            dstrain = engine.flat(engine.strain_param(du))                     # (npts, 6, npar)
            dstress = out["dstress_dp"] + np.einsum("iab,ibm->iam", out["ddsdde"], dstrain)
            dstate = out["dstate_dp"] + np.einsum("isb,ibm->ism", out["dstate_ddstran"], dstrain)
            timings["assembly"] += _time.perf_counter() - tic
            record.du, record.dreaction = du, dreaction
            record.dstress, record.dstate = engine.unflat(dstress), engine.unflat(dstate)
            du_prev = du
        stress, state, u_prev = out["stress"], out["state"], u
        start_stiffness = stiffness
        records.append(record)
        if progress:
            progress(record)

    timings["total"] = _time.perf_counter() - started
    return HistoryResult(
        mode="replay" + ("+reequilibrate" if reequilibrate else "") if replay else "solve",
        parameters=list(material.params), parameter_values=material.parameter_values(props),
        props=props, times=step_times, increments=records, integration=engine.integration,
        material_path=material_path, timings=timings, tolerances=_tolerances_text(),
        sensitivities=sensitivities)


def _check_replay(record, fields, number, engine, out, stiffness_norm, b_norm, u_max,
                  reequilibrate, rtol, correction):
    """Compare the replay with the recording; return the list of excesses."""
    eps = PRECISION_FLOOR
    stress_odb = fields.S[number]
    state_odb = fields.SDV[number]
    d_norm = float(np.abs(out["ddsdde"]).sum(axis=-1).max())
    u_rounding = 8 * eps * max(u_max, np.abs(fields.U[number]).max())
    moved = u_rounding + correction          # how far the kinematics may sit from the recording
    stress_limit = 8 * eps * np.abs(stress_odb).max() + d_norm * b_norm * moved
    stress_error = float(np.abs(record.stress - stress_odb).max())
    parity = {"stress_max_abs": stress_error, "stress_limit": float(stress_limit),
              "stress_max_rel": stress_error / max(np.abs(stress_odb).max(), 1e-300)}
    if state_odb.size:
        tangent = np.abs(out["dstate_ddstran"]).sum(axis=-1).max() if "dstate_ddstran" in out else 0.0
        state_limit = 8 * eps * np.abs(state_odb).max() + tangent * b_norm * moved + 1e-30
        state_error = float(np.abs(record.state - state_odb).max())
        parity.update(state_max_abs=state_error, state_limit=float(state_limit))
    rf_odb = fields.RF[number].reshape(-1)[engine.constrained]
    rf_replay = record.reaction[engine.constrained]
    rf_error = float(np.abs(rf_replay - rf_odb).max())
    rf_limit = 8 * eps * np.abs(rf_odb).max() + stiffness_norm * moved
    parity.update(reaction_max_abs=rf_error, reaction_limit=float(rf_limit))
    nodal = np.abs(record.reaction).reshape(-1, 3)
    mean_force = float(np.linalg.norm(nodal, axis=1).mean())
    if reequilibrate:
        residual_limit = rtol * record.residual_scale
    else:
        residual_limit = ABAQUS_RESIDUAL_RATIO * mean_force + stiffness_norm * u_rounding
    record.residual_limit = float(residual_limit)
    parity["free_residual_limit"] = float(residual_limit)
    record.parity = parity
    failures = []
    if stress_error > stress_limit:
        failures.append("stress %.3e > %.3e" % (stress_error, stress_limit))
    if state_odb.size and parity["state_max_abs"] > parity["state_limit"]:
        failures.append("state %.3e > %.3e" % (parity["state_max_abs"], parity["state_limit"]))
    if rf_error > rf_limit:
        failures.append("reaction %.3e > %.3e" % (rf_error, rf_limit))
    if record.residual_free_max > residual_limit:
        failures.append("free residual %.3e > %.3e" % (record.residual_free_max, residual_limit))
    return failures
