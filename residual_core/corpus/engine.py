"""Assemble, solve and differentiate a C3D8 boundary-value problem driven by a
corpus UMAT. Equations: ``corpus_campaign/batches/B1/noether/DESIGN.md`` (small
strain, B1) and ``.../B2/noether/DESIGN.md`` (finite strain WITH history under
the Abaqus/Standard NLGEOM contract). Short form:

Residual (no external load; prescribed DOFs carry the reaction)::

    R(u_n; p, xi_{n-1}) = A_e f_e(u_n; sigma_q)
    sigma_q = UMAT(STRESS_in, STATEV_{n-1}, STRAN_in, DSTRAN, DFGRD0, DFGRD1, DROT; p)

Kinematics presets (``KINEMATICS``):

``small``        F0 = F1 = DROT = I, DSTRAN = B (u_n - u_{n-1}), reference B (B1).
``b1_legacy``    finite, full 2x2x2 integration, DSTRAN = sym((F1-F0) M^-1),
                 DROT = I, STRESS/STRAN NOT rotated (the B1 convention).
``nlgeom_full``  finite, full integration, Abaqus rotation: DROT Hughes-Winget,
                 STRESS_in = DROT sigma_{n-1} DROT^T, STRAN likewise.
``nlgeom``       Abaqus C3D8: mean-dilatation Fbar to the UMAT, DSTRAN with the
                 element-averaged volumetric increment, DROT from the unmodified
                 spin, residual sum w0 Jbar Bbar^T sigma -- matched to an Abaqus
                 run (B2 probe). ``nlgeom_centroid``: the centroid operators of
                 the RA replay (commit 30a5ab4), for comparison only.

Local parameter sensitivity, u_n and every incoming quantity held fixed::

    dR/dp |_(u_n, sigma_{n-1}, xi_{n-1}, STRAN_{n-1}, u_{n-1}) = A_e f_e(u_n; (d sigma/d p)_local)

Total sensitivity through the history: with seeds carrying d(.)_{n-1}/dp
(through DFGRD0, DSTRAN, DROT, the pre-rotated STRESS and STRAN),
solve ``K_ff du_f/dp = -dR_f/dp|_(u_n)`` (``du_c/dp = 0``), then re-evaluate the
update with the full seeds to obtain d sigma_n/dp, d xi_n/dp.

``K`` is the EXACT Jacobian of R: small strain ``sum B^T DDSDDE B detJ w``
with DDSDDE = d sigma/d DSTRAN from the OTI directions; finite strain the
fixed-stress derivative of f (``C3D8Nlgeom.dforce_fixed_stress``) plus
``sum Bbar^T (d sigma/d Fbar1 : d Fbar1/du)``, d sigma/d Fbar1 from OTI seeds on
every UMAT input that moves with Fbar1 (DFGRD1, DSTRAN, pre-rotated STRESS and
STRAN). DROT itself has no OTI seed in the provider: a source that READS DROT
is refused derivative claims (``provider.reads_drot``).

Frame rotation (objectivity test): a problem may carry ``rotations`` -- one
proper orthogonal Q per increment. The unknowns V are then displacements in the
problem frame and the body is placed at ``x = Q (X + V)``; the residual is
returned to the problem frame node by node (``Q^T R``). An objective material
and kinematics give the same V, the same frame reactions, and sigma' = Q sigma Q^T.
"""
from __future__ import annotations

import time as _time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from residual_core.formulations import c3d8_kernel as kernel
from residual_core.formulations import c3d8_nlgeom as nl
from residual_core.formulations.c3d8_sensitivity import element_dR_dp

from .mesh import Problem, material_coordinates
from .provider import CorpusProvider, MaterialCallError

__all__ = ["Assembly", "HistoryRun", "NewtonFailed", "run_history", "Kinematics",
           "KINEMATICS", "default_kinematics"]

_GAUSS = kernel.ABAQUS_C3D8_GAUSS
_EYE = np.eye(3)


_EPS = float(np.finfo(float).eps)
#: multiple of eps * |K_ff|_inf * L accepted as converged (see run_history)
ROUNDOFF_FLOOR_FACTOR = 8.0


class NewtonFailed(RuntimeError):
    pass


voigt_eng = nl.voigt_eng


@dataclass(frozen=True)
class Kinematics:
    name: str
    finite: bool
    integration: str = "full"          # full | mean_dilatation | centroid
    rotation: str = "none"             # none | abaqus

    def describe(self) -> str:
        if not self.finite:
            return "small strain, F0=F1=DROT=I, reference B, full 2x2x2 Gauss"
        return ("finite strain (current configuration), integration=%s, rotation=%s%s"
                % (self.integration, self.rotation,
                   "" if self.rotation == "abaqus" else
                   " (DROT=I, STRESS/STRAN not pre-rotated: B1 convention)"))


KINEMATICS = {
    "small": Kinematics("small", False),
    "b1_legacy": Kinematics("b1_legacy", True, "full", "none"),
    "nlgeom_full": Kinematics("nlgeom_full", True, "full", "abaqus"),
    # Abaqus C3D8 under NLGEOM (B2 probe: UMAT inputs <= 6.4e-15, reactions 4.7e-8)
    "nlgeom": Kinematics("nlgeom", True, "mean_dilatation", "abaqus"),
    # the centroid operators of residual_core.replay.kinematics (commit 30a5ab4):
    # NOT what Abaqus does on a distorted element; kept to measure the difference
    "nlgeom_centroid": Kinematics("nlgeom_centroid", True, "centroid", "abaqus"),
}


def default_kinematics(finite: bool) -> Kinematics:
    return KINEMATICS["nlgeom" if finite else "small"]


def _resolve(kinematics, finite: bool) -> Kinematics:
    if kinematics is None:
        return default_kinematics(finite)
    if isinstance(kinematics, str):
        kinematics = KINEMATICS[kinematics]
    if kinematics.finite != bool(finite):
        raise ValueError("kinematics %s does not match the case (finite=%s)"
                         % (kinematics.name, finite))
    return kinematics


class Assembly:
    """Per-point geometry of a mesh, the material inputs and the element calls."""

    def __init__(self, problem: Problem, finite: bool, kinematics=None):
        self.problem = problem
        self.mesh = problem.mesh
        self.kin = _resolve(kinematics, finite)
        self.finite = self.kin.finite
        self.rotate_inputs = self.kin.rotation == "abaqus"
        mesh = self.mesh
        ne = len(mesh.conn)
        self.ne = ne
        self.n = 8 * ne
        self.ndof = mesh.ndof
        self.Xe = mesh.coords[mesh.conn]                          # (ne, 8, 3)
        self.edofs = (3 * mesh.conn[:, :, None] + np.arange(3)).reshape(ne, 24)
        self.el = nl.C3D8Nlgeom(self.Xe, self.kin.integration if self.finite else "full")
        self.coords = self.el.coordinates().reshape(self.n, 3)
        if getattr(problem, "material_element", None) is not None:
            # the routine caches its point's position: hand it positions inside
            # the author's own element (mechanics unchanged; see mesh.Problem)
            self.coords = material_coordinates(self.coords, mesh, problem.material_element)
        celent = np.array([float(np.prod(x.max(axis=0) - x.min(axis=0))) ** (1.0 / 3.0)
                           for x in self.Xe])
        self.celent = np.repeat(celent, 8)
        self.noel = np.repeat(np.arange(1, ne + 1), 8)
        self.npt = np.tile(np.arange(1, 9), ne)
        self.dNdX = self.el.G[:, :8]

    # ------------------------------------------------------------- helpers
    def _Ue(self, U):
        return np.asarray(U)[self.edofs].reshape(self.ne, 8, 3)

    def _dUe(self, dU):
        """(ndof, m) -> (ne, m, 8, 3)."""
        dU = np.asarray(dU)
        return np.transpose(dU[self.edofs].reshape(self.ne, 8, 3, -1), (0, 3, 1, 2))

    def _assemble_vec(self, fe):
        """(ne, 24[, m]) -> (ndof[, m])."""
        out = np.zeros((self.ndof,) + fe.shape[2:])
        np.add.at(out, self.edofs, fe)
        return out

    def _assemble_mat(self, Ke):
        K = np.zeros((self.ndof, self.ndof))
        for e in range(self.ne):
            K[np.ix_(self.edofs[e], self.edofs[e])] += Ke[e]
        return K

    def grad(self, U) -> np.ndarray:
        """(n, 3, 3) displacement gradient H_iJ = du_i/dX_J at every point."""
        return (self.el.F(self._Ue(U))[:, :8] - _EYE).reshape(self.n, 3, 3)

    def Fbar(self, U) -> np.ndarray:
        return self.el.Fbar(self._Ue(U)).reshape(self.n, 3, 3)

    def dFbar(self, U, dU) -> np.ndarray:
        """(n, m, 3, 3) for dU (ndof, m)."""
        d = self.el.dFbar(self._Ue(U), self._dUe(dU))
        return d.reshape(self.n, d.shape[2], 3, 3)

    # ---------------------------------------------------------- material inputs
    def inputs(self, U_prev, U, incoming) -> Dict[str, np.ndarray]:
        """Everything the UMAT is handed at every point for the increment."""
        n = self.n
        if not self.finite:
            identity = np.tile(_EYE, (n, 1, 1))
            return {"F0": identity, "F1": identity.copy(), "drot": identity.copy(),
                    "dstran": voigt_eng(self.grad(U) - self.grad(U_prev)),
                    "stress_in": incoming.stress, "stran_in": incoming.stran, "kin": None}
        inc = self.el.increment(self._Ue(U_prev), self._Ue(U))
        if self.rotate_inputs:
            R = inc["drot"].reshape(n, 3, 3)
            stress_in = nl.stress_tensor_to_voigt(nl.rotate(R, nl.stress_voigt_to_tensor(incoming.stress)))
            stran_in = nl.strain_tensor_to_voigt(nl.rotate(R, nl.strain_voigt_to_tensor(incoming.stran)))
            drot = R
        else:
            stress_in, stran_in = incoming.stress, incoming.stran
            drot = np.tile(_EYE, (n, 1, 1))
        return {"F0": inc["Fbar0"].reshape(n, 3, 3), "F1": inc["Fbar1"].reshape(n, 3, 3),
                "drot": drot, "dstran": inc["dstran"].reshape(n, 6),
                "stress_in": stress_in, "stran_in": stran_in, "kin": inc}

    def seeds(self, inp, incoming, dU_prev=None, dU=None, dUe0=None, dUe1=None,
              dstress_prev=None, dstran_prev=None, dstate_prev=None) -> Dict[str, Any]:
        """OTI seeds (provider layout) for a set of m directions.

        Displacement directions either global (dU_prev, dU: (ndof, m)) or
        element-local (dUe0, dUe1: (ne, m, 8, 3)); dstress_prev, dstran_prev,
        dstate_prev (n, ., m): derivatives of the UNROTATED outgoing state of the
        previous increment. Returns dstran_dp, dstress_in, stran_dp (n, ., m),
        dstate_in, dF0_dp, dF1_dp (n, m, 3, 3) -- every UMAT input that moves.
        """
        n = self.n
        if dU_prev is not None:
            dUe0 = self._dUe(dU_prev)
        if dU is not None:
            dUe1 = self._dUe(dU)
        moving = [a for a in (dUe0, dUe1) if a is not None]
        m = moving[0].shape[1] if moving else \
            next(a.shape[2] for a in (dstress_prev, dstran_prev, dstate_prev) if a is not None)
        out = {"dstate_in": dstate_prev}
        if not self.finite:
            zero = np.zeros((self.ne, m, 8, 3))
            g0 = _grad_dirs_e(self, zero if dUe0 is None else dUe0)
            g1 = _grad_dirs_e(self, zero if dUe1 is None else dUe1)
            out.update(dstran_dp=np.transpose(voigt_eng(g1 - g0), (0, 2, 1)),
                       dstress_in=dstress_prev, stran_dp=dstran_prev, dF0_dp=None, dF1_dp=None)
            return out
        if not moving:
            d = None
            out.update(dstran_dp=None, dF0_dp=None, dF1_dp=None)
        else:
            d = self.el.increment_direction(inp["kin"], dUe0, dUe1)
            out.update(dstran_dp=np.transpose(d["dstran"].reshape(n, m, 6), (0, 2, 1)),
                       dF0_dp=d["dFbar0"].reshape(n, m, 3, 3),
                       dF1_dp=d["dFbar1"].reshape(n, m, 3, 3))
        if self.rotate_inputs:
            R = inp["kin"]["drot"].reshape(n, 3, 3)
            ddrot = None if d is None else d["ddrot"].reshape(n, m, 3, 3)
            S = nl.stress_voigt_to_tensor(incoming.stress)
            E = nl.strain_voigt_to_tensor(incoming.stran)
            dS = None if dstress_prev is None else nl.stress_voigt_to_tensor(np.transpose(dstress_prev, (0, 2, 1)))
            dE = None if dstran_prev is None else nl.strain_voigt_to_tensor(np.transpose(dstran_prev, (0, 2, 1)))
            if ddrot is None and dS is None:
                out["dstress_in"] = None
            else:
                out["dstress_in"] = np.transpose(nl.stress_tensor_to_voigt(
                    nl.rotate_direction(R, ddrot, S, dS)), (0, 2, 1))
            if ddrot is None and dE is None:
                out["stran_dp"] = None
            else:
                out["stran_dp"] = np.transpose(nl.strain_tensor_to_voigt(
                    nl.rotate_direction(R, ddrot, E, dE)), (0, 2, 1))
        else:
            out["dstress_in"], out["stran_dp"] = dstress_prev, dstran_prev
        return out

    # --------------------------------------------------------------- assembly
    def residual(self, U, stress, *, with_scale: bool = False):
        """Assembled R; with ``with_scale`` also max_dof sum_e |f_e| -- the size of
        what cancels in R, the right yardstick for its round-off."""
        stress = np.asarray(stress).reshape(self.ne, 8, 6)
        if self.finite:
            fe = self.el.force(self._Ue(U), stress)
        else:
            fe = np.array([kernel.element_internal_force_small_strain(self.Xe[e], stress[e])
                           for e in range(self.ne)])
        R = self._assemble_vec(fe)
        if with_scale:
            A = self._assemble_vec(np.abs(fe))
            return R, float(np.max(A)) if A.size else 0.0
        return R

    def residual_scale_entrywise(self, U, stress) -> np.ndarray:
        """Per-DOF magnitude of what cancels in each entry of R: the FD
        round-off model's scale. sum over elements and points of
        w |g_qa|_1 max_k |sigma_q,k| -- the stress's largest component bounds
        the terms that cancel INSIDE the material (a lateral stress that is
        lambda a + (lambda + 2 mu) b ~ 0 carries round-off of size |sigma_xx|),
        which the assembled |f_e| alone does not see."""
        stress = np.asarray(stress).reshape(self.ne, 8, 6)
        smax = np.max(np.abs(stress), axis=2)                          # (ne, 8)
        if self.finite:
            g, wF, _, _ = self.el.geometry(self._Ue(U))
            gq = g[:, :8]
        else:
            gq, wF = self.el.G[:, :8], self.el.w0
        mag = np.einsum("eq,eqaj,eq->ea", wF, np.abs(gq), smax)        # (ne, 8)
        fe = np.repeat(mag[:, :, None], 3, axis=2).reshape(self.ne, 24)
        return self._assemble_vec(fe)

    def tangent_ra(self, U, stress, ddsdde) -> np.ndarray:
        """DDSDDE-based tangent (small: B^T D B; finite: Abaqus-style
        Bbar^T c Bbar + K_geo with c from the Jaumann/Kirchhoff reading of DDSDDE)."""
        stress = np.asarray(stress).reshape(self.ne, 8, 6)
        D = np.asarray(ddsdde).reshape(self.ne, 8, 6, 6)
        if self.finite:
            return self._assemble_mat(self.el.tangent_abaqus(self._Ue(U), stress, D))
        K = np.zeros((self.ndof, self.ndof))
        for e in range(self.ne):
            k = kernel.element_tangent(self.Xe[e], None, D[e], mode="small")
            K[np.ix_(self.edofs[e], self.edofs[e])] += k
        return K

    def tangent_exact(self, U, stress, dsig_du=None, ddsdde=None) -> np.ndarray:
        """Exact Jacobian of ``residual``. small: ``ddsdde`` = d sigma / d DSTRAN;
        finite: ``dsig_du`` (ne, 8, 24, 6) = total d sigma_q / d u_(e, a) along
        the 24 element DOF directions (every UMAT input that moves with u_n)."""
        if not self.finite:
            return self.tangent_ra(U, stress, ddsdde)
        stress = np.asarray(stress).reshape(self.ne, 8, 6)
        Ue = self._Ue(U)
        E = self.el.unit_directions(self.ne)
        Ke = self.el.dforce_fixed_stress(Ue, stress, E) + self.el.force(Ue, np.asarray(dsig_du))
        return self._assemble_mat(Ke)

    def geometric_fixed_sigma(self, U, stress, dU) -> np.ndarray:
        """(ndof, m): d R / d u . dU at fixed Cauchy stress (zero for small strain)."""
        dU = np.asarray(dU)
        if not self.finite:
            return np.zeros(dU.shape)
        stress = np.asarray(stress).reshape(self.ne, 8, 6)
        return self._assemble_vec(self.el.dforce_fixed_stress(self._Ue(U), stress, self._dUe(dU)))

    def dR_dp(self, U, dstress_dp) -> np.ndarray:
        """sum Bbar^T (d sigma/d p) detJ w at fixed geometry: (ndof, npar)."""
        ds = np.asarray(dstress_dp).reshape(self.ne, 8, 6, -1)
        if self.finite:
            return self._assemble_vec(self.el.force(self._Ue(U), np.transpose(ds, (0, 1, 3, 2))))
        npar = ds.shape[-1]
        out = np.zeros((self.ndof, npar))
        for e in range(self.ne):
            for j in range(npar):
                out[self.edofs[e], j] += element_dR_dp(self.Xe[e], ds[e, :, :, j], None, mode="small")
        return out


# ---------------------------------------------------------------------------
class Frame:
    """x = Q (X + V): displacement u = Q (X + V) - X; residual back with Q^T."""

    def __init__(self, coords: np.ndarray, Q: Optional[np.ndarray]):
        self.X = np.asarray(coords)
        self.Q = None if Q is None else np.asarray(Q, dtype=float)
        nn = len(self.X)
        self.T = None if self.Q is None else np.kron(np.eye(nn), self.Q)

    def actual(self, V):
        if self.Q is None:
            return np.asarray(V).copy()
        x = (self.X + np.asarray(V).reshape(-1, 3)) @ self.Q.T
        return (x - self.X).reshape(-1)

    def to_frame(self, R):
        """Q^T per node for vectors (ndof,) or stacks (ndof, m)."""
        return R if self.T is None else self.T.T @ R

    def from_frame(self, dV):
        return dV if self.T is None else self.T @ dV

    def matrix(self, K):
        return K if self.T is None else self.T.T @ K @ self.T


@dataclass
class IncrementState:
    stress: np.ndarray
    state: np.ndarray
    stran: np.ndarray


@dataclass
class HistoryRun:
    """One solve of the whole load path at one PROPS vector.

    ``U`` actual displacements (configuration x = X + U); ``V`` the unknowns in
    the problem frame (= U without a frame rotation); ``reactions`` in the
    problem frame; ``du_dp`` = dV/dp.
    """

    props: np.ndarray
    material: str                         # 'oti' | 'regular'
    kinematics: str = ""
    U: List[np.ndarray] = field(default_factory=list)
    V: List[np.ndarray] = field(default_factory=list)
    reactions: List[np.ndarray] = field(default_factory=list)
    incoming: List[IncrementState] = field(default_factory=list)
    outgoing: List[IncrementState] = field(default_factory=list)
    newton: List[List[float]] = field(default_factory=list)
    newton_matrix: List[str] = field(default_factory=list)
    #: per increment: 'tolerance' (|R_free| <= tolerance * force scale) or
    #: 'roundoff_floor' (|R_free| <= the round-off floor of R, see run_history)
    converged_by: List[str] = field(default_factory=list)
    roundoff_floor: List[float] = field(default_factory=list)
    #: per increment: index into ``newton`` from which the iteration matrix was
    #: the exact one (0 unless newton_matrix='ra_then_exact' switched later)
    exact_from: List[int] = field(default_factory=list)
    backtracks: List[int] = field(default_factory=list)
    #: the ``newton_matrix`` argument the run was made with
    newton_matrix_arg: str = "auto"
    #: the increments actually solved (problem.schedule() unless cut back)
    schedule: List[dict] = field(default_factory=list)
    #: every cut-back: increment, depth, load-factor span, why
    cutbacks: List[dict] = field(default_factory=list)
    force_scale: List[float] = field(default_factory=list)
    du_dp: List[np.ndarray] = field(default_factory=list)
    dstress_dp: List[np.ndarray] = field(default_factory=list)
    dstate_dp: List[np.ndarray] = field(default_factory=list)
    local_dR_dp: List[np.ndarray] = field(default_factory=list)
    total_dR_dp: List[np.ndarray] = field(default_factory=list)
    sensitivity_equilibrium: List[float] = field(default_factory=list)
    primal_parity: List[float] = field(default_factory=list)
    ddsdde: List[np.ndarray] = field(default_factory=list)
    outputs: List[Dict[str, np.ndarray]] = field(default_factory=list)   # raw material outputs
    seconds: float = 0.0

    def qoi(self, problem: Problem) -> np.ndarray:
        rows = []
        for V, R in zip(self.V, self.reactions):
            row = []
            for q in problem.qois:
                row.append(float(np.sum(R[q["dofs"]])) if q["kind"] == "reaction" else float(V[q["dof"]]))
            rows.append(row)
        return np.array(rows)

    def qoi_dp(self, problem: Problem) -> np.ndarray:
        rows = []
        for du, dR in zip(self.du_dp, self.total_dR_dp):
            row = []
            for q in problem.qois:
                row.append(dR[q["dofs"]].sum(axis=0) if q["kind"] == "reaction" else du[q["dof"]])
            rows.append(row)
        return np.array(rows)


class Driver:
    """Material calls for one problem, one provider."""

    def __init__(self, provider: CorpusProvider, assembly: Assembly):
        self.provider = provider
        self.asm = assembly

    def _kw(self, inc: dict, kinc: int):
        return dict(time=np.array([inc["time_start"], inc["time_start"]]), dtime=inc["dtime"],
                    coords=self.asm.coords, celent=self.asm.celent, noel=self.asm.noel,
                    npt=self.asm.npt, kstep=1, kinc=kinc)

    def regular(self, props, U_prev, U, incoming, inc, kinc, inp=None):
        inp = inp or self.asm.inputs(U_prev, U, incoming)
        return self.provider.regular(props, inp["stress_in"], incoming.state, inp["stran_in"],
                                     inp["dstran"], inp["F0"], inp["F1"], drot=inp["drot"],
                                     **self._kw(inc, kinc))

    def oti(self, props, U_prev, U, incoming, inc, kinc, inp=None, **seeds):
        inp = inp or self.asm.inputs(U_prev, U, incoming)
        seeds = {k: v for k, v in seeds.items() if v is not None}
        return self.provider.total(props, inp["stress_in"], incoming.state, inp["stran_in"],
                                   inp["dstran"], inp["F0"], inp["F1"], drot=inp["drot"],
                                   **self._kw(inc, kinc), **seeds)

    def dsig_du(self, props, U_prev, U, incoming, inc, kinc, base, inp=None):
        """(ne, 8, 24, 6) total d sigma / d u along the 24 element-DOF directions
        (DFGRD1, DSTRAN, DROT-driven pre-rotation of STRESS and STRAN) by OTI
        seeds, all elements at once; the unit-PROPS response is removed."""
        asm = self.asm
        ne, npar = asm.ne, self.provider.nparam
        inp = inp or asm.inputs(U_prev, U, incoming)
        E = asm.el.unit_directions(ne)
        out_d = np.zeros((ne, 8, 24, 6))
        for start in range(0, 24, npar):
            chunk = list(range(start, min(24, start + npar)))
            dUe1 = np.zeros((ne, npar, 8, 3))
            dUe1[:, :len(chunk)] = E[:, chunk]
            s = asm.seeds(inp, incoming, dUe1=dUe1)
            out = self.oti(props, U_prev, U, incoming, inc, kinc, inp=inp,
                           dstran_dp=s["dstran_dp"], dF0_dp=s["dF0_dp"], dF1_dp=s["dF1_dp"],
                           dstress_in=s["dstress_in"], stran_dp=s["stran_dp"])
            d = (out["dstress_dp"] - base["dstress_dp"]).reshape(ne, 8, 6, npar)
            out_d[:, :, chunk] = np.transpose(d[..., :len(chunk)], (0, 1, 3, 2))
        return out_d

    def exact_tangent(self, props, U_prev, U, incoming, inc, kinc, base, inp=None):
        """Exact dR/du (actual frame) at the trial state; ``base`` an OTI output
        with no seeds at the same state."""
        if self.asm.finite:
            return self.asm.tangent_exact(U, base["stress"], dsig_du=self.dsig_du(
                props, U_prev, U, incoming, inc, kinc, base, inp))
        return self.asm.tangent_exact(U, base["stress"], ddsdde=base["ddsdde"])


def run_history(provider: CorpusProvider, problem: Problem, props, *, material: str = "oti",
                sensitivities: bool = False, newton_matrix: str = "auto",
                tolerance: float = 1e-12, max_iterations: int = 40,
                kinematics=None, keep_outputs: bool = False, cutbacks: int = 0,
                schedule: Optional[List[dict]] = None) -> HistoryRun:
    """Solve the load path; optionally carry the total parameter sensitivities.

    Convergence: ``|R_free|_inf <= tolerance * scale`` (``scale`` = max_dof
    sum_e |f_e|), OR ``|R_free|_inf <= ROUNDOFF_FLOOR_FACTOR * eps * |K_ff|_inf * L``:
    the residual cannot be evaluated more accurately than that. Under finite
    strain the routine sees F = I + grad u, whose entries carry an absolute
    round-off of eps, i.e. a nodal position uncertainty of eps * L (L the mesh
    extent); the iteration matrix maps it to a residual uncertainty of
    |K_ff|_inf * eps * L. Small strain: L is max |u| (no identity is added).
    A nearly incompressible material (bulk/shear ~ 5e3) sits above the 1e-12
    relative tolerance at that floor: without it, Newton stalls there and fails
    with a relative residual of ~4e-12 (Jeff97, B2). The floor uses the last
    iteration matrix assembled (none yet: the floor is not used); which
    criterion ended each increment is recorded in ``converged_by``. A floor
    acceptance needs one confirming full Newton step: it must not halve the
    residual (else the iterate was passing through the floor and Newton goes
    on); the better of the two states is kept and both appear in ``newton``.

    ``newton_matrix='ra_then_exact'``: the DDSDDE matrix until the relative
    free residual is below 1e-6 (or 8 iterations), then the exact one -- the
    fallback when the exact matrix alone does not converge from the increment's
    start (nearly incompressible material between clamps: Jeff97, B8). The
    matrix only steers the iteration; the converged state solves R = 0.

    Increment cut-back (``cutbacks`` > 0), as Abaqus/Standard's automatic
    incrementation does: an increment whose Newton iteration fails, or whose
    routine cannot be evaluated or asks for it (PNEWDT < 1), is split into two
    halves in load factor AND time, down to ``cutbacks`` halvings. The increments
    actually solved are ``run.schedule``; a re-solve that must be comparable
    increment by increment (finite differences at p +/- h) passes that list as
    ``schedule`` with ``cutbacks=0``, so both sides use the same increments.

    ``newton_matrix='ra_stalled_exact'``: the DDSDDE matrix, the exact one from
    the 8th iteration on -- what ``'auto'`` does for the original, for either build.

    ``material='regular'`` drives the ORIGINAL routine (reference runs). Its
    Newton iteration matrix is the DDSDDE tangent of the original's own
    DDSDDE; if that has not converged in 8 iterations it switches to the OTI
    exact tangent at the same PROPS -- an iteration matrix only: the converged
    state solves R_original = 0 whichever matrix led there; recorded.
    """
    asm = Assembly(problem, provider.finite, kinematics)
    drv = Driver(provider, asm)
    props = np.asarray(props, dtype=float)
    run = HistoryRun(props=props.copy(), material=material, kinematics=asm.kin.name)
    started = _time.time()
    n, nt, ns, npar = asm.n, provider.ntens, provider.nstatv, provider.nparam
    free, cons = problem.free, problem.constrained
    amp = problem.amplitudes()
    rotations = getattr(problem, "rotations", None)
    incoming = IncrementState(np.zeros((n, nt)), np.zeros((n, ns)), np.zeros((n, nt)))
    initial = provider.case.extra.get("initial_statev") or []
    if initial:
        incoming.state[:, :len(initial)] = np.asarray(initial, dtype=float)[:ns]
    if provider.case.extra.get("initial_state_from_user_subroutine"):
        if not getattr(provider, "has_sdvini", False):
            raise MaterialCallError("the deck initialises STATEV through SDVINI and the source "
                                    "defines none")
        incoming.state = provider.sdvini(incoming.state, asm.coords, asm.noel, asm.npt)
    V_prev = np.zeros(asm.ndof)
    frame_prev = Frame(asm.mesh.coords, None)
    d_in = dict(dstress=np.zeros((n, nt, npar)), dstate=np.zeros((n, ns, npar)),
                dstran=np.zeros((n, nt, npar)), du=np.zeros((asm.ndof, npar)))
    extent = float(np.max(np.ptp(asm.mesh.coords, axis=0))) if asm.finite else 0.0
    k_inf = None                         # |K_ff|_inf of the last iteration matrix
    run.newton_matrix_arg = newton_matrix
    plan = [dict(inc) for inc in (schedule if schedule is not None else problem.schedule())]
    if cutbacks and rotations is not None:
        raise ValueError("increment cut-back with a superposed rotation is not supported")
    while plan:
        inc = plan.pop(0)
        index = len(run.schedule)
        kinc = index + 1
        frame = Frame(asm.mesh.coords, None if rotations is None else rotations[index])
        U_prev = frame_prev.actual(V_prev)
        try:
            (U, V, R, out, inp, scale, floor, history, matrix, exact_from, backtracks,
             converged_by, k_inf) = _solve_increment(
                asm, drv, material, props, problem, inc, kinc, frame, U_prev, V_prev, incoming,
                amp, free, cons, newton_matrix, tolerance, max_iterations, extent, k_inf)
        except (NewtonFailed, MaterialCallError, np.linalg.LinAlgError) as error:
            if int(inc.get("cutback_depth", 0)) >= cutbacks:
                raise
            plan[:0] = _halves(inc)
            run.cutbacks.append({"increment": kinc, "depth": int(inc.get("cutback_depth", 0)) + 1,
                                 "lambda": [inc["lambda_start"], inc["lambda"]],
                                 "reason": "%s: %s" % (type(error).__name__, error)})
            continue
        run.schedule.append(inc)
        run.backtracks.append(backtracks)
        run.converged_by.append(converged_by)
        run.exact_from.append(exact_from if matrix != "ra" else len(history))
        run.roundoff_floor.append(floor)
        run.force_scale.append(scale)
        run.newton.append(history)
        run.newton_matrix.append(matrix)
        run.U.append(U.copy())
        run.V.append(V.copy())
        run.reactions.append(R.copy())
        run.incoming.append(incoming)
        run.ddsdde.append(out["ddsdde"].copy())
        if keep_outputs:
            run.outputs.append({k: np.array(v, copy=True) for k, v in out.items()})
        new_state = IncrementState(out["stress"].copy(), out["state"].copy(),
                                   inp["stran_in"] + inp["dstran"])
        run.outgoing.append(new_state)
        if material == "oti":
            check = drv.regular(props, U_prev, U, incoming, inc, kinc, inp=inp)
            sc = max(float(np.max(np.abs(out["stress"]))), 1e-300)
            run.primal_parity.append(float(np.max(np.abs(check["stress"] - out["stress"]))) / sc)
        if sensitivities:
            _sensitivity_step(run, asm, drv, props, U_prev, U, incoming, inc, kinc, d_in, free,
                              inp, frame)
        incoming = new_state
        V_prev, frame_prev = V, frame
    run.seconds = _time.time() - started
    return run



def _halves(inc: dict) -> List[dict]:
    """Two increments covering ``inc``: half the load-factor change, half the time."""
    mid = 0.5 * (inc["lambda_start"] + inc["lambda"])
    half = 0.5 * inc["dtime"]
    depth = int(inc.get("cutback_depth", 0)) + 1
    first = dict(inc, **{"lambda": mid, "dtime": half, "cutback_depth": depth,
                         "segment_end": False})
    second = dict(inc, **{"lambda_start": mid, "dtime": half, "cutback_depth": depth,
                          "time_start": inc["time_start"] + half})
    return [first, second]


def _solve_increment(asm, drv, material, props, problem, inc, kinc, frame, U_prev, V_prev,
                     incoming, amp, free, cons, newton_matrix, tolerance, max_iterations,
                     extent, k_inf):
    """Newton iteration of one increment (see run_history); nothing is committed."""
    V = V_prev.copy()
    V[cons] = inc["lambda"] * amp
    history, matrix, exact_from = [], newton_matrix, 0
    if matrix == "auto":
        matrix = "exact" if material == "oti" else "ra"
    elif matrix in ("ra_then_exact", "ra_stalled_exact"):
        matrix = "ra"            # the exact matrix once close / stalled (see below)
    out = inp = None
    last = None                      # (V, delta, rfree) of the last accepted Newton step
    alpha, backtracks, iteration = 1.0, 0, 0
    candidate = None                 # the state first found at the round-off floor
    while True:
        U = frame.actual(V)
        try:
            inp = asm.inputs(U_prev, U, incoming)
            if material == "oti":
                out = drv.oti(props, U_prev, U, incoming, inc, kinc, inp=inp)
            else:
                out = drv.regular(props, U_prev, U, incoming, inc, kinc, inp=inp)
            R_act, scale = asm.residual(U, out["stress"], with_scale=True)
            R = frame.to_frame(R_act)
            scale = max(scale, 1e-300)
            rfree = float(np.max(np.abs(R[free]))) if free.size else 0.0
        except (MaterialCallError, nl.InvertedElement) as error:
            if candidate is not None:
                # the confirming step left what the routine can evaluate: the
                # floor state stands (no further progress is possible from it)
                U, V, R, out, inp, scale, rfree, floor = candidate
                history.append(rfree / scale)
                converged_by = "roundoff_floor"
                break
            # a trial state the routine (or the element: det F <= 0) cannot
            # evaluate: a rejected step, halved
            if last is None or alpha <= 1.0 / 256.0:
                if isinstance(error, nl.InvertedElement):
                    raise MaterialCallError(str(error))
                raise
            alpha *= 0.5
            backtracks += 1
            V = last[0].copy()
            V[free] += alpha * last[1]
            continue
        length = extent if asm.finite else float(np.max(np.abs(V))) if V.size else 0.0
        floor = None if k_inf is None else ROUNDOFF_FLOOR_FACTOR * _EPS * k_inf * length
        at_tolerance = rfree <= tolerance * scale or rfree == 0.0
        at_floor = floor is not None and rfree <= floor
        if at_tolerance:
            history.append(rfree / scale)
            converged_by = "tolerance"
            break
        if candidate is not None:
            # the confirming step (Vera, B8 review): the floor is accepted only
            # if one further full Newton step does not halve the residual, i.e.
            # the residual has stopped falling; the better of the two states is
            # kept. A residual that still falls was passing through the floor.
            if rfree >= 0.5 * candidate[6] or at_floor:
                if candidate[6] < rfree:
                    history.append(rfree / scale)
                    U, V, R, out, inp, scale, rfree, floor = candidate
                history.append(rfree / scale)
                converged_by = "roundoff_floor"
                break
            candidate = None
        if at_floor:
            candidate = (U, V.copy(), R, out, inp, scale, rfree, floor)
        elif last is not None and rfree > last[2] and alpha > 1.0 / 256.0:
            alpha *= 0.5
            backtracks += 1
            V = last[0].copy()
            V[free] += alpha * last[1]
            continue
        history.append(rfree / scale)
        if iteration >= max_iterations and candidate is None:   # the confirming step is extra
            raise NewtonFailed("increment %d: free residual %.3e of %.3e after %d iterations "
                               "(%s matrix, %d backtracks)"
                               % (kinc, rfree, scale, iteration, matrix, backtracks))
        iteration += 1
        if matrix == "ra":
            K = asm.tangent_ra(U, out["stress"], out["ddsdde"])
        else:
            base = out if material == "oti" else drv.oti(props, U_prev, U, incoming, inc,
                                                          kinc, inp=inp)
            K = drv.exact_tangent(props, U_prev, U, incoming, inc, kinc, base, inp)
        K = frame.matrix(K)
        Kff = K[np.ix_(free, free)]
        k_inf = float(np.max(np.sum(np.abs(Kff), axis=1))) if free.size else 0.0
        delta = np.linalg.solve(Kff, -R[free])
        last, alpha = (V.copy(), delta, rfree), 1.0
        V[free] += delta
        if matrix == "ra" and iteration >= 8 and (
                newton_matrix == "ra_stalled_exact"
                or (material == "regular" and newton_matrix in ("auto", "ra_then_exact"))):
            matrix = "exact+ra_stalled"
        elif newton_matrix == "ra_then_exact" and matrix == "ra" and \
                (rfree <= 1e-6 * scale or iteration >= 8):
            matrix, exact_from = "ra_then_exact", len(history)
    return (U, V, R, out, inp, scale, floor, history, matrix, exact_from, backtracks,
            converged_by, k_inf)

def _sensitivity_step(run, asm, drv, props, U_prev, U, incoming, inc, kinc, d_in, free, inp,
                      frame):
    """Local dR/dp, the total du/dp solve, and the carried derivatives."""
    npar = drv.provider.nparam
    # local: every input held fixed, unit PROPS seed only
    local = drv.oti(props, U_prev, U, incoming, inc, kinc, inp=inp)
    run.local_dR_dp.append(asm.dR_dp(U, local["dstress_dp"]))
    # phase A: u_n held fixed, history derivatives carried (through Fbar0, DSTRAN,
    # DROT and the pre-rotated STRESS/STRAN)
    sA = asm.seeds(inp, incoming, dU_prev=d_in["du"], dstress_prev=d_in["dstress"],
                   dstran_prev=d_in["dstran"], dstate_prev=d_in["dstate"])
    A = drv.oti(props, U_prev, U, incoming, inc, kinc, inp=inp, **sA)
    dRdp_u = frame.to_frame(asm.dR_dp(U, A["dstress_dp"]))
    K = frame.matrix(drv.exact_tangent(props, U_prev, U, incoming, inc, kinc, local, inp))
    dV = np.zeros((asm.ndof, npar))
    if free.size:
        dV[free] = np.linalg.solve(K[np.ix_(free, free)], -dRdp_u[free])
    du = frame.from_frame(dV)
    # full seeds: the converged update's total derivatives
    sF = asm.seeds(inp, incoming, dU_prev=d_in["du"], dU=du, dstress_prev=d_in["dstress"],
                   dstran_prev=d_in["dstran"], dstate_prev=d_in["dstate"])
    full = drv.oti(props, U_prev, U, incoming, inc, kinc, inp=inp, **sF)
    total = asm.dR_dp(U, full["dstress_dp"]) + asm.geometric_fixed_sigma(U, full["stress"], du)
    total = frame.to_frame(total)
    scale = max(float(np.max(np.abs(dRdp_u))), float(np.max(np.abs(total))), 1e-300)
    run.sensitivity_equilibrium.append(float(np.max(np.abs(total[free]))) / scale if free.size else 0.0)
    run.du_dp.append(dV)
    run.dstress_dp.append(full["dstress_dp"])
    dstate = full["dstate_dp"]
    undefined = [i - 1 for i in (drv.provider.case.extra or {}).get("undefined_statev", ())
                 if 1 <= i <= dstate.shape[1]]
    if undefined:
        # declared undefined in the original and (hidden-state probe) never read:
        # nothing downstream depends on them, so neither does their derivative
        dstate = np.array(dstate, copy=True)
        dstate[:, undefined] = 0.0
    run.dstate_dp.append(dstate)
    run.total_dR_dp.append(total)
    d_in["dstress"] = full["dstress_dp"]
    d_in["dstate"] = dstate
    d_in["dstran"] = (0.0 if sF["stran_dp"] is None else sF["stran_dp"]) + sF["dstran_dp"]
    d_in["du"] = du


def _grad_dirs_e(asm: Assembly, dUe) -> np.ndarray:
    """(n, m, 3, 3) gradient directions of element directions dUe (ne, m, 8, 3)."""
    d = np.einsum("emai,eqaj->eqmij", dUe, asm.dNdX)
    return d.reshape(asm.n, d.shape[2], 3, 3)
