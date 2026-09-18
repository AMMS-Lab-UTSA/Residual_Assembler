"""Output manager: propagate solution sensitivities du/dp to the response
sensitivities dq/dp the collaborator actually asked for.

For any response q(u, p),   dq/dp = @q/@p + (@q/@u) . du/dp.

Supported responses (v1):
  displacement  q = u at (node, component)   -> @q/@p = 0, @q/@u = e_dof
  reaction      q = reaction at a fixed DOF  -> @q/@p = R_,p[dof], @q/@u = K[dof,:]
  stress        q = sigma_component (or von Mises) at (element, ip)
                    @q/@u via D.B at that IP, @q/@p from the explicit d sigma/d p
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping

import numpy as np

from ..formulations import c3d8_kernel as kern


class OutputError(Exception):
    pass


@dataclass
class OutputContext:
    model: Any
    dof_manager: Any
    params: List[str]
    u: np.ndarray                     # converged displacement (ndof,)
    K: np.ndarray                     # global tangent (ndof, ndof)
    R_da: np.ndarray                  # residual derivatives (ndof, m)
    du_dp: np.ndarray                 # displacement sensitivities (ndof, m)
    f_ext: np.ndarray                 # external load (ndof,)
    f_int: np.ndarray                 # assembled internal force at converged u (ndof,)
    tangent_fields: Mapping[int, np.ndarray]           # {eid -> (8,6,6)}
    stress_derivative_fields: Mapping[str, Mapping[int, np.ndarray]]  # {p->{eid->(8,6)}}
    stress_replay: Mapping[int, np.ndarray]            # {eid -> (8,6)}


def _dof_index(dm, node: int, comp: int) -> int:
    return dm.node_dofs(int(node))[int(comp) - 1]


def _vm(sig: np.ndarray) -> float:
    s11, s22, s33, s12, s13, s23 = sig
    return float(np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                         + 3.0 * (s12 * s12 + s13 * s13 + s23 * s23)))


def _dvm_dsig(sig: np.ndarray) -> np.ndarray:
    s11, s22, s33, s12, s13, s23 = sig
    vm = _vm(sig)
    if vm < 1e-30:
        return np.zeros(6)
    return np.array([2 * s11 - s22 - s33, 2 * s22 - s11 - s33, 2 * s33 - s11 - s22,
                     6 * s12, 6 * s13, 6 * s23]) / (2.0 * vm)


def compute_outputs(requests: List[Mapping[str, Any]], ctx: OutputContext
                    ) -> List[Dict[str, Any]]:
    """Return one result per request: {request, value, dq_dp:{param:val}}."""
    out = []
    m = len(ctx.params)
    for req in requests:
        kind = str(req.get("type", "")).lower()
        if kind == "displacement":
            d = _dof_index(ctx.dof_manager, req["node"], req["dof"])
            value = float(ctx.u[d])
            grad = ctx.du_dp[d, :]
        elif kind == "reaction":
            d = _dof_index(ctx.dof_manager, req["node"], req["dof"])
            # reaction = internal force - external at a constrained dof;
            # dReaction/dp = R_,p[d] + K[d,:] . du/dp   (chain rule, @q/@u = K[d,:])
            value = float(ctx.f_int[d] - ctx.f_ext[d])
            grad = ctx.R_da[d, :] + ctx.K[d, :] @ ctx.du_dp
        elif kind == "stress":
            eid, ip = int(req["element"]), int(req["ip"])
            el = ctx.model.elements[eid]
            Xe = ctx.model.coords_of(el.connectivity)
            edofs = np.asarray(ctx.dof_manager.element_dofs(
                el.connectivity, ("UX", "UY", "UZ")), int)
            B = kern.b_matrix_reference(Xe, kern.ABAQUS_C3D8_GAUSS.points[ip])[0]
            D = ctx.tangent_fields[eid][ip]                  # (6,6)
            dsig_du = D @ B                                  # (6,24) = d sigma/d u_e
            du_e = ctx.du_dp[edofs, :]                       # (24, m)
            implicit = dsig_du @ du_e                        # (6, m)
            explicit = np.zeros((6, m))
            for c, p in enumerate(ctx.params):
                explicit[:, c] = ctx.stress_derivative_fields[p][eid][ip]
            dsig_total = explicit + implicit                # (6, m)
            sig0 = ctx.stress_replay[eid][ip]
            if req.get("von_mises"):
                value = _vm(sig0)
                grad = _dvm_dsig(sig0) @ dsig_total
            else:
                comp = int(req["component"])
                value = float(sig0[comp])
                grad = dsig_total[comp, :]
        else:
            raise OutputError("unknown output type %r (want displacement|reaction|stress)"
                              % req.get("type"))
        out.append({
            "request": dict(req),
            "value": None if value is None else float(value),
            "dq_dp": {ctx.params[j]: float(grad[j]) for j in range(m)},
        })
    return out
