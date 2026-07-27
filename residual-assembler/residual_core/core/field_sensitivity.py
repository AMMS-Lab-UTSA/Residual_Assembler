"""Field-driven residual-method parameter sensitivities (Milestone M1).

Given a converged **small-strain** state and, at every integration point, the
material tangent ``D = dsigma/deps`` (the Abaqus ``DDSDDE``) and the stress
parameter-derivatives ``dsigma/da_i``, this module assembles the global tangent

    K = sum_e sum_q  B_eq^T D_eq B_eq  J_eq w_q

and, for every parameter ``a_i``, the residual parameter-derivative column

    R_,a_i = sum_e sum_q  B_eq^T (dsigma_eq/da_i)  J_eq w_q

then solves all parameters simultaneously on the free partition

    K_ff [ U_,a1 ... U_,am ] = -[ R_,a1 ... R_,am ]

for the global displacement sensitivities ``du/da_i``. Prescribed-DOF
sensitivities are zero (parameter-independent boundary conditions -- the
supported v1 scope); this is enforced by only solving the free block.

Relationship to the rest of the framework
------------------------------------------
This is the production generalization of the single-element demo in
``examples/residual_sensitivity_c3d8``. It reuses the audited generic pieces and
does NOT re-implement any element math or DOF bookkeeping:

  * per-IP element tangent / internal force  -> ``formulations.c3d8_kernel``
  * global DOF numbering (node-major)        -> ``core.dof_manager.DofManager``
  * free/prescribed partition (incl. XSYMM/ENCASTRE) -> ``core.constraints.partition``

Nothing here evaluates a material model: it only *consumes* exported derivative
fields. The only material-specific inputs are ``D`` (per IP) and ``dsigma/da``
(per IP, per parameter). For a real UMAT these come from the OTIS-transformed
UMAT's SDVs; for verification they may be synthetic (elastic closed form).

Honesty gate
------------
Every requested parameter must be backed by a real derivative field, every
element in the model must be present in the fields, and every field must have
the right shape. A missing/misspelled parameter or a missing element raises
:class:`FieldSensitivityError` rather than silently producing an all-zero (and
therefore plausible-looking but wrong) sensitivity.

Conventions (binding -- see ``docs`` and ``c3d8_kernel``)
---------------------------------------------------------
* Voigt order (stress/strain): ``(11, 22, 33, 12, 13, 23)`` -- Abaqus order.
* ``DDSDDE`` packs a 6x6 matrix, row-major -> 36 values.
* Integration points: 8 per C3D8, in the kernel's Abaqus Gauss order
  (``ABAQUS_C3D8_GAUSS``); the field arrays MUST be in that IP order.
* SDV layout metadata uses **1-based inclusive Abaqus SDV numbers**; they are
  converted to 0-based Python slices internally (``SDV_k <-> statev[k-1]``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..formulations import c3d8_kernel as _k
from .constraints import partition as _partition
from .dof_manager import DofManager

# --------------------------------------------------------------------------- #
# constants
# --------------------------------------------------------------------------- #
N_IP = 8                                  # C3D8 full integration
N_VOIGT = 6
DDSDDE_LEN = N_VOIGT * N_VOIGT             # 36
SUPPORTED_ETYPES = ("C3D8",)
SOLID_DOF_TYPES = ("UX", "UY", "UZ")

# Volumetric integration of the C3D8 element:
#   'full'              -- 2x2x2 Gauss on every term (the default engine).
#   'selective_reduced' -- mean-dilatation (B-bar) volumetric term; this is what
#                          Abaqus's fully-integrated C3D8 uses. Required to match
#                          fields exported from an Abaqus C3D8 analysis.
INTEGRATION_MODES = ("full", "selective_reduced")


class FieldSensitivityError(ValueError):
    """A user-facing problem with the derivative fields, the SDV layout, or the
    model that must stop the run rather than yield a plausible wrong number."""


# --------------------------------------------------------------------------- #
# result container
# --------------------------------------------------------------------------- #
@dataclass
class FieldSensitivityResult:
    """Output of :func:`solve_field_sensitivities`.

    Attributes
    ----------
    parameters : ordered parameter names, matching the columns below.
    K : (ndof, ndof) assembled global tangent.
    residual_derivatives : (ndof, m) columns R_,a_i (the sensitivity RHS, un-negated).
    displacement_sensitivities : (ndof, m) columns du/da_i (prescribed rows are 0).
    free_mask : (ndof,) bool, True on free DOFs.
    prescribed_idx : (k,) int indices of prescribed DOFs.
    ndof : total number of global DOFs.
    dof_manager : the DofManager used (node-major UX,UY,UZ).
    solution : the converged displacement passed in (or None); unused for
        small-strain assembly, retained for provenance / output chain-rules.
    """
    parameters: List[str]
    K: np.ndarray
    residual_derivatives: np.ndarray
    displacement_sensitivities: np.ndarray
    free_mask: np.ndarray
    prescribed_idx: np.ndarray
    ndof: int
    dof_manager: Any = None
    solution: Optional[np.ndarray] = None
    diagnostics: Dict[str, Any] = field(default_factory=dict)

    def du_da(self, name: str) -> np.ndarray:
        """Displacement-sensitivity column du/da for one parameter name."""
        try:
            j = self.parameters.index(name)
        except ValueError:
            raise FieldSensitivityError(
                "no sensitivity for parameter %r; solved parameters are %r"
                % (name, self.parameters))
        return self.displacement_sensitivities[:, j]

    def R_da(self, name: str) -> np.ndarray:
        """Residual-derivative column R_,a for one parameter name."""
        j = self.parameters.index(name)
        return self.residual_derivatives[:, j]


# --------------------------------------------------------------------------- #
# field validation helpers
# --------------------------------------------------------------------------- #
def _as_tangent_ip(raw, eid) -> np.ndarray:
    """Validate and return a per-IP tangent field of shape (8, 6, 6)."""
    if raw is None:
        raise FieldSensitivityError("element %r has no tangent (DDSDDE) field" % (eid,))
    D = np.asarray(raw, dtype=float)
    if D.shape != (N_IP, N_VOIGT, N_VOIGT):
        raise FieldSensitivityError(
            "element %r tangent field must have shape (%d, 6, 6), got %r; each of "
            "the %d integration points needs a 6x6 DDSDDE"
            % (eid, N_IP, tuple(D.shape), N_IP))
    return D


def _as_dsigma_ip(raw, eid, pname) -> np.ndarray:
    """Validate and return a per-IP stress-derivative field of shape (8, 6)."""
    if raw is None:
        raise FieldSensitivityError(
            "element %r has no dsigma/d%s field" % (eid, pname))
    d = np.asarray(raw, dtype=float)
    if d.shape != (N_IP, N_VOIGT):
        raise FieldSensitivityError(
            "element %r dsigma/d%s field must have shape (%d, 6), got %r; each of "
            "the %d integration points needs a length-6 Voigt stress derivative"
            % (eid, pname, N_IP, tuple(d.shape), N_IP))
    return d


def _resolve_parameters(stress_derivative_fields, parameters) -> List[str]:
    """Ordered parameter list, validated against the supplied fields."""
    available = list(stress_derivative_fields.keys())
    if parameters is None:
        params = list(available)
        if not params:
            raise FieldSensitivityError(
                "no parameters given and stress_derivative_fields is empty; "
                "supply at least one dsigma/da field")
        return params
    params = list(parameters)
    if not params:
        raise FieldSensitivityError("parameters list is empty")
    missing = [p for p in params if p not in stress_derivative_fields]
    if missing:
        raise FieldSensitivityError(
            "requested parameter(s) %r have no dsigma/da field (available: %r); a "
            "misspelled parameter must never become a silent zero sensitivity"
            % (missing, available))
    return params


def _assembled_c3d8_elements(model) -> List[Tuple[int, Any]]:
    """Return [(eid, element)] to assemble, rejecting unsupported element types.

    v1 supports C3D8 only. An unsupported element type is a hard error, not a
    silent skip -- silently dropping an element would produce a wrong global
    system with no warning.
    """
    unsupported: Dict[str, List[int]] = {}
    c3d8: List[Tuple[int, Any]] = []
    for eid, el in sorted(model.elements.items()):
        et = el.etype.upper()
        if et in SUPPORTED_ETYPES:
            c3d8.append((eid, el))
        else:
            unsupported.setdefault(et, []).append(eid)
    if unsupported:
        detail = ", ".join("%s (%d elem)" % (t, len(v)) for t, v in unsupported.items())
        raise FieldSensitivityError(
            "unsupported element type(s) present: %s. The field-driven M1 solver "
            "supports %s only; refusing to assemble a partial system."
            % (detail, "/".join(SUPPORTED_ETYPES)))
    if not c3d8:
        raise FieldSensitivityError("model has no C3D8 elements to assemble")
    return c3d8


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #
def solve_field_sensitivities(
    model,
    tangent_fields: Mapping[int, Any],
    stress_derivative_fields: Mapping[str, Mapping[int, Any]],
    solution: Optional[Sequence[float]] = None,
    parameters: Optional[Sequence[str]] = None,
    dof_manager: Optional[DofManager] = None,
    mode: str = "small",
    integration: str = "full",
) -> FieldSensitivityResult:
    """Assemble ``K`` and ``R_,a_i`` from per-IP derivative fields and solve for
    the global displacement sensitivities ``du/da_i``.

    Parameters
    ----------
    model : a ``core.model.Model`` (mesh + bindings + boundaries).
    tangent_fields : ``{eid -> (8, 6, 6)}`` per-IP DDSDDE (dsigma/deps). Builds K.
    stress_derivative_fields : ``{param_name -> {eid -> (8, 6)}}`` per-IP
        dsigma/da at the converged state. Builds one R_,a column per parameter.
    solution : converged nodal displacement (ndof,). Optional and unused for
        ``mode='small'`` (small-strain K and R_,a do not depend on u); retained
        on the result for provenance and downstream output chain-rules.
    parameters : optional explicit ordered parameter list. If given, every name
        must have a field (a misspelling raises). Defaults to the field keys.
    dof_manager : optional prebuilt DofManager; by default a node-major
        UX,UY,UZ manager is built from the model's nodes.
    mode : 'small' (v1). 'finite' is intentionally not implemented here.
    integration : 'full' (default; 2x2x2 Gauss) or 'selective_reduced'
        (mean-dilatation B-bar volumetric term -- required to match fields
        exported from an Abaqus C3D8 analysis, whose fully-integrated C3D8 uses
        selective reduced integration).

    Returns
    -------
    :class:`FieldSensitivityResult`.
    """
    if mode != "small":
        raise NotImplementedError(
            "field_sensitivity M1 supports mode='small' only; finite-strain "
            "assembly (spatial B, geometric term, objective-rate DDSDDE) is a "
            "later milestone. Got mode=%r." % (mode,))
    if integration not in INTEGRATION_MODES:
        raise FieldSensitivityError(
            "integration must be one of %r, got %r" % (INTEGRATION_MODES, integration))

    params = _resolve_parameters(stress_derivative_fields, parameters)
    elements = _assembled_c3d8_elements(model)

    dm = dof_manager if dof_manager is not None else DofManager(model.nodes.keys())
    ndof = dm.ndof
    m = len(params)

    K = np.zeros((ndof, ndof), dtype=float)
    R_da = np.zeros((ndof, m), dtype=float)

    # every parameter must cover every assembled element -- check up front so a
    # missing element is reported as such, not as a mysterious zero row.
    eids = {eid for eid, _ in elements}
    for pname in params:
        pf = stress_derivative_fields[pname]
        missing = sorted(eids - set(pf.keys()))
        if missing:
            raise FieldSensitivityError(
                "dsigma/d%s is missing %d element(s) present in the model: %r%s"
                % (pname, len(missing), missing[:8],
                   " ..." if len(missing) > 8 else ""))
    missing_tan = sorted(eids - set(tangent_fields.keys()))
    if missing_tan:
        raise FieldSensitivityError(
            "tangent (DDSDDE) field is missing %d element(s) present in the "
            "model: %r%s" % (len(missing_tan), missing_tan[:8],
                             " ..." if len(missing_tan) > 8 else ""))

    bbar = integration == "selective_reduced"
    for eid, el in elements:
        Xe = model.coords_of(el.connectivity)                 # (8,3), conn order
        edofs = np.asarray(dm.element_dofs(el.connectivity, SOLID_DOF_TYPES), int)
        D_ip = _as_tangent_ip(tangent_fields.get(eid), eid)   # (8,6,6)

        if bbar:
            k_e = _k.element_tangent_bbar(Xe, D_ip)
        else:
            k_e = _k.element_tangent(Xe=Xe, Ue=None, Dmat_ip=D_ip, mode="small")
        K[np.ix_(edofs, edofs)] += k_e

        for col, pname in enumerate(params):
            dsig = _as_dsigma_ip(stress_derivative_fields[pname].get(eid), eid, pname)
            if bbar:
                r_e = _k.element_internal_force_bbar(Xe, dsig)      # (24,)
            else:
                r_e = _k.element_internal_force_small_strain(Xe, dsig)
            R_da[edofs, col] += r_e

    free_mask, pres_idx, _pres_vals = _partition(model, dm)
    if not free_mask.any():
        raise FieldSensitivityError(
            "every DOF is prescribed; there is no free partition to solve on")

    U = np.zeros((ndof, m), dtype=float)
    Kff = K[np.ix_(free_mask, free_mask)]
    try:
        # solve all parameters together: one factorization, m right-hand sides
        U[free_mask, :] = np.linalg.solve(Kff, -R_da[free_mask, :])
    except np.linalg.LinAlgError as exc:
        raise FieldSensitivityError(
            "global tangent is singular on the free partition (%s); the model has "
            "rigid-body / zero-stiffness modes -- check that the boundary "
            "conditions fully constrain it (a converged, properly-constrained FE "
            "problem yields a non-singular K)." % exc)

    diagnostics = {
        "n_elements": len(elements),
        "n_parameters": m,
        "n_free": int(free_mask.sum()),
        "n_prescribed": int(pres_idx.size),
        "mode": mode,
        "integration": integration,
        "tangent_source": "exported_ddsdde",
    }
    return FieldSensitivityResult(
        parameters=params, K=K, residual_derivatives=R_da,
        displacement_sensitivities=U, free_mask=free_mask,
        prescribed_idx=pres_idx, ndof=ndof, dof_manager=dm,
        solution=None if solution is None else np.asarray(solution, float).ravel(),
        diagnostics=diagnostics)


# --------------------------------------------------------------------------- #
# SDV layout parsing: raw per-IP statev  ->  tangent / stress-derivative fields
# --------------------------------------------------------------------------- #
def parse_sdv_layout(sdv_layout: Mapping[str, Any]) -> Dict[str, Any]:
    """Convert a metadata SDV layout into 0-based Python slices, with checks.

    Input (Abaqus 1-based inclusive ranges)::

        {
          "indexing": "abaqus_1_based",   # optional, this is the only mode
          "ddsdde": [1, 36],              # -> statev[0:36], row-major 6x6
          "parameters": {"E": [37, 42], "nu": [43, 48]},
        }

    Returns ``{"ddsdde": slice, "parameters": {name: slice}, "n_min": int}``
    where ``n_min`` is the smallest statev length that satisfies the layout.

    Raises :class:`FieldSensitivityError` on a bad indexing mode, wrong DDSDDE
    length, a parameter range that is not length 6, or any overlap between the
    DDSDDE block and a parameter block (or between two parameter blocks).
    """
    indexing = str(sdv_layout.get("indexing", "abaqus_1_based"))
    if indexing != "abaqus_1_based":
        raise FieldSensitivityError(
            "unsupported SDV indexing %r; only 'abaqus_1_based' is supported "
            "(1-based inclusive Abaqus SDV numbers)." % indexing)

    def _to_slice(rng, label, expect_len):
        try:
            lo, hi = int(rng[0]), int(rng[1])
        except (TypeError, ValueError, IndexError):
            raise FieldSensitivityError(
                "SDV range for %s must be [start, end] (1-based inclusive), got %r"
                % (label, rng))
        if lo < 1 or hi < lo:
            raise FieldSensitivityError(
                "SDV range for %s is invalid: [%r, %r] (need 1 <= start <= end)"
                % (label, rng[0], rng[1]))
        length = hi - lo + 1
        if expect_len is not None and length != expect_len:
            raise FieldSensitivityError(
                "SDV range for %s spans %d slots but %d were expected "
                "(range %r)" % (label, length, expect_len, rng))
        return slice(lo - 1, hi), (lo, hi)

    if "ddsdde" not in sdv_layout:
        raise FieldSensitivityError("SDV layout is missing the 'ddsdde' range")
    ddsdde_slice, ddsdde_lohi = _to_slice(sdv_layout["ddsdde"], "ddsdde", DDSDDE_LEN)

    occupied: Dict[int, str] = {}

    def _mark(lohi, label):
        lo, hi = lohi
        for s in range(lo, hi + 1):
            if s in occupied:
                raise FieldSensitivityError(
                    "SDV slot %d is claimed by both %r and %r; ranges must not "
                    "overlap" % (s, occupied[s], label))
            occupied[s] = label

    _mark(ddsdde_lohi, "ddsdde")

    param_slices: Dict[str, slice] = {}
    for name, rng in dict(sdv_layout.get("parameters", {})).items():
        sl, lohi = _to_slice(rng, "parameter %r" % name, N_VOIGT)
        _mark(lohi, "parameter %r" % name)
        param_slices[name] = sl

    n_min = max(occupied) if occupied else 0
    return {"ddsdde": ddsdde_slice, "parameters": param_slices, "n_min": n_min}


def fields_from_statev(
    statev: Mapping[Any, Any],
    sdv_layout: Mapping[str, Any],
    parameters: Optional[Sequence[str]] = None,
) -> Tuple[Dict[int, np.ndarray], Dict[str, Dict[int, np.ndarray]]]:
    """Split a raw per-IP statev export into tangent and stress-derivative fields.

    Parameters
    ----------
    statev : ``{eid -> (nIP, nSDV)}`` per-IP state-variable vectors (as produced
        by ``io.abaqus_odb_export.sdv_fields`` -> element -> IP -> [SDVs]).
    sdv_layout : the metadata layout (see :func:`parse_sdv_layout`).
    parameters : optional subset/order of parameter names to extract. Each must
        be present in the layout (a misspelling raises).

    Returns
    -------
    ``(tangent_fields, stress_derivative_fields)`` ready for
    :func:`solve_field_sensitivities`.
    """
    layout = parse_sdv_layout(sdv_layout)
    all_params = list(layout["parameters"].keys())
    if parameters is None:
        want = all_params
    else:
        want = list(parameters)
        unknown = [p for p in want if p not in layout["parameters"]]
        if unknown:
            raise FieldSensitivityError(
                "parameter(s) %r are not in the SDV layout (layout has %r); a "
                "misspelled parameter must never become a silent zero sensitivity"
                % (unknown, all_params))

    ddsdde_slice = layout["ddsdde"]
    n_min = layout["n_min"]

    tangent_fields: Dict[int, np.ndarray] = {}
    stress_fields: Dict[str, Dict[int, np.ndarray]] = {p: {} for p in want}

    for eid, per_ip in statev.items():
        try:
            eid_i = int(eid)
        except (TypeError, ValueError):
            eid_i = eid
        arr = np.asarray(per_ip, dtype=float)
        if arr.ndim != 2 or arr.shape[0] != N_IP:
            raise FieldSensitivityError(
                "element %r statev must have shape (%d, nSDV), got %r"
                % (eid, N_IP, tuple(arr.shape)))
        if arr.shape[1] < n_min:
            raise FieldSensitivityError(
                "element %r statev has %d SDVs but the layout needs at least %d"
                % (eid, arr.shape[1], n_min))
        D = np.array([arr[ip, ddsdde_slice].reshape(N_VOIGT, N_VOIGT)
                      for ip in range(N_IP)])
        tangent_fields[eid_i] = D
        for p in want:
            sl = layout["parameters"][p]
            stress_fields[p][eid_i] = arr[:, sl].copy()

    return tangent_fields, stress_fields
