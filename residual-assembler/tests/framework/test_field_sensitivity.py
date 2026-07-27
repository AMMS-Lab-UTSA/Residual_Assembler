"""M1 acceptance: field-driven residual-method parameter sensitivities.

Proves the Python side can take ARBITRARY per-integration-point derivative
fields (DDSDDE and dsigma/da) and produce CORRECT global displacement
sensitivities du/da, verified against full central finite differences of the
elastic solve -- on a genuinely multi-element mesh, for two parameters solved
simultaneously.

Nothing here needs Abaqus: the "exported" derivative fields are synthesized in
closed form from linear elasticity (dD/dE = D/E, dsigma/dE = sigma/E; dD/dnu by
a material-level central difference), then handed to the SAME engine a real UMAT
export would feed. The engine never sees the closed forms -- only the numbers.

Checks
------
1. Two-element C3D8, params E and nu: residual-method du/dE, du/dnu match the
   full-solve central FD to < 1e-6 (relative, inf-norm). Both parameter columns
   are solved together; prescribed-DOF sensitivities are exactly zero.
2. Single-element parity: the engine reproduces the analytic single-cube
   reference in examples/residual_sensitivity_c3d8/sensitivity_engine.py.
3. SDV-layout round-trip: fields_from_statev(parse_sdv_layout(...)) reproduces
   the direct-field result.
4. Honest failures: unknown/misspelled parameter, malformed field, missing
   element, overlapping SDV ranges all raise instead of returning zeros.

Run:  pytest tests/framework/test_field_sensitivity.py
  or:  python tests/framework/test_field_sensitivity.py
"""

import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

try:
    import pytest
except ImportError:                        # script mode without pytest installed
    class _Raises:
        def __init__(self, exc):
            self.exc = exc
        def __enter__(self):
            return self
        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError("did not raise %s" % self.exc.__name__)
            return issubclass(et, self.exc)

    class _PytestShim:
        @staticmethod
        def raises(exc):
            return _Raises(exc)
    pytest = _PytestShim()

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.core.model import Element, Model
from residual_core.core.dof_manager import DofManager
from residual_core.core.constraints import partition
from residual_core.formulations import c3d8_kernel as kern
from residual_core.core.field_sensitivity import (
    solve_field_sensitivities, fields_from_statev, parse_sdv_layout,
    FieldSensitivityError, N_IP, DDSDDE_LEN)

E0, NU0 = 210000.0, 0.3


# --------------------------------------------------------------------------- #
# minimal boundary object (duck-typed for core.constraints)
# --------------------------------------------------------------------------- #
@dataclass
class _BC:
    target: int
    dof_start: int = 0
    dof_end: int = 0
    value: float = 0.0
    kind: str = "value"
    amplitude: Optional[str] = None


def _fix(target, comp, value=0.0):
    """Prescribe a single translational component (1=x,2=y,3=z) at a node."""
    return _BC(target=target, dof_start=comp, dof_end=comp, value=value, kind="value")


def _encastre(target):
    return _BC(target=target, kind="ENCASTRE")


# --------------------------------------------------------------------------- #
# two stacked unit cubes (12 nodes, 2 C3D8), base clamped, top pulled in +z
# --------------------------------------------------------------------------- #
def _two_element_model():
    nodes = {
        1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
        5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
        9: (0.0, 0.0, 2.0), 10: (1.0, 0.0, 2.0), 11: (1.0, 1.0, 2.0), 12: (0.0, 1.0, 2.0),
    }
    elements = {
        1: Element(1, "C3D8", [1, 2, 3, 4, 5, 6, 7, 8]),
        2: Element(2, "C3D8", [5, 6, 7, 8, 9, 10, 11, 12]),
    }
    m = Model(nodes=nodes, elements=elements)
    m.element_formulation = {1: "solid_c3d8_small_strain", 2: "solid_c3d8_small_strain"}
    m.element_material = {1: "MAT", 2: "MAT"}
    m.boundaries = [_encastre(n) for n in (1, 2, 3, 4)]     # clamp base face z=0
    return m


def _load_vector(model, dm, total_fz=1000.0):
    """Total +z force spread over the four top nodes (z=2)."""
    f = np.zeros(dm.ndof)
    for n in (9, 10, 11, 12):
        f[dm.node_dofs(n)[2]] += total_fz / 4.0
    return f


# --------------------------------------------------------------------------- #
# an independent full elastic solve (the FD "truth" generator)
# --------------------------------------------------------------------------- #
def _assemble_K(model, dm, E, nu):
    D = kern.isotropic_D(E, nu)
    K = np.zeros((dm.ndof, dm.ndof))
    for el in model.elements.values():
        Xe = model.coords_of(el.connectivity)
        edofs = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
        k_e = kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=D, mode="small")
        K[np.ix_(edofs, edofs)] += k_e
    return K


def _full_solve(model, dm, E, nu, fext):
    """Linear-elastic global solve with the model's Dirichlet BCs (prescribed
    values are 0 here). Returns the global displacement (ndof,)."""
    K = _assemble_K(model, dm, E, nu)
    free_mask, pres_idx, pres = partition(model, dm)
    u = np.zeros(dm.ndof)
    for d, v in pres.items():
        u[d] = v
    Kff = K[np.ix_(free_mask, free_mask)]
    rhs = fext[free_mask] - K[np.ix_(free_mask, pres_idx)] @ u[pres_idx] \
        if pres_idx.size else fext[free_mask]
    u[free_mask] = np.linalg.solve(Kff, rhs)
    return u


def _synthetic_fields(model, dm, u_star, E, nu, hnu=1e-6):
    """Build per-IP DDSDDE and dsigma/dE, dsigma/dnu at the converged state.

    For linear elasticity  sigma = D(E,nu) . eps ,  eps = B u:
        dsigma/dE|_u  = (dD/dE) eps = (D/E) eps = sigma/E     (D linear in E)
        dsigma/dnu|_u = (dD/dnu) eps                          (dD/dnu by central diff)
    These are exactly the residual-method ingredients ( fixed-u derivatives ).
    """
    D = kern.isotropic_D(E, nu)
    dDdnu = (kern.isotropic_D(E, nu + hnu) - kern.isotropic_D(E, nu - hnu)) / (2 * hnu)
    pts = kern.ABAQUS_C3D8_GAUSS.points

    tangent_fields = {}
    dsig_dE = {}
    dsig_dnu = {}
    for el in model.elements.values():
        Xe = model.coords_of(el.connectivity)
        edofs = np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
        ue = u_star[edofs]
        dE = np.zeros((N_IP, 6))
        dnu = np.zeros((N_IP, 6))
        for k in range(N_IP):
            B = kern.b_matrix_reference(Xe, pts[k])[0]
            eps = B @ ue
            sig = D @ eps
            dE[k] = sig / E
            dnu[k] = dDdnu @ eps
        tangent_fields[el.eid] = np.array([D] * N_IP)
        dsig_dE[el.eid] = dE
        dsig_dnu[el.eid] = dnu
    return tangent_fields, {"E": dsig_dE, "nu": dsig_dnu}


def _rel_inf(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30)


# --------------------------------------------------------------------------- #
# 1) the headline acceptance test
# --------------------------------------------------------------------------- #
def test_two_element_E_nu_vs_finite_difference():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)

    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)

    res = solve_field_sensitivities(
        model=model, solution=u_star,
        tangent_fields=tangent_fields, stress_derivative_fields=dsig,
        parameters=["E", "nu"])

    # both columns solved together
    assert res.displacement_sensitivities.shape == (dm.ndof, 2)
    assert res.parameters == ["E", "nu"]

    # prescribed-DOF sensitivities are exactly zero
    assert np.all(res.displacement_sensitivities[res.prescribed_idx] == 0.0)
    assert res.prescribed_idx.size == 12   # 4 clamped base nodes x 3 dofs

    # full central finite differences of the elastic solve
    hE = 1e-6 * E0
    du_dE_fd = (_full_solve(model, dm, E0 + hE, NU0, fext)
                - _full_solve(model, dm, E0 - hE, NU0, fext)) / (2 * hE)
    hnu = 1e-6
    du_dnu_fd = (_full_solve(model, dm, E0, NU0 + hnu, fext)
                 - _full_solve(model, dm, E0, NU0 - hnu, fext)) / (2 * hnu)

    rel_E = _rel_inf(res.du_da("E"), du_dE_fd)
    rel_nu = _rel_inf(res.du_da("nu"), du_dnu_fd)
    print("  two-element  du/dE  rel(inf)=%.3e   du/dnu rel(inf)=%.3e   (tol 1e-6)"
          % (rel_E, rel_nu))
    # sensitivities are genuinely non-trivial (not a zero-vs-zero pass)
    assert np.max(np.abs(du_dE_fd)) > 1e-12
    assert np.max(np.abs(du_dnu_fd)) > 1e-9
    assert rel_E < 1e-6, "du/dE rel=%.3e" % rel_E
    assert rel_nu < 1e-6, "du/dnu rel=%.3e" % rel_nu


# --------------------------------------------------------------------------- #
# 2) single-element parity with the analytic example engine
# --------------------------------------------------------------------------- #
def test_single_element_matches_example_engine():
    from examples.residual_sensitivity_c3d8 import sensitivity_engine as eng

    # one unit cube, node ids 1..8 with the engine's coords; global dof(node-1,comp)
    # matches the engine's dof(node0,comp) exactly, so results are directly comparable.
    nodes = {i + 1: tuple(eng.XE[i]) for i in range(8)}
    model = Model(nodes=nodes,
                  elements={1: Element(1, "C3D8", [1, 2, 3, 4, 5, 6, 7, 8])})
    model.element_formulation = {1: "solid_c3d8_small_strain"}
    model.element_material = {1: "MAT"}

    # engine Case A boundary + load: symmetry BCs, total Fx=100 on the x=1 face
    fixed = eng.symmetry_bcs()                       # {dof0: 0.0}
    model.boundaries = [_fix(d // 3 + 1, d % 3 + 1, v) for d, v in fixed.items()]
    fext = np.zeros(24)
    for n in eng.X1:
        fext[eng.dof(n, 0)] = 100.0 / 4

    ref = eng.sensitivity_dE(E0, fixed, fext)        # analytic du/dE (24,)

    # synthesize the SAME fields the engine's overloaded UMAT stand-in would export
    dm = DofManager(model.nodes.keys())
    u_star = _full_solve(model, dm, E0, eng.NU, fext)
    D = kern.isotropic_D(E0, eng.NU)
    pts = kern.ABAQUS_C3D8_GAUSS.points
    Xe = model.coords_of([1, 2, 3, 4, 5, 6, 7, 8])
    edofs = np.asarray(dm.element_dofs([1, 2, 3, 4, 5, 6, 7, 8], ("UX", "UY", "UZ")), int)
    sig = np.array([D @ (kern.b_matrix_reference(Xe, pts[k])[0] @ u_star[edofs])
                    for k in range(8)])
    res = solve_field_sensitivities(
        model=model, solution=u_star,
        tangent_fields={1: np.array([D] * 8)},
        stress_derivative_fields={"E": {1: sig / E0}}, parameters=["E"])

    rel = _rel_inf(res.du_da("E"), ref["dudE"])
    print("  single-element vs example engine  du/dE rel(inf)=%.3e   (tol 1e-9)" % rel)
    assert rel < 1e-9, "single-element parity rel=%.3e" % rel


# --------------------------------------------------------------------------- #
# 3) SDV-layout round-trip
# --------------------------------------------------------------------------- #
def test_sdv_layout_roundtrip():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)
    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)

    # pack a raw statev array per element: [ DDSDDE(36) | dE(6) | dnu(6) ] per IP
    layout = {"indexing": "abaqus_1_based", "ddsdde": [1, 36],
              "parameters": {"E": [37, 42], "nu": [43, 48]}}
    statev = {}
    for eid in (1, 2):
        rows = []
        for ip in range(8):
            rows.append(list(tangent_fields[eid][ip].ravel())
                        + list(dsig["E"][eid][ip]) + list(dsig["nu"][eid][ip]))
        statev[eid] = rows

    tf2, sf2 = fields_from_statev(statev, layout, parameters=["E", "nu"])
    res_direct = solve_field_sensitivities(
        model=model, tangent_fields=tangent_fields,
        stress_derivative_fields=dsig, parameters=["E", "nu"])
    res_statev = solve_field_sensitivities(
        model=model, tangent_fields=tf2,
        stress_derivative_fields=sf2, parameters=["E", "nu"])

    rel = _rel_inf(res_statev.displacement_sensitivities,
                   res_direct.displacement_sensitivities)
    print("  SDV round-trip  du/da rel(inf)=%.3e   (tol 1e-12)" % rel)
    assert rel < 1e-12
    # parsed slices are correct: SDV37..42 -> python statev[36:42]
    parsed = parse_sdv_layout(layout)
    assert parsed["parameters"]["E"] == slice(36, 42)
    assert parsed["parameters"]["nu"] == slice(42, 48)
    assert parsed["ddsdde"] == slice(0, 36)


# --------------------------------------------------------------------------- #
# 4) honest failure conditions
# --------------------------------------------------------------------------- #
def test_unknown_parameter_raises():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)
    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)
    with pytest.raises(FieldSensitivityError):
        solve_field_sensitivities(model=model, tangent_fields=tangent_fields,
                                  stress_derivative_fields=dsig,
                                  parameters=["E", "YOUNGS_TYPO"])


def test_malformed_tangent_raises():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)
    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)
    tangent_fields[1] = tangent_fields[1][:, :5, :]        # (8,5,6): not 6x6
    with pytest.raises(FieldSensitivityError):
        solve_field_sensitivities(model=model, tangent_fields=tangent_fields,
                                  stress_derivative_fields=dsig, parameters=["E"])


def test_missing_element_raises():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)
    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)
    del dsig["E"][2]                                       # element 2 missing
    with pytest.raises(FieldSensitivityError):
        solve_field_sensitivities(model=model, tangent_fields=tangent_fields,
                                  stress_derivative_fields=dsig, parameters=["E"])


def test_unsupported_element_type_raises():
    model = _two_element_model()
    dm = DofManager(model.nodes.keys())
    fext = _load_vector(model, dm)
    u_star = _full_solve(model, dm, E0, NU0, fext)
    tangent_fields, dsig = _synthetic_fields(model, dm, u_star, E0, NU0)
    model.elements[3] = Element(3, "C3D4", [1, 2, 3, 5])   # a tet sneaks in
    with pytest.raises(FieldSensitivityError):
        solve_field_sensitivities(model=model, tangent_fields=tangent_fields,
                                  stress_derivative_fields=dsig, parameters=["E"])


def test_overlapping_sdv_ranges_raise():
    with pytest.raises(FieldSensitivityError):
        parse_sdv_layout({"ddsdde": [1, 36], "parameters": {"E": [36, 41]}})  # 36 overlaps


def test_misspelled_parameter_in_statev_raises():
    layout = {"ddsdde": [1, 36], "parameters": {"E": [37, 42]}}
    statev = {1: [[0.0] * 42 for _ in range(8)]}
    with pytest.raises(FieldSensitivityError):
        fields_from_statev(statev, layout, parameters=["nu"])   # nu not in layout


# --------------------------------------------------------------------------- #
# script mode
# --------------------------------------------------------------------------- #
def main():
    tests = [
        test_two_element_E_nu_vs_finite_difference,
        test_single_element_matches_example_engine,
        test_sdv_layout_roundtrip,
        test_unknown_parameter_raises,
        test_malformed_tangent_raises,
        test_missing_element_raises,
        test_unsupported_element_type_raises,
        test_overlapping_sdv_ranges_raise,
        test_misspelled_parameter_in_statev_raises,
    ]
    ok = True
    for t in tests:
        try:
            t()
            print("  PASS: %s" % t.__name__)
        except Exception as exc:            # noqa: BLE001 - script-mode reporting
            ok = False
            print("  FAIL: %s -> %s" % (t.__name__, exc))
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
