"""Framework-level verification of the formulation-agnostic assembler.

Proves that the generic core (Model + DofManager + Assembler + formulations)
reproduces the ALREADY-VERIFIED c3d8_kernel results, for both supported offline
modes, WITHOUT the assembler knowing anything about the physics:

  Test A (Mode 1, stress-driven): assembler + StressDrivenC3D8 formulation ==
          the verified kernel assemble_global_internal_force, to machine eps.
  Test B (Mode 2, material-update): assembler + SolidC3D8SmallStrain +
          IsotropicElastic material == kernel small-strain assembly of the same
          elastic stress; and its Level-4 finite-difference tangent passes.
  Test C (Level 1): zero-field residual vanishes.

No Abaqus needed. Run: python tests/framework/test_assembler.py
"""

import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from residual_core.io import abaqus_inp_parser
from residual_core.core.model import from_abaqus
from residual_core.core.dof_manager import DofManager
from residual_core.core.state_manager import StateManager
from residual_core.core.assembler import Assembler
from residual_core.core import verification
from residual_core.formulations import default_formulations
from residual_core.formulations import c3d8_kernel as kern
from residual_core.materials import IsotropicElastic
from residual_core.materials.base import MaterialBinding

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from external_sources import require_external_file  # noqa: E402

INP = os.path.join(_ROOT, "sources", "permissive", "ngrilli_Oxford_Crystal_Plasticity",
                   "ExampleInputFiles", "HCPnoTwin", "Compression111.inp")


@pytest.fixture(scope="module")
def am():
    """The parsed Abaqus model shared by the tests (pytest injects it)."""
    require_external_file(
        INP, "a real C3D8 crystal-plasticity mesh to assemble the global residual against")
    return abaqus_inp_parser.parse_inp(INP)


def _script_run(fn, *args):
    """Script-mode runner: a test passes unless it raises AssertionError."""
    try:
        fn(*args)
    except AssertionError as exc:
        print("  FAIL: %s" % exc)
        return False
    except pytest.skip.Exception as exc:
        print("  SKIP: %s" % exc)
    return True


def _mesh(model):
    node_ids = sorted(model.nodes.keys())
    coords = np.array([model.nodes[n] for n in node_ids], float)
    conn = [(e.eid, list(e.connectivity)) for e in model.elements.values()
            if e.etype.upper() == "C3D8"]
    return node_ids, coords, conn


def test_A_stress_driven(am):
    """Generic assembler + StressDrivenC3D8 == verified kernel assembly."""
    model = from_abaqus(am, formulation_policy=lambda et: "stress_driven_c3d8")
    dm = DofManager(model.nodes.keys())
    asm = Assembler(model, dm, default_formulations(), StateManager(8, 0))

    rng = np.random.default_rng(1)
    S0 = rng.uniform(-50, 50, (3, 3)); S0 = 0.5 * (S0 + S0.T)
    grads = [0.5 * (G + G.T) for G in rng.uniform(-8, 8, (3, 3, 3))]
    pts = kern.ABAQUS_C3D8_GAUSS.points

    def sig(x):
        M = S0 + grads[0]*x[0] + grads[1]*x[1] + grads[2]*x[2]
        return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])

    sigma_all = {}
    for e in model.elements.values():
        Xe = model.coords_of(e.connectivity)
        sigma_all[e.eid] = np.array([sig(kern.shape_functions(pts[k]) @ Xe) for k in range(8)])

    # reference: verified kernel path
    node_ids, coords, conn = _mesh(model)
    F_ref, _ = kern.assemble_global_internal_force(
        node_ids, coords, conn, np.zeros(dm.ndof), sigma_all, mode="small")
    # generic path (no external loads -> R == F_int)
    R, _, diag = asm.assemble(np.zeros(dm.ndof),
                              fields={"stress_ip": sigma_all}, options={"config": "small"})
    rel = np.linalg.norm(R - F_ref) / max(np.linalg.norm(F_ref), 1e-30)
    ok = rel < 1e-12 and diag["elements"] == 125
    print("  A stress-driven vs kernel: rel=%.3e  elems=%d  -> %s"
          % (rel, diag["elements"], "PASS" if ok else "FAIL"))
    assert ok, ("stress-driven vs kernel: rel=%.3e (need < 1e-12), elems=%d (need 125)"
                % (rel, diag["elements"]))


def test_B_material_and_tangent(am):
    """Generic assembler + elastic material == kernel elastic assembly; Level-4 FD."""
    E, nu = 200000.0, 0.3
    bindings = {name: MaterialBinding(IsotropicElastic(), [E, nu], 0, name)
                for name in am.materials}
    model = from_abaqus(am, formulation_policy=lambda et: "solid_c3d8_small_strain",
                        material_bindings=bindings)
    dm = DofManager(model.nodes.keys())
    asm = Assembler(model, dm, default_formulations(), StateManager(8, 0))

    rng = np.random.default_rng(3)
    U = rng.uniform(-1e-3, 1e-3, dm.ndof)      # small strains

    # reference: kernel small-strain assembly of the elastic stress at this U
    D = kern.isotropic_D(E, nu)
    node_ids, coords, conn = _mesh(model)
    id_row = {n: i for i, n in enumerate(node_ids)}
    pts = kern.ABAQUS_C3D8_GAUSS.points
    sigma_all = {}
    for eid, c in conn:
        Xe = coords[[id_row[n] for n in c]]
        ue = np.concatenate([U[dm.node_dofs(n)] for n in c])
        sigma_all[eid] = np.array([D @ (kern.b_matrix_reference(Xe, pts[k])[0] @ ue)
                                   for k in range(8)])
    F_ref, _ = kern.assemble_global_internal_force(node_ids, coords, conn, U, sigma_all, mode="small")
    R, K, _ = asm.assemble(U, compute_tangent=True)
    rel_force = np.linalg.norm(R - F_ref) / max(np.linalg.norm(F_ref), 1e-30)

    # Level-4 finite-difference tangent on a DOF subset (fast, exact for linear elastic)
    subset = list(range(300, 312))
    abs_e, rel_fro, _, _, _ = verification.finite_difference_tangent(
        asm, U, dof_subset=subset)
    ok = rel_force < 1e-12 and rel_fro < 1e-6
    print("  B material path vs kernel: rel_force=%.3e   Level-4 FD tangent rel=%.3e  -> %s"
          % (rel_force, rel_fro, "PASS" if ok else "FAIL"))
    assert ok, ("material path vs kernel: rel_force=%.3e (need < 1e-12), "
                "Level-4 FD tangent rel=%.3e (need < 1e-6)" % (rel_force, rel_fro))


def test_C_zero_field(am):
    bindings = {name: MaterialBinding(IsotropicElastic(), [200000.0, 0.3], 0, name)
                for name in am.materials}
    model = from_abaqus(am, formulation_policy=lambda et: "solid_c3d8_small_strain",
                        material_bindings=bindings)
    dm = DofManager(model.nodes.keys())
    asm = Assembler(model, dm, default_formulations(), StateManager(8, 0))
    r_inf, _ = verification.zero_field_residual(asm)
    ok = r_inf < 1e-9
    print("  C zero-field residual (U=0): max|R|=%.3e  -> %s" % (r_inf, "PASS" if ok else "FAIL"))
    assert ok, "zero-field residual: max|R|=%.3e (need < 1e-9)" % r_inf


def main():
    print("Framework assembler verification (formulation-agnostic core)")
    print("parsing:", os.path.relpath(INP, _ROOT))
    model = abaqus_inp_parser.parse_inp(INP)
    results = [_script_run(fn, model) for fn in (test_A_stress_driven,
                                                 test_B_material_and_tangent,
                                                 test_C_zero_field)]
    ok = all(results)
    print("OVERALL: %s" % ("ALL PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
