"""A verified UMAT's own deck, through the whole assembler, end to end.

Every other check in this repository starts somewhere inside it: a stress is
handed to a kernel, or a tangent to an element. This one starts where the
evidence starts. A fixture carries the numbers Abaqus computed AND the deck
that produced them, so the chain

    deck text -> parser -> steps and boundary conditions -> DOF numbering ->
    element connectivity -> integration points -> B matrix -> strain

can be walked here and compared against a fact: the strain Abaqus handed the
UMAT at that integration point, recorded in the fixture. If any link is
wrong the reconstructed strain is wrong, and the disagreement is a number
rather than an argument.

Then the numbers go the other way. The verified stress and tangent drive a
MULTI-ELEMENT model through ``Assembler``, ``DofManager`` and the registered
C3D8 formulation -- the real objects, not a kernel call -- and the assembled
``dR/du`` is differenced against the assembled ``R`` over a plateau of step
sizes, with the componentwise errors recorded. A single element cannot test a
scatter: every local row is a global row and a wrong map is invisible. Eight
elements sharing nodes can.

What this does NOT establish: that the UMAT is right. Abaqus and the
verification pipeline settled that upstream. What is on trial here is
everything this repository does with the answer.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from residual_core.core import constraints  # noqa: E402
from residual_core.core.assembler import Assembler  # noqa: E402
from residual_core.core.deck_replay import replay_deck  # noqa: E402
from residual_core.core.dof_manager import DofManager  # noqa: E402
from residual_core.core.finite_difference import (  # noqa: E402
    DEFAULT_STEPS, centred_jacobian, componentwise_error, sweep_steps)
from residual_core.core.model import Element, Model  # noqa: E402
from residual_core.core.state_manager import StateManager  # noqa: E402
from residual_core.formulations.solid_c3d8_small_strain import (  # noqa: E402
    SolidC3D8SmallStrain)
from residual_core.materials.fixture_playback import binding_for  # noqa: E402
from verified_fixtures import hex_fixtures  # noqa: E402

#: How closely a differenced column has to sit to the assembled one inside the
#: plateau. Measured on the committed set: the global small-strain assembly
#: reaches 1e-9 to 1e-11 across eight elements, so this leaves two decades of
#: headroom and fails on a wrong assembly rather than on rounding.
AGREES = 1e-7


def ids(fixtures):
    return [f.path.name.split("--")[0] for f in fixtures]


# --------------------------------------------------------------------------- #
# the deck, walked forwards
# --------------------------------------------------------------------------- #
@pytest.mark.integration
@pytest.mark.parametrize("fixture", hex_fixtures(), ids=ids(hex_fixtures()))
def test_the_strain_reconstructed_from_the_deck_is_the_strain_abaqus_recorded(
        fixture):
    """The one fact in this repository that is not self-referential.

    ``STRAN`` at the integration point is what Abaqus handed the UMAT. This
    reconstructs it from the deck alone -- parse, pick the step, resolve the
    boundary conditions, number the DOFs, build B, evaluate at every
    integration point -- and compares. Nothing in the comparison comes from
    this repository on both sides.
    """
    # Each record is replayed in ITS OWN step. Abaqus numbers increments from
    # 1 again in every step, so a four-step J2 cycle has four increment 1s, and
    # replaying all of them against step 1 compares the reversal and the reload
    # against the first loading -- wrong by 3.0 relative, which is not a
    # tolerance question but a different deformation. The fixture carries the
    # step for exactly this; a fixture frozen before it did says None, and
    # those are the single-step ones, where 1 is right.
    replays: dict = {}

    def _replay(step: int):
        if step not in replays:
            replays[step] = replay_deck(fixture.deck, step=step)
        return replays[step]

    first = _replay(1)
    assert first.fully_prescribed, (
        f"{fixture.path.name}: the verification deck drives every degree of "
        f"freedom, so the strain follows from the deck without solving; this "
        f"one prescribes {len(first.prescribed)} of {first.dof_manager.ndof}")
    assert first.increments, "the step must fix its own increment count"

    worst = 0.0
    for record in fixture.original:
        replay = _replay(int(getattr(record, "step", None) or 1))
        strains = replay.strain_at_points(1, record.increment)
        reference = np.asarray(record.strain, dtype=float)
        scale = max(float(np.max(np.abs(reference))), 1e-30)
        for ip, strain in enumerate(strains):
            worst = max(worst, float(np.max(np.abs(strain - reference))) / scale)
    assert worst < 1e-9, (
        f"{fixture.path.name}: the strain this repository reconstructs from "
        f"the deck differs from the strain Abaqus recorded by {worst:.3e} "
        f"relative. Between the two sit the parser, the step selection, the "
        f"boundary conditions, the DOF map, the connectivity, the "
        f"integration points and the Voigt convention.")


@pytest.mark.integration
def test_the_kinematic_reading_is_read_off_the_step_and_not_assumed():
    """Three readings of the same displacement field, and the deck decides.

    The wrong readings are not nearly right. Measured here, not asserted from
    a book: on a small-strain deck the linear reading is exact and the
    logarithmic one is out by 3.0e-03; on a finite-strain deck the mid-point
    accumulation is exact, ``ln V`` is out by 8.3e-08 and the linear reading
    by 3.0e-03. The 8.3e-08 is the one worth having a test for -- it is small
    enough to be called round-off by anything that does not look.
    """
    measured = {}
    for fixture in hex_fixtures():
        record = fixture.original[-1]
        # The record's OWN step. The last record of a multi-step fixture is
        # not in step 1, and replaying it there compares a reload against the
        # first loading -- which made all three readings equally wrong at 0.5
        # and the test read that as the deck choosing none of them.
        replay = replay_deck(fixture.deck,
                             step=int(getattr(record, "step", None) or 1))
        reference = np.asarray(record.strain, dtype=float)
        scale = max(float(np.max(np.abs(reference))), 1e-30)
        errors = {}
        for reading in ("linear", "logarithmic", "hughes_winget"):
            strains = replay.strain_at_points(1, record.increment,
                                              kinematics=reading)
            errors[reading] = float(np.max(np.abs(strains[0] - reference))) / scale
        measured[fixture.path.name] = (replay.nlgeom, errors)
        best = min(errors, key=errors.get)
        expected = "hughes_winget" if replay.nlgeom else "linear"
        assert best == expected, (
            f"{fixture.path.name}: the deck says NLGEOM={replay.nlgeom} so "
            f"the strain should be the {expected} reading, and the numbers "
            f"say {best}: {errors}")
        assert errors[best] < 1e-9, f"{fixture.path.name}: {errors}"
        others = [e for name, e in errors.items() if name != best]
        assert min(others) > 100.0 * errors[best], (
            f"{fixture.path.name}: the other readings are not distinguishable "
            f"from the right one, so this test is not measuring anything: "
            f"{errors}")

    finite = [e for nlgeom, e in measured.values() if nlgeom]
    assert finite, "no finite-strain fixture is committed"
    assert all(1e-9 < e["logarithmic"] < 1e-6 for e in finite), (
        f"ln V should miss Abaqus's accumulated STRAN by a small but "
        f"resolvable amount; measured {[e['logarithmic'] for e in finite]}")


# --------------------------------------------------------------------------- #
# a mesh, so the scatter has something to get wrong
# --------------------------------------------------------------------------- #
def brick_mesh(nx=2, ny=2, nz=2, size=1.0):
    """An ``nx x ny x nz`` grid of C3D8 elements, Abaqus node ordering.

    Deliberately not a cube of unit elements at the origin: the spacing is
    uneven in every direction so that detJ differs from element to element and
    a scatter that lands in the wrong rows changes the answer.
    """
    xs = np.cumsum([0.0] + [size * (1.0 + 0.3 * i) for i in range(nx)])
    ys = np.cumsum([0.0] + [size * (1.0 + 0.2 * i) for i in range(ny)])
    zs = np.cumsum([0.0] + [size * (1.0 + 0.1 * i) for i in range(nz)])
    nodes, ids_of = {}, {}
    nid = 1
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                nodes[nid] = (float(x), float(y), float(z))
                ids_of[(i, j, k)] = nid
                nid += 1
    model = Model()
    model.nodes = nodes
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                connectivity = [ids_of[(i, j, k)], ids_of[(i + 1, j, k)],
                                ids_of[(i + 1, j + 1, k)], ids_of[(i, j + 1, k)],
                                ids_of[(i, j, k + 1)], ids_of[(i + 1, j, k + 1)],
                                ids_of[(i + 1, j + 1, k + 1)],
                                ids_of[(i, j + 1, k + 1)]]
                model.elements[eid] = Element(eid, "C3D8", connectivity)
                model.element_formulation[eid] = "solid_c3d8_small_strain"
                model.element_material[eid] = "VERIFIED"
                eid += 1
    return model


def assembler_for(fixture, nx=2, ny=2, nz=2):
    model = brick_mesh(nx, ny, nz)
    binding = binding_for(fixture)
    model.materials["VERIFIED"] = binding
    formulations = {"solid_c3d8_small_strain": SolidC3D8SmallStrain()}
    dm = DofManager(model.nodes.keys())
    state = StateManager(n_ip=8, n_state=binding.n_state_vars or 0)
    for eid in model.elements:
        state.init_element(eid, binding.material.init_state(binding, n_ip=8))
    return Assembler(model, dm, formulations, state), model, dm


def a_displacement(ndof: int, seed: int = 3) -> np.ndarray:
    """A displacement that is not a special case: no symmetry, no zeros."""
    rng = np.random.default_rng(seed)
    return 1e-3 * rng.standard_normal(ndof)


@pytest.mark.integration
@pytest.mark.parametrize("fixture", hex_fixtures(), ids=ids(hex_fixtures()))
def test_the_global_stiffness_is_the_derivative_of_the_global_residual(fixture):
    """dR/du, assembled across eight elements, against a difference of R.

    Over a plateau of step sizes, with the componentwise error recorded, and
    with the flatness of the error checked: a disagreement that does not move
    when the step moves is a wrong formula rather than a step that was too
    large, and it is the failure mode that reads exactly like success to a
    single-step check.
    """
    assembler, model, dm = assembler_for(fixture)
    U = a_displacement(dm.ndof)
    R, K, diag = assembler.assemble(U, compute_tangent=True)
    assert diag["elements"] == len(model.elements) == 8
    assert K is not None and K.shape == (dm.ndof, dm.ndof)

    def residual(displacement):
        out, _K, _d = assembler.assemble(displacement)
        return out

    columns = list(range(0, dm.ndof, 7))          # a spread of DOFs, not all
    analytic = K[:, columns]
    sweep = sweep_steps(
        analytic,
        lambda h: centred_jacobian(residual, U, h, columns=columns),
        steps=DEFAULT_STEPS, tolerance=AGREES, scale=1.0,
        notes={"fixture": fixture.path.name, "columns": len(columns)})
    assert not sweep.flat, (
        f"{fixture.path.name}: {sweep.verdict()}\n"
        f"worst components (row, column-of-subset): {sweep.worst_components()}")
    assert sweep.converged, (
        f"{fixture.path.name}: the assembled global stiffness is not the "
        f"derivative of the assembled global residual -- {sweep.verdict()}\n"
        f"worst components: {sweep.worst_components()}")


@pytest.mark.integration
@pytest.mark.parametrize("fixture", hex_fixtures(), ids=ids(hex_fixtures()))
def test_the_shear_step_of_the_deck_is_read_as_engineering_shear(fixture):
    """The factor of two, caught where the frozen window cannot catch it.

    Every committed fixture's window sits inside its deck's UNIAXIAL step, so
    its numbers are silent about engineering versus tensorial shear: the two
    differ by a factor of two in components that are zero there, and the
    mapping layer reports that as not established rather than as agreement.

    The deck's SECOND step is not silent. It prescribes ``ux = c Y`` and
    ``uy = c X`` on the corners of the unit cell -- an affine motion whose
    strain is known in closed form without Abaqus and without this
    repository: the tensor component is ``c`` and the ENGINEERING shear the
    UMAT interface is defined in is ``2c``. Reconstructing ``c`` here instead
    of ``2c`` would halve every shear stress the assembler ever integrates.
    """
    replay = replay_deck(fixture.deck, step=2)
    assert replay.fully_prescribed
    strains = replay.strain_at_points(1, replay.increments,
                                      kinematics="linear")
    prescribed = max(abs(v) for v in replay.prescribed.values())
    assert prescribed > 0.0
    for ip, strain in enumerate(strains):
        assert abs(strain[3] - 2.0 * prescribed) < 1e-12 * prescribed, (
            f"{fixture.path.name} point {ip}: the deck imposes a tensor shear "
            f"of {prescribed:g} and the interface is defined in engineering "
            f"shear, so STRAN[3] must be {2 * prescribed:g}; this reads "
            f"{strain[3]:g}")
        assert abs(strain[0]) < 1e-12 * prescribed, "a pure shear has no extension"


@pytest.mark.integration
def test_a_uniform_stress_over_the_mesh_integrates_to_its_surface_tractions():
    """The divergence theorem across eight elements, not one.

    For a constant stress field the internal force at every INTERIOR node
    cancels exactly -- the two elements either side of a shared face pull the
    same way and opposite -- and what is left on the boundary is the
    consistent nodal force of the tractions on the block's outer surface. The
    reference is built from face quadrature alone, so it shares neither the
    volume weights nor the Jacobians nor the scatter with the thing it is
    compared against: a quadrature weight that is 2% wrong, or a detJ taken in
    the wrong configuration, moves one side and not the other.
    """
    from residual_core.formulations import c3d8_kernel as kernel

    model = brick_mesh()
    dm = DofManager(model.nodes.keys())
    uniform = np.array([120.0, -45.0, 33.0, 17.0, -8.0, 5.0])

    assembled = np.zeros(dm.ndof)
    for eid, element in model.elements.items():
        Xe = model.coords_of(element.connectivity)
        force = kernel.element_internal_force_small_strain(
            Xe, np.tile(uniform, (8, 1)))
        rows = dm.element_dofs(element.connectivity, ("UX", "UY", "UZ"))
        np.add.at(assembled, np.asarray(rows, dtype=int), force)

    faces = np.zeros(dm.ndof)
    for eid, element in model.elements.items():
        Xe = model.coords_of(element.connectivity)
        rows = dm.element_dofs(element.connectivity, ("UX", "UY", "UZ"))
        np.add.at(faces, np.asarray(rows, dtype=int),
                  kernel.surface_traction_nodal_forces(Xe, uniform))

    scale = float(np.max(np.abs(faces)))
    assert scale > 0.0
    assert np.max(np.abs(assembled - faces)) / scale < 1e-12, (
        f"the assembled internal force of a constant stress differs from the "
        f"surface tractions it implies by "
        f"{np.max(np.abs(assembled - faces)) / scale:.3e} of their size. The "
        f"divergence theorem is not an approximation.")

    interior = [nid for nid, xyz in model.nodes.items()
                if all(0.0 < xyz[i] < max(c[i] for c in model.nodes.values())
                       for i in range(3))]
    assert interior, "a mesh with no interior node cannot test cancellation"
    for nid in interior:
        rows = dm.node_dofs(nid)
        assert np.max(np.abs(assembled[rows])) / scale < 1e-12, (
            f"node {nid} is interior to the mesh and carries a net force of "
            f"{assembled[rows]}")


@pytest.mark.integration
def test_a_stress_that_varies_in_space_pairs_with_the_point_it_came_from():
    """The check a uniform stress is provably blind to.

    When every integration point carries the same stress, the pairing of a
    stress with the point it was computed at cannot be wrong -- so the patch
    test above passes whatever order the loop uses. With a stress that varies
    linearly in space the pairing matters, and integration by parts gives an
    independent target built from the surface term and the constant
    divergence, sharing none of the volume pairing it is compared against.

    This settles the pairing INSIDE this repository. It does not settle the
    match to the index Abaqus exports its integration points under: that
    needs a real ODB of a job with a spatially varying field, which is filed
    as abaqus_queue/requests/A3_integration_point_ordering.json.
    """
    from residual_core.formulations import c3d8_kernel as kernel

    model = brick_mesh()
    dm = DofManager(model.nodes.keys())

    def field(x):
        return np.array([100.0 + 7.0 * x[0] - 3.0 * x[1], -40.0 + 2.0 * x[2],
                         11.0 * x[0], 5.0 * x[1], -6.0 * x[2],
                         2.0 * x[0] + x[1]])

    assembled, target = np.zeros(dm.ndof), np.zeros(dm.ndof)
    for eid, element in model.elements.items():
        Xe = model.coords_of(element.connectivity)
        at_points = np.array([field(kernel.shape_functions(p) @ Xe)
                              for p in kernel.ABAQUS_C3D8_GAUSS.points])
        rows = np.asarray(dm.element_dofs(element.connectivity,
                                          ("UX", "UY", "UZ")), dtype=int)
        np.add.at(assembled, rows,
                  kernel.element_internal_force_small_strain(Xe, at_points))
        np.add.at(target, rows, kernel.linear_stress_target(Xe, field))

    scale = float(np.max(np.abs(target)))
    error = float(np.max(np.abs(assembled - target))) / scale
    assert error < 1e-10, (
        f"a stress varying linearly in space integrates to {error:.3e} away "
        f"from an independently built target: the quadrature weights, the "
        f"point positions, or the pairing of a stress with the point it was "
        f"computed at.")


@pytest.mark.integration
def test_the_scatter_puts_each_element_where_the_map_says():
    """Assembled from the shared mesh against a scatter written out by hand.

    Both use the same element residuals, so what is being compared is only
    where they were put. An element loop that lands in the wrong global rows
    on a mesh whose elements share nodes changes the sum; on a single element
    it cannot.
    """
    fixture = hex_fixtures()[0]
    assembler, model, dm = assembler_for(fixture)
    U = a_displacement(dm.ndof)
    R, _K, diag = assembler.assemble(U, collect_elements=True)

    by_hand = np.zeros(dm.ndof)
    for eid, (edofs, r_e) in diag["element_residuals"].items():
        expected = [dm.dof_index(nid, component)
                    for nid in model.elements[eid].connectivity
                    for component in (1, 2, 3)]
        assert list(edofs) == expected, f"element {eid} gathered the wrong rows"
        for local, row in enumerate(expected):
            by_hand[row] += r_e[local]
    assert np.allclose(R, by_hand, rtol=0, atol=1e-12 * max(
        float(np.max(np.abs(R))), 1.0))

    shared = [nid for nid in model.nodes
              if sum(nid in e.connectivity for e in model.elements.values()) > 1]
    assert len(shared) >= 7, (
        "a mesh whose elements share no node cannot test a scatter; this one "
        f"shares {len(shared)} nodes")


@pytest.mark.integration
def test_a_prescribed_face_takes_the_reaction_and_the_free_rows_go_to_zero():
    """The boundary treatment, on an assembled problem rather than in theory.

    A block held on one face and pulled on the opposite one: the free rows are
    what a solver drives to zero, the held rows carry the reaction, and the
    reaction is the negative of the applied load to the precision of the
    solve. A partition that leaks -- a DOF counted free and prescribed, or a
    prescribed DOF left in the free set -- breaks the second statement while
    leaving the first looking fine.
    """
    fixture = hex_fixtures()[0]
    assembler, model, dm = assembler_for(fixture)

    class Bound:
        def __init__(self, target, dof, value=0.0):
            self.target, self.dof_start, self.dof_end = target, dof, dof
            self.value, self.kind, self.amplitude = value, "value", None
            self.step, self.op = 1, "NEW"

    held = [nid for nid, xyz in model.nodes.items() if xyz[0] == 0.0]
    pulled = [nid for nid, xyz in model.nodes.items()
              if abs(xyz[0] - max(c[0] for c in model.nodes.values())) < 1e-12]
    model.boundaries = [Bound(nid, d) for nid in held for d in (1, 2, 3)]
    model.boundaries += [Bound(nid, 1, 5e-4) for nid in pulled]

    free_mask, prescribed, values = constraints.partition(model, dm, step=1)
    assert prescribed.size == len(held) * 3 + len(pulled)
    assert int(np.count_nonzero(free_mask)) + prescribed.size == dm.ndof
    assert not np.any(free_mask[prescribed]), "a DOF is both free and held"

    U = np.zeros(dm.ndof)
    for dof, value in values.items():
        U[dof] = value
    # Newton from the prescribed configuration: the free rows go to zero.
    for _ in range(3):
        R, K, _diag = assembler.assemble(U, compute_tangent=True)
        free = np.where(free_mask)[0]
        step = np.linalg.solve(K[np.ix_(free, free)], -R[free])
        U[free] += step
        if np.max(np.abs(R[free])) < 1e-8:
            break
    R, _K, _diag = assembler.assemble(U, compute_tangent=True)
    free_residual, prescribed_idx, reaction = assembler.split(R)
    scale = max(float(np.max(np.abs(R))), 1.0)
    assert float(np.max(np.abs(free_residual))) / scale < 1e-12, (
        "the free rows are not in equilibrium after the solve")
    net = float(abs(np.sum(reaction[np.isin(prescribed_idx, [
        dm.dof_index(nid, 1) for nid in held])])
        + np.sum(reaction[np.isin(prescribed_idx, [
            dm.dof_index(nid, 1) for nid in pulled])])))
    assert net / scale < 1e-12, (
        f"the x reactions on the held and pulled faces do not cancel: {net:.3e} "
        f"against forces of {scale:.3e}. Newton's third law is not a tolerance.")


@pytest.mark.integration
def test_the_fixtures_numbers_reach_the_assembly_unaltered():
    """The stress the element integrates is the stress the fixture carries.

    At the fixture's own strain the playback material must return the
    fixture's own stress and tangent, bit for bit -- otherwise every number
    downstream is about a material this repository invented.
    """
    for fixture in hex_fixtures():
        binding = binding_for(fixture)
        material = binding.material
        stress, tangent, state, _info = material.evaluate(
            {"strain": material.strain0, "dstrain": np.zeros(fixture.ntens)},
            None, binding, (0.0, 0.0), 0.0, None, None)
        # Matched on the STRAIN the playback was built about, not on an
        # increment number. Abaqus numbers increments from 1 again in every
        # step, so a four-step fixture has four increment 1s, and taking the
        # first match compared a state from step 1 against one from step 4 --
        # 1868.09 against 1002.97, two real numbers from two real increments
        # of the same run. The strain is what "at the fixture's own strain"
        # means, so it is what identifies the record.
        record = [r for r in fixture.converted
                  if np.array_equal(np.asarray(r.strain, dtype=float),
                                    np.asarray(material.strain0, dtype=float))][0]
        assert np.array_equal(stress, record.stress), fixture.path.name
        assert np.array_equal(tangent, record.tangent), fixture.path.name
        assert list(binding.constants) == list(fixture.props), fixture.path.name
