"""Which layer a bad residual came from, asked in the order the layers run.

"The residual is wrong by 12%" tells nobody which of nine people to ask. A
residual assembled from a converted UMAT passes through nine places where it
can go wrong, and they need nine different fixes:

1.  **the UMAT** -- the author's routine returned something that is not a
    number, or a stress that does not move when the strain does;
2.  **the OTI transformation** -- the converted build computes a different
    stress from the original, so the two are not the same model and no
    derivative taken from one applies to the other;
3.  **the constitutive tangent** -- the tangent the converted build reports is
    not the derivative of the stress it reports;
4.  **the state sensitivity** -- the derivative with respect to the internal
    state does not differentiate the residual that uses that state;
5.  **the tensor mapping** -- stress and tangent are right and the assembler
    reads them in a different convention: a transposed Voigt order, or
    tensorial shear where the UMAT meant engineering;
6.  **the element integration** -- the quadrature, the integration-point
    ordering, or the pairing of a stress with the point it belongs to;
7.  **the residual assembly** -- everything handed in is right and the
    integration of ``B^T sigma`` into an element force is not;
8.  **the boundary conditions** -- the free/prescribed partition, or which
    step's displacements are in force;
9.  **the global degree-of-freedom mapping** -- the element's force scattered
    to the wrong global rows.

Three rules hold this together.

**Evidence, not a guess.** Every verdict below carries the numbers it was
reached on, in ``Finding.measured``. A layer named without evidence is worse
than no attribution: it sends somebody to read code that was never implicated.

**Unknown is a third answer.** A layer this diagnosis was given nothing to
check reports ``not_established`` and says what would establish it. It is not
reported as holding, because "nothing was checked" and "checked and clean"
are different facts and a consumer that cannot tell them apart is treating an
absence as a pass.

**The first failure stops the run.** Everything downstream of a NaN is
meaningless, and reporting all nine would bury the one that matters.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

from residual_core.core.finite_difference import Sweep, relative_error

#: The layers, in the order they run. A failure is attributed to the first one
#: that does not hold, because everything after it is downstream of it.
LAYERS = ("umat", "transformation", "constitutive_derivative",
          "state_sensitivity", "mapping", "element_integration",
          "residual_assembly", "boundary_conditions", "global_dof_mapping")

#: Kept under its old name because the five it used to have are five of these
#: nine and callers index by name.
STAGES = LAYERS

#: Who fixes each layer, in one line. Carried in the failure message because
#: the point of attributing a failure is to route it.
OWNER = {
    "umat": "the author's own routine, as published",
    "transformation": "the OTI source transformation",
    "constitutive_derivative": "the derivative the converted build extracts",
    "state_sensitivity": "the internal-state derivative and how it is carried",
    "mapping": "the convention this assembler reads the UMAT's arrays in",
    "element_integration": "the quadrature and integration-point order here",
    "residual_assembly": "the integration of B^T sigma in this repository",
    "boundary_conditions": "the prescribed/free partition in this repository",
    "global_dof_mapping": "the element-to-global scatter in this repository",
}

HOLDS, FAILS, NOT_ESTABLISHED = "holds", "fails", "not_established"


@dataclass
class Finding:
    """One layer's verdict, and the numbers behind it."""

    stage: str
    status: str = HOLDS
    detail: str = ""
    measured: dict = field(default_factory=dict)
    #: What would settle a layer that could not be checked. Empty otherwise.
    would_establish: str = ""

    def __post_init__(self):
        if isinstance(self.status, bool):        # historical ok=True/False
            self.status = HOLDS if self.status else FAILS

    @property
    def ok(self) -> bool:
        """True only for a layer that was checked and held."""
        return self.status == HOLDS

    @property
    def failed(self) -> bool:
        return self.status == FAILS

    def as_dict(self) -> dict:
        return {"stage": self.stage, "layer": self.stage, "status": self.status,
                "ok": self.ok, "detail": self.detail,
                "measured": dict(self.measured),
                "would_establish": self.would_establish}

    def message(self) -> str:
        if self.status == HOLDS:
            return f"{self.stage}: holds. {self.detail}"
        if self.status == NOT_ESTABLISHED:
            return (f"{self.stage}: NOT ESTABLISHED -- {self.detail}"
                    + (f"\n  what would establish it: {self.would_establish}"
                       if self.would_establish else ""))
        return (f"{self.stage} is where this goes wrong -- {OWNER[self.stage]}."
                f"\n  {self.detail}"
                + ("\n  measured: " + ", ".join(
                    f"{name}={value!r}" for name, value in self.measured.items())
                   if self.measured else ""))


@dataclass
class Diagnosis:
    """Every layer checked, the first that failed, and what went unchecked."""

    findings: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """Every layer that was CHECKED held.

        A diagnosis with unestablished layers is not a failure and it is not a
        clean bill either; :attr:`complete` is the one that says so.
        """
        return not any(f.failed for f in self.findings)

    @property
    def complete(self) -> bool:
        return bool(self.findings) and all(f.ok for f in self.findings)

    @property
    def blame(self) -> Optional[str]:
        for finding in self.findings:
            if finding.failed:
                return finding.stage
        return None

    @property
    def not_established(self) -> list:
        return [f.stage for f in self.findings if f.status == NOT_ESTABLISHED]

    def finding(self, stage: str) -> Optional[Finding]:
        for found in self.findings:
            if found.stage == stage:
                return found
        return None

    def report(self) -> str:
        lines = [finding.message() for finding in self.findings]
        if self.blame:
            return (f"the failure is in {self.blame} -- {OWNER[self.blame]}.\n"
                    + "\n".join(lines))
        head = ("every layer that was checked holds"
                if not self.complete else
                "every layer holds: the UMAT computed numbers, the conversion "
                "agrees with it, the tangent is the derivative of the stress, "
                "the state derivative differentiates the residual, the "
                "conventions match, the quadrature integrates, the assembly "
                "sums, the partition is a partition, and the scatter lands "
                "where the map says")
        if self.not_established:
            head += (f"; NOT established: {', '.join(self.not_established)}")
        return head + ".\n  " + "\n  ".join(lines)

    def as_dict(self) -> dict:
        return {"ok": self.ok, "complete": self.complete, "blame": self.blame,
                "not_established": list(self.not_established),
                "findings": [f.as_dict() for f in self.findings]}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _relative(a, b) -> float:
    return relative_error(a, b)


#: The element geometry the corpus verification decks are generated on.
UNIT_CUBE = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                      [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])


@dataclass
class StateSensitivityProbe:
    """A dR/dq claim and a way to difference it, supplied by the caller.

    ``analytic`` is what the assembler produced; ``difference_at(h)`` returns
    the same object differenced at an absolute step ``h`` in the state.
    ``scale`` turns the relative sweep steps into absolute ones.
    """

    analytic: np.ndarray
    difference_at: Callable[[float], np.ndarray]
    scale: float = 1.0
    label: str = "dR/dq"
    tolerance: float = 1e-6


# --------------------------------------------------------------------------- #
# the layers
# --------------------------------------------------------------------------- #
def _check_umat(fixture) -> Finding:
    stresses = [record.stress for record in fixture.original]
    strains = [record.strain for record in fixture.original]
    non_finite = [record.increment for record in fixture.original
                  if not _finite(record.stress) or not _finite(record.state)]
    if non_finite:
        return Finding(
            "umat", FAILS,
            f"the ORIGINAL routine returned values that are not numbers at "
            f"increment(s) {non_finite}. Nothing downstream of a NaN is a "
            f"measurement.",
            {"increments": non_finite})
    moved = _relative(stresses[0], stresses[-1]) if len(stresses) > 1 else 0.0
    if len(stresses) > 1 and moved == 0.0 and _relative(strains[0], strains[-1]) > 0:
        return Finding(
            "umat", FAILS,
            "the strain changed over the history and the stress did not, so "
            "the routine is not responding to what it is given.",
            {"strain_change": _relative(strains[0], strains[-1])})
    return Finding(
        "umat", HOLDS,
        f"{len(stresses)} increments, all finite, stress responding to strain",
        {"stress_change": moved, "increments": len(stresses)})


def _check_transformation(fixture, tolerance: float) -> Finding:
    worst = max((_relative(a.stress, b.stress)
                 for a, b in zip(fixture.original, fixture.converted)),
                default=0.0)
    worst_state = max((_relative(a.state, b.state)
                       for a, b in zip(fixture.original, fixture.converted)
                       if a.state.size), default=0.0)
    if worst > tolerance or worst_state > tolerance:
        return Finding(
            "transformation", FAILS,
            f"the converted build's stress differs from the original's by "
            f"{worst:.3e} and its state by {worst_state:.3e}. They are not the "
            f"same model, so no derivative from one describes the other.",
            {"stress": worst, "state": worst_state})
    return Finding(
        "transformation", HOLDS,
        f"original and converted agree to {worst:.3e} in stress and "
        f"{worst_state:.3e} in state over {fixture.increments()} increments",
        {"stress": worst, "state": worst_state,
         "gates_not_true": list(getattr(fixture, "gates_not_true", ()))})


def _check_constitutive_derivative(fixture, tolerance: float) -> Finding:
    from residual_core.materials.verified_fixture import tangent_convention

    convention = tangent_convention(fixture, tolerance=tolerance)
    measured = {"material_reading": convention.material_error,
                "jaumann_reading": convention.jaumann_error,
                "reading": convention.reading,
                "separation": convention.separation,
                "increments_used": convention.increments_used,
                "volumetric_strain_increment": convention.volumetric,
                "prediction_overshoot": convention.prediction_overshoot}
    if convention.increments_used == 0:
        return Finding(
            "constitutive_derivative", NOT_ESTABLISHED, convention.why,
            measured,
            would_establish="a fixture window carrying a moving strain and a "
                            "tangent at two consecutive increments")
    if convention.reading is None and convention.volumetric == 0.0:
        return Finding(
            "constitutive_derivative", NOT_ESTABLISHED, convention.why,
            measured,
            would_establish="a window whose strain increment changes the "
                            "volume, so sigma tr(Deps) is not zero")
    if convention.reading is None and convention.prediction_overshoot > 10.0:
        return Finding(
            "constitutive_derivative", NOT_ESTABLISHED, convention.why,
            measured,
            would_establish=("the ELASTIC strain increment this material's "
                             "stress answers -- the total strain minus "
                             "whatever its state removes -- which means "
                             "reading its state layout, not its STRAN"))
    if convention.reading is None:
        worst, at, rows = _derivative_check(fixture)
        direct, shear = rows
        detail = (f"the tangent the converted build reports does not predict "
                  f"the stress increment it reports under either reading: "
                  f"{convention.why} (worst {worst:.3e} at increment {at}). A "
                  f"tangent that does not differentiate its own stress is not "
                  f"a tangent.")
        if shear > tolerance >= direct:
            detail += (f" The disagreement is confined to the SHEAR rows "
                       f"({shear:.3e} against {direct:.3e} on the direct "
                       f"rows), which is the signature of a shear convention "
                       f"rather than of a wrong tangent -- see the mapping "
                       f"layer.")
        measured.update({"worst": worst, "increment": at,
                         "direct_rows": direct, "shear_rows": shear})
        return Finding("constitutive_derivative", FAILS, detail, measured)
    return Finding(
        "constitutive_derivative", HOLDS,
        f"the tangent differentiates its own stress: {convention.why}. The "
        f"assembler must build this material's stiffness in the "
        f"{convention.reading!r} reading", measured)


def _check_state_sensitivity(fixture,
                             probe: Optional[StateSensitivityProbe],
                             steps) -> Finding:
    if probe is None:
        moved = _state_movement(fixture)
        return Finding(
            "state_sensitivity", NOT_ESTABLISHED,
            f"no dR/dq was offered to difference. The fixture's own state "
            f"moves by {moved:.3e} over its window, which says the material "
            f"has a state -- it does not say anything about the derivative "
            f"with respect to it.",
            {"state_movement": moved,
             "state_components": int(fixture.original[-1].state.size)},
            would_establish=("a dR/dq assembled at a frozen state, with a way "
                             "to re-assemble R at a perturbed state, so the "
                             "two can be differenced over a plateau"))
    from residual_core.core.finite_difference import sweep_steps

    sweep = sweep_steps(probe.analytic, probe.difference_at, steps=steps,
                        tolerance=probe.tolerance, scale=probe.scale)
    if not sweep.converged:
        return Finding(
            "state_sensitivity", FAILS,
            f"{probe.label} does not differentiate the residual it belongs "
            f"to: {sweep.verdict()}.",
            {"best": sweep.best, "plateau": sweep.plateau, "flat": sweep.flat,
             "errors": list(sweep.errors),
             "worst_components": sweep.worst_components()})
    return Finding(
        "state_sensitivity", HOLDS,
        f"{probe.label} {sweep.verdict()}",
        {"best": sweep.best, "plateau": sweep.plateau,
         "plateau_span": list(sweep.plateau_span)})


def _check_mapping(fixture) -> Finding:
    from residual_core.materials.verified_fixture import check_conventions

    problems = check_conventions(fixture)
    shear = _shear_exercised(fixture)
    if problems:
        return Finding(
            "mapping", FAILS,
            "the fixture's own numbers contradict the convention it states: "
            + "; ".join(problems), {"shear_exercised": shear})
    detail = (f"{fixture.ntens} components as {fixture.ndi} direct and "
              f"{fixture.nshr} engineering shear, and the tangent is shaped "
              f"and signed accordingly")
    measured = {"shear_exercised": shear, "ntens": fixture.ntens,
                "ndi": fixture.ndi, "nshr": fixture.nshr,
                "exercised_at": SHEAR_IS_EXERCISED_AT}
    if shear <= SHEAR_IS_EXERCISED_AT:
        return Finding(
            "mapping", NOT_ESTABLISHED,
            detail + f". But the carried window applies effectively NO shear "
                     f"strain -- the largest shear component is {shear:.2e} "
                     f"of the largest extension -- so these numbers cannot "
                     f"tell engineering shear from tensorial: the two differ "
                     f"by a factor of two in components that are all zero "
                     f"here.",
            measured,
            would_establish=("a fixture window inside the deck's shear step, "
                             "where a factor of two in the shear convention "
                             "changes the stress"))
    return Finding("mapping", HOLDS,
                   detail + f", over a window whose largest shear strain is "
                            f"{shear:.3e} of its largest extension -- enough "
                            f"that a factor of two would show", measured)


def _check_element_integration(coordinates, gauss) -> Finding:
    """Quadrature, point ordering and the stress-to-point pairing.

    Three measurements, each independent of the material:

    * a CONSTANT stress integrates to the consistent nodal forces of the
      surface tractions it implies (the divergence theorem);
    * a stress VARYING linearly in space integrates to an independent target
      built from the surface term and the constant divergence -- which the
      uniform-stress check is provably blind to, because every point carries
      the same stress there and the pairing cannot be wrong;
    * a rigid translation produces no strain at any point.
    """
    from residual_core.formulations import c3d8_kernel as kernel

    Xe = np.asarray(coordinates, dtype=float)
    volumes = []
    for point in gauss.points:
        try:
            _B, detJ, _dNdx = kernel._b_from_coords(Xe, point)
        except np.linalg.LinAlgError:
            return Finding(
                "element_integration", FAILS,
                "the element map is not invertible at an integration point, "
                "so there is no B matrix there and nothing can be integrated. "
                "The geometry is degenerate, not the quadrature.",
                {"points": int(len(gauss.points))})
        volumes.append(float(detJ))
    if min(volumes) <= 0.0:
        return Finding(
            "element_integration", FAILS,
            f"the Jacobian determinant is {min(volumes):.3e} at an "
            f"integration point. An element that is inverted or flat "
            f"integrates a negative volume, and every force taken out of it "
            f"has the wrong sign or none.",
            {"min_detJ": min(volumes), "max_detJ": max(volumes)})

    uniform = np.array([120.0, -45.0, 33.0, 17.0, -8.0, 5.0])
    mine = kernel.element_internal_force_small_strain(
        Xe, np.tile(uniform, (len(gauss.weights), 1)), gauss=gauss)
    tractions = kernel.surface_traction_nodal_forces(Xe, uniform)
    uniform_error = _relative(mine, tractions)

    def field(x):
        return np.array([100.0 + 7.0 * x[0] - 3.0 * x[1],
                         -40.0 + 2.0 * x[2], 11.0 * x[0],
                         5.0 * x[1], -6.0 * x[2], 2.0 * x[0] + x[1]])

    varying = np.array([field(kernel.shape_functions(p) @ Xe)
                        for p in gauss.points])
    mine_varying = kernel.element_internal_force_small_strain(Xe, varying,
                                                              gauss=gauss)
    target = kernel.linear_stress_target(Xe, field, gauss=gauss)
    varying_error = _relative(mine_varying, target)

    translation = np.tile(np.array([0.3, -0.2, 0.1]), len(Xe))
    strains = [kernel.b_matrix_reference(Xe, p)[0] @ translation
               for p in gauss.points]
    translation_error = float(np.max(np.abs(strains)))

    measured = {"uniform_stress_vs_tractions": uniform_error,
                "linear_stress_vs_independent_target": varying_error,
                "rigid_translation_strain": translation_error,
                "points": int(len(gauss.points))}
    if uniform_error > 1e-10:
        return Finding("element_integration", FAILS,
                       f"a constant stress does not integrate to the surface "
                       f"tractions it implies: {uniform_error:.3e} relative. "
                       f"The divergence theorem is not an approximation.",
                       measured)
    if varying_error > 1e-10:
        return Finding("element_integration", FAILS,
                       f"a stress varying linearly in space integrates to "
                       f"{varying_error:.3e} away from an independently built "
                       f"target. A uniform stress cannot see this: it is the "
                       f"quadrature weights, the point positions or the "
                       f"pairing of a stress with the point it was computed "
                       f"at.", measured)
    if translation_error > 1e-12:
        return Finding("element_integration", FAILS,
                       f"a rigid translation produces a strain of "
                       f"{translation_error:.3e} at an integration point.",
                       measured)
    return Finding(
        "element_integration", HOLDS,
        f"a constant stress integrates to its surface tractions "
        f"({uniform_error:.3e}), a linearly varying one to an independent "
        f"target ({varying_error:.3e}) -- which is what tests the "
        f"stress-to-point pairing -- and a rigid translation strains nothing "
        f"({translation_error:.3e})", measured)


def _check_residual_assembly(fixture, assemble, reference, gauss,
                             tolerance: float) -> Finding:
    if assemble is None or reference is None:
        return Finding(
            "residual_assembly", NOT_ESTABLISHED,
            "no assembly and no independent integration were offered, so "
            "nothing here compares one against the other",
            would_establish=("an ``assemble(stress_at_ip)`` and a second, "
                             "independently written integration of the same "
                             "quantity"))
    stress_at_ip = np.tile(fixture.original[-1].stress,
                           (len(gauss.weights), 1))
    mine = np.asarray(assemble(stress_at_ip), dtype=float)
    theirs = np.asarray(reference(stress_at_ip), dtype=float)
    difference = _relative(mine, theirs)
    if difference > tolerance:
        return Finding(
            "residual_assembly", FAILS,
            f"the assembled internal force differs from the independent "
            f"reference by {difference:.3e}. Everything handed to it agreed, "
            f"so this is the integration.", {"relative": difference})
    equilibrium = float(np.max(np.abs(
        mine.reshape(-1, 3).sum(axis=0)))) / max(float(np.max(np.abs(mine))), 1.0)
    return Finding(
        "residual_assembly", HOLDS,
        f"B^T sigma integrates to the independent reference within "
        f"{difference:.3e}, and the nodal forces of a uniform stress cancel "
        f"to {equilibrium:.3e} of their own size",
        {"relative": difference, "net_force_relative": equilibrium})


def _check_boundary_conditions(replay) -> Finding:
    if replay is None:
        return Finding(
            "boundary_conditions", NOT_ESTABLISHED,
            "no model with boundary conditions was offered, so nothing here "
            "says which degrees of freedom are prescribed",
            would_establish=("a parsed deck, so the free/prescribed partition "
                             "can be compared against what the deck's step "
                             "says"))
    from residual_core.core import constraints

    dm = replay.dof_manager
    free_mask, prescribed, values = constraints.partition(
        replay.model, dm, step=replay.step)
    n_free = int(np.count_nonzero(free_mask))
    overlap = bool(np.any(free_mask[prescribed])) if prescribed.size else False
    covers = (n_free + prescribed.size) == dm.ndof
    flat = constraints.dirichlet_dofs(replay.model, dm, step=None)
    stepwise = constraints.dirichlet_dofs(replay.model, dm, step=replay.step)
    disagreeing = sum(1 for dof, value in stepwise.items()
                      if abs(flat.get(dof, 0.0) - value) > 1e-15)
    measured = {"ndof": dm.ndof, "prescribed": int(prescribed.size),
                "free": n_free, "step": replay.step,
                "dofs_a_flat_read_would_get_wrong": disagreeing}
    if overlap or not covers:
        return Finding(
            "boundary_conditions", FAILS,
            f"the free and prescribed sets are not a partition of the "
            f"{dm.ndof} degrees of freedom: {n_free} free and "
            f"{prescribed.size} prescribed"
            + (", and they overlap" if overlap else ""), measured)
    if not replay.fully_prescribed and prescribed.size == 0:
        return Finding(
            "boundary_conditions", FAILS,
            "the deck's step prescribes no degree of freedom at all, so the "
            "residual it drives has a null space", measured)
    return Finding(
        "boundary_conditions", HOLDS,
        f"step {replay.step} prescribes {prescribed.size} of {dm.ndof} "
        f"degrees of freedom and leaves {n_free} free, a partition with no "
        f"overlap; reading the deck's steps as one flat list would get "
        f"{disagreeing} of them wrong", measured)


def _check_global_dof_mapping(replay) -> Finding:
    if replay is None:
        return Finding(
            "global_dof_mapping", NOT_ESTABLISHED,
            "no model was offered, so nothing here says where an element's "
            "force lands in the global vector",
            would_establish=("a parsed deck, so the element-to-global map can "
                             "be checked for being a dense bijection and the "
                             "scatter for following it"))
    dm = replay.dof_manager
    seen = {}
    for nid in dm.node_ids:
        for local in range(1, dm.n_node_dofs(nid) + 1):
            index = dm.dof_index(nid, local)
            if index in seen:
                return Finding(
                    "global_dof_mapping", FAILS,
                    f"node {nid} degree of freedom {local} and node "
                    f"{seen[index][0]} degree of freedom {seen[index][1]} "
                    f"both map to global row {index}; two different "
                    f"unknowns cannot share a row",
                    {"collision_row": index})
            seen[index] = (nid, local)
    if sorted(seen) != list(range(dm.ndof)):
        missing = sorted(set(range(dm.ndof)) - set(seen))
        return Finding(
            "global_dof_mapping", FAILS,
            f"the map is not dense: {len(missing)} global row(s) belong to no "
            f"node, the first being {missing[:3]}", {"missing": missing[:10]})

    # The scatter follows the map rather than the order the nodes happen to be
    # written in: gather an element's DOFs, and check each one is the row its
    # own (node, component) resolves to.
    worst = 0
    for eid, element in replay.model.elements.items():
        dofs = dm.element_dofs(element.connectivity, ("UX", "UY", "UZ"))
        expected = [dm.dof_index(nid, component)
                    for nid in element.connectivity
                    for component in (1, 2, 3)]
        if list(dofs) != expected:
            return Finding(
                "global_dof_mapping", FAILS,
                f"element {eid} gathers global rows {list(dofs)[:6]}... and "
                f"its connectivity resolves to {expected[:6]}...; the scatter "
                f"is not following the map", {"element": eid})
        worst = max(worst, max(dofs))
    return Finding(
        "global_dof_mapping", HOLDS,
        f"{dm.ndof} global rows, one per (node, component), dense and without "
        f"collision, and every element gathers the rows its own connectivity "
        f"resolves to (highest row used: {worst})",
        {"ndof": dm.ndof, "elements": len(replay.model.elements)})


# --------------------------------------------------------------------------- #
# the run
# --------------------------------------------------------------------------- #
def diagnose(fixture, *, assemble: Optional[Callable] = None,
             reference: Optional[Callable] = None,
             state_probe: Optional[StateSensitivityProbe] = None,
             replay=None, use_deck: bool = True,
             coordinates=None, gauss=None,
             steps: Sequence[float] = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8),
             primal_tolerance: float = 1e-9,
             derivative_tolerance: float = 1e-5,
             assembly_tolerance: float = 1e-10) -> Diagnosis:
    """Run every layer's check against one verified fixture.

    ``assemble(stress_at_ip)`` returns this repository's element internal
    force and ``reference(stress_at_ip)`` the same quantity computed some
    other way, so the assembly layer is checked against something rather than
    against itself. ``state_probe`` supplies a dR/dq and a way to difference
    it. ``replay`` is a parsed deck for the boundary-condition and DOF-map
    layers; with ``use_deck`` the fixture's OWN deck is parsed for that, which
    is the deck the numbers were produced under.

    Layers with nothing to check report ``not_established``. They are not
    reported as holding.
    """
    from residual_core.formulations.c3d8_kernel import ABAQUS_C3D8_GAUSS

    gauss = gauss or ABAQUS_C3D8_GAUSS
    findings: list = []

    findings.append(_check_umat(fixture))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_transformation(fixture, primal_tolerance))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_constitutive_derivative(fixture,
                                                   derivative_tolerance))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_state_sensitivity(fixture, state_probe, steps))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_mapping(fixture))
    if findings[-1].failed:
        return Diagnosis(findings)

    if replay is None and use_deck:
        replay = _replay_of(fixture)
    if coordinates is None:
        coordinates = _coordinates_of(replay, fixture)

    findings.append(_check_element_integration(coordinates, gauss))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_residual_assembly(fixture, assemble, reference,
                                             gauss, assembly_tolerance))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_boundary_conditions(replay))
    if findings[-1].failed:
        return Diagnosis(findings)

    findings.append(_check_global_dof_mapping(replay))
    return Diagnosis(findings)


def _replay_of(fixture):
    """The fixture's own deck, parsed, or None with the reason swallowed.

    A deck that cannot be replayed leaves the last two layers unestablished,
    which is reported; it is not a reason to fail a diagnosis of the first
    seven.
    """
    from residual_core.core.deck_replay import DeckReplayError, replay_deck

    if not getattr(fixture, "deck", ""):
        return None
    try:
        return replay_deck(fixture.deck, step=1)
    except (DeckReplayError, ValueError, KeyError):
        return None


def _coordinates_of(replay, fixture):
    if replay is not None and replay.model.elements:
        eid = sorted(replay.model.elements)[0]
        coords = replay.element_coordinates(eid)
        if coords.shape == (8, 3):
            return coords
    return UNIT_CUBE


def _state_movement(fixture) -> float:
    states = [record.state for record in fixture.original if record.state.size]
    if len(states) < 2:
        return 0.0
    return _relative(states[0], states[-1])


#: How large a shear strain has to be, as a fraction of the window's largest
#: direct strain, before a factor of two in the shear convention would change
#: a number in the fixture. Abaqus writes -0.0 and values around 1e-20 into
#: components a deck never drives; doubling 1e-20 is still 1e-20, and reading
#: that as "the shear convention was exercised" would report a check that
#: cannot fail as one that passed.
SHEAR_IS_EXERCISED_AT = 1e-6


def _shear_exercised(fixture) -> float:
    """The largest shear strain the window applies, relative to its extension.

    Zero -- or float noise, which is the same thing here -- means the window
    says nothing about the shear convention: a factor of two on components
    that are all zero changes no number in it.
    """
    if fixture.ndi >= fixture.ntens:
        return 0.0
    shear, direct = 0.0, 0.0
    for record in fixture.original:
        if record.strain.size >= fixture.ntens:
            shear = max(shear, float(np.max(np.abs(
                record.strain[fixture.ndi:fixture.ntens]))))
            direct = max(direct, float(np.max(np.abs(
                record.strain[:fixture.ndi]))))
    if direct <= 0.0:
        return shear
    return shear / direct


def _derivative_check(fixture) -> tuple:
    """Does ``D dstrain`` predict the stress increment the same build reported?

    A weaker question than the one the UMAT pipeline asks -- it perturbs the
    routine and this only reads what it wrote -- and a different one: it asks
    whether the tangent and the stress in THIS fixture belong together. That
    is the question the assembler needs answered, because it uses both.

    Returns the worst relative disagreement, the increment it was at, and the
    split of that disagreement between the direct rows and the shear rows --
    because a disagreement confined to the shear rows is a convention and not
    a wrong tangent.
    """
    worst, at, rows = None, None, (0.0, 0.0)
    for previous, current in zip(fixture.converted, fixture.converted[1:]):
        if current.tangent is None or current.dstrain.size != fixture.ntens:
            continue
        change = np.asarray(current.stress) - np.asarray(previous.stress)
        scale = float(np.max(np.abs(change)))
        if scale <= 0.0:
            continue
        predicted = current.tangent @ np.asarray(current.dstrain, dtype=float)
        relative = _relative(predicted, change)
        if worst is None or relative > worst:
            worst, at = relative, current.increment
            residue = np.abs(predicted - change) / scale
            ndi = fixture.ndi or fixture.ntens
            rows = (float(np.max(residue[:ndi])) if ndi else 0.0,
                    float(np.max(residue[ndi:])) if residue.size > ndi else 0.0)
    return worst, at, rows
