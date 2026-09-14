"""Drive the assembler from the deck a verified fixture was produced by.

A fixture carries the numbers Abaqus computed AND the deck that produced
them. That pairing is the only thing in this project that can check the
assembler's inputs rather than its arithmetic: the strain a UMAT was handed
at an integration point is a fact recorded by Abaqus, and this repository
claims to be able to reconstruct it from the deck. Between those two numbers
sit the parser, the boundary conditions, the step the boundary block belongs
to, the global degree-of-freedom numbering, the element connectivity, the
integration-point positions, the B matrix and the Voigt convention. If any of
them is wrong the reconstructed strain is wrong, and the difference is
measurable rather than argued about.

What it does NOT check: the material. The stress is the UMAT's answer to the
strain, and whether that answer is right was settled in Abaqus, upstream.

Two kinematic readings, because the decks come in both flavours and reading
one as the other is exactly the mistake worth catching:

* ``NLGEOM=NO``  -- Abaqus hands the UMAT the linear strain ``B u``;
* ``NLGEOM=YES`` -- it hands the UMAT an INCREMENTALLY ACCUMULATED strain,
  not the logarithmic strain of the total motion. Measured against
  ``irfancn__Abaqus-UMAT-elastic`` at 0.9% extension: ``ln V`` reproduces
  Abaqus's ``STRAN`` to 8.3e-08 relative and the mid-point (Hughes-Winget)
  accumulation reproduces it to 2.2e-16, with the ``ln V`` discrepancy
  growing linearly in the increment count -- the signature of a per-increment
  integration error rather than of a different tensor.

All three are computed here from the same displacement field, and which one a
deck means is read off its own ``*STEP`` rather than assumed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from residual_core.core import constraints
from residual_core.core.dof_manager import DofManager
from residual_core.core.model import Model, from_abaqus
from residual_core.formulations import c3d8_kernel as kernel
from residual_core.io.abaqus_inp_parser import parse_inp

__all__ = ["ReplayedDeck", "replay_deck", "logarithmic_strain_at_points",
           "linear_strain_at_points", "hughes_winget_strain_at_points",
           "deformation_gradients_at_points", "DeckReplayError"]


class DeckReplayError(ValueError):
    """A deck that cannot be replayed, said in terms of what is missing."""


@dataclass
class ReplayedDeck:
    """One deck, parsed, with its prescribed motion resolved."""

    model: Model
    dof_manager: DofManager
    step: int
    nlgeom: bool
    increments: Optional[int]
    prescribed: dict = field(default_factory=dict)   # global dof -> value
    #: Where each prescribed degree of freedom STOOD when this step began,
    #: which is what the earlier steps left it at. Empty for step 1.
    starts_from: dict = field(default_factory=dict)   # global dof -> value
    warnings: tuple = ()

    @property
    def fully_prescribed(self) -> bool:
        """Every degree of freedom in the model is driven by the deck.

        The single-element verification decks are: every node of the one
        element has all three displacements prescribed, so the displacement
        field is known without solving anything, and a strain computed from it
        is a statement about the deck rather than about a solver.
        """
        return len(self.prescribed) == self.dof_manager.ndof

    def displacement_at(self, increment: int) -> np.ndarray:
        """The nodal displacement vector at the END of ``increment``.

        A fixed-increment step ramps its prescribed values linearly over
        ``increments`` steps, so increment ``i`` ends at ``i / increments`` of
        the way. A step this repository cannot read an increment count off is
        refused rather than guessed at.

        From WHERE THE PREVIOUS STEP LEFT IT, not from zero. Abaqus ramps a
        prescribed displacement from its current value to the one the new step
        names; ``OP=NEW`` replaces which conditions are in force, not the
        displacement the model already has. Ramping from zero instead is right
        for step 1 and wrong for every step after it, and on the four-step J2
        cycle it was wrong by 3.0 relative -- a different deformation, not a
        tolerance. Reading it correctly reproduces Abaqus exactly: step 3
        increment 6 takes e11 from 0 to -0.005 and gamma12 from 0.010 back to
        0, giving -0.003 and 0.004, which is what the run recorded.
        """
        if not self.increments:
            raise DeckReplayError(
                f"step {self.step} does not fix its own increment count "
                f"(initial increment and period do not divide), so which "
                f"displacement increment {increment} stands at is a solver "
                f"decision this deck does not record")
        if not 1 <= int(increment) <= int(self.increments):
            raise DeckReplayError(
                f"increment {increment} is outside step {self.step}, which "
                f"takes {self.increments} increment(s)")
        fraction = float(increment) / float(self.increments)
        U = np.zeros(self.dof_manager.ndof)
        for dof, value in self.prescribed.items():
            began = float(self.starts_from.get(dof, 0.0))
            U[dof] = began + fraction * (float(value) - began)
        return U

    def element_displacement(self, eid: int, U: np.ndarray) -> np.ndarray:
        element = self.model.elements[eid]
        dofs = self.dof_manager.element_dofs(element.connectivity,
                                             ("UX", "UY", "UZ"))
        return np.asarray(U, dtype=float)[np.asarray(dofs, dtype=int)]

    def element_coordinates(self, eid: int) -> np.ndarray:
        return self.model.coords_of(self.model.elements[eid].connectivity)

    def strain_at_points(self, eid: int, increment: int,
                         kinematics: Optional[str] = None) -> np.ndarray:
        """The strain at every integration point, in the deck's own kinematics.

        ``kinematics`` overrides what the step says, for the one purpose of
        measuring how far a different reading would be: ``"linear"``,
        ``"logarithmic"`` or ``"hughes_winget"``.
        """
        Xe = self.element_coordinates(eid)
        reading = kinematics or ("hughes_winget" if self.nlgeom else "linear")
        if reading == "linear":
            Ue = self.element_displacement(eid, self.displacement_at(increment))
            return linear_strain_at_points(Xe, Ue)
        if reading == "logarithmic":
            Ue = self.element_displacement(eid, self.displacement_at(increment))
            return logarithmic_strain_at_points(Xe, Ue)
        if reading == "hughes_winget":
            history = [self.element_displacement(eid, self.displacement_at(i))
                       for i in range(1, int(increment) + 1)]
            return hughes_winget_strain_at_points(Xe, history)
        raise DeckReplayError(
            f"{reading!r} is not a kinematic reading this knows: linear, "
            f"logarithmic or hughes_winget")


def linear_strain_at_points(Xe, Ue, gauss=kernel.ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """``B0 u`` at each integration point: (n_ip, 6) engineering-shear Voigt."""
    Xe = np.asarray(Xe, dtype=float)
    Ue = np.asarray(Ue, dtype=float).reshape(-1)
    out = []
    for point in gauss.points:
        B, _detJ = kernel.b_matrix_reference(Xe, point)
        out.append(B @ Ue)
    return np.asarray(out, dtype=float)


def logarithmic_strain_at_points(Xe, Ue,
                                 gauss=kernel.ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """``ln V`` at each integration point: (n_ip, 6) engineering-shear Voigt.

    ``V`` is the left stretch of ``F = I + du/dX``; the logarithm is taken on
    ``B = F F^T`` by eigen-decomposition, which is exact for the symmetric
    positive-definite tensor it is applied to, and the off-diagonal Voigt
    entries are doubled because Abaqus's ``STRAN`` is engineering shear.
    """
    Xe = np.asarray(Xe, dtype=float)
    ue = np.asarray(Ue, dtype=float).reshape(-1, 3)
    out = []
    for point in gauss.points:
        _N, _detJ, dNdX = kernel._b_from_coords(Xe, point)
        F = np.eye(3) + ue.T @ dNdX
        left = F @ F.T
        values, vectors = np.linalg.eigh(left)
        if np.any(values <= 0.0):
            raise DeckReplayError(
                "the deformation gradient is not invertible at an integration "
                "point, so no logarithmic strain exists there")
        log_v = vectors @ np.diag(0.5 * np.log(values)) @ vectors.T
        out.append(np.array([log_v[0, 0], log_v[1, 1], log_v[2, 2],
                             2.0 * log_v[0, 1], 2.0 * log_v[0, 2],
                             2.0 * log_v[1, 2]]))
    return np.asarray(out, dtype=float)


def deformation_gradients_at_points(Xe, Ue,
                                    gauss=kernel.ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """``F = I + du/dX`` at each integration point: (n_ip, 3, 3)."""
    Xe = np.asarray(Xe, dtype=float)
    ue = np.asarray(Ue, dtype=float).reshape(-1, 3)
    out = []
    for point in gauss.points:
        _N, _detJ, dNdX = kernel._b_from_coords(Xe, point)
        out.append(np.eye(3) + ue.T @ dNdX)
    return np.asarray(out, dtype=float)


def hughes_winget_strain_at_points(Xe, displacement_history,
                                   gauss=kernel.ABAQUS_C3D8_GAUSS) -> np.ndarray:
    """Abaqus's NLGEOM ``STRAN``: the mid-point increment, accumulated.

    Per increment, with ``F_mid = (F_n + F_(n+1)) / 2``::

        dL  = (F_(n+1) - F_n) F_mid^-1
        de  = sym(dL),  dW = skew(dL)
        dR  = (I - dW/2)^-1 (I + dW/2)
        e_(n+1) = dR e_n dR^T + de

    This is not ``ln V`` and the difference is not rounding: at 0.9%
    extension over nine increments ``ln V`` misses Abaqus's own ``STRAN`` by
    8.3e-08 relative while this reproduces it to 2.2e-16. Reading one as the
    other puts an error into every strain a finite-strain deck supplies, in
    the eighth significant figure, where a single-step check would call it
    round-off.

    ``displacement_history`` is the element displacement at the end of each
    increment from the first, starting from an undeformed state.
    """
    Xe = np.asarray(Xe, dtype=float)
    history = [np.zeros(3 * len(Xe))] + [np.asarray(u, dtype=float).reshape(-1)
                                         for u in displacement_history]
    n_ip = len(gauss.points)
    accumulated = [np.zeros((3, 3)) for _ in range(n_ip)]
    identity = np.eye(3)
    for previous, current in zip(history, history[1:]):
        F_prev = deformation_gradients_at_points(Xe, previous, gauss)
        F_curr = deformation_gradients_at_points(Xe, current, gauss)
        for ip in range(n_ip):
            F_mid = 0.5 * (F_prev[ip] + F_curr[ip])
            dL = (F_curr[ip] - F_prev[ip]) @ np.linalg.inv(F_mid)
            de = 0.5 * (dL + dL.T)
            dW = 0.5 * (dL - dL.T)
            dR = np.linalg.solve(identity - 0.5 * dW, identity + 0.5 * dW)
            accumulated[ip] = dR @ accumulated[ip] @ dR.T + de
    return np.asarray([[e[0, 0], e[1, 1], e[2, 2],
                        2.0 * e[0, 1], 2.0 * e[0, 2], 2.0 * e[1, 2]]
                       for e in accumulated], dtype=float)


def replay_deck(deck: str, *, step: int = 1,
                path: Optional[Path] = None,
                formulation_policy=None) -> ReplayedDeck:
    """Parse a deck's text and resolve the motion its ``step`` prescribes."""
    import tempfile

    if path is None:
        handle = tempfile.NamedTemporaryFile("w", suffix=".inp", delete=False,
                                             encoding="utf-8")
        handle.write(deck)
        handle.close()
        path = Path(handle.name)
        written = True
    else:
        written = False
    try:
        parsed = parse_inp(str(path))
    finally:
        if written:
            try:
                Path(path).unlink()
            except OSError:                                # pragma: no cover
                pass

    if not parsed.nodes or not parsed.elements:
        raise DeckReplayError("the deck carries no mesh to replay")
    policy = formulation_policy or (
        lambda etype: "solid_c3d8_small_strain"
        if etype.upper() == "C3D8" else None)
    model = from_abaqus(parsed, formulation_policy=policy)
    node_ids = sorted(model.nodes)
    dm = DofManager(node_ids)

    steps = list(getattr(parsed, "steps", ()) or ())
    if not steps:
        raise DeckReplayError(
            "the deck has no *STEP, so nothing says which boundary conditions "
            "are in force or how long they take to be applied")
    if not 1 <= step <= len(steps):
        raise DeckReplayError(
            f"the deck has {len(steps)} step(s); step {step} is not one of them")
    chosen = steps[step - 1]
    prescribed = constraints.dirichlet_dofs(model, dm, step=step)
    # Where the earlier steps left each degree of freedom. Resolved by reading
    # them in order and letting the latest one that names a dof win, which is
    # what OP=NEW does: it replaces the set of conditions in force, and the
    # displacement carries on from where it was.
    starts_from: dict = {}
    for earlier in range(1, step):
        starts_from.update(constraints.dirichlet_dofs(model, dm, step=earlier))
    return ReplayedDeck(model=model, dof_manager=dm, step=step,
                        nlgeom=bool(chosen.nlgeom),
                        increments=chosen.fixed_increments,
                        prescribed=dict(prescribed),
                        starts_from=starts_from,
                        warnings=tuple(parsed.warnings))
