"""A verified fixture's own numbers, presented to the assembler as a material.

The assembler's job is to turn a stress and a tangent at integration points
into a residual and a Jacobian. To test THAT, the stress and the tangent have
to come from somewhere, and where they come from decides what the test is
worth. A material written for the purpose tests the assembly against a model
invented to make it pass. This one carries the stress and the tangent a real
UMAT produced in Abaqus, at an increment the pipeline verified, and hands them
to the assembler through the same ``Material`` interface every backend uses --
so what gets exercised is the DOF manager, the element loop, the quadrature,
the scatter and the boundary partition, driven by numbers nobody in this
repository chose.

**What it is.** The verified state linearised about itself::

    sigma(eps) = sigma_k + D_k (eps - eps_k)

with ``sigma_k``, ``D_k`` and ``eps_k`` read off one increment of one fixture.
Its tangent is exactly ``D_k`` -- by construction, which is the point: the
assembled ``K`` must then be ``sum B^T D_k B detJ w`` and the finite
difference of the assembled ``R`` must reproduce it, with no constitutive
approximation left in between to hide behind.

**What it is not.** It is not the UMAT. Away from ``eps_k`` it is a linear
extrapolation and the real routine is not, so a residual assembled from it at
a distant strain is not a claim about the material. That is deliberate: this
object exists to make the ASSEMBLY checkable against a difference, and the
constitutive question -- is ``D_k`` the derivative of ``sigma_k`` -- is
answered upstream in Abaqus and again by ``tangent_convention``.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from residual_core.materials.base import Material, MaterialBinding

__all__ = ["FixturePlayback", "binding_for"]


class FixturePlayback(Material):
    """One verified increment of one fixture, as a small-strain material."""

    name = "verified_fixture_playback"
    stress_measure = "cauchy"
    tangent_measure = "ddsdde"
    kinematic_input = "small_strain"
    constitutive_kind = "stress_strain"
    input_variables = ("strain", "dstrain")
    output_variables = ("stress", "tangent")
    supported_formulations = ("material-replay",)
    limitations = (
        "the verified state linearised about itself: exact at the frozen "
        "increment, a linear extrapolation away from it",
        "small-strain interface only; the fixture's kinematics are recorded "
        "on the object rather than applied",
    )

    def __init__(self, fixture, increment: Optional[int] = None):
        record = self._pick(fixture, increment)
        if record.tangent is None:
            raise ValueError(
                f"{fixture.source_id}: the increment chosen carries no "
                f"tangent, so there is nothing for an element stiffness to be "
                f"built out of")
        self.fixture = fixture
        self.increment = record.increment
        self.strain0 = np.asarray(record.strain, dtype=float).copy()
        self.stress0 = np.asarray(record.stress, dtype=float).copy()
        self.tangent = np.asarray(record.tangent, dtype=float).copy()
        self.n_state_vars = int(record.state.size)
        self.state0 = np.asarray(record.state, dtype=float).copy()
        self.notes = (
            f"{fixture.source_id} increment {record.increment}, transform "
            f"fingerprint {fixture.transform_fingerprint}")

    @staticmethod
    def _pick(fixture, increment):
        records = [r for r in fixture.converted if r.tangent is not None]
        if not records:
            raise ValueError(f"{fixture.source_id} carries no tangent at all")
        if increment is None:
            return records[-1]
        for record in records:
            if record.increment == increment:
                return record
        raise ValueError(
            f"{fixture.source_id} carries no tangent at increment "
            f"{increment}; it has {[r.increment for r in records]}")

    def evaluate(self, kinematics: Dict[str, Any], state_prev, binding,
                 time, dtime: float, fields: Optional[dict],
                 options: Optional[dict]):
        strain = np.asarray(kinematics["strain"], dtype=float)
        if strain.size != self.stress0.size:
            raise ValueError(
                f"this fixture's tensor has {self.stress0.size} components "
                f"and the element handed {strain.size}")
        stress = self.stress0 + self.tangent @ (strain - self.strain0)
        state_new = (np.asarray(state_prev, dtype=float)
                     if state_prev is not None and np.size(state_prev)
                     else self.state0.copy())
        return stress, self.tangent.copy(), state_new, {
            "material": self.name, "increment": self.increment}

    def init_state(self, binding, coords=None, n_ip: int = 8) -> np.ndarray:
        n = binding.n_state_vars or self.n_state_vars
        if not n:
            return np.zeros((n_ip, 0), dtype=float)
        return np.tile(self.state0[:n], (n_ip, 1))


def binding_for(fixture, increment: Optional[int] = None,
                name: str = "") -> MaterialBinding:
    """A ``MaterialBinding`` carrying the fixture's own PROPS and NSTATV."""
    material = FixturePlayback(fixture, increment)
    return MaterialBinding(
        material=material,
        constants=list(fixture.props),
        n_state_vars=int(fixture.nstatv or material.n_state_vars or 0),
        name=name or fixture.source_id)
