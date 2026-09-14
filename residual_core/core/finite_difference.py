"""Differencing an assembled derivative, and reading the answer honestly.

A centred difference has two errors moving in opposite directions: truncation
falling as ``h^2`` and cancellation rising as ``eps/h``. One step size sits
somewhere on that curve and cannot be told apart from the one step where two
different matrices happen to cross. What settles it is a run of step sizes
over which the disagreement stops moving -- a PLATEAU -- and the plateau is
the evidence, not the single best number.

The distinction that earns this module its place is the other one. A
derivative that disagrees by the same amount at every step size is not
suffering truncation error: truncation error shrinks when the step shrinks. An
error FLAT in ``h`` is a wrong formula, and it reads exactly like a
converged plateau to anything that only checks "is the best error small".
This project has already paid for that once -- a finite-strain tangent wrong
by 6.2e-02 flat over six decades, correct at 1.88e-11 with a plateau -- so
flatness is measured here and reported as its own verdict rather than left to
be inferred.

Nothing here knows what is being differenced. It is given a number and a way
to produce a difference at a step size, and it reports what the sweep showed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

__all__ = ["DEFAULT_STEPS", "Sweep", "sweep_steps", "centred_jacobian",
           "centred_column", "relative_error", "componentwise_error"]

#: Relative step sizes spanning six decades. Truncation dominates one end and
#: cancellation the other, so a plateau between them is a plateau rather than
#: the whole range being flat.
DEFAULT_STEPS = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)


def relative_error(a, b) -> float:
    """Max absolute difference over the larger of the two scales."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = max(float(np.max(np.abs(a))) if a.size else 0.0,
                float(np.max(np.abs(b))) if b.size else 0.0)
    if scale == 0.0:
        return 0.0
    return float(np.max(np.abs(a - b)) / scale)


def componentwise_error(analytic, differenced) -> np.ndarray:
    """|analytic - differenced| entry by entry, relative to the whole scale.

    Entrywise relative error divides by entries that are legitimately zero, so
    the scale is the largest entry of the object rather than each entry's own
    size: what this answers is "which components carry the disagreement", and
    a component that is zero in both cannot.
    """
    analytic = np.asarray(analytic, dtype=float)
    differenced = np.asarray(differenced, dtype=float)
    scale = max(float(np.max(np.abs(analytic))) if analytic.size else 0.0,
                float(np.max(np.abs(differenced))) if differenced.size else 0.0,
                1e-300)
    return np.abs(analytic - differenced) / scale


@dataclass
class Sweep:
    """What a sweep of step sizes showed, and what it means."""

    steps: tuple
    errors: tuple
    tolerance: float
    #: Componentwise error at the best step, for naming WHICH entries are out.
    componentwise: Optional[np.ndarray] = None
    notes: dict = field(default_factory=dict)

    @property
    def best(self) -> float:
        return min(self.errors) if self.errors else float("inf")

    @property
    def best_step(self) -> float:
        return self.steps[int(np.argmin(self.errors))] if self.errors else float("nan")

    @property
    def agreeing(self) -> tuple:
        return tuple(s for s, e in zip(self.steps, self.errors)
                     if e <= self.tolerance)

    @property
    def plateau(self) -> int:
        """The longest run of CONSECUTIVE steps that all agree."""
        longest = run = 0
        for error in self.errors:
            run = run + 1 if error <= self.tolerance else 0
            longest = max(longest, run)
        return longest

    @property
    def plateau_span(self) -> tuple:
        """The (largest, smallest) step of the longest agreeing run."""
        best_run, best_end, run = 0, -1, 0
        for index, error in enumerate(self.errors):
            run = run + 1 if error <= self.tolerance else 0
            if run > best_run:
                best_run, best_end = run, index
        if best_run == 0:
            return ()
        start = best_end - best_run + 1
        return (self.steps[start], self.steps[best_end])

    @property
    def spread(self) -> float:
        """How much the error moved across the sweep, as a ratio.

        1.0 means it did not move at all.
        """
        finite = [e for e in self.errors if np.isfinite(e) and e > 0.0]
        if len(finite) < 2:
            return float("inf")
        return max(finite) / min(finite)

    @property
    def flat(self) -> bool:
        """The error is the same at every step size AND it is not small.

        Truncation error shrinks with the step; an error that does not is a
        wrong formula. Below tolerance this says nothing -- a derivative that
        is right at 1e-15 everywhere is flat too, and that is not a fault.
        """
        return self.best > self.tolerance and self.spread < 10.0

    @property
    def converged(self) -> bool:
        """A real plateau: at least three consecutive steps inside tolerance."""
        return self.plateau >= 3

    def verdict(self) -> str:
        if self.converged:
            span = self.plateau_span
            return (f"agrees to {self.best:.3e} at h={self.best_step:g}, over a "
                    f"plateau of {self.plateau} step sizes ({span[0]:g} to "
                    f"{span[1]:g})")
        if self.flat:
            return (f"disagrees by {self.best:.3e} and the disagreement does "
                    f"not move with the step size (it varies by a factor of "
                    f"{self.spread:.2f} over {len(self.steps)} step sizes, "
                    f"{max(self.steps):g} to {min(self.steps):g}). Truncation "
                    f"error shrinks when the step shrinks; this does not, so "
                    f"it is a wrong formula rather than a step that was too "
                    f"large")
        if self.plateau:
            return (f"agrees to {self.best:.3e} at h={self.best_step:g} but "
                    f"over only {self.plateau} consecutive step size(s); one "
                    f"step cannot separate truncation error from cancellation")
        return (f"does not agree at any of {len(self.steps)} step sizes; best "
                f"{self.best:.3e} at h={self.best_step:g}, errors "
                f"{['%.2e' % e for e in self.errors]}")

    def worst_components(self, limit: int = 5) -> list:
        """The entries carrying the disagreement, largest first."""
        if self.componentwise is None:
            return []
        flat_errors = np.asarray(self.componentwise, dtype=float).reshape(-1)
        order = np.argsort(flat_errors)[::-1][:limit]
        shape = np.asarray(self.componentwise).shape
        return [(tuple(int(i) for i in np.unravel_index(int(index), shape)),
                 float(flat_errors[int(index)])) for index in order]

    def as_dict(self) -> dict:
        return {"steps": list(self.steps), "errors": list(self.errors),
                "tolerance": self.tolerance, "best": self.best,
                "best_step": self.best_step, "plateau": self.plateau,
                "plateau_span": list(self.plateau_span), "flat": self.flat,
                "spread": self.spread, "converged": self.converged,
                "worst_components": self.worst_components(),
                **self.notes}


def sweep_steps(analytic, difference_at: Callable[[float], np.ndarray], *,
                steps: Sequence[float] = DEFAULT_STEPS,
                tolerance: float = 1e-7,
                scale: float = 1.0,
                notes: Optional[dict] = None) -> Sweep:
    """Difference at every step size and report what the sweep showed.

    ``difference_at(h)`` returns the differenced object at an ABSOLUTE step
    ``h``; ``scale`` multiplies the relative steps into absolute ones, so the
    caller states the size of the quantity being perturbed once.
    """
    analytic = np.asarray(analytic, dtype=float)
    errors, best_error, best_components = [], float("inf"), None
    for relative in steps:
        differenced = np.asarray(
            difference_at(float(relative) * float(scale)), dtype=float)
        error = relative_error(analytic, differenced)
        errors.append(error)
        if error < best_error:
            best_error = error
            best_components = componentwise_error(analytic, differenced)
    return Sweep(steps=tuple(steps), errors=tuple(errors),
                 tolerance=float(tolerance), componentwise=best_components,
                 notes=dict(notes or {}))


def centred_column(force: Callable[[np.ndarray], np.ndarray], base: np.ndarray,
                   index: int, step: float) -> np.ndarray:
    """One column of d(force)/d(base) by a centred difference."""
    plus = np.array(base, dtype=float)
    minus = np.array(base, dtype=float)
    plus[index] += step
    minus[index] -= step
    return (np.asarray(force(plus), dtype=float)
            - np.asarray(force(minus), dtype=float)) / (2.0 * step)


def centred_jacobian(force: Callable[[np.ndarray], np.ndarray],
                     base: np.ndarray, step: float,
                     columns: Optional[Sequence[int]] = None) -> np.ndarray:
    """d(force)/d(base), column by column, by centred differences."""
    base = np.asarray(base, dtype=float)
    indices = list(range(base.size)) if columns is None else list(columns)
    first = centred_column(force, base, indices[0], step)
    jacobian = np.zeros((first.size, len(indices)))
    jacobian[:, 0] = first
    for position, index in enumerate(indices[1:], start=1):
        jacobian[:, position] = centred_column(force, base, index, step)
    return jacobian
