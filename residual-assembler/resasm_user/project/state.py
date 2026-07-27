"""Pipeline state machine for a sensitivity project (M5-A).

A project advances through an ordered set of states. Each action requires a
minimum state to run and, on success, advances the project to the state it
produces. A failed action leaves the project at its last successful state, so
the user can fix the cause and resume from exactly where it stopped.
"""

from __future__ import annotations

from enum import IntEnum


class State(IntEnum):
    NEW = 0
    INSPECTED = 1
    CONFIGURED = 2
    PREPARED = 3
    COMPILED = 4
    SOLVED = 5
    EXTRACTED = 6
    SENSITIVITIES_COMPLETE = 7
    VALIDATED = 8
    REPORTED = 9

    @classmethod
    def parse(cls, name: str) -> "State":
        try:
            return cls[str(name).upper()]
        except KeyError:
            raise ValueError("unknown pipeline state %r" % (name,))

    def __str__(self) -> str:              # human-readable name, not "State.NEW"
        return self.name
