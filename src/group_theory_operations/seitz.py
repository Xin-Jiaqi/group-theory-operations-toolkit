"""Seitz operations (R | t) for crystallographic space groups.

A Seitz operation acts on fractional coordinates as ``x' = R x + t`` with
``R`` an integer rotation matrix and ``t`` a fractional translation taken
modulo the lattice (mod 1 per axis).  The package stores space-group
generators in this form; this module provides the multiplication, inverse,
and closure machinery used to verify that the stored generators really
generate the declared operation count.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np

RotationMatrix = np.ndarray
Translation = np.ndarray

_TRANSLATION_TOLERANCE = 1.0e-9


@dataclass(frozen=True, slots=True)
class SeitzOp:
    """A space-group operation in fractional coordinates.

    ``rotation`` is an exact integer 3x3 matrix; ``translation`` is a
    float triplet normalized into [0, 1) modulo the lattice.
    """

    rotation: np.ndarray
    translation: np.ndarray

    def __post_init__(self) -> None:
        rotation = np.asarray(self.rotation)
        translation = np.asarray(self.translation, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError("rotation must be a 3x3 matrix")
        if translation.shape != (3,):
            raise ValueError("translation must be a length-3 vector")
        if not np.all(np.equal(np.rint(rotation), rotation)):
            raise ValueError("rotation must contain exact integers")
        rotation = np.asarray(rotation, dtype=np.int64)
        determinant = int(round(float(np.linalg.det(rotation))))
        if determinant not in {-1, 1}:
            raise ValueError("rotation must be unimodular (determinant +1 or -1)")
        normalized = translation % 1.0
        normalized[np.abs(normalized) < _TRANSLATION_TOLERANCE] = 0.0
        normalized[np.abs(normalized - 1.0) < _TRANSLATION_TOLERANCE] = 0.0
        object.__setattr__(self, "rotation", rotation)
        object.__setattr__(self, "translation", normalized)

    def apply(self, position: np.ndarray) -> np.ndarray:
        """Return the image of a fractional-coordinate position."""
        return (self.rotation @ position + self.translation) % 1.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "rotation": [list(map(int, row)) for row in self.rotation],
            "translation": [float(value) for value in self.translation],
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> SeitzOp:
        return cls(
            rotation=np.asarray(value["rotation"]),
            translation=np.asarray(value["translation"], dtype=np.float64),
        )

    @classmethod
    def identity(cls) -> SeitzOp:
        return cls(rotation=np.eye(3, dtype=np.int64), translation=np.zeros(3))

    def is_identity(self) -> bool:
        return (
            bool(np.array_equal(self.rotation, np.eye(3, dtype=np.int64)))
            and float(np.max(np.abs(self.translation))) < _TRANSLATION_TOLERANCE
        )


def multiply(left: SeitzOp, right: SeitzOp) -> SeitzOp:
    """Return the product ``left * right`` acting as left-after-right."""
    rotation = left.rotation @ right.rotation
    translation = left.rotation @ right.translation + left.translation
    return SeitzOp(rotation=rotation, translation=translation)


def inverse(operation: SeitzOp) -> SeitzOp:
    # Fractional-coordinate rotations preserve the lattice metric, not the
    # Euclidean coordinate metric.  In a non-orthogonal basis (for example a
    # hexagonal conventional cell), R.T is therefore not generally R^-1.
    rotation = np.rint(np.linalg.inv(operation.rotation)).astype(np.int64)
    if not np.array_equal(operation.rotation @ rotation, np.eye(3, dtype=np.int64)):
        raise ValueError("rotation has no exact integer inverse")
    translation = -rotation @ operation.translation
    return SeitzOp(rotation=rotation, translation=translation)


def equivalent(left: SeitzOp, right: SeitzOp) -> bool:
    if not np.array_equal(left.rotation, right.rotation):
        return False
    difference = left.translation - right.translation
    difference -= np.rint(difference)
    return bool(np.max(np.abs(difference)) < _TRANSLATION_TOLERANCE)


def closure(generators: Iterable[SeitzOp]) -> list[SeitzOp]:
    """Generate the finite group from ``generators`` modulo lattice translations.

    The input operations must generate a finite space group; otherwise the
    loop grows without bound, so an internal membership cap guards against
    unbounded growth.
    """
    identity = SeitzOp.identity()
    members: list[SeitzOp] = [identity]
    worklist = [identity]
    while worklist:
        current = worklist.pop()
        for generator in generators:
            product = multiply(current, generator)
            if not any(equivalent(product, member) for member in members):
                members.append(product)
                worklist.append(product)
                if len(members) > 8192:
                    raise ValueError(
                        "closure exceeded 8192 members; generators do not "
                        "define a finite space group"
                    )
    return members
