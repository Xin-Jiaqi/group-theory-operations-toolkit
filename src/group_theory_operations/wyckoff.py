"""Wyckoff-orbit and site-stabilizer analysis for periodic structures."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .catalog import GroupDataError
from .seitz import SeitzOp
from .structure_symmetry import (
    StructureSymmetryContext,
    Vector3,
    _validated_structure_arrays,
    classify_structure_symmetry,
)


def _periodic_cartesian_distance(left: Any, right: Any, lattice: Any) -> float:
    difference = np.asarray(left, dtype=np.float64) - np.asarray(
        right, dtype=np.float64
    )
    difference -= np.rint(difference)
    return float(np.linalg.norm(difference @ np.asarray(lattice, dtype=np.float64)))


def _nearest_standard_site(
    position: Any,
    species: str,
    context: StructureSymmetryContext,
) -> int:
    candidates = [
        index
        for index, symbol in enumerate(context.standardized_species)
        if symbol == species
    ]
    if not candidates:  # pragma: no cover - protected by structure classification
        raise GroupDataError("standardized structure has lost an atomic species")
    return min(
        candidates,
        key=lambda index: _periodic_cartesian_distance(
            position,
            context.standardized_fractional_coordinates[index],
            context.standardized_lattice,
        ),
    )


def _standard_site_stabilizer(
    position: Any,
    context: StructureSymmetryContext,
    tolerance: float,
) -> tuple[SeitzOp, ...]:
    return tuple(
        operation
        for operation in context.standard_operations
        if _periodic_cartesian_distance(
            operation.apply(np.asarray(position, dtype=np.float64)),
            position,
            context.standardized_lattice,
        )
        <= tolerance
    )


def _standard_orbit_size(
    position: Any,
    context: StructureSymmetryContext,
    tolerance: float,
) -> int:
    images: list[np.ndarray] = []
    for operation in context.standard_operations:
        image = operation.apply(np.asarray(position, dtype=np.float64))
        if not any(
            _periodic_cartesian_distance(
                image, known, context.standardized_lattice
            )
            <= tolerance
            for known in images
        ):
            images.append(image)
    return len(images)


def _cartesian_rotation(rotation: Any, lattice: Any) -> np.ndarray:
    basis = np.asarray(lattice, dtype=np.float64).T
    return basis @ np.asarray(rotation, dtype=np.float64) @ np.linalg.inv(basis)


def _canonical_invariant_basis(
    stabilizer: tuple[SeitzOp, ...],
    lattice: Any,
    tolerance: float = 1.0e-9,
) -> tuple[Vector3, ...]:
    rotations = tuple(
        _cartesian_rotation(operation.rotation, lattice)
        for operation in stabilizer
    )
    projector = sum(rotations, np.zeros((3, 3), dtype=np.float64)) / len(rotations)
    projector = (projector + projector.T) / 2.0
    basis: list[np.ndarray] = []
    for axis in np.eye(3, dtype=np.float64):
        vector = projector @ axis
        for known in basis:
            vector -= np.dot(known, vector) * known
        norm = float(np.linalg.norm(vector))
        if norm <= tolerance:
            continue
        vector /= norm
        first = next(
            (value for value in vector if abs(float(value)) > tolerance), 1.0
        )
        if first < 0.0:
            vector *= -1.0
        basis.append(vector)
    for rotation in rotations:
        for vector in basis:
            if not np.allclose(rotation @ vector, vector, atol=1.0e-8, rtol=0.0):
                raise GroupDataError(
                    "site-stabilizer projector did not yield invariant displacements"
                )
    return tuple(
        (float(vector[0]), float(vector[1]), float(vector[2]))
        for vector in basis
    )


@dataclass(frozen=True, slots=True)
class WyckoffOrbit:
    """One occupied crystallographic orbit in a concrete structure.

    ``site_indices`` refer to the supplied cell.  The Wyckoff label,
    multiplicity, representative coordinate, stabilizer and displacement basis
    refer to the detected conventional standard setting.
    """

    species: str
    site_indices: tuple[int, ...]
    representative_site: int
    wyckoff_letter: str
    multiplicity: int
    site_symmetry_symbol: str
    standard_fractional_coordinate: Vector3
    stabilizer_operations: tuple[SeitzOp, ...]
    allowed_displacement_basis_cartesian: tuple[Vector3, ...]

    @property
    def label(self) -> str:
        return f"{self.multiplicity}{self.wyckoff_letter}"

    @property
    def input_orbit_size(self) -> int:
        """Return the occupied-orbit site count in the supplied cell."""

        return len(self.site_indices)

    @property
    def stabilizer_order(self) -> int:
        return len(self.stabilizer_operations)

    @property
    def allowed_displacement_dimension(self) -> int:
        return len(self.allowed_displacement_basis_cartesian)

    def to_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "site_indices": list(self.site_indices),
            "representative_site": self.representative_site,
            "wyckoff_letter": self.wyckoff_letter,
            "multiplicity": self.multiplicity,
            "label": self.label,
            "input_orbit_size": self.input_orbit_size,
            "site_symmetry_symbol": self.site_symmetry_symbol,
            "standard_fractional_coordinate": list(
                self.standard_fractional_coordinate
            ),
            "stabilizer_order": self.stabilizer_order,
            "stabilizer_operations": [
                operation.to_dict() for operation in self.stabilizer_operations
            ],
            "allowed_displacement_dimension": self.allowed_displacement_dimension,
            "allowed_displacement_basis_cartesian": [
                list(vector)
                for vector in self.allowed_displacement_basis_cartesian
            ],
        }


@dataclass(frozen=True, slots=True)
class WyckoffOrbitAnalysis:
    """Occupied Wyckoff orbits together with their structure classification."""

    symmetry: StructureSymmetryContext
    orbits: tuple[WyckoffOrbit, ...]

    def to_dict(self) -> dict[str, Any]:
        space_group = self.symmetry.space_group
        hall_setting = self.symmetry.hall_setting
        return {
            "space_group": {
                "ita_number": space_group.ita_number,
                "international_short": space_group.international_short,
                "international_full": space_group.international_full,
                "crystal_system": space_group.crystal_system,
            },
            "hall_setting": {
                "hall_number": hall_setting.hall_number,
                "hall_symbol": hall_setting.hall_symbol,
                "choice": hall_setting.choice,
                "centering": hall_setting.centering,
            },
            "standardized_lattice": [
                list(row) for row in self.symmetry.standardized_lattice
            ],
            "displacement_coordinate_system": (
                "Cartesian axes of the standardized lattice orientation"
            ),
            "orbits": [orbit.to_dict() for orbit in self.orbits],
            "symprec": self.symmetry.symprec,
            "angle_tolerance": self.symmetry.angle_tolerance,
            "backend": {
                "name": self.symmetry.backend_name,
                "version": self.symmetry.backend_version,
            },
        }


def analyze_wyckoff_orbits(
    structure: Any,
    *,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
) -> WyckoffOrbitAnalysis:
    """Resolve occupied Wyckoff orbits and symmetry-preserving displacements.

    The displacement basis is the fixed-vector space of the site stabilizer.
    It describes infinitesimal displacements that preserve that stabilizer; it
    is not a phonon eigenvector, force constant or energetic prediction.
    """

    _, species, coordinates, _ = _validated_structure_arrays(structure)
    context = classify_structure_symmetry(
        structure,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    tolerance = max(10.0 * float(symprec), 1.0e-8)
    # ``crystallographic_orbits`` follows the primitive crystallographic
    # symmetry and is therefore the correct partition for Wyckoff orbits.
    # ``equivalent_atoms`` can differ for unusual supercell presentations.
    labels = context.crystallographic_orbits
    orbit_members: dict[int, list[int]] = {}
    for index, label in enumerate(labels):
        orbit_members.setdefault(int(label), []).append(index)

    orbits: list[WyckoffOrbit] = []
    for indices in orbit_members.values():
        representative = min(indices)
        letters = {context.wyckoff_letters[index] for index in indices}
        symbols = {context.site_symmetry_symbols[index] for index in indices}
        orbit_species = {species[index] for index in indices}
        if len(letters) != 1 or len(symbols) != 1 or len(orbit_species) != 1:
            raise GroupDataError(
                "one occupied orbit has inconsistent Wyckoff metadata or species"
            )

        input_equivalent_indices = tuple(
            index
            for index, label in enumerate(context.equivalent_atoms)
            if label == context.equivalent_atoms[representative]
        )
        input_stabilizer_order = sum(
            mapping[representative] == representative
            for mapping in context.site_mappings
        )
        if (
            len(input_equivalent_indices) * input_stabilizer_order
            != len(context.input_operations)
        ):
            raise GroupDataError(
                "input-cell orbit and site stabilizer violate orbit-stabilizer"
            )

        transformed = context.to_standard_fractional([coordinates[representative]])[0]
        standard_index = _nearest_standard_site(
            transformed, species[representative], context
        )
        standard_position = context.standardized_fractional_coordinates[
            standard_index
        ]
        stabilizer = _standard_site_stabilizer(
            standard_position, context, tolerance
        )
        multiplicity = _standard_orbit_size(
            standard_position, context, tolerance
        )
        if not stabilizer or (
            multiplicity * len(stabilizer) != len(context.standard_operations)
        ):
            raise GroupDataError(
                "standard-cell orbit and site stabilizer violate orbit-stabilizer"
            )
        orbits.append(
            WyckoffOrbit(
                species=next(iter(orbit_species)),
                site_indices=tuple(indices),
                representative_site=representative,
                wyckoff_letter=next(iter(letters)),
                multiplicity=multiplicity,
                site_symmetry_symbol=next(iter(symbols)),
                standard_fractional_coordinate=standard_position,
                stabilizer_operations=stabilizer,
                allowed_displacement_basis_cartesian=_canonical_invariant_basis(
                    stabilizer, context.standardized_lattice
                ),
            )
        )
    return WyckoffOrbitAnalysis(symmetry=context, orbits=tuple(orbits))
