"""Typed real-structure symmetry classification and setting conversion."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
from typing import Any, Iterable

import numpy as np

from .catalog import GroupDataError
from .seitz import SeitzOp, closure, equivalent, transform_seitz_coordinates
from .space_groups import (
    CrystallographicSpaceGroup,
    HallSetting,
    get_crystallographic_space_group,
)
from .structure import (
    _seitz_site_mapping,
    _site_orbits,
    _transform_fractional_coordinates,
)


Vector3 = tuple[float, float, float]
Matrix3 = tuple[Vector3, Vector3, Vector3]


def _vector3(values: Any) -> Vector3:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise GroupDataError("expected a finite length-3 vector")
    return (float(array[0]), float(array[1]), float(array[2]))


def _matrix3(values: Any) -> Matrix3:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (3, 3) or not np.all(np.isfinite(array)):
        raise GroupDataError("expected a finite 3x3 matrix")
    return (
        (float(array[0, 0]), float(array[0, 1]), float(array[0, 2])),
        (float(array[1, 0]), float(array[1, 1]), float(array[1, 2])),
        (float(array[2, 0]), float(array[2, 1]), float(array[2, 2])),
    )


def _unique_operations(operations: Iterable[SeitzOp]) -> tuple[SeitzOp, ...]:
    unique: list[SeitzOp] = []
    for operation in operations:
        if not any(equivalent(operation, known) for known in unique):
            unique.append(operation)
    return tuple(unique)


def _operation_subset(left: Iterable[SeitzOp], right: Iterable[SeitzOp]) -> bool:
    left_unique = _unique_operations(left)
    right_unique = _unique_operations(right)
    return all(
        any(equivalent(operation, candidate) for candidate in right_unique)
        for operation in left_unique
    )


def _partition(labels: Iterable[int]) -> set[frozenset[int]]:
    groups: dict[int, set[int]] = {}
    for index, label in enumerate(labels):
        groups.setdefault(int(label), set()).add(index)
    return {frozenset(indices) for indices in groups.values()}


def _maximum_same_type_periodic_distance(
    source_positions: Any,
    source_types: Any,
    target_positions: Any,
    target_types: Any,
    lattice: Any,
) -> float:
    sources = np.asarray(source_positions, dtype=np.float64)
    source_labels = np.asarray(source_types)
    targets = np.asarray(target_positions, dtype=np.float64)
    target_labels = np.asarray(target_types)
    metric_lattice = np.asarray(lattice, dtype=np.float64)
    maximum = 0.0
    for position, type_number in zip(sources, source_labels):
        candidates = targets[target_labels == type_number]
        if not len(candidates):
            return float("inf")
        difference = position - candidates
        difference -= np.rint(difference)
        maximum = max(
            maximum,
            float(np.min(np.linalg.norm(difference @ metric_lattice, axis=1))),
        )
    return maximum


@dataclass(frozen=True, slots=True)
class StructureSymmetryContext:
    """Immutable classification tying input and standard settings together.

    ``input_operations`` act in the supplied structure coordinates, while
    ``standard_operations`` are the checked-in operations of the detected Hall
    setting.  They need not have equal cardinality for a supercell or primitive
    input.  The coordinate map follows ``x_std = P x_input + p``.
    """

    space_group: CrystallographicSpaceGroup
    hall_setting: HallSetting
    input_operations: tuple[SeitzOp, ...]
    standard_operations: tuple[SeitzOp, ...]
    site_mappings: tuple[tuple[int, ...], ...]
    equivalent_atoms: tuple[int, ...]
    wyckoff_letters: tuple[str, ...]
    site_symmetry_symbols: tuple[str, ...]
    transformation_matrix: Matrix3
    origin_shift: Vector3
    standardization_rotation: Matrix3
    standardized_lattice: Matrix3
    standardized_species: tuple[str, ...]
    standardized_fractional_coordinates: tuple[Vector3, ...]
    symprec: float
    angle_tolerance: float
    backend_name: str
    backend_version: str

    @property
    def cell_volume_ratio(self) -> float:
        """Return the pre-idealization basis-volume ratio ``|det(P)|``."""

        return abs(float(np.linalg.det(np.asarray(self.transformation_matrix))))

    @property
    def operation_count_ratio(self) -> float:
        """Return input-operation count divided by standard-operation count."""

        return len(self.input_operations) / len(self.standard_operations)

    def to_standard_fractional(
        self,
        coordinates: Any,
        *,
        wrap: bool = True,
    ) -> tuple[Vector3, ...]:
        """Express input fractional coordinates in the standard setting.

        With ``wrap=True``, a supercell map is a quotient and is therefore not
        invertible site by site after reduction modulo the standard lattice.
        """

        return _transform_fractional_coordinates(
            coordinates,
            self.transformation_matrix,
            self.origin_shift,
            wrap=wrap,
        )

    def operation_to_standard(self, operation: SeitzOp) -> SeitzOp:
        """Express an input-setting Seitz operation in standard coordinates."""

        return transform_seitz_coordinates(
            operation,
            np.asarray(self.transformation_matrix),
            np.asarray(self.origin_shift),
        )


def _validated_structure_arrays(
    structure: Any,
) -> tuple[np.ndarray, tuple[str, ...], np.ndarray, np.ndarray]:
    required = ("lattice", "species", "fractional_coordinates", "pbc")
    if any(not hasattr(structure, attribute) for attribute in required):
        raise TypeError("structure does not implement the StructureRecord contract")
    try:
        lattice = np.asarray(structure.lattice, dtype=np.float64)
        coordinates = np.asarray(structure.fractional_coordinates, dtype=np.float64)
        species = tuple(structure.species)
        pbc = tuple(structure.pbc)
    except (TypeError, ValueError) as exc:
        raise GroupDataError("structure arrays cannot be interpreted") from exc
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise GroupDataError("structure lattice must be a finite 3x3 matrix")
    if abs(float(np.linalg.det(lattice))) < 1.0e-12:
        raise GroupDataError("structure lattice must be invertible")
    if coordinates.shape != (len(species), 3) or not species:
        raise GroupDataError(
            "structure needs matching non-empty species and fractional coordinates"
        )
    if not np.all(np.isfinite(coordinates)):
        raise GroupDataError("structure fractional coordinates must be finite")
    if any(not isinstance(symbol, str) or not symbol for symbol in species):
        raise GroupDataError("structure species must be non-empty strings")
    if len(pbc) != 3 or any(type(value) is not bool for value in pbc):
        raise GroupDataError("structure pbc must contain three boolean flags")
    if pbc != (True, True, True):
        raise GroupDataError(
            "crystallographic space-group classification requires three periodic axes"
        )

    unique_species: list[str] = []
    type_numbers: list[int] = []
    for symbol in species:
        if symbol not in unique_species:
            unique_species.append(symbol)
        type_numbers.append(unique_species.index(symbol) + 1)
    return lattice, species, coordinates, np.asarray(type_numbers, dtype=np.intc)


def classify_structure_symmetry(
    structure: Any,
    *,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
) -> StructureSymmetryContext:
    """Classify a three-periodic structure and return a validated context.

    spglib performs the numerical symmetry search.  The result is then checked
    against this package's Hall registry, species-aware site permutations, and
    the verified input-to-standard coordinate convention before it is exposed.
    Install the optional ``structure`` dependencies to use this function.
    """

    if not np.isfinite(symprec) or symprec <= 0.0:
        raise ValueError("symprec must be a positive finite length")
    if not np.isfinite(angle_tolerance) or (
        angle_tolerance != -1.0 and angle_tolerance <= 0.0
    ):
        raise ValueError("angle_tolerance must be -1 or a positive finite angle")
    lattice, species, coordinates, type_numbers = _validated_structure_arrays(
        structure
    )
    try:
        spglib = importlib.import_module("spglib")
    except ImportError as exc:
        raise ImportError(
            "structure symmetry classification requires spglib; install "
            "group-theory-operations-toolkit[structure]"
        ) from exc

    dataset = spglib.get_symmetry_dataset(
        (lattice, coordinates, type_numbers),
        symprec=float(symprec),
        angle_tolerance=float(angle_tolerance),
    )
    if dataset is None:
        raise GroupDataError("spglib could not classify the supplied structure")

    space_group = get_crystallographic_space_group(int(dataset.number))
    try:
        hall_setting = next(
            setting
            for setting in space_group.hall_settings
            if setting.hall_number == int(dataset.hall_number)
        )
    except StopIteration as exc:  # pragma: no cover - registry contract guard
        raise GroupDataError("detected Hall setting is absent from the registry") from exc

    input_operations = tuple(
        SeitzOp(rotation, translation)
        for rotation, translation in zip(dataset.rotations, dataset.translations)
    )
    mappings: list[tuple[int, ...]] = []
    mapping_tolerance = float(symprec) * (1.0 + 1.0e-9)
    for operation in input_operations:
        mapping = _seitz_site_mapping(
            structure,
            operation,
            tolerance=mapping_tolerance,
        )
        if mapping is None:
            raise GroupDataError(
                "a detected operation is not a species-aware structure automorphism"
            )
        mappings.append(mapping)
    equivalent_atoms = _site_orbits(mappings)
    if _partition(equivalent_atoms) != _partition(dataset.equivalent_atoms):
        raise GroupDataError(
            "site permutations disagree with the detected equivalent-atom partition"
        )

    matrix = np.asarray(dataset.transformation_matrix, dtype=np.float64)
    shift = np.asarray(dataset.origin_shift, dtype=np.float64)
    transformed_operations = _unique_operations(
        transform_seitz_coordinates(operation, matrix, shift)
        for operation in input_operations
    )
    standard_operations = tuple(
        closure(
            SeitzOp.from_dict(generator)
            for generator in hall_setting.generators
        )
    )
    inverse_matrix = np.linalg.inv(matrix)
    inverse_shift = -inverse_matrix @ shift
    folded_standard_operations = _unique_operations(
        transform_seitz_coordinates(
            operation,
            inverse_matrix,
            inverse_shift,
        )
        for operation in standard_operations
    )
    if not (
        _operation_subset(transformed_operations, standard_operations)
        and _operation_subset(folded_standard_operations, input_operations)
    ):
        raise GroupDataError(
            "input and registered Hall operations disagree under setting conversion"
        )
    volume_ratio = abs(float(np.linalg.det(matrix)))
    operation_ratio = len(input_operations) / len(standard_operations)
    if not np.isclose(volume_ratio, operation_ratio, atol=1.0e-8, rtol=0.0):
        raise GroupDataError(
            "input and standard operation counts disagree with the cell-volume ratio"
        )

    standard_lattice = np.asarray(dataset.std_lattice, dtype=np.float64)
    standard_rotation = np.asarray(dataset.std_rotation_matrix, dtype=np.float64)
    if not np.allclose(
        standard_rotation @ standard_rotation.T,
        np.eye(3),
        atol=1.0e-10,
        rtol=0.0,
    ):
        raise GroupDataError("standardization rotation is not orthonormal")
    predicted_standard_lattice = (
        np.linalg.inv(matrix).T @ lattice @ standard_rotation.T
    )
    validation_tolerance = 10.0 * float(symprec)
    if not np.allclose(
        predicted_standard_lattice,
        standard_lattice,
        atol=validation_tolerance,
        rtol=1.0e-10,
    ):
        raise GroupDataError(
            "standardized lattice is inconsistent with the basis and rigid rotation"
        )

    standardized_input_positions = np.asarray(
        _transform_fractional_coordinates(coordinates, matrix, shift)
    )
    standard_positions = np.asarray(dataset.std_positions, dtype=np.float64)
    if _maximum_same_type_periodic_distance(
        standardized_input_positions,
        type_numbers,
        standard_positions,
        dataset.std_types,
        standard_lattice,
    ) > validation_tolerance:
        raise GroupDataError(
            "input sites do not map to the idealized standardized structure"
        )
    orbit_positions: list[np.ndarray] = []
    orbit_types: list[int] = []
    for operation in standard_operations:
        for position, type_number in zip(
            standardized_input_positions,
            type_numbers,
        ):
            orbit_positions.append(operation.apply(position))
            orbit_types.append(int(type_number))
    if _maximum_same_type_periodic_distance(
        standard_positions,
        dataset.std_types,
        orbit_positions,
        orbit_types,
        standard_lattice,
    ) > validation_tolerance:
        raise GroupDataError(
            "standard Hall operations do not cover the standardized structure"
        )

    unique_species = tuple(dict.fromkeys(species))
    try:
        standard_species = tuple(
            unique_species[int(type_number) - 1]
            for type_number in dataset.std_types
        )
    except IndexError as exc:  # pragma: no cover - backend contract guard
        raise GroupDataError("standardized atom types are inconsistent") from exc

    return StructureSymmetryContext(
        space_group=space_group,
        hall_setting=hall_setting,
        input_operations=input_operations,
        standard_operations=standard_operations,
        site_mappings=tuple(mappings),
        equivalent_atoms=equivalent_atoms,
        wyckoff_letters=tuple(str(value) for value in dataset.wyckoffs),
        site_symmetry_symbols=tuple(
            str(value).strip() for value in dataset.site_symmetry_symbols
        ),
        transformation_matrix=_matrix3(matrix),
        origin_shift=_vector3(shift),
        standardization_rotation=_matrix3(standard_rotation),
        standardized_lattice=_matrix3(standard_lattice),
        standardized_species=standard_species,
        standardized_fractional_coordinates=tuple(
            _vector3(position) for position in standard_positions
        ),
        symprec=float(symprec),
        angle_tolerance=float(angle_tolerance),
        backend_name="spglib",
        backend_version=str(spglib.__version__),
    )
