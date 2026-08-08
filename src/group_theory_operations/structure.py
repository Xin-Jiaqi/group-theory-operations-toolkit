"""Bridge between catalog operations and materials-structure-core records."""

from __future__ import annotations

from typing import Any

import numpy as np

from .catalog import GroupDataError, OperationRecord
from .seitz import SeitzOp


def _vector(matrix: tuple[tuple[float, ...], ...], vector: tuple[float, ...]) -> list[float]:
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def _is_signed_permutation(matrix: tuple[tuple[float, ...], ...]) -> bool:
    rows = [sum(not _is_zero(value) for value in row) for row in matrix]
    columns = [sum(not _is_zero(matrix[i][j]) for i in range(3)) for j in range(3)]
    return rows == [1, 1, 1] and columns == [1, 1, 1] and all(
        _is_zero(value) or abs(abs(value) - 1.0) < 1.0e-12
        for row in matrix
        for value in row
    )


def _is_zero(value: float) -> bool:
    return abs(value) < 1.0e-12


def _transpose(matrix: tuple[tuple[float, ...], ...]) -> tuple[tuple[float, ...], ...]:
    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))


def _matmul(
    left: tuple[tuple[float, ...], ...],
    right: tuple[tuple[float, ...], ...],
) -> tuple[tuple[float, ...], ...]:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def _preserves_lattice_metric(
    lattice: tuple[tuple[float, ...], ...],
    operation: tuple[tuple[float, ...], ...],
) -> bool:
    """Return whether ``D.T @ (L @ L.T) @ D == L @ L.T``."""

    metric = _matmul(lattice, _transpose(lattice))
    transformed = _matmul(_matmul(_transpose(operation), metric), operation)
    scale = max(1.0, *(abs(value) for row in metric for value in row))
    return all(
        abs(transformed[i][j] - metric[i][j]) <= 1.0e-9 * scale
        for i in range(3)
        for j in range(3)
    )


def apply_fractional_operation(
    structure: Any,
    operation: OperationRecord,
    *,
    wrap: bool = True,
) -> Any:
    """Apply ``f' = D f`` about the origin while keeping the lattice fixed.

    ``structure`` must implement the public ``materials-structure-core``
    `StructureRecord` contract. Operations incompatible with PBC semantics or
    component-wise Selective dynamics are rejected rather than approximated.
    """

    return _apply_affine_fractional_operation(
        structure,
        matrix=operation.matrix_fractional,
        translation=(0.0, 0.0, 0.0),
        operation_name=operation.name,
        wrap=wrap,
    )


def _apply_seitz_operation(
    structure: Any,
    operation: SeitzOp,
    *,
    wrap: bool = True,
) -> Any:
    """Apply ``f' = R f + t`` without exposing a stable structure API yet."""

    matrix = tuple(
        tuple(float(value) for value in row) for row in operation.rotation
    )
    translation = tuple(float(value) for value in operation.translation)
    return _apply_affine_fractional_operation(
        structure,
        matrix=matrix,
        translation=translation,
        operation_name="Seitz operation",
        wrap=wrap,
    )


def _transform_fractional_coordinates(
    coordinates: Any,
    transformation_matrix: Any,
    origin_shift: Any,
    *,
    wrap: bool = True,
) -> tuple[tuple[float, float, float], ...]:
    """Transform points by ``x_new = P x_old + p``.

    This internal helper fixes the same modern spglib/ITA convention used by
    :func:`transform_seitz_coordinates`.  It intentionally transforms only
    coordinate representations; idealization and Cartesian rigid rotations
    of a standardized structure are separate operations.
    """

    points = np.asarray(coordinates, dtype=np.float64)
    matrix = np.asarray(transformation_matrix, dtype=np.float64)
    shift = np.asarray(origin_shift, dtype=np.float64)
    if points.ndim != 2 or points.shape[1:] != (3,):
        raise ValueError("coordinates must have shape (n, 3)")
    if matrix.shape != (3, 3):
        raise ValueError("transformation_matrix must be a 3x3 matrix")
    if shift.shape != (3,):
        raise ValueError("origin_shift must be a length-3 vector")
    if not (
        np.all(np.isfinite(points))
        and np.all(np.isfinite(matrix))
        and np.all(np.isfinite(shift))
    ):
        raise ValueError("coordinate transformation must contain finite values")
    if abs(float(np.linalg.det(matrix))) < 1.0e-12:
        raise ValueError("transformation_matrix must be invertible")
    transformed = points @ matrix.T + shift
    if wrap:
        transformed %= 1.0
        transformed[np.abs(transformed) < 1.0e-12] = 0.0
        transformed[np.abs(transformed - 1.0) < 1.0e-12] = 0.0
    return tuple(
        (float(row[0]), float(row[1]), float(row[2]))
        for row in transformed
    )


def _apply_affine_fractional_operation(
    structure: Any,
    *,
    matrix: tuple[tuple[float, ...], ...],
    translation: tuple[float, ...],
    operation_name: str,
    wrap: bool,
) -> Any:
    required = (
        "fractional_coordinates",
        "lattice",
        "species",
        "pbc",
        "selective_dynamics",
        "length_unit",
        "from_fractional",
        "wrapped",
    )
    if any(not hasattr(structure, attribute) for attribute in required):
        raise TypeError("structure does not implement the StructureRecord contract")
    lattice = tuple(tuple(float(value) for value in row) for row in structure.lattice)
    if not _preserves_lattice_metric(lattice, matrix):
        raise GroupDataError(
            f"{operation_name} is incompatible with the structure lattice metric; "
            "the fractional rotation requires a compatible lattice basis"
        )

    transformed_pbc: list[bool] = []
    for row in matrix:
        dependencies = {bool(structure.pbc[j]) for j, value in enumerate(row) if not _is_zero(value)}
        if len(dependencies) != 1:
            raise GroupDataError(
                f"{operation_name} mixes periodic and non-periodic coordinate directions"
            )
        transformed_pbc.append(dependencies.pop())
    if tuple(transformed_pbc) != tuple(structure.pbc):
        raise GroupDataError(
            f"{operation_name} changes the declared PBC axis semantics"
        )

    fractional = []
    for coordinate in structure.fractional_coordinates:
        rotated = _vector(matrix, tuple(float(value) for value in coordinate))
        fractional.append(
            [rotated[index] + translation[index] for index in range(3)]
        )
    dynamics = structure.selective_dynamics
    if dynamics is not None:
        if not _is_signed_permutation(matrix):
            raise GroupDataError(
                f"{operation_name} cannot exactly map component-wise Selective dynamics"
            )
        dynamics = [
            tuple(
                any(not _is_zero(matrix[i][j]) and allowed[j] for j in range(3))
                for i in range(3)
            )
            for allowed in dynamics
        ]
    result = structure.from_fractional(
        lattice=structure.lattice,
        species=structure.species,
        fractional_coordinates=fractional,
        pbc=structure.pbc,
        selective_dynamics=dynamics,
        length_unit=structure.length_unit,
    )
    return result.wrapped() if wrap else result


def _seitz_site_mapping(
    structure: Any,
    operation: SeitzOp,
    *,
    tolerance: float = 1.0e-5,
) -> tuple[int, ...] | None:
    """Return the same-species site permutation induced by ``operation``.

    Periodic images are compared with the full Cartesian lattice metric.
    ``None`` means that the geometrically valid operation does not map the
    decorated structure onto itself.  A perfect bipartite matching is used so
    that several sites inside the tolerance cannot collapse onto one target.
    """

    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be a positive finite length")
    required = ("fractional_coordinates", "lattice", "species", "pbc")
    if any(not hasattr(structure, attribute) for attribute in required):
        raise TypeError("structure does not implement the StructureRecord contract")

    species = tuple(structure.species)
    lattice = np.asarray(structure.lattice, dtype=np.float64)
    coordinates = np.asarray(structure.fractional_coordinates, dtype=np.float64)
    pbc = tuple(bool(value) for value in structure.pbc)
    if lattice.shape != (3, 3):
        raise GroupDataError("structure lattice must be a 3x3 matrix")
    if coordinates.shape != (len(species), 3) or not species:
        raise GroupDataError(
            "structure needs matching non-empty species and fractional coordinates"
        )
    if len(pbc) != 3:
        raise GroupDataError("structure pbc must contain three axis flags")
    if not np.all(np.isfinite(lattice)) or not np.all(np.isfinite(coordinates)):
        raise GroupDataError("structure lattice and coordinates must be finite")

    images = coordinates @ operation.rotation.T + operation.translation
    candidates: list[list[int]] = []
    for source_index, image in enumerate(images):
        distances: list[tuple[float, int]] = []
        for target_index, target in enumerate(coordinates):
            if species[source_index] != species[target_index]:
                continue
            difference = image - target
            for axis, periodic in enumerate(pbc):
                if periodic:
                    difference[axis] -= np.rint(difference[axis])
            distance = float(np.linalg.norm(difference @ lattice))
            if distance <= tolerance:
                distances.append((distance, target_index))
        candidates.append([target for _, target in sorted(distances)])
    if any(not options for options in candidates):
        return None

    target_to_source: list[int | None] = [None] * len(species)
    source_to_target = [-1] * len(species)

    def augment(source: int, visited: set[int]) -> bool:
        for target in candidates[source]:
            if target in visited:
                continue
            visited.add(target)
            previous = target_to_source[target]
            if previous is None or augment(previous, visited):
                target_to_source[target] = source
                source_to_target[source] = target
                return True
        return False

    for source in range(len(species)):
        if not augment(source, set()):
            return None
    return tuple(source_to_target)


def _site_orbits(site_mappings: Any) -> tuple[int, ...]:
    """Return minimum-index orbit representatives for site permutations."""

    mappings = tuple(tuple(mapping) for mapping in site_mappings)
    if not mappings:
        raise ValueError("at least one site mapping is required")
    site_count = len(mappings[0])
    expected = set(range(site_count))
    if any(
        len(mapping) != site_count
        or any(type(value) is not int for value in mapping)
        or set(mapping) != expected
        for mapping in mappings
    ):
        raise ValueError("site mappings must be equal-length permutations")

    parent = list(range(site_count))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for mapping in mappings:
        for source, target in enumerate(mapping):
            union(source, target)
    members: dict[int, list[int]] = {}
    for index in range(site_count):
        members.setdefault(find(index), []).append(index)
    representative = {
        root: min(indices) for root, indices in members.items()
    }
    return tuple(representative[find(index)] for index in range(site_count))
