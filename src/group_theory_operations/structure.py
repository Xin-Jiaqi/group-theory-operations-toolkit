"""Bridge between catalog operations and materials-structure-core records."""

from __future__ import annotations

from typing import Any

from .catalog import GroupDataError, OperationRecord


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
    matrix = operation.matrix_fractional
    lattice = tuple(tuple(float(value) for value in row) for row in structure.lattice)
    if not _preserves_lattice_metric(lattice, matrix):
        raise GroupDataError(
            f"{operation.name} is incompatible with the structure lattice metric; "
            "the catalog operation is only physical in its documented family basis"
        )

    transformed_pbc: list[bool] = []
    for row in matrix:
        dependencies = {bool(structure.pbc[j]) for j, value in enumerate(row) if not _is_zero(value)}
        if len(dependencies) != 1:
            raise GroupDataError(
                f"{operation.name} mixes periodic and non-periodic coordinate directions"
            )
        transformed_pbc.append(dependencies.pop())
    if tuple(transformed_pbc) != tuple(structure.pbc):
        raise GroupDataError(f"{operation.name} changes the declared PBC axis semantics")

    fractional = [
        _vector(matrix, tuple(float(value) for value in coordinate))
        for coordinate in structure.fractional_coordinates
    ]
    dynamics = structure.selective_dynamics
    if dynamics is not None:
        if not _is_signed_permutation(matrix):
            raise GroupDataError(
                f"{operation.name} cannot exactly map component-wise Selective dynamics"
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
