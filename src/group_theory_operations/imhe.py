"""Planar intrinsic/extrinsic magnetoelectric Hall tensor selection rules."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from .catalog import GroupDataError, load_database
from .magnetic_layer_groups import (
    get_magnetic_layer_group,
    load_magnetic_layer_group_registry,
)
from .magnetic_point_groups import (
    get_magnetic_point_group,
    load_magnetic_point_group_registry,
    magnetic_point_group_operations,
)
from .representations import determinant3

Matrix2 = tuple[tuple[float, float], tuple[float, float]]
IMHE_SECTORS = ("intrinsic", "extrinsic")

_IMHE_ALIASES = {
    "intrinsic": "intrinsic",
    "in": "intrinsic",
    "imhe_intrinsic": "intrinsic",
    "extrinsic": "extrinsic",
    "ex": "extrinsic",
    "imhe_extrinsic": "extrinsic",
}


def canonical_imhe_sector(value: str) -> str:
    """Normalize the intrinsic/extrinsic IMHE sector name."""
    if not isinstance(value, str):
        raise GroupDataError("IMHE sector must be a string")
    key = value.strip().lower().replace("-", "_").replace(" ", "_")
    try:
        return _IMHE_ALIASES[key]
    except KeyError as exc:
        raise GroupDataError(
            "unknown IMHE sector; choices: intrinsic, extrinsic"
        ) from exc


def _clean(value: float, tolerance: float) -> float:
    if abs(value) <= tolerance:
        return 0.0
    nearest = round(value)
    if abs(value - nearest) <= tolerance:
        return float(nearest)
    return float(value)


def _nullspace(rows: list[list[float]], columns: int, tolerance: float) -> tuple[tuple[float, ...], ...]:
    if not rows:
        return tuple(
            tuple(1.0 if index == column else 0.0 for index in range(columns))
            for column in range(columns)
        )
    matrix = [list(row) for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        if pivot_row >= len(matrix):
            break
        candidate = max(
            range(pivot_row, len(matrix)),
            key=lambda row: abs(matrix[row][column]),
        )
        if abs(matrix[candidate][column]) <= tolerance:
            continue
        matrix[pivot_row], matrix[candidate] = matrix[candidate], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            factor = matrix[row][column]
            if abs(factor) <= tolerance:
                continue
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_columns.append(column)
        pivot_row += 1

    free_columns = [column for column in range(columns) if column not in pivot_columns]
    basis = []
    for free_column in free_columns:
        vector = [0.0] * columns
        vector[free_column] = 1.0
        for row, pivot_column in reversed(list(enumerate(pivot_columns))):
            vector[pivot_column] = -sum(
                matrix[row][column] * vector[column]
                for column in free_columns
            )
        basis.append(tuple(_clean(value, tolerance * 10.0) for value in vector))
    return tuple(basis)


def _matmul(left: Matrix2, right: Matrix2) -> Matrix2:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def _subtract(left: Matrix2, right: Matrix2) -> Matrix2:
    return (
        (left[0][0] - right[0][0], left[0][1] - right[0][1]),
        (left[1][0] - right[1][0], left[1][1] - right[1][1]),
    )


def _scale(matrix: Matrix2, factor: float) -> Matrix2:
    return (
        (factor * matrix[0][0], factor * matrix[0][1]),
        (factor * matrix[1][0], factor * matrix[1][1]),
    )


def _linear_combination(coefficients: Sequence[float], basis: Sequence[Matrix2]) -> Matrix2:
    values = [[0.0, 0.0], [0.0, 0.0]]
    for coefficient, matrix in zip(coefficients, basis, strict=True):
        for row in range(2):
            for column in range(2):
                values[row][column] += coefficient * matrix[row][column]
    return (
        (values[0][0], values[0][1]),
        (values[1][0], values[1][1]),
    )


def _planar_block(
    matrix: Sequence[Sequence[Any]], *, tolerance: float
) -> tuple[Matrix2, float]:
    if (
        not isinstance(matrix, Sequence)
        or len(matrix) != 3
        or any(not isinstance(row, Sequence) or len(row) != 3 for row in matrix)
    ):
        raise GroupDataError("IMHE spatial operation must be a 3x3 matrix")
    values = tuple(tuple(float(item) for item in row) for row in matrix)
    if any(not math.isfinite(item) for row in values for item in row):
        raise GroupDataError("IMHE spatial operation contains non-finite values")

    mixing = (
        values[0][2],
        values[1][2],
        values[2][0],
        values[2][1],
    )
    if any(abs(value) > tolerance for value in mixing) or abs(abs(values[2][2]) - 1.0) > tolerance:
        raise GroupDataError(
            "fixed-zz planar IMHE is not closed for operations that mix the chosen z axis with x/y"
        )
    planar: Matrix2 = (
        (values[0][0], values[0][1]),
        (values[1][0], values[1][1]),
    )
    return planar, determinant3(values)


def _candidate_basis(sector: str) -> tuple[Matrix2, ...]:
    if sector == "intrinsic":
        return (((0.0, -1.0), (1.0, 0.0)),)
    return (
        ((1.0, 0.0), (0.0, 0.0)),
        ((0.0, 1.0), (1.0, 0.0)),
        ((0.0, 0.0), (0.0, 1.0)),
    )


def _solve_imhe_operations(
    operations: Sequence[tuple[Sequence[Sequence[Any]], bool]],
    sector: str,
    *,
    tolerance: float,
) -> tuple[Matrix2, ...]:
    resolved_sector = canonical_imhe_sector(sector)
    if not operations:
        raise ValueError("IMHE operation list cannot be empty")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be a positive finite number")

    eta_t = 1.0 if resolved_sector == "intrinsic" else -1.0
    candidates = _candidate_basis(resolved_sector)
    constraints: list[list[float]] = []

    for matrix, time_reversal in operations:
        if type(time_reversal) is not bool:
            raise ValueError("time-reversal labels must be boolean")
        planar, determinant = _planar_block(matrix, tolerance=tolerance)
        temporal = eta_t if time_reversal else 1.0
        factor = temporal * determinant
        residuals = tuple(
            _subtract(_matmul(planar, basis), _scale(_matmul(basis, planar), factor))
            for basis in candidates
        )
        for row in range(2):
            for column in range(2):
                equation = [residual[row][column] for residual in residuals]
                if any(abs(value) > tolerance for value in equation):
                    constraints.append(equation)

    coefficient_basis = _nullspace(constraints, len(candidates), tolerance)
    result = []
    for coefficients in coefficient_basis:
        matrix = _linear_combination(coefficients, candidates)
        result.append(
            tuple(
                tuple(_clean(value, tolerance * 10.0) for value in row)
                for row in matrix
            )
        )
    return tuple(result)


@dataclass(frozen=True, slots=True)
class IMHETensorBasis:
    """Allowed planar IMHE coefficient matrices for one magnetic group."""

    symmetry_class: str
    group_number: int
    group_identifier: str
    group_symbol: str
    magnetic_type: str
    sector: str
    time_character: str
    coefficient_symmetry: str
    basis: tuple[Matrix2, ...]

    @property
    def dimension(self) -> int:
        return len(self.basis)

    @property
    def shape(self) -> tuple[int, int]:
        return (2, 2)

    @property
    def allowed(self) -> bool:
        return self.dimension > 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "symmetry_class": self.symmetry_class,
            "group_number": self.group_number,
            "group_identifier": self.group_identifier,
            "group_symbol": self.group_symbol,
            "magnetic_type": self.magnetic_type,
            "response": f"imhe_{self.sector}",
            "sector": self.sector,
            "time_character": self.time_character,
            "coefficient_symmetry": self.coefficient_symmetry,
            "equation": "j_a = chi_ab;zz E_b E_z B_z; a,b in {x,y}",
            "shape": [2, 2],
            "dimension": self.dimension,
            "allowed": self.allowed,
            "row_basis": ["j_x", "j_y"],
            "column_basis": ["E_x", "E_y"],
            "basis": [[list(row) for row in matrix] for matrix in self.basis],
        }


def magnetic_imhe_tensor_basis(
    magnetic_point_group: str | int,
    sector: str,
    *,
    database: Mapping[str, Any] | None = None,
    registry: Mapping[str, Any] | None = None,
    tolerance: float = 1e-10,
) -> IMHETensorBasis:
    """Return planar IMHE selection rules for one magnetic point group."""
    resolved_sector = canonical_imhe_sector(sector)
    source_database = load_database() if database is None else database
    source_registry = (
        load_magnetic_point_group_registry() if registry is None else registry
    )
    group = get_magnetic_point_group(magnetic_point_group, source_registry)
    operations = magnetic_point_group_operations(
        group.number,
        database=source_database,
        registry=source_registry,
    )
    basis = _solve_imhe_operations(
        tuple(
            (operation.spatial.matrix_cartesian, operation.time_reversal)
            for operation in operations
        ),
        resolved_sector,
        tolerance=tolerance,
    )
    return IMHETensorBasis(
        symmetry_class="magnetic_point_group",
        group_number=group.number,
        group_identifier=group.magnetic_number,
        group_symbol=group.hm_symbol,
        magnetic_type=group.category,
        sector=resolved_sector,
        time_character="even" if resolved_sector == "intrinsic" else "odd",
        coefficient_symmetry=(
            "antisymmetric" if resolved_sector == "intrinsic" else "symmetric"
        ),
        basis=basis,
    )


def magnetic_layer_imhe_tensor_basis(
    magnetic_layer_group: str | int,
    sector: str,
    *,
    registry: Mapping[str, Any] | None = None,
    tolerance: float = 1e-10,
) -> IMHETensorBasis:
    """Return planar IMHE selection rules for one of the 528 magnetic layer groups."""
    resolved_sector = canonical_imhe_sector(sector)
    source_registry = (
        load_magnetic_layer_group_registry() if registry is None else registry
    )
    group = get_magnetic_layer_group(magnetic_layer_group, source_registry)
    basis = _solve_imhe_operations(
        tuple(
            (operation.matrix_cartesian, operation.time_reversal)
            for operation in group.point_operations
        ),
        resolved_sector,
        tolerance=tolerance,
    )
    return IMHETensorBasis(
        symmetry_class="magnetic_layer_group",
        group_number=group.global_number,
        group_identifier=group.og_number,
        group_symbol=group.litvin_og_symbol_ascii,
        magnetic_type=group.magnetic_type,
        sector=resolved_sector,
        time_character="even" if resolved_sector == "intrinsic" else "odd",
        coefficient_symmetry=(
            "antisymmetric" if resolved_sector == "intrinsic" else "symmetric"
        ),
        basis=basis,
    )
