"""Induced representations on quadratic optical-field spaces."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


Matrix3 = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]
Matrix6 = tuple[
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
]


def _matrix3(values: Sequence[Sequence[Any]]) -> Matrix3:
    if (
        isinstance(values, (str, bytes))
        or not isinstance(values, Sequence)
        or len(values) != 3
        or any(
            isinstance(row, (str, bytes))
            or not isinstance(row, Sequence)
            or len(row) != 3
            for row in values
        )
    ):
        raise ValueError("matrix must have shape (3, 3)")
    if any(type(value) not in (int, float) for row in values for value in row):
        raise ValueError("matrix must contain finite real numbers")
    matrix = tuple(tuple(float(value) for value in row) for row in values)
    if not all(math.isfinite(value) for row in matrix for value in row):
        raise ValueError("matrix must contain finite real numbers")
    return matrix  # type: ignore[return-value]


def determinant3(matrix: Sequence[Sequence[Any]]) -> float:
    """Return the determinant of a finite real 3 x 3 matrix."""

    rows = _matrix3(matrix)
    return (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )


def _clean(value: float, *, tolerance: float = 1e-12) -> float:
    if abs(value) < tolerance:
        return 0.0
    nearest = round(value)
    if abs(value - nearest) < tolerance:
        return float(nearest)
    return value


def symmetric_field_matrix(
    matrix: Sequence[Sequence[Any]],
    *,
    normalization: str = "sum",
) -> Matrix6:
    r"""Return :math:`M_+(R)=\operatorname{Sym}^2D(R)`.

    ``normalization="sum"`` uses the basis

    ``(|Ex|², |Ey|², |Ez|², Ex Ey* + Ey Ex*, Ex Ez* + Ez Ex*, Ey Ez* + Ez Ey*)``.

    ``normalization="half"`` uses one half of each off-diagonal sum and is
    retained for the historical human-readable tables.
    """

    if normalization not in {"sum", "half"}:
        raise ValueError("normalization must be 'sum' or 'half'")
    rows = _matrix3(matrix)
    result: list[tuple[float, ...]] = []
    diagonal_to_off_diagonal = 1.0 if normalization == "sum" else 2.0
    off_diagonal_to_diagonal = 2.0 if normalization == "sum" else 1.0

    for row in rows:
        x, y, z = row
        result.append(
            (
                x * x,
                y * y,
                z * z,
                diagonal_to_off_diagonal * x * y,
                diagonal_to_off_diagonal * x * z,
                diagonal_to_off_diagonal * y * z,
            )
        )
    for first, second in ((rows[0], rows[1]), (rows[0], rows[2]), (rows[1], rows[2])):
        result.append(
            (
                off_diagonal_to_diagonal * first[0] * second[0],
                off_diagonal_to_diagonal * first[1] * second[1],
                off_diagonal_to_diagonal * first[2] * second[2],
                first[0] * second[1] + first[1] * second[0],
                first[0] * second[2] + first[2] * second[0],
                first[1] * second[2] + first[2] * second[1],
            )
        )
    return tuple(
        tuple(_clean(value) for value in row) for row in result
    )  # type: ignore[return-value]


def antisymmetric_field_matrix(matrix: Sequence[Sequence[Any]]) -> Matrix3:
    r"""Return the axial-vector representation :math:`M_-(R)=\det(D)D`."""

    rows = _matrix3(matrix)
    determinant = determinant3(rows)
    return tuple(
        tuple(_clean(determinant * value) for value in row) for row in rows
    )  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class QuadraticFieldRepresentation:
    """The symmetric and antisymmetric quadratic-field actions of one operation."""

    determinant: float
    matrix_symmetric: Matrix6
    matrix_antisymmetric: Matrix3

    def to_dict(self) -> dict[str, Any]:
        return {
            "determinant": _clean(self.determinant),
            "matrix_symmetric": [list(row) for row in self.matrix_symmetric],
            "matrix_antisymmetric": [list(row) for row in self.matrix_antisymmetric],
        }


def quadratic_field_representation(
    operation_or_matrix: Mapping[str, Any] | Sequence[Sequence[Any]],
    *,
    normalization: str = "sum",
) -> QuadraticFieldRepresentation:
    """Build both induced representations from an operation or Cartesian matrix."""

    matrix = (
        operation_or_matrix["matrix_cartesian"]
        if isinstance(operation_or_matrix, Mapping)
        else operation_or_matrix
    )
    rows = _matrix3(matrix)
    return QuadraticFieldRepresentation(
        determinant=_clean(determinant3(rows)),
        matrix_symmetric=symmetric_field_matrix(rows, normalization=normalization),
        matrix_antisymmetric=antisymmetric_field_matrix(rows),
    )


def load_quadratic_field_catalog(path: str | Path | None = None) -> dict[str, Any]:
    """Load the packaged, reproducibly generated all-operation field catalog."""

    if path is None:
        repository_copy = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "quadratic_field_representations.json"
        )
        text = (
            repository_copy.read_text(encoding="utf-8")
            if repository_copy.is_file()
            else resources.files("group_theory_operations")
            .joinpath("data/quadratic_field_representations.json")
            .read_text(encoding="utf-8")
        )
    else:
        text = Path(path).read_text(encoding="utf-8")
    data = json.loads(text)
    if data.get("schema_version") != 1:
        raise ValueError("quadratic-field catalog schema_version must be 1")
    return data
