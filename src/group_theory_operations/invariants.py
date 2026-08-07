"""Spatial point-group invariant bases for nonlinear optical response tensors."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from .catalog import GroupDataError, load_database
from .point_groups import (
    get_crystallographic_point_group,
    load_point_group_registry,
    point_group_operations,
)
from .magnetic_point_groups import (
    get_magnetic_point_group,
    load_magnetic_point_group_registry,
    magnetic_point_group_operations,
)
from .representations import determinant3, quadratic_field_representation


Matrix = tuple[tuple[float, ...], ...]

POLAR_BASIS = ("x", "y", "z")
SYMMETRIC_FIELD_BASIS = ("xx", "yy", "zz", "xy", "xz", "yz")
ANTISYMMETRIC_FIELD_BASIS = ("h_x", "h_y", "h_z")

TENSOR_SPACE_BASES = {
    "scalar": ("1",),
    "pseudoscalar": ("p",),
    "polar_vector": POLAR_BASIS,
    "axial_vector": ("a_x", "a_y", "a_z"),
    "symmetric_quadratic": SYMMETRIC_FIELD_BASIS,
    "antisymmetric_quadratic": ANTISYMMETRIC_FIELD_BASIS,
}

_TENSOR_SPACE_ALIASES = {
    "scalar": "scalar",
    "pseudoscalar": "pseudoscalar",
    "polar": "polar_vector",
    "polar_vector": "polar_vector",
    "axial": "axial_vector",
    "axial_vector": "axial_vector",
    "symmetric": "symmetric_quadratic",
    "symmetric_quadratic": "symmetric_quadratic",
    "antisymmetric": "antisymmetric_quadratic",
    "antisymmetric_quadratic": "antisymmetric_quadratic",
}

RESPONSE_SPECS = {
    "shift_current": {
        "input_space": "symmetric",
        "output_basis": POLAR_BASIS,
        "input_basis": SYMMETRIC_FIELD_BASIS,
        "equation": "J^a(0) = sigma^{a;bc} E^b E^{c*}; sigma^{a;bc}=sigma^{a;cb}",
    },
    "shg": {
        "input_space": "symmetric",
        "output_basis": POLAR_BASIS,
        "input_basis": SYMMETRIC_FIELD_BASIS,
        "equation": "P^a(2omega) = chi^{a;bc} E^b(omega) E^c(omega); chi^{a;bc}=chi^{a;cb}",
    },
    "circular_injection_current": {
        "input_space": "antisymmetric",
        "output_basis": POLAR_BASIS,
        "input_basis": ANTISYMMETRIC_FIELD_BASIS,
        "equation": "dJ^a/dt = beta^{a;j} h_j; h = i E x E*",
    },
}

_RESPONSE_ALIASES = {
    "shift": "shift_current",
    "shiftcurrent": "shift_current",
    "shift_current": "shift_current",
    "shg": "shg",
    "secondharmonicgeneration": "shg",
    "circularinjection": "circular_injection_current",
    "circularinjectioncurrent": "circular_injection_current",
    "cpge": "circular_injection_current",
    "circular_injection_current": "circular_injection_current",
}


def canonical_response_name(value: str) -> str:
    """Normalize supported response aliases."""

    if not isinstance(value, str):
        raise GroupDataError("response name must be a string")
    key = value.strip().lower().replace("-", "_").replace(" ", "_")
    compact = key.replace("_", "")
    for candidate in (key, compact):
        if candidate in _RESPONSE_ALIASES:
            return _RESPONSE_ALIASES[candidate]
    choices = ", ".join(RESPONSE_SPECS)
    raise GroupDataError(f"unknown response {value!r}; choices: {choices}")


def _matrix(values: Sequence[Sequence[Any]], *, square: bool = False) -> Matrix:
    if (
        isinstance(values, (str, bytes))
        or not isinstance(values, Sequence)
        or not values
        or any(
            isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row
            for row in values
        )
    ):
        raise ValueError("matrix must be a non-empty rectangular sequence")
    width = len(values[0])
    if any(len(row) != width for row in values):
        raise ValueError("matrix must be rectangular")
    if square and len(values) != width:
        raise ValueError("representation matrices must be square")
    if any(type(item) not in (int, float) for row in values for item in row):
        raise ValueError("matrix must contain finite real numbers")
    result = tuple(tuple(float(item) for item in row) for row in values)
    if not all(math.isfinite(item) for row in result for item in row):
        raise ValueError("matrix must contain finite real numbers")
    return result


def _clean(value: float, tolerance: float = 1e-10) -> float:
    if abs(value) < tolerance:
        return 0.0
    nearest = round(value)
    if abs(value - nearest) < tolerance:
        return float(nearest)
    return float(value)


def _json_number(value: float) -> int | float:
    cleaned = _clean(value)
    nearest = round(cleaned)
    if cleaned == nearest:
        return int(nearest)
    return round(cleaned, 12)


def _nullspace(rows: list[list[float]], columns: int, tolerance: float) -> tuple[tuple[float, ...], ...]:
    if not rows:
        return tuple(
            tuple(1.0 if row == column else 0.0 for row in range(columns))
            for column in range(columns)
        )
    matrix = [list(row) for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        candidate = max(
            range(pivot_row, len(matrix)),
            key=lambda row: abs(matrix[row][column]),
            default=pivot_row,
        )
        if pivot_row >= len(matrix) or abs(matrix[candidate][column]) <= tolerance:
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
        if pivot_row == len(matrix):
            break
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
        basis.append(tuple(_clean(value, tolerance * 10) for value in vector))
    return tuple(basis)


def equivariant_map_basis(
    output_representations: Sequence[Sequence[Sequence[Any]]],
    input_representations: Sequence[Sequence[Sequence[Any]]],
    *,
    tolerance: float = 1e-10,
) -> tuple[Matrix, ...]:
    r"""Solve all maps :math:`T` satisfying :math:`A_gT=T B_g`."""

    if type(tolerance) not in (int, float) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be a positive finite number")
    if len(output_representations) != len(input_representations) or not output_representations:
        raise ValueError("input and output representation lists must have equal non-zero length")
    outputs = tuple(_matrix(value, square=True) for value in output_representations)
    inputs = tuple(_matrix(value, square=True) for value in input_representations)
    output_dimension = len(outputs[0])
    input_dimension = len(inputs[0])
    if any(len(value) != output_dimension for value in outputs):
        raise ValueError("all output representations must have the same shape")
    if any(len(value) != input_dimension for value in inputs):
        raise ValueError("all input representations must have the same shape")
    columns = output_dimension * input_dimension
    constraints: list[list[float]] = []
    for output, source in zip(outputs, inputs, strict=True):
        for row in range(output_dimension):
            for column in range(input_dimension):
                equation = [0.0] * columns
                for index in range(output_dimension):
                    equation[index * input_dimension + column] += output[row][index]
                for index in range(input_dimension):
                    equation[row * input_dimension + index] -= source[index][column]
                if any(abs(value) > tolerance for value in equation):
                    constraints.append(equation)
    vectors = _nullspace(constraints, columns, tolerance)
    return tuple(
        tuple(
            tuple(vector[row * input_dimension + column] for column in range(input_dimension))
            for row in range(output_dimension)
        )
        for vector in vectors
    )


@dataclass(frozen=True, slots=True)
class InvariantTensorBasis:
    """A deterministic spatial invariant basis for one optical response."""

    point_group_number: int
    point_group: str
    schoenflies_symbol: str
    response: str
    output_basis: tuple[str, ...]
    input_basis: tuple[str, ...]
    equation: str
    basis: tuple[Matrix, ...]

    @property
    def dimension(self) -> int:
        return len(self.basis)

    @property
    def shape(self) -> tuple[int, int]:
        return len(self.output_basis), len(self.input_basis)

    def to_dict(self) -> dict[str, Any]:
        return {
            "response": self.response,
            "shape": list(self.shape),
            "dimension": self.dimension,
            "output_basis": list(self.output_basis),
            "input_basis": list(self.input_basis),
            "equation": self.equation,
            "basis": [
                [[_json_number(value) for value in row] for row in matrix]
                for matrix in self.basis
            ],
        }


@dataclass(frozen=True, slots=True)
class MagneticInvariantTensorBasis:
    """A real tensor-map basis with explicit spatial and temporal parities.

    This representation applies to static or otherwise real response objects.
    Frequency-domain antiunitary constraints may additionally exchange
    frequencies and complex-conjugate coefficients, so they must be formulated
    separately rather than inferred from temporal parity alone.
    """

    magnetic_point_group_number: int
    magnetic_number: str
    magnetic_point_group: str
    category: str
    output_space: str
    input_space: str
    output_time_parity: str
    input_time_parity: str
    output_basis: tuple[str, ...]
    input_basis: tuple[str, ...]
    basis: tuple[Matrix, ...]

    @property
    def dimension(self) -> int:
        return len(self.basis)

    @property
    def shape(self) -> tuple[int, int]:
        return len(self.output_basis), len(self.input_basis)

    def to_dict(self) -> dict[str, Any]:
        return {
            "magnetic_point_group_number": self.magnetic_point_group_number,
            "magnetic_number": self.magnetic_number,
            "magnetic_point_group": self.magnetic_point_group,
            "category": self.category,
            "output_space": self.output_space,
            "input_space": self.input_space,
            "output_time_parity": self.output_time_parity,
            "input_time_parity": self.input_time_parity,
            "shape": list(self.shape),
            "dimension": self.dimension,
            "output_basis": list(self.output_basis),
            "input_basis": list(self.input_basis),
            "basis": [
                [[_json_number(value) for value in row] for row in matrix]
                for matrix in self.basis
            ],
        }


def canonical_tensor_space(value: str) -> str:
    """Normalize a supported real tensor representation name."""

    if not isinstance(value, str):
        raise GroupDataError("tensor space must be a string")
    key = value.strip().lower().replace("-", "_").replace(" ", "_")
    try:
        return _TENSOR_SPACE_ALIASES[key]
    except KeyError as exc:
        choices = ", ".join(TENSOR_SPACE_BASES)
        raise GroupDataError(f"unknown tensor space {value!r}; choices: {choices}") from exc


def canonical_time_parity(value: str | int) -> str:
    """Normalize temporal parity to ``even`` or ``odd``."""

    if type(value) is bool:
        raise GroupDataError("time parity must be 'even' (+1) or 'odd' (-1)")
    if value in (1, "+1", "+", "even"):
        return "even"
    if value in (-1, "-1", "-", "odd"):
        return "odd"
    raise GroupDataError("time parity must be 'even' (+1) or 'odd' (-1)")


def _scale_matrix(matrix: Matrix, factor: float) -> Matrix:
    return tuple(tuple(factor * value for value in row) for row in matrix)


def _spatial_representation(matrix: Sequence[Sequence[Any]], space: str) -> Matrix:
    resolved = _matrix(matrix, square=True)
    if len(resolved) != 3:
        raise ValueError("spatial point-operation matrices must be 3x3")
    if space == "scalar":
        return ((1.0,),)
    determinant = determinant3(resolved)
    if space == "pseudoscalar":
        return ((determinant,),)
    if space == "polar_vector":
        return resolved
    fields = quadratic_field_representation(resolved)
    if space == "axial_vector":
        return fields.matrix_antisymmetric
    if space == "symmetric_quadratic":
        return fields.matrix_symmetric
    if space == "antisymmetric_quadratic":
        return fields.matrix_antisymmetric
    raise AssertionError(f"unhandled tensor space {space}")


def magnetic_equivariant_map_basis(
    output_spatial_representations: Sequence[Sequence[Sequence[Any]]],
    input_spatial_representations: Sequence[Sequence[Sequence[Any]]],
    time_reversals: Sequence[bool],
    *,
    output_time_parity: str | int = "even",
    input_time_parity: str | int = "even",
    tolerance: float = 1e-10,
) -> tuple[Matrix, ...]:
    r"""Solve magnetic equivariance with explicit temporal parity.

    For operation :math:`(R,\theta)`, a time-odd object gains a factor
    :math:`(-1)^\theta`; a time-even object does not.  Spatial polar/axial
    behavior must already be represented by the supplied matrices.
    """

    if (
        len(output_spatial_representations) != len(time_reversals)
        or len(input_spatial_representations) != len(time_reversals)
        or not time_reversals
        or any(type(value) is not bool for value in time_reversals)
    ):
        raise ValueError("representations and boolean time-reversal labels must align")
    output_parity = canonical_time_parity(output_time_parity)
    input_parity = canonical_time_parity(input_time_parity)
    outputs = tuple(
        _scale_matrix(_matrix(matrix, square=True), -1.0)
        if time_reversal and output_parity == "odd"
        else _matrix(matrix, square=True)
        for matrix, time_reversal in zip(
            output_spatial_representations, time_reversals, strict=True
        )
    )
    inputs = tuple(
        _scale_matrix(_matrix(matrix, square=True), -1.0)
        if time_reversal and input_parity == "odd"
        else _matrix(matrix, square=True)
        for matrix, time_reversal in zip(
            input_spatial_representations, time_reversals, strict=True
        )
    )
    return equivariant_map_basis(outputs, inputs, tolerance=tolerance)


def magnetic_tensor_basis(
    magnetic_point_group: str | int,
    output_space: str,
    input_space: str,
    *,
    output_time_parity: str | int = "even",
    input_time_parity: str | int = "even",
    database: Mapping[str, Any] | None = None,
    registry: Mapping[str, Any] | None = None,
    tolerance: float = 1e-10,
) -> MagneticInvariantTensorBasis:
    """Return a real magnetic tensor-map basis for one magnetic point group."""

    source_database = load_database() if database is None else database
    source_registry = (
        load_magnetic_point_group_registry() if registry is None else registry
    )
    group = get_magnetic_point_group(magnetic_point_group, source_registry)
    output_name = canonical_tensor_space(output_space)
    input_name = canonical_tensor_space(input_space)
    output_parity = canonical_time_parity(output_time_parity)
    input_parity = canonical_time_parity(input_time_parity)
    operations = magnetic_point_group_operations(
        group.number,
        database=source_database,
        registry=source_registry,
    )
    output_representations = [
        _spatial_representation(operation.spatial.matrix_cartesian, output_name)
        for operation in operations
    ]
    input_representations = [
        _spatial_representation(operation.spatial.matrix_cartesian, input_name)
        for operation in operations
    ]
    basis = magnetic_equivariant_map_basis(
        output_representations,
        input_representations,
        [operation.time_reversal for operation in operations],
        output_time_parity=output_parity,
        input_time_parity=input_parity,
        tolerance=tolerance,
    )
    return MagneticInvariantTensorBasis(
        magnetic_point_group_number=group.number,
        magnetic_number=group.magnetic_number,
        magnetic_point_group=group.hm_symbol,
        category=group.category,
        output_space=output_name,
        input_space=input_name,
        output_time_parity=output_parity,
        input_time_parity=input_parity,
        output_basis=TENSOR_SPACE_BASES[output_name],
        input_basis=TENSOR_SPACE_BASES[input_name],
        basis=basis,
    )


def response_tensor_basis(
    point_group: str | int,
    response: str,
    *,
    database: Mapping[str, Any] | None = None,
    registry: Mapping[str, Any] | None = None,
    tolerance: float = 1e-10,
) -> InvariantTensorBasis:
    """Return spatially allowed tensor combinations for a standard point group."""

    source_database = load_database() if database is None else database
    source_registry = load_point_group_registry() if registry is None else registry
    group = get_crystallographic_point_group(point_group, source_registry)
    response_name = canonical_response_name(response)
    specification = RESPONSE_SPECS[response_name]
    operations = point_group_operations(
        group.number,
        database=source_database,
        registry=source_registry,
    )
    output_representations = [operation.matrix_cartesian for operation in operations]
    input_representations: list[tuple[tuple[float, ...], ...]]
    if specification["input_space"] == "symmetric":
        input_representations = [
            quadratic_field_representation(operation.matrix_cartesian).matrix_symmetric
            for operation in operations
        ]
    else:
        input_representations = [
            quadratic_field_representation(operation.matrix_cartesian).matrix_antisymmetric
            for operation in operations
        ]
    basis = equivariant_map_basis(
        output_representations,
        input_representations,
        tolerance=tolerance,
    )
    return InvariantTensorBasis(
        point_group_number=group.number,
        point_group=group.hm_symbol,
        schoenflies_symbol=group.schoenflies_symbol,
        response=response_name,
        output_basis=tuple(specification["output_basis"]),
        input_basis=tuple(specification["input_basis"]),
        equation=str(specification["equation"]),
        basis=basis,
    )


def load_optical_response_catalog(path: str | Path | None = None) -> dict[str, Any]:
    """Load the packaged all-point-group optical invariant catalog."""

    if path is None:
        repository_copy = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "optical_response_invariants.json"
        )
        text = (
            repository_copy.read_text(encoding="utf-8")
            if repository_copy.is_file()
            else resources.files("group_theory_operations")
            .joinpath("data/optical_response_invariants.json")
            .read_text(encoding="utf-8")
        )
    else:
        text = Path(path).read_text(encoding="utf-8")
    data = json.loads(text)
    groups = data.get("point_groups")
    if (
        data.get("schema_version") != 1
        or not isinstance(groups, list)
        or len(groups) != 32
        or [item.get("number") for item in groups if isinstance(item, dict)] != list(range(1, 33))
    ):
        raise GroupDataError("optical-response catalog schema_version must be 1")
    return data
