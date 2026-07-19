"""Typed access and scientific validation for the operation catalog."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


Matrix3 = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]


class GroupDataError(ValueError):
    """Raised when catalog data or a requested selection is inconsistent."""


@dataclass(frozen=True, slots=True)
class OperationRecord:
    """Stable machine interface for one point operation."""

    family: str
    index: int
    name: str
    matrix_fractional: Matrix3
    matrix_cartesian: Matrix3
    xyz_fractional: str | None = None
    ita: str | None = None

    def matrix(self, coordinate: str = "cartesian") -> Matrix3:
        if coordinate == "cartesian":
            return self.matrix_cartesian
        if coordinate == "fractional":
            return self.matrix_fractional
        raise GroupDataError("coordinate must be 'fractional' or 'cartesian'")

    def to_dict(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "index": self.index,
            "name": self.name,
            "xyz_fractional": self.xyz_fractional,
            "ita": self.ita,
            "matrix_fractional": [list(row) for row in self.matrix_fractional],
            "matrix_cartesian": [list(row) for row in self.matrix_cartesian],
        }


def canonical_name(value: str) -> str:
    """Normalize supported historical and LaTeX-like operation labels."""

    if not isinstance(value, str):
        raise GroupDataError("operation name must be a string")
    name = (
        value.strip()
        .replace("−", "-")
        .replace("^", "")
        .replace("{", "")
        .replace("}", "")
        .replace(" ", "")
    )
    if name.lower() in {"e", "identity"}:
        return "1"
    if name.lower() in {"i", "inversion"}:
        return "-1"
    if name in {"1", "-1"} or "_" in name:
        return name
    if name.startswith("m"):
        return "m_" + name[1:]
    for prefix in (
        "-6+", "-6-", "-4+", "-4-", "-3+", "-3-",
        "6+", "6-", "4+", "4-", "3+", "3-", "2",
    ):
        if name.startswith(prefix) and len(name) > len(prefix):
            return prefix + "_" + name[len(prefix):]
    return name


def _default_text() -> str:
    repository_copy = Path(__file__).resolve().parents[2] / "data" / "group_operations.json"
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/group_operations.json")
        .read_text(encoding="utf-8")
    )


def _default_schema_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2] / "schema" / "group-operations-v1.schema.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("schema/group-operations-v1.schema.json")
        .read_text(encoding="utf-8")
    )


def load_schema(path: str | Path | None = None) -> dict[str, Any]:
    """Load the public JSON Schema for catalog schema version 1."""

    try:
        text = _default_schema_text() if path is None else Path(path).read_text(encoding="utf-8")
        schema = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read schema file: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid schema JSON: {exc}") from exc
    if schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema":
        raise GroupDataError("unsupported catalog schema declaration")
    return schema


def load_database(path: str | Path | None = None, *, validate: bool = True) -> dict[str, Any]:
    """Load schema v1 from a path or the packaged canonical catalog."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read data file: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid JSON: {exc}") from exc
    if validate:
        errors = validate_database(data)
        if errors:
            raise GroupDataError("catalog validation failed:\n- " + "\n- ".join(errors))
    return data


def family_data(database: Mapping[str, Any], family: str) -> dict[str, Any]:
    try:
        return database["families"][family]
    except (KeyError, TypeError) as exc:
        choices = ", ".join(database.get("families", {}))
        raise GroupDataError(f"unknown family {family!r}; choices: {choices}") from exc


def iter_operations(
    database: Mapping[str, Any], family: str | None = None
) -> Iterable[tuple[str, dict[str, Any]]]:
    families = [family] if family else database["families"]
    for family_name in families:
        for operation in family_data(database, family_name)["operations"]:
            yield family_name, operation


def find_operations(
    database: Mapping[str, Any], name: str, family: str | None = None
) -> list[tuple[str, dict[str, Any]]]:
    target = canonical_name(name)
    matches: list[tuple[str, dict[str, Any]]] = []
    for family_name, operation in iter_operations(database, family):
        aliases = {operation["name"]}
        if "source_name" in operation:
            aliases.add(canonical_name(operation["source_name"]))
        if target in aliases:
            matches.append((family_name, operation))
    return matches


def get_operation(database: Mapping[str, Any], name: str, family: str) -> dict[str, Any]:
    matches = find_operations(database, name, family)
    if not matches:
        raise GroupDataError(f"operation {name!r} not found in {family}")
    if len(matches) != 1:
        raise GroupDataError(f"operation {name!r} is not unique in {family}")
    return matches[0][1]


def operation_record(
    database: Mapping[str, Any], name: str, family: str
) -> OperationRecord:
    operation = get_operation(database, name, family)
    return OperationRecord(
        family=family,
        index=int(operation["index"]),
        name=str(operation["name"]),
        matrix_fractional=_matrix3(operation["matrix_fractional"]),
        matrix_cartesian=_matrix3(operation["matrix_cartesian"]),
        xyz_fractional=operation.get("xyz_fractional", operation.get("xyz")),
        ita=operation.get("ita"),
    )


def get_layer_group(database: Mapping[str, Any], family: str, lg_number: int) -> dict[str, Any]:
    for entry in family_data(database, family)["layer_groups"]:
        if entry["LG"] == lg_number:
            return entry
    raise GroupDataError(f"LG{lg_number} not found in {family}")


def get_point_group(database: Mapping[str, Any], family: str, name: str) -> dict[str, Any]:
    target = name.lower()
    for entry in family_data(database, family)["point_groups"]:
        if entry["name"].lower() == target:
            return entry
    raise GroupDataError(f"point group {name!r} not found in {family}")


def multiply_operations(
    database: Mapping[str, Any], family: str, left: str, right: str
) -> str:
    """Return ``left * right``; with column vectors, ``right`` acts first."""

    multiplication = family_data(database, family).get("multiplication")
    if multiplication is None:
        raise GroupDataError(f"{family} has no multiplication table")
    left_name = get_operation(database, left, family)["name"]
    right_name = get_operation(database, right, family)["name"]
    try:
        return str(multiplication["table"][left_name][right_name])
    except KeyError as exc:
        raise GroupDataError(f"table lacks {left_name} * {right_name}") from exc


def matrix_for(operation: Mapping[str, Any], coordinate: str) -> list[list[float]]:
    key = "matrix_fractional" if coordinate == "fractional" else "matrix_cartesian"
    if coordinate not in {"fractional", "cartesian"}:
        raise GroupDataError("coordinate must be 'fractional' or 'cartesian'")
    try:
        return [list(row) for row in _matrix3(operation[key])]
    except (KeyError, TypeError, ValueError) as exc:
        raise GroupDataError(f"operation {operation.get('name')} lacks {coordinate} matrix") from exc


def _matrix3(values: Sequence[Sequence[Any]]) -> Matrix3:
    if (
        not isinstance(values, Sequence)
        or isinstance(values, (str, bytes))
        or len(values) != 3
        or any(
            not isinstance(row, Sequence)
            or isinstance(row, (str, bytes))
            or len(row) != 3
            for row in values
        )
    ):
        raise ValueError("matrix must have shape (3, 3)")
    if any(type(item) not in (int, float) for row in values for item in row):
        raise ValueError("matrix must contain JSON numbers, not strings or booleans")
    try:
        matrix = tuple(tuple(float(item) for item in row) for row in values)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("matrix must contain finite numbers") from exc
    if not all(math.isfinite(item) for row in matrix for item in row):
        raise ValueError("matrix must contain finite numbers")
    return matrix  # type: ignore[return-value]


def _matmul(left: Matrix3, right: Matrix3) -> Matrix3:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def _matrix_key(matrix: Matrix3, digits: int = 9) -> tuple[float, ...]:
    return tuple(round(item, digits) for row in matrix for item in row)


_IDENTITY: Matrix3 = (
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
)


def _transpose(matrix: Matrix3) -> Matrix3:
    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _basis_matrices(crystal_system: str) -> tuple[Matrix3, Matrix3] | None:
    if crystal_system in {"cubic", "tetragonal"}:
        return _IDENTITY, _IDENTITY
    if crystal_system == "hexagonal":
        root3 = math.sqrt(3.0)
        basis: Matrix3 = (
            (1.0, -0.5, 0.0),
            (0.0, root3 / 2.0, 0.0),
            (0.0, 0.0, 1.0),
        )
        inverse: Matrix3 = (
            (1.0, root3 / 3.0, 0.0),
            (0.0, 2.0 * root3 / 3.0, 0.0),
            (0.0, 0.0, 1.0),
        )
        return basis, inverse
    return None


def _validate_subgroup(
    context: str,
    members: Sequence[Any],
    matrices: Mapping[str, Matrix3],
    errors: list[str],
) -> None:
    if (
        isinstance(members, (str, bytes))
        or not isinstance(members, Sequence)
        or any(not isinstance(member, str) for member in members)
    ):
        errors.append(f"{context}: operation names must be a list of strings")
        return
    names = list(members)
    if not names or len(set(names)) != len(names):
        errors.append(f"{context}: operation names must be non-empty and unique")
        return
    if any(name not in matrices for name in names):
        errors.append(f"{context}: references an unknown operation")
        return
    identity_names = [
        name for name, matrix in matrices.items() if _matrix_key(matrix) == _matrix_key(_IDENTITY)
    ]
    if len(identity_names) != 1 or identity_names[0] not in names:
        errors.append(f"{context}: subgroup must contain the unique identity")
        return
    identity_name = identity_names[0]
    matrix_to_name = {_matrix_key(matrix): name for name, matrix in matrices.items()}
    member_set = set(names)
    for left in names:
        has_inverse = False
        for right in names:
            product_name = matrix_to_name.get(
                _matrix_key(_matmul(matrices[left], matrices[right]))
            )
            if product_name not in member_set:
                errors.append(f"{context}: not closed under {left} * {right}")
                return
            reverse_name = matrix_to_name.get(
                _matrix_key(_matmul(matrices[right], matrices[left]))
            )
            if product_name == identity_name and reverse_name == identity_name:
                has_inverse = True
        if not has_inverse:
            errors.append(f"{context}: {left} has no two-sided inverse")
            return


def validate_database(database: Any) -> tuple[str, ...]:
    """Return all detected schema/group-consistency errors without mutating data."""

    errors: list[str] = []
    if (
        not isinstance(database, dict)
        or type(database.get("schema_version")) is not int
        or database.get("schema_version") != 1
    ):
        return ("schema_version must be 1",)
    families = database.get("families")
    if not isinstance(families, dict) or not families:
        return ("families must be a non-empty object",)
    for family_name, family in families.items():
        if not isinstance(family_name, str) or not isinstance(family, dict):
            errors.append(f"{family_name}: family must be an object")
            continue
        operations = family.get("operations")
        if not isinstance(operations, list) or not operations:
            errors.append(f"{family_name}: operations must be non-empty")
            continue
        names: set[str] = set()
        indices: set[int] = set()
        matrices: dict[str, Matrix3] = {}
        matrix_owners: dict[tuple[float, ...], str] = {}
        crystal_system = family.get("crystal_system")
        basis_pair = (
            _basis_matrices(crystal_system) if isinstance(crystal_system, str) else None
        )
        if basis_pair is None:
            errors.append(f"{family_name}: unsupported or missing crystal_system")
        parent_point_group = family.get("parent_point_group")
        if not isinstance(parent_point_group, str) or not parent_point_group:
            errors.append(f"{family_name}: parent_point_group must be a non-empty string")
        for position, operation in enumerate(operations):
            context = f"{family_name}.operations[{position}]"
            if not isinstance(operation, dict):
                errors.append(f"{context}: operation must be an object")
                continue
            try:
                name = operation["name"]
                index = operation["index"]
                if not isinstance(name, str) or not name:
                    raise ValueError("name must be non-empty")
                if not isinstance(index, int) or isinstance(index, bool):
                    raise ValueError("index must be integer")
                fractional = _matrix3(operation["matrix_fractional"])
                cartesian = _matrix3(operation["matrix_cartesian"])
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"{context}: {exc}")
                continue
            if name in names:
                errors.append(f"{family_name}: duplicate operation name {name}")
            if index in indices:
                errors.append(f"{family_name}: duplicate operation index {index}")
            names.add(name)
            indices.add(index)
            matrices[name] = fractional
            for field_name in ("xyz", "xyz_fractional", "ita", "source_name", "parity_block"):
                if field_name in operation and not isinstance(operation[field_name], str):
                    errors.append(f"{context}: {field_name} must be a string")
            matrix_key = _matrix_key(fractional)
            if matrix_key in matrix_owners:
                errors.append(
                    f"{family_name}: {name} duplicates matrix of {matrix_owners[matrix_key]}"
                )
            matrix_owners[matrix_key] = name
            product = _matmul(_transpose(cartesian), cartesian)
            if _matrix_key(product) != _matrix_key(_IDENTITY):
                errors.append(f"{family_name}.{name}: Cartesian matrix is not orthogonal")
            if basis_pair is not None:
                basis, basis_inverse = basis_pair
                expected_cartesian = _matmul(_matmul(basis, fractional), basis_inverse)
                if _matrix_key(expected_cartesian) != _matrix_key(cartesian):
                    errors.append(
                        f"{family_name}.{name}: fractional/Cartesian basis mapping disagrees"
                    )

        by_index = {
            op.get("index"): op.get("name")
            for op in operations
            if isinstance(op, dict) and type(op.get("index")) is int
        }
        for group_kind in ("point_groups", "layer_groups"):
            entries = family.get(group_kind, [])
            if not isinstance(entries, list):
                errors.append(f"{family_name}.{group_kind}: must be a list")
                continue
            point_group_names: set[str] = set()
            layer_group_numbers: set[int] = set()
            for entry_index, entry in enumerate(entries):
                context = f"{family_name}.{group_kind}[{entry_index}]"
                if not isinstance(entry, dict):
                    errors.append(f"{context}: entry must be an object")
                    continue
                if group_kind == "point_groups":
                    group_name = entry.get("name")
                    group_order = entry.get("order")
                    referenced_operations = entry.get("operations")
                    if not isinstance(group_name, str) or not group_name:
                        errors.append(f"{context}: name must be a non-empty string")
                    elif group_name in point_group_names:
                        errors.append(f"{context}: duplicate point-group name {group_name}")
                    else:
                        point_group_names.add(group_name)
                    if (
                        type(group_order) is not int
                        or not isinstance(referenced_operations, list)
                        or group_order != len(referenced_operations)
                    ):
                        errors.append(f"{context}: order must equal the operation count")
                    pairs = ((entry.get("operation_indices", []), entry.get("operations", [])),)
                else:
                    layer_number = entry.get("LG")
                    if type(layer_number) is not int or not 1 <= layer_number <= 80:
                        errors.append(f"{context}: LG must be an integer from 1 to 80")
                    elif layer_number in layer_group_numbers:
                        errors.append(f"{context}: duplicate LG number {layer_number}")
                    else:
                        layer_group_numbers.add(layer_number)
                    for field_name in ("point_group", "point_group_base"):
                        value = entry.get(field_name)
                        if not isinstance(value, str) or not value:
                            errors.append(f"{context}: {field_name} must be a non-empty string")
                    embedding = entry.get("point_group_embedding")
                    if embedding is not None and not isinstance(embedding, str):
                        errors.append(
                            f"{context}: point_group_embedding must be a string or null"
                        )
                    pairs = (
                        (entry.get("R+_indices", []), entry.get("R+", [])),
                        (entry.get("R-_indices", []), entry.get("R-", [])),
                    )
                for referenced_indices, referenced_names in pairs:
                    if not isinstance(referenced_indices, list) or not isinstance(
                        referenced_names, list
                    ):
                        errors.append(f"{context}: index/name references must be lists")
                        continue
                    if any(type(index) is not int for index in referenced_indices):
                        errors.append(f"{context}: operation indices must be integers")
                        continue
                    if any(not isinstance(name, str) for name in referenced_names):
                        errors.append(f"{context}: operation names must be strings")
                        continue
                    resolved = [by_index.get(index) for index in referenced_indices]
                    if resolved != referenced_names or any(
                        name not in names for name in referenced_names
                    ):
                        errors.append(f"{context}: index/name references disagree")
                members = (
                    entry.get("operations", [])
                    if group_kind == "point_groups"
                    else list(entry.get("R+", [])) + list(entry.get("R-", []))
                    if isinstance(entry.get("R+", []), list)
                    and isinstance(entry.get("R-", []), list)
                    else ()
                )
                _validate_subgroup(context, members, matrices, errors)

        multiplication = family.get("multiplication")
        if multiplication is not None:
            if not isinstance(multiplication, dict):
                errors.append(f"{family_name}.multiplication: must be an object")
                continue
            order = multiplication.get("element_order", [])
            table = multiplication.get("table", {})
            if (
                not isinstance(order, list)
                or any(not isinstance(item, str) for item in order)
                or not isinstance(table, dict)
                or set(order) != names
                or list(table) != order
            ):
                errors.append(f"{family_name}.multiplication: element order disagrees with operations")
            else:
                identity = multiplication.get("identity")
                identity_names = [
                    name
                    for name, matrix in matrices.items()
                    if _matrix_key(matrix) == _matrix_key(_IDENTITY)
                ]
                if len(identity_names) != 1 or identity != identity_names[0]:
                    errors.append(f"{family_name}.multiplication: wrong identity")
                inverse = multiplication.get("inverse")
                if not isinstance(inverse, dict) or set(inverse) != set(order):
                    errors.append(f"{family_name}.multiplication: inverse map is incomplete")
                for left in order:
                    row = table.get(left)
                    if not isinstance(row, dict) or list(row) != order:
                        errors.append(f"{family_name}.multiplication.{left}: incomplete row")
                        continue
                    for right in order:
                        result = row[right]
                        if (
                            not isinstance(result, str)
                            or result not in matrices
                            or _matrix_key(_matmul(matrices[left], matrices[right]))
                            != _matrix_key(matrices[result])
                        ):
                            errors.append(f"{family_name}.multiplication: wrong {left} * {right}")
                    row_results = list(row.values())
                    if any(not isinstance(result, str) for result in row_results) or set(
                        result for result in row_results if isinstance(result, str)
                    ) != set(order):
                        errors.append(f"{family_name}.multiplication.{left}: row is not Latin")
                for right in order:
                    column = [table[left].get(right) for left in order if isinstance(table.get(left), dict)]
                    if (
                        any(not isinstance(result, str) for result in column)
                        or set(result for result in column if isinstance(result, str))
                        != set(order)
                        or len(column) != len(order)
                    ):
                        errors.append(f"{family_name}.multiplication.{right}: column is not Latin")
                if isinstance(inverse, dict) and isinstance(identity, str):
                    for name in order:
                        inverse_name = inverse.get(name)
                        if inverse_name not in order:
                            errors.append(f"{family_name}.multiplication: invalid inverse for {name}")
                            continue
                        left_row = table.get(name)
                        inverse_row = table.get(inverse_name)
                        if not isinstance(left_row, dict) or not isinstance(
                            inverse_row, dict
                        ):
                            errors.append(
                                f"{family_name}.multiplication: invalid inverse rows for {name}"
                            )
                            continue
                        if (
                            left_row.get(inverse_name) != identity
                            or inverse_row.get(name) != identity
                        ):
                            errors.append(f"{family_name}.multiplication: wrong inverse for {name}")
    return tuple(dict.fromkeys(errors))
