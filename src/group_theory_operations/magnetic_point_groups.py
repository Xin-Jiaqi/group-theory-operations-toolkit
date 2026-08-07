"""Registry and typed access for the 122 magnetic point groups."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .catalog import GroupDataError, OperationRecord, load_database, operation_record
from .point_groups import iter_crystallographic_point_groups, load_point_group_registry


MAGNETIC_CATEGORIES = (
    "type_I",
    "type_II_gray",
    "type_III_black_white",
)


def _matrix_key(matrix: tuple[tuple[float, ...], ...]) -> tuple[float, ...]:
    return tuple(round(float(value), 9) for row in matrix for value in row)


def _matmul(
    left: tuple[tuple[float, ...], ...],
    right: tuple[tuple[float, ...], ...],
) -> tuple[tuple[float, ...], ...]:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


@dataclass(frozen=True, slots=True)
class MagneticPointOperation:
    """A spatial point operation with an optional time-reversal factor."""

    name: str
    time_reversal: bool

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MagneticPointOperation":
        if (
            not isinstance(value, Mapping)
            or not isinstance(value.get("name"), str)
            or not value["name"]
            or type(value.get("time_reversal")) is not bool
        ):
            raise GroupDataError("invalid magnetic point operation")
        return cls(name=value["name"], time_reversal=value["time_reversal"])

    @property
    def label(self) -> str:
        return f"{self.name}'" if self.time_reversal else self.name

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "time_reversal": self.time_reversal}


@dataclass(frozen=True, slots=True)
class MagneticPointGroup:
    """One standard magnetic point-group embedding."""

    number: int
    magnetic_number: str
    hm_symbol: str
    category: str
    parent_point_group_number: int
    parent_point_group_hm: str
    crystal_system: str
    host_family: str
    order: int
    unitary_subgroup_order: int
    generators: tuple[MagneticPointOperation, ...]
    operations: tuple[MagneticPointOperation, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MagneticPointGroup":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid magnetic point-group record")
        integer_fields = (
            "number",
            "parent_point_group_number",
            "order",
            "unitary_subgroup_order",
        )
        string_fields = (
            "magnetic_number",
            "hm_symbol",
            "category",
            "parent_point_group_hm",
            "crystal_system",
            "host_family",
        )
        if any(type(value.get(field)) is not int for field in integer_fields):
            raise GroupDataError("magnetic point-group counts and numbers must be integers")
        if any(not isinstance(value.get(field), str) or not value[field] for field in string_fields):
            raise GroupDataError("magnetic point-group symbols and metadata must be strings")
        if value["category"] not in MAGNETIC_CATEGORIES:
            raise GroupDataError("unknown magnetic point-group category")
        try:
            generators = tuple(
                MagneticPointOperation.from_mapping(item) for item in value["generators"]
            )
            operations = tuple(
                MagneticPointOperation.from_mapping(item) for item in value["operations"]
            )
        except (KeyError, TypeError) as exc:
            raise GroupDataError("magnetic point-group operations must be lists") from exc
        return cls(
            number=value["number"],
            magnetic_number=value["magnetic_number"],
            hm_symbol=value["hm_symbol"],
            category=value["category"],
            parent_point_group_number=value["parent_point_group_number"],
            parent_point_group_hm=value["parent_point_group_hm"],
            crystal_system=value["crystal_system"],
            host_family=value["host_family"],
            order=value["order"],
            unitary_subgroup_order=value["unitary_subgroup_order"],
            generators=generators,
            operations=operations,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "number": self.number,
            "magnetic_number": self.magnetic_number,
            "hm_symbol": self.hm_symbol,
            "category": self.category,
            "parent_point_group_number": self.parent_point_group_number,
            "parent_point_group_hm": self.parent_point_group_hm,
            "crystal_system": self.crystal_system,
            "host_family": self.host_family,
            "order": self.order,
            "unitary_subgroup_order": self.unitary_subgroup_order,
            "generators": [item.to_dict() for item in self.generators],
            "operations": [item.to_dict() for item in self.operations],
        }


@dataclass(frozen=True, slots=True)
class ResolvedMagneticPointOperation:
    """A magnetic operation paired with its resolved spatial matrix record."""

    spatial: OperationRecord
    time_reversal: bool

    @property
    def name(self) -> str:
        return self.spatial.name

    @property
    def label(self) -> str:
        return f"{self.name}'" if self.time_reversal else self.name


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2] / "data" / "magnetic_point_groups.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/magnetic_point_groups.json")
        .read_text(encoding="utf-8")
    )


def load_magnetic_point_group_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the versioned 122-magnetic-point-group registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read magnetic point-group registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid magnetic point-group registry JSON: {exc}") from exc
    values = data.get("magnetic_point_groups")
    if data.get("schema_version") != 1 or not isinstance(values, list):
        raise GroupDataError("magnetic point-group registry must use schema_version 1")
    records = [MagneticPointGroup.from_mapping(item) for item in values]
    if [record.number for record in records] != list(range(1, 123)):
        raise GroupDataError("magnetic point-group registry must contain ordered numbers 1-122")
    if len({record.magnetic_number for record in records}) != 122:
        raise GroupDataError("magnetic point-group numbers must be unique")
    if len({record.hm_symbol for record in records}) != 122:
        raise GroupDataError("magnetic point-group symbols must be unique")
    if Counter(record.category for record in records) != {
        "type_I": 32,
        "type_II_gray": 32,
        "type_III_black_white": 58,
    }:
        raise GroupDataError("magnetic point-group category counts must be 32 + 32 + 58")
    parents = {
        item.number: item
        for item in iter_crystallographic_point_groups(load_point_group_registry())
    }
    database = load_database()
    local_numbers: Counter[int] = Counter()
    for record in records:
        parent = parents.get(record.parent_point_group_number)
        pairs = {(item.name, item.time_reversal) for item in record.operations}
        try:
            parent_number, local_number, global_number = (
                int(value) for value in record.magnetic_number.split(".")
            )
        except (TypeError, ValueError) as exc:
            raise GroupDataError("invalid three-part magnetic point-group number") from exc
        local_numbers[parent_number] += 1
        if (
            parent is None
            or (parent_number, local_number, global_number)
            != (
                record.parent_point_group_number,
                local_numbers[parent_number],
                record.number,
            )
            or parent.hm_symbol != record.parent_point_group_hm
            or record.order != len(record.operations)
            or len(pairs) != record.order
            or ("1", False) not in pairs
            or not {(item.name, item.time_reversal) for item in record.generators}.issubset(pairs)
            or record.unitary_subgroup_order
            != sum(not item.time_reversal for item in record.operations)
            or any(item.name not in parent.operations for item in record.operations)
        ):
            raise GroupDataError("magnetic point-group operations are inconsistent")
        spatial_counts = Counter(item.name for item in record.operations)
        expected_multiplicity = 2 if record.category == "type_II_gray" else 1
        if spatial_counts != Counter(
            {name: expected_multiplicity for name in parent.operations}
        ):
            raise GroupDataError("magnetic point group must cover every parent operation")
        pure_time_reversal = ("1", True) in pairs
        if record.category == "type_I" and (
            record.order != parent.order or any(item.time_reversal for item in record.operations)
        ):
            raise GroupDataError("type-I magnetic point group is inconsistent")
        if record.category == "type_II_gray" and (
            record.order != 2 * parent.order
            or not pure_time_reversal
            or record.unitary_subgroup_order != parent.order
        ):
            raise GroupDataError("type-II gray magnetic point group is inconsistent")
        if record.category == "type_III_black_white" and (
            record.order != parent.order
            or pure_time_reversal
            or record.unitary_subgroup_order * 2 != parent.order
        ):
            raise GroupDataError("type-III magnetic point group is inconsistent")
        spatial_records = {
            name: operation_record(database, name, parent.host_family)
            for name in parent.operations
        }
        by_matrix = {
            _matrix_key(item.matrix_cartesian): item.name
            for item in spatial_records.values()
        }
        for left in record.operations:
            for right in record.operations:
                product_name = by_matrix[
                    _matrix_key(
                        _matmul(
                            spatial_records[left.name].matrix_cartesian,
                            spatial_records[right.name].matrix_cartesian,
                        )
                    )
                ]
                if (product_name, left.time_reversal ^ right.time_reversal) not in pairs:
                    raise GroupDataError("magnetic point group is not closed")
    return data


def iter_magnetic_point_groups(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[MagneticPointGroup]:
    """Yield all 122 magnetic point groups in standard order."""

    source = load_magnetic_point_group_registry() if registry is None else registry
    for item in source["magnetic_point_groups"]:
        yield MagneticPointGroup.from_mapping(item)


def _lookup_key(value: str) -> str:
    return (
        value.strip()
        .replace(" ", "")
        .replace("−", "-")
        .replace("′", "'")
        .lower()
    )


def get_magnetic_point_group(
    identifier: str | int,
    registry: Mapping[str, Any] | None = None,
) -> MagneticPointGroup:
    """Resolve a global number, three-part number, or magnetic HM symbol."""

    groups = tuple(iter_magnetic_point_groups(registry))
    if isinstance(identifier, int):
        matches = [group for group in groups if group.number == identifier]
    elif isinstance(identifier, str):
        key = _lookup_key(identifier)
        matches = [group for group in groups if key == _lookup_key(group.hm_symbol)]
        if not matches:
            matches = [
                group for group in groups if key == _lookup_key(group.magnetic_number)
            ]
        if not matches and key.isdecimal():
            matches = [group for group in groups if int(key) == group.number]
    else:
        raise GroupDataError("magnetic point group must be a number or symbol")
    if len(matches) != 1:
        raise GroupDataError(f"unknown magnetic point group {identifier!r}")
    return matches[0]


def magnetic_point_group_operations(
    identifier: str | int,
    *,
    database: Mapping[str, Any] | None = None,
    registry: Mapping[str, Any] | None = None,
) -> tuple[ResolvedMagneticPointOperation, ...]:
    """Resolve spatial matrices and time-reversal labels for one group."""

    group = get_magnetic_point_group(identifier, registry)
    source = load_database() if database is None else database
    return tuple(
        ResolvedMagneticPointOperation(
            spatial=operation_record(source, item.name, group.host_family),
            time_reversal=item.time_reversal,
        )
        for item in group.operations
    )
