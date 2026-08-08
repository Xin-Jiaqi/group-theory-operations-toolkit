"""Typed access to all 528 magnetic layer-group point co-groups."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, cast

from .catalog import GroupDataError
from .layer_groups import (
    get_crystallographic_layer_group,
    load_layer_group_registry,
)


MAGNETIC_LAYER_TYPES = ("I", "II", "III", "IV")
_EXPECTED_TYPES = {"I": 80, "II": 80, "III": 246, "IV": 122}


def _matrix(value: Any, field: str) -> tuple[tuple[float, ...], ...]:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(not isinstance(row, list) or len(row) != 3 for row in value)
        or any(
            type(item) not in {int, float} for row in value for item in row
        )
    ):
        raise GroupDataError(f"{field} must be a numeric 3x3 matrix")
    return tuple(tuple(float(item) for item in row) for row in value)


@dataclass(frozen=True, slots=True)
class MagneticLayerPointOperation:
    name: str
    label: str
    time_reversal: bool
    matrix_fractional: tuple[tuple[float, ...], ...]
    matrix_cartesian: tuple[tuple[float, ...], ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MagneticLayerPointOperation":
        if (
            not isinstance(value, Mapping)
            or not isinstance(value.get("name"), str)
            or not value["name"]
            or not isinstance(value.get("label"), str)
            or type(value.get("time_reversal")) is not bool
        ):
            raise GroupDataError("invalid magnetic layer point operation")
        expected_label = (
            f"{value['name']}'" if value["time_reversal"] else value["name"]
        )
        if value["label"] != expected_label:
            raise GroupDataError("magnetic operation label and time reversal disagree")
        return cls(
            name=value["name"],
            label=value["label"],
            time_reversal=value["time_reversal"],
            matrix_fractional=_matrix(value.get("matrix_fractional"), "matrix_fractional"),
            matrix_cartesian=_matrix(value.get("matrix_cartesian"), "matrix_cartesian"),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "label": self.label,
            "time_reversal": self.time_reversal,
            "matrix_fractional": [list(row) for row in self.matrix_fractional],
            "matrix_cartesian": [list(row) for row in self.matrix_cartesian],
        }


@dataclass(frozen=True, slots=True)
class CorrespondingMagneticSpaceGroup:
    bns_number: str
    uni_number: int
    og_number: str
    basis_transform_key: str
    basis_transform: tuple[tuple[float, ...], ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "CorrespondingMagneticSpaceGroup":
        if (
            not isinstance(value, Mapping)
            or not isinstance(value.get("bns_number"), str)
            or type(value.get("uni_number")) is not int
            or not 1 <= value["uni_number"] <= 1651
            or not isinstance(value.get("og_number"), str)
            or value.get("basis_transform_key") not in {"I", "P1", "P2", "P3", "P4", "P5"}
        ):
            raise GroupDataError("invalid corresponding magnetic space group")
        return cls(
            bns_number=value["bns_number"],
            uni_number=value["uni_number"],
            og_number=value["og_number"],
            basis_transform_key=value["basis_transform_key"],
            basis_transform=_matrix(value.get("basis_transform"), "basis_transform"),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "bns_number": self.bns_number,
            "uni_number": self.uni_number,
            "og_number": self.og_number,
            "basis_transform_key": self.basis_transform_key,
            "basis_transform": [list(row) for row in self.basis_transform],
        }


@dataclass(frozen=True, slots=True)
class MagneticLayerGroup:
    global_number: int
    og_number: str
    family_number: int
    magnetic_type: str
    litvin_og_symbol_ascii: str
    magnetic_point_group_symbol: str
    parent_layer_group_number: int
    parent_layer_group_symbol: str
    crystal_system: str
    host_family: str
    point_operation_count: int
    unitary_subgroup_order: int
    anti_translation_fractional: tuple[float, float, float] | None
    corresponding_magnetic_space_group: CorrespondingMagneticSpaceGroup
    source_pdf_page: int
    point_operations: tuple[MagneticLayerPointOperation, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "MagneticLayerGroup":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid magnetic layer-group record")
        integer_fields = (
            "global_number",
            "family_number",
            "parent_layer_group_number",
            "point_operation_count",
            "unitary_subgroup_order",
            "source_pdf_page",
        )
        string_fields = (
            "og_number",
            "magnetic_type",
            "litvin_og_symbol_ascii",
            "magnetic_point_group_symbol",
            "parent_layer_group_symbol",
            "crystal_system",
            "host_family",
        )
        if any(type(value.get(field)) is not int for field in integer_fields):
            raise GroupDataError("magnetic layer-group numbers and counts must be integers")
        if any(
            not isinstance(value.get(field), str) or not value[field]
            for field in string_fields
        ):
            raise GroupDataError("magnetic layer-group symbols and metadata must be strings")
        if value["magnetic_type"] not in MAGNETIC_LAYER_TYPES:
            raise GroupDataError("unknown magnetic layer-group type")
        anti = value.get("anti_translation_fractional")
        anti_translation: tuple[float, float, float] | None
        if anti is None:
            anti_translation = None
        elif (
            isinstance(anti, list)
            and len(anti) == 3
            and all(type(item) in {int, float} for item in anti)
        ):
            anti_translation = cast(
                tuple[float, float, float], tuple(float(item) for item in anti)
            )
        else:
            raise GroupDataError("anti-translation must be null or a numeric triplet")
        try:
            operations = tuple(
                MagneticLayerPointOperation.from_mapping(item)
                for item in value["point_operations"]
            )
            corresponding = CorrespondingMagneticSpaceGroup.from_mapping(
                value["corresponding_magnetic_space_group"]
            )
        except (KeyError, TypeError) as exc:
            raise GroupDataError("invalid magnetic layer-group operation list") from exc
        return cls(
            global_number=value["global_number"],
            og_number=value["og_number"],
            family_number=value["family_number"],
            magnetic_type=value["magnetic_type"],
            litvin_og_symbol_ascii=value["litvin_og_symbol_ascii"],
            magnetic_point_group_symbol=value["magnetic_point_group_symbol"],
            parent_layer_group_number=value["parent_layer_group_number"],
            parent_layer_group_symbol=value["parent_layer_group_symbol"],
            crystal_system=value["crystal_system"],
            host_family=value["host_family"],
            point_operation_count=value["point_operation_count"],
            unitary_subgroup_order=value["unitary_subgroup_order"],
            anti_translation_fractional=anti_translation,
            corresponding_magnetic_space_group=corresponding,
            source_pdf_page=value["source_pdf_page"],
            point_operations=operations,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "global_number": self.global_number,
            "og_number": self.og_number,
            "family_number": self.family_number,
            "magnetic_type": self.magnetic_type,
            "litvin_og_symbol_ascii": self.litvin_og_symbol_ascii,
            "magnetic_point_group_symbol": self.magnetic_point_group_symbol,
            "parent_layer_group_number": self.parent_layer_group_number,
            "parent_layer_group_symbol": self.parent_layer_group_symbol,
            "crystal_system": self.crystal_system,
            "host_family": self.host_family,
            "point_operation_count": self.point_operation_count,
            "unitary_subgroup_order": self.unitary_subgroup_order,
            "anti_translation_fractional": (
                list(self.anti_translation_fractional)
                if self.anti_translation_fractional is not None
                else None
            ),
            "corresponding_magnetic_space_group": (
                self.corresponding_magnetic_space_group.to_dict()
            ),
            "source_pdf_page": self.source_pdf_page,
            "point_operations": [item.to_dict() for item in self.point_operations],
        }


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2] / "data" / "magnetic_layer_groups.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/magnetic_layer_groups.json")
        .read_text(encoding="utf-8")
    )


def load_magnetic_layer_group_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the 528-group point-co-group registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read magnetic layer-group registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid magnetic layer-group registry JSON: {exc}") from exc
    values = data.get("magnetic_layer_groups")
    if data.get("schema_version") != 1 or not isinstance(values, list):
        raise GroupDataError("magnetic layer-group registry must use schema_version 1")
    records = [MagneticLayerGroup.from_mapping(item) for item in values]
    if [item.global_number for item in records] != list(range(1, 529)):
        raise GroupDataError("magnetic layer groups must contain ordered numbers 1-528")
    if Counter(item.magnetic_type for item in records) != _EXPECTED_TYPES:
        raise GroupDataError("magnetic layer-group type counts must be 80 + 80 + 246 + 122")
    if len({item.og_number for item in records}) != 528 or len(
        {item.litvin_og_symbol_ascii for item in records}
    ) != 528:
        raise GroupDataError("magnetic layer-group numbers and symbols must be unique")
    layer_registry = load_layer_group_registry()
    family_counts: Counter[int] = Counter()
    for record in records:
        parent_number, family_number, global_number = map(int, record.og_number.split("."))
        family_counts[parent_number] += 1
        parent = get_crystallographic_layer_group(parent_number, layer_registry)
        spatial_counts = Counter(item.name for item in record.point_operations)
        expected_multiplicity = 2 if record.magnetic_type in {"II", "IV"} else 1
        if (
            (parent_number, family_number, global_number)
            != (record.parent_layer_group_number, record.family_number, record.global_number)
            or family_number != family_counts[parent_number]
            or record.parent_layer_group_symbol != parent.international_short
            or record.crystal_system != parent.crystal_system
            or record.point_operation_count != len(record.point_operations)
            or record.unitary_subgroup_order
            != sum(not item.time_reversal for item in record.point_operations)
            or len(set(spatial_counts.values())) != 1
            or next(iter(spatial_counts.values())) != expected_multiplicity
            or (record.magnetic_type == "IV")
            != (record.anti_translation_fractional is not None)
        ):
            raise GroupDataError("magnetic layer-group record is inconsistent")
        pure_time_reversal = any(
            item.name == "1" and item.time_reversal for item in record.point_operations
        )
        if record.magnetic_type == "I" and any(
            item.time_reversal for item in record.point_operations
        ):
            raise GroupDataError("type-I magnetic layer group must be unitary")
        if record.magnetic_type == "III" and (
            pure_time_reversal
            or record.unitary_subgroup_order * 2 != record.point_operation_count
        ):
            raise GroupDataError("type-III magnetic layer group is inconsistent")
        if record.magnetic_type in {"II", "IV"} and not pure_time_reversal:
            raise GroupDataError("type-II/IV point co-groups must contain time reversal")
    return data


def iter_magnetic_layer_groups(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[MagneticLayerGroup]:
    data = load_magnetic_layer_group_registry() if registry is None else registry
    for item in data["magnetic_layer_groups"]:
        yield MagneticLayerGroup.from_mapping(item)


def _normalize_symbol(value: str) -> str:
    return "".join(value.lower().replace("_", "").split())


def get_magnetic_layer_group(
    identifier: str | int,
    registry: Mapping[str, Any] | None = None,
) -> MagneticLayerGroup:
    data = load_magnetic_layer_group_registry() if registry is None else registry
    values = data.get("magnetic_layer_groups")
    if not isinstance(values, list):
        raise GroupDataError("invalid magnetic layer-group registry")
    if type(identifier) is int:
        matches = [item for item in values if item.get("global_number") == identifier]
    elif isinstance(identifier, str):
        token = identifier.strip()
        if token.upper().startswith("MLG") and token[3:].isdigit():
            matches = [
                item for item in values if item.get("global_number") == int(token[3:])
            ]
        elif token.isdigit():
            matches = [
                item for item in values if item.get("global_number") == int(token)
            ]
        elif _is_full_number(token):
            matches = [item for item in values if item.get("og_number") == token]
        else:
            normalized = _normalize_symbol(token)
            matches = [
                item
                for item in values
                if isinstance(item.get("litvin_og_symbol_ascii"), str)
                and _normalize_symbol(item["litvin_og_symbol_ascii"]) == normalized
            ]
    else:
        raise GroupDataError("magnetic layer-group identifier must be a string or integer")
    if len(matches) != 1:
        raise GroupDataError(f"unknown or ambiguous magnetic layer group {identifier!r}")
    return MagneticLayerGroup.from_mapping(matches[0])


def _is_full_number(value: str) -> bool:
    parts = value.split(".")
    return len(parts) == 3 and all(part.isdigit() for part in parts)
