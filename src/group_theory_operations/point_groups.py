"""Registry and typed access for the 32 crystallographic point groups."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .catalog import GroupDataError, OperationRecord, load_database, operation_record


@dataclass(frozen=True, slots=True)
class CrystallographicPointGroup:
    """One standard crystallographic point-group embedding."""

    number: int
    hm_symbol: str
    schoenflies_symbol: str
    aliases: tuple[str, ...]
    crystal_system: str
    host_family: str
    order: int
    generators: tuple[str, ...]
    operations: tuple[str, ...]
    centrosymmetric: bool
    polar: bool
    chiral: bool

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "CrystallographicPointGroup":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid crystallographic point-group record")
        integer_fields = ("number", "order")
        string_fields = (
            "hm_symbol",
            "schoenflies_symbol",
            "crystal_system",
            "host_family",
        )
        sequence_fields = ("aliases", "generators", "operations")
        boolean_fields = ("centrosymmetric", "polar", "chiral")
        if any(type(value.get(field)) is not int for field in integer_fields):
            raise GroupDataError("point-group number and order must be integers")
        if any(not isinstance(value.get(field), str) or not value[field] for field in string_fields):
            raise GroupDataError("point-group symbols, crystal system and host family must be strings")
        if any(
            not isinstance(value.get(field), list)
            or any(not isinstance(item, str) or not item for item in value[field])
            for field in sequence_fields
        ):
            raise GroupDataError("point-group aliases, generators and operations must be string lists")
        if any(type(value.get(field)) is not bool for field in boolean_fields):
            raise GroupDataError("point-group classification flags must be booleans")
        try:
            return cls(
                number=value["number"],
                hm_symbol=value["hm_symbol"],
                schoenflies_symbol=value["schoenflies_symbol"],
                aliases=tuple(value["aliases"]),
                crystal_system=value["crystal_system"],
                host_family=value["host_family"],
                order=value["order"],
                generators=tuple(value["generators"]),
                operations=tuple(value["operations"]),
                centrosymmetric=value["centrosymmetric"],
                polar=value["polar"],
                chiral=value["chiral"],
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GroupDataError("invalid crystallographic point-group record") from exc

    def to_dict(self) -> dict[str, Any]:
        return {
            "number": self.number,
            "hm_symbol": self.hm_symbol,
            "schoenflies_symbol": self.schoenflies_symbol,
            "aliases": list(self.aliases),
            "crystal_system": self.crystal_system,
            "host_family": self.host_family,
            "order": self.order,
            "generators": list(self.generators),
            "operations": list(self.operations),
            "centrosymmetric": self.centrosymmetric,
            "polar": self.polar,
            "chiral": self.chiral,
        }


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "crystallographic_point_groups.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/crystallographic_point_groups.json")
        .read_text(encoding="utf-8")
    )


def load_point_group_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load the versioned 32-point-group registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read point-group registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid point-group registry JSON: {exc}") from exc
    groups = data.get("point_groups")
    if data.get("schema_version") != 1 or not isinstance(groups, list):
        raise GroupDataError("point-group registry must use schema_version 1")
    records = [CrystallographicPointGroup.from_mapping(item) for item in groups]
    if [record.number for record in records] != list(range(1, 33)):
        raise GroupDataError("point-group registry must contain the ordered numbers 1-32")
    if len({record.hm_symbol for record in records}) != 32:
        raise GroupDataError("point-group Hermann-Mauguin symbols must be unique")
    if any(record.order != len(record.operations) for record in records):
        raise GroupDataError("point-group order must equal its operation count")
    if any(
        not record.operations
        or record.operations[0] != "1"
        or len(set(record.operations)) != record.order
        or not set(record.generators).issubset(record.operations)
        or record.centrosymmetric != ("-1" in record.operations)
        for record in records
    ):
        raise GroupDataError("point-group operations or classification flags are inconsistent")
    return data


def iter_crystallographic_point_groups(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[CrystallographicPointGroup]:
    """Yield the 32 point groups in their standard crystallographic order."""

    source = load_point_group_registry() if registry is None else registry
    for item in source["point_groups"]:
        yield CrystallographicPointGroup.from_mapping(item)


def _lookup_key(value: str) -> str:
    return value.strip().replace(" ", "").replace("−", "-").lower()


def get_crystallographic_point_group(
    name_or_number: str | int,
    registry: Mapping[str, Any] | None = None,
) -> CrystallographicPointGroup:
    """Resolve a standard number, Hermann-Mauguin symbol, or Schoenflies symbol."""

    groups = tuple(iter_crystallographic_point_groups(registry))
    if isinstance(name_or_number, int):
        matches = [group for group in groups if group.number == name_or_number]
    elif isinstance(name_or_number, str):
        key = _lookup_key(name_or_number)
        matches = [
            group
            for group in groups
            if key
            in {
                str(group.number),
                _lookup_key(group.hm_symbol),
                _lookup_key(group.schoenflies_symbol),
                *(_lookup_key(alias) for alias in group.aliases),
            }
        ]
    else:
        raise GroupDataError("point group must be a number or symbol")
    if len(matches) != 1:
        raise GroupDataError(f"unknown crystallographic point group {name_or_number!r}")
    return matches[0]


def point_group_operations(
    name_or_number: str | int,
    *,
    database: Mapping[str, Any] | None = None,
    registry: Mapping[str, Any] | None = None,
) -> tuple[OperationRecord, ...]:
    """Resolve all operations of one standard point-group embedding."""

    group = get_crystallographic_point_group(name_or_number, registry)
    source = load_database() if database is None else database
    return tuple(
        operation_record(source, name, group.host_family) for name in group.operations
    )
