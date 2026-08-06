"""Registry and typed access for the 230 crystallographic space groups."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .catalog import GroupDataError


@dataclass(frozen=True, slots=True)
class HallSetting:
    """One Hall setting of a space group, with its Seitz generators."""

    hall_number: int
    hall_symbol: str
    choice: str
    centering: str
    symmorphic: bool
    operation_count: int
    generators: tuple[dict[str, Any], ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "HallSetting":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid Hall-setting record")
        if type(value.get("hall_number")) is not int or value["hall_number"] < 1:
            raise GroupDataError("hall_number must be a positive integer")
        if not isinstance(value.get("hall_symbol"), str) or not value["hall_symbol"]:
            raise GroupDataError("hall_symbol must be a non-empty string")
        if not isinstance(value.get("choice"), str):
            raise GroupDataError("choice must be a string")
        if value.get("centering") not in {"P", "A", "B", "C", "F", "I", "R"}:
            raise GroupDataError("centering must be one of P, A, B, C, F, I, R")
        if type(value.get("symmorphic")) is not bool:
            raise GroupDataError("symmorphic must be a boolean")
        if type(value.get("operation_count")) is not int or value["operation_count"] < 1:
            raise GroupDataError("operation_count must be a positive integer")
        generators = value.get("generators")
        if not isinstance(generators, list):
            raise GroupDataError("generators must be a list")
        if generators and any(
            not isinstance(generator, Mapping)
            or not isinstance(generator.get("rotation"), list)
            or len(generator.get("rotation", [])) != 3
            or any(len(row) != 3 for row in generator.get("rotation", []))
            or not isinstance(generator.get("translation"), list)
            or len(generator.get("translation", [])) != 3
            for generator in generators
        ):
            raise GroupDataError("each generator needs a 3x3 rotation and a length-3 translation")
        return cls(
            hall_number=value["hall_number"],
            hall_symbol=value["hall_symbol"],
            choice=value["choice"],
            centering=value["centering"],
            symmorphic=value["symmorphic"],
            operation_count=value["operation_count"],
            generators=tuple(generators),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "hall_number": self.hall_number,
            "hall_symbol": self.hall_symbol,
            "choice": self.choice,
            "centering": self.centering,
            "symmorphic": self.symmorphic,
            "operation_count": self.operation_count,
            "generators": list(self.generators),
        }


@dataclass(frozen=True, slots=True)
class CrystallographicSpaceGroup:
    """One of the 230 crystallographic space groups."""

    ita_number: int
    international_short: str
    international_full: str
    schoenflies: str
    point_group_number: int
    point_group_hm: str
    point_group_schoenflies: str
    crystal_system: str
    centering: str
    symmorphic: bool
    primary_hall_number: int
    hall_settings: tuple[HallSetting, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "CrystallographicSpaceGroup":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid crystallographic space-group record")
        if type(value.get("ita_number")) is not int or not 1 <= value["ita_number"] <= 230:
            raise GroupDataError("ita_number must be an integer from 1 to 230")
        string_fields = (
            "international_short",
            "international_full",
            "schoenflies",
            "point_group_hm",
            "point_group_schoenflies",
            "crystal_system",
        )
        if any(
            not isinstance(value.get(field), str) or not value[field]
            for field in string_fields
        ):
            raise GroupDataError("space-group symbols and crystal system must be strings")
        if value.get("crystal_system") not in {
            "triclinic",
            "monoclinic",
            "orthorhombic",
            "tetragonal",
            "trigonal",
            "hexagonal",
            "cubic",
        }:
            raise GroupDataError("unknown crystal system")
        if value.get("centering") not in {"P", "A", "B", "C", "F", "I", "R"}:
            raise GroupDataError("centering must be one of P, A, B, C, F, I, R")
        if type(value.get("point_group_number")) is not int:
            raise GroupDataError("point_group_number must be an integer")
        if type(value.get("symmorphic")) is not bool:
            raise GroupDataError("symmorphic must be a boolean")
        if type(value.get("primary_hall_number")) is not int:
            raise GroupDataError("primary_hall_number must be an integer")
        settings = value.get("hall_settings")
        if not isinstance(settings, list) or not settings:
            raise GroupDataError("hall_settings must be a non-empty list")
        parsed_settings = tuple(HallSetting.from_mapping(item) for item in settings)
        if any(
            setting.hall_number != value["primary_hall_number"]
            for setting in parsed_settings[:1]
        ) and not any(
            setting.hall_number == value["primary_hall_number"]
            for setting in parsed_settings
        ):
            raise GroupDataError("primary_hall_number must name one hall setting")
        if value["symmorphic"] != any(setting.symmorphic for setting in parsed_settings):
            raise GroupDataError("symmorphic must match the hall settings")
        try:
            return cls(
                ita_number=value["ita_number"],
                international_short=value["international_short"],
                international_full=value["international_full"],
                schoenflies=value["schoenflies"],
                point_group_number=value["point_group_number"],
                point_group_hm=value["point_group_hm"],
                point_group_schoenflies=value["point_group_schoenflies"],
                crystal_system=value["crystal_system"],
                centering=value["centering"],
                symmorphic=value["symmorphic"],
                primary_hall_number=value["primary_hall_number"],
                hall_settings=parsed_settings,
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GroupDataError("invalid crystallographic space-group record") from exc

    def to_dict(self) -> dict[str, Any]:
        return {
            "ita_number": self.ita_number,
            "international_short": self.international_short,
            "international_full": self.international_full,
            "schoenflies": self.schoenflies,
            "point_group_number": self.point_group_number,
            "point_group_hm": self.point_group_hm,
            "point_group_schoenflies": self.point_group_schoenflies,
            "crystal_system": self.crystal_system,
            "centering": self.centering,
            "symmorphic": self.symmorphic,
            "primary_hall_number": self.primary_hall_number,
            "hall_settings": [setting.to_dict() for setting in self.hall_settings],
        }


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "crystallographic_space_groups.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/crystallographic_space_groups.json")
        .read_text(encoding="utf-8")
    )


def load_space_group_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load the versioned 230-space-group registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read space-group registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid space-group registry JSON: {exc}") from exc
    groups = data.get("space_groups")
    if data.get("schema_version") != 1 or not isinstance(groups, list):
        raise GroupDataError("space-group registry must use schema_version 1")
    records = [CrystallographicSpaceGroup.from_mapping(item) for item in groups]
    if [record.ita_number for record in records] != list(range(1, 231)):
        raise GroupDataError("space-group registry must contain the ordered numbers 1-230")
    if len({record.international_short for record in records}) != 230:
        raise GroupDataError("space-group symbols must be unique")
    if len({record.point_group_hm for record in records}) != 32:
        raise GroupDataError("space-group registry must cover all 32 point groups")
    expected_counts = {
        "triclinic": 2,
        "monoclinic": 13,
        "orthorhombic": 59,
        "tetragonal": 68,
        "trigonal": 25,
        "hexagonal": 27,
        "cubic": 36,
    }
    if {system: sum(r.crystal_system == system for r in records) for system in expected_counts} != expected_counts:
        raise GroupDataError("space-group crystal-system counts do not match ITA")
    if sum(record.symmorphic for record in records) != 73:
        raise GroupDataError("expected 73 symmorphic space groups")
    if any(
        record.primary_hall_number not in {s.hall_number for s in record.hall_settings}
        for record in records
    ):
        raise GroupDataError("primary_hall_number must name a stored hall setting")
    return data


def iter_crystallographic_space_groups(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[CrystallographicSpaceGroup]:
    """Iterate the 230 space groups in ITA order."""

    data = load_space_group_registry() if registry is None else registry
    for item in data["space_groups"]:
        yield CrystallographicSpaceGroup.from_mapping(item)


def get_crystallographic_space_group(
    ita_number: int,
    registry: Mapping[str, Any] | None = None,
) -> CrystallographicSpaceGroup:
    """Return one space group by ITA number."""

    if type(ita_number) is not int or not 1 <= ita_number <= 230:
        raise GroupDataError("ita_number must be an integer from 1 to 230")
    data = load_space_group_registry() if registry is None else registry
    for item in data["space_groups"]:
        if item.get("ita_number") == ita_number:
            return CrystallographicSpaceGroup.from_mapping(item)
    raise GroupDataError(f"space group {ita_number} not found in registry")
