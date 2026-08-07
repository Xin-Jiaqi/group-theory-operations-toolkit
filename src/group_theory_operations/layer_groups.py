"""Typed access to the 80 crystallographic layer groups."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .catalog import GroupDataError


_CRYSTAL_SYSTEM_COUNTS = {
    "triclinic": 2,
    "monoclinic": 16,
    "orthorhombic": 30,
    "tetragonal": 16,
    "trigonal": 8,
    "hexagonal": 8,
}


@dataclass(frozen=True, slots=True)
class LayerHallSetting:
    layer_hall_number: int
    setting_number: int
    choice: str
    standard: bool
    hall_symbol: str
    international: str
    international_full: str
    international_short: str
    centering: str
    operation_count: int
    generators: tuple[dict[str, Any], ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "LayerHallSetting":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid layer Hall-setting record")
        if type(value.get("layer_hall_number")) is not int or not 1 <= value[
            "layer_hall_number"
        ] <= 116:
            raise GroupDataError("layer_hall_number must be an integer from 1 to 116")
        if type(value.get("setting_number")) is not int or value["setting_number"] < 1:
            raise GroupDataError("setting_number must be a positive integer")
        if type(value.get("standard")) is not bool:
            raise GroupDataError("standard must be a boolean")
        string_fields = (
            "choice",
            "hall_symbol",
            "international",
            "international_full",
            "international_short",
        )
        if not isinstance(value.get("choice"), str) or any(
            not isinstance(value.get(field), str) or not value[field]
            for field in string_fields[1:]
        ):
            raise GroupDataError("layer-setting symbols must be strings")
        if value.get("centering") not in {"P", "C"}:
            raise GroupDataError("layer centering must be P or C")
        if type(value.get("operation_count")) is not int or not 1 <= value[
            "operation_count"
        ] <= 24:
            raise GroupDataError("operation_count must be an integer from 1 to 24")
        generators = value.get("generators")
        if not isinstance(generators, list):
            raise GroupDataError("generators must be a list")
        for generator in generators:
            if not isinstance(generator, Mapping):
                raise GroupDataError("invalid layer-group generator")
            rotation = generator.get("rotation")
            translation = generator.get("translation")
            if (
                not isinstance(rotation, list)
                or len(rotation) != 3
                or any(not isinstance(row, list) or len(row) != 3 for row in rotation)
                or not isinstance(translation, list)
                or len(translation) != 3
            ):
                raise GroupDataError("layer generators need 3x3 rotations and triplets")
        return cls(
            layer_hall_number=value["layer_hall_number"],
            setting_number=value["setting_number"],
            choice=value["choice"],
            standard=value["standard"],
            hall_symbol=value["hall_symbol"],
            international=value["international"],
            international_full=value["international_full"],
            international_short=value["international_short"],
            centering=value["centering"],
            operation_count=value["operation_count"],
            generators=tuple(dict(generator) for generator in generators),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "layer_hall_number": self.layer_hall_number,
            "setting_number": self.setting_number,
            "choice": self.choice,
            "standard": self.standard,
            "hall_symbol": self.hall_symbol,
            "international": self.international,
            "international_full": self.international_full,
            "international_short": self.international_short,
            "centering": self.centering,
            "operation_count": self.operation_count,
            "generators": [dict(generator) for generator in self.generators],
        }


@dataclass(frozen=True, slots=True)
class CrystallographicLayerGroup:
    number: int
    international_short: str
    international_full: str
    schoenflies: str
    point_group_number: int
    point_group_hm: str
    point_group_schoenflies: str
    point_group_embedding: str | None
    crystal_system: str
    centering: str
    primary_layer_hall_number: int
    hall_settings: tuple[LayerHallSetting, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "CrystallographicLayerGroup":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid crystallographic layer-group record")
        if type(value.get("number")) is not int or not 1 <= value["number"] <= 80:
            raise GroupDataError("layer-group number must be an integer from 1 to 80")
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
            raise GroupDataError("layer-group symbols and crystal system must be strings")
        if value["crystal_system"] not in _CRYSTAL_SYSTEM_COUNTS:
            raise GroupDataError("unknown layer-group crystal system")
        if value.get("centering") not in {"P", "C"}:
            raise GroupDataError("layer centering must be P or C")
        if type(value.get("point_group_number")) is not int or not 1 <= value[
            "point_group_number"
        ] <= 27:
            raise GroupDataError("point_group_number must be an integer from 1 to 27")
        embedding = value.get("point_group_embedding")
        if embedding is not None and (not isinstance(embedding, str) or not embedding):
            raise GroupDataError("point_group_embedding must be a non-empty string")
        settings = value.get("hall_settings")
        if not isinstance(settings, list) or not settings:
            raise GroupDataError("hall_settings must be a non-empty list")
        parsed_settings = tuple(LayerHallSetting.from_mapping(item) for item in settings)
        primary = value.get("primary_layer_hall_number")
        if type(primary) is not int or primary not in {
            setting.layer_hall_number for setting in parsed_settings
        }:
            raise GroupDataError("primary_layer_hall_number must name one setting")
        if sum(setting.standard for setting in parsed_settings) != 1:
            raise GroupDataError("each layer group must have exactly one standard setting")
        return cls(
            number=value["number"],
            international_short=value["international_short"],
            international_full=value["international_full"],
            schoenflies=value["schoenflies"],
            point_group_number=value["point_group_number"],
            point_group_hm=value["point_group_hm"],
            point_group_schoenflies=value["point_group_schoenflies"],
            point_group_embedding=embedding,
            crystal_system=value["crystal_system"],
            centering=value["centering"],
            primary_layer_hall_number=primary,
            hall_settings=parsed_settings,
        )

    def to_dict(self) -> dict[str, Any]:
        result = {
            "number": self.number,
            "international_short": self.international_short,
            "international_full": self.international_full,
            "schoenflies": self.schoenflies,
            "point_group_number": self.point_group_number,
            "point_group_hm": self.point_group_hm,
            "point_group_schoenflies": self.point_group_schoenflies,
            "crystal_system": self.crystal_system,
            "centering": self.centering,
            "primary_layer_hall_number": self.primary_layer_hall_number,
            "hall_settings": [setting.to_dict() for setting in self.hall_settings],
        }
        if self.point_group_embedding is not None:
            result["point_group_embedding"] = self.point_group_embedding
        return result


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "crystallographic_layer_groups.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/crystallographic_layer_groups.json")
        .read_text(encoding="utf-8")
    )


def load_layer_group_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the versioned 80-layer-group registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read layer-group registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid layer-group registry JSON: {exc}") from exc
    groups = data.get("layer_groups")
    if data.get("schema_version") != 1 or not isinstance(groups, list):
        raise GroupDataError("layer-group registry must use schema_version 1")
    records = [CrystallographicLayerGroup.from_mapping(item) for item in groups]
    if [record.number for record in records] != list(range(1, 81)):
        raise GroupDataError("layer-group registry must contain ordered LG1-LG80")
    if len({record.international_short for record in records}) != 80:
        raise GroupDataError("standard layer-group symbols must be unique")
    settings = [setting for record in records for setting in record.hall_settings]
    if [setting.layer_hall_number for setting in settings] != list(range(1, 117)):
        raise GroupDataError("layer Hall settings must contain ordered numbers 1-116")
    actual_counts = {
        system: sum(record.crystal_system == system for record in records)
        for system in _CRYSTAL_SYSTEM_COUNTS
    }
    if actual_counts != _CRYSTAL_SYSTEM_COUNTS:
        raise GroupDataError("layer-group crystal-system counts are invalid")
    if len({record.point_group_hm for record in records}) != 27:
        raise GroupDataError("layer groups must cover the 27 axial point groups")
    return data


def iter_crystallographic_layer_groups(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[CrystallographicLayerGroup]:
    data = load_layer_group_registry() if registry is None else registry
    for item in data["layer_groups"]:
        yield CrystallographicLayerGroup.from_mapping(item)


def _normalize_symbol(value: str) -> str:
    return "".join(value.lower().replace("_", "").split())


def get_crystallographic_layer_group(
    identifier: int | str,
    registry: Mapping[str, Any] | None = None,
) -> CrystallographicLayerGroup:
    """Return one layer group by LG number or international symbol."""

    records = tuple(iter_crystallographic_layer_groups(registry))
    if type(identifier) is int:
        if not 1 <= identifier <= 80:
            raise GroupDataError("layer-group number must be from 1 to 80")
        for record in records:
            if record.number == identifier:
                return record
    elif isinstance(identifier, str) and identifier.strip():
        cleaned = identifier.strip()
        if cleaned.isdigit():
            return get_crystallographic_layer_group(int(cleaned), registry)
        if cleaned.upper().startswith("LG") and cleaned[2:].isdigit():
            return get_crystallographic_layer_group(int(cleaned[2:]), registry)
        target = _normalize_symbol(identifier)
        matches = [
            record
            for record in records
            if target
            in {
                _normalize_symbol(record.international_short),
                _normalize_symbol(record.international_full),
                *(
                    _normalize_symbol(setting.international_short)
                    for setting in record.hall_settings
                ),
            }
        ]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise GroupDataError(f"layer-group symbol {identifier!r} is ambiguous")
    else:
        raise GroupDataError("layer-group identifier must be a number or symbol")
    raise GroupDataError(f"layer group {identifier!r} not found")
