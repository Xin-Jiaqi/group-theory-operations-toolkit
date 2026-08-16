"""Complete Hall-setting Wyckoff registry and orbit splitting."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from itertools import product
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .catalog import GroupDataError
from .seitz import SeitzOp, closure, equivalent, transform_seitz_coordinates
from .space_groups import get_crystallographic_space_group


Vector3 = tuple[float, float, float]
IntegerVector3 = tuple[int, int, int]
IntegerMatrix3 = tuple[tuple[int, int, int], ...]
FloatVector3 = tuple[float, float, float]
FloatMatrix3 = tuple[tuple[float, float, float], ...]
_MAX_CONVENTIONAL_CELL_INDEX = 4096
_MAX_EXPANDED_ORBIT_SIZE = 4096


def _vector3(value: Iterable[Any], *, label: str) -> IntegerVector3:
    values = tuple(value)
    if len(values) != 3 or any(type(item) is not int for item in values):
        raise GroupDataError(f"{label} must contain three integers")
    return int(values[0]), int(values[1]), int(values[2])


def _matrix3(value: Iterable[Any]) -> IntegerMatrix3:
    rows = tuple(_vector3(row, label="parameter-matrix row") for row in value)
    if len(rows) != 3:
        raise GroupDataError("parameter_matrix must contain three rows")
    return rows


def _parameters(value: Iterable[Any]) -> np.ndarray:
    try:
        result = np.asarray(tuple(value), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise GroupDataError("Wyckoff parameters must be three finite numbers") from exc
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise GroupDataError("Wyckoff parameters must be three finite numbers")
    return result


def _integer_determinant_and_adjugate(
    matrix: np.ndarray,
) -> tuple[int, IntegerMatrix3]:
    rows = tuple(tuple(int(value) for value in row) for row in matrix)
    a, b, c = rows[0]
    d, e, f = rows[1]
    g, h, i = rows[2]
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (
        d * h - e * g
    )
    adjugate: IntegerMatrix3 = (
        (e * i - f * h, c * h - b * i, b * f - c * e),
        (f * g - d * i, a * i - c * g, c * d - a * f),
        (d * h - e * g, b * g - a * h, a * e - b * d),
    )
    return determinant, adjugate


def _coordinate_transformation(
    transformation_matrix: Iterable[Iterable[Any]] | None,
    origin_shift: Iterable[Any] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[IntegerVector3, ...]]:
    try:
        matrix = (
            np.eye(3, dtype=np.float64)
            if transformation_matrix is None
            else np.asarray(tuple(tuple(row) for row in transformation_matrix), dtype=float)
        )
        shift = (
            np.zeros(3, dtype=np.float64)
            if origin_shift is None
            else np.asarray(tuple(origin_shift), dtype=float)
        )
    except (TypeError, ValueError) as exc:
        raise GroupDataError("invalid subgroup coordinate transformation") from exc
    if matrix.shape != (3, 3) or shift.shape != (3,):
        raise GroupDataError(
            "subgroup transformation must be a 3x3 matrix and a length-3 shift"
        )
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(shift)):
        raise GroupDataError("subgroup coordinate transformation must be finite")
    try:
        inverse_matrix = np.linalg.inv(matrix)
    except np.linalg.LinAlgError as exc:
        raise GroupDataError("subgroup transformation matrix must be invertible") from exc
    rounded_inverse = np.rint(inverse_matrix)
    if not np.allclose(inverse_matrix, rounded_inverse, atol=1.0e-9, rtol=0.0):
        raise GroupDataError(
            "P^-1 must be an integer parent-cell supercell matrix for "
            "x_subgroup = P x_parent + p"
        )
    if float(np.max(np.abs(rounded_inverse))) > float(np.iinfo(np.int32).max):
        raise GroupDataError("subgroup supercell-matrix entries are too large")
    supercell_matrix = rounded_inverse.astype(np.int64)
    determinant, adjugate = _integer_determinant_and_adjugate(supercell_matrix)
    if determinant == 0:  # pragma: no cover - guarded by invertibility above
        raise GroupDataError("subgroup supercell matrix must be invertible")
    if abs(determinant) > _MAX_CONVENTIONAL_CELL_INDEX:
        raise GroupDataError(
            "subgroup conventional-cell index exceeds the supported maximum "
            f"of {_MAX_CONVENTIONAL_CELL_INDEX}"
        )
    canonical_matrix = np.asarray(adjugate, dtype=np.float64) / determinant
    if not np.allclose(matrix, canonical_matrix, atol=1.0e-9, rtol=0.0):
        raise GroupDataError("subgroup transformation is not an exact supercell map")
    cosets = _translation_coset_representatives(supercell_matrix)
    return canonical_matrix, shift, supercell_matrix, cosets


def _translation_coset_representatives(
    supercell_matrix: np.ndarray,
) -> tuple[IntegerVector3, ...]:
    """Return representatives of ``Z^3 / A Z^3`` for integer ``A``.

    The representatives are parent-cell lattice translations.  Their images
    ``P n`` with ``P = A^-1`` place every copy of a parent conventional cell
    inside the subgroup conventional cell.  Residues are evaluated with the
    exact integer adjugate, avoiding a tolerance-dependent coset partition.
    """

    matrix = np.asarray(supercell_matrix, dtype=np.int64)
    determinant, adjugate = _integer_determinant_and_adjugate(matrix)
    cell_index = abs(determinant)

    def residue(vector: np.ndarray) -> IntegerVector3:
        values = tuple(
            sum(row[column] * int(vector[column]) for column in range(3))
            % cell_index
            for row in adjugate
        )
        return int(values[0]), int(values[1]), int(values[2])

    zero = np.zeros(3, dtype=np.int64)
    representatives = [zero]
    seen = {residue(zero)}
    cursor = 0
    generators = np.eye(3, dtype=np.int64)
    while cursor < len(representatives) and len(representatives) < cell_index:
        current = representatives[cursor]
        cursor += 1
        for generator in generators:
            candidate = current + generator
            key = residue(candidate)
            if key in seen:
                continue
            seen.add(key)
            representatives.append(candidate)
            if len(representatives) == cell_index:
                break
    if len(representatives) != cell_index:  # pragma: no cover - algebraic guard
        raise GroupDataError("could not construct all supercell translation cosets")
    return tuple(
        (int(item[0]), int(item[1]), int(item[2]))
        for item in representatives
    )


def _periodic_close(left: Any, right: Any, tolerance: float) -> bool:
    difference = np.asarray(left, dtype=np.float64) - np.asarray(
        right, dtype=np.float64
    )
    difference -= np.rint(difference)
    return bool(np.max(np.abs(difference)) <= tolerance)


def _unique_positions(values: Iterable[Any], tolerance: float) -> list[np.ndarray]:
    unique: list[np.ndarray] = []
    for value in values:
        position = np.asarray(value, dtype=np.float64) % 1.0
        if not any(_periodic_close(position, known, tolerance) for known in unique):
            unique.append(position)
    return unique


@dataclass(frozen=True, slots=True)
class WyckoffCoordinateMap:
    """One affine coordinate map ``x = A q + t/24``."""

    parameter_matrix: IntegerMatrix3
    translation_numerators: IntegerVector3

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "WyckoffCoordinateMap":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid Wyckoff coordinate map")
        try:
            matrix = _matrix3(value["parameter_matrix"])
            translation = _vector3(
                value["translation_numerators"],
                label="translation_numerators",
            )
        except (KeyError, TypeError) as exc:
            raise GroupDataError("invalid Wyckoff coordinate map") from exc
        if any(not 0 <= item < 24 for item in translation):
            raise GroupDataError("Wyckoff translation numerators must be in 0..23")
        return cls(matrix, translation)

    def evaluate(self, parameters: Iterable[Any]) -> Vector3:
        q = _parameters(parameters)
        position = (
            np.asarray(self.parameter_matrix, dtype=np.float64) @ q
            + np.asarray(self.translation_numerators, dtype=np.float64) / 24.0
        ) % 1.0
        return float(position[0]), float(position[1]), float(position[2])

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameter_matrix": [list(row) for row in self.parameter_matrix],
            "translation_numerators": list(self.translation_numerators),
        }


@dataclass(frozen=True, slots=True)
class WyckoffPositionRecord:
    """One Wyckoff position in a specified Hall setting."""

    multiplicity: int
    letter: str
    site_symmetry: str
    parameter_dimension: int
    representative_maps: tuple[WyckoffCoordinateMap, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "WyckoffPositionRecord":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid Wyckoff-position record")
        multiplicity = value.get("multiplicity")
        dimension = value.get("parameter_dimension")
        if type(multiplicity) is not int or not 1 <= multiplicity <= 192:
            raise GroupDataError("Wyckoff multiplicity must be an integer in 1..192")
        if type(dimension) is not int or dimension not in {0, 1, 2, 3}:
            raise GroupDataError("Wyckoff parameter dimension must be 0, 1, 2, or 3")
        letter = value.get("letter")
        site_symmetry = value.get("site_symmetry")
        if not isinstance(letter, str) or len(letter) != 1 or not letter.isalpha():
            raise GroupDataError("Wyckoff letter must be one alphabetic character")
        if not isinstance(site_symmetry, str) or not site_symmetry:
            raise GroupDataError("site-symmetry symbol must be non-empty")
        raw_maps = value.get("representative_maps")
        if not isinstance(raw_maps, list) or not raw_maps:
            raise GroupDataError("Wyckoff position needs representative maps")
        maps = tuple(WyckoffCoordinateMap.from_mapping(item) for item in raw_maps)
        ranks = {
            int(np.linalg.matrix_rank(np.asarray(item.parameter_matrix, dtype=float)))
            for item in maps
        }
        if ranks != {dimension}:
            raise GroupDataError("parameter dimension disagrees with coordinate maps")
        return cls(multiplicity, letter, site_symmetry, dimension, maps)

    @property
    def label(self) -> str:
        return f"{self.multiplicity}{self.letter}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "multiplicity": self.multiplicity,
            "letter": self.letter,
            "label": self.label,
            "site_symmetry": self.site_symmetry,
            "parameter_dimension": self.parameter_dimension,
            "representative_maps": [item.to_dict() for item in self.representative_maps],
        }


@dataclass(frozen=True, slots=True)
class HallWyckoffSetting:
    """All Wyckoff positions for one origin- and axis-specific Hall setting."""

    hall_number: int
    ita_number: int
    hall_symbol: str
    choice: str
    source_symbol: str
    centering: str
    centering_translation_numerators: tuple[IntegerVector3, ...]
    positions: tuple[WyckoffPositionRecord, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "HallWyckoffSetting":
        if not isinstance(value, Mapping):
            raise GroupDataError("invalid Hall Wyckoff setting")
        hall_number = value.get("hall_number")
        ita_number = value.get("ita_number")
        if type(hall_number) is not int or not 1 <= hall_number <= 530:
            raise GroupDataError("Hall number must be an integer in 1..530")
        if type(ita_number) is not int or not 1 <= ita_number <= 230:
            raise GroupDataError("ITA number must be an integer in 1..230")
        string_fields = ("hall_symbol", "choice", "source_symbol", "centering")
        if any(not isinstance(value.get(field), str) for field in string_fields):
            raise GroupDataError("Hall-setting symbols must be strings")
        if not value["hall_symbol"] or not value["source_symbol"]:
            raise GroupDataError("Hall-setting symbols must be non-empty")
        if value["centering"] not in {"P", "A", "B", "C", "F", "I", "R", "H"}:
            raise GroupDataError("invalid Wyckoff centering symbol")
        raw_centering = value.get("centering_translation_numerators")
        raw_positions = value.get("positions")
        if not isinstance(raw_centering, list) or not raw_centering:
            raise GroupDataError("Hall setting needs centering translations")
        if not isinstance(raw_positions, list) or not raw_positions:
            raise GroupDataError("Hall setting needs Wyckoff positions")
        centerings = tuple(
            _vector3(item, label="centering translation") for item in raw_centering
        )
        if any(any(not 0 <= component < 24 for component in item) for item in centerings):
            raise GroupDataError("centering translation numerators must be in 0..23")
        positions = tuple(WyckoffPositionRecord.from_mapping(item) for item in raw_positions)
        if len({item.letter for item in positions}) != len(positions):
            raise GroupDataError("Wyckoff letters must be unique within a Hall setting")
        if any(
            item.multiplicity != len(item.representative_maps) * len(centerings)
            for item in positions
        ):
            raise GroupDataError("Wyckoff multiplicity disagrees with centering expansion")
        return cls(
            hall_number,
            ita_number,
            value["hall_symbol"],
            value["choice"],
            value["source_symbol"],
            value["centering"],
            centerings,
            positions,
        )

    def position(self, letter: str) -> WyckoffPositionRecord:
        if not isinstance(letter, str) or len(letter) != 1:
            raise GroupDataError("Wyckoff letter must be one character")
        for position in self.positions:
            if position.letter == letter:
                return position
        raise GroupDataError(
            f"Hall setting {self.hall_number} has no Wyckoff position {letter!r}"
        )

    def coordinates(
        self,
        letter: str,
        parameters: Iterable[Any] = (0.173, 0.287, 0.419),
        *,
        tolerance: float = 1.0e-9,
    ) -> tuple[Vector3, ...]:
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise GroupDataError("tolerance must be positive and finite")
        record = self.position(letter)
        q = _parameters(parameters)
        values = []
        for coordinate_map in record.representative_maps:
            base = np.asarray(coordinate_map.evaluate(q), dtype=np.float64)
            for translation in self.centering_translation_numerators:
                values.append(base + np.asarray(translation, dtype=float) / 24.0)
        return tuple(
            (float(item[0]), float(item[1]), float(item[2]))
            for item in _unique_positions(values, tolerance)
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "hall_number": self.hall_number,
            "ita_number": self.ita_number,
            "hall_symbol": self.hall_symbol,
            "choice": self.choice,
            "source_symbol": self.source_symbol,
            "centering": self.centering,
            "centering_translation_numerators": [
                list(item) for item in self.centering_translation_numerators
            ],
            "positions": [item.to_dict() for item in self.positions],
        }


def _default_text() -> str:
    repository_copy = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "crystallographic_wyckoff_positions.json"
    )
    if repository_copy.is_file():
        return repository_copy.read_text(encoding="utf-8")
    return (
        resources.files("group_theory_operations")
        .joinpath("data/crystallographic_wyckoff_positions.json")
        .read_text(encoding="utf-8")
    )


def load_wyckoff_registry(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the complete 530-Hall-setting Wyckoff registry."""

    try:
        text = _default_text() if path is None else Path(path).read_text(encoding="utf-8")
        data = json.loads(text)
    except OSError as exc:
        raise GroupDataError(f"cannot read Wyckoff registry: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"invalid Wyckoff registry JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise GroupDataError("Wyckoff registry root must be a JSON object")
    settings = data.get("hall_settings")
    expected_counts = {
        "hall_settings": 530,
        "wyckoff_positions": 3467,
        "representative_maps": 15117,
        "expanded_coordinate_maps": 24295,
    }
    if data.get("schema_version") != 1 or not isinstance(settings, list):
        raise GroupDataError("Wyckoff registry must use schema_version 1")
    source = data.get("source")
    expected_source = {
        "name": "spglib",
        "version": "2.5.0",
        "tag": "v2.5.0",
        "commit": "e4531bb49371dce3e807c2095a4d9d9b7245c524",
        "license": "BSD-3-Clause",
        "file": "database/Wyckoff.csv",
        "sha256": "d3d786a1f0187e5c6d69a3ade35648ffab34fd1b977d61ad84d8b0434b8b7ca0",
        "url": (
            "https://github.com/spglib/spglib/blob/"
            "e4531bb49371dce3e807c2095a4d9d9b7245c524/database/Wyckoff.csv"
        ),
    }
    if not isinstance(source, Mapping) or dict(source) != expected_source:
        raise GroupDataError("Wyckoff registry source provenance is invalid")
    if data.get("counts") != expected_counts:
        raise GroupDataError("Wyckoff registry dimensions do not match v1")
    records = tuple(HallWyckoffSetting.from_mapping(item) for item in settings)
    if tuple(item.hall_number for item in records) != tuple(range(1, 531)):
        raise GroupDataError("Wyckoff registry must contain ordered Hall numbers 1-530")
    if sum(len(item.positions) for item in records) != 3467:
        raise GroupDataError("Wyckoff-position count disagrees with registry metadata")
    if sum(
        len(position.representative_maps)
        for setting in records
        for position in setting.positions
    ) != 15117:
        raise GroupDataError("representative-map count disagrees with registry metadata")
    if sum(
        position.multiplicity
        for setting in records
        for position in setting.positions
    ) != 24295:
        raise GroupDataError("expanded-coordinate count disagrees with registry metadata")
    return data


@lru_cache(maxsize=1)
def _default_settings() -> tuple[HallWyckoffSetting, ...]:
    data = load_wyckoff_registry()
    return tuple(
        HallWyckoffSetting.from_mapping(item) for item in data["hall_settings"]
    )


def iter_wyckoff_settings(
    registry: Mapping[str, Any] | None = None,
) -> Iterable[HallWyckoffSetting]:
    if registry is None:
        yield from _default_settings()
        return
    data = registry
    for item in data["hall_settings"]:
        yield HallWyckoffSetting.from_mapping(item)


def get_wyckoff_setting(
    hall_number: int,
    registry: Mapping[str, Any] | None = None,
) -> HallWyckoffSetting:
    if type(hall_number) is not int or not 1 <= hall_number <= 530:
        raise GroupDataError("Hall number must be an integer in 1..530")
    if registry is None:
        return _default_settings()[hall_number - 1]
    data = registry
    return HallWyckoffSetting.from_mapping(data["hall_settings"][hall_number - 1])


@lru_cache(maxsize=530)
def _hall_operations(hall_number: int) -> tuple[SeitzOp, ...]:
    setting = get_wyckoff_setting(hall_number)
    space_group = get_crystallographic_space_group(setting.ita_number)
    hall = next(
        (item for item in space_group.hall_settings if item.hall_number == hall_number),
        None,
    )
    if hall is None:  # pragma: no cover - cross-registry contract guard
        raise GroupDataError("Wyckoff and space-group registries disagree")
    return tuple(closure(SeitzOp.from_dict(item) for item in hall.generators))


def _operation_subset(subgroup: Iterable[SeitzOp], parent: Iterable[SeitzOp]) -> bool:
    parent_operations = tuple(parent)
    return all(
        any(equivalent(operation, candidate) for candidate in parent_operations)
        for operation in subgroup
    )


def _coordinate_orbits(
    coordinates: list[np.ndarray], operations: tuple[SeitzOp, ...], tolerance: float
) -> list[list[np.ndarray]]:
    remaining = set(range(len(coordinates)))
    orbits: list[list[np.ndarray]] = []
    while remaining:
        seed = min(remaining)
        members: set[int] = set()
        for operation in operations:
            image = operation.apply(coordinates[seed])
            matches = [
                index
                for index, coordinate in enumerate(coordinates)
                if _periodic_close(image, coordinate, tolerance)
            ]
            if len(matches) != 1:
                raise GroupDataError(
                    "subgroup operation does not map uniquely within the parent orbit"
                )
            members.add(matches[0])
        if not members or not members.issubset(remaining | {seed}):
            raise GroupDataError("subgroup operations do not permute the parent orbit")
        remaining -= members
        orbits.append([coordinates[index] for index in sorted(members)])
    return orbits


def _fit_parameters(
    coordinate_map: WyckoffCoordinateMap,
    target: np.ndarray,
    tolerance: float,
) -> np.ndarray | None:
    matrix = np.asarray(coordinate_map.parameter_matrix, dtype=np.float64)
    translation = np.asarray(
        coordinate_map.translation_numerators, dtype=np.float64
    ) / 24.0
    for lattice_shift in product(range(-2, 3), repeat=3):
        right = target - translation + np.asarray(lattice_shift, dtype=float)
        parameters, _, _, _ = np.linalg.lstsq(matrix, right, rcond=None)
        if np.max(np.abs(matrix @ parameters - right)) <= tolerance:
            return parameters
    return None


def _same_coordinate_set(
    left: Iterable[Any], right: Iterable[Any], tolerance: float
) -> bool:
    left_values = tuple(left)
    right_values = tuple(right)
    return len(left_values) == len(right_values) and all(
        any(_periodic_close(item, candidate, tolerance) for candidate in right_values)
        for item in left_values
    )


@dataclass(frozen=True, slots=True)
class SubgroupWyckoffOrbit:
    """One child Wyckoff orbit produced by symmetry reduction."""

    letter: str
    multiplicity: int
    site_symmetry: str
    parameter_dimension: int
    coordinates: tuple[Vector3, ...]

    @property
    def label(self) -> str:
        return f"{self.multiplicity}{self.letter}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "letter": self.letter,
            "multiplicity": self.multiplicity,
            "label": self.label,
            "site_symmetry": self.site_symmetry,
            "parameter_dimension": self.parameter_dimension,
            "coordinates": [list(item) for item in self.coordinates],
        }


@dataclass(frozen=True, slots=True)
class WyckoffOrbitSplitting:
    """Decomposition of one parent orbit under an embedded subgroup."""

    parent_hall_number: int
    parent_label: str
    subgroup_hall_number: int
    subgroup_index: int
    conventional_cell_index: int
    translation_subgroup_index: int
    point_group_index: int
    subgroup_supercell_matrix: IntegerMatrix3
    parent_translation_cosets: tuple[IntegerVector3, ...]
    subgroup_transformation_matrix: FloatMatrix3
    subgroup_origin_shift: FloatVector3
    child_orbits: tuple[SubgroupWyckoffOrbit, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent_hall_number": self.parent_hall_number,
            "parent_label": self.parent_label,
            "subgroup_hall_number": self.subgroup_hall_number,
            "subgroup_index": self.subgroup_index,
            "conventional_cell_index": self.conventional_cell_index,
            "translation_subgroup_index": self.translation_subgroup_index,
            "point_group_index": self.point_group_index,
            "coordinate_convention": "x_subgroup = P @ x_parent + p",
            "subgroup_supercell_matrix": [
                list(row) for row in self.subgroup_supercell_matrix
            ],
            "parent_translation_cosets": [
                list(item) for item in self.parent_translation_cosets
            ],
            "subgroup_transformation_matrix": [
                list(row) for row in self.subgroup_transformation_matrix
            ],
            "subgroup_origin_shift": list(self.subgroup_origin_shift),
            "child_orbits": [item.to_dict() for item in self.child_orbits],
        }


def split_wyckoff_orbit(
    parent_hall_number: int,
    parent_letter: str,
    subgroup_hall_number: int,
    parameters: Iterable[Any] = (0.173, 0.287, 0.419),
    *,
    subgroup_transformation_matrix: Iterable[Iterable[Any]] | None = None,
    subgroup_origin_shift: Iterable[Any] | None = None,
    tolerance: float = 1.0e-8,
) -> WyckoffOrbitSplitting:
    """Split a parent Wyckoff orbit under an explicitly embedded subgroup.

    The optional coordinate relation follows ``x_subgroup = P x_parent + p``.
    Identity ``P`` and zero ``p`` are used by default.  ``P`` must be an
    matrix whose inverse ``A = P^-1`` is an integer supercell matrix in the
    parent conventional basis.  The translation cosets ``Z^3 / A Z^3`` are
    constructed exactly and used to expand the parent orbit into the subgroup
    conventional cell.  If the transformed subgroup operations are not a
    subset of the parent operations, the pair is rejected; the function never
    infers an unstated change of basis, origin, or group--subgroup relation.
    """

    if not np.isfinite(tolerance) or tolerance <= 0:
        raise GroupDataError("tolerance must be positive and finite")
    parent = get_wyckoff_setting(parent_hall_number)
    subgroup = get_wyckoff_setting(subgroup_hall_number)
    transformation, origin_shift, supercell_matrix, translation_cosets = (
        _coordinate_transformation(
            subgroup_transformation_matrix,
            subgroup_origin_shift,
        )
    )
    parent_position = parent.position(parent_letter)
    parent_coordinates = [
        np.asarray(item, dtype=float)
        for item in parent.coordinates(parent_letter, parameters, tolerance=tolerance)
    ]
    if len(parent_coordinates) != parent_position.multiplicity:
        raise GroupDataError(
            "parameters specialize the parent position; query the resulting special "
            "Wyckoff letter before splitting"
        )
    parent_operations = _hall_operations(parent_hall_number)
    native_subgroup_operations = _hall_operations(subgroup_hall_number)
    inverse_transformation = np.linalg.inv(transformation)
    inverse_origin_shift = -inverse_transformation @ origin_shift
    try:
        subgroup_operations = tuple(
            transform_seitz_coordinates(
                operation,
                inverse_transformation,
                inverse_origin_shift,
            )
            for operation in native_subgroup_operations
        )
    except ValueError as exc:
        raise GroupDataError(
            "subgroup coordinate transformation is incompatible with its Hall operations"
        ) from exc
    unique_subgroup_operations: list[SeitzOp] = []
    for operation in subgroup_operations:
        if not any(
            equivalent(operation, known) for known in unique_subgroup_operations
        ):
            unique_subgroup_operations.append(operation)
    if len(unique_subgroup_operations) != len(subgroup_operations):
        raise GroupDataError(
            "subgroup coordinate transformation collapses distinct Hall operations"
        )
    if not _operation_subset(subgroup_operations, parent_operations):
        raise GroupDataError(
            "subgroup Hall operations are not embedded in the parent setting; "
            "an explicit basis/origin transformation is required"
        )
    conventional_cell_index = len(translation_cosets)
    index_numerator = conventional_cell_index * len(parent_operations)
    if index_numerator % len(subgroup_operations):
        raise GroupDataError("parent and subgroup operation orders have noninteger index")
    centering_numerator = (
        conventional_cell_index * len(parent.centering_translation_numerators)
    )
    if centering_numerator % len(subgroup.centering_translation_numerators):
        raise GroupDataError(
            "parent and subgroup translation lattices have noninteger index"
        )
    translation_subgroup_index = (
        centering_numerator // len(subgroup.centering_translation_numerators)
    )
    parent_point_order = (
        len(parent_operations) // len(parent.centering_translation_numerators)
    )
    subgroup_point_order = (
        len(subgroup_operations) // len(subgroup.centering_translation_numerators)
    )
    if parent_point_order % subgroup_point_order:
        raise GroupDataError("parent and subgroup point groups have noninteger index")
    point_group_index = parent_point_order // subgroup_point_order
    subgroup_index = index_numerator // len(subgroup_operations)
    if subgroup_index != translation_subgroup_index * point_group_index:
        raise GroupDataError(
            "space-group index does not factor into lattice and point indices"
        )

    expected_coordinate_count = (
        parent_position.multiplicity * conventional_cell_index
    )
    if expected_coordinate_count > _MAX_EXPANDED_ORBIT_SIZE:
        raise GroupDataError(
            "expanded parent Wyckoff orbit exceeds the supported maximum "
            f"of {_MAX_EXPANDED_ORBIT_SIZE} coordinates"
        )
    expanded_subgroup_coordinates = _unique_positions(
        (
            transformation @ (coordinate + np.asarray(coset, dtype=float))
            + origin_shift
            for coset in translation_cosets
            for coordinate in parent_coordinates
        ),
        tolerance,
    )
    if len(expanded_subgroup_coordinates) != expected_coordinate_count:
        raise GroupDataError(
            "supercell translation cosets do not produce the expected parent-orbit copies"
        )
    coordinate_orbits = _coordinate_orbits(
        expanded_subgroup_coordinates,
        native_subgroup_operations,
        tolerance,
    )
    children: list[SubgroupWyckoffOrbit] = []
    for orbit in coordinate_orbits:
        subgroup_orbit = orbit
        matched: tuple[WyckoffPositionRecord, tuple[Vector3, ...]] | None = None
        for candidate in subgroup.positions:
            if candidate.multiplicity != len(subgroup_orbit):
                continue
            for coordinate_map in candidate.representative_maps:
                for centering in subgroup.centering_translation_numerators:
                    target = (
                        subgroup_orbit[0]
                        - np.asarray(centering, dtype=float) / 24.0
                    )
                    fitted = _fit_parameters(coordinate_map, target, tolerance)
                    if fitted is None:
                        continue
                    generated = subgroup.coordinates(
                        candidate.letter, fitted, tolerance=tolerance
                    )
                    if _same_coordinate_set(generated, subgroup_orbit, tolerance):
                        matched = candidate, generated
                        break
                if matched is not None:
                    break
            if matched is not None:
                break
        if matched is None:
            raise GroupDataError(
                "subgroup orbit could not be matched to its Wyckoff registry"
            )
        candidate, generated = matched
        children.append(
            SubgroupWyckoffOrbit(
                candidate.letter,
                candidate.multiplicity,
                candidate.site_symmetry,
                candidate.parameter_dimension,
                generated,
            )
        )
    transformation_tuple: FloatMatrix3 = (
        (
            float(transformation[0, 0]),
            float(transformation[0, 1]),
            float(transformation[0, 2]),
        ),
        (
            float(transformation[1, 0]),
            float(transformation[1, 1]),
            float(transformation[1, 2]),
        ),
        (
            float(transformation[2, 0]),
            float(transformation[2, 1]),
            float(transformation[2, 2]),
        ),
    )
    origin_shift_tuple: FloatVector3 = (
        float(origin_shift[0]),
        float(origin_shift[1]),
        float(origin_shift[2]),
    )
    supercell_matrix_tuple: IntegerMatrix3 = (
        (
            int(supercell_matrix[0, 0]),
            int(supercell_matrix[0, 1]),
            int(supercell_matrix[0, 2]),
        ),
        (
            int(supercell_matrix[1, 0]),
            int(supercell_matrix[1, 1]),
            int(supercell_matrix[1, 2]),
        ),
        (
            int(supercell_matrix[2, 0]),
            int(supercell_matrix[2, 1]),
            int(supercell_matrix[2, 2]),
        ),
    )
    return WyckoffOrbitSplitting(
        parent_hall_number,
        parent_position.label,
        subgroup_hall_number,
        subgroup_index,
        conventional_cell_index,
        translation_subgroup_index,
        point_group_index,
        supercell_matrix_tuple,
        translation_cosets,
        transformation_tuple,
        origin_shift_tuple,
        tuple(children),
    )
