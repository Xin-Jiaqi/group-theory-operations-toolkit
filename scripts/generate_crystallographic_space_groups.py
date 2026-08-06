#!/usr/bin/env python3
"""Generate the canonical registry of the 230 crystallographic space groups.

The registry is machine-generated from the spglib database (BSD-3-Clause,
same license as this repository).  For every ITA number 1..230:

1. Every Hall setting is listed.  The primary setting is the lowest Hall
   number whose full Seitz operation set round-trips: applying the
   operations to three generic positions on a conventional lattice of the
   right crystal system and asking spglib to identify the resulting
   structure returns the same ITA number.
2. A deterministic greedy reduction derives a compact generating set whose
   closure reproduces the setting's full operation count.
3. The symmorphic flag is computed by the fixed-point criterion: the group
   is symmorphic iff some origin shift makes every non-identity rotation
   pure, modulo the cell's lattice translations (identity-rotation
   translations, including centering vectors).

Run with spglib installed:

    python3 scripts/generate_crystallographic_space_groups.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from spglib import (
    get_symmetry,
    get_symmetry_dataset,
    get_symmetry_from_database,
    get_spacegroup_type,
)


ROOT = Path(__file__).resolve().parents[1]

DATA_PATH = ROOT / "data" / "crystallographic_space_groups.json"
POINT_GROUP_PATH = ROOT / "data" / "crystallographic_point_groups.json"

CRYSTAL_SYSTEMS = (
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
)

EXPECTED_COUNTS = {
    "triclinic": 2,
    "monoclinic": 13,
    "orthorhombic": 59,
    "tetragonal": 68,
    "trigonal": 25,
    "hexagonal": 27,
    "cubic": 36,
}

EXPECTED_SYMMORPHIC = 73

_POINT_GROUP_CRYSTAL_SYSTEM = {
    "1": "triclinic",
    "-1": "triclinic",
    "2": "monoclinic",
    "m": "monoclinic",
    "2/m": "monoclinic",
    "222": "orthorhombic",
    "mm2": "orthorhombic",
    "mmm": "orthorhombic",
    "4": "tetragonal",
    "-4": "tetragonal",
    "4/m": "tetragonal",
    "422": "tetragonal",
    "4mm": "tetragonal",
    "-42m": "tetragonal",
    "4/mmm": "tetragonal",
    "3": "trigonal",
    "-3": "trigonal",
    "32": "trigonal",
    "3m": "trigonal",
    "-3m": "trigonal",
    "6": "hexagonal",
    "-6": "hexagonal",
    "6/m": "hexagonal",
    "622": "hexagonal",
    "6mm": "hexagonal",
    "-6m2": "hexagonal",
    "6/mmm": "hexagonal",
    "23": "cubic",
    "m-3": "cubic",
    "432": "cubic",
    "-43m": "cubic",
    "m-3m": "cubic",
}

_GENERIC_POSITIONS = (
    np.array([0.113, 0.237, 0.341]),
    np.array([0.154, 0.304, 0.424]),
    np.array([0.234, 0.391, 0.439]),
)


def _lattice_for(crystal_system: str) -> np.ndarray:
    if crystal_system == "cubic":
        return np.eye(3) * 5.2
    if crystal_system == "tetragonal":
        return np.array([[5.2, 0.0, 0.0], [0.0, 5.2, 0.0], [0.0, 0.0, 6.3]])
    if crystal_system in ("hexagonal", "trigonal"):
        a, c = 5.2, 6.3
        return np.array([[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]])
    if crystal_system == "orthorhombic":
        return np.array([[5.2, 0.0, 0.0], [0.0, 6.1, 0.0], [0.0, 0.0, 4.4]])
    if crystal_system == "monoclinic":
        return np.array([[5.1, 0.0, 0.0], [0.0, 6.3, 0.0], [0.7, 0.0, 4.2]])
    return np.array([[4.1, 0.0, 0.0], [0.3, 5.3, 0.0], [0.2, 0.4, 6.7]])


def _round_trip_structure(hall_number: int, crystal_system: str):
    """Build a conventional-cell structure whose symmetry is the space group."""
    lattice = _lattice_for(crystal_system)
    operations = get_symmetry_from_database(hall_number)
    rotations = np.asarray(operations["rotations"], dtype=np.int64)
    translations = np.asarray(operations["translations"], dtype=np.float64)
    positions = []
    for base in _GENERIC_POSITIONS:
        for rotation, translation in zip(rotations, translations):
            positions.append((rotation @ base + translation) % 1.0)
    return lattice, positions


def _round_trips(hall_number: int, ita_number: int, crystal_system: str) -> bool:
    lattice, positions = _round_trip_structure(hall_number, crystal_system)
    dataset = get_symmetry_dataset(
        (lattice, positions, [1] * len(positions)), symprec=1.0e-4
    )
    return int(dataset.number) == ita_number


def _lattice_translations(rotations: np.ndarray, translations: np.ndarray) -> list[np.ndarray]:
    identity = np.eye(3, dtype=np.int64)
    return [
        np.asarray(translation, dtype=np.float64) % 1.0
        for rotation, translation in zip(rotations, translations)
        if np.array_equal(rotation, identity)
    ]


def _in_lattice(vector: np.ndarray, generators: list[np.ndarray]) -> bool:
    from itertools import product

    residue = np.asarray(vector, dtype=np.float64) % 1.0
    for coefficients in product(range(-3, 4), repeat=len(generators)):
        combination = sum(c * g for c, g in zip(coefficients, generators)) % 1.0
        difference = combination - residue
        difference -= np.rint(difference)
        if np.max(np.abs(difference)) < 1.0e-9:
            return True
    return False


def _symmorphic(rotations: np.ndarray, translations: np.ndarray) -> bool:
    """Fixed-point criterion: an origin shift makes every rotation pure."""
    identity = np.eye(3, dtype=np.int64)
    lattice = _lattice_translations(rotations, translations)
    candidates = [np.zeros(3)]
    for rotation, translation in zip(rotations, translations):
        mat = np.eye(3) - rotation
        if abs(np.linalg.det(mat)) > 1.0e-9:
            candidates.append((np.linalg.solve(mat, -translation)) % 1.0)
    for origin in candidates:
        consistent = True
        for rotation, translation in zip(rotations, translations):
            if np.array_equal(rotation, identity):
                continue
            residual = (np.eye(3) - rotation) @ origin + translation
            if not _in_lattice(residual, lattice):
                consistent = False
                break
        if consistent:
            return True
    return False


def _normalize_hm(symbol: str) -> str:
    return " ".join(str(symbol).replace("_", " ").split())


def _point_group_map() -> dict[str, dict[str, object]]:
    registry = json.loads(POINT_GROUP_PATH.read_text(encoding="utf-8"))
    by_symbol: dict[str, dict[str, object]] = {}
    for entry in registry["point_groups"]:
        normalized = " ".join(str(entry["hm_symbol"]).split())
        by_symbol[normalized] = {
            "number": int(entry["number"]),
            "crystal_system": str(entry["crystal_system"]),
        }
    return by_symbol


def _generating_set(full: list) -> list:
    from group_theory_operations.seitz import SeitzOp, closure, equivalent

    identity = SeitzOp.identity()
    generators: list[SeitzOp] = []
    generated: list[SeitzOp] = [identity]
    for candidate in full:
        if any(equivalent(candidate, operation) for operation in generated):
            continue
        generators.append(candidate)
        generated = closure(generators)
        if len(generated) == len(full):
            break
    if len(generated) != len(full):
        raise ValueError("generating set does not reproduce the full operation set")
    return generators


def _seitz_ops(hall_number: int) -> list:
    from group_theory_operations.seitz import SeitzOp

    operations = get_symmetry_from_database(hall_number)
    rotations = np.asarray(operations["rotations"], dtype=np.int64)
    translations = np.asarray(operations["translations"], dtype=np.float64)
    return [
        SeitzOp(rotation=rotation, translation=translation)
        for rotation, translation in zip(rotations, translations)
    ]


def _centering(hall_symbol: str) -> str:
    centering = next(
        (token for token in str(hall_symbol).split()[0] if token in "PABC FIR"),
        "",
    )
    if centering not in {"P", "A", "B", "C", "F", "I", "R"}:
        raise ValueError(
            f"unexpected centering {centering!r} in hall symbol {hall_symbol!r}"
        )
    return centering


def main() -> int:
    point_groups = _point_group_map()
    by_ita: dict[int, list] = {}
    for hall_number in range(1, 531):
        info = get_spacegroup_type(hall_number)
        by_ita.setdefault(int(info.number), []).append(
            (hall_number, str(info.hall_symbol), str(info.choice))
        )
    if set(by_ita) != set(range(1, 231)):
        raise ValueError(
            "missing space group numbers: "
            f"{sorted(set(range(1, 231)) - set(by_ita))}"
        )

    entries = []
    counts = {system: 0 for system in CRYSTAL_SYSTEMS}
    symmorphic_count = 0
    for ita_number in range(1, 231):
        hall_rows = sorted(by_ita[ita_number])
        first = get_spacegroup_type(hall_rows[0][0])
        point_group = point_groups[_normalize_hm(first.pointgroup_international)]
        crystal_system = str(point_group["crystal_system"])
        counts[crystal_system] += 1

        primary_hall = None
        for hall_number, _, _ in hall_rows:
            if _round_trips(hall_number, ita_number, crystal_system):
                primary_hall = hall_number
                break
        if primary_hall is None:
            raise ValueError(f"no Hall setting round-trips to space group {ita_number}")
        primary = get_spacegroup_type(primary_hall)
        point_group = point_groups[_normalize_hm(primary.pointgroup_international)]

        settings = []
        setting_symmorphic = []
        for hall_number, hall_symbol, choice in hall_rows:
            operations = _seitz_ops(hall_number)
            rotations = np.asarray([op.rotation for op in operations], dtype=np.int64)
            translations = np.asarray([op.translation for op in operations])
            setting_symmorphic.append(_symmorphic(rotations, translations))
            settings.append(
                {
                    "hall_number": hall_number,
                    "hall_symbol": hall_symbol,
                    "choice": choice,
                    "centering": _centering(hall_symbol),
                    "symmorphic": setting_symmorphic[-1],
                    "operation_count": len(operations),
                    "generators": [generator.to_dict() for generator in _generating_set(operations)],
                }
            )
        symmorphic = any(setting_symmorphic)
        if symmorphic:
            symmorphic_count += 1

        entries.append(
            {
                "ita_number": ita_number,
                "international_short": str(primary.international_short),
                "international_full": str(primary.international_full),
                "schoenflies": str(primary.schoenflies),
                "point_group_number": int(point_group["number"]),
                "point_group_hm": _normalize_hm(primary.pointgroup_international),
                "point_group_schoenflies": str(primary.pointgroup_schoenflies),
                "crystal_system": crystal_system,
                "centering": _centering(primary.hall_symbol),
                "symmorphic": symmorphic,
                "primary_hall_number": primary_hall,
                "hall_settings": settings,
            }
        )

    if len(entries) != 230:
        raise ValueError(f"expected 230 space groups, got {len(entries)}")
    for system, expected in EXPECTED_COUNTS.items():
        if counts[system] != expected:
            raise ValueError(
                f"crystal system {system}: expected {expected}, got {counts[system]}"
            )
    if symmorphic_count != EXPECTED_SYMMORPHIC:
        raise ValueError(
            f"expected {EXPECTED_SYMMORPHIC} symmorphic space groups, got {symmorphic_count}"
        )
    point_group_numbers = {entry["point_group_number"] for entry in entries}
    if len(point_group_numbers) != 32:
        raise ValueError(f"expected all 32 point groups, got {len(point_group_numbers)}")

    document = {
        "schema_version": 1,
        "source": {
            "name": "spglib",
            "version": __import__("spglib").__version__,
            "license": "BSD-3-Clause",
            "note": "international and Schoenflies labels, Hall symbols, and "
            "Seitz operations generated from the spglib database; every Hall "
            "setting is stored, with the round-trip-verified setting marked "
            "primary.",
        },
        "generated_by": "scripts/generate_crystallographic_space_groups.py",
        "point_group_registry_sha256": hashlib.sha256(
            POINT_GROUP_PATH.read_bytes()
        ).hexdigest(),
        "conventions": {
            "seitz": "x' = R x + t in fractional coordinates of the "
            "conventional cell, translations modulo the lattice",
            "primary": "lowest Hall number whose operations round-trip: "
            "three generic positions on a conventional lattice identified "
            "by spglib as the same ITA number",
            "generators": "deterministic greedy reduction whose closure "
            "reproduces the full operation set of each Hall setting",
            "symmorphic": "fixed-point criterion modulo the cell lattice: "
            "an origin shift makes every non-identity rotation pure",
        },
        "space_groups": entries,
    }
    DATA_PATH.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {DATA_PATH}: {len(entries)} space groups, "
        f"{symmorphic_count} symmorphic, "
        f"{sum(len(e['hall_settings']) for e in entries)} hall settings"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
