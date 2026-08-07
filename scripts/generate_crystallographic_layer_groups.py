#!/usr/bin/env python3
"""Generate the registry of all 80 crystallographic layer groups.

Metadata comes from spglib v2.5.0 ``database/layer_spg.csv``.  That source
file is kept verbatim under ``scripts/sources``; its copyright notice and
BSD-3-Clause license are retained in ``NOTICE.spglib``. Full operations are
read from the same spglib layer
Hall database through its private extension entry point because spglib does
not yet expose a public operation-retrieval wrapper for negative (layer)
Hall numbers.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import spglib
from spglib import _spglib

from group_theory_operations.seitz import SeitzOp, closure, equivalent


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "scripts" / "sources" / "spglib-layer-groups-v2.5.0.csv"
POINT_GROUP_PATH = ROOT / "data" / "crystallographic_point_groups.json"
POINT_OPERATION_PATH = ROOT / "data" / "group_operations.json"
DATA_PATH = ROOT / "data" / "crystallographic_layer_groups.json"

SOURCE_VERSION = "2.5.0"
EXPECTED_SETTINGS = 116


def _crystal_system(number: int) -> str:
    if number <= 2:
        return "triclinic"
    if number <= 18:
        return "monoclinic"
    if number <= 48:
        return "orthorhombic"
    if number <= 64:
        return "tetragonal"
    if number <= 72:
        return "trigonal"
    return "hexagonal"


def _read_source_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with SOURCE_PATH.open(encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle):
            if len(values) != 9:
                raise ValueError(f"layer metadata row needs 9 fields: {values!r}")
            rows.append(
                {
                    "layer_hall_number": int(values[0]),
                    "setting_number": int(values[1]),
                    "choice": values[2],
                    "layer_group_number": int(values[4]),
                    "schoenflies": values[5],
                    "hall_symbol": values[6],
                    "international": values[7],
                    "international_full": values[8],
                    "international_short": values[7].split("=")[0].replace(" ", ""),
                }
            )
    if len(rows) != EXPECTED_SETTINGS:
        raise ValueError(f"expected {EXPECTED_SETTINGS} layer Hall settings")
    if [row["layer_hall_number"] for row in rows] != list(
        range(1, EXPECTED_SETTINGS + 1)
    ):
        raise ValueError("layer Hall settings must be ordered 1-116")
    if {row["layer_group_number"] for row in rows} != set(range(1, 81)):
        raise ValueError("layer metadata must cover LG1-LG80")
    return rows


def _source_operations(layer_hall_number: int) -> list[SeitzOp]:
    # The extension requires the same 192-row buffers used by spglib's public
    # space-group wrapper, even though layer groups contain at most 24 entries.
    rotations = np.zeros((192, 3, 3), dtype=np.intc)
    translations = np.zeros((192, 3), dtype=np.float64)
    count = _spglib.symmetry_from_database(
        rotations, translations, -layer_hall_number
    )
    if type(count) is not int or count < 1:
        raise ValueError(f"spglib has no layer Hall setting {layer_hall_number}")
    operations = [
        SeitzOp(rotations[index], translations[index]) for index in range(count)
    ]
    for operation in operations:
        if np.any(operation.rotation[:2, 2]) or np.any(operation.rotation[2, :2]):
            raise ValueError("layer operation mixes periodic and aperiodic axes")
        if abs(float(operation.translation[2])) > 1.0e-12:
            raise ValueError("layer operation translates along the aperiodic axis")
    return operations


def _generating_set(operations: list[SeitzOp]) -> list[SeitzOp]:
    generators: list[SeitzOp] = []
    generated = [SeitzOp.identity()]
    for candidate in operations:
        if any(equivalent(candidate, member) for member in generated):
            continue
        generators.append(candidate)
        generated = closure(generators)
        if len(generated) == len(operations):
            break
    if len(generated) != len(operations) or any(
        not any(equivalent(operation, member) for member in generated)
        for operation in operations
    ):
        raise ValueError("layer generators do not reproduce the source operations")
    return generators


def _point_group_by_schoenflies() -> dict[str, dict[str, Any]]:
    data = json.loads(POINT_GROUP_PATH.read_text(encoding="utf-8"))
    return {
        entry["schoenflies_symbol"]: entry for entry in data["point_groups"]
    }


def _legacy_layer_map() -> dict[int, dict[str, Any]]:
    data = json.loads(POINT_OPERATION_PATH.read_text(encoding="utf-8"))
    result: dict[int, dict[str, Any]] = {}
    for family in data["families"].values():
        for entry in family["layer_groups"]:
            result[int(entry["LG"])] = entry
    if set(result) != set(range(1, 81)):
        raise ValueError("point-operation catalog must cover LG1-LG80")
    return result


def main() -> int:
    rows = _read_source_rows()
    point_groups = _point_group_by_schoenflies()
    legacy = _legacy_layer_map()
    by_group: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        by_group.setdefault(row["layer_group_number"], []).append(row)

    entries = []
    for number in range(1, 81):
        source_rows = by_group[number]
        standard = next(
            (row for row in source_rows if row["setting_number"] == 1), None
        )
        if standard is None:
            raise ValueError(f"LG{number} has no standard setting")
        point_group_symbol = str(standard["schoenflies"]).split("^")[0]
        point_group = point_groups[point_group_symbol]
        if legacy[number]["point_group_base"] != point_group_symbol:
            raise ValueError(f"LG{number} point-group mapping disagrees")

        settings = []
        for row in source_rows:
            operations = _source_operations(row["layer_hall_number"])
            generators = _generating_set(operations)
            settings.append(
                {
                    "layer_hall_number": row["layer_hall_number"],
                    "setting_number": row["setting_number"],
                    "choice": row["choice"],
                    "standard": row["setting_number"] == 1,
                    "hall_symbol": row["hall_symbol"],
                    "international": row["international"],
                    "international_full": row["international_full"],
                    "international_short": row["international_short"],
                    "centering": row["international_short"][0].upper(),
                    "operation_count": len(operations),
                    "generators": [generator.to_dict() for generator in generators],
                }
            )

        entry = {
            "number": number,
            "international_short": standard["international_short"],
            "international_full": standard["international_full"],
            "schoenflies": standard["schoenflies"],
            "point_group_number": point_group["number"],
            "point_group_hm": point_group["hm_symbol"],
            "point_group_schoenflies": point_group_symbol,
            "crystal_system": _crystal_system(number),
            "centering": standard["international_short"][0].upper(),
            "primary_layer_hall_number": standard["layer_hall_number"],
            "hall_settings": settings,
        }
        if "point_group_embedding" in legacy[number]:
            entry["point_group_embedding"] = legacy[number][
                "point_group_embedding"
            ]
        entries.append(entry)

    document = {
        "schema_version": 1,
        "source": {
            "name": "spglib layer-group database",
            "metadata_version": SOURCE_VERSION,
            "operations_version": spglib.__version__,
            "license": "BSD-3-Clause",
            "url": f"https://github.com/spglib/spglib/blob/v{SOURCE_VERSION}/database/layer_spg.csv",
        },
        "generated_by": "scripts/generate_crystallographic_layer_groups.py",
        "source_metadata_sha256": hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
        "point_operation_catalog_sha256": hashlib.sha256(
            POINT_OPERATION_PATH.read_bytes()
        ).hexdigest(),
        "conventions": {
            "periodic_axes": ["x", "y"],
            "aperiodic_axis": "z",
            "seitz": "x' = R x + t in fractional coordinates; translations are modulo the periodic x/y lattice only",
            "settings": "all 116 spglib layer Hall settings are retained; setting_number 1 is the standard setting",
            "generators": "deterministic greedy reduction whose closure reproduces every source operation",
        },
        "layer_groups": entries,
    }
    DATA_PATH.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {DATA_PATH}: {len(entries)} layer groups, "
        f"{sum(len(entry['hall_settings']) for entry in entries)} Hall settings"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
