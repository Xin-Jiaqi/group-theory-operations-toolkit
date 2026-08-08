#!/usr/bin/env python3
"""Generate the 528 magnetic-layer-group point-co-group registry."""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = (
    ROOT / "scripts" / "sources" / "magnetic-layer-groups-528-v2023.tsv"
)
OPERATION_PATH = ROOT / "data" / "group_operations.json"
LAYER_PATH = ROOT / "data" / "crystallographic_layer_groups.json"
DATA_PATH = ROOT / "data" / "magnetic_layer_groups.json"

EXPECTED_TYPES = {"I": 80, "II": 80, "III": 246, "IV": 122}
TRANSFORMS = {
    "I": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "P1": [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
    "P2": [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
    "P3": [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    "P4": [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    "P5": [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_rows() -> list[dict[str, str]]:
    with SOURCE_PATH.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda item: int(item["global_number"]))
    if [int(item["global_number"]) for item in rows] != list(range(1, 529)):
        raise ValueError("source must contain ordered magnetic groups 1-528")
    if Counter(item["magnetic_type"] for item in rows) != EXPECTED_TYPES:
        raise ValueError("magnetic group type counts must be 80 + 80 + 246 + 122")
    return rows


def _operation_catalog() -> tuple[
    dict[int, str], dict[tuple[str, str], dict[str, Any]]
]:
    data = json.loads(OPERATION_PATH.read_text(encoding="utf-8"))
    layer_family: dict[int, str] = {}
    operations: dict[tuple[str, str], dict[str, Any]] = {}
    for family_name, family in data["families"].items():
        for item in family["operations"]:
            operations[(family_name, item["name"])] = item
        for group in family["layer_groups"]:
            layer_family[group["LG"]] = family_name
    return layer_family, operations


def _anti_translation(value: str) -> list[int | float] | None:
    if not value:
        return None
    result: list[int | float] = []
    for token in value.split(","):
        if "/" in token:
            numerator, denominator = token.split("/", maxsplit=1)
            result.append(int(numerator) / int(denominator))
        else:
            result.append(int(token))
    if len(result) != 3:
        raise ValueError("anti-translation must have three fractional coordinates")
    return result


def _point_operations(
    row: dict[str, str],
    family: str,
    catalog: dict[tuple[str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    names = row["point_operation_names"].split(",")
    primed = {
        int(value)
        for value in row["primed_operation_indices"].split(",")
        if value
    }
    values: list[tuple[str, bool]]
    if row["magnetic_type"] in {"II", "IV"}:
        values = [(name, False) for name in names] + [(name, True) for name in names]
    else:
        values = [
            (name, index in primed) for index, name in enumerate(names, start=1)
        ]
    operations: list[dict[str, Any]] = []
    for name, time_reversal in values:
        source = catalog[(family, name)]
        operations.append(
            {
                "name": name,
                "label": f"{name}'" if time_reversal else name,
                "time_reversal": time_reversal,
                "matrix_fractional": source["matrix_fractional"],
                "matrix_cartesian": source["matrix_cartesian"],
            }
        )
    return operations


def _validate_closure(operations: list[dict[str, Any]]) -> None:
    lookup = {
        (
            tuple(np.asarray(item["matrix_fractional"], dtype=np.int64).ravel()),
            item["time_reversal"],
        )
        for item in operations
    }
    if len(lookup) != len(operations):
        raise ValueError("magnetic point operations must be unique")
    for left in operations:
        left_matrix = np.asarray(left["matrix_fractional"], dtype=np.int64)
        for right in operations:
            right_matrix = np.asarray(right["matrix_fractional"], dtype=np.int64)
            product = (
                tuple((left_matrix @ right_matrix).ravel()),
                bool(left["time_reversal"]) ^ bool(right["time_reversal"]),
            )
            if product not in lookup:
                raise ValueError("magnetic point operations are not closed")


def generate() -> dict[str, Any]:
    rows = _source_rows()
    layer_family, catalog = _operation_catalog()
    layer_groups = {
        item["number"]: item
        for item in json.loads(LAYER_PATH.read_text(encoding="utf-8"))["layer_groups"]
    }
    records: list[dict[str, Any]] = []
    for row in rows:
        parent_number = int(row["parent_layer_group_number"])
        family = layer_family[parent_number]
        operations = _point_operations(row, family, catalog)
        _validate_closure(operations)
        anti_translation = _anti_translation(row["anti_translation"])
        if (row["magnetic_type"] == "IV") != (anti_translation is not None):
            raise ValueError("only type-IV groups may carry an anti-translation")
        records.append(
            {
                "global_number": int(row["global_number"]),
                "og_number": row["og_number"],
                "family_number": int(row["family_number"]),
                "magnetic_type": row["magnetic_type"],
                "litvin_og_symbol_ascii": row["litvin_og_symbol_ascii"],
                "magnetic_point_group_symbol": row["magnetic_point_group_symbol"],
                "parent_layer_group_number": parent_number,
                "parent_layer_group_symbol": layer_groups[parent_number][
                    "international_short"
                ],
                "crystal_system": layer_groups[parent_number]["crystal_system"],
                "host_family": family,
                "point_operation_count": len(operations),
                "unitary_subgroup_order": sum(
                    not item["time_reversal"] for item in operations
                ),
                "anti_translation_fractional": anti_translation,
                "corresponding_magnetic_space_group": {
                    "bns_number": row["corresponding_msg_bns_number"],
                    "uni_number": int(row["corresponding_msg_uni_number"]),
                    "og_number": row["corresponding_msg_og_number"],
                    "basis_transform_key": row["basis_transform"],
                    "basis_transform": TRANSFORMS[row["basis_transform"]],
                },
                "source_pdf_page": int(row["source_pdf_page"]),
                "point_operations": operations,
            }
        )
    return {
        "schema_version": 1,
        "registry_version": "2026.08",
        "scope": (
            "Finite magnetic point-co-group representatives for homogeneous "
            "tensor constraints. Spatial translations are omitted; the "
            "distinguished type-IV anti-translation is retained separately."
        ),
        "source": {
            "classification_and_operation_labels": {
                "citation": (
                    "D. B. Litvin, Magnetic Group Tables, Part 3: Magnetic "
                    "Layer Groups (IUCr, 2005)"
                ),
                "url": "https://journals.iucr.org/a/issues/2005/03/00/sh5024/sh5024sup3.pdf",
                "sha256": "f4d7a37354830b5f52d0478fa3a8b7272a1a1f48791aa6cc6ebe8bc174e4c294",
            },
            "magnetic_space_group_correspondence": {
                "citation": (
                    "Z.-M. Zhang et al., Encyclopedia of emergent particles "
                    "in 528 magnetic layer groups and 394 magnetic rod groups, "
                    "Phys. Rev. B 107, 075405 (2023), Supplementary Table S1/S2"
                ),
                "url": "https://arxiv.org/abs/2210.11080",
                "supplement_sha256": "db054fc9467e919cc03dcd72937f16ddba95769bf1185251aa94b63e9c49f8f8",
            },
        },
        "source_table_sha256": _sha256(SOURCE_PATH),
        "point_operation_catalog_sha256": _sha256(OPERATION_PATH),
        "layer_group_registry_sha256": _sha256(LAYER_PATH),
        "magnetic_layer_groups": records,
    }


def main() -> int:
    document = generate()
    DATA_PATH.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {len(document['magnetic_layer_groups'])} groups to {DATA_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
