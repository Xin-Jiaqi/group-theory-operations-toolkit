#!/usr/bin/env python3
"""Generate the Hall-setting Wyckoff-position registry from spglib v2.5.0.

The source file is spglib's BSD-3-Clause ``database/Wyckoff.csv``.  Pass its
path explicitly so registry regeneration never depends on an unpinned package
installation:

    python3 scripts/generate_wyckoff_positions.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
import spglib


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "crystallographic_wyckoff_positions.json"
DEFAULT_SOURCE = ROOT / "scripts" / "sources" / "spglib-Wyckoff-v2.5.0.csv"
SOURCE_VERSION = "2.5.0"
SOURCE_TAG = "v2.5.0"
SOURCE_COMMIT = "e4531bb49371dce3e807c2095a4d9d9b7245c524"
EXPECTED_SOURCE_SHA256 = (
    "d3d786a1f0187e5c6d69a3ade35648ffab34fd1b977d61ad84d8b0434b8b7ca0"
)

_CENTERING_TRANSLATIONS = {
    "P": ((0, 0, 0),),
    "A": ((0, 0, 0), (0, 12, 12)),
    "B": ((0, 0, 0), (12, 0, 12)),
    "C": ((0, 0, 0), (12, 12, 0)),
    "I": ((0, 0, 0), (12, 12, 12)),
    "R": ((0, 0, 0),),
    "H": ((0, 0, 0), (16, 8, 8), (8, 16, 16)),
    "F": ((0, 0, 0), (0, 12, 12), (12, 0, 12), (12, 12, 0)),
}

_TERM = re.compile(r"([+-]?)(?:(\d+))?([xyz])")
_FRACTION = re.compile(r"([+-]?)(\d+)/(\d+)")


def _source_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_source(path: Path) -> list[dict[str, Any]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        lines = lines[: lines.index("end of data")]
    except ValueError as exc:
        raise ValueError("Wyckoff.csv is missing its end-of-data marker") from exc
    headers = [index for index, line in enumerate(lines) if line.split(":")[0].isdigit()]
    headers.append(len(lines))
    settings: list[dict[str, Any]] = []
    for start, stop in zip(headers, headers[1:]):
        header = lines[start].split(":")
        hall_number = int(header[0])
        source_symbol = header[1].strip()
        centering = source_symbol[0]
        if hall_number in {433, 436, 444, 450, 452, 458, 460}:
            centering = "H"
        positions: list[dict[str, Any]] = []
        current: dict[str, Any] | None = None
        for line in lines[start + 1 : stop]:
            fields = line.split(":")
            if len(fields) < 9:
                fields.extend([""] * (9 - len(fields)))
            if fields[2].strip().isdigit():
                current = {
                    "multiplicity": int(fields[2]),
                    "letter": fields[3].strip(),
                    "site_symmetry": fields[4].strip(),
                    "coordinate_expressions": [],
                }
                positions.append(current)
            if current is not None:
                current["coordinate_expressions"].extend(
                    value.strip()
                    for value in fields[5:9]
                    if value.strip()
                )
        settings.append(
            {
                "hall_number": hall_number,
                "source_symbol": source_symbol,
                "centering": centering,
                "positions": positions,
            }
        )
    return settings


def _expression_component(component: str) -> tuple[list[int], int]:
    coefficients = [0, 0, 0]
    for sign, magnitude, variable in _TERM.findall(component):
        value = int(magnitude) if magnitude else 1
        if sign == "-":
            value *= -1
        coefficients["xyz".index(variable)] += value
    translation = 0
    for sign, numerator, denominator in _FRACTION.findall(component):
        denominator_int = int(denominator)
        if 24 % denominator_int:
            raise ValueError(f"translation denominator {denominator_int} does not divide 24")
        value = int(numerator) * (24 // denominator_int)
        translation += -value if sign == "-" else value
    return coefficients, translation % 24


def _affine_map(expression: str) -> dict[str, Any]:
    components = expression.strip().strip("()").split(",")
    if len(components) != 3:
        raise ValueError(f"invalid Wyckoff coordinate expression: {expression}")
    matrix: list[list[int]] = []
    translation: list[int] = []
    for component in components:
        row, shift = _expression_component(component)
        matrix.append(row)
        translation.append(shift)
    return {
        "parameter_matrix": matrix,
        "translation_numerators": translation,
    }


def _parameter_dimension(maps: list[dict[str, Any]]) -> int:
    ranks = {
        int(np.linalg.matrix_rank(np.asarray(item["parameter_matrix"], dtype=float)))
        for item in maps
    }
    if len(ranks) != 1:
        raise ValueError("one Wyckoff position has inconsistent parameter ranks")
    return ranks.pop()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, nargs="?", default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    source_hash = _source_sha256(args.source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            "Wyckoff.csv does not match pinned spglib v2.5.0: " + source_hash
        )
    if spglib.__version__ != SOURCE_VERSION:
        raise ValueError(
            f"generator requires spglib {SOURCE_VERSION}, found {spglib.__version__}"
        )

    settings = _parse_source(args.source)
    if [item["hall_number"] for item in settings] != list(range(1, 531)):
        raise ValueError("source does not contain ordered Hall numbers 1-530")
    position_count = 0
    representative_count = 0
    expanded_count = 0
    for setting in settings:
        info = spglib.get_spacegroup_type(setting["hall_number"])
        setting["ita_number"] = int(info.number)
        setting["hall_symbol"] = str(info.hall_symbol)
        setting["choice"] = str(info.choice)
        setting["centering_translation_numerators"] = [
            list(translation)
            for translation in _CENTERING_TRANSLATIONS[setting["centering"]]
        ]
        position_count += len(setting["positions"])
        for position in setting["positions"]:
            maps = [
                _affine_map(expression)
                for expression in position.pop("coordinate_expressions")
            ]
            position["parameter_dimension"] = _parameter_dimension(maps)
            position["representative_maps"] = maps
            representative_count += len(maps)
            expanded_count += position["multiplicity"]

    if (position_count, representative_count, expanded_count) != (
        3467,
        15117,
        24295,
    ):
        raise ValueError("unexpected Wyckoff registry dimensions")
    payload = {
        "schema_version": 1,
        "source": {
            "name": "spglib",
            "version": SOURCE_VERSION,
            "tag": SOURCE_TAG,
            "commit": SOURCE_COMMIT,
            "license": "BSD-3-Clause",
            "file": "database/Wyckoff.csv",
            "sha256": source_hash,
            "url": f"https://github.com/spglib/spglib/blob/{SOURCE_COMMIT}/database/Wyckoff.csv",
        },
        "generated_by": "scripts/generate_wyckoff_positions.py",
        "conventions": {
            "setting_key": "spglib Hall number 1-530",
            "coordinate_action": "fractional column coordinates x' = A q + t/24",
            "parameters": "q = (x, y, z); repeated or absent symbols encode constrained coordinates",
            "centering": "representative maps are expanded by the separately stored centering translations",
            "rhombohedral": "Hall choices H and R follow spglib; H uses the obverse hexagonal cell and R the primitive rhombohedral cell",
            "origin_and_axes": "the Hall symbol and spglib choice define the origin and axis setting; no implicit conversion is applied",
        },
        "counts": {
            "hall_settings": len(settings),
            "wyckoff_positions": position_count,
            "representative_maps": representative_count,
            "expanded_coordinate_maps": expanded_count,
        },
        "hall_settings": settings,
    }
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
