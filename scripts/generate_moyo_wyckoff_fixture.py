#!/usr/bin/env python3
"""Regenerate the pinned independent moyopy Wyckoff cross-check."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import moyopy


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "tests" / "fixtures" / "spglib_real_structures_v2.5.0.json"
DEFAULT_OUTPUT = ROOT / "tests" / "fixtures" / "moyo_wyckoff_crosscheck_v0.10.0.json"
EXPECTED_MOYO_VERSION = "0.10.0"
MOYO_TAG_COMMIT = "43ab761cdb2ff1936c569a2c348456caf49ebc73"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    if moyopy.__version__ != EXPECTED_MOYO_VERSION:
        parser.error(
            f"requires moyopy {EXPECTED_MOYO_VERSION}, found {moyopy.__version__}"
        )

    source_bytes = input_path.read_bytes()
    source = json.loads(source_bytes)
    records = []
    for item in source["structures"]:
        unique_species = list(dict.fromkeys(item["species"]))
        numbers = [unique_species.index(symbol) + 1 for symbol in item["species"]]
        expected = item["expected"]
        dataset = moyopy.MoyoDataset(
            moyopy.Cell(
                basis=item["lattice"],
                positions=item["fractional_coordinates"],
                numbers=numbers,
            ),
            symprec=1.0e-5,
            setting=moyopy.Setting.hall_number(expected["hall_number"]),
        )
        records.append(
            {
                "id": item["id"],
                "ita_number": dataset.number,
                "hall_number": dataset.hall_number,
                "wyckoff_letters": dataset.wyckoffs,
                "site_symmetry_symbols": dataset.site_symmetry_symbols,
                "orbit_labels": dataset.orbits,
            }
        )
    payload = {
        "schema_version": 1,
        "fixture_kind": "independent_wyckoff_crosscheck",
        "generator": "scripts/generate_moyo_wyckoff_fixture.py",
        "upstream": {
            "name": "moyopy",
            "version": EXPECTED_MOYO_VERSION,
            "git_tag_commit": MOYO_TAG_COMMIT,
            "license": "MIT OR Apache-2.0",
            "url": "https://github.com/spglib/moyo",
        },
        "input_fixture": {
            "path": str(input_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
        },
        "classification": {
            "symprec": 1.0e-5,
            "setting": "explicit fixture Hall number",
        },
        "structures": records,
    }
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
