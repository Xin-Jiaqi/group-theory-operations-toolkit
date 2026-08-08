#!/usr/bin/env python3
"""Generate real-structure classification fixtures from spglib v2.5.0.

The source POSCAR files remain in the upstream spglib repository.  This
generator verifies their pinned SHA-256 values before converting the selected
cells into a compact, self-contained JSON fixture for this repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import spglib


SPGLIB_VERSION = "2.5.0"
SPGLIB_COMMIT = "e4531bb49371dce3e807c2095a4d9d9b7245c524"
SYMPREC = 1.0e-5
ANGLE_TOLERANCE = -1.0

SOURCES = (
    {
        "id": "sio2_triclinic",
        "formula": "SiO2",
        "crystal_system": "triclinic",
        "path": "test/functional/python/data/triclinic/POSCAR-001",
        "sha256": "9400b4f33dc874e52a28767e71e9106dee13d46c7655d4a14725bc90163a7526",
    },
    {
        "id": "sio2_monoclinic",
        "formula": "SiO2",
        "crystal_system": "monoclinic",
        "path": "test/functional/python/data/monoclinic/POSCAR-004",
        "sha256": "a5c00727d039b5eb9f6c7b976f37205a35b26b134f554b5fa738c379a99fc35f",
    },
    {
        "id": "bates3_orthorhombic",
        "formula": "BaTeS3",
        "crystal_system": "orthorhombic",
        "path": "test/functional/python/data/orthorhombic/POSCAR-062",
        "sha256": "a208fbf236ef58f1fc0b66da62bfe5a85f6cba7b98ec7506436054e619f7b9f9",
    },
    {
        "id": "mno2_tetragonal",
        "formula": "MnO2",
        "crystal_system": "tetragonal",
        "path": "test/functional/python/data/tetragonal/POSCAR-136",
        "sha256": "969491b355cc65de1182cbb26086a7d3209ecca50a905fe475eca9e9cb9e4b3c",
    },
    {
        "id": "rucl3_trigonal",
        "formula": "RuCl3",
        "crystal_system": "trigonal",
        "path": "test/functional/python/data/trigonal/POSCAR-158",
        "sha256": "7684e9510fb0c29bc076bb5b896b9b4786e83e061179b065960755dd30d598a3",
    },
    {
        "id": "aucn_hexagonal",
        "formula": "AuCN",
        "crystal_system": "hexagonal",
        "path": "test/functional/python/data/hexagonal/POSCAR-183-2",
        "sha256": "74b77a9ee8cfda258f4b3bb78627ada7994d165497f23e9351974f2b0f8c1801",
    },
    {
        "id": "cssnbr3_cubic",
        "formula": "CsSnBr3",
        "crystal_system": "cubic",
        "path": "test/functional/python/data/cubic/POSCAR-221-2",
        "sha256": "bfacab3c23563283fc5d4aff9637823da7788bb5db2afb33de476505481dd7ff",
    },
)


def _read_poscar(path: Path) -> tuple[list[list[float]], list[str], list[list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    scale = float(lines[1])
    if scale <= 0:
        raise ValueError(f"{path}: only a positive POSCAR scale is supported")
    lattice = [
        [scale * float(value) for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]

    try:
        counts = [int(value) for value in lines[5].split()]
        coordinate_header = 6
    except ValueError:
        counts = [int(value) for value in lines[6].split()]
        coordinate_header = 7
    if lines[coordinate_header].lower().startswith("s"):
        coordinate_header += 1
    if not lines[coordinate_header].lower().startswith("d"):
        raise ValueError(f"{path}: fixture sources must use direct coordinates")

    atom_count = sum(counts)
    coordinate_lines = lines[coordinate_header + 1 : coordinate_header + 1 + atom_count]
    if len(coordinate_lines) != atom_count:
        raise ValueError(f"{path}: coordinate count does not match POSCAR counts")

    positions: list[list[float]] = []
    species: list[str] = []
    for index, line in enumerate(coordinate_lines):
        positions.append([float(value) for value in line.split()[:3]])
        match = re.search(r"#\s*([A-Z][a-z]?)\d+\s*$", line)
        if match is None:
            raise ValueError(f"{path}: atom {index + 1} has no element label")
        species.append(match.group(1))
    if [species.count(symbol) for symbol in dict.fromkeys(species)] != counts:
        raise ValueError(f"{path}: element labels do not match POSCAR counts")
    return lattice, species, positions


def _type_numbers(species: list[str]) -> list[int]:
    mapping: dict[str, int] = {}
    return [mapping.setdefault(symbol, len(mapping) + 1) for symbol in species]


def _structure_fixture(root: Path, source: dict[str, str]) -> dict[str, Any]:
    path = root / source["path"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != source["sha256"]:
        raise ValueError(f"{source['path']}: expected SHA-256 {source['sha256']}, got {digest}")
    lattice, species, positions = _read_poscar(path)
    dataset = spglib.get_symmetry_dataset(
        (lattice, positions, _type_numbers(species)),
        symprec=SYMPREC,
        angle_tolerance=ANGLE_TOLERANCE,
    )
    if dataset is None:
        raise RuntimeError(f"spglib could not classify {source['path']}")
    database_symmetry = spglib.get_symmetry_from_database(dataset.hall_number)
    return {
        "id": source["id"],
        "formula": source["formula"],
        "source": {
            "path": source["path"],
            "url": f"https://github.com/spglib/spglib/blob/{SPGLIB_COMMIT}/{source['path']}",
            "sha256": digest,
        },
        "lattice": lattice,
        "species": species,
        "fractional_coordinates": positions,
        "expected": {
            "ita_number": int(dataset.number),
            "hall_number": int(dataset.hall_number),
            "international_short": str(dataset.international),
            "hall_symbol": str(dataset.hall),
            "setting_choice": str(dataset.choice),
            "point_group_hm": str(dataset.pointgroup),
            "crystal_system": source["crystal_system"],
            "input_cell_operation_count": len(dataset.rotations),
            "hall_setting_operation_count": len(database_symmetry["rotations"]),
            "equivalent_atom_orbit_count": len(set(dataset.equivalent_atoms.tolist())),
            "wyckoff_letters": [str(value) for value in dataset.wyckoffs],
        },
    }


def generate(spglib_root: Path) -> dict[str, Any]:
    if spglib.__version__ != SPGLIB_VERSION:
        raise RuntimeError(
            f"fixture generation requires spglib {SPGLIB_VERSION}, got {spglib.__version__}"
        )
    return {
        "schema_version": 1,
        "fixture_kind": "real_crystal_structure_classification",
        "generator": "scripts/generate_spglib_structure_fixtures.py",
        "upstream": {
            "name": "spglib",
            "version": SPGLIB_VERSION,
            "git_commit": SPGLIB_COMMIT,
            "license": "BSD-3-Clause",
            "license_path": "COPYING",
        },
        "classification": {
            "reference_spglib_version": SPGLIB_VERSION,
            "symprec": SYMPREC,
            "angle_tolerance": ANGLE_TOLERANCE,
        },
        "structures": [_structure_fixture(spglib_root, source) for source in SOURCES],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spglib-root", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/spglib_real_structures_v2.5.0.json"),
    )
    args = parser.parse_args()
    data = generate(args.spglib_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
