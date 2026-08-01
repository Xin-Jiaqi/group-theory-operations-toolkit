#!/usr/bin/env python3
"""Generate the canonical registry of the 32 crystallographic point groups."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from group_theory_operations import load_database  # noqa: E402


DATA_PATH = ROOT / "data" / "crystallographic_point_groups.json"


@dataclass(frozen=True, slots=True)
class Definition:
    number: int
    hm_symbol: str
    schoenflies_symbol: str
    crystal_system: str
    host_family: str
    generators: tuple[str, ...]
    aliases: tuple[str, ...] = ()


DEFINITIONS = (
    Definition(1, "1", "C1", "triclinic", "cubic_Oh", ()),
    Definition(2, "-1", "Ci", "triclinic", "cubic_Oh", ("-1",), ("S2",)),
    Definition(3, "2", "C2", "monoclinic", "tetragonal_D4h", ("2_001",)),
    Definition(4, "m", "Cs", "monoclinic", "tetragonal_D4h", ("m_001",)),
    Definition(5, "2/m", "C2h", "monoclinic", "tetragonal_D4h", ("2_001", "-1")),
    Definition(6, "222", "D2", "orthorhombic", "tetragonal_D4h", ("2_001", "2_100")),
    Definition(7, "mm2", "C2v", "orthorhombic", "tetragonal_D4h", ("2_001", "m_100")),
    Definition(8, "mmm", "D2h", "orthorhombic", "tetragonal_D4h", ("2_001", "2_100", "-1")),
    Definition(9, "4", "C4", "tetragonal", "tetragonal_D4h", ("4+_001",)),
    Definition(10, "-4", "S4", "tetragonal", "tetragonal_D4h", ("-4+_001",)),
    Definition(11, "4/m", "C4h", "tetragonal", "tetragonal_D4h", ("4+_001", "-1")),
    Definition(12, "422", "D4", "tetragonal", "tetragonal_D4h", ("4+_001", "2_100")),
    Definition(13, "4mm", "C4v", "tetragonal", "tetragonal_D4h", ("4+_001", "m_100")),
    Definition(14, "-42m", "D2d", "tetragonal", "tetragonal_D4h", ("-4+_001", "2_100")),
    Definition(15, "4/mmm", "D4h", "tetragonal", "tetragonal_D4h", ("4+_001", "2_100", "-1")),
    Definition(16, "3", "C3", "trigonal", "hexagonal_D6h", ("3+_001",)),
    Definition(17, "-3", "C3i", "trigonal", "hexagonal_D6h", ("-3+_001",), ("S6",)),
    Definition(18, "32", "D3", "trigonal", "hexagonal_D6h", ("3+_001", "2_100")),
    Definition(19, "3m", "C3v", "trigonal", "hexagonal_D6h", ("3+_001", "m_100")),
    Definition(20, "-3m", "D3d", "trigonal", "hexagonal_D6h", ("-3+_001", "m_100")),
    Definition(21, "6", "C6", "hexagonal", "hexagonal_D6h", ("6+_001",)),
    Definition(22, "-6", "C3h", "hexagonal", "hexagonal_D6h", ("-6+_001",)),
    Definition(23, "6/m", "C6h", "hexagonal", "hexagonal_D6h", ("6+_001", "-1")),
    Definition(24, "622", "D6", "hexagonal", "hexagonal_D6h", ("6+_001", "2_100")),
    Definition(25, "6mm", "C6v", "hexagonal", "hexagonal_D6h", ("6+_001", "m_100")),
    Definition(26, "-6m2", "D3h", "hexagonal", "hexagonal_D6h", ("-6+_001", "m_100")),
    Definition(27, "6/mmm", "D6h", "hexagonal", "hexagonal_D6h", ("6+_001", "2_100", "-1")),
    Definition(28, "23", "T", "cubic", "cubic_Oh", ("3+_111", "2_001")),
    Definition(29, "m-3", "Th", "cubic", "cubic_Oh", ("3+_111", "2_001", "-1")),
    Definition(30, "432", "O", "cubic", "cubic_Oh", ("4+_001", "3+_111")),
    Definition(31, "-43m", "Td", "cubic", "cubic_Oh", ("-4+_001", "3+_111")),
    Definition(32, "m-3m", "Oh", "cubic", "cubic_Oh", ("4+_001", "3+_111", "-1")),
)


POLAR_GROUPS = {"1", "2", "m", "mm2", "4", "4mm", "3", "3m", "6", "6mm"}
CHIRAL_GROUPS = {"1", "2", "222", "4", "422", "3", "32", "6", "622", "23", "432"}


def _matmul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def _matrix_key(matrix, digits: int = 9):
    return tuple(round(float(value), digits) for row in matrix for value in row)


def _closure(family: dict, generators: tuple[str, ...]) -> list[str]:
    operations = {operation["name"]: operation for operation in family["operations"]}
    by_matrix = {
        _matrix_key(operation["matrix_cartesian"]): operation["name"]
        for operation in family["operations"]
    }
    if "1" not in operations or any(name not in operations for name in generators):
        raise ValueError(f"invalid generators {generators!r}")
    members = {"1", *generators}
    while True:
        expanded = set(members)
        for left_name in members:
            for right_name in members:
                product = _matrix_key(
                    _matmul(
                        operations[left_name]["matrix_cartesian"],
                        operations[right_name]["matrix_cartesian"],
                    )
                )
                try:
                    expanded.add(by_matrix[product])
                except KeyError as exc:
                    raise ValueError(
                        f"{left_name} * {right_name} is absent from host family"
                    ) from exc
        if expanded == members:
            break
        members = expanded
    index = {operation["name"]: operation["index"] for operation in family["operations"]}
    return sorted(members, key=index.__getitem__)


def build_payload(database: dict) -> dict:
    groups = []
    for definition in DEFINITIONS:
        family = database["families"][definition.host_family]
        operations = _closure(family, definition.generators)
        groups.append(
            {
                "number": definition.number,
                "hm_symbol": definition.hm_symbol,
                "schoenflies_symbol": definition.schoenflies_symbol,
                "aliases": list(definition.aliases),
                "crystal_system": definition.crystal_system,
                "host_family": definition.host_family,
                "order": len(operations),
                "generators": list(definition.generators),
                "operations": operations,
                "centrosymmetric": "-1" in operations,
                "polar": definition.hm_symbol in POLAR_GROUPS,
                "chiral": definition.hm_symbol in CHIRAL_GROUPS,
            }
        )
    source_bytes = (ROOT / "data" / "group_operations.json").read_bytes()
    return {
        "schema_version": 1,
        "source_catalog_schema_version": database["schema_version"],
        "source_catalog_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "generated_by": "scripts/generate_crystallographic_point_groups.py",
        "conventions": {
            "point_group_order": "standard crystallographic point-group sequence 1-32",
            "setting": "toolkit canonical orthonormal Cartesian embedding; principal rotation axis along z; explicit operations are authoritative",
            "operation_source": "operation names resolve through data/group_operations.json",
        },
        "point_groups": groups,
    }


def main() -> int:
    payload = build_payload(load_database())
    DATA_PATH.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    for group in payload["point_groups"]:
        print(group["number"], group["hm_symbol"], group["order"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
