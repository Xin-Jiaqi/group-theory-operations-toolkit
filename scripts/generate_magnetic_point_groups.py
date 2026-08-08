#!/usr/bin/env python3
"""Generate the canonical registry of the 122 magnetic point groups.

The standard three-part numbering and Hermann--Mauguin symbols are frozen in
the versioned table below. Spatial operations are generated from the toolkit's
validated 32 crystallographic point-group embeddings. Type-III groups are
constructed from an index-two unitary subgroup, encoded by the time-reversal
parity of the parent generators.
"""

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

from group_theory_operations import load_database, load_point_group_registry  # noqa: E402


DATA_PATH = ROOT / "data" / "magnetic_point_groups.json"


@dataclass(frozen=True, slots=True)
class BlackWhiteDefinition:
    """One standard type-III symbol and generator parity assignment."""

    hm_symbol: str
    generator_time_reversals: tuple[bool, ...]


# Order within every parent follows the standard MPOINT numbering.  The bits
# refer, in order, to the generators in crystallographic_point_groups.json.
BLACK_WHITE_DEFINITIONS: dict[int, tuple[BlackWhiteDefinition, ...]] = {
    2: (BlackWhiteDefinition("-1'", (True,)),),
    3: (BlackWhiteDefinition("2'", (True,)),),
    4: (BlackWhiteDefinition("m'", (True,)),),
    5: (
        BlackWhiteDefinition("2'/m", (True, True)),
        BlackWhiteDefinition("2/m'", (False, True)),
        BlackWhiteDefinition("2'/m'", (True, False)),
    ),
    6: (BlackWhiteDefinition("2'2'2", (False, True)),),
    7: (
        BlackWhiteDefinition("m'm2'", (True, True)),
        BlackWhiteDefinition("m'm'2", (False, True)),
    ),
    8: (
        BlackWhiteDefinition("m'mm", (True, False, True)),
        BlackWhiteDefinition("m'm'm", (False, True, False)),
        BlackWhiteDefinition("m'm'm'", (False, False, True)),
    ),
    9: (BlackWhiteDefinition("4'", (True,)),),
    10: (BlackWhiteDefinition("-4'", (True,)),),
    11: (
        BlackWhiteDefinition("4'/m", (True, False)),
        BlackWhiteDefinition("4/m'", (False, True)),
        BlackWhiteDefinition("4'/m'", (True, True)),
    ),
    12: (
        BlackWhiteDefinition("4'22'", (True, False)),
        BlackWhiteDefinition("42'2'", (False, True)),
    ),
    13: (
        BlackWhiteDefinition("4'm'm", (True, True)),
        BlackWhiteDefinition("4m'm'", (False, True)),
    ),
    14: (
        BlackWhiteDefinition("-4'2'm", (True, True)),
        BlackWhiteDefinition("-4'2m'", (True, False)),
        BlackWhiteDefinition("-42'm'", (False, True)),
    ),
    15: (
        BlackWhiteDefinition("4/m'mm", (False, True, True)),
        BlackWhiteDefinition("4'/mm'm", (True, True, False)),
        BlackWhiteDefinition("4'/m'm'm", (True, False, True)),
        BlackWhiteDefinition("4/mm'm'", (False, True, False)),
        BlackWhiteDefinition("4/m'm'm'", (False, False, True)),
    ),
    17: (BlackWhiteDefinition("-3'", (True,)),),
    18: (BlackWhiteDefinition("32'", (False, True)),),
    19: (BlackWhiteDefinition("3m'", (False, True)),),
    20: (
        BlackWhiteDefinition("-3'm", (True, False)),
        BlackWhiteDefinition("-3'm'", (True, True)),
        BlackWhiteDefinition("-3m'", (False, True)),
    ),
    21: (BlackWhiteDefinition("6'", (True,)),),
    22: (BlackWhiteDefinition("-6'", (True,)),),
    23: (
        BlackWhiteDefinition("6'/m", (True, True)),
        BlackWhiteDefinition("6/m'", (False, True)),
        BlackWhiteDefinition("6'/m'", (True, False)),
    ),
    24: (
        BlackWhiteDefinition("6'22'", (True, False)),
        BlackWhiteDefinition("62'2'", (False, True)),
    ),
    25: (
        BlackWhiteDefinition("6'mm'", (True, False)),
        BlackWhiteDefinition("6m'm'", (False, True)),
    ),
    26: (
        BlackWhiteDefinition("-6'm'2", (True, True)),
        BlackWhiteDefinition("-6'm2'", (True, False)),
        BlackWhiteDefinition("-6m'2'", (False, True)),
    ),
    27: (
        BlackWhiteDefinition("6/m'mm", (False, True, True)),
        BlackWhiteDefinition("6'/mmm'", (True, True, True)),
        BlackWhiteDefinition("6'/m'mm'", (True, False, False)),
        BlackWhiteDefinition("6/mm'm'", (False, True, False)),
        BlackWhiteDefinition("6/m'm'm'", (False, False, True)),
    ),
    29: (BlackWhiteDefinition("m'-3'", (False, False, True)),),
    30: (BlackWhiteDefinition("4'32'", (True, False)),),
    31: (BlackWhiteDefinition("-4'3m'", (True, False)),),
    32: (
        BlackWhiteDefinition("m'-3'm", (True, False, True)),
        BlackWhiteDefinition("m-3m'", (True, False, False)),
        BlackWhiteDefinition("m'-3'm'", (False, False, True)),
    ),
}


def _matrix_key(matrix: list[list[float]]) -> tuple[float, ...]:
    return tuple(round(float(value), 9) for row in matrix for value in row)


def _matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def _multiplication_table(group: dict, database: dict) -> dict[tuple[str, str], str]:
    family = database["families"][group["host_family"]]
    records = {item["name"]: item for item in family["operations"]}
    by_matrix = {
        _matrix_key(item["matrix_cartesian"]): item["name"]
        for item in family["operations"]
    }
    return {
        (left, right): by_matrix[
            _matrix_key(
                _matmul(
                    records[left]["matrix_cartesian"],
                    records[right]["matrix_cartesian"],
                )
            )
        ]
        for left in group["operations"]
        for right in group["operations"]
    }


def _parity_map(
    group: dict,
    generator_time_reversals: tuple[bool, ...],
    multiplication: dict[tuple[str, str], str],
) -> dict[str, bool]:
    generators = group["generators"]
    if len(generators) != len(generator_time_reversals):
        raise ValueError(f"wrong generator parity count for {group['hm_symbol']}")
    parity: dict[str, bool] = {"1": False}
    for name, value in zip(generators, generator_time_reversals, strict=True):
        parity[name] = value
    while len(parity) < len(group["operations"]):
        changed = False
        for left in tuple(parity):
            for right in tuple(parity):
                product = multiplication[left, right]
                value = parity[left] ^ parity[right]
                if product in parity and parity[product] != value:
                    raise ValueError(f"inconsistent time parity for {group['hm_symbol']}")
                if product not in parity:
                    parity[product] = value
                    changed = True
        if not changed:
            raise ValueError(f"generators do not span {group['hm_symbol']}")
    for left in group["operations"]:
        for right in group["operations"]:
            if parity[multiplication[left, right]] != parity[left] ^ parity[right]:
                raise ValueError(f"time parity is not a homomorphism for {group['hm_symbol']}")
    if sum(parity.values()) * 2 != len(parity):
        raise ValueError(f"unitary subgroup is not index two for {group['hm_symbol']}")
    return parity


def _operation(name: str, time_reversal: bool) -> dict[str, object]:
    return {"name": name, "time_reversal": time_reversal}


def build_payload(database: dict, point_registry: dict) -> dict:
    records: list[dict[str, object]] = []
    global_number = 0
    for parent in point_registry["point_groups"]:
        multiplication = _multiplication_table(parent, database)
        definitions: list[tuple[str, str, tuple[bool, ...] | None]] = [
            ("type_I", parent["hm_symbol"], tuple(False for _ in parent["generators"])),
            ("type_II_gray", f"{parent['hm_symbol']}1'", None),
        ]
        definitions.extend(
            ("type_III_black_white", item.hm_symbol, item.generator_time_reversals)
            for item in BLACK_WHITE_DEFINITIONS.get(parent["number"], ())
        )
        for local_number, (category, symbol, generator_bits) in enumerate(definitions, start=1):
            global_number += 1
            if category == "type_II_gray":
                operations = [
                    _operation(name, time_reversal)
                    for time_reversal in (False, True)
                    for name in parent["operations"]
                ]
                generators = [
                    *(_operation(name, False) for name in parent["generators"]),
                    _operation("1", True),
                ]
            else:
                assert generator_bits is not None
                parity = (
                    {name: False for name in parent["operations"]}
                    if category == "type_I"
                    else _parity_map(parent, generator_bits, multiplication)
                )
                operations = [_operation(name, parity[name]) for name in parent["operations"]]
                generators = [
                    _operation(name, value)
                    for name, value in zip(parent["generators"], generator_bits, strict=True)
                ]
            records.append(
                {
                    "number": global_number,
                    "magnetic_number": f"{parent['number']}.{local_number}.{global_number}",
                    "hm_symbol": symbol,
                    "category": category,
                    "parent_point_group_number": parent["number"],
                    "parent_point_group_hm": parent["hm_symbol"],
                    "crystal_system": parent["crystal_system"],
                    "host_family": parent["host_family"],
                    "order": len(operations),
                    "unitary_subgroup_order": sum(
                        not bool(item["time_reversal"]) for item in operations
                    ),
                    "generators": generators,
                    "operations": operations,
                }
            )
    if global_number != 122:
        raise ValueError(f"expected 122 magnetic point groups, generated {global_number}")
    if sum(len(items) for items in BLACK_WHITE_DEFINITIONS.values()) != 58:
        raise ValueError("expected 58 type-III magnetic point groups")
    point_bytes = (ROOT / "data" / "crystallographic_point_groups.json").read_bytes()
    operation_bytes = (ROOT / "data" / "group_operations.json").read_bytes()
    return {
        "schema_version": 1,
        "point_group_registry_sha256": hashlib.sha256(point_bytes).hexdigest(),
        "operation_catalog_sha256": hashlib.sha256(operation_bytes).hexdigest(),
        "generated_by": "scripts/generate_magnetic_point_groups.py",
        "sources": [
            {
                "title": "MPOINT Magnetic Point Group Tables",
                "url": "https://www.cryst.ehu.es/cryst/mpoint.html",
                "role": "standard magnetic point-group numbers and traditional Hermann-Mauguin symbols",
            },
            {
                "title": "Teaching crystallographic and magnetic point group symmetry using three-dimensional rendered visualizations",
                "url": "https://www.iucr.org/education/pamphlets/23/full-text",
                "role": "32 + 32 + 58 classification and index-two subgroup construction",
            },
        ],
        "conventions": {
            "magnetic_number": "parent crystallographic point-group number.local magnetic type.global number (1-122)",
            "time_reversal": "false for a unitary spatial operation and true for the same operation composed with time reversal",
            "construction": "type I: G; type II: G + theta G; type III: H + theta(G-H), where H is an index-two subgroup",
            "setting": "inherits the toolkit canonical Cartesian embedding of the parent crystallographic point group",
        },
        "magnetic_point_groups": records,
    }


def main() -> int:
    payload = build_payload(load_database(), load_point_group_registry())
    DATA_PATH.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    counts: dict[str, int] = {}
    for group in payload["magnetic_point_groups"]:
        category = str(group["category"])
        counts[category] = counts.get(category, 0) + 1
    print(f"generated 122 magnetic point groups: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
