#!/usr/bin/env python3
"""Generate invariant optical-response bases for all crystallographic point groups."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from group_theory_operations.invariants import (  # noqa: E402
    RESPONSE_SPECS,
    response_tensor_basis,
)
from group_theory_operations import load_database  # noqa: E402
from group_theory_operations.point_groups import (  # noqa: E402
    iter_crystallographic_point_groups,
    load_point_group_registry,
)


DATA_PATH = ROOT / "data" / "optical_response_invariants.json"


def build_payload() -> dict:
    registry = load_point_group_registry()
    database = load_database()
    groups = []
    for group in iter_crystallographic_point_groups(registry):
        responses = {
            name: response_tensor_basis(
                group.number,
                name,
                database=database,
                registry=registry,
            ).to_dict()
            for name in RESPONSE_SPECS
        }
        groups.append(
            {
                "number": group.number,
                "hm_symbol": group.hm_symbol,
                "schoenflies_symbol": group.schoenflies_symbol,
                "centrosymmetric": group.centrosymmetric,
                "responses": responses,
            }
        )
    return {
        "schema_version": 1,
        "point_group_registry_sha256": hashlib.sha256(
            (ROOT / "data" / "crystallographic_point_groups.json").read_bytes()
        ).hexdigest(),
        "operation_catalog_sha256": hashlib.sha256(
            (ROOT / "data" / "group_operations.json").read_bytes()
        ).hexdigest(),
        "generated_by": "scripts/generate_optical_response_invariants.py",
        "conventions": {
            "equivariance": "A(R) T = T B(R)",
            "output": "polar Cartesian vector with A(R)=D(R)",
            "symmetric_input": "M_+(R) on (xx, yy, zz, xy, xz, yz); off-diagonal entries are full sums",
            "antisymmetric_input": "M_-(R) on h=i E x E*",
            "scope": "spatial point-group selection rules only; no microscopic, time-reversal, frequency, or unit constraints",
        },
        "point_groups": groups,
    }


def main() -> int:
    payload = build_payload()
    DATA_PATH.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
