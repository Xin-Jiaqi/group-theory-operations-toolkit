"""Verification of the complete 80-layer-group registry."""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
from io import StringIO
import json
from pathlib import Path
import unittest

import numpy as np
from jsonschema import Draft202012Validator

from group_theory_operations.cli import main
from group_theory_operations.layer_groups import (
    get_crystallographic_layer_group,
    iter_crystallographic_layer_groups,
    load_layer_group_registry,
)
from group_theory_operations.seitz import SeitzOp, closure, equivalent, inverse, multiply

try:
    from spglib import _spglib
except ImportError:  # pragma: no cover
    _spglib = None


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "crystallographic_layer_groups.json"
SCHEMA_PATH = ROOT / "schema" / "crystallographic-layer-groups-v1.schema.json"
SOURCE_PATH = ROOT / "scripts" / "sources" / "spglib-layer-groups-v2.5.0.csv"
POINT_OPERATION_PATH = ROOT / "data" / "group_operations.json"

EXPECTED_COUNTS = {
    "triclinic": 2,
    "monoclinic": 16,
    "orthorhombic": 30,
    "tetragonal": 16,
    "trigonal": 8,
    "hexagonal": 8,
}


def _source_operations(layer_hall_number: int) -> list[SeitzOp]:
    if _spglib is None:  # pragma: no cover
        raise RuntimeError("spglib is not installed")
    rotations = np.zeros((192, 3, 3), dtype=np.intc)
    translations = np.zeros((192, 3), dtype=np.float64)
    count = _spglib.symmetry_from_database(
        rotations, translations, -layer_hall_number
    )
    return [SeitzOp(rotations[index], translations[index]) for index in range(count)]


class LayerGroupRegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = load_layer_group_registry(DATA_PATH)
        cls.records = tuple(iter_crystallographic_layer_groups(cls.registry))

    def test_registry_satisfies_schema(self) -> None:
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        errors = sorted(
            Draft202012Validator(schema).iter_errors(self.registry), key=str
        )
        self.assertEqual(errors, [])

    def test_source_hashes_are_bound_to_checked_in_inputs(self) -> None:
        self.assertEqual(
            self.registry["source_metadata_sha256"],
            hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.registry["point_operation_catalog_sha256"],
            hashlib.sha256(POINT_OPERATION_PATH.read_bytes()).hexdigest(),
        )

    def test_complete_numbering_settings_and_crystal_systems(self) -> None:
        self.assertEqual([record.number for record in self.records], list(range(1, 81)))
        settings = [setting for record in self.records for setting in record.hall_settings]
        self.assertEqual(
            [setting.layer_hall_number for setting in settings], list(range(1, 117))
        )
        self.assertEqual(
            {
                system: sum(record.crystal_system == system for record in self.records)
                for system in EXPECTED_COUNTS
            },
            EXPECTED_COUNTS,
        )
        self.assertEqual(len({record.point_group_hm for record in self.records}), 27)

    def test_known_groups_and_multiple_settings(self) -> None:
        checks = {
            1: ("p1", "1", 1, 1),
            2: ("p-1", "-1", 1, 2),
            5: ("p11a", "m", 3, 2),
            7: ("p112/a", "2/m", 3, 4),
            49: ("p4", "4", 1, 4),
            52: ("p4/n", "4/m", 2, 8),
            64: ("p4/nmm", "4/mmm", 2, 16),
            65: ("p3", "3", 1, 3),
            74: ("p-6", "-6", 1, 6),
            80: ("p6/mmm", "6/mmm", 1, 24),
        }
        for number, (symbol, point_group, setting_count, operation_count) in checks.items():
            record = get_crystallographic_layer_group(number, self.registry)
            self.assertEqual(record.international_short, symbol)
            self.assertEqual(record.point_group_hm, point_group)
            self.assertEqual(len(record.hall_settings), setting_count)
            standard = next(setting for setting in record.hall_settings if setting.standard)
            self.assertEqual(standard.operation_count, operation_count)

    def test_number_and_symbol_queries(self) -> None:
        expected = get_crystallographic_layer_group(80, self.registry)
        for identifier in ("80", "LG80", "p6/mmm", "p 6/m m m"):
            self.assertEqual(
                get_crystallographic_layer_group(identifier, self.registry), expected
            )

    def test_generators_close_and_preserve_the_layer_plane(self) -> None:
        for record in self.records:
            for setting in record.hall_settings:
                generators = [
                    SeitzOp.from_dict(generator) for generator in setting.generators
                ]
                generated = closure(generators)
                self.assertEqual(
                    len(generated),
                    setting.operation_count,
                    f"LG{record.number} layer Hall {setting.layer_hall_number}",
                )
                for operation in generated:
                    self.assertTrue(np.all(operation.rotation[:2, 2] == 0))
                    self.assertTrue(np.all(operation.rotation[2, :2] == 0))
                    self.assertEqual(float(operation.translation[2]), 0.0)
                    candidate = inverse(operation)
                    self.assertTrue(multiply(operation, candidate).is_identity())
                    self.assertTrue(multiply(candidate, operation).is_identity())


@unittest.skipUnless(_spglib is not None, "spglib not installed")
class LayerGroupSpglibCrossCheckTests(unittest.TestCase):
    def test_all_generators_reproduce_spglib_operations(self) -> None:
        for record in iter_crystallographic_layer_groups():
            for setting in record.hall_settings:
                generators = [
                    SeitzOp.from_dict(generator) for generator in setting.generators
                ]
                generated = closure(generators)
                source = _source_operations(setting.layer_hall_number)
                self.assertEqual(len(generated), len(source))
                for operation in generated:
                    self.assertTrue(
                        any(equivalent(operation, candidate) for candidate in source),
                        f"LG{record.number} layer Hall {setting.layer_hall_number}",
                    )


class LayerGroupCliTests(unittest.TestCase):
    def test_cli_queries_numeric_identifier_as_json(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["layer-groups", "80", "--json"]), 0)
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["number"], 80)
        self.assertEqual(payload["international_short"], "p6/mmm")

    def test_cli_queries_symbol(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["layer-groups", "p4/nmm"]), 0)
        self.assertIn("64 p4/nmm", output.getvalue())
