"""Real-structure classification checks against the space-group registry."""

from __future__ import annotations

import json
import math
from pathlib import Path
import unittest

import numpy as np

from group_theory_operations.space_groups import get_crystallographic_space_group

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "spglib_real_structures_v2.5.0.json"
SPGLIB_COMMIT = "e4531bb49371dce3e807c2095a4d9d9b7245c524"
CRYSTAL_SYSTEMS = {
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
}


def _load_fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _type_numbers(species: list[str]) -> list[int]:
    mapping: dict[str, int] = {}
    return [mapping.setdefault(symbol, len(mapping) + 1) for symbol in species]


class RealStructureFixtureContractTests(unittest.TestCase):
    def test_fixture_provenance_shape_and_crystal_system_coverage(self) -> None:
        fixture = _load_fixture()
        self.assertEqual(fixture["schema_version"], 1)
        self.assertEqual(fixture["generator"], "scripts/generate_spglib_structure_fixtures.py")
        self.assertEqual(fixture["upstream"]["version"], "2.5.0")
        self.assertEqual(fixture["upstream"]["git_commit"], SPGLIB_COMMIT)
        self.assertEqual(fixture["upstream"]["license"], "BSD-3-Clause")
        self.assertEqual(len(fixture["structures"]), 7)
        self.assertEqual(
            {item["expected"]["crystal_system"] for item in fixture["structures"]},
            CRYSTAL_SYSTEMS,
        )
        self.assertEqual(len({item["id"] for item in fixture["structures"]}), 7)

        for item in fixture["structures"]:
            with self.subTest(item=item["id"]):
                self.assertRegex(item["source"]["sha256"], r"^[0-9a-f]{64}$")
                self.assertIn(SPGLIB_COMMIT, item["source"]["url"])
                self.assertEqual(len(item["lattice"]), 3)
                self.assertTrue(all(len(row) == 3 for row in item["lattice"]))
                self.assertEqual(len(item["species"]), len(item["fractional_coordinates"]))
                self.assertTrue(item["species"])
                self.assertTrue(
                    all(
                        len(position) == 3 and all(math.isfinite(value) for value in position)
                        for position in item["fractional_coordinates"]
                    )
                )
                self.assertGreater(abs(float(np.linalg.det(item["lattice"]))), 1.0e-9)
                self.assertEqual(
                    len(item["expected"]["wyckoff_letters"]), len(item["species"])
                )


@unittest.skipUnless(spglib is not None, "spglib not installed")
class RealStructureClassificationTests(unittest.TestCase):
    def test_spglib_classification_matches_fixture_and_registry(self) -> None:
        fixture = _load_fixture()
        settings = fixture["classification"]
        for item in fixture["structures"]:
            with self.subTest(item=item["id"], formula=item["formula"]):
                expected = item["expected"]
                dataset = spglib.get_symmetry_dataset(
                    (
                        item["lattice"],
                        item["fractional_coordinates"],
                        _type_numbers(item["species"]),
                    ),
                    symprec=settings["symprec"],
                    angle_tolerance=settings["angle_tolerance"],
                )
                self.assertIsNotNone(dataset)
                assert dataset is not None
                self.assertEqual(dataset.number, expected["ita_number"])
                self.assertEqual(dataset.hall_number, expected["hall_number"])
                self.assertEqual(dataset.international, expected["international_short"])
                self.assertEqual(dataset.hall, expected["hall_symbol"])
                self.assertEqual(dataset.choice, expected["setting_choice"])
                self.assertEqual(dataset.pointgroup, expected["point_group_hm"])
                self.assertEqual(
                    len(dataset.rotations), expected["input_cell_operation_count"]
                )
                self.assertEqual(
                    len(set(dataset.equivalent_atoms.tolist())),
                    expected["equivalent_atom_orbit_count"],
                )
                self.assertEqual(
                    [str(value) for value in dataset.wyckoffs],
                    expected["wyckoff_letters"],
                )

                record = get_crystallographic_space_group(dataset.number)
                self.assertEqual(record.international_short, dataset.international)
                self.assertEqual(record.point_group_hm, dataset.pointgroup)
                self.assertEqual(record.crystal_system, expected["crystal_system"])
                hall_setting = next(
                    setting
                    for setting in record.hall_settings
                    if setting.hall_number == dataset.hall_number
                )
                self.assertEqual(hall_setting.hall_symbol, dataset.hall)
                self.assertEqual(hall_setting.choice, dataset.choice)
                self.assertEqual(
                    hall_setting.operation_count,
                    expected["hall_setting_operation_count"],
                )
