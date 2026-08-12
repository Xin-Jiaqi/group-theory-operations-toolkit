"""Public real-structure symmetry-context tests."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from group_theory_operations import (
    GroupDataError,
    StructureSymmetryContext,
    classify_structure_symmetry,
)
from group_theory_operations.seitz import equivalent

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "spglib_real_structures_v2.5.0.json"


class StructureDouble:
    def __init__(
        self,
        *,
        lattice,
        species,
        fractional_coordinates,
        pbc=(True, True, True),
    ):
        self.lattice = tuple(tuple(row) for row in lattice)
        self.species = tuple(species)
        self.fractional_coordinates = tuple(
            tuple(row) for row in fractional_coordinates
        )
        self.pbc = tuple(pbc)


def _load_fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _structure(item, *, lattice=None, positions=None, species=None, pbc=None):
    return StructureDouble(
        lattice=item["lattice"] if lattice is None else lattice,
        species=item["species"] if species is None else species,
        fractional_coordinates=(
            item["fractional_coordinates"] if positions is None else positions
        ),
        pbc=(True, True, True) if pbc is None else pbc,
    )


@unittest.skipUnless(spglib is not None, "spglib not installed")
class PublicStructureSymmetryTests(unittest.TestCase):
    def test_all_real_fixtures_return_validated_typed_contexts(self) -> None:
        fixture = _load_fixture()
        settings = fixture["classification"]
        for item in fixture["structures"]:
            with self.subTest(item=item["id"]):
                context = classify_structure_symmetry(
                    _structure(item),
                    symprec=settings["symprec"],
                    angle_tolerance=settings["angle_tolerance"],
                )
                expected = item["expected"]
                self.assertIsInstance(context, StructureSymmetryContext)
                self.assertEqual(context.space_group.ita_number, expected["ita_number"])
                self.assertEqual(context.hall_setting.hall_number, expected["hall_number"])
                self.assertEqual(
                    len(context.input_operations),
                    expected["input_cell_operation_count"],
                )
                self.assertEqual(
                    len(context.standard_operations),
                    expected["hall_setting_operation_count"],
                )
                self.assertEqual(context.wyckoff_letters, tuple(expected["wyckoff_letters"]))
                self.assertEqual(
                    len(set(context.equivalent_atoms)),
                    expected["equivalent_atom_orbit_count"],
                )
                self.assertEqual(
                    len(context.crystallographic_orbits),
                    len(item["species"]),
                )
                self.assertEqual(
                    len(context.site_mappings),
                    len(context.input_operations),
                )
                self.assertEqual(
                    len(context.standardized_species),
                    len(context.standardized_fractional_coordinates),
                )
                self.assertAlmostEqual(
                    context.operation_count_ratio,
                    context.cell_volume_ratio,
                    places=8,
                )
                transformed = tuple(
                    context.operation_to_standard(operation)
                    for operation in context.input_operations
                )
                self.assertTrue(
                    all(
                        any(
                            equivalent(operation, standard)
                            for standard in context.standard_operations
                        )
                        for operation in transformed
                    )
                )
                self.assertEqual(context.backend_name, "spglib")
                self.assertEqual(context.backend_version, spglib.__version__)

                if item["id"] == "sio2_monoclinic":
                    self.assertAlmostEqual(context.cell_volume_ratio, 2.0)
                    self.assertEqual(
                        len(context.input_operations),
                        2 * len(context.standard_operations),
                    )

        with self.assertRaises(FrozenInstanceError):
            context.symprec = 1.0  # type: ignore[misc]

    def test_changed_basis_and_origin_recover_the_same_hall_setting(self) -> None:
        item = _load_fixture()["structures"][2]
        basis_change = np.array(
            [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
        )
        origin_shift = np.array([0.137, 0.271, 0.419])
        lattice = (
            np.linalg.inv(basis_change).T
            @ np.asarray(item["lattice"], dtype=np.float64)
        )
        positions = (
            np.asarray(item["fractional_coordinates"]) @ basis_change.T
            + origin_shift
        ) % 1.0
        context = classify_structure_symmetry(_structure(item, lattice=lattice, positions=positions))
        self.assertEqual(context.space_group.ita_number, item["expected"]["ita_number"])
        self.assertEqual(context.hall_setting.hall_number, item["expected"]["hall_number"])
        self.assertGreater(max(abs(value) for value in context.origin_shift), 1.0e-3)
        standardized = context.to_standard_fractional(positions)
        self.assertEqual(len(standardized), len(positions))

    def test_small_metric_distortion_is_handled_at_declared_tolerance(self) -> None:
        item = _load_fixture()["structures"][-1]
        lattice = np.asarray(item["lattice"], dtype=np.float64)
        lattice[0, 0] += 2.0e-6
        context = classify_structure_symmetry(
            _structure(item, lattice=lattice),
            symprec=1.0e-5,
        )
        self.assertEqual(context.space_group.ita_number, 221)
        self.assertEqual(len(context.input_operations), 48)

    def test_primitive_input_expands_to_centered_standard_operations(self) -> None:
        rock_salt_primitive = StructureDouble(
            lattice=[
                [0.0, 2.82, 2.82],
                [2.82, 0.0, 2.82],
                [2.82, 2.82, 0.0],
            ],
            species=["Na", "Cl"],
            fractional_coordinates=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        )
        context = classify_structure_symmetry(rock_salt_primitive)
        self.assertEqual(context.space_group.ita_number, 225)
        self.assertEqual(context.hall_setting.centering, "F")
        self.assertEqual(len(context.input_operations), 48)
        self.assertEqual(len(context.standard_operations), 192)
        self.assertAlmostEqual(context.cell_volume_ratio, 0.25)
        self.assertAlmostEqual(context.operation_count_ratio, 0.25)
        self.assertEqual(len(context.standardized_species), 8)
        self.assertEqual(context.standardized_species.count("Na"), 4)
        self.assertEqual(context.standardized_species.count("Cl"), 4)

    def test_invalid_contracts_and_tolerances_are_rejected(self) -> None:
        item = _load_fixture()["structures"][-1]
        with self.assertRaisesRegex(ValueError, "positive finite"):
            classify_structure_symmetry(_structure(item), symprec=0.0)
        with self.assertRaisesRegex(ValueError, "positive finite angle"):
            classify_structure_symmetry(_structure(item), angle_tolerance=0.0)
        with self.assertRaisesRegex(GroupDataError, "three periodic axes"):
            classify_structure_symmetry(
                _structure(item, pbc=(True, True, False))
            )
        with self.assertRaisesRegex(GroupDataError, "non-empty strings"):
            classify_structure_symmetry(
                _structure(item, species=[1] * len(item["species"]))
            )
        with self.assertRaisesRegex(GroupDataError, "finite 3x3"):
            classify_structure_symmetry(
                _structure(item, lattice=np.eye(2))
            )

    def test_spglib_is_an_optional_lazy_dependency(self) -> None:
        item = _load_fixture()["structures"][-1]
        with patch(
            "group_theory_operations.structure_symmetry.importlib.import_module",
            side_effect=ModuleNotFoundError,
        ):
            with self.assertRaisesRegex(ImportError, r"\[structure\]"):
                classify_structure_symmetry(_structure(item))


if __name__ == "__main__":
    unittest.main()
