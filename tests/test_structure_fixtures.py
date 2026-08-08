"""Real-structure classification checks against the space-group registry."""

from __future__ import annotations

import json
import math
from pathlib import Path
import unittest

import numpy as np

from group_theory_operations.space_groups import get_crystallographic_space_group
from group_theory_operations.seitz import (
    SeitzOp,
    closure,
    equivalent,
    transform_seitz_coordinates,
)
from group_theory_operations.structure import (
    _seitz_site_mapping,
    _site_orbits,
    _transform_fractional_coordinates,
)

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


class FixtureStructureRecord:
    def __init__(
        self,
        *,
        lattice,
        species,
        fractional_coordinates,
        pbc=(True, True, True),
        selective_dynamics=None,
        length_unit="angstrom",
    ):
        self.lattice = tuple(tuple(row) for row in lattice)
        self.species = tuple(species)
        self.fractional_coordinates = tuple(
            tuple(row) for row in fractional_coordinates
        )
        self.pbc = tuple(pbc)
        self.selective_dynamics = selective_dynamics
        self.length_unit = length_unit

    @classmethod
    def from_fractional(cls, **kwargs):
        return cls(**kwargs)

    def wrapped(self):
        return self.from_fractional(
            lattice=self.lattice,
            species=self.species,
            fractional_coordinates=[
                [value % 1.0 for value in position]
                for position in self.fractional_coordinates
            ],
            pbc=self.pbc,
            selective_dynamics=self.selective_dynamics,
            length_unit=self.length_unit,
        )


def _partition(labels) -> set[frozenset[int]]:
    groups: dict[int, set[int]] = {}
    for index, label in enumerate(labels):
        groups.setdefault(int(label), set()).add(index)
    return {frozenset(indices) for indices in groups.values()}


def _unique_operations(operations) -> tuple[SeitzOp, ...]:
    unique: list[SeitzOp] = []
    for operation in operations:
        if not any(equivalent(operation, known) for known in unique):
            unique.append(operation)
    return tuple(unique)


def _same_operation_set(left, right) -> bool:
    left_unique = _unique_operations(left)
    right_unique = _unique_operations(right)
    return len(left_unique) == len(right_unique) and all(
        any(equivalent(operation, candidate) for candidate in right_unique)
        for operation in left_unique
    )


def _maximum_same_type_periodic_distance(
    source_positions,
    source_types,
    target_positions,
    target_types,
    lattice,
) -> float:
    target_positions = np.asarray(target_positions, dtype=np.float64)
    target_types = np.asarray(target_types)
    lattice = np.asarray(lattice, dtype=np.float64)
    maximum = 0.0
    for position, type_number in zip(source_positions, source_types):
        candidates = target_positions[target_types == type_number]
        if not len(candidates):
            return math.inf
        difference = np.asarray(position) - candidates
        difference -= np.rint(difference)
        maximum = max(
            maximum,
            float(np.min(np.linalg.norm(difference @ lattice, axis=1))),
        )
    return maximum


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

    def test_all_detected_operations_are_species_aware_automorphisms(self) -> None:
        fixture = _load_fixture()
        settings = fixture["classification"]
        for item in fixture["structures"]:
            with self.subTest(item=item["id"], formula=item["formula"]):
                structure = FixtureStructureRecord(
                    lattice=item["lattice"],
                    species=item["species"],
                    fractional_coordinates=item["fractional_coordinates"],
                )
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
                mappings = []
                for rotation, translation in zip(
                    dataset.rotations, dataset.translations
                ):
                    mapping = _seitz_site_mapping(
                        structure,
                        SeitzOp(rotation, translation),
                        tolerance=settings["symprec"],
                    )
                    self.assertIsNotNone(mapping)
                    assert mapping is not None
                    self.assertTrue(
                        all(
                            item["species"][source] == item["species"][target]
                            for source, target in enumerate(mapping)
                        )
                    )
                    mappings.append(mapping)

                orbit_labels = _site_orbits(mappings)
                self.assertEqual(
                    _partition(orbit_labels),
                    _partition(dataset.equivalent_atoms),
                )
                self.assertEqual(
                    len(set(orbit_labels)),
                    item["expected"]["equivalent_atom_orbit_count"],
                )
                for orbit in set(orbit_labels):
                    letters = {
                        item["expected"]["wyckoff_letters"][index]
                        for index, label in enumerate(orbit_labels)
                        if label == orbit
                    }
                    self.assertEqual(len(letters), 1)

    def test_input_to_standard_setting_round_trips_end_to_end(self) -> None:
        fixture = _load_fixture()
        settings = fixture["classification"]
        basis_change = np.array(
            [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
        )
        shifted_origin = np.array([0.137, 0.271, 0.419])

        for item in fixture["structures"]:
            with self.subTest(item=item["id"], coordinates="source"):
                self._assert_standard_setting_round_trip(
                    item,
                    np.asarray(item["lattice"], dtype=np.float64),
                    np.asarray(item["fractional_coordinates"], dtype=np.float64),
                    settings,
                )

            transformed_lattice = (
                np.linalg.inv(basis_change).T
                @ np.asarray(item["lattice"], dtype=np.float64)
            )
            transformed_positions = _transform_fractional_coordinates(
                item["fractional_coordinates"],
                basis_change,
                shifted_origin,
            )
            with self.subTest(item=item["id"], coordinates="changed-basis-origin"):
                dataset = self._assert_standard_setting_round_trip(
                    item,
                    transformed_lattice,
                    np.asarray(transformed_positions),
                    settings,
                )
                if dataset.number != 1:
                    self.assertGreater(
                        float(np.max(np.abs(dataset.origin_shift))),
                        1.0e-3,
                    )

    def _assert_standard_setting_round_trip(
        self,
        item,
        lattice,
        positions,
        settings,
    ):
        type_numbers = np.asarray(_type_numbers(item["species"]))
        dataset = spglib.get_symmetry_dataset(
            (lattice, positions, type_numbers),
            symprec=settings["symprec"],
            angle_tolerance=settings["angle_tolerance"],
        )
        self.assertIsNotNone(dataset)
        assert dataset is not None
        self.assertEqual(dataset.number, item["expected"]["ita_number"])
        self.assertEqual(dataset.hall_number, item["expected"]["hall_number"])

        matrix = np.asarray(dataset.transformation_matrix, dtype=np.float64)
        shift = np.asarray(dataset.origin_shift, dtype=np.float64)
        rotation = np.asarray(dataset.std_rotation_matrix, dtype=np.float64)
        np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=1.0e-12)
        predicted_standard_lattice = np.linalg.inv(matrix).T @ lattice @ rotation.T
        np.testing.assert_allclose(
            predicted_standard_lattice,
            dataset.std_lattice,
            atol=10.0 * settings["symprec"],
            rtol=1.0e-10,
        )

        unwrapped_standardized_positions = _transform_fractional_coordinates(
            positions,
            matrix,
            shift,
            wrap=False,
        )
        inverse_matrix = np.linalg.inv(matrix)
        restored_positions = _transform_fractional_coordinates(
            unwrapped_standardized_positions,
            inverse_matrix,
            -inverse_matrix @ shift,
            wrap=False,
        )
        np.testing.assert_allclose(restored_positions, positions, atol=1.0e-12)
        standardized_positions = _transform_fractional_coordinates(
            positions,
            matrix,
            shift,
        )

        tolerance = 10.0 * settings["symprec"]
        self.assertLessEqual(
            _maximum_same_type_periodic_distance(
                standardized_positions,
                type_numbers,
                dataset.std_positions,
                dataset.std_types,
                dataset.std_lattice,
            ),
            tolerance,
        )
        self.assertLessEqual(
            _maximum_same_type_periodic_distance(
                dataset.std_positions,
                dataset.std_types,
                standardized_positions,
                type_numbers,
                dataset.std_lattice,
            ),
            tolerance,
        )

        transformed_operations = _unique_operations(
            transform_seitz_coordinates(
                SeitzOp(operation, translation),
                matrix,
                shift,
            )
            for operation, translation in zip(
                dataset.rotations,
                dataset.translations,
            )
        )
        record = get_crystallographic_space_group(dataset.number)
        hall_setting = next(
            setting
            for setting in record.hall_settings
            if setting.hall_number == dataset.hall_number
        )
        registry_operations = tuple(
            closure(
                SeitzOp.from_dict(generator)
                for generator in hall_setting.generators
            )
        )
        self.assertTrue(
            _same_operation_set(transformed_operations, registry_operations)
        )
        self.assertAlmostEqual(
            len(dataset.rotations) / len(transformed_operations),
            abs(float(np.linalg.det(matrix))),
            places=8,
        )
        return dataset
