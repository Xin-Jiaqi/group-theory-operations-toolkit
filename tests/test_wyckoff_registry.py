"""Complete Wyckoff registry and embedded orbit-splitting tests."""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import unittest

import numpy as np
from jsonschema import Draft202012Validator

from group_theory_operations import (
    GroupDataError,
    get_wyckoff_setting,
    iter_wyckoff_settings,
    load_wyckoff_registry,
    split_wyckoff_orbit,
)
from group_theory_operations.cli import main
from group_theory_operations.seitz import SeitzOp, closure
from group_theory_operations.space_groups import get_crystallographic_space_group

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "crystallographic_wyckoff_positions.json"
SCHEMA_PATH = ROOT / "schema" / "crystallographic-wyckoff-positions-v1.schema.json"
SOURCE_PATH = ROOT / "scripts" / "sources" / "spglib-Wyckoff-v2.5.0.csv"


def _operation_orbit(hall_number: int, coordinate: np.ndarray) -> list[np.ndarray]:
    space_group = get_crystallographic_space_group(
        get_wyckoff_setting(hall_number).ita_number
    )
    setting = next(
        item for item in space_group.hall_settings if item.hall_number == hall_number
    )
    operations = closure(SeitzOp.from_dict(item) for item in setting.generators)
    unique: list[np.ndarray] = []
    for operation in operations:
        image = operation.apply(coordinate)
        if not any(
            np.max(np.abs((image - known) - np.rint(image - known))) < 1.0e-9
            for known in unique
        ):
            unique.append(image)
    return unique


class WyckoffRegistryContractTests(unittest.TestCase):
    def test_schema_provenance_and_complete_counts(self) -> None:
        registry = load_wyckoff_registry()
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        errors = sorted(
            Draft202012Validator(schema).iter_errors(registry), key=str
        )
        self.assertEqual(errors, [])
        self.assertEqual(registry["counts"]["hall_settings"], 530)
        self.assertEqual(registry["counts"]["wyckoff_positions"], 3467)
        self.assertEqual(registry["counts"]["representative_maps"], 15117)
        self.assertEqual(registry["counts"]["expanded_coordinate_maps"], 24295)
        self.assertEqual(
            hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
            registry["source"]["sha256"],
        )

    def test_all_hall_settings_match_space_group_registry(self) -> None:
        settings = tuple(iter_wyckoff_settings())
        self.assertEqual(tuple(item.hall_number for item in settings), tuple(range(1, 531)))
        for setting in settings:
            space_group = get_crystallographic_space_group(setting.ita_number)
            hall = next(
                item
                for item in space_group.hall_settings
                if item.hall_number == setting.hall_number
            )
            self.assertEqual(setting.hall_symbol, hall.hall_symbol)
            self.assertEqual(setting.choice, hall.choice)
            self.assertEqual(
                len({position.letter for position in setting.positions}),
                len(setting.positions),
            )
            for position in setting.positions:
                self.assertEqual(
                    position.multiplicity,
                    len(position.representative_maps)
                    * len(setting.centering_translation_numerators),
                )

    def test_all_generated_coordinates_equal_registered_operation_orbits(self) -> None:
        parameters = np.asarray([0.173, 0.287, 0.419])
        for setting in iter_wyckoff_settings():
            hall_number = setting.hall_number
            for position in setting.positions:
                coordinates = setting.coordinates(position.letter, parameters)
                self.assertEqual(len(coordinates), position.multiplicity)
                orbit = _operation_orbit(hall_number, np.asarray(coordinates[0]))
                self.assertEqual(len(orbit), position.multiplicity)
                self.assertTrue(
                    all(
                        any(
                            np.max(
                                np.abs(
                                    (np.asarray(coordinate) - image)
                                    - np.rint(np.asarray(coordinate) - image)
                                )
                            )
                            < 1.0e-9
                            for image in orbit
                        )
                        for coordinate in coordinates
                    )
                )

    @unittest.skipUnless(spglib is not None, "spglib not installed")
    def test_labels_and_site_symmetries_round_trip_through_spglib(self) -> None:
        for hall_number in (1, 2, 227):
            setting = get_wyckoff_setting(hall_number)
            lattice = (
                np.asarray([[4.7, 0.0, 0.0], [0.4, 5.3, 0.0], [0.2, 0.7, 6.1]])
                if hall_number in {1, 2}
                else np.diag([4.7, 5.3, 6.1])
            )
            general = setting.positions[0]
            anchor_one = setting.coordinates(general.letter, (0.137, 0.263, 0.389))
            anchor_two = setting.coordinates(general.letter, (0.191, 0.317, 0.443))
            for position in setting.positions:
                target = setting.coordinates(position.letter)
                coordinates = (*anchor_one, *anchor_two, *target)
                species = (
                    [1] * len(anchor_one)
                    + [2] * len(anchor_two)
                    + [3] * len(target)
                )
                dataset = spglib.get_symmetry_dataset(
                    (lattice, coordinates, species),
                    symprec=1.0e-5,
                )
                self.assertEqual(int(dataset.hall_number), hall_number)
                target_slice = slice(len(anchor_one) + len(anchor_two), None)
                self.assertEqual(
                    set(dataset.wyckoffs[target_slice]),
                    {position.letter},
                )
                self.assertEqual(
                    {
                        str(item).strip()
                        for item in dataset.site_symmetry_symbols[target_slice]
                    },
                    {position.site_symmetry},
                )

    def test_invalid_provenance_and_corrupt_multiplicity_are_rejected(self) -> None:
        registry = load_wyckoff_registry()
        wrong_source = json.loads(json.dumps(registry))
        wrong_source["source"]["sha256"] = "0" * 64
        wrong_multiplicity = json.loads(json.dumps(registry))
        wrong_multiplicity["hall_settings"][0]["positions"][0]["multiplicity"] = 2
        with self.assertRaises(GroupDataError):
            load_wyckoff_registry_from_data(wrong_source)
        with self.assertRaises(GroupDataError):
            list(iter_wyckoff_settings(wrong_multiplicity))

    def test_malformed_registry_root_and_parameters_are_rejected(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registry.json"
            path.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(GroupDataError, "JSON object"):
                load_wyckoff_registry(path)
        with self.assertRaisesRegex(GroupDataError, "three finite numbers"):
            get_wyckoff_setting(2).coordinates("i", ("x", "y", "z"))


def load_wyckoff_registry_from_data(data: dict) -> dict:
    """Exercise the public file loader without adding a second validation path."""

    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "registry.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return load_wyckoff_registry(path)


class WyckoffSplittingTests(unittest.TestCase):
    def test_p_minus_one_general_orbit_splits_into_two_p_one_orbits(self) -> None:
        result = split_wyckoff_orbit(2, "i", 1)
        self.assertEqual(result.parent_label, "2i")
        self.assertEqual(result.subgroup_index, 2)
        self.assertEqual([item.label for item in result.child_orbits], ["1a", "1a"])
        self.assertEqual(sum(item.multiplicity for item in result.child_orbits), 2)
        self.assertEqual(result.translation_subgroup_index, 1)
        self.assertEqual(result.point_group_index, 2)

    def test_pmmm_general_orbit_splitting_conserves_multiplicity(self) -> None:
        result = split_wyckoff_orbit(227, "A", 2)
        self.assertEqual(result.subgroup_index, 4)
        self.assertEqual([item.label for item in result.child_orbits], ["2i"] * 4)
        self.assertEqual(sum(item.multiplicity for item in result.child_orbits), 8)

    def test_centered_subgroup_matching_includes_centering_translations(self) -> None:
        result = split_wyckoff_orbit(400, "u", 349)
        self.assertEqual(result.subgroup_index, 4)
        self.assertEqual([item.label for item in result.child_orbits], ["4d"] * 4)
        self.assertEqual(sum(item.multiplicity for item in result.child_orbits), 16)

    def test_centering_loss_is_counted_in_translation_subgroup_index(self) -> None:
        result = split_wyckoff_orbit(523, "l", 517)
        self.assertEqual(result.conventional_cell_index, 1)
        self.assertEqual(result.translation_subgroup_index, 4)
        self.assertEqual(result.point_group_index, 1)
        self.assertEqual(result.subgroup_index, 4)
        self.assertEqual([item.label for item in result.child_orbits], ["48n"] * 4)

    def test_published_i23_to_p23_wyckoff_splitting(self) -> None:
        result = split_wyckoff_orbit(491, "a", 489)
        self.assertEqual(result.parent_label, "2a")
        self.assertEqual(result.translation_subgroup_index, 2)
        self.assertEqual(result.point_group_index, 1)
        self.assertEqual([item.label for item in result.child_orbits], ["1a", "1b"])

    def test_unembedded_setting_is_rejected(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "basis/origin transformation"):
            split_wyckoff_orbit(3, "e", 4)

    def test_explicit_axis_transformation_embeds_alternative_settings(self) -> None:
        transformation = (
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0),
        )
        result = split_wyckoff_orbit(
            3,
            "e",
            4,
            subgroup_transformation_matrix=transformation,
        )
        self.assertEqual(result.subgroup_index, 1)
        self.assertEqual([item.label for item in result.child_orbits], ["2e"])
        self.assertEqual(result.subgroup_transformation_matrix, transformation)

    def test_explicit_origin_shift_embeds_alternative_settings(self) -> None:
        result = split_wyckoff_orbit(
            39,
            "a",
            43,
            subgroup_origin_shift=(0.0, 0.25, 0.0),
        )
        self.assertEqual(result.subgroup_index, 1)
        self.assertEqual([item.label for item in result.child_orbits], ["4a"])
        self.assertEqual(result.subgroup_origin_shift, (0.0, 0.25, 0.0))

    def test_specialized_parent_parameters_are_rejected(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "special Wyckoff letter"):
            split_wyckoff_orbit(2, "i", 1, parameters=(0.0, 0.0, 0.0))

    def test_doubled_cell_expands_parent_orbit_and_constructs_cosets(self) -> None:
        result = split_wyckoff_orbit(
            2,
            "i",
            1,
            subgroup_transformation_matrix=(
                (0.5, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ),
        )
        self.assertEqual(result.conventional_cell_index, 2)
        self.assertEqual(result.translation_subgroup_index, 2)
        self.assertEqual(result.point_group_index, 2)
        self.assertEqual(result.subgroup_index, 4)
        self.assertEqual(
            result.subgroup_supercell_matrix,
            ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        )
        self.assertEqual(
            result.parent_translation_cosets,
            ((0, 0, 0), (1, 0, 0)),
        )
        self.assertEqual([item.label for item in result.child_orbits], ["1a"] * 4)
        self.assertEqual(sum(item.multiplicity for item in result.child_orbits), 4)

    def test_non_diagonal_supercell_cosets_are_complete(self) -> None:
        result = split_wyckoff_orbit(
            1,
            "a",
            1,
            subgroup_transformation_matrix=(
                (0.5, -0.5, 0.0),
                (0.5, 0.5, 0.0),
                (0.0, 0.0, 1.0),
            ),
            subgroup_origin_shift=(0.125, 0.25, 0.0),
        )
        self.assertEqual(result.conventional_cell_index, 2)
        self.assertEqual(result.subgroup_index, 2)
        self.assertEqual(len(result.parent_translation_cosets), 2)
        self.assertEqual([item.label for item in result.child_orbits], ["1a", "1a"])

    def test_noncrystallographic_basis_change_is_rejected(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "integer parent-cell supercell"):
            split_wyckoff_orbit(
                2,
                "i",
                1,
                subgroup_transformation_matrix=(
                    (0.5, 0.0, 0.0),
                    (0.0, 2.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
            )

    def test_cell_contraction_is_rejected_as_non_subgroup_embedding(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "integer parent-cell supercell"):
            split_wyckoff_orbit(
                2,
                "i",
                1,
                subgroup_transformation_matrix=(
                    (2.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
            )

    def test_unbounded_supercell_expansion_is_rejected(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "supported maximum"):
            split_wyckoff_orbit(
                1,
                "a",
                1,
                subgroup_transformation_matrix=(
                    (1.0 / 4097.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
            )

    def test_unbounded_expanded_orbit_is_rejected_before_construction(self) -> None:
        with self.assertRaisesRegex(GroupDataError, "expanded parent Wyckoff orbit"):
            split_wyckoff_orbit(
                227,
                "A",
                1,
                subgroup_transformation_matrix=(
                    (1.0 / 1024.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
            )

    def test_cli_queries_registry_and_splitting(self) -> None:
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["wyckoff-positions", "2", "i", "--json"]), 0)
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["positions"][0]["label"], "2i")
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["wyckoff-split", "2", "i", "1", "--json"]), 0)
        payload = json.loads(output.getvalue())
        self.assertEqual([item["label"] for item in payload["child_orbits"]], ["1a", "1a"])
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(
                    [
                        "wyckoff-split",
                        "3",
                        "e",
                        "4",
                        "--transformation-matrix",
                        "1",
                        "0",
                        "0",
                        "0",
                        "0",
                        "1",
                        "0",
                        "1",
                        "0",
                        "--json",
                    ]
                ),
                0,
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["child_orbits"][0]["label"], "2e")


if __name__ == "__main__":
    unittest.main()
