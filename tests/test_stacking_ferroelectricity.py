"""Scientific regression tests for stacking-ferroelectricity symmetry kernels."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import json
import unittest

import numpy as np

from group_theory_operations import (
    iter_crystallographic_layer_groups,
    point_group_operations,
)
from group_theory_operations.cli import main
from group_theory_operations.stacking import (
    bravais_lattice_operations,
    equivalent_interface_orbit,
    equivalent_interface_translation,
    layer_group_polarization,
    partition_left_cosets,
    polarization_space,
    polarization_switch,
    preserved_recursive_stacking_operations,
    preserves_recursive_stacking_step,
    stacking_rotation_cosets,
)


class LayerGroupPolarizationTests(unittest.TestCase):
    # Independent scientific answer set: Table S2 of
    # Phys. Rev. Lett. 130, 146801 (2023), Supplemental Material.
    EXPECTED = {
        "IP": [4, 5, *range(8, 11), *range(27, 37)],
        "OP": [3, *range(23, 27), 49, 55, 56, 65, 69, 70, 73, 77],
        "CP": [1, 11, 12, 13],
        "NP": [
            2,
            6,
            7,
            *range(14, 23),
            *range(37, 49),
            *range(50, 55),
            *range(57, 65),
            66,
            67,
            68,
            71,
            72,
            74,
            75,
            76,
            78,
            79,
            80,
        ],
    }

    def test_all_80_layer_groups_match_published_polar_types(self) -> None:
        actual = {key: [] for key in self.EXPECTED}
        for number in range(1, 81):
            actual[layer_group_polarization(number).polar_type].append(number)
        self.assertEqual(actual, self.EXPECTED)

    def test_representative_fixed_space_dimensions(self) -> None:
        self.assertEqual(layer_group_polarization(1).dimension, 3)
        self.assertEqual(layer_group_polarization(8).in_plane_dimension, 1)
        self.assertEqual(layer_group_polarization(3).out_of_plane_dimension, 1)
        self.assertEqual(layer_group_polarization(68).dimension, 0)

    def test_all_layer_hall_settings_keep_the_same_polar_type(self) -> None:
        for group in iter_crystallographic_layer_groups():
            expected = layer_group_polarization(group.number).polar_type
            for setting in group.hall_settings:
                self.assertEqual(
                    layer_group_polarization(
                        group.number,
                        layer_hall_number=setting.layer_hall_number,
                    ).polar_type,
                    expected,
                )


class RotationCosetTests(unittest.TestCase):
    def test_ci_rotation_counts_grow_with_lattice_symmetry(self) -> None:
        # Table SI of Phys. Rev. B 111, 224102 (2025), Supplemental
        # Material: LG2 (Ci) has 2, 4, 4, 8, and 12 orientation classes.
        expected = {"oP": 2, "rP": 4, "rC": 4, "sP": 8, "hP": 12}
        self.assertEqual(
            {
                lattice: len(stacking_rotation_cosets("-1", lattice))
                for lattice in expected
            },
            expected,
        )

    def test_left_cosets_support_nonnormal_subgroups(self) -> None:
        lattice = bravais_lattice_operations("sP")
        by_name = {operation.name: operation for operation in lattice}
        # C2h with its twofold axis along x, the MoH3 embedding used to
        # motivate the improved left-coset construction.
        subgroup = [by_name[name] for name in ("1", "m_100", "-1", "2_100")]
        cosets = partition_left_cosets(
            (operation.matrix_fractional for operation in lattice),
            (operation.matrix_fractional for operation in subgroup),
        )
        self.assertEqual(len(cosets), 4)
        self.assertTrue(all(len(coset.members) == 4 for coset in cosets))
        self.assertEqual(
            len({matrix for coset in cosets for matrix in coset.members}), 16
        )

    def test_incompatible_coordinate_embedding_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "not embedded"):
            stacking_rotation_cosets("mmm", "hP")

    def test_nonclosed_subgroup_is_rejected(self) -> None:
        lattice = bravais_lattice_operations("sP")
        by_name = {operation.name: operation for operation in lattice}
        with self.assertRaisesRegex(ValueError, "subgroup operations are not closed"):
            partition_left_cosets(
                (operation.matrix_fractional for operation in lattice),
                (
                    by_name[name].matrix_fractional
                    for name in ("1", "4+_001")
                ),
            )


class InterfaceEquivalenceTests(unittest.TestCase):
    def test_bn_ab_and_ba_are_related_by_layer_exchange(self) -> None:
        mirror_z = next(
            operation
            for operation in point_group_operations("6/mmm")
            if operation.name == "m_001"
        )
        self.assertTrue(
            np.allclose(
                equivalent_interface_translation((1 / 3, 2 / 3), mirror_z),
                (2 / 3, 1 / 3),
            )
        )

    def test_equivalent_interface_orbit_is_unique_and_periodic(self) -> None:
        operations = point_group_operations("6/mmm")
        orbit = equivalent_interface_orbit((1 / 3, 2 / 3), operations)
        expected = ((1 / 3, 2 / 3), (2 / 3, 1 / 3))
        self.assertEqual(len(orbit), 2)
        self.assertTrue(np.allclose(orbit, expected))
        self.assertTrue(
            np.allclose(
                equivalent_interface_orbit((4 / 3, -1 / 3), operations),
                orbit,
            )
        )

    def test_centered_translation_is_modded_out(self) -> None:
        identity = np.eye(3)
        self.assertEqual(
            equivalent_interface_translation(
                (0.75, 0.25), identity, centering="C"
            ),
            (0.25, 0.75),
        )


class MultilayerCriterionTests(unittest.TestCase):
    def test_layer_preserving_and_exchanging_cases(self) -> None:
        identity = np.eye(3)
        inversion = -np.eye(3)
        first = (1 / 3, 2 / 3)
        second = (2 / 3, 1 / 3)
        self.assertTrue(preserves_recursive_stacking_step(identity, first, second))
        self.assertTrue(
            preserves_recursive_stacking_step(inversion, first, second)
        )
        self.assertFalse(
            preserves_recursive_stacking_step(inversion, first, first)
        )

    def test_graphene_abc_and_aba_reproduce_expected_point_groups(self) -> None:
        operations = point_group_operations("6/mmm")
        abc = preserved_recursive_stacking_operations(
            operations, (1 / 3, 2 / 3), (2 / 3, 1 / 3)
        )
        aba = preserved_recursive_stacking_operations(
            operations, (1 / 3, 2 / 3), (1 / 3, 2 / 3)
        )
        self.assertEqual(len(abc), 12)  # D3d
        self.assertEqual(len(aba), 12)  # D3h
        self.assertIn("-1", {operation.name for operation in abc})
        self.assertNotIn("m_001", {operation.name for operation in abc})
        self.assertIn("m_001", {operation.name for operation in aba})
        self.assertNotIn("-1", {operation.name for operation in aba})
        self.assertEqual(
            polarization_space(
                operation.matrix_fractional for operation in abc
            ).polar_type,
            "NP",
        )
        self.assertEqual(
            polarization_space(
                operation.matrix_fractional for operation in aba
            ).polar_type,
            "NP",
        )


class SwitchingModeTests(unittest.TestCase):
    def test_op_reversal_with_ip_unchanged(self) -> None:
        result = polarization_switch(np.diag((1, 1, -1)), (1, 0, 1))
        self.assertEqual(result.out_of_plane, "reversed")
        self.assertEqual(result.in_plane_angle_degrees, 0.0)

    def test_synchronous_120_degree_switch(self) -> None:
        angle = np.deg2rad(120)
        operation = np.asarray(
            (
                (np.cos(angle), -np.sin(angle), 0),
                (np.sin(angle), np.cos(angle), 0),
                (0, 0, -1),
            )
        )
        result = polarization_switch(operation, (1, 0, 1))
        self.assertEqual(result.out_of_plane, "reversed")
        self.assertEqual(result.in_plane_angle_degrees, 120.0)


class StackingCliTests(unittest.TestCase):
    def test_layer_polarity_query(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["layer-polarity", "LG68", "--json"]), 0)
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["layer_group_number"], 68)
        self.assertEqual(payload["polar_type"], "NP")

    def test_stacking_rotation_query(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(
                    [
                        "stacking-rotations",
                        "-1",
                        "--lattice",
                        "rP",
                        "--json",
                    ]
                ),
                0,
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["bravais_point_group"], "mmm")
        self.assertEqual(payload["rotation_class_count"], 4)


if __name__ == "__main__":
    unittest.main()
