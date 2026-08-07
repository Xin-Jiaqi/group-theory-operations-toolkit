"""Tests for complex antiunitary maps and magnetic optical sectors."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
import unittest

from group_theory_operations import (
    GroupDataError,
    antiunitary_equivariant_map_basis,
    iter_magnetic_point_groups,
    load_database,
    load_magnetic_point_group_registry,
    magnetic_response_tensor_basis,
)
from group_theory_operations.cli import main


I1 = (((1.0,),),)


class AntiunitaryMapTests(unittest.TestCase):
    def test_pure_time_reversal_separates_real_and_imaginary_parts(self):
        even = antiunitary_equivariant_map_basis(
            I1, I1, [True], antiunitary_character="even"
        )
        odd = antiunitary_equivariant_map_basis(
            I1, I1, [True], antiunitary_character="odd"
        )
        self.assertEqual((even.real_dimension, even.imaginary_dimension), (1, 0))
        self.assertEqual((odd.real_dimension, odd.imaginary_dimension), (0, 1))
        self.assertEqual(even.dimension, 1)
        self.assertEqual(odd.to_dict()["antiunitary_character"], "odd")

    def test_unitary_operations_do_not_conjugate(self):
        basis = antiunitary_equivariant_map_basis(I1, I1, [False])
        self.assertEqual((basis.real_dimension, basis.imaginary_dimension), (1, 1))

    def test_rejects_misaligned_or_non_boolean_labels(self):
        with self.assertRaises(ValueError):
            antiunitary_equivariant_map_basis(I1, I1, [])
        with self.assertRaises(ValueError):
            antiunitary_equivariant_map_basis(I1, I1, [1])


class MagneticOpticalResponseTests(unittest.TestCase):
    responses = (
        "normal_shift_current",
        "magnetic_shift_current",
        "normal_injection_current",
        "magnetic_injection_current",
        "shg_even",
        "shg_odd",
    )

    def test_all_122_groups_and_six_sectors_are_solved(self):
        database = load_database()
        registry = load_magnetic_point_group_registry()
        for group in iter_magnetic_point_groups(registry):
            for response in self.responses:
                result = magnetic_response_tensor_basis(
                    group.number,
                    response,
                    database=database,
                    registry=registry,
                )
                self.assertEqual(result.magnetic_point_group_number, group.number)
                self.assertEqual(result.shape[0], 3)
                self.assertIn(result.shape[1], (3, 6))

    def test_gray_group_keeps_only_time_even_sectors(self):
        self.assertEqual(
            magnetic_response_tensor_basis("11'", "normal_shift_current").dimension,
            18,
        )
        self.assertEqual(
            magnetic_response_tensor_basis("11'", "magnetic_injection_current").dimension,
            0,
        )
        self.assertEqual(
            magnetic_response_tensor_basis("11'", "shg_odd").dimension,
            0,
        )

    def test_pt_group_keeps_only_time_odd_polar_rank_three_sectors(self):
        self.assertEqual(
            magnetic_response_tensor_basis("-1'", "normal_shift_current").dimension,
            0,
        )
        self.assertEqual(
            magnetic_response_tensor_basis("-1'", "magnetic_injection_current").dimension,
            18,
        )
        self.assertEqual(
            magnetic_response_tensor_basis("-1'", "magnetic_shift_current").dimension,
            9,
        )

    def test_published_pt_shg_example_has_two_odd_components(self):
        even = magnetic_response_tensor_basis("-3'm'", "shg_even")
        odd = magnetic_response_tensor_basis("-3'm'", "shg_odd")
        self.assertEqual((even.dimension, odd.dimension), (0, 2))
        self.assertEqual(
            odd.basis,
            (
                (
                    (-1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
                (
                    (0.0, 0.0, 0.0, 0.0, 0.0, -1.0),
                    (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ),
        )

    def test_alias_validation_and_cli(self):
        self.assertEqual(
            magnetic_response_tensor_basis("-1'", "MIC").response,
            "magnetic_injection_current",
        )
        with self.assertRaises(GroupDataError):
            magnetic_response_tensor_basis("1", "unknown")
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(["magnetic-responses", "74", "shg_odd", "--json"]),
                0,
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["dimension"], 2)
        self.assertEqual(payload["time_character"], "odd")


if __name__ == "__main__":
    unittest.main()
