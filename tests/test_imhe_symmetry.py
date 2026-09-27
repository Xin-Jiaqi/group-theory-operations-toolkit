"""Regression tests for planar intrinsic/extrinsic IMHE tensor selection rules."""

from __future__ import annotations

import unittest

from group_theory_operations import (
    GroupDataError,
    iter_magnetic_layer_groups,
    iter_magnetic_point_groups,
    load_database,
    load_magnetic_layer_group_registry,
    load_magnetic_point_group_registry,
    magnetic_imhe_tensor_basis,
    magnetic_layer_imhe_tensor_basis,
)


class IMHESymmetryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.point_registry = load_magnetic_point_group_registry()
        cls.layer_registry = load_magnetic_layer_group_registry()

    def point(self, group, sector):
        return magnetic_imhe_tensor_basis(
            group,
            sector,
            database=self.database,
            registry=self.point_registry,
        )

    def test_time_reversal_and_inversion(self):
        self.assertEqual(self.point("11'", "intrinsic").dimension, 1)
        self.assertEqual(self.point("11'", "extrinsic").dimension, 0)
        self.assertEqual(self.point("-1", "intrinsic").dimension, 0)
        self.assertEqual(self.point("-1", "extrinsic").dimension, 0)
        self.assertEqual(self.point("m", "intrinsic").dimension, 0)
        self.assertEqual(self.point("m", "extrinsic").dimension, 0)

    def test_c3z(self):
        intrinsic = self.point("3", "intrinsic")
        extrinsic = self.point("3", "extrinsic")
        self.assertEqual(
            intrinsic.basis,
            (((0.0, -1.0), (1.0, 0.0)),),
        )
        self.assertEqual(
            extrinsic.basis,
            (((1.0, 0.0), (0.0, 1.0)),),
        )

    def test_gray_c3z(self):
        self.assertEqual(self.point("31'", "intrinsic").dimension, 1)
        self.assertEqual(self.point("31'", "extrinsic").dimension, 0)

    def test_c4z_time_reversal_extrinsic(self):
        result = self.point("4'", "extrinsic")
        self.assertEqual(result.dimension, 2)
        for matrix in result.basis:
            self.assertAlmostEqual(matrix[0][0] + matrix[1][1], 0.0)
            self.assertAlmostEqual(matrix[0][1], matrix[1][0])

    def test_pt_independent_prediction(self):
        # Das & Agarwal arXiv v1 Table I lists extrinsic PT as forbidden.
        # The stated T-odd pseudotensor transformation instead allows the
        # complete planar symmetric tensor, so the independent derivation is
        # deliberately frozen here for later comparison with the missing SM.
        self.assertEqual(self.point("-1'", "intrinsic").dimension, 0)
        self.assertEqual(self.point("-1'", "extrinsic").dimension, 3)

    def test_axis_preserving_magnetic_point_groups(self):
        applicable = 0
        rejected = 0
        for group in iter_magnetic_point_groups(self.point_registry):
            try:
                self.point(group.number, "intrinsic")
            except GroupDataError as exc:
                self.assertIn("mix the chosen z axis", str(exc))
                rejected += 1
            else:
                applicable += 1
        self.assertEqual((applicable, rejected), (106, 16))

    def test_all_magnetic_layer_groups(self):
        groups = tuple(iter_magnetic_layer_groups(self.layer_registry))
        self.assertEqual(len(groups), 528)
        for group in groups:
            for sector in ("intrinsic", "extrinsic"):
                result = magnetic_layer_imhe_tensor_basis(
                    group.global_number,
                    sector,
                    registry=self.layer_registry,
                )
                self.assertEqual(result.shape, (2, 2))
                for matrix in result.basis:
                    if sector == "intrinsic":
                        self.assertAlmostEqual(matrix[0][0], 0.0)
                        self.assertAlmostEqual(matrix[1][1], 0.0)
                        self.assertAlmostEqual(matrix[0][1], -matrix[1][0])
                    else:
                        self.assertAlmostEqual(matrix[0][1], matrix[1][0])


if __name__ == "__main__":
    unittest.main()
