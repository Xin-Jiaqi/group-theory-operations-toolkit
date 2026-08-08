"""Coordinate-contract tests for affine Seitz operations on structures."""

from __future__ import annotations

import unittest

import numpy as np

from group_theory_operations import GroupDataError
from group_theory_operations.seitz import (
    SeitzOp,
    equivalent,
    inverse,
    multiply,
    transform_seitz_coordinates,
)
from group_theory_operations.structure import _apply_seitz_operation


class StructureDouble:
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
        self.selective_dynamics = (
            None
            if selective_dynamics is None
            else tuple(tuple(row) for row in selective_dynamics)
        )
        self.length_unit = length_unit

    @classmethod
    def from_fractional(cls, **kwargs):
        return cls(**kwargs)

    def wrapped(self):
        coordinates = [
            tuple(
                value % 1.0 if self.pbc[axis] else value
                for axis, value in enumerate(row)
            )
            for row in self.fractional_coordinates
        ]
        return self.from_fractional(
            lattice=self.lattice,
            species=self.species,
            fractional_coordinates=coordinates,
            pbc=self.pbc,
            selective_dynamics=self.selective_dynamics,
            length_unit=self.length_unit,
        )


class SeitzCoordinateConventionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.left = SeitzOp(
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([0.25, 0.0, 0.0]),
        )
        self.right = SeitzOp(-np.eye(3, dtype=np.int64), np.array([0.0, 0.5, 0.0]))
        self.transformation = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        self.origin_shift = np.array([0.125, 0.25, 0.375])

    def test_transformed_operation_commutes_with_point_coordinate_change(self) -> None:
        operation = transform_seitz_coordinates(
            self.left, self.transformation, self.origin_shift
        )
        old_point = np.array([0.13, 0.27, 0.41])
        new_point = (self.transformation @ old_point + self.origin_shift) % 1.0
        expected = (
            self.transformation @ self.left.apply(old_point) + self.origin_shift
        ) % 1.0
        self.assertTrue(np.allclose(operation.apply(new_point), expected, atol=1.0e-12))

    def test_coordinate_change_preserves_products_and_inverses(self) -> None:
        transformed_product = transform_seitz_coordinates(
            multiply(self.left, self.right),
            self.transformation,
            self.origin_shift,
        )
        product_of_transformed = multiply(
            transform_seitz_coordinates(
                self.left, self.transformation, self.origin_shift
            ),
            transform_seitz_coordinates(
                self.right, self.transformation, self.origin_shift
            ),
        )
        self.assertTrue(equivalent(transformed_product, product_of_transformed))
        self.assertTrue(
            equivalent(
                transform_seitz_coordinates(
                    inverse(self.left), self.transformation, self.origin_shift
                ),
                inverse(
                    transform_seitz_coordinates(
                        self.left, self.transformation, self.origin_shift
                    )
                ),
            )
        )

    def test_origin_only_change_uses_declared_sign_convention(self) -> None:
        operation = SeitzOp(
            np.diag([-1, 1, -1]),
            np.array([0.0, 0.5, 0.0]),
        )
        transformed = transform_seitz_coordinates(
            operation,
            np.eye(3),
            np.array([0.25, 0.125, 0.5]),
        )
        self.assertTrue(np.allclose(transformed.translation, [0.5, 0.5, 0.0]))

    def test_invalid_coordinate_transformations_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "3x3"):
            transform_seitz_coordinates(self.left, np.eye(2), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "length-3"):
            transform_seitz_coordinates(self.left, np.eye(3), np.zeros(2))
        with self.assertRaisesRegex(ValueError, "invertible"):
            transform_seitz_coordinates(self.left, np.zeros((3, 3)), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "integer fractional rotation"):
            transform_seitz_coordinates(
                self.left,
                np.diag([2.0, 1.0, 1.0]),
                np.zeros(3),
            )


class AffineStructureOperationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.structure = StructureDouble(
            lattice=np.eye(3),
            species=["B", "N"],
            fractional_coordinates=[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            selective_dynamics=[[True, False, True], [True, True, True]],
        )
        self.left = SeitzOp(
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([0.25, 0.0, 0.0]),
        )
        self.right = SeitzOp(-np.eye(3, dtype=np.int64), np.array([0.0, 0.5, 0.0]))

    def test_affine_operation_applies_rotation_and_translation(self) -> None:
        result = _apply_seitz_operation(self.structure, self.left)
        self.assertTrue(
            np.allclose(result.fractional_coordinates[0], [0.05, 0.1, 0.3])
        )
        self.assertEqual(result.selective_dynamics[0], (False, True, True))
        self.assertEqual(result.species, self.structure.species)
        self.assertEqual(result.lattice, self.structure.lattice)

    def test_structure_action_follows_seitz_product_order(self) -> None:
        sequential = _apply_seitz_operation(
            _apply_seitz_operation(self.structure, self.right), self.left
        )
        direct = _apply_seitz_operation(
            self.structure, multiply(self.left, self.right)
        )
        self.assertTrue(
            np.allclose(
                sequential.fractional_coordinates,
                direct.fractional_coordinates,
                atol=1.0e-12,
            )
        )

    def test_inverse_restores_wrapped_structure_coordinates(self) -> None:
        transformed = _apply_seitz_operation(self.structure, self.left)
        restored = _apply_seitz_operation(transformed, inverse(self.left))
        self.assertTrue(
            np.allclose(
                restored.fractional_coordinates,
                self.structure.fractional_coordinates,
                atol=1.0e-12,
            )
        )

    def test_nonorthogonal_hexagonal_screw_operation(self) -> None:
        hexagonal = StructureDouble(
            lattice=[[2.0, 0.0, 0.0], [-1.0, np.sqrt(3.0), 0.0], [0.0, 0.0, 7.0]],
            species=["H"],
            fractional_coordinates=[[0.25, 0.5, 0.1]],
        )
        screw = SeitzOp(
            np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]]),
            np.array([0.0, 0.0, 0.5]),
        )
        result = _apply_seitz_operation(hexagonal, screw)
        self.assertTrue(
            np.allclose(result.fractional_coordinates[0], [0.5, 0.75, 0.6])
        )

    def test_wrap_false_keeps_unreduced_affine_coordinates(self) -> None:
        result = _apply_seitz_operation(self.structure, self.right, wrap=False)
        self.assertTrue(
            np.allclose(result.fractional_coordinates[0], [-0.1, 0.3, -0.3])
        )

    def test_structure_semantic_guards_apply_to_seitz_operations(self) -> None:
        anisotropic = StructureDouble(
            lattice=np.diag([2.0, 3.0, 4.0]),
            species=["H"],
            fractional_coordinates=[[0.1, 0.2, 0.3]],
        )
        with self.assertRaisesRegex(GroupDataError, "lattice metric"):
            _apply_seitz_operation(anisotropic, self.left)

        layer = StructureDouble(
            lattice=np.eye(3),
            species=["H"],
            fractional_coordinates=[[0.1, 0.2, 0.3]],
            pbc=(True, True, False),
        )
        cyclic_axes = SeitzOp(
            np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
            np.zeros(3),
        )
        with self.assertRaisesRegex(GroupDataError, "PBC axis semantics"):
            _apply_seitz_operation(layer, cyclic_axes)


if __name__ == "__main__":
    unittest.main()
