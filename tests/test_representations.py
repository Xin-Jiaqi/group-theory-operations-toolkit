from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import unittest

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from group_theory_operations import (  # noqa: E402
    antisymmetric_field_matrix,
    determinant3,
    iter_operations,
    load_database,
    load_quadratic_field_catalog,
    quadratic_field_representation,
    symmetric_field_matrix,
)
from group_theory_operations.cli import main  # noqa: E402


def matmul(left, right):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def matvec(matrix, vector):
    return tuple(sum(row[index] * vector[index] for index in range(len(vector))) for row in matrix)


def transpose(matrix):
    return tuple(tuple(row[index] for row in matrix) for index in range(len(matrix[0])))


def matrix_key(matrix, digits=9):
    return tuple(round(float(value), digits) for row in matrix for value in row)


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def symmetric_components(field):
    x, y, z = field
    return (
        abs(x) ** 2,
        abs(y) ** 2,
        abs(z) ** 2,
        x * y.conjugate() + y * x.conjugate(),
        x * z.conjugate() + z * x.conjugate(),
        y * z.conjugate() + z * y.conjugate(),
    )


def antisymmetric_components(field):
    conjugate = tuple(value.conjugate() for value in field)
    return tuple(1j * value for value in cross(field, conjugate))


class QuadraticFieldRepresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.catalog = load_quadratic_field_catalog()

    def assertVectorAlmostEqual(self, left, right, places=9):
        self.assertEqual(len(left), len(right))
        for actual, expected in zip(left, right, strict=True):
            self.assertAlmostEqual(actual.real, expected.real, places=places)
            self.assertAlmostEqual(actual.imag, expected.imag, places=places)

    def test_all_operation_artifact_matches_public_api(self):
        self.assertEqual(
            self.catalog["source_catalog_sha256"],
            hashlib.sha256((ROOT / "data" / "group_operations.json").read_bytes()).hexdigest(),
        )
        total = 0
        for family_name, family in self.database["families"].items():
            generated = self.catalog["families"][family_name]
            self.assertEqual(generated["operation_count"], len(family["operations"]))
            self.assertEqual(
                [item["name"] for item in generated["operations"]],
                [item["name"] for item in family["operations"]],
            )
            for operation, item in zip(family["operations"], generated["operations"], strict=True):
                representation = quadratic_field_representation(operation)
                self.assertEqual(item["determinant"], representation.determinant)
                self.assertEqual(
                    matrix_key(item["matrix_symmetric"]),
                    matrix_key(representation.matrix_symmetric),
                )
                self.assertEqual(
                    matrix_key(item["matrix_antisymmetric"]),
                    matrix_key(representation.matrix_antisymmetric),
                )
                total += 1
        self.assertEqual(total, 88)

    def test_generated_catalog_satisfies_schema(self):
        schema = json.loads(
            (ROOT / "schema" / "quadratic-field-representations-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        errors = sorted(Draft202012Validator(schema).iter_errors(self.catalog), key=str)
        self.assertEqual(errors, [])

    def test_defining_field_actions_for_every_operation(self):
        field = (0.31 + 0.27j, -0.43 + 0.19j, 0.23 - 0.37j)
        symmetric = symmetric_components(field)
        antisymmetric = antisymmetric_components(field)
        for family_name, operation in iter_operations(self.database):
            matrix = operation["matrix_cartesian"]
            transformed_field = matvec(matrix, field)
            representation = quadratic_field_representation(operation)
            self.assertVectorAlmostEqual(
                matvec(representation.matrix_symmetric, symmetric),
                symmetric_components(transformed_field),
            )
            self.assertVectorAlmostEqual(
                matvec(representation.matrix_antisymmetric, antisymmetric),
                antisymmetric_components(transformed_field),
            )

    def test_cross_product_identity_for_every_operation(self):
        left = (0.31, -0.43, 0.23)
        right = (-0.17, 0.29, 0.47)
        source_cross = cross(left, right)
        for _, operation in iter_operations(self.database):
            matrix = operation["matrix_cartesian"]
            self.assertVectorAlmostEqual(
                cross(matvec(matrix, left), matvec(matrix, right)),
                matvec(antisymmetric_field_matrix(matrix), source_cross),
            )

    def test_representations_are_homomorphisms(self):
        for family_name, family in self.database["families"].items():
            operations = family["operations"]
            for left in operations:
                left_representation = quadratic_field_representation(left)
                for right in operations:
                    product = matmul(left["matrix_cartesian"], right["matrix_cartesian"])
                    self.assertEqual(
                        matrix_key(symmetric_field_matrix(product)),
                        matrix_key(
                            matmul(
                                left_representation.matrix_symmetric,
                                quadratic_field_representation(right).matrix_symmetric,
                            )
                        ),
                        family_name,
                    )
                    self.assertEqual(
                        matrix_key(antisymmetric_field_matrix(product)),
                        matrix_key(
                            matmul(
                                left_representation.matrix_antisymmetric,
                                quadratic_field_representation(right).matrix_antisymmetric,
                            )
                        ),
                        family_name,
                    )

    def test_antisymmetric_matrices_are_proper_orthogonal(self):
        identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        for _, operation in iter_operations(self.database):
            matrix = antisymmetric_field_matrix(operation["matrix_cartesian"])
            self.assertEqual(matrix_key(matmul(transpose(matrix), matrix)), matrix_key(identity))
            self.assertAlmostEqual(determinant3(matrix), 1.0, places=9)

    def test_reference_examples_and_inversion(self):
        expected = {
            "1": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            "-1": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            "2_001": ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
            "m_100": ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
        }
        operations = {
            operation["name"]: operation
            for operation in self.database["families"]["tetragonal_D4h"]["operations"]
        }
        for name, matrix in expected.items():
            self.assertEqual(
                matrix_key(antisymmetric_field_matrix(operations[name]["matrix_cartesian"])),
                matrix_key(matrix),
            )

    def test_symmetric_normalizations_are_similarity_related(self):
        scale = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
        scale_matrix = tuple(
            tuple(scale[row] if row == column else 0.0 for column in range(6))
            for row in range(6)
        )
        inverse_scale = tuple(
            tuple((1.0 / scale[row]) if row == column else 0.0 for column in range(6))
            for row in range(6)
        )
        for _, operation in iter_operations(self.database):
            matrix = operation["matrix_cartesian"]
            converted = matmul(
                matmul(scale_matrix, symmetric_field_matrix(matrix, normalization="half")),
                inverse_scale,
            )
            self.assertEqual(
                matrix_key(converted),
                matrix_key(symmetric_field_matrix(matrix, normalization="sum")),
            )

    def test_cli_returns_machine_readable_antisymmetric_matrix(self):
        output = io.StringIO()
        with redirect_stdout(output):
            code = main(
                [
                    "field-representation",
                    "m_100",
                    "--family",
                    "tetragonal_D4h",
                    "--space",
                    "antisymmetric",
                    "--json",
                ]
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(code, 0)
        self.assertEqual(payload[0]["matrix_antisymmetric"], [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        self.assertNotIn("matrix_symmetric", payload[0])

    def test_generated_markdown_covers_all_operations(self):
        document = (ROOT / "docs" / "quadratic_field_representations.md").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("<details", document)
        for family_name, family in self.catalog["families"].items():
            self.assertIn(f"### `{family_name}`：{family['operation_count']} 个操作", document)
            for operation in family["operations"]:
                self.assertIn(f"| {operation['index']} | `{operation['name']}` |", document)

    def test_invalid_input_is_rejected(self):
        for matrix in ([[1, 0], [0, 1]], [[True, 0, 0], [0, 1, 0], [0, 0, 1]]):
            with self.assertRaises(ValueError):
                quadratic_field_representation(matrix)
        with self.assertRaises(ValueError):
            symmetric_field_matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)), normalization="bad")


if __name__ == "__main__":
    unittest.main()
