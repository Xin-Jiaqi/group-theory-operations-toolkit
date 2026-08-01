from __future__ import annotations

from contextlib import redirect_stdout
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

from jsonschema import Draft202012Validator
import spglib


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from group_theory_operations import (  # noqa: E402
    GroupDataError,
    equivariant_map_basis,
    get_crystallographic_point_group,
    iter_crystallographic_point_groups,
    load_database,
    load_optical_response_catalog,
    load_point_group_registry,
    point_group_operations,
    quadratic_field_representation,
    response_tensor_basis,
)
from group_theory_operations.cli import main  # noqa: E402


EXPECTED_ORDERS = {
    "1": 1,
    "-1": 2,
    "2": 2,
    "m": 2,
    "2/m": 4,
    "222": 4,
    "mm2": 4,
    "mmm": 8,
    "4": 4,
    "-4": 4,
    "4/m": 8,
    "422": 8,
    "4mm": 8,
    "-42m": 8,
    "4/mmm": 16,
    "3": 3,
    "-3": 6,
    "32": 6,
    "3m": 6,
    "-3m": 12,
    "6": 6,
    "-6": 6,
    "6/m": 12,
    "622": 12,
    "6mm": 12,
    "-6m2": 12,
    "6/mmm": 24,
    "23": 12,
    "m-3": 24,
    "432": 24,
    "-43m": 24,
    "m-3m": 48,
}

EXPECTED_DIMENSIONS = {
    "1": (18, 9),
    "-1": (0, 0),
    "2": (8, 5),
    "m": (10, 4),
    "2/m": (0, 0),
    "222": (3, 3),
    "mm2": (5, 2),
    "mmm": (0, 0),
    "4": (4, 3),
    "-4": (4, 2),
    "4/m": (0, 0),
    "422": (1, 2),
    "4mm": (3, 1),
    "-42m": (2, 1),
    "4/mmm": (0, 0),
    "3": (6, 3),
    "-3": (0, 0),
    "32": (2, 2),
    "3m": (4, 1),
    "-3m": (0, 0),
    "6": (4, 3),
    "-6": (2, 0),
    "6/m": (0, 0),
    "622": (1, 2),
    "6mm": (3, 1),
    "-6m2": (1, 0),
    "6/mmm": (0, 0),
    "23": (1, 1),
    "m-3": (0, 0),
    "432": (0, 1),
    "-43m": (1, 0),
    "m-3m": (0, 0),
}


def matmul(left, right):
    return tuple(
        tuple(
            sum(left[row][index] * right[index][column] for index in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def matrix_key(matrix, digits=8):
    return tuple(round(float(value), digits) for row in matrix for value in row)


def matrix_signature(matrix):
    determinant = round(
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    trace = round(sum(matrix[index][index] for index in range(3)))
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    power = identity
    order = None
    for candidate in range(1, 13):
        power = matmul(power, matrix)
        if matrix_key(power) == matrix_key(identity):
            order = candidate
            break
    if order is None:
        raise AssertionError("crystallographic rotation order exceeds 12")
    return determinant, trace, order


class CrystallographicPointGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.registry = load_point_group_registry()
        cls.groups = tuple(iter_crystallographic_point_groups(cls.registry))

    def test_registry_schema_hash_order_and_metadata(self):
        schema = json.loads(
            (ROOT / "schema" / "crystallographic-point-groups-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        errors = sorted(Draft202012Validator(schema).iter_errors(self.registry), key=str)
        self.assertEqual(errors, [])
        self.assertEqual(
            self.registry["source_catalog_sha256"],
            hashlib.sha256((ROOT / "data" / "group_operations.json").read_bytes()).hexdigest(),
        )
        self.assertEqual([group.number for group in self.groups], list(range(1, 33)))
        self.assertEqual({group.hm_symbol: group.order for group in self.groups}, EXPECTED_ORDERS)

    def test_every_embedding_is_a_generated_closed_group(self):
        for group in self.groups:
            operations = point_group_operations(
                group.number,
                database=self.database,
                registry=self.registry,
            )
            by_matrix = {matrix_key(operation.matrix_cartesian): operation.name for operation in operations}
            self.assertEqual(len(by_matrix), group.order, group.hm_symbol)
            self.assertIn("1", group.operations)
            self.assertTrue(set(group.generators).issubset(group.operations))
            for left in operations:
                for right in operations:
                    product = matrix_key(matmul(left.matrix_cartesian, right.matrix_cartesian))
                    self.assertIn(product, by_matrix, group.hm_symbol)
            generated = {"1", *group.generators}
            while True:
                expanded = set(generated)
                selected = [operation for operation in operations if operation.name in generated]
                for left in selected:
                    for right in selected:
                        expanded.add(by_matrix[matrix_key(matmul(left.matrix_cartesian, right.matrix_cartesian))])
                if expanded == generated:
                    break
                generated = expanded
            self.assertEqual(generated, set(group.operations), group.hm_symbol)

    def test_lookup_accepts_number_hm_schoenflies_and_alias(self):
        self.assertEqual(get_crystallographic_point_group(15).hm_symbol, "4/mmm")
        self.assertEqual(get_crystallographic_point_group("4/mmm").schoenflies_symbol, "D4h")
        self.assertEqual(get_crystallographic_point_group("d4H").number, 15)
        self.assertEqual(get_crystallographic_point_group("S6").hm_symbol, "-3")
        self.assertEqual(get_crystallographic_point_group("S2").hm_symbol, "-1")
        with self.assertRaises(GroupDataError):
            get_crystallographic_point_group("not-a-group")

    def test_loader_rejects_non_boolean_flags(self):
        malformed = json.loads(json.dumps(self.registry))
        malformed["point_groups"][0]["polar"] = "true"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registry.json"
            path.write_text(json.dumps(malformed), encoding="utf-8")
            with self.assertRaises(GroupDataError):
                load_point_group_registry(path)

    def test_polar_chiral_and_centrosymmetric_flags(self):
        polar = {group.hm_symbol for group in self.groups if group.polar}
        chiral = {group.hm_symbol for group in self.groups if group.chiral}
        centrosymmetric = {group.hm_symbol for group in self.groups if group.centrosymmetric}
        self.assertEqual(polar, {"1", "2", "m", "mm2", "4", "4mm", "3", "3m", "6", "6mm"})
        self.assertEqual(chiral, {"1", "2", "222", "4", "422", "3", "32", "6", "622", "23", "432"})
        self.assertEqual(len(centrosymmetric), 11)

    def test_readme_lists_all_registered_groups(self):
        document = (ROOT / "README.md").read_text(encoding="utf-8")
        for group in self.groups:
            self.assertIn(
                f"| {group.number} | `{group.hm_symbol}` | `{group.schoenflies_symbol}` |",
                document,
            )

    def test_independent_spglib_point_group_signatures(self):
        spglib_rotations = {}
        for hall_number in range(1, 531):
            symbol = spglib.get_spacegroup_type(hall_number).pointgroup_international
            if symbol in spglib_rotations:
                continue
            symmetry = spglib.get_symmetry_from_database(hall_number)
            unique = {
                tuple(int(value) for row in rotation for value in row): rotation
                for rotation in symmetry["rotations"]
            }
            spglib_rotations[symbol] = tuple(unique.values())
        self.assertEqual(set(spglib_rotations), set(EXPECTED_ORDERS))
        for group in self.groups:
            internal = point_group_operations(
                group.number,
                database=self.database,
                registry=self.registry,
            )
            self.assertEqual(
                Counter(matrix_signature(operation.matrix_cartesian) for operation in internal),
                Counter(matrix_signature(rotation) for rotation in spglib_rotations[group.hm_symbol]),
                group.hm_symbol,
            )


class OpticalInvariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.registry = load_point_group_registry()
        cls.catalog = load_optical_response_catalog()

    def test_catalog_schema_hashes_and_public_solver_match(self):
        schema = json.loads(
            (ROOT / "schema" / "optical-response-invariants-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        errors = sorted(Draft202012Validator(schema).iter_errors(self.catalog), key=str)
        self.assertEqual(errors, [])
        self.assertEqual(
            self.catalog["point_group_registry_sha256"],
            hashlib.sha256((ROOT / "data" / "crystallographic_point_groups.json").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.catalog["operation_catalog_sha256"],
            hashlib.sha256((ROOT / "data" / "group_operations.json").read_bytes()).hexdigest(),
        )
        for item in self.catalog["point_groups"]:
            for response_name, stored in item["responses"].items():
                solved = response_tensor_basis(
                    item["number"],
                    response_name,
                    database=self.database,
                    registry=self.registry,
                )
                self.assertEqual(stored, solved.to_dict(), (item["hm_symbol"], response_name))

    def test_frozen_spatial_selection_rule_dimensions(self):
        for group in iter_crystallographic_point_groups(self.registry):
            symmetric, circular = EXPECTED_DIMENSIONS[group.hm_symbol]
            shift = response_tensor_basis(group.number, "shift", database=self.database, registry=self.registry)
            shg = response_tensor_basis(group.number, "shg", database=self.database, registry=self.registry)
            cpge = response_tensor_basis(group.number, "cpge", database=self.database, registry=self.registry)
            self.assertEqual((shift.dimension, cpge.dimension), (symmetric, circular), group.hm_symbol)
            self.assertEqual(shift.basis, shg.basis, group.hm_symbol)

    def test_every_basis_vector_is_equivariant_under_every_operation(self):
        for item in self.catalog["point_groups"]:
            operations = point_group_operations(
                item["number"],
                database=self.database,
                registry=self.registry,
            )
            for response_name in ("shift_current", "circular_injection_current"):
                solved = response_tensor_basis(
                    item["number"],
                    response_name,
                    database=self.database,
                    registry=self.registry,
                )
                for operation in operations:
                    fields = quadratic_field_representation(operation.matrix_cartesian)
                    source = (
                        fields.matrix_symmetric
                        if response_name == "shift_current"
                        else fields.matrix_antisymmetric
                    )
                    for basis in solved.basis:
                        self.assertEqual(
                            matrix_key(matmul(operation.matrix_cartesian, basis)),
                            matrix_key(matmul(basis, source)),
                            (item["hm_symbol"], response_name, operation.name),
                        )

    def test_dimensions_equal_independent_character_inner_products(self):
        for group in iter_crystallographic_point_groups(self.registry):
            operations = point_group_operations(
                group.number,
                database=self.database,
                registry=self.registry,
            )
            for response_name in ("shift_current", "circular_injection_current"):
                character_sum = 0.0
                for operation in operations:
                    fields = quadratic_field_representation(operation.matrix_cartesian)
                    source = (
                        fields.matrix_symmetric
                        if response_name == "shift_current"
                        else fields.matrix_antisymmetric
                    )
                    character_sum += sum(operation.matrix_cartesian[i][i] for i in range(3)) * sum(
                        source[i][i] for i in range(len(source))
                    )
                expected = round(character_sum / group.order)
                solved = response_tensor_basis(
                    group.number,
                    response_name,
                    database=self.database,
                    registry=self.registry,
                )
                self.assertAlmostEqual(character_sum / group.order, expected, places=8)
                self.assertEqual(solved.dimension, expected, (group.hm_symbol, response_name))

    def test_generic_solver_and_invalid_inputs(self):
        identity2 = ((1, 0), (0, 1))
        basis = equivariant_map_basis((identity2,), (identity2,))
        self.assertEqual(len(basis), 4)
        with self.assertRaises(ValueError):
            equivariant_map_basis((), ())
        with self.assertRaises(ValueError):
            equivariant_map_basis((identity2,), (identity2,), tolerance=0)
        with self.assertRaises(ValueError):
            equivariant_map_basis((identity2,), (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (identity2)))

    def test_readme_lists_all_dimensions(self):
        document = (ROOT / "README.md").read_text(encoding="utf-8")
        for item in self.catalog["point_groups"]:
            responses = item["responses"]
            self.assertIn(
                f"| {item['number']} | `{item['hm_symbol']}` | `{item['schoenflies_symbol']}` |",
                document,
            )
            self.assertEqual(
                responses["shift_current"]["dimension"],
                responses["shg"]["dimension"],
            )

    def test_cli_lists_registry_and_returns_json_basis(self):
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["point-groups", "4mm", "--json"]), 0)
        point_group = json.loads(output.getvalue())
        self.assertEqual(point_group["number"], 13)
        self.assertEqual(len(point_group["operation_records"]), 8)
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["invariants", "4mm", "shift_current", "--json"]), 0)
        response = json.loads(output.getvalue())
        self.assertEqual(response["shape"], [3, 6])
        self.assertEqual(response["dimension"], 3)


if __name__ == "__main__":
    unittest.main()
