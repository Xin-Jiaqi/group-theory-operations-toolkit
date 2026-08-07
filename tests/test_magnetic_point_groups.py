"""Verification of the 122 magnetic point groups and time-parity solver."""

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import redirect_stdout
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from jsonschema import Draft202012Validator
import spglib

from group_theory_operations import (
    GroupDataError,
    get_magnetic_point_group,
    iter_magnetic_point_groups,
    load_database,
    load_magnetic_point_group_registry,
    magnetic_point_group_operations,
    magnetic_tensor_basis,
)
from group_theory_operations.cli import main


ROOT = Path(__file__).resolve().parents[1]

# Frozen transcription of the standard traditional-symbol table.  Each row is
# ordered type I, type II, then the type-III groups for one parent point group.
EXPECTED_BY_PARENT = {
    1: ("1", "11'"),
    2: ("-1", "-11'", "-1'"),
    3: ("2", "21'", "2'"),
    4: ("m", "m1'", "m'"),
    5: ("2/m", "2/m1'", "2'/m", "2/m'", "2'/m'"),
    6: ("222", "2221'", "2'2'2"),
    7: ("mm2", "mm21'", "m'm2'", "m'm'2"),
    8: ("mmm", "mmm1'", "m'mm", "m'm'm", "m'm'm'"),
    9: ("4", "41'", "4'"),
    10: ("-4", "-41'", "-4'"),
    11: ("4/m", "4/m1'", "4'/m", "4/m'", "4'/m'"),
    12: ("422", "4221'", "4'22'", "42'2'"),
    13: ("4mm", "4mm1'", "4'm'm", "4m'm'"),
    14: ("-42m", "-42m1'", "-4'2'm", "-4'2m'", "-42'm'"),
    15: (
        "4/mmm", "4/mmm1'", "4/m'mm", "4'/mm'm", "4'/m'm'm",
        "4/mm'm'", "4/m'm'm'",
    ),
    16: ("3", "31'"),
    17: ("-3", "-31'", "-3'"),
    18: ("32", "321'", "32'"),
    19: ("3m", "3m1'", "3m'"),
    20: ("-3m", "-3m1'", "-3'm", "-3'm'", "-3m'"),
    21: ("6", "61'", "6'"),
    22: ("-6", "-61'", "-6'"),
    23: ("6/m", "6/m1'", "6'/m", "6/m'", "6'/m'"),
    24: ("622", "6221'", "6'22'", "62'2'"),
    25: ("6mm", "6mm1'", "6'mm'", "6m'm'"),
    26: ("-6m2", "-6m21'", "-6'm'2", "-6'm2'", "-6m'2'"),
    27: (
        "6/mmm", "6/mmm1'", "6/m'mm", "6'/mmm'", "6'/m'mm'",
        "6/mm'm'", "6/m'm'm'",
    ),
    28: ("23", "231'"),
    29: ("m-3", "m-31'", "m'-3'"),
    30: ("432", "4321'", "4'32'"),
    31: ("-43m", "-43m1'", "-4'3m'"),
    32: ("m-3m", "m-3m1'", "m'-3'm", "m-3m'", "m'-3'm'"),
}


def matrix_key(matrix, digits: int = 8):
    return tuple(round(float(value), digits) for row in matrix for value in row)


def matmul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def rotation_signature(matrix):
    determinant = round(
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    trace = round(sum(matrix[index][index] for index in range(3)))
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    power = identity
    for order in range(1, 13):
        power = matmul(power, matrix)
        if matrix_key(power) == matrix_key(identity):
            return determinant, trace, order
    raise AssertionError("rotation order exceeds 12")


def colored_signature(pairs):
    return tuple(
        sorted(
            Counter(
                (rotation_signature(matrix), int(time_reversal))
                for matrix, time_reversal in pairs
            ).items()
        )
    )


class MagneticPointGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.registry = load_magnetic_point_group_registry()
        cls.groups = tuple(iter_magnetic_point_groups(cls.registry))

    def test_schema_hashes_numbers_symbols_and_category_counts(self):
        schema = json.loads(
            (ROOT / "schema" / "magnetic-point-groups-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        errors = sorted(Draft202012Validator(schema).iter_errors(self.registry), key=str)
        self.assertEqual(errors, [])
        self.assertEqual(
            self.registry["point_group_registry_sha256"],
            hashlib.sha256((ROOT / "data" / "crystallographic_point_groups.json").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.registry["operation_catalog_sha256"],
            hashlib.sha256((ROOT / "data" / "group_operations.json").read_bytes()).hexdigest(),
        )
        expected_symbols = [symbol for values in EXPECTED_BY_PARENT.values() for symbol in values]
        self.assertEqual([group.number for group in self.groups], list(range(1, 123)))
        self.assertEqual([group.hm_symbol for group in self.groups], expected_symbols)
        self.assertEqual(
            Counter(group.category for group in self.groups),
            {"type_I": 32, "type_II_gray": 32, "type_III_black_white": 58},
        )
        for group in self.groups:
            parent, local, global_number = map(int, group.magnetic_number.split("."))
            self.assertEqual((parent, global_number), (group.parent_point_group_number, group.number))
            self.assertEqual(
                local,
                EXPECTED_BY_PARENT[parent].index(group.hm_symbol) + 1,
            )

    def test_every_colored_operation_set_is_closed_and_generated(self):
        for group in self.groups:
            operations = magnetic_point_group_operations(
                group.number, database=self.database, registry=self.registry
            )
            by_pair = {(item.name, item.time_reversal) for item in operations}
            spatial_by_matrix = {
                matrix_key(item.spatial.matrix_cartesian): item.name for item in operations
            }
            spatial_by_name = {item.name: item.spatial.matrix_cartesian for item in operations}
            multiplication = {
                (left, right): spatial_by_matrix[
                    matrix_key(matmul(spatial_by_name[left], spatial_by_name[right]))
                ]
                for left in spatial_by_name
                for right in spatial_by_name
            }
            self.assertEqual(len(by_pair), group.order, group.hm_symbol)
            for left in operations:
                for right in operations:
                    product = (
                        multiplication[left.name, right.name],
                        left.time_reversal ^ right.time_reversal,
                    )
                    self.assertIn(product, by_pair, group.hm_symbol)
            generated = {
                ("1", False),
                *((generator.name, generator.time_reversal) for generator in group.generators),
            }
            while True:
                expanded = set(generated)
                for left_name, left_time in generated:
                    for right_name, right_time in generated:
                        expanded.add(
                            (
                                multiplication[left_name, right_name],
                                left_time ^ right_time,
                            )
                        )
                if expanded == generated:
                    break
                generated = expanded
            self.assertEqual(generated, set(by_pair), group.hm_symbol)

    def test_all_122_types_match_independent_spglib_magnetic_database(self):
        internal = defaultdict(set)
        category_names = {
            "type_I": "ordinary",
            "type_II_gray": "gray",
            "type_III_black_white": "black_white",
        }
        for group in self.groups:
            operations = magnetic_point_group_operations(
                group.number, database=self.database, registry=self.registry
            )
            pairs = [
                (item.spatial.matrix_cartesian, item.time_reversal) for item in operations
            ]
            internal[(group.parent_point_group_hm, category_names[group.category])].add(
                colored_signature(pairs)
            )

        space_group_point_groups = {}
        for hall_number in range(1, 531):
            item = spglib.get_spacegroup_type(hall_number)
            space_group_point_groups.setdefault(item.number, item.pointgroup_international)
        external = defaultdict(set)
        identity_key = matrix_key(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        for uni_number in range(1, 1652):
            item = spglib.get_magnetic_spacegroup_type(uni_number)
            symmetry = spglib.get_magnetic_symmetry_from_database(uni_number)
            pairs_by_key = {}
            for rotation, time_reversal in zip(
                symmetry["rotations"], symmetry["time_reversals"], strict=True
            ):
                matrix = tuple(tuple(float(value) for value in row) for row in rotation)
                pairs_by_key[(matrix_key(matrix), bool(time_reversal))] = matrix
            pure_time_reversal = (identity_key, True) in pairs_by_key
            if pure_time_reversal:
                category = "gray"
            elif all(not key[1] for key in pairs_by_key):
                category = "ordinary"
            else:
                category = "black_white"
            external[(space_group_point_groups[item.number], category)].add(
                colored_signature(
                    (matrix, key[1]) for key, matrix in pairs_by_key.items()
                )
            )
        self.assertEqual(internal, external)

    def test_lookup_and_loader_rejection(self):
        self.assertEqual(get_magnetic_point_group(104).hm_symbol, "6'/m'mm'")
        self.assertEqual(get_magnetic_point_group("27.5.104").number, 104)
        self.assertEqual(get_magnetic_point_group("6′/m′mm′").number, 104)
        with self.assertRaises(GroupDataError):
            get_magnetic_point_group("not-a-group")
        malformed = json.loads(json.dumps(self.registry))
        malformed["magnetic_point_groups"][0]["operations"][0]["time_reversal"] = 0
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registry.json"
            path.write_text(json.dumps(malformed), encoding="utf-8")
            with self.assertRaises(GroupDataError):
                load_magnetic_point_group_registry(path)

    def test_cli_list_and_show(self):
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(["magnetic-point-groups", "27.5.104", "--json"]), 0
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["hm_symbol"], "6'/m'mm'")
        self.assertEqual(len(payload["operation_records"]), 24)


class MagneticTensorTests(unittest.TestCase):
    def test_gray_group_enforces_time_parity(self):
        forbidden = magnetic_tensor_basis(
            "11'",
            "polar_vector",
            "polar_vector",
            output_time_parity="odd",
            input_time_parity="even",
        )
        allowed = magnetic_tensor_basis(
            "11'",
            "polar_vector",
            "polar_vector",
            output_time_parity="odd",
            input_time_parity="odd",
        )
        self.assertEqual(forbidden.dimension, 0)
        self.assertEqual(allowed.dimension, 9)

    def test_type_I_is_independent_of_declared_time_parity(self):
        even = magnetic_tensor_basis("2", "axial", "scalar")
        odd = magnetic_tensor_basis(
            "2", "axial", "scalar", output_time_parity="odd"
        )
        self.assertEqual(even.basis, odd.basis)
        self.assertEqual(even.dimension, 1)

    def test_black_white_group_changes_allowed_axial_direction(self):
        basis = magnetic_tensor_basis(
            "2'",
            "axial_vector",
            "scalar",
            output_time_parity="odd",
            input_time_parity="even",
        )
        self.assertEqual(basis.dimension, 2)
        self.assertEqual(basis.shape, (3, 1))

    def test_real_tensor_spaces_and_cli(self):
        spaces = (
            "scalar",
            "pseudoscalar",
            "polar_vector",
            "axial_vector",
            "symmetric_quadratic",
            "antisymmetric_quadratic",
        )
        for space in spaces:
            result = magnetic_tensor_basis("1", space, "scalar")
            self.assertGreaterEqual(result.dimension, 1)
            self.assertTrue(all(math.isfinite(value) for matrix in result.basis for row in matrix for value in row))
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(
                    [
                        "magnetic-invariants", "2'", "axial_vector", "scalar",
                        "--output-time", "odd", "--json",
                    ]
                ),
                0,
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["dimension"], 2)
        self.assertEqual(payload["shape"], [3, 1])


if __name__ == "__main__":
    unittest.main()
