"""Verification of all 528 magnetic layer-group point co-groups."""

from __future__ import annotations

from collections import Counter
from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import unittest

from jsonschema import Draft202012Validator
import spglib

from group_theory_operations import (
    GroupDataError,
    get_magnetic_layer_group,
    iter_magnetic_layer_groups,
    load_database,
    load_magnetic_layer_group_registry,
    magnetic_layer_response_tensor_basis,
    magnetic_layer_tensor_basis,
    screen_response_symmetry,
)
from group_theory_operations.cli import main


ROOT = Path(__file__).resolve().parents[1]
TYPE_TO_SPGLIB = {"I": 1, "II": 2, "III": 3, "IV": 4}


def matrix_key(matrix):
    return tuple(round(float(value), 8) for row in matrix for value in row)


def matmul(left, right):
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(3))
            for column in range(3)
        )
        for row in range(3)
    )


def rotation_signature(matrix):
    determinant = round(
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
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
    return Counter(
        (rotation_signature(matrix), bool(time_reversal))
        for matrix, time_reversal in pairs
    )


class MagneticLayerGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = load_database()
        cls.registry = load_magnetic_layer_group_registry()
        cls.groups = tuple(iter_magnetic_layer_groups(cls.registry))

    def test_schema_hashes_numbering_and_type_counts(self):
        schema = json.loads(
            (ROOT / "schema" / "magnetic-layer-groups-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        errors = sorted(
            Draft202012Validator(schema).iter_errors(self.registry), key=str
        )
        self.assertEqual(errors, [])
        self.assertEqual(
            self.registry["source_table_sha256"],
            hashlib.sha256(
                (
                    ROOT
                    / "scripts"
                    / "sources"
                    / "magnetic-layer-groups-528-v2023.tsv"
                ).read_bytes()
            ).hexdigest(),
        )
        self.assertEqual(
            self.registry["point_operation_catalog_sha256"],
            hashlib.sha256((ROOT / "data" / "group_operations.json").read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.registry["layer_group_registry_sha256"],
            hashlib.sha256(
                (ROOT / "data" / "crystallographic_layer_groups.json").read_bytes()
            ).hexdigest(),
        )
        self.assertEqual(
            [group.global_number for group in self.groups], list(range(1, 529))
        )
        self.assertEqual(
            Counter(group.magnetic_type for group in self.groups),
            {"I": 80, "II": 80, "III": 246, "IV": 122},
        )

    def test_source_table_round_trips_into_registry(self):
        with (
            ROOT / "scripts" / "sources" / "magnetic-layer-groups-528-v2023.tsv"
        ).open(encoding="utf-8", newline="") as handle:
            source = {row["og_number"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(len(source), 528)
        for group in self.groups:
            row = source[group.og_number]
            self.assertEqual(row["magnetic_type"], group.magnetic_type)
            self.assertEqual(row["litvin_og_symbol_ascii"], group.litvin_og_symbol_ascii)
            self.assertEqual(
                row["corresponding_msg_bns_number"],
                group.corresponding_magnetic_space_group.bns_number,
            )
            self.assertEqual(
                int(row["corresponding_msg_uni_number"]),
                group.corresponding_magnetic_space_group.uni_number,
            )

    def test_point_operation_sets_match_parent_layer_groups_and_are_closed(self):
        parent_operations = {}
        for family_name, family in self.database["families"].items():
            by_index = {item["index"]: item["name"] for item in family["operations"]}
            for layer_group in family["layer_groups"]:
                indices = set(layer_group["R+_indices"] + layer_group["R-_indices"])
                parent_operations[layer_group["LG"]] = {
                    by_index[index] for index in indices
                }
        self.assertEqual(set(parent_operations), set(range(1, 81)))
        for group in self.groups:
            operations = group.point_operations
            spatial_names = {operation.name for operation in operations}
            self.assertEqual(
                spatial_names,
                parent_operations[group.parent_layer_group_number],
                group.og_number,
            )
            lookup = {
                (matrix_key(operation.matrix_fractional), operation.time_reversal)
                for operation in operations
            }
            self.assertEqual(len(lookup), len(operations), group.og_number)
            for left in operations:
                for right in operations:
                    self.assertIn(
                        (
                            matrix_key(
                                matmul(left.matrix_fractional, right.matrix_fractional)
                            ),
                            left.time_reversal ^ right.time_reversal,
                        ),
                        lookup,
                        group.og_number,
                    )

    def test_all_corresponding_msg_records_match_spglib(self):
        for group in self.groups:
            correspondence = group.corresponding_magnetic_space_group
            external_type = spglib.get_magnetic_spacegroup_type(
                correspondence.uni_number
            )
            self.assertEqual(external_type.bns_number, correspondence.bns_number)
            self.assertEqual(external_type.og_number, correspondence.og_number)
            self.assertEqual(external_type.type, TYPE_TO_SPGLIB[group.magnetic_type])
            symmetry = spglib.get_magnetic_symmetry_from_database(
                correspondence.uni_number
            )
            unique = {}
            for rotation, time_reversal in zip(
                symmetry["rotations"], symmetry["time_reversals"], strict=True
            ):
                matrix = tuple(
                    tuple(float(value) for value in row) for row in rotation
                )
                unique[(matrix_key(matrix), bool(time_reversal))] = matrix
            external_signature = colored_signature(
                (matrix, key[1]) for key, matrix in unique.items()
            )
            internal_signature = colored_signature(
                (operation.matrix_fractional, operation.time_reversal)
                for operation in group.point_operations
            )
            self.assertEqual(external_signature, internal_signature, group.og_number)

    def test_type_specific_time_reversal_and_anti_translation(self):
        anti_translations = Counter()
        for group in self.groups:
            pairs = Counter(
                (operation.name, operation.time_reversal)
                for operation in group.point_operations
            )
            pure_time_reversal = ("1", True) in pairs
            if group.magnetic_type == "I":
                self.assertFalse(any(key[1] for key in pairs))
            elif group.magnetic_type == "III":
                self.assertFalse(pure_time_reversal)
                self.assertEqual(group.unitary_subgroup_order * 2, len(pairs))
            else:
                self.assertTrue(pure_time_reversal)
                for name in {key[0] for key in pairs}:
                    self.assertIn((name, False), pairs)
                    self.assertIn((name, True), pairs)
            if group.anti_translation_fractional is not None:
                anti_translations[group.anti_translation_fractional] += 1
        self.assertEqual(
            anti_translations,
            {(1.0, 0.0, 0.0): 47, (0.0, 1.0, 0.0): 43, (0.5, 0.5, 0.0): 32},
        )

    def test_representative_groups_and_identifiers(self):
        expected = {
            1: ("1.1.1", "I", "p1", None),
            2: ("1.2.2", "II", "p11'", None),
            3: ("1.3.3", "IV", "p2a1", (1.0, 0.0, 0.0)),
            25: ("6.5.25", "III", "p112'/m'", None),
            520: ("80.1.520", "I", "p6/mmm", None),
            528: ("80.9.528", "III", "p6'/mm'm", None),
        }
        for number, values in expected.items():
            group = get_magnetic_layer_group(number)
            self.assertEqual(
                (
                    group.og_number,
                    group.magnetic_type,
                    group.litvin_og_symbol_ascii,
                    group.anti_translation_fractional,
                ),
                values,
            )
            self.assertEqual(get_magnetic_layer_group(f"MLG{number}"), group)
            self.assertEqual(get_magnetic_layer_group(group.og_number), group)
        with self.assertRaises(GroupDataError):
            get_magnetic_layer_group("not-a-group")

    def test_all_groups_and_six_optical_sectors_are_solved(self):
        responses = (
            "normal_shift_current",
            "magnetic_shift_current",
            "normal_injection_current",
            "magnetic_injection_current",
            "shg_even",
            "shg_odd",
        )
        screened = {
            (result.group_number, result.response): result.dimension
            for result in screen_response_symmetry(
                "magnetic_layer_group",
                registry=self.registry,
            )
        }
        self.assertEqual(len(screened), 528 * 6)
        for group in self.groups:
            for response in responses:
                result = magnetic_layer_response_tensor_basis(
                    group.global_number,
                    response,
                    registry=self.registry,
                )
                self.assertEqual(result.magnetic_layer_group_number, group.global_number)
                self.assertEqual(result.shape[0], 3)
                self.assertIn(result.shape[1], (3, 6))
                self.assertEqual(
                    screened[(group.global_number, response)],
                    result.dimension,
                )

    def test_gray_point_cogroups_remove_time_odd_homogeneous_tensors(self):
        for identifier in (2, 3):
            self.assertEqual(
                magnetic_layer_tensor_basis(
                    identifier, "polar_vector", "symmetric_quadratic"
                ).dimension,
                18,
            )
            self.assertEqual(
                magnetic_layer_response_tensor_basis(
                    identifier, "magnetic_injection_current"
                ).dimension,
                0,
            )
            self.assertEqual(
                magnetic_layer_response_tensor_basis(identifier, "shg_odd").dimension,
                0,
            )

    def test_cli_query_filter_solver_and_validation(self):
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["magnetic-layer-groups", "6.5.25", "--json"]), 0)
        self.assertEqual(json.loads(output.getvalue())["global_number"], 25)
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(
                    [
                        "magnetic-layer-responses",
                        "6.5.25",
                        "shg_odd",
                        "--json",
                    ]
                ),
                0,
            )
        self.assertEqual(json.loads(output.getvalue())["response"], "shg_odd")
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["magnetic-layer-groups", "--type", "IV"]), 0)
        self.assertEqual(len(output.getvalue().splitlines()), 122)


if __name__ == "__main__":
    unittest.main()
