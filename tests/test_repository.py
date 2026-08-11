from __future__ import annotations

import json
import math
import re
import unittest
from pathlib import Path

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import group_tools  # noqa: E402
import group_theory_operations  # noqa: E402
from jsonschema import Draft202012Validator  # noqa: E402


class MinimalStructureRecord:
    """Small contract double; real cross-repository integration is tested separately."""

    def __init__(self, lattice, species, fractional_coordinates, pbc=(True, True, True), selective_dynamics=None, length_unit="angstrom"):
        self.lattice = tuple(tuple(row) for row in lattice)
        self.species = tuple(species)
        self.fractional_coordinates = tuple(tuple(row) for row in fractional_coordinates)
        self.pbc = tuple(pbc)
        self.selective_dynamics = None if selective_dynamics is None else tuple(tuple(row) for row in selective_dynamics)
        self.length_unit = length_unit

    @classmethod
    def from_fractional(cls, **kwargs):
        return cls(**kwargs)

    def wrapped(self):
        coordinates = [
            tuple(value % 1.0 if self.pbc[axis] else value for axis, value in enumerate(row))
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


def matmul(left, right):
    return [
        [sum(float(left[i][k]) * float(right[k][j]) for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def matrix_key(matrix, digits=9):
    return tuple(round(float(value), digits) for row in matrix for value in row)


class RepositoryDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = group_tools.load_database(ROOT / "data" / "group_operations.json")

    def test_expected_dataset_sizes(self):
        families = self.database["families"]
        self.assertEqual(len(families["cubic_Oh"]["operations"]), 48)
        self.assertEqual(len(families["cubic_Oh"]["point_groups"]), 5)
        self.assertEqual(len(families["tetragonal_D4h"]["operations"]), 16)
        self.assertEqual(len(families["tetragonal_D4h"]["layer_groups"]), 64)
        self.assertEqual(len(families["hexagonal_D6h"]["operations"]), 24)
        self.assertEqual(len(families["hexagonal_D6h"]["layer_groups"]), 16)

    def test_release_version_and_numpy_dependency_are_consistent(self):
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        version = re.search(r'^version = "([^"]+)"$', pyproject, re.MULTILINE)
        self.assertIsNotNone(version)
        self.assertEqual(version.group(1), group_theory_operations.__version__)
        self.assertIn(f"version: {group_theory_operations.__version__}", citation)
        self.assertRegex(pyproject, r'dependencies = \[[^\]]*"numpy>=1\.24,<3"')
        self.assertRegex(
            pyproject,
            r'(?m)^structure = \[.*"spglib>=2\.5,<3"\]$',
        )

    def test_public_validator_accepts_canonical_catalog(self):
        self.assertEqual(group_tools.validate_database(self.database), ())

    def test_catalog_satisfies_versioned_json_schema(self):
        schema = group_tools.load_schema(ROOT / "schema" / "group-operations-v1.schema.json")
        errors = sorted(Draft202012Validator(schema).iter_errors(self.database), key=str)
        self.assertEqual(errors, [])

    def test_public_validator_reports_matrix_corruption(self):
        corrupted = json.loads(json.dumps(self.database))
        corrupted["families"]["tetragonal_D4h"]["operations"][0]["matrix_cartesian"][0][0] = 2
        errors = group_tools.validate_database(corrupted)
        self.assertTrue(any("not orthogonal" in error for error in errors))

    def test_public_validator_is_total_for_malformed_json(self):
        mutations = []
        for path, value in (
            (("families", "tetragonal_D4h"), []),
            (("families", "tetragonal_D4h", "operations", 0), "bad"),
            (("families", "cubic_Oh", "point_groups", 0), "bad"),
        ):
            corrupted = json.loads(json.dumps(self.database))
            target = corrupted
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            mutations.append(corrupted)
        for corrupted in mutations:
            errors = group_tools.validate_database(corrupted)
            self.assertTrue(errors)

    def test_public_validator_requires_integer_schema_version(self):
        for invalid in (True, 1.0):
            corrupted = json.loads(json.dumps(self.database))
            corrupted["schema_version"] = invalid
            self.assertEqual(
                group_tools.validate_database(corrupted), ("schema_version must be 1",)
            )

    def test_public_validator_rejects_nested_or_boolean_operation_indices(self):
        for invalid in ([[]], [True]):
            corrupted = json.loads(json.dumps(self.database))
            corrupted["families"]["cubic_Oh"]["point_groups"][0][
                "operation_indices"
            ] = invalid
            errors = group_tools.validate_database(corrupted)
            self.assertTrue(any("indices must be integers" in error for error in errors))

    def test_public_validator_is_total_for_extreme_numbers_and_table_values(self):
        mutations = []
        huge_matrix = json.loads(json.dumps(self.database))
        huge_matrix["families"]["cubic_Oh"]["operations"][0]["matrix_fractional"][0][
            0
        ] = 10**400
        mutations.append(huge_matrix)
        for invalid in ([], {}):
            corrupted = json.loads(json.dumps(self.database))
            corrupted["families"]["tetragonal_D4h"]["multiplication"]["table"]["1"][
                "1"
            ] = invalid
            mutations.append(corrupted)
        invalid_row = json.loads(json.dumps(self.database))
        invalid_row["families"]["hexagonal_D6h"]["multiplication"]["table"][
            "-3+_001"
        ] = True
        mutations.append(invalid_row)
        for corrupted in mutations:
            errors = group_tools.validate_database(corrupted)
            self.assertTrue(errors)

    def test_public_validator_rejects_boolean_and_string_matrix_scalars(self):
        for invalid in (True, "1"):
            corrupted = json.loads(json.dumps(self.database))
            corrupted["families"]["tetragonal_D4h"]["operations"][0][
                "matrix_fractional"
            ][0][0] = invalid
            errors = group_tools.validate_database(corrupted)
            self.assertTrue(any("JSON numbers" in error for error in errors))

    def test_public_validator_rejects_corrupt_group_metadata(self):
        mutations = []
        for family, path, value in (
            ("cubic_Oh", ("point_groups", 0, "name"), []),
            ("cubic_Oh", ("point_groups", 0, "order"), True),
            ("tetragonal_D4h", ("layer_groups", 0, "LG"), "bad"),
            ("tetragonal_D4h", ("layer_groups", 0, "point_group"), {}),
            ("tetragonal_D4h", ("layer_groups", 0, "point_group_base"), ""),
            ("hexagonal_D6h", ("layer_groups", 0, "point_group_embedding"), []),
        ):
            corrupted = json.loads(json.dumps(self.database))
            target = corrupted["families"][family]
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            mutations.append(corrupted)
        parent = json.loads(json.dumps(self.database))
        parent["families"]["cubic_Oh"]["parent_point_group"] = {}
        mutations.append(parent)
        source_name = json.loads(json.dumps(self.database))
        source_name["families"]["cubic_Oh"]["operations"][0]["source_name"] = []
        mutations.append(source_name)
        for corrupted in mutations:
            self.assertTrue(group_tools.validate_database(corrupted))

    def test_public_validator_checks_basis_subgroups_and_multiplication_metadata(self):
        corrupted = json.loads(json.dumps(self.database))
        family = corrupted["families"]["tetragonal_D4h"]
        family["multiplication"]["identity"] = "-1"
        family["multiplication"]["inverse"]["4+_001"] = "4+_001"
        family["layer_groups"][0]["R+_indices"] = [0, 0]
        family["layer_groups"][0]["R+"] = ["1", "1"]
        family["operations"][2]["matrix_cartesian"] = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
        errors = group_tools.validate_database(corrupted)
        self.assertTrue(any("wrong identity" in error for error in errors))
        self.assertTrue(any("wrong inverse" in error for error in errors))
        self.assertTrue(any("unique" in error or "not closed" in error for error in errors))
        self.assertTrue(any("basis mapping" in error for error in errors))

    def test_typed_operation_record_is_json_safe(self):
        record = group_tools.operation_record(self.database, "4^+_{001}", "tetragonal_D4h")
        self.assertEqual(record.name, "4+_001")
        self.assertEqual(record.family, "tetragonal_D4h")
        self.assertEqual(record.matrix("fractional")[0], (0.0, -1.0, 0.0))
        json.dumps(record.to_dict())

    def test_operation_indices_names_and_matrix_shapes(self):
        for family in self.database["families"].values():
            operations = family["operations"]
            self.assertEqual(len({op["index"] for op in operations}), len(operations))
            self.assertEqual(len({op["name"] for op in operations}), len(operations))
            for operation in operations:
                for key in ("matrix_fractional", "matrix_cartesian"):
                    matrix = operation[key]
                    self.assertEqual(len(matrix), 3)
                    self.assertTrue(all(len(row) == 3 for row in matrix))

    def test_cartesian_matrices_are_orthogonal(self):
        for _, operation in group_tools.iter_operations(self.database):
            matrix = operation["matrix_cartesian"]
            transpose = [list(row) for row in zip(*matrix)]
            product = matmul(transpose, matrix)
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(product[i][j], 1.0 if i == j else 0.0, places=9)

    def test_group_references_and_closure(self):
        for family_name, family in self.database["families"].items():
            by_index = {op["index"]: op for op in family["operations"]}
            groups = []
            for point_group in family["point_groups"]:
                names = [by_index[index]["name"] for index in point_group["operation_indices"]]
                self.assertEqual(names, point_group["operations"])
                groups.append(point_group["operation_indices"])
            for layer_group in family["layer_groups"]:
                indices = layer_group["R+_indices"] + layer_group["R-_indices"]
                names = [by_index[index]["name"] for index in indices]
                self.assertEqual(names, layer_group["R+"] + layer_group["R-"])
                groups.append(indices)

            for indices in groups:
                matrices = [by_index[index]["matrix_cartesian"] for index in indices]
                keys = {matrix_key(matrix) for matrix in matrices}
                for left in matrices:
                    for right in matrices:
                        self.assertIn(matrix_key(matmul(left, right)), keys, family_name)

    def test_hexagonal_basis_conversion(self):
        family = self.database["families"]["hexagonal_D6h"]
        root_three = math.sqrt(3.0)
        basis = [[1.0, -0.5, 0.0], [0.0, root_three / 2.0, 0.0], [0.0, 0.0, 1.0]]
        inverse = [[1.0, root_three / 3.0, 0.0], [0.0, 2.0 * root_three / 3.0, 0.0], [0.0, 0.0, 1.0]]
        for operation in family["operations"]:
            converted = matmul(matmul(basis, operation["matrix_fractional"]), inverse)
            self.assertEqual(matrix_key(converted), matrix_key(operation["matrix_cartesian"]))

    def test_historical_aliases(self):
        operation = group_tools.get_operation(self.database, "4^+_{001}", "tetragonal_D4h")
        self.assertEqual(operation["name"], "4+_001")
        operation = group_tools.get_operation(self.database, "m100", "hexagonal_D6h")
        self.assertEqual(operation["name"], "m_100")

    def test_operation_order_and_coordinate_actions(self):
        expected = {
            "tetragonal_D4h": [
                ("1", "x,y,z"), ("2001", "-x,-y,z"), ("4+001", "-y,x,z"),
                ("4-001", "y,-x,z"), ("2010", "-x,y,-z"), ("2100", "x,-y,-z"),
                ("2110", "y,x,-z"), ("2_1-10", "-y,-x,-z"), ("-1", "-x,-y,-z"),
                ("m001", "x,y,-z"), ("-4+001", "y,-x,-z"), ("-4-001", "-y,x,-z"),
                ("m010", "x,-y,z"), ("m100", "-x,y,z"), ("m110", "-y,-x,z"),
                ("m1-10", "y,x,z"),
            ],
            "hexagonal_D6h": [
                ("1", "x,y,z"), ("3+001", "-y,x-y,z"), ("3-001", "-x+y,-x,z"),
                ("2001", "-x,-y,z"), ("6-001", "y,-x+y,z"), ("6+001", "x-y,x,z"),
                ("2110", "y,x,-z"), ("2100", "x-y,-y,-z"), ("2010", "-x,-x+y,-z"),
                ("2_1-10", "-y,-x,-z"), ("2120", "-x+y,y,-z"), ("2210", "x,x-y,-z"),
                ("-1", "-x,-y,-z"), ("-3+001", "y,-x+y,-z"), ("-3-001", "x-y,x,-z"),
                ("m001", "x,y,-z"), ("-6-001", "-y,x-y,-z"), ("-6+001", "-x+y,-x,-z"),
                ("m110", "-y,-x,z"), ("m100", "-x+y,y,z"), ("m010", "x,x-y,z"),
                ("m1-10", "y,x,z"), ("m120", "x-y,-y,z"), ("m210", "-x,-x+y,z"),
            ],
        }
        for family_name, reference_rows in expected.items():
            operations = self.database["families"][family_name]["operations"]
            actual = [(op.get("source_name", op["name"]), op["xyz_fractional"]) for op in operations]
            self.assertEqual(actual, reference_rows)

    def test_layer_group_embedding_labels(self):
        family = self.database["families"]["hexagonal_D6h"]
        actual = {
            entry["LG"]: entry.get("point_group_embedding")
            for entry in family["layer_groups"]
            if entry["LG"] in {67, 68, 69, 70, 71, 72, 78, 79}
        }
        self.assertEqual(actual, {
            67: "120", 68: "100", 69: "100", 70: "120",
            71: "120", 72: "100", 78: "100", 79: "120",
        })

    def test_layer_group_point_group_catalog(self):
        expected = {
            **{1: "C1", 2: "Ci", 3: "C2"},
            **{lg: "Cs" for lg in range(4, 6)},
            **{lg: "C2h" for lg in range(6, 8)},
            **{lg: "C2" for lg in range(8, 11)},
            **{lg: "Cs" for lg in range(11, 14)},
            **{lg: "C2h" for lg in range(14, 19)},
            **{lg: "D2" for lg in range(19, 23)},
            **{lg: "C2v" for lg in range(23, 37)},
            **{lg: "D2h" for lg in range(37, 49)},
            49: "C4", 50: "S4",
            **{lg: "C4h" for lg in range(51, 53)},
            **{lg: "D4" for lg in range(53, 55)},
            **{lg: "C4v" for lg in range(55, 57)},
            **{lg: "D2d" for lg in range(57, 61)},
            **{lg: "D4h" for lg in range(61, 65)},
            65: "C3", 66: "C3i", 67: "D3", 68: "D3",
            69: "C3v", 70: "C3v", 71: "D3d", 72: "D3d",
            73: "C6", 74: "C3h", 75: "C6h", 76: "D6",
            77: "C6v", 78: "D3h", 79: "D3h", 80: "D6h",
        }
        actual = {}
        for family_name in ("tetragonal_D4h", "hexagonal_D6h"):
            for entry in self.database["families"][family_name]["layer_groups"]:
                actual[entry["LG"]] = entry["point_group_base"]
        self.assertEqual(actual, expected)

    def test_d4h_and_d6h_multiplication_tables(self):
        expected_orders = {
            "tetragonal_D4h": [
                "1", "-1", "2_001", "m_001", "2_100", "2_010", "m_100", "m_010",
                "4+_001", "4-_001", "-4+_001", "-4-_001", "2_110", "2_1-10", "m_110", "m_1-10",
            ],
            "hexagonal_D6h": [
                "1", "-1", "2_001", "m_001", "2_100", "2_010", "m_100", "m_010",
                "3+_001", "3-_001", "-3+_001", "-3-_001", "6+_001", "6-_001",
                "-6+_001", "-6-_001", "2_110", "2_120", "2_210", "2_1-10",
                "m_110", "m_120", "m_210", "m_1-10",
            ],
        }
        for family_name, expected_order in expected_orders.items():
            family = self.database["families"][family_name]
            multiplication = family["multiplication"]
            order = multiplication["element_order"]
            table = multiplication["table"]
            self.assertEqual(order, expected_order)
            self.assertEqual(set(order), {op["name"] for op in family["operations"]})
            self.assertEqual(list(table), order)

            by_name = {op["name"]: op["matrix_fractional"] for op in family["operations"]}
            for left in order:
                self.assertEqual(list(table[left]), order)
                self.assertEqual(set(table[left].values()), set(order))
                for right in order:
                    result = table[left][right]
                    self.assertEqual(
                        matrix_key(matmul(by_name[left], by_name[right])),
                        matrix_key(by_name[result]),
                        f"{family_name}: {left} * {right}",
                    )
            for right in order:
                self.assertEqual({table[left][right] for left in order}, set(order))
            self.assertEqual([table["1"][right] for right in order], order)
            self.assertEqual([table[left]["1"] for left in order], order)
            for element, inverse in multiplication["inverse"].items():
                self.assertEqual(table[element][inverse], "1")
                self.assertEqual(table[inverse][element], "1")
            for left in order:
                for middle in order:
                    for right in order:
                        self.assertEqual(
                            table[table[left][middle]][right],
                            table[left][table[middle][right]],
                        )

    def test_multiplication_query_accepts_aliases(self):
        result = group_tools.multiply_operations(
            self.database, "tetragonal_D4h", "4^+_{001}", "2_100"
        )
        self.assertEqual(result, "2_110")

    def test_markdown_multiplication_tables_match_json(self):
        document = (ROOT / "docs" / "group_theory.md").read_text(encoding="utf-8")
        self.assertNotIn("<details", document)

        def display(name):
            return {"1": "E", "-1": "I"}.get(name, name)

        for family_name in ("tetragonal_D4h", "hexagonal_D6h"):
            multiplication = self.database["families"][family_name]["multiplication"]
            order = multiplication["element_order"]
            table = multiplication["table"]
            for start in range(0, len(order), 8):
                columns = order[start:start + 8]
                header = "| 行 $\\times$ 列 | " + " | ".join(
                    f"`{display(name)}`" for name in columns
                ) + " |"
                self.assertIn(header, document)
                for left in order:
                    cells = " | ".join(
                        f"`{display(table[left][right])}`" for right in columns
                    )
                    row = f"| **`{display(left)}`** | {cells} |"
                    self.assertIn(row, document)

    def test_markdown_table_rows_have_consistent_columns(self):
        for path in (ROOT / "README.md", ROOT / "docs" / "group_theory.md"):
            lines = path.read_text(encoding="utf-8").splitlines()
            index = 0
            while index < len(lines):
                if not lines[index].startswith("|"):
                    index += 1
                    continue
                block = []
                while index < len(lines) and lines[index].startswith("|"):
                    block.append(lines[index])
                    index += 1
                self.assertGreaterEqual(len(block), 2, path)
                pipe_count = block[0].count("|")
                self.assertTrue(all(line.count("|") == pipe_count for line in block), path)

    def test_rendered_markdown_uses_github_compatible_math(self):
        markdown_paths = (
            ROOT / "README.md",
            ROOT / "ROADMAP.md",
            *(ROOT / "docs").glob("*.md"),
        )
        for path in markdown_paths:
            document = path.read_text(encoding="utf-8")
            self.assertNotIn(r"\operatorname", document, path)

        readme_lines = (ROOT / "README.md").read_text(encoding="utf-8").splitlines()
        for line in readme_lines:
            if line.startswith("|"):
                self.assertNotIn("$", line, line)

    def test_structure_contract_applies_fractional_operation(self):
        structure = MinimalStructureRecord(
            lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 12]],
            species=["B", "N"],
            fractional_coordinates=[[0.25, 0.5, 0.75], [0.1, 0.2, 0.3]],
            pbc=[True, True, False],
            selective_dynamics=[[True, False, True], [True, True, True]],
        )
        operation = group_tools.operation_record(self.database, "4+_001", "tetragonal_D4h")
        result = group_tools.apply_fractional_operation(structure, operation)
        self.assertEqual(result.fractional_coordinates[0], (0.5, 0.25, 0.75))
        self.assertEqual(result.selective_dynamics[0], (False, True, True))
        self.assertEqual(result.species, structure.species)

    def test_structure_contract_rejects_pbc_axis_mixing(self):
        structure = MinimalStructureRecord(
            lattice=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            species=["H"],
            fractional_coordinates=[[0.1, 0.2, 0.3]],
            pbc=[True, True, False],
        )
        operation = group_tools.operation_record(self.database, "3+_111", "cubic_Oh")
        with self.assertRaisesRegex(group_tools.GroupDataError, "PBC"):
            group_tools.apply_fractional_operation(structure, operation)

    def test_structure_contract_rejects_unrepresentable_selective_flags(self):
        structure = MinimalStructureRecord(
            lattice=[[1, 0, 0], [-0.5, math.sqrt(3) / 2, 0], [0, 0, 1]],
            species=["H"],
            fractional_coordinates=[[0.1, 0.2, 0.3]],
            selective_dynamics=[[True, False, True]],
        )
        operation = group_tools.operation_record(self.database, "3+_001", "hexagonal_D6h")
        with self.assertRaisesRegex(group_tools.GroupDataError, "Selective"):
            group_tools.apply_fractional_operation(structure, operation)

    def test_structure_contract_rejects_incompatible_lattice_metric(self):
        structure = MinimalStructureRecord(
            lattice=[[2, 0, 0], [0, 3, 0], [0, 0, 4]],
            species=["H"],
            fractional_coordinates=[[0.25, 0.0, 0.0]],
        )
        operation = group_tools.operation_record(
            self.database, "4+_001", "tetragonal_D4h"
        )
        with self.assertRaisesRegex(group_tools.GroupDataError, "lattice metric"):
            group_tools.apply_fractional_operation(structure, operation)

    def test_structure_contract_accepts_hexagonal_metric(self):
        structure = MinimalStructureRecord(
            lattice=[[2, 0, 0], [-1, math.sqrt(3), 0], [0, 0, 7]],
            species=["H"],
            fractional_coordinates=[[0.25, 0.5, 0.0]],
        )
        operation = group_tools.operation_record(
            self.database, "3+_001", "hexagonal_D6h"
        )
        result = group_tools.apply_fractional_operation(structure, operation)
        self.assertEqual(result.species, structure.species)

    def test_structure_application_composes_in_table_order(self):
        structure = MinimalStructureRecord(
            lattice=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            species=["H"],
            fractional_coordinates=[[0.13, 0.27, 0.41]],
        )
        family = "tetragonal_D4h"
        left = group_tools.operation_record(self.database, "4+_001", family)
        right = group_tools.operation_record(self.database, "2_100", family)
        result_name = group_tools.multiply_operations(self.database, family, left.name, right.name)
        composed = group_tools.operation_record(self.database, result_name, family)
        sequential = group_tools.apply_fractional_operation(
            group_tools.apply_fractional_operation(structure, right, wrap=False), left, wrap=False
        )
        direct = group_tools.apply_fractional_operation(structure, composed, wrap=False)
        for actual, expected in zip(sequential.fractional_coordinates[0], direct.fractional_coordinates[0]):
            self.assertAlmostEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
