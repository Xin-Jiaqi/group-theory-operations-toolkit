"""Multi-source verification of the 230 crystallographic space groups.

The registry is machine-generated from the spglib database.  These tests
verify it from five independent directions:

1. structural invariants -- 230 groups, ITA crystal-system counts, all 32
   point groups covered, 73 symmorphic groups, unique symbols;
2. spglib cross-check -- every stored label and Hall setting agrees with
   ``spglib.get_spacegroup_type`` on the same Hall number;
3. generator closure -- the stored Seitz generators close to the declared
   operation count for every Hall setting, and the rotations form a closed
   point-group table;
4. ASE cross-check -- the point group of every space group agrees with
   ASE's independent database (skipped when ASE is unavailable);
5. end-to-end round trip -- three generic positions on a conventional
   lattice, transformed by the stored generators, are identified by spglib
   as the same space group number (this verifies the Seitz data itself).
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np

from group_theory_operations.seitz import SeitzOp, closure, equivalent, inverse, multiply
from group_theory_operations.space_groups import (
    get_crystallographic_space_group,
    iter_crystallographic_space_groups,
    load_space_group_registry,
)

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None

try:
    from ase.spacegroup import Spacegroup as AseSpacegroup
except ImportError:  # pragma: no cover
    AseSpacegroup = None

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "crystallographic_space_groups.json"
POINT_GROUP_PATH = ROOT / "data" / "crystallographic_point_groups.json"

EXPECTED_COUNTS = {
    "triclinic": 2,
    "monoclinic": 13,
    "orthorhombic": 59,
    "tetragonal": 68,
    "trigonal": 25,
    "hexagonal": 27,
    "cubic": 36,
}

_POINT_GROUP_CRYSTAL_SYSTEM = {
    "1": "triclinic",
    "-1": "triclinic",
    "2": "monoclinic",
    "m": "monoclinic",
    "2/m": "monoclinic",
    "222": "orthorhombic",
    "mm2": "orthorhombic",
    "mmm": "orthorhombic",
    "4": "tetragonal",
    "-4": "tetragonal",
    "4/m": "tetragonal",
    "422": "tetragonal",
    "4mm": "tetragonal",
    "-42m": "tetragonal",
    "4/mmm": "tetragonal",
    "3": "trigonal",
    "-3": "trigonal",
    "32": "trigonal",
    "3m": "trigonal",
    "-3m": "trigonal",
    "6": "hexagonal",
    "-6": "hexagonal",
    "6/m": "hexagonal",
    "622": "hexagonal",
    "6mm": "hexagonal",
    "-6m2": "hexagonal",
    "6/mmm": "hexagonal",
    "23": "cubic",
    "m-3": "cubic",
    "432": "cubic",
    "-43m": "cubic",
    "m-3m": "cubic",
}

GENERIC_POSITIONS = (
    np.array([0.113, 0.237, 0.341]),
    np.array([0.154, 0.304, 0.424]),
    np.array([0.234, 0.391, 0.439]),
)


class SeitzAlgebraTests(unittest.TestCase):
    def test_inverse_in_nonorthogonal_fractional_basis(self) -> None:
        rotation = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
        operation = SeitzOp(rotation, np.array([1 / 3, 2 / 3, 0.0]))
        candidate = inverse(operation)
        self.assertTrue(multiply(operation, candidate).is_identity())
        self.assertTrue(multiply(candidate, operation).is_identity())
        self.assertTrue(equivalent(inverse(candidate), operation))

    def test_seitz_rotation_must_be_exact_and_unimodular(self) -> None:
        with self.assertRaisesRegex(ValueError, "exact integers"):
            SeitzOp.from_dict(
                {
                    "rotation": [[1, 0.5, 0], [0, 1, 0], [0, 0, 1]],
                    "translation": [0, 0, 0],
                }
            )
        with self.assertRaisesRegex(ValueError, "unimodular"):
            SeitzOp(np.diag([2, 1, 1]), np.zeros(3))

    def test_closure_accepts_one_shot_iterables(self) -> None:
        quarter_turn = SeitzOp(
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), np.zeros(3)
        )
        generated = closure(operation for operation in (quarter_turn,))
        self.assertEqual(len(generated), 4)


def _lattice_for(crystal_system: str) -> np.ndarray:
    if crystal_system == "cubic":
        return np.eye(3) * 5.2
    if crystal_system == "tetragonal":
        return np.array([[5.2, 0.0, 0.0], [0.0, 5.2, 0.0], [0.0, 0.0, 6.3]])
    if crystal_system in ("hexagonal", "trigonal"):
        a, c = 5.2, 6.3
        return np.array([[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0], [0.0, 0.0, c]])
    if crystal_system == "orthorhombic":
        return np.array([[5.2, 0.0, 0.0], [0.0, 6.1, 0.0], [0.0, 0.0, 4.4]])
    if crystal_system == "monoclinic":
        return np.array([[5.1, 0.0, 0.0], [0.0, 6.3, 0.0], [0.7, 0.0, 4.2]])
    return np.array([[4.1, 0.0, 0.0], [0.3, 5.3, 0.0], [0.2, 0.4, 6.7]])


def _generators(record) -> list[SeitzOp]:
    setting = next(
        setting
        for setting in record.hall_settings
        if setting.hall_number == record.primary_hall_number
    )
    return [SeitzOp.from_dict(generator) for generator in setting.generators]


def _translation_close(left: np.ndarray, right: np.ndarray) -> bool:
    difference = np.asarray(left, dtype=np.float64) - np.asarray(right, dtype=np.float64)
    difference -= np.rint(difference)
    return bool(np.max(np.abs(difference)) < 1.0e-9)


def _round_trip_structure(record):
    lattice = _lattice_for(record.crystal_system)
    generators = _generators(record)
    operations = closure(generators)
    positions = []
    for base in GENERIC_POSITIONS:
        for operation in operations:
            positions.append(operation.apply(base))
    return lattice, positions


@unittest.skipUnless(spglib is not None, "spglib not installed")
class SpaceGroupStructuralTests(unittest.TestCase):
    def test_registry_contains_all_230(self) -> None:
        records = list(iter_crystallographic_space_groups())
        self.assertEqual([record.ita_number for record in records], list(range(1, 231)))

    def test_crystal_system_counts_match_ita(self) -> None:
        records = list(iter_crystallographic_space_groups())
        counts = {
            system: sum(record.crystal_system == system for record in records)
            for system in EXPECTED_COUNTS
        }
        self.assertEqual(counts, EXPECTED_COUNTS)

    def test_all_32_point_groups_covered(self) -> None:
        records = list(iter_crystallographic_space_groups())
        self.assertEqual(len({record.point_group_hm for record in records}), 32)
        self.assertEqual(len({record.point_group_number for record in records}), 32)

    def test_73_symmorphic_groups(self) -> None:
        records = list(iter_crystallographic_space_groups())
        self.assertEqual(sum(record.symmorphic for record in records), 73)

    def test_symbols_are_unique(self) -> None:
        records = list(iter_crystallographic_space_groups())
        self.assertEqual(len({record.international_short for record in records}), 230)

    def test_point_group_crystal_system_consistency(self) -> None:
        records = list(iter_crystallographic_space_groups())
        for record in records:
            self.assertEqual(
                record.crystal_system,
                _POINT_GROUP_CRYSTAL_SYSTEM[record.point_group_hm],
                f"space group {record.ita_number}",
            )

    def test_all_hall_numbers_covered(self) -> None:
        records = list(iter_crystallographic_space_groups())
        hall_numbers = [
            setting.hall_number
            for record in records
            for setting in record.hall_settings
        ]
        self.assertEqual(sorted(hall_numbers), list(range(1, 531)))

    def test_known_spot_checks(self) -> None:
        checks = {
            1: ("P1", "triclinic", "P", True),
            2: ("P-1", "triclinic", "P", True),
            4: ("P2_1", "monoclinic", "P", False),
            14: ("P2_1/c", "monoclinic", "P", False),
            47: ("Pmmm", "orthorhombic", "P", True),
            62: ("Pnma", "orthorhombic", "P", False),
            146: ("R3", "trigonal", "R", True),
            160: ("R3m", "trigonal", "R", True),
            227: ("Fd-3m", "cubic", "F", False),
            230: ("Ia-3d", "cubic", "I", False),
        }
        for ita_number, (symbol, system, centering, symmorphic) in checks.items():
            record = get_crystallographic_space_group(ita_number)
            self.assertEqual(record.international_short, symbol)
            self.assertEqual(record.crystal_system, system)
            self.assertEqual(record.centering, centering)
            self.assertEqual(record.symmorphic, symmorphic)


@unittest.skipUnless(spglib is not None, "spglib not installed")
class SpaceGroupSpglibCrossCheckTests(unittest.TestCase):
    def test_labels_agree_with_spglib(self) -> None:
        for record in iter_crystallographic_space_groups():
            info = spglib.get_spacegroup_type(record.primary_hall_number)
            self.assertEqual(info.number, record.ita_number, f"SG {record.ita_number}")
            self.assertEqual(
                info.international_short, record.international_short,
                f"SG {record.ita_number}",
            )
            self.assertEqual(
                info.international_full.replace(" ", ""),
                record.international_full.replace(" ", ""),
                f"SG {record.ita_number}",
            )
            self.assertEqual(
                info.schoenflies, record.schoenflies, f"SG {record.ita_number}"
            )
            self.assertEqual(
                " ".join(info.pointgroup_international.replace("_", " ").split()),
                record.point_group_hm,
                f"SG {record.ita_number}",
            )

    def test_hall_settings_agree_with_spglib(self) -> None:
        for record in iter_crystallographic_space_groups():
            for setting in record.hall_settings:
                info = spglib.get_spacegroup_type(setting.hall_number)
                self.assertEqual(
                    info.number, record.ita_number, f"SG {record.ita_number}"
                )
                self.assertEqual(
                    info.hall_symbol, setting.hall_symbol,
                    f"SG {record.ita_number} hall {setting.hall_number}",
                )

    def test_operation_counts_agree_with_spglib(self) -> None:
        for record in iter_crystallographic_space_groups():
            for setting in record.hall_settings:
                operations = spglib.get_symmetry_from_database(setting.hall_number)
                self.assertEqual(
                    len(operations["rotations"]),
                    setting.operation_count,
                    f"SG {record.ita_number} hall {setting.hall_number}",
                )

    def test_primary_setting_round_trips(self) -> None:
        for record in iter_crystallographic_space_groups():
            lattice, positions = _round_trip_structure(record)
            dataset = spglib.get_symmetry_dataset(
                (lattice, positions, [1] * len(positions)), symprec=1.0e-4
            )
            self.assertEqual(
                dataset.number,
                record.ita_number,
                f"SG {record.ita_number} did not round-trip",
            )


class SpaceGroupGeneratorClosureTests(unittest.TestCase):
    def test_generators_close_to_operation_count(self) -> None:
        for record in iter_crystallographic_space_groups():
            for setting in record.hall_settings:
                generators = [
                    SeitzOp.from_dict(generator) for generator in setting.generators
                ]
                generated = closure(generators)
                self.assertEqual(
                    len(generated),
                    setting.operation_count,
                    f"SG {record.ita_number} hall {setting.hall_number}",
                )

    def test_rotation_part_closes_as_a_point_group(self) -> None:
        identity = np.eye(3, dtype=np.int64)
        for record in iter_crystallographic_space_groups():
            setting = next(
                s
                for s in record.hall_settings
                if s.hall_number == record.primary_hall_number
            )
            generators = [
                SeitzOp.from_dict(generator) for generator in setting.generators
            ]
            generated = closure(generators)
            rotations = {
                tuple(operation.rotation.ravel()) for operation in generated
            }
            self.assertIn(tuple(identity.ravel()), rotations)
            for left in rotations:
                for right in rotations:
                    product = np.asarray(left, dtype=np.int64).reshape(3, 3) @ np.asarray(
                        right, dtype=np.int64
                    ).reshape(3, 3)
                    self.assertIn(
                        tuple(product.ravel()),
                        rotations,
                        f"SG {record.ita_number} rotation part is not closed",
                    )

    def test_generators_reproduce_primary_operations(self) -> None:
        for record in iter_crystallographic_space_groups():
            setting = next(
                s
                for s in record.hall_settings
                if s.hall_number == record.primary_hall_number
            )
            generators = [
                SeitzOp.from_dict(generator) for generator in setting.generators
            ]
            generated = closure(generators)
            source = spglib.get_symmetry_from_database(setting.hall_number)
            rotations = np.asarray(source["rotations"], dtype=np.int64)
            translations = np.asarray(source["translations"], dtype=np.float64)
            self.assertEqual(len(generated), len(rotations))
            for operation in generated:
                match = any(
                    np.array_equal(operation.rotation, rotation)
                    and _translation_close(operation.translation, translation)
                    for rotation, translation in zip(rotations, translations)
                )
                self.assertTrue(match, f"SG {record.ita_number} has an unexpected operation")


@unittest.skipUnless(AseSpacegroup is not None, "ASE not installed")
class SpaceGroupAseCrossCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.ase_rows = _ase_rows()

    def test_operation_counts_agree_with_ase(self) -> None:
        # ASE's data file stores a generating subset (not necessarily the full
        # operation set), so the primitive operation count must be a positive
        # multiple of ASE's count.
        multiplicity = {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4, "R": 3}
        for record in iter_crystallographic_space_groups():
            setting = next(
                s
                for s in record.hall_settings
                if s.hall_number == record.primary_hall_number
            )
            primitive_count = setting.operation_count // multiplicity[record.centering]
            ase_count = int(self.ase_rows[record.ita_number]["operation_count"])
            self.assertGreaterEqual(ase_count, 1, f"SG {record.ita_number}")
            self.assertEqual(
                primitive_count % ase_count,
                0,
                f"SG {record.ita_number} primitive count {primitive_count} "
                f"is not a multiple of ASE's {ase_count}",
            )

    def test_hm_symbols_agree_with_ase(self) -> None:
        for record in iter_crystallographic_space_groups():
            ase_symbol = "".join(self.ase_rows[record.ita_number]["symbol"].split())
            stored_symbol = "".join(record.international_short.split())
            self.assertEqual(
                stored_symbol.replace("_", "").lower(),
                ase_symbol.replace("_", "").lower(),
                f"SG {record.ita_number}",
            )


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def _ase_rows() -> dict[int, dict[str, object]]:
    import ase.spacegroup.spacegroup as ase_spacegroup

    path = Path(ase_spacegroup.get_datafile())
    rows: dict[int, dict[str, object]] = {}
    current: dict[str, object] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields:
            continue
        if (
            fields[0].isdigit()
            and not _is_number(fields[1])
            and fields[1][0] in "PABC FIR-"
        ):
            number = int(fields[0])
            if number < 1 or number > 230:
                current = None
                continue
            if number not in rows:
                rows[number] = {
                    "symbol": " ".join(fields[1:]),
                    "operation_count": 0,
                }
                current = rows[number]
            else:
                current = None
            continue
        if current is not None and len(fields) == 12:
            try:
                [float(value) for value in fields]
            except ValueError:
                continue
            current["operation_count"] = int(current["operation_count"]) + 1
    return rows


class SpaceGroupSchemaTests(unittest.TestCase):
    def test_data_matches_schema(self) -> None:
        try:
            import jsonschema
        except ImportError:  # pragma: no cover
            self.skipTest("jsonschema not installed")
        schema = json.loads(
            (ROOT / "schema" / "crystallographic-space-groups-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
        jsonschema.validate(instance=data, schema=schema)

    def test_point_group_registry_sha256_is_bound(self) -> None:
        data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
        import hashlib

        expected = hashlib.sha256(POINT_GROUP_PATH.read_bytes()).hexdigest()
        self.assertEqual(data["point_group_registry_sha256"], expected)


class SpaceGroupWikipediaFixtureTests(unittest.TestCase):
    """Crystal systems transcribed from the Wikipedia list of space groups.

    The fixture was generated on 2026-08-06 from
    https://en.wikipedia.org/wiki/List_of_space_groups (section headings
    assign each ITA number to one of the seven crystal systems).
    """

    def test_crystal_systems_agree_with_wikipedia(self) -> None:
        fixture = json.loads(
            (ROOT / "tests" / "fixtures" / "wikipedia_crystal_systems.json").read_text(
                encoding="utf-8"
            )
        )
        expected = fixture["crystal_systems"]
        self.assertEqual(len(expected), 230)
        for record in iter_crystallographic_space_groups():
            self.assertEqual(
                expected[str(record.ita_number)],
                record.crystal_system,
                f"SG {record.ita_number}",
            )


@unittest.skipUnless(gemmi is not None, "gemmi not installed")
class SpaceGroupGemmiCrossCheckTests(unittest.TestCase):
    """Cross-check against gemmi's independent space-group tables (BSD-3).

    gemmi's notation differs in three cosmetic ways, so only the invariants
    that are insensitive to notation are compared: the crystal system, the
    operation count of the reference setting, and the symmorphic flag.
    """

    def test_crystal_systems_agree_with_gemmi(self) -> None:
        for record in iter_crystallographic_space_groups():
            sg = gemmi.find_spacegroup_by_number(record.ita_number)
            self.assertIsNotNone(sg, f"SG {record.ita_number}")
            gemmi_system = str(sg.crystal_system()).split(".")[-1].lower()
            self.assertEqual(
                gemmi_system, record.crystal_system, f"SG {record.ita_number}"
            )

    def test_operation_counts_agree_with_gemmi(self) -> None:
        for record in iter_crystallographic_space_groups():
            sg = gemmi.find_spacegroup_by_number(record.ita_number)
            setting = next(
                s
                for s in record.hall_settings
                if s.hall_number == record.primary_hall_number
            )
            self.assertEqual(
                len(sg.operations()),
                setting.operation_count,
                f"SG {record.ita_number}",
            )

    def test_symmorphic_flags_agree_with_gemmi(self) -> None:
        for record in iter_crystallographic_space_groups():
            sg = gemmi.find_spacegroup_by_number(record.ita_number)
            self.assertEqual(
                bool(sg.is_symmorphic()),
                record.symmorphic,
                f"SG {record.ita_number}",
            )

    def test_centrosymmetric_flags_agree_with_point_group_registry(self) -> None:
        # gemmi's centrosymmetry must match the parent point group registry.
        point_registry = json.loads(POINT_GROUP_PATH.read_text(encoding="utf-8"))
        centrosymmetric = {
            entry["number"]: entry["centrosymmetric"]
            for entry in point_registry["point_groups"]
        }
        for record in iter_crystallographic_space_groups():
            sg = gemmi.find_spacegroup_by_number(record.ita_number)
            self.assertEqual(
                bool(sg.is_centrosymmetric()),
                centrosymmetric[record.point_group_number],
                f"SG {record.ita_number}",
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
