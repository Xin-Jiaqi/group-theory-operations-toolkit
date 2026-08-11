"""Concrete-structure screening of nonmagnetic nonlinear responses."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

from group_theory_operations import (
    StructureResponseAnalysis,
    analyze_structure_responses,
)
from group_theory_operations.cli import main

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "spglib_real_structures_v2.5.0.json"

EXPECTED_RESPONSE_DIMENSIONS = {
    "sio2_triclinic": (18, 18, 9),
    "sio2_monoclinic": (8, 8, 5),
    "bates3_orthorhombic": (0, 0, 0),
    "mno2_tetragonal": (0, 0, 0),
    "rucl3_trigonal": (4, 4, 1),
    "aucn_hexagonal": (3, 3, 1),
    "cssnbr3_cubic": (0, 0, 0),
}


class StructureDouble:
    def __init__(self, item):
        self.lattice = tuple(tuple(row) for row in item["lattice"])
        self.species = tuple(item["species"])
        self.fractional_coordinates = tuple(
            tuple(row) for row in item["fractional_coordinates"]
        )
        self.pbc = (True, True, True)


def _load_fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _vasp_text(item: dict) -> str:
    species_order = tuple(dict.fromkeys(item["species"]))
    counts = [item["species"].count(symbol) for symbol in species_order]
    if list(item["species"]) != [
        symbol for symbol, count in zip(species_order, counts) for _ in range(count)
    ]:
        raise AssertionError("test POSCAR species must be grouped")
    lines = [item["formula"], "1.0"]
    lines.extend(" ".join(str(value) for value in row) for row in item["lattice"])
    lines.append(" ".join(species_order))
    lines.append(" ".join(str(count) for count in counts))
    lines.append("Direct")
    lines.extend(
        " ".join(str(value) for value in row)
        for row in item["fractional_coordinates"]
    )
    return "\n".join(lines) + "\n"


@unittest.skipUnless(spglib is not None, "spglib not installed")
class StructureResponseAnalysisTests(unittest.TestCase):
    def test_all_seven_crystal_system_fixtures_match_expected_dimensions(self):
        for item in _load_fixture()["structures"]:
            with self.subTest(item=item["id"]):
                analysis = analyze_structure_responses(StructureDouble(item))
                expected = item["expected"]
                self.assertIsInstance(analysis, StructureResponseAnalysis)
                self.assertEqual(
                    analysis.symmetry.space_group.ita_number,
                    expected["ita_number"],
                )
                self.assertEqual(
                    analysis.symmetry.space_group.point_group_hm,
                    expected["point_group_hm"],
                )
                self.assertEqual(
                    tuple(result.dimension for result in analysis.responses),
                    EXPECTED_RESPONSE_DIMENSIONS[item["id"]],
                )
                self.assertTrue(
                    all(
                        result.group_number
                        == analysis.symmetry.space_group.point_group_number
                        for result in analysis.responses
                    )
                )

    def test_centrosymmetric_structure_has_no_selected_second_order_response(self):
        item = _load_fixture()["structures"][-1]
        analysis = analyze_structure_responses(
            StructureDouble(item),
            allowed_only=True,
        )
        self.assertEqual(analysis.symmetry.space_group.point_group_hm, "m-3m")
        self.assertEqual(analysis.responses, ())
        self.assertEqual(analysis.allowed_responses, ())

    def test_selected_response_and_json_summary_preserve_physical_context(self):
        item = _load_fixture()["structures"][-2]
        analysis = analyze_structure_responses(
            StructureDouble(item),
            responses=("shift", "cpge"),
        )
        self.assertEqual(
            [(result.response, result.dimension) for result in analysis.responses],
            [("shift_current", 3), ("circular_injection_current", 1)],
        )
        payload = analysis.to_dict()
        self.assertEqual(payload["space_group"]["ita_number"], 183)
        self.assertEqual(payload["space_group"]["international_short"], "P6mm")
        self.assertEqual(payload["point_group"]["hm_symbol"], "6mm")
        self.assertEqual(payload["classification"]["backend"], "spglib")
        self.assertEqual(len(payload["responses"]), 2)

    def test_cli_reads_poscar_and_returns_structure_responses(self):
        item = _load_fixture()["structures"][-2]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(_vasp_text(item), encoding="utf-8")
            output = io.StringIO()
            with redirect_stdout(output):
                self.assertEqual(
                    main(
                        [
                            "structure-responses",
                            str(path),
                            "--input-format",
                            "vasp",
                            "--response",
                            "shift_current",
                            "--response",
                            "circular_injection_current",
                            "--json",
                        ]
                    ),
                    0,
                )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["space_group"]["ita_number"], 183)
        self.assertEqual(payload["point_group"]["hm_symbol"], "6mm")
        self.assertEqual(
            [(result["response"], result["dimension"]) for result in payload["responses"]],
            [("shift_current", 3), ("circular_injection_current", 1)],
        )


if __name__ == "__main__":
    unittest.main()
