"""Occupied Wyckoff-orbit and site-stabilizer tests."""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from group_theory_operations import (
    WyckoffOrbitAnalysis,
    analyze_wyckoff_orbits,
)
from group_theory_operations.cli import main

try:
    import spglib
except ImportError:  # pragma: no cover
    spglib = None


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "spglib_real_structures_v2.5.0.json"
MOYO_FIXTURE_PATH = (
    ROOT / "tests" / "fixtures" / "moyo_wyckoff_crosscheck_v0.10.0.json"
)


def _partition(labels) -> set[frozenset[int]]:
    groups: dict[int, set[int]] = {}
    for index, label in enumerate(labels):
        groups.setdefault(int(label), set()).add(index)
    return {frozenset(indices) for indices in groups.values()}


class StructureDouble:
    def __init__(self, lattice, species, fractional_coordinates):
        self.lattice = tuple(tuple(row) for row in lattice)
        self.species = tuple(species)
        self.fractional_coordinates = tuple(
            tuple(row) for row in fractional_coordinates
        )
        self.pbc = (True, True, True)


def _structure(item) -> StructureDouble:
    return StructureDouble(
        item["lattice"], item["species"], item["fractional_coordinates"]
    )


def _vasp_text(item: dict) -> str:
    species_order = tuple(dict.fromkeys(item["species"]))
    counts = [item["species"].count(symbol) for symbol in species_order]
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
class WyckoffOrbitTests(unittest.TestCase):
    def test_matches_pinned_independent_moyo_results(self) -> None:
        source = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
        independent = json.loads(MOYO_FIXTURE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(independent["upstream"]["version"], "0.10.0")
        self.assertEqual(
            hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest(),
            independent["input_fixture"]["sha256"],
        )
        expected_by_id = {item["id"]: item for item in independent["structures"]}
        for item in source["structures"]:
            with self.subTest(item=item["id"]):
                expected = expected_by_id[item["id"]]
                analysis = analyze_wyckoff_orbits(
                    _structure(item),
                    symprec=independent["classification"]["symprec"],
                )
                self.assertEqual(
                    analysis.symmetry.space_group.ita_number,
                    expected["ita_number"],
                )
                self.assertEqual(
                    analysis.symmetry.hall_setting.hall_number,
                    expected["hall_number"],
                )
                self.assertEqual(
                    analysis.symmetry.wyckoff_letters,
                    tuple(expected["wyckoff_letters"]),
                )
                self.assertEqual(
                    analysis.symmetry.site_symmetry_symbols,
                    tuple(expected["site_symmetry_symbols"]),
                )
                self.assertEqual(
                    _partition(analysis.symmetry.crystallographic_orbits),
                    _partition(expected["orbit_labels"]),
                )

    def test_all_real_fixtures_obey_both_orbit_stabilizer_relations(self) -> None:
        fixture = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
        settings = fixture["classification"]
        for item in fixture["structures"]:
            with self.subTest(item=item["id"]):
                analysis = analyze_wyckoff_orbits(
                    _structure(item),
                    symprec=settings["symprec"],
                    angle_tolerance=settings["angle_tolerance"],
                )
                self.assertIsInstance(analysis, WyckoffOrbitAnalysis)
                self.assertEqual(
                    len(analysis.orbits),
                    item["expected"]["equivalent_atom_orbit_count"],
                )
                covered = sorted(
                    index
                    for orbit in analysis.orbits
                    for index in orbit.site_indices
                )
                self.assertEqual(covered, list(range(len(item["species"]))))
                for orbit in analysis.orbits:
                    self.assertEqual(
                        orbit.multiplicity * orbit.stabilizer_order,
                        len(analysis.symmetry.standard_operations),
                    )
                    for operation in orbit.stabilizer_operations:
                        lattice = np.asarray(analysis.symmetry.standardized_lattice)
                        basis = lattice.T
                        rotation = basis @ operation.rotation @ np.linalg.inv(basis)
                        for vector in orbit.allowed_displacement_basis_cartesian:
                            self.assertTrue(
                                np.allclose(
                                    rotation @ np.asarray(vector),
                                    vector,
                                    atol=1.0e-8,
                                    rtol=0.0,
                                )
                            )

    def test_nacl_special_sites_have_zero_symmetry_preserving_displacement(self) -> None:
        structure = StructureDouble(
            [[0.0, 2.82, 2.82], [2.82, 0.0, 2.82], [2.82, 2.82, 0.0]],
            ["Na", "Cl"],
            [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        )
        analysis = analyze_wyckoff_orbits(structure)
        by_species = {orbit.species: orbit for orbit in analysis.orbits}
        self.assertEqual(by_species["Na"].label, "4a")
        self.assertEqual(by_species["Cl"].label, "4b")
        for orbit in by_species.values():
            self.assertEqual(orbit.site_symmetry_symbol, "m-3m")
            self.assertEqual(orbit.multiplicity, 4)
            self.assertEqual(orbit.stabilizer_order, 48)
            self.assertEqual(orbit.allowed_displacement_dimension, 0)

    def test_general_position_allows_three_cartesian_displacements(self) -> None:
        item = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))["structures"][0]
        analysis = analyze_wyckoff_orbits(_structure(item))
        self.assertTrue(
            all(orbit.allowed_displacement_dimension == 3 for orbit in analysis.orbits)
        )
        payload = analysis.to_dict()
        self.assertEqual(len(payload["orbits"]), len(analysis.orbits))
        self.assertIn("displacement_coordinate_system", payload)

    def test_mirror_site_allows_only_displacements_in_the_mirror_plane(self) -> None:
        item = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))["structures"][2]
        analysis = analyze_wyckoff_orbits(_structure(item))
        mirror_orbits = [
            orbit for orbit in analysis.orbits if orbit.site_symmetry_symbol == ".m."
        ]
        self.assertTrue(mirror_orbits)
        for orbit in mirror_orbits:
            self.assertEqual(orbit.allowed_displacement_dimension, 2)
            self.assertTrue(
                all(abs(vector[1]) < 1.0e-10 for vector in orbit.allowed_displacement_basis_cartesian)
            )

    def test_cli_reads_poscar_and_returns_wyckoff_orbits(self) -> None:
        item = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))["structures"][-2]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(_vasp_text(item), encoding="utf-8")
            output = io.StringIO()
            with redirect_stdout(output):
                self.assertEqual(
                    main(
                        [
                            "wyckoff-orbits",
                            str(path),
                            "--input-format",
                            "vasp",
                            "--json",
                        ]
                    ),
                    0,
                )
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["space_group"]["ita_number"], 183)
        self.assertEqual([orbit["label"] for orbit in payload["orbits"]], ["1a"] * 3)
        self.assertEqual(
            [orbit["allowed_displacement_dimension"] for orbit in payload["orbits"]],
            [1, 1, 1],
        )


if __name__ == "__main__":
    unittest.main()
