"""High-throughput symmetry screening for nonlinear optical responses."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
import unittest

from group_theory_operations import GroupDataError, screen_response_symmetry
from group_theory_operations.cli import main


class ResponseSymmetryScreeningTests(unittest.TestCase):
    def test_all_registry_rows_are_reported_in_standard_order(self):
        point = screen_response_symmetry("point_group")
        magnetic_point = screen_response_symmetry("magnetic_point_group")
        magnetic_layer = screen_response_symmetry("magnetic_layer_group")

        self.assertEqual(len(point), 32 * 3)
        self.assertEqual(len(magnetic_point), 122 * 6)
        self.assertEqual(len(magnetic_layer), 528 * 6)
        self.assertEqual(
            [(item.group_number, item.response) for item in point[:3]],
            [
                (1, "shift_current"),
                (1, "shg"),
                (1, "circular_injection_current"),
            ],
        )
        self.assertEqual(magnetic_point[-1].group_number, 122)
        self.assertEqual(magnetic_layer[-1].group_number, 528)

    def test_filters_preserve_physical_group_and_response_order(self):
        results = screen_response_symmetry(
            "point",
            groups=("4mm", "-1", "4mm"),
            responses=("cpge", "shift"),
        )
        self.assertEqual(
            [
                (item.group_symbol, item.response, item.dimension)
                for item in results
            ],
            [
                ("-1", "shift_current", 0),
                ("-1", "circular_injection_current", 0),
                ("4mm", "shift_current", 3),
                ("4mm", "circular_injection_current", 1),
            ],
        )
        allowed = screen_response_symmetry(
            "point_group",
            groups=("4mm", "-1"),
            responses=("shift_current", "circular_injection_current"),
            allowed_only=True,
        )
        self.assertEqual([item.group_symbol for item in allowed], ["4mm", "4mm"])
        self.assertTrue(all(item.allowed for item in allowed))

    def test_gray_group_removes_time_odd_responses(self):
        results = screen_response_symmetry(
            "magnetic_point",
            groups=("11'",),
            responses=("NSC", "MIC", "shg_c"),
        )
        dimensions = {item.response: item.dimension for item in results}
        self.assertEqual(dimensions["normal_shift_current"], 18)
        self.assertEqual(dimensions["magnetic_injection_current"], 0)
        self.assertEqual(dimensions["shg_odd"], 0)
        self.assertEqual(
            {item.time_character for item in results if item.dimension == 0},
            {"odd"},
        )

    def test_invalid_selections_are_rejected(self):
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("unknown")
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("point_group", groups="4mm")
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("point_group", groups=(True,))
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("point_group", responses=())
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("point_group", responses=("unknown",))
        with self.assertRaises(GroupDataError):
            screen_response_symmetry("point_group", allowed_only=1)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            screen_response_symmetry("point_group", tolerance=0.0)

    def test_cli_returns_compact_json_rows(self):
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(
                main(
                    [
                        "screen-responses",
                        "--symmetry-class",
                        "point_group",
                        "--group",
                        "4mm",
                        "--response",
                        "shift_current",
                        "--response",
                        "circular_injection_current",
                        "--allowed-only",
                        "--json",
                    ]
                ),
                0,
            )
        payload = json.loads(output.getvalue())
        self.assertEqual(len(payload), 2)
        self.assertEqual(
            [(item["response"], item["dimension"]) for item in payload],
            [("shift_current", 3), ("circular_injection_current", 1)],
        )
        self.assertTrue(all(item["allowed"] for item in payload))


if __name__ == "__main__":
    unittest.main()
