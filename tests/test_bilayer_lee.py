"""Regression tests for the symmetry-only bilayer layer Edelstein effect."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import json
import unittest

from group_theory_operations import (
    bilayer_edelstein_analysis,
    bilayer_edelstein_point_groups,
    partition_bilayer_operations,
    point_group_operations,
)
from group_theory_operations.cli import main


class BilayerEdelsteinTests(unittest.TestCase):
    def test_mosse_inversion_stacking_is_type_i(self) -> None:
        result = bilayer_edelstein_point_groups("C3v", "D3d")
        self.assertEqual(result.classification, "type_i")
        self.assertEqual(result.type_i_components, ("xy", "yx"))
        self.assertEqual(result.type_ii_components, ())
        self.assertEqual(result.zero_field_global_basis, ())

    def test_mote2_direct_stacking_is_type_ii(self) -> None:
        result = bilayer_edelstein_point_groups("D3h", "D3h")
        self.assertEqual(result.classification, "type_ii")
        self.assertEqual(result.type_i_components, ())
        self.assertEqual(result.type_ii_components, ("xy", "yx"))
        self.assertEqual(result.monolayer_basis, ())
        self.assertEqual(result.zero_field_global_basis, ())

    def test_d3_inversion_route_has_both_component_classes(self) -> None:
        result = bilayer_edelstein_point_groups("D3", "Ci")
        self.assertEqual(result.classification, "both")
        self.assertEqual(result.type_i_components, ("xx", "yy", "zz"))
        self.assertEqual(
            result.type_ii_components, ("xy", "xz", "yx", "yz", "zx", "zy")
        )

    def test_no_layer_exchange_operation_is_forbidden(self) -> None:
        result = bilayer_edelstein_point_groups("C3v", "C3v")
        self.assertEqual(result.classification, "forbidden")
        self.assertEqual(result.interlayer_operation_count, 0)

    def test_generic_operation_api_matches_standard_partition(self) -> None:
        monolayer = point_group_operations("3m")
        intra, inter = partition_bilayer_operations(point_group_operations("-3m"))
        result = bilayer_edelstein_analysis(monolayer, intra, inter)
        self.assertEqual(result.classification, "type_i")

    def test_cli_returns_componentwise_json(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["bilayer-edelstein", "C3v", "D3d", "--json"]), 0)
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["classification"], "type_i")
        self.assertEqual(payload["type_i_components"], ["xy", "yx"])
