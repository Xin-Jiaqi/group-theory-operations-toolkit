"""Command-line interface for the operation catalog."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from . import __version__

from .catalog import (
    GroupDataError,
    find_operations,
    get_layer_group,
    get_operation,
    get_point_group,
    load_database,
    multiply_operations,
    operation_record,
    validate_database,
)
from .representations import quadratic_field_representation
from .point_groups import (
    get_crystallographic_point_group,
    iter_crystallographic_point_groups,
    point_group_operations,
)
from .invariants import RESPONSE_SPECS, response_tensor_basis
from .structure import apply_fractional_operation


def _list(args: argparse.Namespace, database: dict) -> int:
    families = [args.family] if args.family else list(database["families"])
    payload = []
    for name in families:
        family = database["families"].get(name)
        if family is None:
            raise GroupDataError(f"unknown family {name!r}")
        payload.append(
            {
                "family": name,
                "parent_point_group": family["parent_point_group"],
                "operation_count": len(family["operations"]),
                "point_group_count": len(family["point_groups"]),
                "layer_group_count": len(family["layer_groups"]),
                "operations": [op["name"] for op in family["operations"]],
            }
        )
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        for item in payload:
            print(
                f"{item['family']}: parent={item['parent_point_group']}, "
                f"operations={item['operation_count']}, point_groups={item['point_group_count']}, "
                f"layer_groups={item['layer_group_count']}"
            )
    return 0


def _show(args: argparse.Namespace, database: dict) -> int:
    matches = find_operations(database, args.name, args.family)
    if not matches:
        raise GroupDataError(f"operation {args.name!r} not found")
    records = [operation_record(database, operation["name"], family) for family, operation in matches]
    if args.json:
        print(json.dumps([record.to_dict() for record in records], ensure_ascii=False, indent=2))
    else:
        for record in records:
            print(f"{record.family} / {record.name}")
            for row in record.matrix(args.coordinate):
                print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _group(args: argparse.Namespace, database: dict) -> int:
    result = (
        get_layer_group(database, args.family, args.lg)
        if args.lg is not None
        else get_point_group(database, args.family, args.point_group)
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


def _multiply(args: argparse.Namespace, database: dict) -> int:
    left = get_operation(database, args.left, args.family)["name"]
    right = get_operation(database, args.right, args.family)["name"]
    result = multiply_operations(database, args.family, left, right)
    print(json.dumps({"family": args.family, "left": left, "right": right, "result": result}) if args.json else f"{left} * {right} = {result}")
    return 0


def _field(args: argparse.Namespace, database: dict) -> int:
    matches = find_operations(database, args.name, args.family)
    if not matches:
        raise GroupDataError(f"operation {args.name!r} not found")
    payload = []
    for family, operation in matches:
        representation = quadratic_field_representation(operation)
        item = {
            "family": family,
            "index": operation["index"],
            "name": operation["name"],
            "determinant": representation.determinant,
        }
        if args.space in {"symmetric", "both"}:
            item["matrix_symmetric"] = [
                list(row) for row in representation.matrix_symmetric
            ]
        if args.space in {"antisymmetric", "both"}:
            item["matrix_antisymmetric"] = [
                list(row) for row in representation.matrix_antisymmetric
            ]
        payload.append(item)
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    for item in payload:
        print(f"{item['family']} / {item['name']} / det(D)={item['determinant']:.12g}")
        for key in ("matrix_symmetric", "matrix_antisymmetric"):
            if key not in item:
                continue
            print(key)
            for row in item[key]:
                print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _point_groups(args: argparse.Namespace, database: dict) -> int:
    if args.name is None:
        records = list(iter_crystallographic_point_groups())
        payload: dict[str, Any] | list[dict[str, Any]] = [record.to_dict() for record in records]
    else:
        record = get_crystallographic_point_group(args.name)
        payload = record.to_dict()
        payload["operation_records"] = [
            operation.to_dict()
            for operation in point_group_operations(record.number, database=database)
        ]
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    records = (
        list(iter_crystallographic_point_groups())
        if args.name is None
        else [get_crystallographic_point_group(args.name)]
    )
    for record in records:
        print(
            f"{record.number:2d} {record.hm_symbol:6s} {record.schoenflies_symbol:4s} "
            f"order={record.order:2d} family={record.host_family}"
        )
    return 0


def _invariants(args: argparse.Namespace, database: dict) -> int:
    result = response_tensor_basis(args.point_group, args.response, database=database)
    payload = {
        "point_group_number": result.point_group_number,
        "point_group": result.point_group,
        "schoenflies_symbol": result.schoenflies_symbol,
        **result.to_dict(),
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    print(
        f"{result.point_group} ({result.schoenflies_symbol}) / {result.response} / "
        f"shape={result.shape[0]}x{result.shape[1]} / dimension={result.dimension}"
    )
    for index, matrix in enumerate(result.basis, start=1):
        print(f"basis[{index}]")
        for row in matrix:
            print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _validate(_: argparse.Namespace, database: dict) -> int:
    errors = validate_database(database)
    if errors:
        raise GroupDataError("\n".join(errors))
    print("catalog valid")
    return 0


def _apply(args: argparse.Namespace, database: dict) -> int:
    try:
        from materials_structure_core import (  # type: ignore[import-not-found]
            StructureIOError,
            read_structure,
            write_structure,
        )
    except ImportError as exc:
        raise GroupDataError(
            "apply-structure requires materials-structure-core[io]>=0.0.2"
        ) from exc
    try:
        parsed = read_structure(args.input, format=args.input_format)
        operation = operation_record(database, args.operation, args.family)
        transformed = apply_fractional_operation(
            parsed.structure, operation, wrap=not args.no_wrap
        )
        output = Path(args.output)
        write_structure(
            transformed,
            output,
            format=args.output_format,
            direct=not args.cartesian,
            overwrite=args.force,
        )
    except StructureIOError as exc:
        raise GroupDataError(str(exc)) from exc
    print(output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="group-ops")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--database", help="override the packaged schema-v1 JSON")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list")
    list_parser.add_argument("--family")
    list_parser.add_argument("--json", action="store_true")
    list_parser.set_defaults(handler=_list)

    show_parser = subparsers.add_parser("show")
    show_parser.add_argument("name")
    show_parser.add_argument("--family")
    show_parser.add_argument("--coordinate", choices=("fractional", "cartesian"), default="cartesian")
    show_parser.add_argument("--json", action="store_true")
    show_parser.set_defaults(handler=_show)

    group_parser = subparsers.add_parser("group")
    group_parser.add_argument("--family", required=True)
    selection = group_parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--lg", type=int)
    selection.add_argument("--point-group")
    group_parser.set_defaults(handler=_group)

    multiply_parser = subparsers.add_parser("multiply")
    multiply_parser.add_argument("left")
    multiply_parser.add_argument("right")
    multiply_parser.add_argument("--family", required=True)
    multiply_parser.add_argument("--json", action="store_true")
    multiply_parser.set_defaults(handler=_multiply)

    field_parser = subparsers.add_parser("field-representation")
    field_parser.add_argument("name")
    field_parser.add_argument("--family", required=True)
    field_parser.add_argument(
        "--space",
        choices=("symmetric", "antisymmetric", "both"),
        default="both",
    )
    field_parser.add_argument("--json", action="store_true")
    field_parser.set_defaults(handler=_field)

    point_groups_parser = subparsers.add_parser("point-groups")
    point_groups_parser.add_argument("name", nargs="?")
    point_groups_parser.add_argument("--json", action="store_true")
    point_groups_parser.set_defaults(handler=_point_groups)

    invariants_parser = subparsers.add_parser("invariants")
    invariants_parser.add_argument("point_group")
    invariants_parser.add_argument("response", choices=tuple(RESPONSE_SPECS))
    invariants_parser.add_argument("--json", action="store_true")
    invariants_parser.set_defaults(handler=_invariants)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.set_defaults(handler=_validate)

    apply_parser = subparsers.add_parser("apply-structure")
    apply_parser.add_argument("input")
    apply_parser.add_argument("operation")
    apply_parser.add_argument("--family", required=True)
    apply_parser.add_argument("--output", required=True)
    apply_parser.add_argument("--input-format", choices=("vasp", "cif"))
    apply_parser.add_argument("--output-format", choices=("vasp", "cif"))
    apply_parser.add_argument("--cartesian", action="store_true")
    apply_parser.add_argument("--no-wrap", action="store_true")
    apply_parser.add_argument("--force", action="store_true")
    apply_parser.set_defaults(handler=_apply)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        database = load_database(args.database)
        return int(args.handler(args, database))
    except GroupDataError as exc:
        parser.error(str(exc))
    return 2
