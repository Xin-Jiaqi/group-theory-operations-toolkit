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
from .magnetic_point_groups import (
    MAGNETIC_CATEGORIES,
    get_magnetic_point_group,
    iter_magnetic_point_groups,
    load_magnetic_point_group_registry,
    magnetic_point_group_operations,
)
from .magnetic_layer_groups import (
    MAGNETIC_LAYER_TYPES,
    get_magnetic_layer_group,
    iter_magnetic_layer_groups,
    load_magnetic_layer_group_registry,
)
from .space_groups import (
    get_crystallographic_space_group,
    iter_crystallographic_space_groups,
    load_space_group_registry,
)
from .layer_groups import (
    get_crystallographic_layer_group,
    iter_crystallographic_layer_groups,
    load_layer_group_registry,
)
from .stacking import (
    BRAVAIS_LATTICE_POINT_GROUPS,
    layer_group_polarization,
    stacking_rotation_cosets,
)
from .invariants import (
    MAGNETIC_RESPONSE_SPECS,
    RESPONSE_SYMMETRY_CLASSES,
    RESPONSE_SPECS,
    TENSOR_SPACE_BASES,
    magnetic_layer_response_tensor_basis,
    magnetic_layer_tensor_basis,
    magnetic_tensor_basis,
    magnetic_response_tensor_basis,
    response_tensor_basis,
    screen_response_symmetry,
)
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


def _magnetic_point_groups(args: argparse.Namespace, database: dict) -> int:
    if args.identifier is None:
        records = [
            record
            for record in iter_magnetic_point_groups()
            if args.category is None or record.category == args.category
        ]
        payload: dict[str, Any] | list[dict[str, Any]] = [
            record.to_dict() for record in records
        ]
    else:
        record = get_magnetic_point_group(args.identifier)
        if args.category is not None and record.category != args.category:
            raise GroupDataError(
                f"magnetic point group {record.hm_symbol} is {record.category}, not {args.category}"
            )
        payload = record.to_dict()
        payload["operation_records"] = [
            {
                "name": operation.name,
                "label": operation.label,
                "time_reversal": operation.time_reversal,
                "matrix_cartesian": [list(row) for row in operation.spatial.matrix_cartesian],
            }
            for operation in magnetic_point_group_operations(
                record.number, database=database
            )
        ]
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    records = (
        [
            record
            for record in iter_magnetic_point_groups()
            if args.category is None or record.category == args.category
        ]
        if args.identifier is None
        else [get_magnetic_point_group(args.identifier)]
    )
    for record in records:
        print(
            f"{record.magnetic_number:9s} {record.hm_symbol:12s} "
            f"{record.category:21s} order={record.order:2d} "
            f"parent={record.parent_point_group_hm}"
        )
    return 0


def _space_groups(args: argparse.Namespace, database: dict) -> int:
    if args.ita is None:
        records = list(iter_crystallographic_space_groups())
        payload: dict[str, Any] | list[dict[str, Any]] = [record.to_dict() for record in records]
    else:
        record = get_crystallographic_space_group(args.ita)
        payload = record.to_dict()
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    records = (
        list(iter_crystallographic_space_groups())
        if args.ita is None
        else [get_crystallographic_space_group(args.ita)]
    )
    for record in records:
        print(
            f"{record.ita_number:3d} {record.international_short:8s} "
            f"{record.schoenflies:7s} {record.crystal_system:12s} "
            f"center={record.centering} point={record.point_group_hm:5s} "
            f"{'symmorphic' if record.symmorphic else 'non-symmorphic'}"
        )
    return 0


def _layer_groups(args: argparse.Namespace, database: dict) -> int:
    if args.identifier is None:
        records = list(iter_crystallographic_layer_groups())
        payload: dict[str, Any] | list[dict[str, Any]] = [
            record.to_dict() for record in records
        ]
    else:
        record = get_crystallographic_layer_group(args.identifier)
        payload = record.to_dict()
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    records = (
        list(iter_crystallographic_layer_groups())
        if args.identifier is None
        else [get_crystallographic_layer_group(args.identifier)]
    )
    for record in records:
        print(
            f"{record.number:2d} {record.international_short:9s} "
            f"{record.schoenflies:7s} {record.crystal_system:12s} "
            f"center={record.centering} point={record.point_group_hm:5s} "
            f"settings={len(record.hall_settings)}"
        )
    return 0


def _layer_polarity(args: argparse.Namespace, database: dict) -> int:
    del database
    group = get_crystallographic_layer_group(args.identifier)
    result = layer_group_polarization(
        group.number, layer_hall_number=args.layer_hall_number
    )
    payload = {
        "layer_group_number": group.number,
        "layer_group": group.international_short,
        **result.to_dict(),
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(
            f"LG{group.number} {group.international_short} / "
            f"{result.polar_type} / dimension={result.dimension}"
        )
    return 0


def _stacking_rotations(args: argparse.Namespace, database: dict) -> int:
    del database
    cosets = stacking_rotation_cosets(args.point_group, args.lattice)
    payload = {
        "monolayer_point_group": args.point_group,
        "bravais_lattice": args.lattice,
        "bravais_point_group": BRAVAIS_LATTICE_POINT_GROUPS[args.lattice],
        "rotation_class_count": len(cosets),
        "left_cosets": [coset.to_dict() for coset in cosets],
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(
            f"{args.point_group} in {args.lattice} "
            f"({payload['bravais_point_group']}) / "
            f"rotation classes={len(cosets)}"
        )
        for index, coset in enumerate(cosets, start=1):
            print(
                f"C{index}: order={len(coset.members)} "
                f"representative={coset.representative}"
            )
    return 0


def _magnetic_layer_groups(args: argparse.Namespace, database: dict) -> int:
    del database
    if args.identifier is None:
        records = [
            record
            for record in iter_magnetic_layer_groups()
            if args.type is None or record.magnetic_type == args.type
        ]
        payload: dict[str, Any] | list[dict[str, Any]] = [
            record.to_dict() for record in records
        ]
    else:
        record = get_magnetic_layer_group(args.identifier)
        if args.type is not None and record.magnetic_type != args.type:
            raise GroupDataError(
                f"magnetic layer group {record.og_number} is type "
                f"{record.magnetic_type}, not {args.type}"
            )
        payload = record.to_dict()
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    records = (
        [
            record
            for record in iter_magnetic_layer_groups()
            if args.type is None or record.magnetic_type == args.type
        ]
        if args.identifier is None
        else [get_magnetic_layer_group(args.identifier)]
    )
    for record in records:
        anti = (
            "-"
            if record.anti_translation_fractional is None
            else ",".join(f"{value:g}" for value in record.anti_translation_fractional)
        )
        print(
            f"{record.og_number:10s} {record.litvin_og_symbol_ascii:18s} "
            f"type={record.magnetic_type:3s} parent=LG{record.parent_layer_group_number:02d} "
            f"point={record.magnetic_point_group_symbol:10s} anti-t={anti}"
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


def _screen_responses(args: argparse.Namespace, database: dict) -> int:
    results = screen_response_symmetry(
        args.symmetry_class,
        groups=args.groups,
        responses=args.responses,
        allowed_only=args.allowed_only,
        database=database,
    )
    if args.json:
        print(
            json.dumps(
                [result.to_dict() for result in results],
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    for result in results:
        time_character = (
            "spatial" if result.time_character is None else result.time_character
        )
        print(
            f"{result.group_identifier:10s} {result.group_symbol:18s} "
            f"{result.response:28s} {time_character:7s} "
            f"shape={result.shape[0]}x{result.shape[1]} "
            f"dimension={result.dimension:2d} "
            f"allowed={'yes' if result.allowed else 'no'}"
        )
    return 0


def _magnetic_invariants(args: argparse.Namespace, database: dict) -> int:
    result = magnetic_tensor_basis(
        args.magnetic_point_group,
        args.output_space,
        args.input_space,
        output_time_parity=args.output_time,
        input_time_parity=args.input_time,
        database=database,
    )
    if args.json:
        print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
        return 0
    print(
        f"{result.magnetic_point_group} ({result.magnetic_number}) / "
        f"{result.output_space}[T-{result.output_time_parity}] <- "
        f"{result.input_space}[T-{result.input_time_parity}] / "
        f"shape={result.shape[0]}x{result.shape[1]} / dimension={result.dimension}"
    )
    for index, matrix in enumerate(result.basis, start=1):
        print(f"basis[{index}]")
        for row in matrix:
            print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _magnetic_responses(args: argparse.Namespace, database: dict) -> int:
    result = magnetic_response_tensor_basis(
        args.magnetic_point_group,
        args.response,
        database=database,
    )
    if args.json:
        print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
        return 0
    print(
        f"{result.magnetic_point_group} ({result.magnetic_number}) / "
        f"{result.response}[T-{result.time_character}] / "
        f"shape={result.shape[0]}x{result.shape[1]} / dimension={result.dimension}"
    )
    for index, matrix in enumerate(result.basis, start=1):
        print(f"basis[{index}]")
        for row in matrix:
            print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _magnetic_layer_invariants(args: argparse.Namespace, database: dict) -> int:
    del database
    result = magnetic_layer_tensor_basis(
        args.magnetic_layer_group,
        args.output_space,
        args.input_space,
        output_time_parity=args.output_time,
        input_time_parity=args.input_time,
    )
    if args.json:
        print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
        return 0
    print(
        f"{result.magnetic_layer_group} ({result.og_number}) / "
        f"{result.output_space}[T-{result.output_time_parity}] <- "
        f"{result.input_space}[T-{result.input_time_parity}] / "
        f"shape={result.shape[0]}x{result.shape[1]} / dimension={result.dimension}"
    )
    for index, matrix in enumerate(result.basis, start=1):
        print(f"basis[{index}]")
        for row in matrix:
            print("[" + ", ".join(f"{value:.12g}" for value in row) + "]")
    return 0


def _magnetic_layer_responses(args: argparse.Namespace, database: dict) -> int:
    del database
    result = magnetic_layer_response_tensor_basis(
        args.magnetic_layer_group,
        args.response,
    )
    if args.json:
        print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
        return 0
    print(
        f"{result.magnetic_layer_group} ({result.og_number}) / "
        f"{result.response}[T-{result.time_character}] / "
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
    load_space_group_registry()
    load_layer_group_registry()
    load_magnetic_point_group_registry()
    load_magnetic_layer_group_registry()
    print("catalog valid")
    return 0


def _apply(args: argparse.Namespace, database: dict) -> int:
    try:
        from materials_structure_core import (
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

    magnetic_point_groups_parser = subparsers.add_parser("magnetic-point-groups")
    magnetic_point_groups_parser.add_argument("identifier", nargs="?")
    magnetic_point_groups_parser.add_argument("--category", choices=MAGNETIC_CATEGORIES)
    magnetic_point_groups_parser.add_argument("--json", action="store_true")
    magnetic_point_groups_parser.set_defaults(handler=_magnetic_point_groups)

    space_groups_parser = subparsers.add_parser("space-groups")
    space_groups_parser.add_argument("ita", type=int, nargs="?")
    space_groups_parser.add_argument("--json", action="store_true")
    space_groups_parser.set_defaults(handler=_space_groups)

    layer_groups_parser = subparsers.add_parser("layer-groups")
    layer_groups_parser.add_argument("identifier", nargs="?")
    layer_groups_parser.add_argument("--json", action="store_true")
    layer_groups_parser.set_defaults(handler=_layer_groups)

    layer_polarity_parser = subparsers.add_parser("layer-polarity")
    layer_polarity_parser.add_argument("identifier")
    layer_polarity_parser.add_argument("--layer-hall-number", type=int)
    layer_polarity_parser.add_argument("--json", action="store_true")
    layer_polarity_parser.set_defaults(handler=_layer_polarity)

    stacking_rotations_parser = subparsers.add_parser("stacking-rotations")
    stacking_rotations_parser.add_argument("point_group")
    stacking_rotations_parser.add_argument(
        "--lattice",
        choices=tuple(BRAVAIS_LATTICE_POINT_GROUPS),
        required=True,
    )
    stacking_rotations_parser.add_argument("--json", action="store_true")
    stacking_rotations_parser.set_defaults(handler=_stacking_rotations)

    magnetic_layer_groups_parser = subparsers.add_parser("magnetic-layer-groups")
    magnetic_layer_groups_parser.add_argument("identifier", nargs="?")
    magnetic_layer_groups_parser.add_argument("--type", choices=MAGNETIC_LAYER_TYPES)
    magnetic_layer_groups_parser.add_argument("--json", action="store_true")
    magnetic_layer_groups_parser.set_defaults(handler=_magnetic_layer_groups)

    invariants_parser = subparsers.add_parser("invariants")
    invariants_parser.add_argument("point_group")
    invariants_parser.add_argument("response", choices=tuple(RESPONSE_SPECS))
    invariants_parser.add_argument("--json", action="store_true")
    invariants_parser.set_defaults(handler=_invariants)

    screen_responses_parser = subparsers.add_parser("screen-responses")
    screen_responses_parser.add_argument(
        "--symmetry-class",
        choices=RESPONSE_SYMMETRY_CLASSES,
        required=True,
    )
    screen_responses_parser.add_argument(
        "--group", dest="groups", action="append"
    )
    screen_responses_parser.add_argument(
        "--response", dest="responses", action="append"
    )
    screen_responses_parser.add_argument("--allowed-only", action="store_true")
    screen_responses_parser.add_argument("--json", action="store_true")
    screen_responses_parser.set_defaults(handler=_screen_responses)

    magnetic_invariants_parser = subparsers.add_parser("magnetic-invariants")
    magnetic_invariants_parser.add_argument("magnetic_point_group")
    magnetic_invariants_parser.add_argument("output_space", choices=tuple(TENSOR_SPACE_BASES))
    magnetic_invariants_parser.add_argument("input_space", choices=tuple(TENSOR_SPACE_BASES))
    magnetic_invariants_parser.add_argument("--output-time", choices=("even", "odd"), default="even")
    magnetic_invariants_parser.add_argument("--input-time", choices=("even", "odd"), default="even")
    magnetic_invariants_parser.add_argument("--json", action="store_true")
    magnetic_invariants_parser.set_defaults(handler=_magnetic_invariants)

    magnetic_responses_parser = subparsers.add_parser("magnetic-responses")
    magnetic_responses_parser.add_argument("magnetic_point_group")
    magnetic_responses_parser.add_argument(
        "response", choices=tuple(MAGNETIC_RESPONSE_SPECS)
    )
    magnetic_responses_parser.add_argument("--json", action="store_true")
    magnetic_responses_parser.set_defaults(handler=_magnetic_responses)

    magnetic_layer_invariants_parser = subparsers.add_parser(
        "magnetic-layer-invariants"
    )
    magnetic_layer_invariants_parser.add_argument("magnetic_layer_group")
    magnetic_layer_invariants_parser.add_argument(
        "output_space", choices=tuple(TENSOR_SPACE_BASES)
    )
    magnetic_layer_invariants_parser.add_argument(
        "input_space", choices=tuple(TENSOR_SPACE_BASES)
    )
    magnetic_layer_invariants_parser.add_argument(
        "--output-time", choices=("even", "odd"), default="even"
    )
    magnetic_layer_invariants_parser.add_argument(
        "--input-time", choices=("even", "odd"), default="even"
    )
    magnetic_layer_invariants_parser.add_argument("--json", action="store_true")
    magnetic_layer_invariants_parser.set_defaults(
        handler=_magnetic_layer_invariants
    )

    magnetic_layer_responses_parser = subparsers.add_parser(
        "magnetic-layer-responses"
    )
    magnetic_layer_responses_parser.add_argument("magnetic_layer_group")
    magnetic_layer_responses_parser.add_argument(
        "response", choices=tuple(MAGNETIC_RESPONSE_SPECS)
    )
    magnetic_layer_responses_parser.add_argument("--json", action="store_true")
    magnetic_layer_responses_parser.set_defaults(
        handler=_magnetic_layer_responses
    )

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
