#!/usr/bin/env python3
"""Query point-group data and apply a point operation to POSCAR coordinates."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
DEFAULT_DATABASE = ROOT / "data" / "group_operations.json"


class GroupDataError(ValueError):
    """Raised when the database or a user selection is inconsistent."""


def canonical_name(value: str) -> str:
    """Normalize common historical and LaTeX-like operation labels."""
    name = (
        value.strip()
        .replace("−", "-")
        .replace("^", "")
        .replace("{", "")
        .replace("}", "")
        .replace(" ", "")
    )
    if name.lower() == "i":
        return "-1"
    if name in {"1", "-1"} or "_" in name:
        return name
    if name.startswith("m"):
        return "m_" + name[1:]
    for prefix in (
        "-6+", "-6-", "-4+", "-4-", "-3+", "-3-",
        "6+", "6-", "4+", "4-", "3+", "3-", "2",
    ):
        if name.startswith(prefix) and len(name) > len(prefix):
            return prefix + "_" + name[len(prefix):]
    return name


def load_database(path: str | Path = DEFAULT_DATABASE) -> dict:
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise GroupDataError(f"无法读取数据文件：{path}") from exc
    except json.JSONDecodeError as exc:
        raise GroupDataError(f"JSON 格式错误：{path}: {exc}") from exc
    if data.get("schema_version") != 1 or not isinstance(data.get("families"), dict):
        raise GroupDataError("不支持的数据结构或 schema_version")
    return data


def family_data(database: dict, family: str) -> dict:
    try:
        return database["families"][family]
    except KeyError as exc:
        choices = ", ".join(database["families"])
        raise GroupDataError(f"未知晶系数据集 {family!r}；可选：{choices}") from exc


def iter_operations(database: dict, family: str | None = None) -> Iterable[tuple[str, dict]]:
    families = [family] if family else database["families"]
    for family_name in families:
        for operation in family_data(database, family_name)["operations"]:
            yield family_name, operation


def find_operations(database: dict, name: str, family: str | None = None) -> list[tuple[str, dict]]:
    target = canonical_name(name)
    matches = []
    for family_name, operation in iter_operations(database, family):
        aliases = {operation["name"]}
        if "source_name" in operation:
            aliases.add(canonical_name(operation["source_name"]))
        if target in aliases:
            matches.append((family_name, operation))
    return matches


def get_operation(database: dict, name: str, family: str) -> dict:
    matches = find_operations(database, name, family)
    if not matches:
        raise GroupDataError(f"在 {family} 中找不到点操作 {name!r}")
    if len(matches) != 1:
        raise GroupDataError(f"点操作 {name!r} 在 {family} 中不是唯一项")
    return matches[0][1]


def get_layer_group(database: dict, family: str, lg_number: int) -> dict:
    for entry in family_data(database, family)["layer_groups"]:
        if entry["LG"] == lg_number:
            return entry
    raise GroupDataError(f"在 {family} 中找不到 LG{lg_number}")


def get_point_group(database: dict, family: str, name: str) -> dict:
    target = name.lower()
    for entry in family_data(database, family)["point_groups"]:
        if entry["name"].lower() == target:
            return entry
    raise GroupDataError(f"在 {family} 中找不到点群 {name!r}")


def multiply_operations(database: dict, family: str, left: str, right: str) -> str:
    """Return the canonical name of left * right from the stored Cayley table."""
    family_entry = family_data(database, family)
    multiplication = family_entry.get("multiplication")
    if multiplication is None:
        raise GroupDataError(f"{family} 尚未收录乘法表")
    left_name = get_operation(database, left, family)["name"]
    right_name = get_operation(database, right, family)["name"]
    try:
        return multiplication["table"][left_name][right_name]
    except KeyError as exc:
        raise GroupDataError(f"乘法表缺少 {left_name} * {right_name}") from exc


def matrix_for(operation: dict, coordinate: str) -> list[list[float]]:
    key = "matrix_fractional" if coordinate == "fractional" else "matrix_cartesian"
    try:
        return [[float(value) for value in row] for row in operation[key]]
    except (KeyError, TypeError, ValueError) as exc:
        raise GroupDataError(f"操作 {operation.get('name')} 缺少可用的 {coordinate} 矩阵") from exc


def multiply_matrix_vector(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def parse_poscar(path: str | Path) -> dict:
    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise GroupDataError(f"无法读取 POSCAR：{path}") from exc
    if len(lines) < 8:
        raise GroupDataError("POSCAR 内容过短")

    try:
        atom_count = sum(int(value) for value in lines[6].split())
    except ValueError as exc:
        raise GroupDataError("POSCAR 第 7 行不是有效的元素数目行（当前仅支持 VASP 5 格式）") from exc

    cursor = 7
    selective = None
    if lines[cursor].strip().lower().startswith("s"):
        selective = lines[cursor]
        cursor += 1
    if cursor >= len(lines):
        raise GroupDataError("POSCAR 缺少坐标类型行")

    coordinate_line = lines[cursor]
    token = coordinate_line.strip().split()[0].lower()
    if token.startswith("d"):
        coordinate_type = "Direct"
    elif token.startswith("c") or token.startswith("k"):
        coordinate_type = "Cartesian"
    else:
        raise GroupDataError(f"未知 POSCAR 坐标类型：{coordinate_line}")
    cursor += 1

    if len(lines) < cursor + atom_count:
        raise GroupDataError("POSCAR 原子坐标行数不足")
    coordinates = []
    flags = []
    for line in lines[cursor:cursor + atom_count]:
        fields = line.split()
        if len(fields) < 3:
            raise GroupDataError(f"无效坐标行：{line}")
        try:
            coordinates.append([float(fields[0]), float(fields[1]), float(fields[2])])
        except ValueError as exc:
            raise GroupDataError(f"无效坐标行：{line}") from exc
        flags.append(fields[3:])

    return {
        "title": lines[0],
        "scale": lines[1],
        "lattice": lines[2:5],
        "species": lines[5],
        "counts": lines[6],
        "selective": selective,
        "coordinate_line": coordinate_line,
        "coordinate_type": coordinate_type,
        "coordinates": coordinates,
        "flags": flags,
        "tail": lines[cursor + atom_count:],
    }


def transform_coordinates(
    coordinates: list[list[float]], matrix: list[list[float]], direct: bool
) -> list[list[float]]:
    transformed = []
    for coordinate in coordinates:
        result = multiply_matrix_vector(matrix, coordinate)
        if direct:
            result = [value % 1.0 for value in result]
            result = [0.0 if math.isclose(value, 0.0, abs_tol=1e-12) or math.isclose(value, 1.0, abs_tol=1e-12) else value for value in result]
        transformed.append(result)
    return transformed


def format_number(value: float) -> str:
    if math.isclose(value, 0.0, abs_tol=1e-14):
        value = 0.0
    return f"{value:.16f}".rstrip("0").rstrip(".") or "0"


def write_poscar(path: str | Path, poscar: dict, coordinates: list[list[float]], force: bool = False) -> Path:
    path = Path(path)
    if path.exists() and not force:
        raise GroupDataError(f"输出文件已存在：{path}；如需覆盖请使用 --force")
    lines = [
        poscar["title"], poscar["scale"], *poscar["lattice"],
        poscar["species"], poscar["counts"],
    ]
    if poscar["selective"] is not None:
        lines.append(poscar["selective"])
    lines.append(poscar["coordinate_line"])
    for coordinate, flags in zip(coordinates, poscar["flags"]):
        line = "  ".join(format_number(value) for value in coordinate)
        if flags:
            line += "  " + "  ".join(flags)
        lines.append(line)
    lines.extend(poscar["tail"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def print_matrix(matrix: list[list[float]]) -> None:
    for row in matrix:
        print("[" + ", ".join(format_number(value) for value in row) + "]")


def command_list(args: argparse.Namespace, database: dict) -> None:
    if args.family:
        family = family_data(database, args.family)
        print(f"{args.family}: {family['parent_point_group']} ({len(family['operations'])} operations)")
        for operation in family["operations"]:
            alias = f" [source: {operation['source_name']}]" if "source_name" in operation else ""
            print(f"  {operation['index']:>2}  {operation['name']}{alias}")
        return
    for name, family in database["families"].items():
        multiplication = family.get("multiplication")
        table_size = len(multiplication["element_order"]) if multiplication else 0
        print(
            f"{name}: parent={family['parent_point_group']}, "
            f"operations={len(family['operations'])}, "
            f"point_groups={len(family['point_groups'])}, "
            f"layer_groups={len(family['layer_groups'])}, "
            f"multiplication_table={table_size}x{table_size if table_size else 0}"
        )


def command_show(args: argparse.Namespace, database: dict) -> None:
    matches = find_operations(database, args.name, args.family)
    if not matches:
        raise GroupDataError(f"找不到点操作 {args.name!r}")
    for family_name, operation in matches:
        print(f"{family_name} / {operation['name']}")
        print_matrix(matrix_for(operation, args.coordinate))


def command_group(args: argparse.Namespace, database: dict) -> None:
    if args.lg is not None:
        print(json.dumps(get_layer_group(database, args.family, args.lg), ensure_ascii=False, indent=2))
    else:
        print(json.dumps(get_point_group(database, args.family, args.point_group), ensure_ascii=False, indent=2))


def command_multiply(args: argparse.Namespace, database: dict) -> None:
    left = get_operation(database, args.left, args.family)["name"]
    right = get_operation(database, args.right, args.family)["name"]
    result = multiply_operations(database, args.family, left, right)
    print(f"{left} * {right} = {result}")


def command_apply(args: argparse.Namespace, database: dict) -> None:
    poscar = parse_poscar(args.input)
    operation = get_operation(database, args.operation, args.family)
    direct = poscar["coordinate_type"] == "Direct"
    matrix = matrix_for(operation, "fractional" if direct else "cartesian")
    coordinates = transform_coordinates(poscar["coordinates"], matrix, direct)
    output = Path(args.output) if args.output else Path(f"{operation['name']}_{Path(args.input).stem}_out.vasp")
    write_poscar(output, poscar, coordinates, force=args.force)
    print(f"family: {args.family}")
    print(f"operation: {operation['name']}")
    print(f"coordinate type: {poscar['coordinate_type']}")
    print(f"output: {output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", default=str(DEFAULT_DATABASE), help="统一 JSON 数据文件路径")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="列出数据集或操作")
    list_parser.add_argument("--family", help="限定数据集")
    list_parser.set_defaults(handler=command_list)

    show_parser = subparsers.add_parser("show", help="显示一个操作的矩阵")
    show_parser.add_argument("name")
    show_parser.add_argument("--family")
    show_parser.add_argument("--coordinate", choices=("fractional", "cartesian"), default="cartesian")
    show_parser.set_defaults(handler=command_show)

    group_parser = subparsers.add_parser("group", help="查询点群或层群的操作集合")
    group_parser.add_argument("--family", required=True)
    selection = group_parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--lg", type=int)
    selection.add_argument("--point-group")
    group_parser.set_defaults(handler=command_group)

    multiply_parser = subparsers.add_parser("multiply", help="查询两个点操作的乘积")
    multiply_parser.add_argument("left", help="左侧（行）操作")
    multiply_parser.add_argument("right", help="右侧（列）操作")
    multiply_parser.add_argument("--family", required=True)
    multiply_parser.set_defaults(handler=command_multiply)

    apply_parser = subparsers.add_parser("apply-poscar", help="对 POSCAR 坐标施加点操作")
    apply_parser.add_argument("input")
    apply_parser.add_argument("operation")
    apply_parser.add_argument("--family", required=True)
    apply_parser.add_argument("--output")
    apply_parser.add_argument("--force", action="store_true")
    apply_parser.set_defaults(handler=command_apply)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        database = load_database(args.database)
        args.handler(args, database)
    except GroupDataError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    sys.exit(main())
