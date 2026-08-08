#!/usr/bin/env python3
"""Extract a compact, reviewable source table for all 528 magnetic layer groups.

The script combines two primary-source supplements without checking the PDFs
into the repository:

* Litvin (2005), complete magnetic-layer-group tables, supplies the OG symbol,
  operation numbering, time-reversal primes, and type-IV anti-translation.
* Zhang et al. (2023), Table S1/S2, supplies the type, corresponding magnetic
  space group, and basis-change key.

Only point-co-group data are extracted.  Spatial matrices are added later from
the already versioned crystallographic layer-group operation catalog.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any

import numpy as np
import spglib


ROOT = Path(__file__).resolve().parents[1]
OPERATION_PATH = ROOT / "data" / "group_operations.json"
DEFAULT_OUTPUT = (
    ROOT / "scripts" / "sources" / "magnetic-layer-groups-528-v2023.tsv"
)

TYPE_TO_SPGLIB = {"I": 1, "II": 2, "III": 3, "IV": 4}
FIELDS = (
    "og_number",
    "global_number",
    "parent_layer_group_number",
    "family_number",
    "magnetic_type",
    "litvin_og_symbol_ascii",
    "magnetic_point_group_symbol",
    "source_pdf_page",
    "point_operation_names",
    "primed_operation_indices",
    "anti_translation",
    "corresponding_msg_bns_number",
    "corresponding_msg_uni_number",
    "corresponding_msg_og_number",
    "basis_transform",
)


def _pdf_text(path: Path, *, first: int | None = None, last: int | None = None) -> str:
    command = ["pdftotext", "-layout"]
    if first is not None:
        command.extend(["-f", str(first)])
    if last is not None:
        command.extend(["-l", str(last)])
    with tempfile.NamedTemporaryFile(suffix=".txt") as handle:
        command.extend([str(path), handle.name])
        subprocess.run(command, check=True)
        return Path(handle.name).read_text(encoding="utf-8", errors="replace")


def _correspondence_rows(text: str) -> dict[str, dict[str, str]]:
    pattern = re.compile(
        r"(?P<og>\d{1,2}\.\d{1,2}\.\d{1,3})\s+"
        r"(?P<type>IV|III|II|I)\s+"
        r"(?P<symbols>.*?)\s+"
        r"(?P<msg>\d{1,3}\.\d{1,3})\s+[PABCIFR]",
        re.DOTALL,
    )
    records = {
        match.group("og"): {
            "magnetic_type": match.group("type"),
            "corresponding_msg_bns_number": match.group("msg"),
        }
        for match in pattern.finditer(text)
    }
    if len(records) != 528:
        raise ValueError(f"expected 528 correspondence rows, found {len(records)}")
    transforms = {
        og: transform
        for og, _msg, transform in re.findall(
            r"(\d{1,2}\.\d{1,2}\.\d{1,3})\s+"
            r"(\d{1,3}\.\d{1,3})\s+(P[1-5])",
            text,
        )
    }
    if len(transforms) != 237:
        raise ValueError(f"expected 237 nonidentity transforms, found {len(transforms)}")
    for og, record in records.items():
        record["basis_transform"] = transforms.get(og, "I")
    return records


def _first_pages(text: str) -> dict[str, tuple[int, str]]:
    pages: dict[str, tuple[int, str]] = {}
    for page_number, page in enumerate(text.split("\f"), start=1):
        match = re.search(
            r"(?:N[0o]|No)\.\s*(\d+)\s*\.\s*(\d+)\s*\.\s*(\d+)", page
        )
        if match and "Symmetry operations" in page:
            pages[".".join(match.groups())] = (page_number, page)
    # The 2005 supplement prints 55.2.287 on its first page; the continuation
    # page and both independent number tables correctly identify 55.2.387.
    pages["55.2.387"] = pages.pop("55.2.287")
    if len(pages) != 528:
        raise ValueError(f"expected 528 first group pages, found {len(pages)}")
    return pages


def _first_column(line: str) -> str:
    return re.split(r"\s{3,}", line.strip())[0].replace(" ", "")


def _litvin_symbol(page: str) -> str:
    lines = [line for line in page.splitlines() if line.strip()]
    symbol = _first_column(lines[0])
    if symbol.startswith("&"):
        base = _first_column(lines[1])
        symbol = f"{base[0]}-{base[1:]}{symbol[1:]}"
    if not symbol.lower().startswith(("p", "c")) or "&" in symbol:
        raise ValueError(f"could not normalize Litvin symbol {symbol!r}")
    return symbol


def _litvin_point_group_symbol(page: str) -> str:
    """Reassemble the middle header column, including overbar glyphs."""

    lines = page.splitlines()
    first = next(line for line in lines if line.strip())
    fields = re.split(r"\s{3,}", first.strip())
    if len(fields) >= 3 and not fields[-2].startswith("&"):
        return fields[-2].replace(" ", "")
    if len(fields) == 2 and not fields[-1].startswith("&"):
        return fields[-1].replace(" ", "")
    number_line = next(
        line for line in lines if re.search(r"(?:N[0o]|No)\.", line)
    )
    number_match = re.search(
        r"(?:N[0o]|No)\.\s*\d+\s*\.\s*\d+\s*\.\s*\d+", number_line
    )
    if number_match is None:
        raise ValueError("missing magnetic-layer-group number")
    remainder = number_line[number_match.end() :]
    repeated_symbol = re.search(r"\S", remainder)
    if repeated_symbol is None:
        raise ValueError("missing repeated magnetic-layer-group symbol")
    center = number_match.end() + repeated_symbol.start()
    header = lines[: lines.index(number_line)]
    crystal_start = min(
        (
            match.start()
            for line in header
            if (
                match := re.search(
                    r"(?:Triclinic|Monoclinic|Orthorhombic|Tetragonal|Trigonal|Hexagonal)/",
                    line,
                )
            )
        ),
        default=max(len(line) for line in header),
    )
    fragments = [
        line[max(20, center - 20) : crystal_start].strip().replace(" ", "")
        for line in header
    ]
    fragments = [fragment for fragment in fragments if fragment]
    normal = next((item for item in fragments if not item.startswith("&")), "")
    overbar = next((item for item in fragments if item.startswith("&")), "")
    if overbar and normal:
        symbol = f"-{normal[0]}{overbar[1:]}{normal[1:]}"
    elif normal:
        symbol = normal
    else:
        raise ValueError(f"could not reconstruct magnetic point group from {fragments!r}")
    return symbol


def _operation_segments(page: str) -> dict[int, str]:
    section = page.split("Symmetry operations", maxsplit=1)[1]
    segments: dict[int, str] = {}
    for line in section.splitlines():
        labels = list(re.finditer(r"\(\s*(\d+)\)\s*", line))
        for position, label in enumerate(labels):
            number = int(label.group(1))
            end = labels[position + 1].start() if position + 1 < len(labels) else len(line)
            segments.setdefault(number, line[label.end() : end].strip())
    return segments


def _axis_suffix(token: str, segment: str, *, hexagonal: bool) -> str | None:
    segment = segment.split("&", maxsplit=1)[0]
    compact = re.sub(
        r"^(?:m|[abcegn]|[1-6](?:[+-])?)", "", segment.lstrip()
    ).replace(" ", "")
    compact = re.sub(r"\([^)]*\)", "", compact)
    mirror_like = token in {"m", "a", "b", "c", "e", "g", "n"}
    if "x,2x" in compact:
        return "100" if mirror_like and hexagonal else "120"
    if "2x,x" in compact:
        return "010" if mirror_like and hexagonal else "210"
    if compact.count("x") >= 2 and "y" not in compact:
        barred = ")" in compact
        if mirror_like:
            return "110" if barred else "1-10"
        return "1-10" if barred else "110"
    variables = {axis for axis in "xyz" if axis in compact}
    mirror_like = token in {"m", "a", "b", "c", "e", "g", "n"}
    if mirror_like:
        if hexagonal:
            return {
                frozenset(("x", "y")): "001",
                frozenset(("x", "z")): "120",
                frozenset(("y", "z")): "210",
            }.get(frozenset(variables))
        return {
            frozenset(("x", "y")): "001",
            frozenset(("x", "z")): "010",
            frozenset(("y", "z")): "100",
        }.get(frozenset(variables))
    return {
        frozenset(("z",)): "001",
        frozenset(("x",)): "100",
        frozenset(("y",)): "010",
    }.get(frozenset(variables))


def _matches_label(
    operation: dict[str, Any], number: int, segment: str, *, hexagonal: bool
) -> bool:
    match = re.match(r"(m|[abcegn]|[1-6](?:[+-])?)", segment.lstrip())
    if match is None:
        return False
    token = match.group(1)
    name = operation["name"]
    if number == 1:
        return name == "1"
    if token in {"m", "a", "b", "c", "e", "g", "n"}:
        if not name.startswith("m"):
            return False
    else:
        base = name.removeprefix("-").split("_", maxsplit=1)[0]
        if not base.startswith(token):
            return False
    suffix = _axis_suffix(token, segment, hexagonal=hexagonal)
    return suffix is None or "_" not in name or name.rsplit("_", maxsplit=1)[1] == suffix


def _table_order(
    operations: list[dict[str, Any]], page: str
) -> list[dict[str, Any]]:
    segments = _operation_segments(page)
    if set(segments) != set(range(1, len(operations) + 1)):
        raise ValueError("type-I operation numbers are incomplete")
    hexagonal = any(item["name"].endswith(("_120", "_210")) for item in operations)
    choices = {
        number: (
            matches
            if (
                matches := [
                    index
                    for index, operation in enumerate(operations)
                    if _matches_label(
                        operation, number, segments[number], hexagonal=hexagonal
                    )
                ]
            )
            else list(range(len(operations)))
        )
        for number in segments
    }
    solutions: list[list[int]] = []

    def visit(number: int, used: set[int], assignment: list[int]) -> None:
        if number > len(operations):
            solutions.append(assignment.copy())
            return
        for index in choices[number]:
            if index not in used:
                used.add(index)
                assignment.append(index)
                visit(number + 1, used, assignment)
                assignment.pop()
                used.remove(index)

    visit(1, set(), [])
    if not solutions:
        raise ValueError(f"could not match type-I operation labels: {choices!r}")
    # Residual ambiguity is limited to inverse +/- pairs.  A Z2 character has
    # the same value on an element and its inverse, so the first deterministic
    # catalog ordering is sufficient for extracting time-reversal labels.
    return [operations[index] for index in solutions[0]]


def _parent_operations(
    first_pages: dict[str, tuple[int, str]],
    correspondence: dict[str, dict[str, str]],
) -> dict[int, list[dict[str, Any]]]:
    database = json.loads(OPERATION_PATH.read_text(encoding="utf-8"))
    result: dict[int, list[dict[str, Any]]] = {}
    for family in database["families"].values():
        by_index = {item["index"]: item for item in family["operations"]}
        for group in family["layer_groups"]:
            indices = sorted(set(group["R+_indices"] + group["R-_indices"]))
            result[group["LG"]] = [by_index[index] for index in indices]
    if set(result) != set(range(1, 81)):
        raise ValueError("point-operation catalog must cover LG1-LG80")
    type_i = {
        int(og.split(".")[0]): og
        for og, source in correspondence.items()
        if source["magnetic_type"] == "I"
    }
    return {
        number: _table_order(operations, first_pages[type_i[number]][1])
        for number, operations in result.items()
    }


def _homomorphisms(operations: list[dict[str, Any]]) -> list[tuple[bool, ...]]:
    matrices = [np.asarray(item["matrix_fractional"], dtype=np.int64) for item in operations]
    lookup = {tuple(matrix.ravel()): index for index, matrix in enumerate(matrices)}
    equation_masks: set[int] = set()
    for left_index, left in enumerate(matrices):
        for right_index, right in enumerate(matrices):
            product_index = lookup[tuple((left @ right).ravel())]
            mask = (1 << left_index) ^ (1 << right_index) ^ (1 << product_index)
            if mask:
                equation_masks.add(mask)

    rows = list(equation_masks)
    pivots: list[int] = []
    for column in range(len(operations)):
        pivot = next(
            (
                index
                for index in range(len(pivots), len(rows))
                if (rows[index] >> column) & 1
            ),
            None,
        )
        if pivot is None:
            continue
        rows[len(pivots)], rows[pivot] = rows[pivot], rows[len(pivots)]
        for index in range(len(rows)):
            if index != len(pivots) and (rows[index] >> column) & 1:
                rows[index] ^= rows[len(pivots)]
        pivots.append(column)

    free = [index for index in range(len(operations)) if index not in pivots]
    results: list[tuple[bool, ...]] = []
    for choice in range(1 << len(free)):
        vector = sum(((choice >> bit) & 1) << column for bit, column in enumerate(free))
        for row_index in range(len(pivots) - 1, -1, -1):
            pivot = pivots[row_index]
            remainder = rows[row_index] & ~(1 << pivot)
            if (remainder & vector).bit_count() % 2:
                vector |= 1 << pivot
        results.append(
            tuple(bool((vector >> index) & 1) for index in range(len(operations)))
        )
    return results


def _visible_prime_flags(page: str) -> dict[int, bool]:
    section = page.split("Symmetry operations", maxsplit=1)[1]
    flags: dict[int, bool] = {}
    for line in section.splitlines():
        labels = list(re.finditer(r"\(\s*(\d+)\)\s*", line))
        for position, label in enumerate(labels):
            number = int(label.group(1))
            end = labels[position + 1].start() if position + 1 < len(labels) else len(line)
            segment = line[label.end() : end].lstrip()
            token = re.match(
                r"(?:m|[abcegn]|[1-6](?:[+-])?)\s*(')?(?:\s|$)", segment
            )
            if token:
                flags[number] = bool(token.group(1))
    # Identity is unitary in every type-III group.  This also shields the
    # parser from overbar glyphs that the PDF text layer moves across columns.
    flags[1] = False
    return flags


def _type_iii_primes(
    page: str,
    operations: list[dict[str, Any]],
) -> tuple[int, ...]:
    visible = _visible_prime_flags(page)
    known = {
        index: visible[index + 1]
        for index, operation in enumerate(operations)
        if index + 1 in visible and not operation["name"].startswith("-")
    }
    candidates = [
        assignment
        for assignment in _homomorphisms(operations)
        if any(assignment)
        and all(assignment[index] == value for index, value in known.items())
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"type-III {_litvin_point_group_symbol(page)} resolved to "
            f"{len(candidates)} homomorphisms"
        )
    return tuple(index + 1 for index, value in enumerate(candidates[0]) if value)


def _anti_translation(page: str) -> str:
    section = page.split("Symmetry operations", maxsplit=1)[1]
    sets = list(re.finditer(r"For\s+\(([^\n]*?)\)'?\s*\+\s*set", section))
    matches: list[str] = []
    for index, item in enumerate(sets):
        end = sets[index + 1].start() if index + 1 < len(sets) else len(section)
        if "t'" in section[item.end() : end]:
            matches.append("".join(item.group(1).split()))
    if len(matches) != 1:
        raise ValueError(f"expected one anti-translation, found {matches!r}")
    return matches[0]


def extract(litvin_pdf: Path, correspondence_pdf: Path) -> list[dict[str, Any]]:
    correspondence = _correspondence_rows(
        _pdf_text(correspondence_pdf, first=2, last=10)
    )
    first_pages = _first_pages(_pdf_text(litvin_pdf))
    parent_operations = _parent_operations(first_pages, correspondence)
    bns_to_uni = {
        spglib.get_magnetic_spacegroup_type(number).bns_number: number
        for number in range(1, 1652)
    }
    records: list[dict[str, Any]] = []
    for og_number, source in correspondence.items():
        parent_number, family_number, global_number = map(int, og_number.split("."))
        page_number, page = first_pages[og_number]
        magnetic_type = source["magnetic_type"]
        operations = parent_operations[parent_number]
        primed = (
            _type_iii_primes(page, operations)
            if magnetic_type == "III"
            else ()
        )
        anti_translation = _anti_translation(page) if magnetic_type == "IV" else ""
        bns_number = source["corresponding_msg_bns_number"]
        uni_number = bns_to_uni[bns_number]
        msg_type = spglib.get_magnetic_spacegroup_type(uni_number)
        if msg_type.type != TYPE_TO_SPGLIB[magnetic_type]:
            raise ValueError(f"{og_number} and MSG {bns_number} have different types")
        records.append(
            {
                "og_number": og_number,
                "global_number": global_number,
                "parent_layer_group_number": parent_number,
                "family_number": family_number,
                "magnetic_type": magnetic_type,
                "litvin_og_symbol_ascii": _litvin_symbol(page),
                "magnetic_point_group_symbol": _litvin_point_group_symbol(page),
                "source_pdf_page": page_number,
                "point_operation_names": ",".join(
                    operation["name"] for operation in operations
                ),
                "primed_operation_indices": ",".join(map(str, primed)),
                "anti_translation": anti_translation,
                "corresponding_msg_bns_number": bns_number,
                "corresponding_msg_uni_number": uni_number,
                "corresponding_msg_og_number": msg_type.og_number,
                "basis_transform": source["basis_transform"],
            }
        )
    records.sort(key=lambda item: item["global_number"])
    if [item["global_number"] for item in records] != list(range(1, 529)):
        raise ValueError("global magnetic-layer-group numbering must be 1-528")
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("litvin_pdf", type=Path)
    parser.add_argument("correspondence_pdf", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    records = extract(args.litvin_pdf, args.correspondence_pdf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"wrote {len(records)} rows to {args.output}")
    print(f"litvin_sha256={hashlib.sha256(args.litvin_pdf.read_bytes()).hexdigest()}")
    print(
        "correspondence_sha256="
        f"{hashlib.sha256(args.correspondence_pdf.read_bytes()).hexdigest()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
