#!/usr/bin/env python3
"""Generate machine and human views of quadratic-field representations."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from group_theory_operations import (  # noqa: E402
    load_database,
    quadratic_field_representation,
)


DATA_PATH = ROOT / "data" / "quadratic_field_representations.json"
DOC_PATH = ROOT / "docs" / "quadratic_field_representations.md"


def _number(value: float) -> int | float:
    if abs(value) < 1e-12:
        return 0
    nearest = round(value)
    if abs(value - nearest) < 1e-12:
        return int(nearest)
    return round(value, 12)


def _matrix(matrix) -> list[list[int | float]]:
    return [[_number(value) for value in row] for row in matrix]


def build_payload(database: dict) -> dict:
    source_bytes = (ROOT / "data" / "group_operations.json").read_bytes()
    families = {}
    for family_name, family in database["families"].items():
        operations = []
        for operation in family["operations"]:
            representation = quadratic_field_representation(operation)
            operations.append(
                {
                    "index": operation["index"],
                    "name": operation["name"],
                    "determinant": _number(representation.determinant),
                    "matrix_symmetric": _matrix(representation.matrix_symmetric),
                    "matrix_antisymmetric": _matrix(representation.matrix_antisymmetric),
                }
            )
        families[family_name] = {
            "parent_point_group": family["parent_point_group"],
            "operation_count": len(operations),
            "operations": operations,
        }
    return {
        "schema_version": 1,
        "source_catalog_schema_version": database["schema_version"],
        "source_catalog_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "generated_by": "scripts/generate_quadratic_field_representations.py",
        "conventions": {
            "polar_action": "r' = D(R) r; Cartesian orthonormal basis; column vectors",
            "symmetric_basis": [
                "|E_x|^2",
                "|E_y|^2",
                "|E_z|^2",
                "E_x E_y^* + E_y E_x^*",
                "E_x E_z^* + E_z E_x^*",
                "E_y E_z^* + E_z E_y^*",
            ],
            "symmetric_action": "M_+(R) = Sym^2 D(R)",
            "antisymmetric_basis": [
                "i(E_y E_z^* - E_z E_y^*)",
                "i(E_z E_x^* - E_x E_z^*)",
                "i(E_x E_y^* - E_y E_x^*)",
            ],
            "antisymmetric_action": "M_-(R) = det(D(R)) D(R)",
        },
        "families": families,
    }


def _tex_number(value: float) -> str:
    candidates = (
        (0.0, "0"),
        (1.0, "1"),
        (-1.0, "-1"),
        (0.5, r"\frac12"),
        (-0.5, r"-\frac12"),
        (math.sqrt(3) / 2, r"\frac{\sqrt3}{2}"),
        (-math.sqrt(3) / 2, r"-\frac{\sqrt3}{2}"),
    )
    for candidate, text in candidates:
        if abs(value - candidate) < 1e-9:
            return text
    return f"{value:.9g}"


def _tex_matrix(matrix: list[list[int | float]]) -> str:
    rows = ["&".join(_tex_number(float(value)) for value in row) for row in matrix]
    return r"$\begin{bmatrix}" + r"\\".join(rows) + r"\end{bmatrix}$"


def build_markdown(payload: dict) -> str:
    lines = [
        "# 二次场空间的诱导表示",
        "",
        "本文档由 `data/group_operations.json` 确定性生成。机器读取请使用 "
        "`data/quadratic_field_representations.json`；不要解析下列表格。",
        "",
        "## 约定与推导",
        "",
        r"令 $Q^{bc}=E^bE^{c*}$。在对称基底",
        "",
        r"$$\boldsymbol{\mathcal E}_+=(|E^x|^2,|E^y|^2,|E^z|^2,E^xE^{y*}+E^yE^{x*},E^xE^{z*}+E^zE^{x*},E^yE^{z*}+E^zE^{y*})^{\mathsf T}$$",
        "",
        r"上，$\boldsymbol{\mathcal E}'_+=M_+(R)\boldsymbol{\mathcal E}_+$，且 $M_+(R)=\operatorname{Sym}^2D(R)$。",
        "",
        r"反对称子空间与实轴矢量 $\mathbf h=i\mathbf E\times\mathbf E^*$ 等价。由",
        "",
        r"$$[D\mathbf u]\times[D\mathbf v]=\det(D)D(\mathbf u\times\mathbf v)$$",
        "",
        r"得到 $\mathbf h'=M_-(R)\mathbf h$，其中",
        "",
        r"$$\boxed{M_-(R)=\det[D(R)]D(R)}.$$",
        "",
        r"因此 $M_-(R_1R_2)=M_-(R_1)M_-(R_2)$、$M_-^{-1}=M_-^{\mathsf T}$。反演满足 $D(I)=-I_3$，但 $M_-(I)=I_3$。所有 $M_-(R)$ 均为行列式 $+1$ 的正交矩阵。",
        "",
        "## 全部操作的 $M_-(R)$",
        "",
    ]
    for family_name, family in payload["families"].items():
        lines.extend(
            [
                f"### `{family_name}`：{family['operation_count']} 个操作",
                "",
                r"| Index | $R$ | $\det D(R)$ | $M_-(R)$ |",
                "|---:|---|---:|---|",
            ]
        )
        for operation in family["operations"]:
            lines.append(
                f"| {operation['index']} | `{operation['name']}` | "
                f"{operation['determinant']:+d} | "
                f"{_tex_matrix(operation['matrix_antisymmetric'])} |"
            )
        lines.append("")
    lines.extend(
        [
            "## 机器验证",
            "",
            "测试逐操作验证定义式、正交性、行列式、叉积恒等式、$M_+$/$M_-$ 同态，"
            "并核对本文件与生成 JSON。88 个母群操作均覆盖；同名操作在不同母群中保留各自记录。",
            "",
            "重新生成：`python scripts/generate_quadratic_field_representations.py`。",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    payload = build_payload(load_database())
    DATA_PATH.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    DOC_PATH.write_text(build_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
