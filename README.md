# Group Theory Operations Toolkit

[![Verify repository](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml/badge.svg)](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml)

这是我为二维材料、层间堆叠、铁电与非线性光学研究维护的群论数据工具。它把常用对称操作、矩阵、群注册表和允许张量基放在同一套经过验证的机器接口中，既适合直接查阅，也适合作为高通量程序的底层依赖。

## 数据与能力

| 数据 | 覆盖范围 | 主要内容 |
|---|---:|---|
| 基础点操作 | 88 | $O_h$、$D_{4h}$、$D_{6h}$ 的分数/笛卡尔矩阵与乘法表 |
| 二次场表示 | 88 | $M_+(R)=\operatorname{Sym}^2D(R)$ 与 $M_-(R)=\det[D(R)]D(R)$ |
| 晶体学点群 | 32 | 标准编号、HM/Schoenflies、生成元、闭包操作 |
| 磁点群 | 122 | type I、灰群、黑白群及时间反演标签 |
| 空间群 | 230 | 530 个 Hall 设置、Seitz 生成元与操作数 |
| 层群 | 80 | 116 个 layer Hall 设置、符号、母点群与 Seitz 生成元 |
| 磁性层群 | 528 | OG 编号、I–IV 型、点余群操作、反平移和对应磁空间群 |
| 响应求解 | 按需 | shift current、injection current、SHG 与通用时间奇偶张量 |

主要机器数据位于 [`data/`](data)，对应 JSON Schema 位于 [`schema/`](schema)。[`data/group_operations.json`](data/group_operations.json) 是基础点操作的事实源；常用入口包括 [`data/crystallographic_point_groups.json`](data/crystallographic_point_groups.json)、[`data/optical_response_invariants.json`](data/optical_response_invariants.json) 和 [`data/magnetic_layer_groups.json`](data/magnetic_layer_groups.json)。其他表由版本化脚本确定性生成并绑定来源文件或上游数据的 SHA-256。

磁性层群注册表面向均匀张量选择定则，保存有限磁性点余群 $(R,\theta)$。普通空间平移不进入 $q=0$ 均匀张量约束；type IV 的特征反平移仍单独保存在 `anti_translation_fractional`。这一区分避免把点余群接口误当成完整仿射磁性层群作用。

## 安装与查询

需要 Python 3.10 或更高版本，核心运行依赖只有 NumPy。

```bash
python -m pip install -e .
group-ops validate

# 点操作、乘法与 M_-/M_+
group-ops show '4+_001' --family tetragonal_D4h
group-ops multiply '4+_001' 2_100 --family tetragonal_D4h
group-ops field-representation m_100 --family tetragonal_D4h --space antisymmetric --json

# 晶体学与磁性群注册表
group-ops point-groups 4mm --json
group-ops space-groups 227 --json
group-ops layer-groups p4/nmm --json
group-ops magnetic-point-groups "6'/m'mm'" --json
group-ops magnetic-layer-groups 80.9.528 --json
group-ops magnetic-layer-groups --type IV

# 非线性光学允许张量基
group-ops invariants 4mm shift_current --json
group-ops magnetic-responses 74 shg_odd --json
group-ops magnetic-layer-responses 6.5.25 shg_odd --json
group-ops magnetic-layer-invariants 6.5.25 polar_vector symmetric_quadratic \
  --output-time even --input-time odd --json
```

Python 接口返回带类型的不可变记录；示例：

```python
from group_theory_operations import (
    get_magnetic_layer_group,
    magnetic_layer_response_tensor_basis,
)

group = get_magnetic_layer_group("6.5.25")
basis = magnetic_layer_response_tensor_basis(group.og_number, "shg_odd")
print(group.point_operations)
print(basis.dimension, basis.basis)
```

## 约定

矩阵采用列向量约定 $\mathbf r'=D(R)\mathbf r$。乘法表使用“行操作 × 列操作”，即 $D(\text{left})D(\text{right})$，右操作先作用。对称二次场基底为 `(xx, yy, zz, xy, xz, yz)`；反对称场以 $\mathbf h=i\mathbf E\times\mathbf E^*$ 表示，基底为 `(h_x, h_y, h_z)`。

空间点群求解器使用 $D(R)T_+=T_+M_+(R)$ 与 $D(R)T_-=T_-M_-(R)$。磁性接口进一步写入 `time_reversal`，并区分 normal/magnetic shift current、normal/magnetic injection current 及 i-type/c-type SHG。输出是对称性允许空间，不包含材料数值、单位、共振、耗散或具体频率约定。

对具体 POSCAR/CIF 施加操作时，文件读写由可选的 `materials-structure-core[io]` 提供：

```bash
group-ops apply-structure structure.vasp 3+_001 \
  --family hexagonal_D6h --output transformed.vasp
```

## 核验与来源

```bash
python3 -m unittest discover -s tests -v
```

测试覆盖矩阵、逆运算、乘法表、群闭包、$M_\pm$ 表示同态、32 个点群、122 个磁点群、230 个空间群、80 个层群、528 个磁性层群及六类磁性光学响应。空间/层群数据与固定版本的 spglib 数据交叉验证；磁性层群的编号、类型、操作标记和磁空间群对应关系分别依据 D. B. Litvin 的 *Magnetic Group Tables, Part 3* 与 [Zhang et al., Phys. Rev. B **107**, 075405 (2023)](https://doi.org/10.1103/PhysRevB.107.075405)，并再次对照 spglib 的磁空间群数据库。外部原始出版物不收入仓库，只保存可复现派生表、正式链接与校验值。

磁性光学扇区约定对照 [Gao et al., npj Computational Materials **7**, 10 (2021)](https://doi.org/10.1038/s41524-020-00462-9) 和 [Xiao et al., npj Quantum Materials **8**, 62 (2023)](https://doi.org/10.1038/s41535-023-00594-3)。更细的接口说明见 [`docs/MACHINE_INTERFACE.md`](docs/MACHINE_INTERFACE.md)。

当前版本为 `0.7.0`，采用 [BSD 3-Clause License](LICENSE)。科研使用请通过 [`CITATION.cff`](CITATION.cff) 引用实际使用的版本。
