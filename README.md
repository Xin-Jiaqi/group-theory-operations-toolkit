# Group Theory Operations Toolkit

[![Verify repository](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml/badge.svg)](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml)

这是我为二维材料、层间堆叠与铁电理论研究整理的群论数据仓库。我希望把后续推导中反复使用的点操作、矩阵表示、层群分类和操作乘法集中在一个可靠的数据源中，便于人工查阅，也便于程序和 AI 直接调用。

## 内容

仓库目前包含六个相互校验的机器数据集。下面三项是基础点操作族；另外三项分别是 32 个点群与光学响应、230 个空间群、80 个层群的完整注册表。

| 数据集 | 母点群 | 基础操作 | 群目录 | 乘法表 |
|---|---|---:|---|---:|
| `cubic_Oh` | $O_h$ | 48 | $O_h$、$O$、$T_d$、$T_h$、$T$ | — |
| `tetragonal_D4h` | $D_{4h}$ | 16 | LG1–LG64 | 16×16 |
| `hexagonal_D6h` | $D_{6h}$ | 24 | LG65–LG80 | 24×24 |

`data/group_operations.json` 是点操作的唯一事实源；`data/quadratic_field_representations.json` 是由它确定性生成的 $M_+$/$M_-$ 全操作机器表。`docs/group_theory.md` 汇总点操作、层群集合和乘法表，[`docs/quadratic_field_representations.md`](docs/quadratic_field_representations.md) 给出二次场推导及 88 个 $M_-(R)$。Python/CLI 负责查询，测试负责矩阵、群闭合、叉积恒等式、诱导表示同态和跨仓库结构契约。

新增的 [`data/crystallographic_point_groups.json`](data/crystallographic_point_groups.json) 固定 32 个晶体学点群的标准序号、HM/Schoenflies 名称、设置、生成元和闭包操作；[`data/optical_response_invariants.json`](data/optical_response_invariants.json) 给出 shift current、SHG 与 circular injection current 的全部空间点群不变量基。两者均绑定源数据 SHA-256，并有独立 JSON Schema。

新增的 [`data/crystallographic_space_groups.json`](data/crystallographic_space_groups.json) 覆盖全部 230 个空间群：每个 ITA 序号的国际符号、Schoenflies、母点群、晶系、心型、symmorphic 判定，以及全部 530 个 Hall 设置各自的 Seitz 生成元 $(R\mid\mathbf t)$ 与操作数。数据由 spglib 数据库（BSD-3-Clause）生成，并经过 spglib、ASE 双源交叉验证与端到端回环验证；详见 [`docs/space_groups.md`](docs/space_groups.md)。

[`data/crystallographic_layer_groups.json`](data/crystallographic_layer_groups.json) 覆盖 LG1–LG80 及全部 116 个 layer Hall 设置，保存标准/备选符号、母点群、晶系、心型、Seitz 生成元和操作数，并与原有层群点操作分类交叉校验。它由仓库内固定的 spglib v2.5.0 BSD-3-Clause 源表确定性生成；读取已生成数据只依赖 NumPy，不依赖 spglib。

## 使用

核心查询工具需要 Python 3.10 或更高版本；NumPy 是唯一的核心运行时依赖。

```bash
# 查看数据集
python -m pip install -e .
group-ops list

# 查询矩阵
group-ops show '4+_001' --family tetragonal_D4h

# 查询旧版层群点操作分类
group-ops group --family hexagonal_D6h --lg 80

# 查询完整 80 层群注册表（符号、设置与 Seitz 生成元）
group-ops layer-groups 80
group-ops layer-groups p4/nmm --json

# 查询点群
group-ops group --family cubic_Oh --point-group Td

# 查询乘积
group-ops multiply '4+_001' 2_100 --family tetragonal_D4h

# 查询对称/反对称二次场表示
group-ops field-representation m_100 \
  --family tetragonal_D4h --space antisymmetric --json

# 查询 32 点群注册表与允许张量基
group-ops point-groups 4mm --json
group-ops invariants 4mm shift_current --json

# 查询 230 空间群（含 Seitz 生成元）
group-ops space-groups 227
group-ops space-groups 227 --json

# 对结构施加点操作；文件 I/O 委托给 materials-structure-core
group-ops apply-structure structure.vasp 3+_001 \
  --family hexagonal_D6h --output transformed.vasp

# 供脚本稳定读取
group-ops show 4+_001 --family tetragonal_D4h --json
group-ops validate
```

## 约定

矩阵采用列向量约定 $\mathbf r'=D(R)\mathbf r$。操作名统一写成 `2_001`、`4+_001`、`m_100`；程序也接受 `2001`、`4+001`、`m100` 和 `4^+_{001}` 等旧写法。六角晶系同时保存晶格分数坐标矩阵、精确根式笛卡尔矩阵和数值笛卡尔矩阵。对称二次场默认基底为 $(|E_x|^2,|E_y|^2,|E_z|^2,E_xE_y^*+E_yE_x^*,E_xE_z^*+E_zE_x^*,E_yE_z^*+E_zE_y^*)$；反对称基底为 $\mathbf h=i\mathbf E\times\mathbf E^*$，因此 $M_-(R)=\det[D(R)]D(R)$。

不变量求解器使用 $D(R)T_+=T_+M_+(R)$ 与 $D(R)T_-=T_-M_-(R)$。shift current/SHG 的基矩阵为 3×6；circular injection current 为 3×3。它们只表示空间点群选择定则，不替代时间反演、频率、单位、共振或微观机制分析。

## 32 点群与响应维数

注册表使用工具包固定的正交 Cartesian 嵌入，主旋转轴沿 `z`；显式操作列表是最终约定。下表给出全部空间点群解的自由度，完整生成元、闭包操作和基矩阵分别位于两个版本化 JSON 中。

| No. | HM | Schoenflies | Order | Centrosymmetric | Shift | SHG | Circular injection |
|---:|---|---|---:|:---:|---:|---:|---:|
| 1 | `1` | `C1` | 1 | no | 18 | 18 | 9 |
| 2 | `-1` | `Ci` | 2 | yes | 0 | 0 | 0 |
| 3 | `2` | `C2` | 2 | no | 8 | 8 | 5 |
| 4 | `m` | `Cs` | 2 | no | 10 | 10 | 4 |
| 5 | `2/m` | `C2h` | 4 | yes | 0 | 0 | 0 |
| 6 | `222` | `D2` | 4 | no | 3 | 3 | 3 |
| 7 | `mm2` | `C2v` | 4 | no | 5 | 5 | 2 |
| 8 | `mmm` | `D2h` | 8 | yes | 0 | 0 | 0 |
| 9 | `4` | `C4` | 4 | no | 4 | 4 | 3 |
| 10 | `-4` | `S4` | 4 | no | 4 | 4 | 2 |
| 11 | `4/m` | `C4h` | 8 | yes | 0 | 0 | 0 |
| 12 | `422` | `D4` | 8 | no | 1 | 1 | 2 |
| 13 | `4mm` | `C4v` | 8 | no | 3 | 3 | 1 |
| 14 | `-42m` | `D2d` | 8 | no | 2 | 2 | 1 |
| 15 | `4/mmm` | `D4h` | 16 | yes | 0 | 0 | 0 |
| 16 | `3` | `C3` | 3 | no | 6 | 6 | 3 |
| 17 | `-3` | `C3i` | 6 | yes | 0 | 0 | 0 |
| 18 | `32` | `D3` | 6 | no | 2 | 2 | 2 |
| 19 | `3m` | `C3v` | 6 | no | 4 | 4 | 1 |
| 20 | `-3m` | `D3d` | 12 | yes | 0 | 0 | 0 |
| 21 | `6` | `C6` | 6 | no | 4 | 4 | 3 |
| 22 | `-6` | `C3h` | 6 | no | 2 | 2 | 0 |
| 23 | `6/m` | `C6h` | 12 | yes | 0 | 0 | 0 |
| 24 | `622` | `D6` | 12 | no | 1 | 1 | 2 |
| 25 | `6mm` | `C6v` | 12 | no | 3 | 3 | 1 |
| 26 | `-6m2` | `D3h` | 12 | no | 1 | 1 | 0 |
| 27 | `6/mmm` | `D6h` | 24 | yes | 0 | 0 | 0 |
| 28 | `23` | `T` | 12 | no | 1 | 1 | 1 |
| 29 | `m-3` | `Th` | 24 | yes | 0 | 0 | 0 |
| 30 | `432` | `O` | 24 | no | 0 | 0 | 1 |
| 31 | `-43m` | `Td` | 24 | no | 1 | 1 | 0 |
| 32 | `m-3m` | `Oh` | 48 | yes | 0 | 0 | 0 |

`basis[k]` 的行按 `(x,y,z)`；shift/SHG 列按 `(xx,yy,zz,xy,xz,yz)`；circular injection 列按 `(h_x,h_y,h_z)`。任意线性组合都是允许张量。重新生成时先运行 `scripts/generate_crystallographic_point_groups.py`，再运行 `scripts/generate_optical_response_invariants.py`。

乘法表使用“行操作 × 列操作”：`multiplication.table[left][right]` 对应 $D(left)D(right)$，对列向量而言右操作先作用。JSON 使用 `1` 和 `-1` 表示恒等与反演，Markdown 表中显示为 $E$ 和 $I$。Markdown 乘法表按 8 个列操作分块，以保证 Obsidian 和 GitHub 都能正常渲染；程序应直接读取 JSON，而不是解析 Markdown 表格。

## 数据结构

每个操作都有稳定的 `index`、规范化 `name`、矩阵和坐标作用。层群条目同时保存操作索引与名称，并用 `point_group_base`、`point_group_embedding` 明确记录 `100/110/120`、IP/OP 等不同嵌入。$D_{4h}$ 和 $D_{6h}$ 的 `multiplication` 字段包含元素顺序、单位元、逆元以及显式的 `table[left][right] -> result` 映射。派生表逐操作保存 `determinant`、`matrix_symmetric` 与 `matrix_antisymmetric`，并绑定源数据 SHA-256。

## 核对与测试

我对 $D_{4h}$、$D_{6h}$ 的操作顺序、坐标作用、矩阵、层群映射和点群嵌入进行了逐项一致性检查。两张乘法表的 832 个乘积不是人工抄写，而是由精确整数矩阵生成并反向验证。

```bash
python3 -m unittest discover -s tests -v
```

当前测试覆盖矩阵正交性、非正交分数坐标中的 Seitz 逆运算、群闭合性、坐标作用、LG1–LG80 与全部 116 个 layer Hall 设置、全部乘法单元、JSON/Markdown 同步、$M_+$/$M_-$ 定义与同态、32 点群生成闭包、230 空间群、三类响应的零空间维数与逐操作等变性、机器接口、PBC 与 Selective-dynamics 拒绝边界。

参见[机器接口与跨仓库契约](docs/MACHINE_INTERFACE.md)和[发布路线图](ROADMAP.md)。POSCAR/CIF 语法不再由本仓库手写解析，而由 `materials-structure-core` 的维护型后端统一负责。

科研使用时请通过 [`CITATION.cff`](CITATION.cff) 引用所使用的准确版本；只有在对应发布或论文确实存在后，才会填写 DOI 或首选论文引用。

## 范围与后续

目前点操作数据描述操作的线性部分；结构变换默认围绕原点。空间群的 Seitz 对 $(R\mid\mathbf t)$（含非零平移）已由 [`crystallographic_space_groups.json`](data/crystallographic_space_groups.json) 覆盖，但逐结构的平移/滑移操作应用于具体 POSCAR 的 affine Seitz 契约仍待固定。张量工具给出允许分量空间，不计算材料响应数值。

当前版本为 `0.4.1`；代码、JSON 数据和文档均采用 [BSD 3-Clause License](LICENSE)。在接受含平移的操作直接变换具体结构前，仍需固定原点选择、Wyckoff 位置和结构级 affine Seitz 契约。
