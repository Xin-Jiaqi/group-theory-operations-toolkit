# Group Theory Operations Toolkit

这是我为二维材料、层间堆叠与铁电理论研究整理的群论数据仓库。我希望把后续推导中反复使用的点操作、矩阵表示、层群分类和操作乘法集中在一个可靠的数据源中，便于人工查阅，也便于程序和 AI 直接调用。

## 内容

仓库目前包含三个数据集：

| 数据集 | 母点群 | 基础操作 | 群目录 | 乘法表 |
|---|---:|---:|---|---:|
| `cubic_Oh` | $O_h$ | 48 | $O_h$、$O$、$T_d$、$T_h$、$T$ | — |
| `tetragonal_D4h` | $D_{4h}$ | 16 | LG1–LG64 | 16×16 |
| `hexagonal_D6h` | $D_{6h}$ | 24 | LG65–LG80 | 24×24 |

`data/group_operations.json` 是唯一机器可读数据源；`docs/group_theory.md` 汇总公式、全部矩阵、层群操作集合和乘法表；`group_theory_operations` 提供可安装的 Python/CLI 接口；`group_tools.py` 仅保留为兼容入口；测试负责数据、群闭合和跨仓库结构契约校验。

## 使用

核心查询工具只需要 Python 3.10 或更高版本，没有第三方运行时依赖。

```bash
# 查看数据集
python -m pip install -e .
group-ops list

# 查询矩阵
group-ops show '4+_001' --family tetragonal_D4h

# 查询层群
group-ops group --family hexagonal_D6h --lg 80

# 查询点群
group-ops group --family cubic_Oh --point-group Td

# 查询乘积
group-ops multiply '4+_001' 2_100 --family tetragonal_D4h

# 对结构施加点操作；文件 I/O 委托给 materials-structure-core
group-ops apply-structure structure.vasp 3+_001 \
  --family hexagonal_D6h --output transformed.vasp

# 供脚本稳定读取
group-ops show 4+_001 --family tetragonal_D4h --json
group-ops validate
```

## 约定

矩阵采用列向量约定 $\mathbf r'=D(R)\mathbf r$。操作名统一写成 `2_001`、`4+_001`、`m_100`；程序也接受 `2001`、`4+001`、`m100` 和 `4^+_{001}` 等旧写法。六角晶系同时保存晶格分数坐标矩阵、精确根式笛卡尔矩阵和数值笛卡尔矩阵。

乘法表使用“行操作 × 列操作”：`multiplication.table[left][right]` 对应 $D(left)D(right)$，对列向量而言右操作先作用。JSON 使用 `1` 和 `-1` 表示恒等与反演，Markdown 表中显示为 $E$ 和 $I$。Markdown 乘法表按 8 个列操作分块，以保证 Obsidian 和 GitHub 都能正常渲染；程序应直接读取 JSON，而不是解析 Markdown 表格。

## 数据结构

每个操作都有稳定的 `index`、规范化 `name`、矩阵和坐标作用。层群条目同时保存操作索引与名称，并用 `point_group_base`、`point_group_embedding` 明确记录 `100/110/120`、IP/OP 等不同嵌入。$D_{4h}$ 和 $D_{6h}$ 的 `multiplication` 字段包含元素顺序、单位元、逆元以及显式的 `table[left][right] -> result` 映射。

## 核对与测试

我逐项对照了 Bilbao Crystallographic Server 的 $D_{4h}(4/mmm)$、$D_{6h}(6/mmm)$ General Positions 截图，以及 Layer Group PDF 中的 LG1–LG80、$R^+/R^-$ 操作集合和点群嵌入。两张乘法表的元素顺序与显示名称也已核对；832 个乘积不是人工抄写，而是由精确整数矩阵生成。

```bash
python3 -m unittest discover -s tests -v
```

当前测试覆盖矩阵正交性、六角基底变换、群闭合性、Bilbao 坐标作用、LG1–LG80 分类、全部乘法单元、单位元、逆元、结合律、JSON/Markdown 同步、机器接口、乘法/结构变换一致性、PBC 与 Selective-dynamics 拒绝边界。

参见[机器接口与跨仓库契约](docs/MACHINE_INTERFACE.md)、[0.1 迁移说明](docs/MIGRATION_0_1.md)、[来源与验证边界](docs/SOURCES.md)和[发布路线图](ROADMAP.md)。POSCAR/CIF 语法不再由本仓库手写解析，而由 `materials-structure-core` 的维护型后端统一负责。

科研使用时请通过 [`CITATION.cff`](CITATION.cff) 引用所使用的准确版本；只有在对应发布或论文确实存在后，才会填写 DOI 或首选论文引用。

## 范围与后续

目前数据描述点操作的线性部分，不含非零平移的完整 Seitz 对 $(R\mid\mathbf t)$；结构变换默认围绕原点。后续扩展会优先明确平移、旋转中心、容差和数据来源，再加入表示分解与不变量构造，避免把尚未验证的应用目标写成现有能力。

参考资料的核对范围记录在 JSON 顶层 `sources` 字段中。`0.1.0` 仍是发布候选；正式稳定发布前还需要确定代码、数据和文档许可证，并补充每个 Bilbao 页面、访问日期及可复现快照策略。
