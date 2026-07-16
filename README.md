# Group Theory Operations Toolkit

这是我为二维材料、层间堆叠与铁电理论研究整理的群论数据仓库。我希望把后续推导中反复使用的点操作、矩阵表示、层群分类和操作乘法集中在一个可靠的数据源中，便于人工查阅，也便于程序和 AI 直接调用。

## 内容

仓库目前包含三个数据集：

| 数据集 | 母点群 | 基础操作 | 群目录 | 乘法表 |
|---|---:|---:|---|---:|
| `cubic_Oh` | $O_h$ | 48 | $O_h$、$O$、$T_d$、$T_h$、$T$ | — |
| `tetragonal_D4h` | $D_{4h}$ | 16 | LG1–LG64 | 16×16 |
| `hexagonal_D6h` | $D_{6h}$ | 24 | LG65–LG80 | 24×24 |

文件结构很简单：`data/group_operations.json` 是唯一机器可读数据源；`docs/group_theory.md` 汇总公式、全部矩阵、层群操作集合和乘法表；`group_tools.py` 用于查询和变换 POSCAR；`tests/test_repository.py` 负责数据校验。

## 使用

核心工具只需要 Python 3.10 或更高版本，没有第三方依赖。

```bash
# 查看数据集
python3 group_tools.py list

# 查询矩阵
python3 group_tools.py show '4+_001' --family tetragonal_D4h

# 查询层群
python3 group_tools.py group --family hexagonal_D6h --lg 80

# 查询点群
python3 group_tools.py group --family cubic_Oh --point-group Td

# 查询乘积
python3 group_tools.py multiply '4+_001' 2_100 --family tetragonal_D4h

# 对 POSCAR 坐标施加点操作
python3 group_tools.py apply-poscar structure.vasp 3+_001 \
  --family hexagonal_D6h --output transformed.vasp
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

当前测试覆盖矩阵正交性、六角基底变换、群闭合性、Bilbao 坐标作用、LG1–LG80 分类、全部乘法单元、单位元、逆元、结合律、JSON/Markdown 同步和 POSCAR 回写。

## 范围与后续

目前数据描述点操作的线性部分，不含非零平移的完整 Seitz 对 $(R\mid\mathbf t)$；POSCAR 变换也默认围绕原点。后续我会在这个基础上加入层间平移、旋转中心、双层堆叠对称性、铁电极化约束、高通量筛选、表示分解与不变量构造，并继续扩充分数量子铁电相关的群论数据。

参考资料的核对范围记录在 JSON 顶层 `sources` 字段中。正式公开发布时还需要确定代码、数据和文档的许可证，并补充 Bilbao 的具体页面与访问日期。
