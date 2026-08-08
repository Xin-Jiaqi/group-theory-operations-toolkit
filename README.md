# Group Theory Operations Toolkit

[![Verify repository](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml/badge.svg)](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/actions/workflows/verify.yml)
[![Release](https://img.shields.io/github/v/release/Xin-Jiaqi/group-theory-operations-toolkit)](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/releases)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

这是我为二维材料、层间堆叠、铁电与非线性光学研究维护的群论数据工具。它把晶体学与磁性群注册表、点操作矩阵、群乘法、光场二次表示和允许张量基放在同一套可验证的机器接口中，可直接查阅，也可作为高通量筛选与材料工作流的底层依赖。

当前版本为 **0.8.0**。仓库回答四类问题：一个操作怎样作用与复合；某个（磁）点群允许哪些非线性响应分量；时间反演怎样区分普通与磁性响应；给定单层对称性与层间平移后，堆叠结构允许何种极化与切换关系。

## 科学能力

| 层次 | 覆盖范围 | 物理用途 |
|---|---:|---|
| 点操作代数 | 88 个基础记录 | 查询 $O_h$、$D_{4h}$、$D_{6h}$ 中的分数/笛卡尔矩阵、逆操作与乘法表 |
| 晶体学群 | 32 点群、230 空间群、80 层群 | 固定标准 setting、生成元与 Seitz 闭包，连接三维周期晶体和二维周期层状体系 |
| 磁性群 | 122 磁点群、528 磁性层群 | 显式保存 $(R,\theta)$ 中的时间反演标签；区分 I–IV 型和 type-IV 反平移 |
| 二次光场与响应 | 88 个 $M_\pm(R)$；全部注册群的允许基 | 求 shift current、circular injection current、SHG 及通用时间奇偶张量的对称性允许空间 |
| 堆叠铁电 | 80 层群与五类二维 Bravais 晶格 | 判断单层极化类型、双层取向类、等价界面、极化切换和递归多层保留对称性 |

机器数据位于 [`data/`](data)，对应 JSON Schema 位于 [`schema/`](schema)。[`data/group_operations.json`](data/group_operations.json) 是基础点操作的唯一事实源；常用入口包括 [`data/crystallographic_point_groups.json`](data/crystallographic_point_groups.json) 与 [`data/optical_response_invariants.json`](data/optical_response_invariants.json)。派生表由版本化脚本确定性生成，并记录上游数据或输入文件的 SHA-256。

## 物理框架

### 从光场表示到允许张量

矩阵采用列向量约定 $\mathbf r'=D(R)\mathbf r$。线偏振相关的对称二次场位于 $\operatorname{Sym}^2D$，在 `(xx, yy, zz, xy, xz, yz)` 基底上表示为

$$M_+(R)=\operatorname{Sym}^2D(R).$$

圆偏振相关的反对称场可写成实轴矢量 $\mathbf h=i\mathbf E\times\mathbf E^*$；极性与轴性在非固有操作下相差一个行列式，因此

$$M_-(R)=\det[D(R)]D(R).$$

对极性输出矢量，Neumann 原理转化为等变映射方程

$$D(R)T_+=T_+M_+(R),\qquad D(R)T_-=T_-M_-(R).$$

求解所有群操作的公共零空间即可得到完整的对称性允许张量基。`shift_current` 与 `shg` 使用 $3\times6$ 的 $D\leftarrow M_+$ 映射；`circular_injection_current` 使用 $3\times3$ 的 $D\leftarrow M_-$ 映射。磁性接口把操作扩展为 $(R,\theta)$，再按输入/输出的时间奇偶性处理反幺正约束，从而区分 normal/magnetic shift current、normal/magnetic injection current 与 i-type/c-type SHG。

### 从层群到堆叠铁电

单层允许极化是其点操作共同不变子空间

$$V_P=\{\mathbf P\mid D(R)\mathbf P=\mathbf P,\ \forall R\in G\},$$

由此得到面内（IP）、面外（OP）、耦合（CP）和非极性（NP）四类。双层相对取向用 Bravais 晶格点群对单层点群的**左陪集**分类，因此不要求单层子群为正规子群。对称等价界面满足

$$R^{\pm}\boldsymbol\tau_p\equiv\pm\boldsymbol\tau_q\pmod{\Lambda},$$

其中 $R^+$ 保持层次，$R^-$ 交换上下层。递归多层中，$R^+$ 必须分别固定相邻界面平移，$R^-$ 必须交换它们；程序同时处理 primitive 与 centered 面内晶格。

这些结论只回答“对称性是否允许”和“两个构型是否由指定操作连接”。能量简并、滑移势垒、极化大小、频率依赖与动力学仍需第一性原理、模型或实验确定。

## 安装与快速查询

需要 Python 3.10 或更高版本，核心运行依赖只有 NumPy。

```bash
python -m pip install -e .
group-ops validate

# 点操作、群注册表与 M_-/M_+
group-ops multiply '4+_001' 2_100 --family tetragonal_D4h
group-ops field-representation m_100 --family tetragonal_D4h --space antisymmetric --json
group-ops point-groups 4mm --json
group-ops layer-groups p4/nmm --json
group-ops magnetic-point-groups "6'/m'mm'" --json
group-ops magnetic-layer-groups 80.9.528 --json

# 非线性光学与堆叠铁电
group-ops invariants 4mm shift_current --json
group-ops magnetic-layer-responses 6.5.25 shg_odd --json
group-ops layer-polarity LG68 --json
group-ops stacking-rotations -1 --lattice rP --json
```

Python 接口返回带类型的不可变记录：

```python
from group_theory_operations import (
    equivalent_interface_orbit,
    layer_group_polarization,
    magnetic_layer_response_tensor_basis,
    point_group_operations,
)

print(layer_group_polarization(68).polar_type)  # NP
print(magnetic_layer_response_tensor_basis("6.5.25", "shg_odd").dimension)
print(equivalent_interface_orbit((1 / 3, 2 / 3), point_group_operations("6/mmm")))
```

POSCAR/CIF 读写由可选依赖 `materials-structure-core[io]` 提供；本仓库只负责经过约定检查的对称操作：

```bash
group-ops apply-structure structure.vasp 3+_001 \
  --family hexagonal_D6h --output transformed.vasp
```

## 版本演进：每一版增加了什么物理能力

| 版本 | 新增内容 | 科学意义与核验 |
|---|---|---|
| **0.1.0** | 统一操作目录、矩阵、逆操作、乘法表、JSON Schema 与 Python/CLI 接口 | 建立后续所有表示与群闭包的代数事实源；逐项核验正交性、逆元和乘法一致性 |
| **0.2.0** | 为全部 88 个操作生成 $M_+(R)$ 与 $M_-(R)$ | 把空间操作提升到二次光场空间，明确线偏振对称场与圆偏振轴矢量的不同变换律；验证直接光场作用、叉积恒等式和表示同态 |
| **0.3.0** | 补齐 32 个晶体学点群，并加入 shift current、SHG、circular injection current 不变量求解器 | 将 Neumann 原理变成可执行的公共不变子空间问题；每个输出基向量均逐群操作验证 |
| **0.4.0** | 注册 230 空间群的 530 个 Hall setting 与 80 层群的 116 个 layer Hall setting | 从点群扩展到含平移的 Seitz 代数，覆盖三维晶体与二维层状体系；修正非正交分数基下的逆运算并与 spglib 交叉验证 |
| **0.4.1** | 修复一次性迭代器的闭包问题，并复核 LG1–LG80 的标准旋转集合 | 这是数值与数据一致性维护版，保证同一群在重复调用和不同容器输入下给出稳定结果 |
| **0.5.0** | 注册全部 122 个磁点群，加入 type I、灰群、黑白群及时间奇偶张量映射 | 对称操作从 $R$ 扩展为 $(R,\theta)$，使磁畴反转和反幺正约束进入张量选择定则 |
| **0.6.0** | 加入六个磁性非线性光学扇区与复反幺正映射求解器 | 区分 NSC/NIC 与 MSC/MIC，并把 SHG 分为时间偶 i-type 和时间奇 c-type；以灰群、$PT$ 群和 $\bar3'm'$ 规则回归测试 |
| **0.7.0** | 注册全部 528 个磁性层群及六类响应基 | 把磁性张量筛选扩展到二维周期体系；保存有限磁点余群、I–IV 型、母层群、对应磁空间群和 type-IV 反平移 |
| **0.8.0** | 加入双层与多层堆叠铁电对称性内核 | 从 80 层群导出 IP/OP/CP/NP；用左陪集处理取向，用 $R^\pm$ 处理等价界面与极化切换，并实现递归多层判据；复现 BN AB/BA、graphene ABC/$D_{3d}$ 与 ABA/$D_{3h}$ 基准 |

完整开发门槛与后续边界见 [`ROADMAP.md`](ROADMAP.md)。版本号描述的是机器接口和科学能力的演进，不代表已经计算任何具体材料的响应强度。

## 理论来源与文献要点

1. **晶体学注册表。** 空间群和层群的 Hall setting、Seitz 操作与标识基于固定版本的 [spglib](https://spglib.readthedocs.io/) 数据库；层群元数据固定为 v2.5.0，并保留许可证、来源 URL 与校验值。spglib 在本仓库中是独立的数据与结构识别基准，不是本工具对称性推导的替代品。

2. **磁性层群。** D. B. Litvin, *Magnetic Group Tables, Part 3: Magnetic Layer Groups* 提供 OG 符号、磁点群和操作表；Z. Zhang *et al.*, [“Encyclopedia of emergent particles in 528 magnetic layer groups and 394 magnetic rod groups,” *Phys. Rev. B* **107**, 075405 (2023)](https://doi.org/10.1103/PhysRevB.107.075405) 提供 528 个磁性层群的类型及其与磁空间群的系统对应。仓库还利用 spglib 磁空间群数据库进行第三方交叉核验。

3. **磁性非线性光电流。** H. Wang and X. Qian, [“Electrically and magnetically switchable nonlinear photocurrent in $PT$-symmetric magnetic topological quantum materials,” *npj Computational Materials* **6**, 199 (2020)](https://doi.org/10.1038/s41524-020-00462-9) 系统区分 NSC、NIC、MSC 与 MIC，并说明这些响应对空间反演、时间反演和磁畴的不同奇偶性。0.6.0 将这种物理分类落实为显式时间奇偶的群论接口；仓库只求允许分量，不实现论文中的微观 Berry 几何或第一性原理响应公式。

4. **磁性 SHG。** R.-C. Xiao *et al.*, [“Classification of second harmonic generation effect in magnetically ordered materials,” *npj Quantum Materials* **8**, 62 (2023)](https://doi.org/10.1038/s41535-023-00594-3) 将 SHG 分解为时间偶 i-type 与时间奇 c-type，并用磁点群连接磁序、反演破缺来源和 SHG 选择定则。0.6.0/0.7.0 采用这一时间奇偶框架求允许张量基，但不复刻论文的七类材料学分类或材料数据库。

5. **双层堆叠铁电。** J. Ji *et al.*, [“General Theory for Bilayer Stacking Ferroelectricity,” *Phys. Rev. Lett.* **130**, 146801 (2023)](https://doi.org/10.1103/PhysRevLett.130.146801) 建立 80 层群上的一般 BSF 框架，说明堆叠可在原本非极性的单层之间产生极化，并给出 IP/OP/CP/NP 与 $R^\pm$ 的基本语言。0.8.0 的层群极化固定空间和双层操作分类以此为基础。

6. **切换与取向的增强理论。** J. Xin, Y. Guo, and Q. Wang, [“Enhanced theoretical framework for bilayer stacking ferroelectricity,” *Phys. Rev. B* **111**, 224102 (2025)](https://doi.org/10.1103/PhysRevB.111.224102) 修正旋转操作的选择：应从晶格系统允许的操作出发，并用左陪集区分物理取向，而非默认商群；同时把“极性构型”与“存在对称等价反向极化态的铁电构型”区分开。0.8.0 据此实现非正规子群也适用的左陪集和等价界面映射。

7. **多层堆叠铁电。** J. Xin, Y. Guo, and Q. Wang, [“Multilayer stacking ferroelectricity in two-dimensional materials with Bravais lattice symmetry: Theory and applications,” *Phys. Rev. B* **113**, 075310 (2026)](https://doi.org/10.1103/9tt5-qm26) 把 BSF 推广到 Bravais-lattice monolayers 与递归多层结构：层保持操作固定相邻界面平移，层交换操作交换它们，由此决定多层堆叠后的剩余对称性和允许极化。0.8.0 实现的是该递归对称性判据，不包含文中的材料设计、BPVE 数值或能量计算。

## 约定、边界与核验

- 乘法表使用“行操作 × 列操作”，即 $D(\text{left})D(\text{right})$，右操作先作用。
- 分数坐标矩阵只用于相应晶格基；轴矢量公式 $M_-(R)=\det(D)D$ 由正交笛卡尔矩阵导出。
- 磁性层群接口面向 $q=0$ 均匀张量，保存有限磁点余群。普通平移不进入均匀张量约束；type-IV 反平移单独保留，不能把该接口误作完整的仿射磁性层群作用。
- 输出是对称性允许空间，不包含材料数值、单位、共振、耗散、频率置换、界面能量或切换动力学。

```bash
python3 -m unittest discover -s tests -v
```

测试覆盖矩阵与逆运算、乘法表、Seitz 闭包、$M_\pm$ 表示同态、32 点群、122 磁点群、230 空间群、80 层群、528 磁性层群、六类磁性光学响应及堆叠铁电基准。外部论文原文与截图不收入仓库；仓库只保存可复现派生数据、正式链接、许可信息和校验值。

本项目采用 [BSD 3-Clause License](LICENSE)。科研使用请通过 [`CITATION.cff`](CITATION.cff) 引用实际使用的软件版本，并同时引用与你调用的理论模块对应的上述原始文献。
