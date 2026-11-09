# Born 有效电荷与双层滑移模式的对称性

新增文档与增量模块，现有 `group-ops` CLI、点群数据库、包初始化及项目配置不变。

- [完整推导](../born-effective-charges.md)：原子位点与原子置换、32 群六类形式、一般模式表示、层交换、孤立双层三类耦合、切换路径与 VASP 约定。
- [速查表](../born-charge-quick-reference.md)：32 群逐项表以及滑移耦合方向。
- [可机读表](../born-charge-reference.json)：显式生成元、张量形式、维数、允许掩码与模式空间。
- [核验记录](../born-charge-validation.md)：文献对应、165 项新增测试与 32 群精确代数检查；明确区分未运行的完整上游集成测试。

使用 `from group_theory_operations.born_effective_charges import ...` 或 `python -m group_theory_operations.born_effective_charges ...`。模块包括原子/结构/模式允许空间、原始 OUTCAR 读取、层模式投影与 ASR 残差。

示例文件位于 `examples/born/`。其中 `OUTCAR.synthetic` 是人工构造的测试数据，不是第一性原理结果。
