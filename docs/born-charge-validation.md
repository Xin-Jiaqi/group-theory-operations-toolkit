# Born 有效电荷推导与实现的核验记录

核验日期：2026-11-09。目标仓库：`Xin-Jiaqi/group-theory-operations-toolkit`。通过 GitHub 读取的 `main` 快照为 `551a0daf8d4ba5a5f24b2811694bfe438616de48`，版本 `0.18.0`。本记录的“通过”仅指下列已经执行的检查，不表示完成新材料的 DFT 验证或全部上游集成测试。

## GitHub 实仓库与远端 CI 的后续复核（2026-11-09）

- **实仓库点群表，32/32 通过：** 通过已连接的 GitHub API 分别读取 `main` 中的 `data/crystallographic_point_groups.json` 与 `data/group_operations.json`，逐群建立九分量 `I-R⊗R` 约束，数值求秩并核对本功能分支 `docs/born-charge-reference.json` 的独立参数数目和张量基。全部 32 群通过；最大操作不变性残差为 `7.60e-13`。这是实际 registry 数据交叉核验，不应等同于执行完整 Python 仓库测试。
- **新增测试，165/165 通过：** 在同一会话的本地容器内重新运行 `python -m unittest discover -s tests -p test_born_effective_charges.py -q`，全部通过；独立精确代数验证脚本亦完成。
- **GitHub Actions 无法作为此次功能测试证据：** PR #33 的 `Verify repository` 工作流运行 `42611757474`，原始与重新运行两次的 Python 3.14 任务均约 2 秒结束，`runner_name` 为空、`steps=[]`；3.10 任务取消。说明未实际进入 runner 执行代码。GitHub job logs 下载接口返回 404，不能断言是账户额度、账单还是 GitHub runner 配置的具体哪一项导致。
- **仍未验证：** 完整仓库现有回归测试、项目 mypy、wheel 构建，以及外部 VASP/DFPT 数据。CI 红灯不得解释为新增代码的测试断言失败，也不得视为全套测试已通过。

## 文献交叉核验

| 来源 | 实际核对的位置或内容 | 本次推导/实现中的对应项 |
|---|---|---|
| VASP 官方 Born effective charges 文档 | 定义、原始 OUTCAR 与 py4vasp 的指标次序、DFPT/有限电场路线 | `Z[i,j]` 统一为 field-first；力取 `Z.T @ E`；读取原始矩阵不自动转置 |
| Octopus 官方 Born effective charges 文档 | 原子 BEC 一般不对称、位点群约束、电中性求和规则 | 九分量求解；区分位点稳定子与原子置换；显式 ASR |
| Ke, Liu, Liu, PRL 135, 046201 (2025) | 正文第 2–3 页及 Fig. 1–2 | 高对称 C3 构型的相关交叉耦合为零；层电荷求和；局部对称性破缺与畴壁驱动 |
| Wang, Dong, PRB 111, L201406 (2025) | 原版 PDF 第 2 页 Eqs. (1)–(2)、Fig. 1，以及第 3 页 Fig. 2 | 原文采用相反指标次序；端点为零与中间态非零不矛盾；区分各层合力 |
| Fachin 等, npj Computational Materials 11, 307 (2025) | Germanene；Symmetries of jacutingaite atoms；Eq. (3) | 3m 位点两参数；m 和 2/m 位点五参数；反演伙伴与 ASR 联合消零；质量加权模式投影 |

论文、DOI 和官方文档地址列于 [完整推导的参考文献](born-effective-charges.md#参考文献与可核查出处)。PDF 中涉及公式和图的部分以原版页面截图核对，而不是仅依赖提取文本。32 群的完整表和孤立双层模式的三类穷尽分类为独立推导，并非声称某篇论文已经逐项给出。

## 已执行的数学与软件检查

环境：Python 3.13.5，NumPy 2.3.5；精确代数核验另使用 SymPy 1.14.0。SymPy 仅用于可选验证脚本，不是新模块的运行时依赖。

| 检查 | 结果 | 覆盖范围 |
|---|---|---|
| 新增单元测试 | **165 项通过，0 失败** | 32 群每群四类检查，共 128 项；另有 37 项物理边界、文献实例、输入校验和接口测试 |
| 精确符号计算 | **32/32 通过** | 在含 sqrt(3) 的精确代数域上求九分量零空间，核对解析矩阵逐操作不变性、参数独立性与完备性 |
| 群平均/特征标交叉检查 | **32/32 通过** | 完整群闭包阶数；Reynolds 投影；投影幂等性；特征标维数；与生成元 SVD 零空间一致 |
| 坐标旋转检验 | **32/32 通过** | 随机正交变换下允许空间协变；比较投影算符而不是任意选取的 SVD 基矢 |
| 原子置换与中性约束 | **通过** | 局域 C3 与三原子轨道的区别；反演、水平镜面交换原子；联合特征标与 ASR 秩—零度关系 |
| 模式选择与虚功关系 | **通过** | 层交换符号；一维允许方向；层内抵消；非正交模式坐标；对角 BEC 驱动面外模式 |
| 索引、量纲及 OUTCAR | **通过** | 故意不对称的张量；V/nm 与 V/Å 换算；最后数据块、不完整块、NIONS、D 指数 |
| 命令行冒烟检查 | **通过** | C2x 滑移空间为 1 维；明确标注为合成数据的 OUTCAR 得到 zeta=(0.05,0)，ASR 为零 |
| Python 语法编译 | **通过** | 新模块、测试与脚本的 compileall |

独立群生成元位于 `tests/born_reference_groups.py`，仅用于验证与生成文档，不作为第二套运行时数据库。正式的 `point_group_born_space()` 延迟调用仓库原有点群 registry；轴向以其操作矩阵为准。

## 尚未执行的检查与外部操作

- **未执行完整仓库的 Python 回归测试、mypy、wheel 构建。** 现已通过已连接的 GitHub API 实际读取上游 32 点群 registry 并进行独立逐群数值/解析交叉核验；但 `scripts/check_born_registry.py` 尚未在真实完整 checkout 下运行，且 GitHub Actions 在 runner 启动前失败。
- 附带 `scripts/check_born_registry.py`，用于在完整仓库中将实际 32 群 registry 与独立解析表逐项比较，包括数字字符串与 Hermann–Mauguin 符号的区分。
- **没有进行新的 VASP/DFPT、有限电场或畴壁动力学计算。** 合成 OUTCAR 只测试解析器和模式投影，不是任何材料的预测数据。
- **远端发布记录：** 本增量包随后通过 GitHub 连接器提交到独立功能分支，未直接修改 `main`。远端 GitHub Actions 状态应以相应提交和 PR 页面为准。

## 复现命令

在完整仓库并具备既有运行/测试依赖的环境中：

```bash
export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"
python -m unittest discover -s tests -p 'test_born_effective_charges.py' -v
python scripts/verify_born_exact.py             # 可选审计，需 SymPy
python scripts/check_born_registry.py           # 需要完整仓库
python -m unittest discover -s tests -v         # 上游完整测试，应用补丁后执行
```

检查依据为明确的输入假设：静态绝缘体、固定晶胞、正交笛卡尔坐标、真实原子映射和正确模式表示。改变电边界条件、进入金属态、模式不封闭、或直接把晶体点群当成每个位点群，都超出当前表格的正确使用范围。
