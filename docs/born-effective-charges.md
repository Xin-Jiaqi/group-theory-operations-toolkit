# Born 有效电荷的对称性与双层滑移电场耦合

本文给出静态 Born 有效电荷的完整局域对称性推导、32 个晶体位点点群的六类张量形式，以及面外电场驱动孤立双层面内相对滑移的判据。对应实现为 `src/group_theory_operations/born_effective_charges.py`；速查表见 [born-charge-quick-reference.md](born-charge-quick-reference.md)，检验记录见 [born-charge-validation.md](born-charge-validation.md)。

**结论的对象必须明确：单原子的张量、多个原子组成的结构、指定集体模式，是三个不同层次。只知道晶体点群，不能唯一给出每个原子的 BEC。对称性给出的是允许的线性响应空间，不给出响应大小、能垒或可切换性的充分条件。**

## 1. 定义、坐标与适用范围

考虑绝缘体系、固定晶胞、零宏观电场附近的静态绝热响应。原子子晶格位移为 $u_{\kappa j}$，晶胞体积为 $\Omega$，$e>0$。统一使用 VASP **原始 OUTCAR** 的指标次序 [1]：

$$
Z^*_{\kappa,ij}
=\frac{\Omega}{e}\left.\frac{\partial P_i}{\partial u_{\kappa j}}\right|_{\mathcal E=0}
=\frac1e\left.\frac{\partial F_{\kappa j}}{\partial\mathcal E_i}\right|_{\mathcal E=0}.
\tag{1}
$$

第一个指标是极化/电场，第二个是位移/受力；数值以 $|e|$ 为单位。星号不是复共轭。所有矩阵用正交笛卡尔坐标和列向量。$\kappa$ 表示原胞中的一个原子及其平移等价子晶格，扰动波矢为 $\mathbf q=0$。

令 $\mathbf p=\Omega\mathbf P$，局域地选择连续的极化分支，则

$$
\delta\mathbf p=e\sum_\kappa Z^*_\kappa\delta\mathbf u_\kappa,
\qquad
\Delta\mathbf F_\kappa=e(Z^*_\kappa)^T\boldsymbol{\mathcal E}+O(\mathcal E^2).
\tag{2}
$$

力公式中的转置不可省略。尤其 $\mathcal E_z\to F_x$ 对应 $Z_{zx}$，在 OUTCAR 中是第 3 行第 1 列。VASP 官方说明 **py4vasp 使用相反的 force-first 约定**，两者不能直接混用 [1]。

### 1.1 为什么没有 $Z_{ij}=Z_{ji}$？

电焓的局域展开为

$$
\mathcal F(\{\mathbf u\},\boldsymbol{\mathcal E})
=U_0(\{\mathbf u\})-\boldsymbol{\mathcal E}\cdot\mathbf p_0(\{\mathbf u\})
-\tfrac12\mathcal E_i\alpha_{ij}(\{\mathbf u\})\mathcal E_j+\cdots.
$$

因此 $eZ_{\kappa,ij}=-\partial^2\mathcal F/(\partial\mathcal E_i\partial u_{\kappa j})$。交换求导顺序，只能交换“对同一个 $\mathcal E_i$ 和同一个 $u_{\kappa j}$ 的先后”；不能把它变成 $-\partial^2\mathcal F/(\partial\mathcal E_j\partial u_{\kappa i})$。BEC 一般是九分量的实矩阵，而不是预先对称的六分量矩阵 [2]。

### 1.2 本文不直接覆盖的情况

金属的非绝热动态有效电荷、有限频率复数响应、有限波矢声子、固定 $D$ 与固定 $E$ 边界条件的转换，不属于本文的静态绝缘体表格。磁性体系的静态电响应可在实际磁结构的对称操作下分析，但不得直接用破坏该磁结构的非磁母相操作。有限外场时应采用该外场下仍保留的对称群。

## 2. 从空间操作推导原子 BEC 的变换

设当前结构的空间对称操作 $g=\{R_g|\mathbf t_g\}$ 将原子 $\kappa$ 映到 $\pi_g(\kappa)$，允许相差晶格平移。微小位移及偶极变化都是极矢量：

$$
\delta\mathbf p'=R_g\delta\mathbf p,
\qquad
\delta\mathbf u'_{\pi_g(\kappa)}=R_g\delta\mathbf u_\kappa.
$$

将其代入响应关系，并利用任意 $\delta\mathbf u_\kappa$ 的独立性，得到

$$
R_gZ^*_\kappa=Z^*_{\pi_g(\kappa)}R_g,
\qquad
\boxed{Z^*_{\pi_g(\kappa)}=R_gZ^*_\kappa R_g^T.}
\tag{3}
$$

没有 $\det R_g$ 因子：两个指标都是极矢量指标。这里平移不直接出现在 $\Gamma$ 点的笛卡尔张量变换中，但**决定原子映射 $\pi_g$，不能在建立映射之前删除**。

### 2.1 位点稳定子与等价原子

位点稳定子为

$$
G_\kappa=\{g:\pi_g(\kappa)=\kappa\}.
$$

仅当 $g\in G_\kappa$ 时，式 (3) 才成为单个原子的自约束

$$
\boxed{Z^*_\kappa=R_gZ^*_\kappa R_g^T.}
\tag{4}
$$

若 $g$ 交换 A、B，约束是 $Z_B=R_gZ_AR_g^T$，不是 $Z_A=R_gZ_AR_g^T$。例如 $C_3$ 循环交换三个一般位置原子时，单个原子可以有非零 $Z_{zx},Z_{zy}$；三原子完整轨道的和才必须符合整体 $C_3$ 不变性。代码提供该反例的回归测试。

### 2.2 声学求和规则

对中性体系整体平移所有原子，极化不变，因此 [2]

$$
\boxed{\sum_\kappa Z^*_{\kappa,ij}=0.}
\tag{5}
$$

这约束的是**全部原子之和**，不是要求每个原子张量为零。对一个只有两个反演相关原子的中性原胞，反演给出 $Z_A=Z_B$，求和规则又给出 $Z_A+Z_B=0$，二者合起来确实迫使二者都为零。germanene 的文献分析明确给出了这一例子 [5]。

因此“晶体整体的 Born 张量”不能简单当作一个一般二阶材料张量来查表：中性原胞的全部原子 BEC 直接求和恒为零。

## 3. 九分量线性方程与两种独立检验

采用按行展开

$$
\mathbf z=(Z_{xx},Z_{xy},Z_{xz},Z_{yx},Z_{yy},Z_{yz},Z_{zx},Z_{zy},Z_{zz})^T.
$$

式 (4) 化为

$$
\boxed{[I_9-R_g\otimes R_g]\mathbf z=0.}
\tag{6}
$$

将生成元的约束纵向堆叠，零空间即允许张量空间。生成元已经足够，因为对每个生成元不变就对其任意乘积不变。实现使用 SVD，保留九个未知量；返回正交基和投影矩阵，而不把 SVD 的任意基选择误当作物理参数的唯一表示。

独立检验可用完整有限群的 Reynolds 投影：

$$
\Pi=\frac1{|G|}\sum_{g\in G}R_g\otimes R_g,
\quad \Pi^2=\Pi,
\quad
\boxed{d=\operatorname{tr}\Pi=\frac1{|G|}\sum_g[\operatorname{tr}R_g]^2.}
\tag{7}
$$

本补丁另外用 $\mathbb Q(\sqrt3)$ 上的精确符号零空间检验全部 32 个点群；不是仅比较浮点小数是否“接近零”。

### 3.1 实际多原子结构

将所有 $\mathbf z_\kappa$ 拼接为 $9N$ 个未知量，对每个原子施加

$$
\mathbf z_{\pi_g(\kappa)}-(R_g\otimes R_g)\mathbf z_\kappa=0,
$$

再附加式 (5)。这样同时处理位点约束、等价原子关系和 ASR，不会把宏观点群错误地强加给每个原子。

当完整有限操作集合及原子置换形成一致表示时，未加 ASR 的维数还可独立检查：

$$
d_{\rm atoms}=\frac1{|G|}\sum_g n_{\rm fixed}(g)[\operatorname{tr}R_g]^2.
\tag{8}
$$

其中 $n_{\rm fixed}(g)$ 是模晶格平移不变的原子数。加中性 ASR 后

$$
d_{\rm neutral}=d_{\rm atoms}-d_{\rm tensor}(G).
\tag{9}
$$

证明：求和映射把满足原子协变关系的张量组映到整体不变二阶张量；该映射满射，因为任意整体不变 $X$ 可由每个原子取 $Z_\kappa=X/N$ 实现。ASR 空间就是该满射的核，故由秩—零度定理得到式 (9)。这只计算对称性及中性约束，不包含额外微观数值关系。

## 4. 为什么全部 32 个点群只得到六种形式

各字母为独立实参数；不同矩阵中的同名字母无联系。本节使用**主轴沿 $z$，单斜镜面为 $xy$，次级二重轴沿 $x$**的参考设置，不能只凭群名忽略实际取向。

| 类别 | 点群 | 允许矩阵 | 维数 |
|---|---|---|---:|
| A | $1,\bar1$ | $\begin{pmatrix}a&b&c\\d&e&f\\g&h&i\end{pmatrix}$ | 9 |
| B | $2,m,2/m$ | $\begin{pmatrix}a&b&0\\c&d&0\\0&0&e\end{pmatrix}$ | 5 |
| C | $222,mm2,mmm$ | $\operatorname{diag}(a,b,c)$ | 3 |
| D | $3,\bar3,4,\bar4,4/m,6,\bar6,6/m$ | $\begin{pmatrix}a&b&0\\-b&a&0\\0&0&c\end{pmatrix}$ | 3 |
| E | $32,3m,\bar3m,422,4mm,\bar42m,4/mmm,622,6mm,\bar6m2,6/mmm$ | $\operatorname{diag}(a,a,c)$ | 2 |
| F | $23,m\bar3,432,\bar43m,m\bar3m$ | $aI_3$ | 1 |

### 4.1 A、B、C 类

反演 $R=-I$ 对二阶极性张量的作用为 $(-I)Z(-I)=Z$，所以 $1$ 与 $\bar1$ 都不给出局域分量限制。

对 $C_{2z}=\operatorname{diag}(-1,-1,1)$，$Z_{ij}=r_ir_jZ_{ij}$。恰有一个指标为 $z$ 的项变号，故 $Z_{xz},Z_{yz},Z_{zx},Z_{zy}=0$；面内 $2\times2$ 块仍一般不对称。$m_z=-C_{2z}$ 对二阶张量作用完全相同，加入反演也不改变结果，得到 B 类。

再加入独立的面内二重轴或相应镜面，$xy,yx$ 也被禁止，只剩三个独立对角元，得到 C 类。

### 4.2 D、E 类

对绕 $z$ 的 $n=3,4,6$ 重旋转，写

$$
R=\begin{pmatrix}A&0\\0&1\end{pmatrix},\quad
Z=\begin{pmatrix}M&v\\w^T&c\end{pmatrix}.
$$

不变性要求 $v=Av,w=Aw,M=AMA^T$。非平凡的这些面内旋转无非零不动矢量，故 $v=w=0$；$M$ 与 $A$ 对易，得到

$$
M=aI_2+b\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
$$

旋转反演 $-R$ 与 $R$ 对该二阶张量作用相同，产生 D 类所有成员。若再有竖直镜面或面内二重轴，面内反对称部分变号，故 $b=0$，得到 E 类。**$C_3$ 只禁止面外—面内交叉分量，并不禁止所有非对角元。**

### 4.3 F 类

立方群中的二重轴先消去非对角元，绕 $[111]$ 的三重旋转循环交换三个轴，再要求三个对角元相同。因此所有五个立方点群均只有 $aI_3$。

另一种一致性解释是 $Z=(\operatorname{tr}Z/3)I+S_0+A$：标量、对称无迹和反对称三部分分别约束；反对称部分等价于轴矢量，只有相应轴矢量被点群允许时才可保留。

### 4.4 坐标变化不能省略

若 $Q$ 将参考坐标旋转到实际笛卡尔坐标，必须同步使用

$$
R'=QRQ^T,\qquad Z'=QZQ^T.
$$

例如 B 类若实际独特轴是 $x$，矩阵变成

$$
Z=\begin{pmatrix}a&0&0\\0&b&c\\0&d&e\end{pmatrix},
$$

此时 $Z_{zy}$ 允许非零。对于非正交分数坐标操作 $W$，若晶格矢量为列组成 $L$，须先转换 $R=LWL^{-1}$；不能把 $W^T$ 当作 $W^{-1}$。

## 5. 任意集体模式的有效电荷

用 $m$ 个广义坐标 $q_a$ 定义位移模式：

$$
\delta u_{\kappa j}=\sum_a U_{\kappa j,a}\delta q_a.
$$

代入式 (2)：

$$
\delta p_i=e\sum_a C_{ia}\delta q_a,
\qquad
\boxed{C_{ia}=\sum_{\kappa j}Z^*_{\kappa,ij}U_{\kappa j,a}.}
\tag{10}
$$

质量加权本征模的 $1/\sqrt{M_\kappa}$ 必须包含在 $U$ 中，不能和刚性滑移坐标的归一化混用。文献的红外振子强度就是这一模式投影的一种具体归一化 [5]。

若模式子空间在对称操作下封闭，$\delta q'=D_g\delta q$，由 $\delta p'=R_g\delta p$ 得

$$
\boxed{R_gC=CD_g,\qquad C=R_gCD_g^{-1}.}
\tag{11}
$$

实际构造 $D_g$ 时，先将 $U$ 排成 $3N\times m$ 矩阵，令 $T_g$ 表示完整原子位移空间的操作：它同时将原子 $\kappa$ 置换到 $\pi_g(\kappa)$，并将该原子的三维位移乘以 $R_g$。模式子空间封闭的条件是

$$
T_gU=UD_g.
$$

对满列秩的模式矩阵，可计算 $D_g=(U^TU)^{-1}U^TT_gU$，但必须检查 $T_gU-UD_g=0$ 的残差；否则这些模式不组成封闭子空间，应补齐被对称操作混入的模式，而不能直接沿用投影所得的矩阵。当前实现接收已经建立的 $D_g$，不自动从任意结构识别模式。

这才是一般的模式选择定则，而非对所有模式都使用 $C=RCR^T$。按行展开后的约束为

$$
[R_g\otimes I_m-I_3\otimes D_g^T]\operatorname{vec}_r C=0.
\tag{12}
$$

相应的广义场致力由虚功关系给出

$$
\boxed{\mathbf f_q^{(E)}=eC^T\boldsymbol{\mathcal E}.}
\tag{13}
$$

若模式包括面外原子运动，对角元 $Z_{zz}$ 也能贡献。因此“所有电场切换都必须有原子 BEC 非对角元”并不成立。

## 6. 双层面内相对滑移

### 6.1 坐标、层电荷和归一化

以下限定：自由、孤立的中性双层，层平面为 $xy$，固定原胞，刚性平移每一层；$z$ 为实际层法向。上下层以空间位置标记。

定义

$$
\mathbf s=\mathbf u_T-\mathbf u_B,
\qquad \delta\mathbf u_T=\tfrac12\delta\mathbf s,
\qquad \delta\mathbf u_B=-\tfrac12\delta\mathbf s.
$$

这里固定的是两层位移的算术平均，不要求两层质量相等。总中性下换用另一种固定整体平移的规范，不改变相对模式的电荷。

面外偶极响应为

$$
\delta p_z=e\boldsymbol\zeta^T\delta\mathbf s,
\qquad
\boxed{\zeta_a=\frac12\left(\sum_{\kappa\in T}Z_{\kappa,za}-\sum_{\kappa\in B}Z_{\kappa,za}\right),\quad a=x,y.}
\tag{14}
$$

对仅含这两层的中性体系，ASR 进一步给出

$$
\zeta_a=\sum_{\kappa\in T}Z_{\kappa,za}=-\sum_{\kappa\in B}Z_{\kappa,za}.
\tag{15}
$$

因此

$$
\boxed{\mathbf f_s^{(E)}=e\mathcal E_z\boldsymbol\zeta.}
\tag{16}
$$

层内原子可以互相抵消；仅见某个原子的非零交叉元，不能断定有净层间滑移驱动力。实际计算应同时报告 ASR 残差，而不是无提示地把它清零。PRL 的上层 B、N 交叉 BEC 求和公式与式 (15)–(16) 一致 [3]。

### 6.2 层交换的负号

将保持层平面的操作写成

$$
R_g=\begin{pmatrix}A_g&0\\0&\eta_g\end{pmatrix},\qquad \eta_g=\pm1.
$$

另外定义层标签的置换奇偶性：$\chi_g=+1$ 表示分别保留两层，$\chi_g=-1$ 表示交换上下层。于是

$$
\delta p_z'=\eta_g\delta p_z,
\qquad
\boxed{\delta\mathbf s'=\chi_g A_g\delta\mathbf s.}
\tag{17}
$$

第二式的 $\chi_g$ 来自 $u_T-u_B$ 在层交换下反号。这不是普通面内位移矢量的变换，漏掉它会得到错误选择定则。

将式 (17) 代入式 (14)，得到一般条件

$$
\boxed{(\chi_gA_g)^T\boldsymbol\zeta=\eta_g\boldsymbol\zeta.}
\tag{18}
$$

对上述**孤立双层**，翻转法向恰好交换两层，因此 $\chi_g=\eta_g$，约束化简为

$$
\boxed{A_g^T\boldsymbol\zeta=\boldsymbol\zeta
\quad\Longleftrightarrow\quad A_g\boldsymbol\zeta=\boldsymbol\zeta.}
\tag{19}
$$

这给出的是当前构型的局域 $\mathcal E_z\to f_{s,x},f_{s,y}$ 一阶耦合。周期多层、带衬底结构、不同标签定义，必须重新建立实际模式表示，不能自动令 $\chi=\eta$。

### 6.3 完整分类只有三种结果

令 $H=\{A_g\}\subset O(2)$。允许的 $\boldsymbol\zeta$ 就是 $H$ 的不动子空间。

若存在非零不动向量 $v$，每个二维正交矩阵既然固定 $v$，在垂直于 $v$ 的一维空间上只能取 $+1$ 或 $-1$。因此 $A_g$ 只能是恒等或以 $v$ 为轴的同一个面内反射。这证明：

| $H$ | 允许空间 | 对应的孤立双层点群及取向 |
|---|---|---|
| $\{I_2\}$ | 整个面内平面，2 维 | $C_1$；仅含水平镜面的 $C_s^{(h)}$ |
| $\{I_2,S_\ell\}$ | 沿面内直线 $\ell$，1 维 | 竖直镜面 $C_s^{(v)}$；二重轴沿 $\ell$ 的 $C_2$；包含 $m_z$ 和相应竖直镜面的、主轴沿 $\ell$ 的 $C_{2v}$ |
| 其他有限面内投影群 | 零空间，0 维 | 其余保持该层平面的点群/取向 |

这是在上述假设下的穷尽分类，不是只列几个常见操作。特别是任意非平凡的面内旋转都使 $\boldsymbol\zeta=0$；两个不同方向的面内反射相乘给出非平凡旋转，也使其为零。

### 6.4 两个容易混淆的反例

**水平镜面 $m_z$：** 若它固定一个原子位点，会禁止该原子的 $Z_{zx},Z_{zy}$；若它交换上下层，则只要求 $Z^B_{za}=-Z^T_{za}$，相对滑移的模式电荷反而可以非零。该高对称参考态的 $p_z$ 为零，与其滑移导数非零不矛盾。

**反演 $I$：** 它不限制单个位点二阶极性张量，但在孤立双层相对模式中 $A=-I_2,\chi=\eta=-1$，故 $\delta s'=\delta s$ 而 $\delta p_z'=-\delta p_z$，迫使 $\boldsymbol\zeta=0$。

## 7. 连接两个构型的操作不是单态约束

若 $g$ 将不同构型 I、II 联系起来，且相应局域模式坐标满足 $\delta q_{II}=D_g\delta q_I$，则

$$
\boxed{C_{II}=R_g C_I D_g^{-1}.}
\tag{20}
$$

除非 I、II 本来就是同一个构型，否则不能将其改写为 $C_I=R_gC_ID_g^{-1}$。对于上述孤立双层，伙伴态的面外滑移电荷满足 $\boldsymbol\zeta_{II}=A_g\boldsymbol\zeta_I$；一般情况则为 $\boldsymbol\zeta_{II}=\eta_g\chi_g A_g\boldsymbol\zeta_I$。常数堆叠平移只改变坐标原点，不改变局域切向表示。

因此在 BSF 工作流中应分开保存：当前结构的 stabilizer，用于允许性；构型间的 connection operator，用于两态间的响应关系。

## 8. 从 BEC 到切换路径：能推断什么

沿连续绝缘路径 $\mathbf u(q)$：

$$
\frac1e\frac{dp_z}{dq}=\sum_{\kappa j}Z^*_{\kappa,zj}\frac{du_{\kappa j}}{dq},
\qquad
\Delta p_z=e\int C_{zq}(q)\,dq.
\tag{21}
$$

对严格面内刚性滑移，这等于 $\boldsymbol\zeta\cdot d\mathbf s/dq$。若两端偶极不同，则路径上不可能处处为零，但端点的导数完全可以为零；计算端点 BEC 不能替代整条路径。Wang–Dong 的 Fig. 1(d) 和 Ke 等的高对称态分析与这一点一致 [3,4]。

**若对每个滑移点还弛豫了层间距、起伏或层内结构，式 (21) 必须包含所有 $j=x,y,z$。** 这时“弛豫后极化曲线的导数”不必等于刚性上层 $Z_{zx},Z_{zy}$ 求和；$Z_{zz}\,du_z/dq$ 等项也参与。应比较相同约束条件的 DFPT 与有限差分结果。

在最低阶场近似下

$$
G(q,\mathcal E_z)=U_0(q)-\mathcal E_zp_z(q),
\quad
f_q=-U'_0(q)+\mathcal E_zp'_z(q),
\quad
G''=U''_0-\mathcal E_zp''_z.
\tag{22}
$$

在精确 $C_{3z}$ 构型，纯 $\mathcal E_z$ 保留该旋转对称性，不能选出非零定向面内力；只要该高对称结构及电子态仍被保持，这一论证不只限于一阶场。可是能量曲率仍可能随场改变；若发生不稳定、热涨落、成核、缺陷或畴壁引起的对称性破缺，应在新构型重新判断。这不是对所有外场切换机制的“不可能性定理”。

例如 $p_z(q)=p_0\cos(\pi q/L)$ 的两端导数为零，中间导数非零，而两端偶极相反。再如 $C_3$ 允许 $p_z(\mathbf s)=p_0+a|\mathbf s|^2+O(s^3)$：梯度在原点为零，但 $-\mathcal E_zp_z$ 改变滑移模刚度。这两个例子是解析示意，不是 h-BN 的拟合模型。

非零 BEC 只说明某种电耦合存在；是否克服恢复力、沿目标路径做功、形成畴壁、出现可观测翻转，仍需要能量和动力学计算。

## 9. VASP 数值核验路线

对固定的绝缘构型，通常可在已收敛的静态设置上使用 [1]：

```ini
LEPSILON = .TRUE.
LRPA     = .FALSE.
NSW      = 0
```

这里不规定通用 ENCUT、k 网格或真空厚度；它们必须分别收敛。中间滑移构型不能不加约束地重新弛豫回稳定堆叠态。依次检查：实际对称性及位点映射、单原子禁戒分量、等价原子协变关系、ASR 残差，以及模式投影。

二维超胞的 $P^{3D}$ 会被真空体积稀释；使用 $p_z=\Omega P_z^{3D}$ 或 $P_z^{2D}=L_zP_z^{3D}$ 统一比较。固定晶胞并保持相同电边界条件时，BEC 无需人为再除以真空厚度；应检查它随真空及偶极处理是否收敛。

独立的数值检查包括：

$$
\zeta_a\simeq\frac{p_z(s+\delta\hat a)-p_z(s-\delta\hat a)}{2e\delta},
\qquad
Z^*_{\kappa,zj}\simeq\frac{F_{\kappa j}(+E_z)-F_{\kappa j}(-E_z)}{2eE_z}.
$$

两者均须固定比较对象的原子坐标/路径约束，检查差分步长收敛和连续极化分支。这里没有开展任何新的 DFT 计算。

## 10. Python 与命令行使用

本模块是增量文件，不修改既有 `group-ops` 命令。完整仓库中可直接导入：

```python
import numpy as np
from group_theory_operations.born_effective_charges import (
    point_group_born_space, atomic_born_space, sliding_born_space,
    structure_born_space, read_outcar_born, bilayer_mode_charge,
    acoustic_sum_residual,
)

site = point_group_born_space("C3v")  # SITE symmetry, not a blind crystal-group substitution.
assert site.dimension == 2

c2x = np.diag([1., -1., -1.])
sliding = sliding_born_space([c2x])  # isolated xy bilayer; one allowed direction.
assert sliding.dimension == 1

z = read_outcar_born("OUTCAR")
labels = [1, 1, -1, -1]            # example for four atoms, in the actual OUTCAR order.
c = bilayer_mode_charge(z, labels)
print("Ez/sliding charge:", c[2])
print("ASR residual:", acoustic_sum_residual(z))
```

```bash
python -m group_theory_operations.born_effective_charges atomic --point-group C3v
python -m group_theory_operations.born_effective_charges sliding --operations examples/born/c2x.json
python -m group_theory_operations.born_effective_charges project --outcar OUTCAR --layers examples/born/four_atom_layers.json
python -m unittest discover -s tests -p 'test_born_effective_charges.py' -v
python scripts/verify_born_exact.py
python scripts/check_born_registry.py
```

字符串 `"3"` 按 Hermann–Mauguin 点群 $3$ 解释；整数 `3` 按点群编号 3，即 $C_2$ 解释。轴向以实际操作矩阵为准。输出的 `allowed_mask` 不是数值预测，基底还编码了分量相等、相反等关系，不能只保留布尔零/非零掩码。

## 参考文献与可核查出处

[1] VASP Wiki, **Born effective charges**, definition, raw OUTCAR/py4vasp index conventions, DFPT and finite-field routes. https://vasp.at/wiki/index.php/Born_effective_charges （访问日期 2026-11-09）

[2] Octopus documentation, **Born effective charges**, nonsymmetric atomic tensors, site symmetry, acoustic sum rule. https://octopus-code.desy.de/documentation/main/manual/calculations/born_effective_charges/ （访问日期 2026-11-09）

[3] C. Ke, F. Liu, S. Liu, **Superlubric Motion of Wavelike Domain Walls in Sliding Ferroelectrics**, Phys. Rev. Lett. **135**, 046201 (2025). DOI: https://doi.org/10.1103/jhlq-2dd7 . 核对正文 pp. 046201-2–3、Fig. 1–2；采用 field-first 的 $Z_{31},Z_{32}$。

[4] Z. Wang, S. Dong, **Polarization switching in sliding ferroelectrics: Roles of fluctuation and domain wall**, Phys. Rev. B **111**, L201406 (2025). DOI: https://doi.org/10.1103/PhysRevB.111.L201406 . 核对原版 PDF p. L201406-2 的 Eqs. (1)–(2)、Fig. 1，以及 p. L201406-3 的 Fig. 2；该文采用 force-first 次序，$Z^{\rm paper}_{xz}=Z^{\rm OUTCAR}_{zx}$。原图的导数为 $\partial P_j/\partial u_i$，不可依据错误的文本提取结果判其公式有误。

[5] P. Fachin, F. Macheda, P. Barone, F. Mauri, **Infrared markers of topological phase transitions in quantum spin Hall insulators**, npj Computational Materials **11**, 307 (2025). DOI: https://doi.org/10.1038/s41524-025-01780-6 . 核对 “Germanene”、 “Symmetries of jacutingaite atoms” 与 Eq. (3)：$3m$ 位点的两参数形式；$2/m$ 和 $m$ 位点的五参数形式；反演关系与电中性导致的额外消零；质量加权模式投影。

32 点群全表、式 (8)–(9) 与双层不动子空间的穷尽分类，是基于上述定义独立推导并由代码验证的结果，不声称这些文献逐项发表了本表。
