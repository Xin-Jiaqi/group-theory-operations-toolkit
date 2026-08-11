# 二次场空间的诱导表示

本文档由 `data/group_operations.json` 确定性生成。机器读取请使用 `data/quadratic_field_representations.json`；不要解析下列表格。

## 约定与推导

令 $Q^{bc}=E^bE^{c*}$。在对称基底

$$\boldsymbol{\mathcal E}_+=(|E^x|^2,|E^y|^2,|E^z|^2,E^xE^{y*}+E^yE^{x*},E^xE^{z*}+E^zE^{x*},E^yE^{z*}+E^zE^{y*})^{\mathsf T}$$

上，$\boldsymbol{\mathcal E}'_+=M_+(R)\boldsymbol{\mathcal E}_+$，且 $M_+(R)=\mathrm{Sym}^2D(R)$。

反对称子空间与实轴矢量 $\mathbf h=i\mathbf E\times\mathbf E^*$ 等价。由

$$[D\mathbf u]\times[D\mathbf v]=\det(D)D(\mathbf u\times\mathbf v)$$

得到 $\mathbf h'=M_-(R)\mathbf h$，其中

$$\boxed{M_-(R)=\det[D(R)]D(R)}.$$

因此 $M_-(R_1R_2)=M_-(R_1)M_-(R_2)$、$M_-^{-1}=M_-^{\mathsf T}$。反演满足 $D(I)=-I_3$，但 $M_-(I)=I_3$。所有 $M_-(R)$ 均为行列式 $+1$ 的正交矩阵。

## 全部操作的 $M_-(R)$

### `cubic_Oh`：48 个操作

| Index | $R$ | $\det D(R)$ | $M_-(R)$ |
|---:|---|---:|---|
| 0 | `1` | +1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 1 | `2_001` | +1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 2 | `2_010` | +1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 3 | `2_100` | +1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 4 | `3+_111` | +1 | $\begin{bmatrix}0&0&1\\1&0&0\\0&1&0\end{bmatrix}$ |
| 5 | `3+_-11-1` | +1 | $\begin{bmatrix}0&0&1\\-1&0&0\\0&-1&0\end{bmatrix}$ |
| 6 | `3+_1-1-1` | +1 | $\begin{bmatrix}0&0&-1\\-1&0&0\\0&1&0\end{bmatrix}$ |
| 7 | `3+_-1-11` | +1 | $\begin{bmatrix}0&0&-1\\1&0&0\\0&-1&0\end{bmatrix}$ |
| 8 | `3-_111` | +1 | $\begin{bmatrix}0&1&0\\0&0&1\\1&0&0\end{bmatrix}$ |
| 9 | `3-_1-1-1` | +1 | $\begin{bmatrix}0&-1&0\\0&0&1\\-1&0&0\end{bmatrix}$ |
| 10 | `3-_-1-11` | +1 | $\begin{bmatrix}0&1&0\\0&0&-1\\-1&0&0\end{bmatrix}$ |
| 11 | `3-_-11-1` | +1 | $\begin{bmatrix}0&-1&0\\0&0&-1\\1&0&0\end{bmatrix}$ |
| 12 | `2_110` | +1 | $\begin{bmatrix}0&1&0\\1&0&0\\0&0&-1\end{bmatrix}$ |
| 13 | `2_1-10` | +1 | $\begin{bmatrix}0&-1&0\\-1&0&0\\0&0&-1\end{bmatrix}$ |
| 14 | `4-_001` | +1 | $\begin{bmatrix}0&1&0\\-1&0&0\\0&0&1\end{bmatrix}$ |
| 15 | `4+_001` | +1 | $\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}$ |
| 16 | `4-_100` | +1 | $\begin{bmatrix}1&0&0\\0&0&1\\0&-1&0\end{bmatrix}$ |
| 17 | `2_011` | +1 | $\begin{bmatrix}-1&0&0\\0&0&1\\0&1&0\end{bmatrix}$ |
| 18 | `2_01-1` | +1 | $\begin{bmatrix}-1&0&0\\0&0&-1\\0&-1&0\end{bmatrix}$ |
| 19 | `4+_100` | +1 | $\begin{bmatrix}1&0&0\\0&0&-1\\0&1&0\end{bmatrix}$ |
| 20 | `4+_010` | +1 | $\begin{bmatrix}0&0&1\\0&1&0\\-1&0&0\end{bmatrix}$ |
| 21 | `2_101` | +1 | $\begin{bmatrix}0&0&1\\0&-1&0\\1&0&0\end{bmatrix}$ |
| 22 | `4-_010` | +1 | $\begin{bmatrix}0&0&-1\\0&1&0\\1&0&0\end{bmatrix}$ |
| 23 | `2_-101` | +1 | $\begin{bmatrix}0&0&-1\\0&-1&0\\-1&0&0\end{bmatrix}$ |
| 24 | `-1` | -1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 25 | `m_001` | -1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 26 | `m_010` | -1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 27 | `m_100` | -1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 28 | `-3+_111` | -1 | $\begin{bmatrix}0&0&1\\1&0&0\\0&1&0\end{bmatrix}$ |
| 29 | `-3+_-11-1` | -1 | $\begin{bmatrix}0&0&1\\-1&0&0\\0&-1&0\end{bmatrix}$ |
| 30 | `-3+_1-1-1` | -1 | $\begin{bmatrix}0&0&-1\\-1&0&0\\0&1&0\end{bmatrix}$ |
| 31 | `-3+_-1-11` | -1 | $\begin{bmatrix}0&0&-1\\1&0&0\\0&-1&0\end{bmatrix}$ |
| 32 | `-3-_111` | -1 | $\begin{bmatrix}0&1&0\\0&0&1\\1&0&0\end{bmatrix}$ |
| 33 | `-3-_1-1-1` | -1 | $\begin{bmatrix}0&-1&0\\0&0&1\\-1&0&0\end{bmatrix}$ |
| 34 | `-3-_-1-11` | -1 | $\begin{bmatrix}0&1&0\\0&0&-1\\-1&0&0\end{bmatrix}$ |
| 35 | `-3-_-11-1` | -1 | $\begin{bmatrix}0&-1&0\\0&0&-1\\1&0&0\end{bmatrix}$ |
| 36 | `m_110` | -1 | $\begin{bmatrix}0&1&0\\1&0&0\\0&0&-1\end{bmatrix}$ |
| 37 | `m_1-10` | -1 | $\begin{bmatrix}0&-1&0\\-1&0&0\\0&0&-1\end{bmatrix}$ |
| 38 | `-4-_001` | -1 | $\begin{bmatrix}0&1&0\\-1&0&0\\0&0&1\end{bmatrix}$ |
| 39 | `-4+_001` | -1 | $\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}$ |
| 40 | `-4-_100` | -1 | $\begin{bmatrix}1&0&0\\0&0&1\\0&-1&0\end{bmatrix}$ |
| 41 | `m_011` | -1 | $\begin{bmatrix}-1&0&0\\0&0&1\\0&1&0\end{bmatrix}$ |
| 42 | `m_01-1` | -1 | $\begin{bmatrix}-1&0&0\\0&0&-1\\0&-1&0\end{bmatrix}$ |
| 43 | `-4+_100` | -1 | $\begin{bmatrix}1&0&0\\0&0&-1\\0&1&0\end{bmatrix}$ |
| 44 | `-4+_010` | -1 | $\begin{bmatrix}0&0&1\\0&1&0\\-1&0&0\end{bmatrix}$ |
| 45 | `m_101` | -1 | $\begin{bmatrix}0&0&1\\0&-1&0\\1&0&0\end{bmatrix}$ |
| 46 | `-4-_010` | -1 | $\begin{bmatrix}0&0&-1\\0&1&0\\1&0&0\end{bmatrix}$ |
| 47 | `m_-101` | -1 | $\begin{bmatrix}0&0&-1\\0&-1&0\\-1&0&0\end{bmatrix}$ |

### `tetragonal_D4h`：16 个操作

| Index | $R$ | $\det D(R)$ | $M_-(R)$ |
|---:|---|---:|---|
| 0 | `1` | +1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 1 | `2_001` | +1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 2 | `4+_001` | +1 | $\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}$ |
| 3 | `4-_001` | +1 | $\begin{bmatrix}0&1&0\\-1&0&0\\0&0&1\end{bmatrix}$ |
| 4 | `2_010` | +1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 5 | `2_100` | +1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 6 | `2_110` | +1 | $\begin{bmatrix}0&1&0\\1&0&0\\0&0&-1\end{bmatrix}$ |
| 7 | `2_1-10` | +1 | $\begin{bmatrix}0&-1&0\\-1&0&0\\0&0&-1\end{bmatrix}$ |
| 8 | `-1` | -1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 9 | `m_001` | -1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 10 | `-4+_001` | -1 | $\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}$ |
| 11 | `-4-_001` | -1 | $\begin{bmatrix}0&1&0\\-1&0&0\\0&0&1\end{bmatrix}$ |
| 12 | `m_010` | -1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 13 | `m_100` | -1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 14 | `m_110` | -1 | $\begin{bmatrix}0&1&0\\1&0&0\\0&0&-1\end{bmatrix}$ |
| 15 | `m_1-10` | -1 | $\begin{bmatrix}0&-1&0\\-1&0&0\\0&0&-1\end{bmatrix}$ |

### `hexagonal_D6h`：24 个操作

| Index | $R$ | $\det D(R)$ | $M_-(R)$ |
|---:|---|---:|---|
| 0 | `1` | +1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 1 | `3+_001` | +1 | $\begin{bmatrix}-\frac12&-\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&-\frac12&0\\0&0&1\end{bmatrix}$ |
| 2 | `3-_001` | +1 | $\begin{bmatrix}-\frac12&\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&-\frac12&0\\0&0&1\end{bmatrix}$ |
| 3 | `2_001` | +1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 4 | `6-_001` | +1 | $\begin{bmatrix}\frac12&\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&\frac12&0\\0&0&1\end{bmatrix}$ |
| 5 | `6+_001` | +1 | $\begin{bmatrix}\frac12&-\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&\frac12&0\\0&0&1\end{bmatrix}$ |
| 6 | `2_110` | +1 | $\begin{bmatrix}-\frac12&\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&\frac12&0\\0&0&-1\end{bmatrix}$ |
| 7 | `2_100` | +1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 8 | `2_010` | +1 | $\begin{bmatrix}-\frac12&-\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&\frac12&0\\0&0&-1\end{bmatrix}$ |
| 9 | `2_1-10` | +1 | $\begin{bmatrix}\frac12&-\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&-\frac12&0\\0&0&-1\end{bmatrix}$ |
| 10 | `2_120` | +1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 11 | `2_210` | +1 | $\begin{bmatrix}\frac12&\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&-\frac12&0\\0&0&-1\end{bmatrix}$ |
| 12 | `-1` | -1 | $\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$ |
| 13 | `-3+_001` | -1 | $\begin{bmatrix}-\frac12&-\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&-\frac12&0\\0&0&1\end{bmatrix}$ |
| 14 | `-3-_001` | -1 | $\begin{bmatrix}-\frac12&\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&-\frac12&0\\0&0&1\end{bmatrix}$ |
| 15 | `m_001` | -1 | $\begin{bmatrix}-1&0&0\\0&-1&0\\0&0&1\end{bmatrix}$ |
| 16 | `-6-_001` | -1 | $\begin{bmatrix}\frac12&\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&\frac12&0\\0&0&1\end{bmatrix}$ |
| 17 | `-6+_001` | -1 | $\begin{bmatrix}\frac12&-\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&\frac12&0\\0&0&1\end{bmatrix}$ |
| 18 | `m_110` | -1 | $\begin{bmatrix}-\frac12&\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&\frac12&0\\0&0&-1\end{bmatrix}$ |
| 19 | `m_100` | -1 | $\begin{bmatrix}1&0&0\\0&-1&0\\0&0&-1\end{bmatrix}$ |
| 20 | `m_010` | -1 | $\begin{bmatrix}-\frac12&-\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&\frac12&0\\0&0&-1\end{bmatrix}$ |
| 21 | `m_1-10` | -1 | $\begin{bmatrix}\frac12&-\frac{\sqrt3}{2}&0\\-\frac{\sqrt3}{2}&-\frac12&0\\0&0&-1\end{bmatrix}$ |
| 22 | `m_120` | -1 | $\begin{bmatrix}-1&0&0\\0&1&0\\0&0&-1\end{bmatrix}$ |
| 23 | `m_210` | -1 | $\begin{bmatrix}\frac12&\frac{\sqrt3}{2}&0\\\frac{\sqrt3}{2}&-\frac12&0\\0&0&-1\end{bmatrix}$ |

## 机器验证

测试逐操作验证定义式、正交性、行列式、叉积恒等式、$M_+$/$M_-$ 同态，并核对本文件与生成 JSON。88 个母群操作均覆盖；同名操作在不同母群中保留各自记录。

重新生成：`python scripts/generate_quadratic_field_representations.py`。
