# 实验报告：线性方程组的直接法与迭代法求解及误差分析

## 第一部分：直接法

直接法是通过有限步算术运算（不考虑舍入误差）得到线性方程组 $Ax = b$ 精确解的一类方法。

---

### 一、Doolittle 分解法

#### 1. 基本原理

Doolittle 分解将非奇异方阵 $A$ 分解为一个**单位下三角矩阵** $L$（对角线元素全为1）和一个**上三角矩阵** $U$ 的乘积
$A = LU$

---
$$
\begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nn}
\end{pmatrix}=\begin{pmatrix}
1 & 0 & \cdots & 0 \\
l_{21} & 1 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
l_{n1} & l_{n2} & \cdots & 1
\end{pmatrix}
\begin{pmatrix}
u_{11} & u_{12} & \cdots & u_{1n} \\
0 & u_{22} & \cdots & u_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & u_{nn}
\end{pmatrix}
$$

#### 2. 分解公式（按行计算 $U$，按列计算 $L$）

对于 $k = 1, 2, \ldots, n$：

**计算 $U$ 的第 $k$ 行元素**（$j = k, k+1, \ldots, n$）：

$u_{kj} = a_{kj} - \sum_{r=1}^{k-1} l_{kr} u_{rj}$

**计算 $L$ 的第 $k$ 列元素**（$i = k+1, k+2, \ldots, n$）：

$l_{ik} = \frac{1}{u_{kk}} \left( a_{ik} - \sum_{r=1}^{k-1} l_{ir} u_{rk} \right)$

#### 3. 求解过程

原方程组 $Ax = b$ 转化为：

$LUx = b$

**前代**：令 $y = Ux$，解 $Ly = b$

$
\begin{cases}
y_1 = b_1 \\
y_i = b_i - \sum_{j=1}^{i-1} l_{ij} y_j, \quad i = 2, 3, \ldots, n
\end{cases}
$

**回代**：解 $Ux = y$

$
\begin{cases}
x_n = y_n / u_{nn} \\
x_i = \left( y_i - \sum_{j=i+1}^n u_{ij} x_j \right) / u_{ii}, \quad i = n-1, n-2, \ldots, 1
\end{cases}
$

---

### 二、Crout 分解法

#### 1. 基本原理

Crout 分解将矩阵 $A$ 分解为一个**下三角矩阵** $L$ 和一个**单位上三角矩阵** $U$（对角线元素全为1）的乘积：

$A = LU$

其中：

$
L = \begin{pmatrix}
l_{11} & 0 & \cdots & 0 \\
l_{21} & l_{22} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
l_{n1} & l_{n2} & \cdots & l_{nn}
\end{pmatrix}, \quad
U = \begin{pmatrix}
1 & u_{12} & \cdots & u_{1n} \\
0 & 1 & \cdots & u_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{pmatrix}
$

#### 2. 分解公式（按列计算 $L$ 和 $U$）

对于 $k = 1, 2, \ldots, n$：

**计算 $L$ 的第 $k$ 列元素**（$i = k, k+1, \ldots, n$）：

$l_{ik} = a_{ik} - \sum_{r=1}^{k-1} l_{ir} u_{rk}$

**计算 $U$ 的第 $k$ 行元素**（$j = k+1, k+2, \ldots, n$）：

$u_{kj} = \frac{1}{l_{kk}} \left( a_{kj} - \sum_{r=1}^{k-1} l_{kr} u_{rj} \right)$

#### 3. 求解过程

与 Doolittle 相同：先前代解 $Ly = b$，再回代解 $Ux = y$。

---

### 三、对称正定矩阵的 Cholesky 分解法

#### 1. 适用条件

矩阵 $A$ 必须满足：
- **对称性**：$A^T = A$
- **正定性**：对所有非零向量 $x$，有 $x^T A x > 0$

#### 2. 基本原理

Cholesky 分解将对称正定矩阵 $A$ 分解为：

$A = LL^T$

其中 $L$ 是**下三角矩阵**，且对角线元素 $l_{ii} > 0$。

#### 3. 分解公式

对于 $j = 1, 2, \ldots, n$：

**计算对角元**：

$l_{jj} = \sqrt{a_{jj} - \sum_{k=1}^{j-1} l_{jk}^2}$

**计算非对角元**（$i = j+1, j+2, \ldots, n$）：

$l_{ij} = \frac{1}{l_{jj}} \left( a_{ij} - \sum_{k=1}^{j-1} l_{ik} l_{jk} \right)$

#### 4. 求解过程

原方程 $Ax = b$ 转化为 $LL^T x = b$。

**前代**：解 $Ly = b$（$y = L^T x$）

$
y_i = \frac{1}{l_{ii}} \left( b_i - \sum_{j=1}^{i-1} l_{ij} y_j \right), \quad i = 1, 2, \ldots, n
$

**回代**：解 $L^T x = y$

$
x_i = \frac{1}{l_{ii}} \left( y_i - \sum_{j=i+1}^n l_{ji} x_j \right), \quad i = n, n-1, \ldots, 1
$

#### 5. 特点

- 计算量约为 Doolittle 的一半（只需计算 $L$）
- 数值稳定，不需要选主元
- 节省存储空间（只需存储 $L$）

---

### 四、三对角矩阵的追赶法

#### 1. 矩阵形式

三对角矩阵形式如下：

$
A = \begin{pmatrix}
b_1 & c_1 & & & 0 \\
a_2 & b_2 & c_2 & & \\
& a_3 & b_3 & \ddots & \\
& & \ddots & \ddots & c_{n-1} \\
0 & & & a_n & b_n
\end{pmatrix}
$

其中 $a_i, b_i, c_i$ 为已知系数，且满足 $a_i \neq 0$（通常可保证）。

#### 2. 基本原理

追赶法是 Doolittle 分解在三对角矩阵上的特化形式。分解 $A = LU$，其中：

$
L = \begin{pmatrix}
1 & & & & 0 \\
l_2 & 1 & & & \\
& l_3 & 1 & & \\
& & \ddots & \ddots & \\
0 & & & l_n & 1
\end{pmatrix}, \quad
U = \begin{pmatrix}
u_1 & c_1 & & & 0 \\
& u_2 & c_2 & & \\
& & u_3 & \ddots & \\
& & & \ddots & c_{n-1} \\
0 & & & & u_n
\end{pmatrix}
$

#### 3. 追的过程（计算 $l_i, u_i$）

$
\begin{cases}
u_1 = b_1 \\
l_i = \dfrac{a_i}{u_{i-1}}, \quad i = 2, 3, \ldots, n \\
u_i = b_i - l_i c_{i-1}, \quad i = 2, 3, \ldots, n
\end{cases}
$

#### 4. 赶的过程（求解）

**前代**：解 $Ly = b$

$
\begin{cases}
y_1 = b_1 \\
y_i = b_i - l_i y_{i-1}, \quad i = 2, 3, \ldots, n
\end{cases}
$

**回代**：解 $Ux = y$

$
\begin{cases}
x_n = y_n / u_n \\
x_i = (y_i - c_i x_{i+1}) / u_i, \quad i = n-1, n-2, \ldots, 1
\end{cases}
$

#### 5. 特点

- 计算量仅 $O(n)$，非常高效
- 只需存储 4 个一维数组（$a_i, b_i, c_i, b_i$）
- 广泛用于求解常微分方程边值问题和偏微分方程差分格式

---

### 五、列主元消去法

#### 1. 基本原理

在 Gauss 消去过程中，每次选取当前列中**绝对值最大的元素**作为主元，并通过**行交换**将其移到主元位置，然后进行消元。

#### 2. 详细步骤

对于 $k = 1, 2, \ldots, n-1$：

**第 1 步：选主元**
在 $A(k:n, k)$ 中寻找绝对值最大的元素：

$|a_{p,k}| = \max_{k \le i \le n} |a_{i,k}|$

**第 2 步：判断**
若 $|a_{p,k}| < \varepsilon$（$\varepsilon$ 为容差），则矩阵奇异，停止计算。

**第 3 步：行交换**
若 $p \neq k$，交换第 $k$ 行与第 $p$ 行（同时交换右端项 $b$ 的对应分量）。

**第 4 步：消元**
对于 $i = k+1, k+2, \ldots, n$：

$l_{ik} = \frac{a_{ik}}{a_{kk}} \quad$（消元因子）

对于 $j = k+1, k+2, \ldots, n$：

$a_{ij} = a_{ij} - l_{ik} a_{kj}$

$b_i = b_i - l_{ik} b_k$

#### 3. 回代求解

消元完成后得到上三角方程组 $Ux = \tilde{b}$，回代求解：

$
\begin{cases}
x_n = \tilde{b}_n / u_{nn} \\
x_i = \left( \tilde{b}_i - \sum_{j=i+1}^n u_{ij} x_j \right) / u_{ii}, \quad i = n-1, n-2, \ldots, 1
\end{cases}
$

#### 4. 特点

- 有效避免小主元导致的数值不稳定
- 计算量约为 $O(n^3/3)$
- 实际计算中最常用的直接法之一

---

## 第二部分：迭代法

迭代法从一个初始向量 $x^{(0)}$ 出发，通过递推格式生成序列 $\{x^{(k)}\}$，使其收敛于真解。

---

### 一、Jacobi 迭代法

#### 1. 基本原理

将系数矩阵 $A$ 分解为：

$A = D - L - U$

其中：
- $D$：对角矩阵，$D = \text{diag}(a_{11}, a_{22}, \ldots, a_{nn})$
- $-L$：严格下三角部分（$L$ 是下三角且对角元为 0）
- $-U$：严格上三角部分（$U$ 是上三角且对角元为 0）

#### 2. 迭代格式

由 $Ax = b$ 得 $(D - L - U)x = b$，整理为 $Dx = (L+U)x + b$，则：

$x^{(k+1)} = D^{-1}(L+U)x^{(k)} + D^{-1}b$

分量形式（更常用）：

$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^n a_{ij} x_j^{(k)} \right), \quad i = 1, 2, \ldots, n$

#### 3. 计算步骤

1. 选取初始向量 $x^{(0)} = (x_1^{(0)}, x_2^{(0)}, \ldots, x_n^{(0)})^T$
2. 对于 $k = 0, 1, 2, \ldots$，依次计算 $x_1^{(k+1)}, x_2^{(k+1)}, \ldots, x_n^{(k+1)}$
3. 检查收敛性：若 $\|x^{(k+1)} - x^{(k)}\| < \varepsilon$，停止

#### 4. 收敛条件

充分条件（满足其一即可）：
- $A$ 严格对角占优：$|a_{ii}| > \sum_{j \neq i} |a_{ij}|, \quad \forall i$
- $A$ 不可约且弱对角占优

#### 5. 特点

- 每步计算独立，**天然并行**
- 收敛速度通常较慢（线性收敛）
- 存储量小，只需保留 $x^{(k)}$ 和 $x^{(k+1)}$

---

### 二、Gauss-Seidel 迭代法

#### 1. 基本原理

与 Jacobi 不同，Gauss-Seidel 在计算 $x_i^{(k+1)}$ 时，**立即使用**已经计算出的最新分量。

#### 2. 迭代格式

矩阵形式：

$x^{(k+1)} = (D-L)^{-1}U x^{(k)} + (D-L)^{-1}b$

分量形式：

$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^n a_{ij} x_j^{(k)} \right)$

#### 3. 计算步骤

1. 选取初始向量 $x^{(0)}$
2. 对于 $k = 0, 1, 2, \ldots$：
   - 计算 $x_1^{(k+1)}$（使用 $x_2^{(k)}, x_3^{(k)}, \ldots$）
   - 计算 $x_2^{(k+1)}$（使用 $x_1^{(k+1)}, x_3^{(k)}, \ldots$）
   - 计算 $x_3^{(k+1)}$（使用 $x_1^{(k+1)}, x_2^{(k+1)}, x_4^{(k)}, \ldots$）
   - 依次类推
3. 检查收敛性

#### 4. 收敛条件

与 Jacobi 类似，且：
- 若 Jacobi 收敛，Gauss-Seidel 收敛更快（通常）
- 对某些矩阵，Gauss-Seidel 收敛但 Jacobi 不收敛

#### 5. 特点

- 收敛速度一般快于 Jacobi
- 只需存储一组向量（可原地更新）
- **天然串行**，无法直接并行

---

### 三、SOR（逐次超松弛）迭代法

#### 1. 基本原理

SOR 在 Gauss-Seidel 的基础上引入**松弛因子** $\omega$，将当前修正量放大或缩小，以加速收敛。

#### 2. 迭代格式

首先计算 Gauss-Seidel 的更新值：

$\tilde{x}_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^n a_{ij} x_j^{(k)} \right)$

然后进行松弛：

$x_i^{(k+1)} = (1-\omega)x_i^{(k)} + \omega \tilde{x}_i^{(k+1)}$

合并为一个公式：

$x_i^{(k+1)} = (1-\omega)x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^n a_{ij} x_j^{(k)} \right)$

#### 3. 松弛因子的选择

- $0 < \omega < 1$：低松弛（用于某些难以收敛的问题）
- $\omega = 1$：退化为 Gauss-Seidel
- $1 < \omega < 2$：超松弛（**最常用**，加速收敛）
- $\omega \ge 2$：一般不收敛

**最优松弛因子**（对对称正定矩阵）：

$\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho(B_{\text{Jacobi}})^2}}$

其中 $\rho(B_{\text{Jacobi}})$ 是 Jacobi 迭代矩阵的谱半径。

#### 4. 计算步骤

1. 选取初始向量 $x^{(0)}$ 和松弛因子 $\omega$（通常取 1.2～1.6）
2. 按分量形式迭代更新
3. 检查收敛性

#### 5. 特点

- 适当选择 $\omega$ 可**显著加速**收敛（比 Gauss-Seidel 快数倍）
- 最优 $\omega$ 的估计需要了解矩阵性质
- 最常用的经典迭代法之一

---

## 第三部分：误差分析

### 一、误差范数的定义

设：
- $x^*$：精确解（或足够精确的参考解）
- $x$：数值方法求得的近似解

**误差向量**：

$e = x - x^*$

### 二、常用向量范数

#### 1. 1-范数（列和范数）

$\|e\|_1 = \sum_{i=1}^n |e_i|$

#### 2. 2-范数（欧几里得范数）

$\|e\|_2 = \sqrt{\sum_{i=1}^n e_i^2}$

#### 3. 无穷范数（最大范数，最常用）

$\|e\|_\infty = \max_{1 \le i \le n} |e_i|$

### 三、误差分析步骤

#### 第 1 步：确定参考解 $x^*$

- **情况 A**：方程组有理论解（如人工构造的例题），直接使用
- **情况 B**：理论解未知，使用高精度直接法（如列主元消去法或 MATLAB 的 `A\b`）计算一个参考解

#### 第 2 步：计算各方法的近似解

分别用 Doolittle、Crout、Cholesky、追赶法、列主元消去法、Jacobi、Gauss-Seidel、SOR 等方法求解，得到 $x_{\text{method}}$。

#### 第 3 步：计算误差向量

$e_{\text{method}} = x_{\text{method}} - x^*$

#### 第 4 步：计算各范数

例如，对于 $e = (e_1, e_2, \ldots, e_n)^T$：

$\|e\|_2 = \sqrt{e_1^2 + e_2^2 + \cdots + e_n^2}$

$\|e\|_\infty = \max(|e_1|, |e_2|, \ldots, |e_n|)$

#### 第 5 步：结果分析

- 范数越小 → 近似解越准确
- 比较不同方法的范数 → 判断哪种方法在该问题上更优
- 对于迭代法，记录迭代次数与误差范数的关系

### 四、数值示例

**设**：

$x^* = \begin{pmatrix} 1 \\ -1 \\ 2 \end{pmatrix}, \quad
x = \begin{pmatrix} 1.01 \\ -0.98 \\ 1.99 \end{pmatrix}$

**误差向量**：

$e = \begin{pmatrix} 0.01 \\ 0.02 \\ -0.01 \end{pmatrix}$

**计算**：

$\|e\|_2 = \sqrt{0.01^2 + 0.02^2 + (-0.01)^2} = \sqrt{0.0001 + 0.0004 + 0.0001} = \sqrt{0.0006} \approx 0.02449$

$\|e\|_\infty = \max(0.01, 0.02, 0.01) = 0.02$

### 五、实验报告结果表格模板

| 方法 | $\|e\|_2$ | $\|e\|_\infty$ | 迭代次数（如适用） |
|------|-------------|------------------|-------------------|
| Doolittle 分解 | 1.23e-14 | 2.45e-14 | — |
| Crout 分解 | 1.25e-14 | 2.50e-14 | — |
| Cholesky 分解 | 1.20e-14 | 2.40e-14 | — |
| 追赶法 | 1.22e-14 | 2.44e-14 | — |
| 列主元消去法 | 1.18e-14 | 2.38e-14 | — |
| Jacobi 迭代 | 3.45e-3 | 5.67e-3 | 100 |
| Gauss-Seidel 迭代 | 1.23e-3 | 2.34e-3 | 68 |
| SOR（$\omega = 1.2$） | 2.31e-4 | 4.12e-4 | 32 |

### 六、结论分析

1. **直接法**：精度高（误差通常在 $10^{-14}$ 量级，仅受机器精度限制），但计算量大（$O(n^3)$），适合中小规模问题。

2. **迭代法**：精度可控（取决于迭代次数），计算量小（每步 $O(n^2)$），适合大规模稀疏矩阵。

3. **收敛速度**：SOR > Gauss-Seidel > Jacobi（当参数选择适当时）。

4. **实际选择**：
   - 小规模稠密矩阵 → 列主元消去法或 Doolittle
   - 对称正定矩阵 → Cholesky
   - 三对角矩阵 → 追赶法
   - 大规模稀疏矩阵 → SOR 或共轭梯度法（CG）