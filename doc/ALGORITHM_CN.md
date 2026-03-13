# PV-LIO-PLUS 算法详解

## 1. 算法概述

PV-LIO-PLUS 的核心是一个**紧耦合LiDAR-惯性里程计**系统。系统融合了激光雷达点云和IMU惯性测量数据，通过**迭代扩展卡尔曼滤波器（iESKF）**在**流形（Manifold）**上进行状态估计。地图使用**自适应概率体素地图（VoxelMap / VoxelMap++）**进行管理，通过**点到平面（Point-to-Plane）**匹配提供观测约束。

### 1.1 算法核心流程

```
传感器数据 → 预处理 → IMU传播 → 点云去畸变 → 降采样 → 迭代EKF更新 → 地图更新 → 输出
```

---

## 2. 状态空间与流形表示

### 2.1 状态定义

系统状态定义在**流形（Manifold）**上，总共包含 **23 维**的状态空间（误差状态为 23 维）：

```
x = [p, R, R_LI, t_LI, v, bg, ba, g]^T
```

| 状态变量 | 符号 | 维度 | 流形 | 说明 |
|----------|------|------|------|------|
| 位置 | p | 3 | ℝ³ | 机体在世界坐标系下的位置 |
| 旋转 | R | 3 | SO(3) | 机体在世界坐标系下的旋转 |
| LiDAR-IMU旋转外参 | R_LI | 3 | SO(3) | LiDAR到IMU的旋转外参 |
| LiDAR-IMU平移外参 | t_LI | 3 | ℝ³ | LiDAR到IMU的平移外参 |
| 速度 | v | 3 | ℝ³ | 机体在世界坐标系下的速度 |
| 陀螺仪偏置 | bg | 3 | ℝ³ | IMU陀螺仪偏置 |
| 加速度计偏置 | ba | 3 | ℝ³ | IMU加速度计偏置 |
| 重力加速度 | g | 2 | S² | 重力向量（在S2球面流形上） |

> **为什么使用流形？** 旋转矩阵属于SO(3)群，不是普通的欧氏空间。直接在欧氏空间中对旋转进行加减运算会破坏其正交性约束。通过在流形上定义状态和运算，可以保证旋转的约束始终被满足。同样，重力向量的模长固定为 g ≈ 9.81 m/s²，使用 S² 流形表示可以确保这一约束。

### 2.2 流形运算

流形上定义了以下关键运算（⊞ 和 ⊟ 运算符）：

```
加法（⊞）: x' = x ⊞ δx    // 流形上的更新操作
减法（⊟）: δx = x' ⊟ x     // 计算流形上的差异

对于 SO(3) 部分:
    R' = R · Exp(δφ)         // 旋转更新
    δφ = Log(R^T · R')       // 旋转差异

对于 S² 部分:
    g' = g ⊞ δu             // S2球面上的更新
    δu = g ⊟ g'             // S2球面上的差异（2维）

对于 ℝ³ 部分:
    v' = v + δv              // 普通向量加法
```

### 2.3 输入和噪声定义

#### 系统输入（IMU测量）

```
u = [a_m, ω_m]^T

其中：
  a_m = a_true + ba + na     // 加速度计测量 = 真值 + 偏置 + 噪声
  ω_m = ω_true + bg + ng     // 陀螺仪测量 = 真值 + 偏置 + 噪声
```

#### 过程噪声

```
w = [ng, na, nbg, nba]^T    // 12维过程噪声

  ng  ~ N(0, σ²_g · I)      // 陀螺仪测量噪声
  na  ~ N(0, σ²_a · I)      // 加速度计测量噪声
  nbg ~ N(0, σ²_bg · I)     // 陀螺仪偏置随机游走噪声
  nba ~ N(0, σ²_ba · I)     // 加速度计偏置随机游走噪声
```

---

## 3. IMU预积分与运动模型

### 3.1 连续时间运动方程

IMU驱动的系统动力学方程如下：

```
ṗ = v                                              // 位置导数 = 速度
Ṙ = R · [ω_m - bg - ng]×                           // 旋转导数（×表示反对称矩阵）
v̇ = R · (a_m - ba - na) + g                        // 速度导数
ḃg = nbg                                           // 陀螺仪偏置随机游走
ḃa = nba                                           // 加速度计偏置随机游走
ġ = 0                                              // 重力不变
Ṙ_LI = 0                                           // 外参不变
ṫ_LI = 0                                           // 外参不变
```

### 3.2 离散时间传播

在离散时间步长 Δt 内，使用中值积分进行状态传播：

```
// 角速度中值
ω_mid = (ω_k + ω_{k+1}) / 2 - bg

// 旋转更新（Rodrigues公式）
R_{k+1} = R_k · Exp(ω_mid · Δt)

// 加速度中值（世界系下）
a_mid = (R_k · (a_k - ba) + R_{k+1} · (a_{k+1} - ba)) / 2 + g

// 位置和速度更新
p_{k+1} = p_k + v_k · Δt + 0.5 · a_mid · Δt²
v_{k+1} = v_k + a_mid · Δt
```

### 3.3 IMU初始化

系统启动时需要进行IMU初始化（`IMU_init()`），步骤如下：

1. **静止检测**：收集初始约200帧IMU数据
2. **均值估计**：计算加速度和陀螺仪的均值
   ```
   mean_acc = (1/N) · Σ a_i
   mean_gyr = (1/N) · Σ ω_i
   ```
3. **噪声估计**：计算加速度和陀螺仪的协方差
   ```
   cov_acc = (1/N) · Σ (a_i - mean_acc)(a_i - mean_acc)^T
   cov_gyr = (1/N) · Σ (ω_i - mean_gyr)(ω_i - mean_gyr)^T
   ```
4. **重力对齐**：根据加速度均值估计初始重力方向
   ```
   g_init = -mean_acc / ||mean_acc|| · G_m_s2
   ```
5. **状态初始化**：设置初始状态估计和协方差矩阵

### 3.4 点云运动去畸变

激光雷达扫描期间，机体在运动中采集点云，导致同一帧点云中的各点对应不同时刻的机体位姿。去畸变（`UndistortPcl()`）的目的是将所有点统一到扫描结束时刻的坐标系下。

**算法步骤**：

1. **前向传播**：利用IMU测量，从扫描开始时刻逐步传播到扫描结束时刻，记录每个IMU时刻的位姿
   ```
   对于每个IMU测量 (ω_k, a_k, dt_k):
     R_{k+1} = R_k · Exp(ω_k · dt_k)
     v_{k+1} = v_k + (R_k · a_k + g) · dt_k
     p_{k+1} = p_k + v_k · dt_k + 0.5 · (R_k · a_k + g) · dt_k²
   ```

2. **后向补偿**：对每个LiDAR点，根据其时间戳找到对应的IMU位姿，将其从采集时刻变换到扫描结束时刻
   ```
   对于每个LiDAR点 p_i（采集时刻 t_i）：
     // 找到 t_i 对应的IMU位姿 (R_i, p_i_imu)
     // 计算从 t_i 到 t_end 的相对变换
     R_rel = R_end^T · R_i
     t_rel = R_end^T · (p_i_imu - p_end)
     // 变换点到扫描结束时刻
     p_corrected = R_rel · (R_LI · p_i_lidar + t_LI) + t_rel + t_LI
   ```

---

## 4. 迭代扩展卡尔曼滤波器（iESKF）

### 4.1 算法概述

系统使用**迭代扩展卡尔曼滤波器（iterated Extended Kalman Filter on Manifold, iESKF）**进行状态估计。与标准EKF不同，iESKF在每次观测更新时进行多次迭代线性化，以提高非线性系统的估计精度。

### 4.2 预测步骤

利用IMU测量进行状态传播和协方差传播：

```
// 状态预测（使用运动模型）
x̂_{k|k-1} = f(x̂_{k-1|k-1}, u_k)

// 状态转移Jacobian
F_x = ∂f/∂x |_{x̂_{k-1|k-1}}    // 24×23 矩阵

// 噪声Jacobian
F_w = ∂f/∂w |_{x̂_{k-1|k-1}}    // 24×12 矩阵

// 协方差传播
P_{k|k-1} = F_x · P_{k-1|k-1} · F_x^T + F_w · Q · F_w^T
```

其中 Q 是过程噪声协方差矩阵（12×12对角矩阵）：

```
Q = diag(σ²_g, σ²_g, σ²_g, σ²_a, σ²_a, σ²_a,
         σ²_bg, σ²_bg, σ²_bg, σ²_ba, σ²_ba, σ²_ba)
```

### 4.3 迭代更新步骤

观测更新在每帧LiDAR数据到来时执行，使用迭代线性化提高精度：

```
初始化: x̂^(0) = x̂_{k|k-1}    // 使用预测值作为初始迭代点

for j = 0 to max_iteration - 1:

  1. 计算观测残差 h 和Jacobian H:
     // 在当前迭代估计 x̂^(j) 处线性化
     h_i = n_i^T · (R^(j) · (R_LI · p_i^L + t_LI) + p^(j) - q_i)   // 点到平面距离
     H_i = ∂h_i/∂δx |_{x̂^(j)}                                       // 观测Jacobian (1×23)

  2. 组装观测方程:
     H = [H_1; H_2; ...; H_m]    // m×23 Jacobian矩阵
     h = [h_1; h_2; ...; h_m]    // m×1 残差向量
     R_obs = diag(r_1, r_2, ..., r_m)  // m×m 观测噪声协方差

  3. 计算卡尔曼增益:
     // 使用对角观测噪声简化计算
     K = P_{k|k-1} · H^T · (H · P_{k|k-1} · H^T + R_obs)^{-1}

  4. 更新状态:
     δx = K · (h - H · (x̂^(j) ⊟ x̂_{k|k-1}))
     x̂^(j+1) = x̂^(j) ⊞ δx

  5. 收敛判断:
     if ||δx|| < ε:
       break

// 最终更新协方差
P_{k|k} = (I - K · H) · P_{k|k-1}
```

### 4.4 观测模型

系统提供两种观测模型，分别对应VoxelMap和VoxelMap++：

#### VoxelMap观测模型 (`observation_model_share`)

```
对于每个匹配的点-平面对 (p_i^L, π_i):

  // 将LiDAR系点变换到世界系
  p_i^W = R · (R_LI · p_i^L + t_LI) + p

  // 计算点到平面的距离（残差）
  h_i = n_i^T · p_i^W + d_i

  // 计算Jacobian矩阵 H_i (1×23)
  // 对位置 p 的偏导: n_i^T
  // 对旋转 R 的偏导: -n_i^T · R · [R_LI · p_i^L + t_LI]×
  // 对外参 R_LI 的偏导: -n_i^T · R · R_LI · [p_i^L]×
  // 对外参 t_LI 的偏导: n_i^T · R

  // 计算观测噪声（包含平面协方差和点协方差）
  r_i = n_i^T · J_p · Cov_point · J_p^T · n_i + n_i^T · Cov_plane · n_i
```

#### VoxelMap++观测模型 (`observation_model_share_plus`)

```
对于每个匹配的点-平面对 (p_i^L, π_i):

  // VoxelMap++使用不同的平面编码方式
  // 根据主方向（0/1/2）选择不同的残差计算
  // 主方向0: z = -(ax + by + d)
  // 主方向1: y = -(ax + bz + d)
  // 主方向2: x = -(ay + bz + d)

  // 残差计算（以主方向0为例）
  h_i = p_i^W[2] + omega[0]*p_i^W[0] + omega[1]*p_i^W[1] + omega[2]

  // Jacobian根据主方向不同而不同
  // 观测噪声同样考虑平面和点的协方差
```

---

## 5. 自适应概率体素地图

### 5.1 VoxelMap算法

#### 5.1.1 基本概念

VoxelMap将三维空间划分为规则的体素网格，每个体素内部使用**八叉树（OctoTree）**进行多层级平面拟合。

```
空间划分示意（2D示意，实际为3D）：

层级0（最粗）     层级1         层级2         层级3（最细）
┌─────────┐    ┌────┬────┐   ┌──┬──┬──┬──┐  每个子格进一步
│         │    │    │    │   │  │  │  │  │  划分为8个子格
│         │    │    │    │   │  │  │  │  │  （3D为8叉）
│         │    ├────┼────┤   ├──┼──┼──┼──┤
│         │    │    │    │   │  │  │  │  │
└─────────┘    └────┴────┘   └──┴──┴──┴──┘

体素大小: voxel_size=0.5m
层级数: max_layer=4（层0到层3）
```

#### 5.1.2 体素定位

给定一个世界坐标系下的三维点 `p = (x, y, z)`，其所在的体素位置计算如下：

```
voxel_x = floor(x / voxel_size)
voxel_y = floor(y / voxel_size)
voxel_z = floor(z / voxel_size)
```

体素位置使用 `VOXEL_LOC` 结构存储，并通过哈希映射 `unordered_map<VOXEL_LOC, OctoTree*>` 快速查找。

#### 5.1.3 平面拟合算法

对于每个体素（或八叉树叶子节点）中的点集，使用**主成分分析（PCA）**拟合平面：

```
输入: 点集 P = {p_1, p_2, ..., p_n}

1. 计算点集均值:
   μ = (1/n) · Σ p_i

2. 计算协方差矩阵:
   Σ = (1/n) · Σ (p_i - μ)(p_i - μ)^T    // 3×3矩阵

3. 特征值分解:
   Σ = V · Λ · V^T
   其中 Λ = diag(λ_1, λ_2, λ_3), λ_1 ≥ λ_2 ≥ λ_3

4. 平面判定:
   if λ_3 < plannar_threshold:
     is_plane = true
     法向量 n = V 中 λ_3 对应的特征向量
     平面中心 c = μ
     平面方程: n^T · (p - c) = 0, 即 n^T · p + d = 0, 其中 d = -n^T · c
   else:
     is_plane = false
     如果当前层 < max_layer，则继续划分为8个子节点
```

**平面判定标准**：最小特征值 λ₃ 反映了点集在法向方向上的离散程度。当 λ₃ 足够小（< plannar_threshold）时，表明点集近似分布在一个平面上。

#### 5.1.4 平面协方差计算

VoxelMap 不仅估计平面参数，还计算平面的不确定性（6×6协方差矩阵），考虑了每个点自身的测量噪声：

```
平面协方差 Σ_plane (6×6) 包含:
  ┌─────────────┬─────────────┐
  │ Σ_nn (3×3) │ Σ_nc (3×3)  │   Σ_nn: 法向量的协方差
  ├─────────────┼─────────────┤   Σ_nc: 法向-中心的交叉协方差
  │ Σ_cn (3×3) │ Σ_cc (3×3)  │   Σ_cc: 平面中心的协方差
  └─────────────┴─────────────┘

计算方法:
  Σ_plane = J · Σ_points · J^T

  其中 J 是平面参数对点坐标的Jacobian矩阵
  Σ_points 是所有点的联合协方差矩阵
```

#### 5.1.5 八叉树操作

```
初始化 (init_octo_tree):
  1. 收集体素内的点
  2. 如果点数 ≥ layer_point_size[layer]:
     尝试拟合平面
     如果拟合成功: 设置为叶子节点（平面）
     如果拟合失败且当前层 < max_layer: 划分为8个子节点

更新 (UpdateOctoTree):
  1. 如果是叶子节点（已检测到平面）:
     将新点加入
     如果新增点数 ≥ update_size_threshold:
       重新拟合平面
  2. 如果是非叶子节点:
     根据新点坐标，递归更新对应的子节点

查询最近平面 (BuildResidualListOMP):
  对于查询点 p:
  1. 计算 p 所在的体素坐标
  2. 在当前体素及其邻域中搜索
  3. 遍历八叉树各层，找到包含平面的叶子节点
  4. 计算点到平面的距离
  5. 如果距离 < 阈值（sigma_num * σ），则记录为有效匹配
```

### 5.2 VoxelMap++算法

#### 5.2.1 主要改进

VoxelMap++ 在VoxelMap的基础上做了以下优化：

1. **单层体素结构**：去掉了八叉树，每个体素只维护一个平面，减少内存消耗
2. **增量平面更新**：使用Welford算法增量计算协方差，避免重复计算
3. **并查集平面融合**：相邻体素中的共面平面自动融合
4. **多方向平面编码**：根据平面法向主方向选择最稳定的编码方式

#### 5.2.2 增量平面计算

VoxelMap++ 使用增量方式维护平面参数，无需存储所有点：

```
维护的统计量:
  xx = Σ x_i²,  yy = Σ y_i²,  zz = Σ z_i²
  xy = Σ x_i·y_i,  xz = Σ x_i·z_i,  yz = Σ y_i·z_i
  x = Σ x_i,  y = Σ y_i,  z = Σ z_i
  n = 点的数量

添加新点 (x_new, y_new, z_new):
  xx += x_new²,  yy += y_new²,  zz += z_new²
  xy += x_new·y_new,  xz += x_new·z_new,  yz += y_new·z_new
  x += x_new,  y += y_new,  z += z_new
  n += 1

从统计量计算平面参数:
  center = (x/n, y/n, z/n)
  cov = | xx/n - (x/n)²    xy/n - (x/n)(y/n)  xz/n - (x/n)(z/n) |
        | xy/n - (x/n)(y/n) yy/n - (y/n)²     yz/n - (y/n)(z/n) |
        | xz/n - (x/n)(z/n) yz/n - (y/n)(z/n) zz/n - (z/n)²     |
```

#### 5.2.3 平面编码方式

VoxelMap++ 根据法向量的主方向，选择不同的平面编码：

```
主方向选择:
  法向量 n = (n_x, n_y, n_z)
  如果 |n_z| 最大 → 主方向 0: z = -(a·x + b·y + d)，参数 ω = (a, b, d)
  如果 |n_y| 最大 → 主方向 1: y = -(a·x + b·z + d)，参数 ω = (a, b, d)
  如果 |n_x| 最大 → 主方向 2: x = -(a·y + b·z + d)，参数 ω = (a, b, d)

这种编码保证了在法向主方向上的数值稳定性，避免了当法向接近坐标轴时的奇异问题。
```

#### 5.2.4 并查集平面融合

```
对于每个体素 V:
  1. 检查V的6个相邻体素
  2. 对每个相邻体素 V_adj:
     如果 V 和 V_adj 都检测到了平面:
       计算两个平面之间的Mahalanobis距离
       如果距离 < 阈值:
         使用并查集将两个平面合并
         合并后的平面参数 = 加权平均

Find(node):  // 并查集查找根节点（带路径压缩）
  if node.parent ≠ node:
    node.parent = Find(node.parent)
  return node.parent

Union(A, B):  // 并查集合并
  root_A = Find(A)
  root_B = Find(B)
  if root_A ≠ root_B:
    合并平面参数
    root_B.parent = root_A
```

---

## 6. 点到平面匹配

### 6.1 匹配流程

```
输入: 降采样后的LiDAR点云, 当前位姿估计, 体素地图
输出: 点-平面匹配列表 ptpl_list

BuildResidualListOMP(voxel_map, points):

  // 使用OpenMP并行处理每个点
  #pragma omp parallel for
  for 每个点 p_i:
    1. 将 p_i 从LiDAR系变换到世界系: p_w = R·(R_LI·p_i + t_LI) + p
    2. 将协方差从LiDAR系传播到世界系

    3. 计算 p_w 所在的体素坐标 (vx, vy, vz)

    4. 在 (vx, vy, vz) 及其邻域中搜索平面:
       // 搜索当前体素及相邻体素
       for 每个候选体素:
         遍历八叉树各层（VoxelMap）或直接查找（VoxelMap++）
         如果找到有效平面 π:
           计算点到平面的距离 d = |n^T · (p_w - c)|
           如果 d < sigma_num × σ:  // σ是平面协方差的对应分量
             记录匹配 (p_i, π, d)

    5. 选择最佳匹配（距离最小的平面）
```

### 6.2 协方差传播

点的协方差需要从LiDAR坐标系传播到世界坐标系：

```
transformLidarCovToWorld(p_lidar, kf_state, Cov_lidar):

  // 获取当前状态
  R = state.rot           // 世界到机体的旋转
  R_LI = state.offset_R_L_I  // LiDAR到IMU的旋转
  t_LI = state.offset_T_L_I  // LiDAR到IMU的平移
  P = state.covariance    // 状态协方差

  // 计算变换后的协方差
  p_body = R_LI · p_lidar + t_LI

  // Jacobian矩阵
  J = [I_3×3, -(R · p_body)×, -(R · R_LI · p_lidar)×, R, ...]
  //  对 p     对 R             对 R_LI              对 t_LI

  // 协方差传播
  Cov_world = J · P · J^T + R · R_LI · Cov_lidar · (R · R_LI)^T

  return Cov_world
```

---

## 7. LiDAR点云协方差模型

### 7.1 测量噪声模型

每个LiDAR点的测量可以表示为极坐标 `(r, θ, φ)`，其中 r 为距离，θ 和 φ 为角度。测量噪声模型为：

```
Cov_spherical = diag(σ²_r, σ²_θ, σ²_φ)

其中:
  σ_r = ranging_cov       // 测距标准差（配置参数）
  σ_θ = angle_cov         // 角度标准差（配置参数）
  σ_φ = angle_cov         // 角度标准差（配置参数）
```

### 7.2 协方差转换到笛卡尔坐标系

```
calcLidarCov(p, ranging_cov, angle_cov):

  // 点的球坐标
  r = ||p||
  θ = arccos(p_z / r)
  φ = arctan2(p_y, p_x)

  // Jacobian: 球坐标到笛卡尔坐标
  J = ∂(x, y, z) / ∂(r, θ, φ)

  // 笛卡尔坐标系下的协方差
  Cov_cartesian = J · Cov_spherical · J^T

  return Cov_cartesian    // 3×3矩阵
```

---

## 8. 点云预处理算法

### 8.1 总体流程

```
原始LiDAR消息
      │
      ▼
┌─────────────────┐
│ 消息解析        │  根据LiDAR类型选择不同的处理函数:
│ (Handler)       │  avia_handler() / velodyne_handler() / oust64_handler()
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ 盲区过滤        │  去除距离 < blind 的点
│ (Blind Filter)  │  默认 blind = 0.1m
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ 降采样          │  每 point_filter_num 个点取一个
│ (Downsample)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ 特征提取        │  (可选) 识别平面和边缘特征
│ (Feature Ext.)  │  give_feature()
└────────┬────────┘
         │
         ▼
输出: 预处理后的点云
```

### 8.2 特征分类

预处理模块可以将点云中的每个点分类为以下特征类型：

| 特征类型 | 枚举值 | 说明 |
|----------|--------|------|
| `Nor` | 0 | 普通点（默认） |
| `Poss_Plane` | 1 | 可能的平面点 |
| `Real_Plane` | 2 | 确认的平面点 |
| `Edge_Jump` | 3 | 边缘跳变点 |
| `Edge_Plane` | 4 | 平面边缘点 |
| `Wire` | 5 | 线状物体点 |
| `ZeroPoint` | 6 | 无效点（距离为零） |

### 8.3 平面判定算法

```
plane_judge(点云 pl, 当前点索引 i_cur):

  1. 从 i_cur 开始，向后搜索 group_size 个点
  2. 对每组相邻点，计算两两之间的向量
  3. 计算法向量（叉积）和各方向的变化
  4. 如果组内所有点共面（法向量变化小）:
     标记为 Real_Plane
     返回平面法向量方向
  5. 否则:
     标记为 Nor 或 Poss_Plane
```

### 8.4 边缘跳变判定算法

```
edge_jump_judge(点云 pl, 类型 types, 点索引 i, 方向 nor_dir):

  1. 获取当前点与相邻点（前/后方向）
  2. 计算距离比值: ratio = d_next / d_curr
  3. 根据比值判断跳变类型:
     - Nr_nor: 正常（无跳变）
     - Nr_zero: 相邻点距离为零
     - Nr_180: 接近180度转向
     - Nr_inf: 距离突变（遮挡/边缘）
     - Nr_blind: 进入盲区

  4. 如果有跳变，标记为 Edge_Jump
```

---

## 9. SO(3) 旋转群数学

### 9.1 反对称矩阵

```
给定向量 v = [v_1, v_2, v_3]^T

反对称矩阵 [v]× = |  0   -v_3   v_2 |
                    |  v_3   0   -v_1 |
                    | -v_2  v_1    0   |

性质: [v]× · u = v × u（叉积）
```

### 9.2 指数映射（Rodrigues公式）

```
Exp: so(3) → SO(3)    // 角速度向量 → 旋转矩阵

给定旋转向量 φ = θ · â (θ为旋转角度, â为旋转轴单位向量):

R = Exp(φ) = I + sin(θ)/θ · [φ]× + (1 - cos(θ))/θ² · [φ]×²

特殊情况: 当 θ → 0 时
R ≈ I + [φ]×
```

### 9.3 对数映射

```
Log: SO(3) → so(3)    // 旋转矩阵 → 角速度向量

给定旋转矩阵 R:

θ = arccos((tr(R) - 1) / 2)

如果 θ ≈ 0:
  φ = [R_32 - R_23, R_13 - R_31, R_21 - R_12]^T / 2

否则:
  φ = θ / (2·sin(θ)) · [R_32 - R_23, R_13 - R_31, R_21 - R_12]^T
```

### 9.4 旋转矩阵到欧拉角转换

```
RotMtoEuler(R):
  // ZYX顺序（Yaw-Pitch-Roll）
  sy = sqrt(R[0][0]² + R[1][0]²)

  if sy > 1e-6:  // 非奇异
    roll  = atan2(R[2][1], R[2][2])
    pitch = atan2(-R[2][0], sy)
    yaw   = atan2(R[1][0], R[0][0])
  else:  // 万向锁
    roll  = atan2(-R[1][2], R[1][1])
    pitch = atan2(-R[2][0], sy)
    yaw   = 0

  return [roll, pitch, yaw]^T
```

---

## 10. 残差权重计算

### 10.1 PV-LIO-PLUS的修正

PV-LIO-PLUS 修正了原始PV-LIO中的误差传播公式。在计算观测噪声时，需要同时考虑：

1. **点的测量噪声**：LiDAR测量的距离和角度噪声
2. **平面的拟合噪声**：由平面拟合的不确定性引起的噪声
3. **状态估计噪声**：由位姿估计的不确定性引起的噪声

```
观测噪声计算:

r_i = σ²_point + σ²_plane

其中:
  // 点的测量噪声传播到观测空间
  σ²_point = n^T · J_transform · Cov_lidar · J_transform^T · n

  // 平面的拟合噪声传播到观测空间
  σ²_plane = J_plane · Cov_plane · J_plane^T

  J_transform: 坐标变换的Jacobian矩阵
  J_plane: 平面参数对观测的Jacobian矩阵
```

### 10.2 VoxelMap++中的权重

```
VoxelMap++ 使用3参数平面表示 ω = (a, b, d)

对于主方向0 (z = -(ax + by + d)):
  残差 h = p_z + a·p_x + b·p_y + d

  观测噪声:
  r = [p_x, p_y, 1] · Cov_ω · [p_x, p_y, 1]^T + 平面法向对应的点协方差分量

  其中 Cov_ω 是平面参数 ω 的3×3协方差矩阵
```
