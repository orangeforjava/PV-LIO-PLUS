# PV-LIO-PLUS 配置与使用指南

## 1. 环境要求

### 1.1 操作系统

| 操作系统 | 推荐版本 | ROS版本 |
|----------|----------|---------|
| Ubuntu 18.04 | ✅ 支持 | ROS Melodic |
| Ubuntu 20.04 | ✅ 推荐 | ROS Noetic |
| Ubuntu 16.04 | ⚠️ 最低支持 | ROS Kinetic |

### 1.2 依赖库

| 依赖 | 最低版本 | 说明 |
|------|----------|------|
| ROS | Melodic | 机器人操作系统 |
| PCL | 1.8 | 点云处理库（Ubuntu 18.04+ 默认满足） |
| Eigen | 3.3.4 | 线性代数库（Ubuntu 18.04+ 默认满足） |
| Boost | - | C++工具库（随ROS安装） |
| OpenMP | - | 并行计算（GCC默认支持） |
| C++ 编译器 | C++14 | GCC 5.0+ 或 Clang 3.4+ |

### 1.3 Livox ROS Driver

PV-LIO-PLUS 使用 Livox 激光雷达时需要安装 [livox_ros_driver](https://github.com/Livox-SDK/livox_ros_driver)。

```bash
# 安装Livox ROS Driver
cd ~/catkin_ws/src
git clone https://github.com/Livox-SDK/livox_ros_driver.git
cd ..
catkin_make

# 将驱动source到环境中（添加到 ~/.bashrc）
echo "source ~/catkin_ws/devel/setup.bash" >> ~/.bashrc
source ~/.bashrc
```

> **注意**：必须在编译和运行 PV-LIO-PLUS 之前 source Livox ROS Driver。

---

## 2. 编译安装

### 2.1 获取源代码

```bash
cd ~/catkin_ws/src
git clone https://github.com/orangeforjava/PV-LIO-PLUS.git
```

### 2.2 编译

```bash
cd ~/catkin_ws
catkin_make
source devel/setup.bash
```

### 2.3 编译选项说明

CMakeLists.txt 中的关键编译选项：

| 选项 | 值 | 说明 |
|------|-----|------|
| `CMAKE_BUILD_TYPE` | Release | 编译模式（Release为优化模式） |
| `CMAKE_CXX_FLAGS` | -std=c++14 -O3 | C++14标准，O3优化级别 |
| `MP_EN` | 宏定义 | 启用OpenMP多线程并行 |
| `MP_PROC_NUM` | 3 | 并行线程数 |
| `ROOT_DIR` | 自动生成 | 输出文件根目录 |

---

## 3. 运行

### 3.1 使用Launch文件启动

```bash
# 基本启动（包含RViz可视化）
roslaunch pv_lio_plus mapping_mid360.launch

# 不启动RViz
roslaunch pv_lio_plus mapping_mid360.launch rviz:=false

# 启用rosbag录制
roslaunch pv_lio_plus mapping_mid360.launch record:=true
```

### 3.2 播放数据集

```bash
# 在另一个终端中播放rosbag数据集
rosbag play your_dataset.bag
```

### 3.3 Launch文件参数

`launch/mapping_mid360.launch` 支持以下参数：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `rviz` | true | 是否启动RViz可视化 |
| `debug` | false | 是否启用调试模式 |
| `record` | false | 是否录制rosbag |

### 3.4 Launch文件结构解析

```xml
<launch>
    <arg name="rviz" default="true"/>
    <arg name="debug" default="false"/>
    <arg name="record" default="false"/>

    <!-- 加载YAML配置参数到ROS参数服务器 -->
    <rosparam command="load"
              file="$(find pv_lio_plus)/config/mid360_indoor.yaml"/>

    <!-- 启动主处理节点 -->
    <group if="$(eval arg('debug') == False)">
        <node pkg="pv_lio_plus" type="pv_lio_plus_node"
              name="pv_lio_plus_node" output="screen"/>
    </group>

    <!-- 启动RViz可视化 -->
    <group if="$(arg rviz)">
        <node pkg="rviz" type="rviz" name="rviz"
              args="-d $(find pv_lio_plus)/config/rviz_cfg/voxel_mapping.rviz"/>
    </group>

    <!-- 录制输出的点云话题 -->
    <group if="$(arg record)">
        <node pkg="rosbag" type="record" name="bag_record" output="screen"
              args="-O $(find pv_lio_plus)/../../output/pv_lio_plus.bag
                     /cloud_registered_lidar"/>
    </group>
</launch>
```

---

## 4. 配置文件详解

配置文件位于 `config/mid360_indoor.yaml`，包含以下部分：

### 4.1 通用配置（common）

```yaml
common:
    lid_topic: "/livox/lidar"    # LiDAR点云话题名
    imu_topic: "/livox/imu"      # IMU数据话题名
    time_sync_en: false          # 是否启用外部时间同步
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `lid_topic` | string | 订阅的LiDAR话题名。根据使用的LiDAR驱动设置 |
| `imu_topic` | string | 订阅的IMU话题名。必须与IMU驱动发布的话题一致 |
| `time_sync_en` | bool | 是否启用LiDAR和IMU的外部时间同步。如果传感器自身已同步，设为false |

### 4.2 预处理配置（preprocess）

```yaml
preprocess:
    lidar_type: 1                # LiDAR类型
    scan_line: 4                 # 扫描线数
    scan_rate: 10                # 扫描频率
    blind: 0.1                   # 盲区距离
    point_filter_num: 1          # 点过滤因子
```

| 参数 | 类型 | 取值范围 | 说明 |
|------|------|----------|------|
| `lidar_type` | int | 1/2/3 | LiDAR类型：1=Livox Avia, 2=Velodyne VLP-16, 3=Ouster OS1-64 |
| `scan_line` | int | 1-128 | 激光雷达的扫描线数。Livox Avia=4, Velodyne=16, Ouster=64 |
| `scan_rate` | int | Hz | 激光雷达的扫描频率。一般为10Hz |
| `blind` | double | 米 | 盲区距离。距离传感器小于此值的点将被过滤，避免噪声。默认0.1m |
| `point_filter_num` | int | ≥1 | 点云降采样因子。每N个点取1个。1=不降采样，4=取25%的点 |

#### LiDAR类型选择指南

| lidar_type | 传感器型号 | 消息类型 | 备注 |
|------------|-----------|----------|------|
| 1 | Livox Avia / Mid360 | `livox_ros_driver/CustomMsg` | 需安装 livox_ros_driver |
| 2 | Velodyne VLP-16 | `sensor_msgs/PointCloud2` | 标准PCL格式 |
| 3 | Ouster OS1-64 | `sensor_msgs/PointCloud2` | 标准PCL格式 |

### 4.3 建图配置（mapping）

```yaml
mapping:
    down_sample_size: 0.1        # 降采样大小
    max_iteration: 10            # 最大迭代次数
    voxel_size: 0.5              # 体素大小
    max_layer: 4                 # 最大八叉树层数
    layer_point_size: [5, 5, 5, 5, 5]  # 各层点数阈值
    plannar_threshold: 0.01      # 平面判定阈值
    max_points_size: 1000        # 每个平面最大点数
    max_cov_points_size: 1000    # 协方差计算最大点数
    update_size_threshold: 5     # 触发平面更新的最少新增点数
    sigma_num: 3                 # 异常值剔除阈值
    b_use_voxelmap_plus: true    # 是否使用VoxelMap++

    fov_degree: 360              # 视场角
    det_range: 100.0             # 最大检测距离
    extrinsic_est_en: false      # 是否在线估计外参
    extrinsic_T: [0.011, 0.02329, -0.04412]  # LiDAR-IMU平移外参
    extrinsic_R: [1, 0, 0, 0, 1, 0, 0, 0, 1]  # LiDAR-IMU旋转外参
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `down_sample_size` | double | 0.1 | VoxelGrid降采样的体素大小（米）。值越大降采样越多 |
| `max_iteration` | int | 10 | iESKF的最大迭代次数。通常5-10次即可收敛 |
| `voxel_size` | double | 0.5 | 体素地图的体素大小（米）。影响地图精度和内存 |
| `max_layer` | int | 4 | 八叉树最大层数（仅VoxelMap）。层数越多精度越高但计算量越大 |
| `layer_point_size` | list[int] | [5,5,5,5,5] | 各层初始化平面所需的最少点数 |
| `plannar_threshold` | double | 0.01 | 平面判定阈值（最小特征值）。值越小要求越严格 |
| `max_points_size` | int | 1000 | 每个平面允许的最大点数。超过后不再添加 |
| `max_cov_points_size` | int | 1000 | 用于协方差计算的最大点数 |
| `update_size_threshold` | int | 5 | 触发平面重新拟合的最少新增点数 |
| `sigma_num` | int | 3 | 异常值剔除阈值（σ的倍数）。3σ约覆盖99.7%的正常值 |
| `b_use_voxelmap_plus` | bool | true | 是否使用VoxelMap++算法。false则使用原始VoxelMap |
| `fov_degree` | double | 360 | 视场角（度）。360°为全向 |
| `det_range` | double | 100.0 | 最大检测距离（米）。超过此距离的点将被忽略 |
| `extrinsic_est_en` | bool | false | 是否在线估计LiDAR-IMU外参。建议先标定好设为false |
| `extrinsic_T` | list[3] | - | LiDAR到IMU的平移向量 [x, y, z]（米） |
| `extrinsic_R` | list[9] | - | LiDAR到IMU的旋转矩阵（行优先，3×3） |

#### 体素大小选择建议

| 场景 | voxel_size | down_sample_size | 说明 |
|------|-----------|-----------------|------|
| 室内小空间 | 0.25-0.5 | 0.05-0.1 | 精细建图 |
| 室内大空间 | 0.5-1.0 | 0.1-0.2 | 平衡精度和效率 |
| 室外场景 | 1.0-2.0 | 0.2-0.5 | 大范围覆盖 |
| 地下隧道 | 0.5 | 0.1 | 结构化环境 |

### 4.4 噪声模型配置（noise_model）

```yaml
noise_model:
    ranging_cov: 0.03            # 测距标准差
    angle_cov: 0.15              # 角度标准差
    acc_cov: 0.1                 # 加速度计噪声
    gyr_cov: 0.1                 # 陀螺仪噪声
    b_acc_cov: 0.0001            # 加速度计偏置漂移噪声
    b_gyr_cov: 0.0001            # 陀螺仪偏置漂移噪声
```

| 参数 | 类型 | 默认值 | 单位 | 说明 |
|------|------|--------|------|------|
| `ranging_cov` | double | 0.03 | 米 | LiDAR测距噪声标准差。值越大表示测距精度越低 |
| `angle_cov` | double | 0.15 | 弧度 | LiDAR角度噪声标准差。值越大表示角度精度越低 |
| `acc_cov` | double | 0.1 | m/s² | 加速度计白噪声标准差 |
| `gyr_cov` | double | 0.1 | rad/s | 陀螺仪白噪声标准差 |
| `b_acc_cov` | double | 0.0001 | m/s²/√s | 加速度计偏置随机游走噪声 |
| `b_gyr_cov` | double | 0.0001 | rad/s/√s | 陀螺仪偏置随机游走噪声 |

#### 噪声参数调节指南

- **测距噪声（ranging_cov）**：根据LiDAR规格书设置。Livox一般为0.02-0.05m
- **角度噪声（angle_cov）**：根据LiDAR规格书设置。一般为0.1-0.3 rad
- **IMU噪声（acc_cov, gyr_cov）**：
  - 高端MEMS IMU（如BMI088）：acc=0.05, gyr=0.05
  - 普通MEMS IMU：acc=0.1, gyr=0.1
  - 低端MEMS IMU：acc=0.5, gyr=0.5
- **偏置噪声（b_acc_cov, b_gyr_cov）**：一般设为0.0001-0.001，值越小表示偏置越稳定

### 4.5 发布配置（publish）

```yaml
publish:
    pub_voxel_map: false         # 是否发布体素地图可视化
    publish_max_voxel_layer: 2   # 发布的最大体素层数
    path_en: true                # 是否发布轨迹
    scan_publish_en: true        # 是否发布扫描点云
    dense_publish_en: false      # 是否发布密集点云
    scan_bodyframe_pub_en: false # 是否发布机体坐标系点云
    scan_lidarframe_pub_en: true # 是否发布LiDAR坐标系点云
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `pub_voxel_map` | bool | false | 发布体素地图平面可视化（MarkerArray）。开启会增加计算和通信负载 |
| `publish_max_voxel_layer` | int | 2 | 发布体素地图时的最大层数。层数越大细节越多 |
| `path_en` | bool | true | 是否发布运动轨迹（`/path` 话题） |
| `scan_publish_en` | bool | true | 是否发布当前帧的配准点云 |
| `dense_publish_en` | bool | false | 是否发布密集点云（未降采样）。开启会增加通信负载 |
| `scan_bodyframe_pub_en` | bool | false | 是否在机体坐标系下发布点云 |
| `scan_lidarframe_pub_en` | bool | true | 是否在LiDAR坐标系下发布点云 |

### 4.6 PCD保存配置（pcd_save）

```yaml
pcd_save:
    pcd_save_en: true            # 是否保存点云到PCD文件
    interval: -1                 # 保存间隔
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `pcd_save_en` | bool | true | 是否将全局地图保存为PCD文件 |
| `interval` | int | -1 | 保存间隔。-1=程序结束时保存全部；N=每N帧保存一次 |

> PCD文件保存路径为 `${ROOT_DIR}/scans.pcd`。ROOT_DIR在CMakeLists.txt中定义，默认为项目上两级目录的output文件夹。

---

## 5. ROS话题说明

### 5.1 订阅话题

| 话题名 | 消息类型 | 频率 | 说明 |
|--------|----------|------|------|
| `lid_topic` (配置) | `livox_ros_driver/CustomMsg` 或 `sensor_msgs/PointCloud2` | 10Hz | LiDAR点云数据 |
| `imu_topic` (配置) | `sensor_msgs/Imu` | 200-400Hz | IMU惯性测量数据 |

### 5.2 发布话题

| 话题名 | 消息类型 | 频率 | 说明 |
|--------|----------|------|------|
| `/cloud_registered` | `sensor_msgs/PointCloud2` | 10Hz | 世界坐标系下的配准点云 |
| `/cloud_registered_body` | `sensor_msgs/PointCloud2` | 10Hz | 机体坐标系下的点云 |
| `/cloud_registered_lidar` | `sensor_msgs/PointCloud2` | 10Hz | LiDAR坐标系下的点云 |
| `/cloud_effected` | `sensor_msgs/PointCloud2` | 10Hz | 有效匹配的点云 |
| `/Laser_map` | `sensor_msgs/PointCloud2` | ~1Hz | 全局地图点云 |
| `/Odometry` | `nav_msgs/Odometry` | 10Hz | 里程计（位置+姿态+协方差） |
| `/path` | `nav_msgs/Path` | 10Hz | 运动轨迹 |
| `/Extrinsic` | `nav_msgs/Odometry` | 10Hz | LiDAR-IMU外参 |
| `/planes` | `visualization_msgs/MarkerArray` | ~1Hz | 体素地图平面可视化（可选） |

### 5.3 TF变换

系统发布以下TF变换：

```
odom → body
  // odom: 里程计坐标系（与世界坐标系对齐）
  // body: 机体坐标系（IMU坐标系）
```

---

## 6. 传感器配置要求

### 6.1 LiDAR与IMU同步

**关键要求**：LiDAR和IMU的数据必须时间同步。

- **硬件同步**（推荐）：通过PPS信号或硬件触发器同步
- **软件同步**：通过`time_sync_en`参数启用ROS时间戳同步

> 如果运行时出现警告 "Failed to find match for field 'time'"，说明LiDAR点云中缺少各点的时间戳，这会影响去畸变的精度。

### 6.2 IMU频率要求

- **最低频率**：100Hz
- **推荐频率**：200-400Hz
- **高频率的好处**：更精确的运动去畸变、更快的IMU初始化

### 6.3 LiDAR-IMU外参标定

外参参数 `extrinsic_T` 和 `extrinsic_R` 定义了LiDAR坐标系到IMU坐标系的刚体变换：

```
p_imu = R_extrinsic * p_lidar + T_extrinsic
```

- **平移外参 extrinsic_T**：[x, y, z]，单位为米
- **旋转外参 extrinsic_R**：行优先的3×3旋转矩阵，以9个元素的列表表示

可以通过CAD模型或外参标定工具（如[lidar_imu_calib](https://github.com/hku-mars/lidar_IMU_calib)）获取。

如果LiDAR和IMU几乎同轴安装，可以设为：
```yaml
extrinsic_T: [0, 0, 0]
extrinsic_R: [1, 0, 0, 0, 1, 0, 0, 0, 1]
```

---

## 7. 常见问题

### 7.1 编译相关

**Q: 编译时Eigen报错**
A: Ubuntu 20.04 已修复此问题（2023.07.18 更新）。如仍有问题，确保 Eigen ≥ 3.3.4。

**Q: 找不到 livox_ros_driver**
A: 确保已安装并 source 了 livox_ros_driver：
```bash
source ~/ws_livox/devel/setup.bash
```

### 7.2 运行相关

**Q: 程序启动后没有输出**
A: 检查以下几点：
1. 确认LiDAR和IMU话题名与配置文件一致
2. 确认rosbag正在播放或传感器正在发布数据
3. 使用 `rostopic list` 确认话题存在
4. 使用 `rostopic hz <话题名>` 确认数据频率

**Q: 出现 "Failed to find match for field 'time'" 警告**
A: 说明LiDAR消息中缺少各点的时间戳字段。这对于点云去畸变很重要。请确保：
- 使用正确的LiDAR类型设置
- rosbag中的点云包含 time 字段

**Q: 里程计漂移严重**
A: 可尝试以下调整：
1. 增大 `max_iteration`（如15-20）
2. 减小 `voxel_size`（如0.25）
3. 检查IMU噪声参数是否过大
4. 确认LiDAR-IMU外参是否准确

**Q: 程序崩溃**
A: PV-LIO-PLUS 已修复了多个稳定性问题。如仍崩溃，请检查：
1. 是否使用了正确的 `lidar_type`
2. 数据集质量是否良好
3. IMU数据频率是否足够高

### 7.3 性能优化

**Q: 实时性不足**
A: 可尝试以下优化：
1. 增大 `down_sample_size`（减少处理点数）
2. 增大 `point_filter_num`（预处理降采样）
3. 减小 `max_iteration`（减少迭代次数）
4. 增大 `voxel_size`（减少体素数量）
5. 关闭 `pub_voxel_map`（减少发布负载）
6. 关闭 `dense_publish_en`（减少点云发布量）

---

## 8. VoxelMap vs VoxelMap++ 选择

通过配置参数 `b_use_voxelmap_plus` 切换两种算法：

| 特性 | VoxelMap | VoxelMap++ |
|------|---------|------------|
| 配置值 | `false` | `true` |
| 数据结构 | 八叉树（多层） | 单层体素 + 并查集 |
| 内存消耗 | 较高 | 较低 |
| 计算效率 | 一般 | 较高 |
| 平面拟合 | PCA全量计算 | 增量计算 |
| 平面合并 | 不支持 | 并查集自动合并 |
| 平面表示 | 法向量 + 中心 | ω向量（3参数） |
| 适用场景 | 通用 | 大规模环境 |

**推荐**：一般情况下建议使用 VoxelMap++（`b_use_voxelmap_plus: true`），其计算效率更高且内存占用更少。

---

## 9. 输出文件

### 9.1 PCD点云文件

当 `pcd_save_en: true` 时，系统将在程序结束时保存全局点云地图到PCD文件。

文件路径：`${ROOT_DIR}/scans.pcd`

> `ROOT_DIR` 在 `CMakeLists.txt` 中定义，默认为编译目录上两级的 `output/` 文件夹。

可使用以下工具查看PCD文件：
```bash
# PCL点云查看器
pcl_viewer scans.pcd

# 或使用CloudCompare
cloudcompare.CloudCompare scans.pcd
```

### 9.2 Rosbag录制

当启动时设置 `record:=true`，系统会录制 `/cloud_registered_lidar` 话题到rosbag文件：

文件路径：`${ROOT_DIR}/pv_lio_plus.bag`

---

## 10. 自定义传感器配置

如果使用不同的传感器组合，需要创建新的配置文件。以下是步骤：

### 10.1 创建新配置文件

```bash
# 复制模板
cp config/mid360_indoor.yaml config/my_sensor.yaml
```

### 10.2 修改关键参数

1. 设置正确的话题名（`lid_topic`, `imu_topic`）
2. 设置正确的LiDAR类型（`lidar_type`）
3. 设置正确的扫描线数（`scan_line`）
4. 设置正确的LiDAR-IMU外参（`extrinsic_T`, `extrinsic_R`）
5. 根据传感器规格调整噪声参数

### 10.3 创建新Launch文件

```bash
cp launch/mapping_mid360.launch launch/mapping_my_sensor.launch
```

修改Launch文件中加载的配置文件路径：
```xml
<rosparam command="load"
          file="$(find pv_lio_plus)/config/my_sensor.yaml"/>
```

### 10.4 运行

```bash
roslaunch pv_lio_plus mapping_my_sensor.launch
```
