# PV-LIO-PLUS 代码结构与模块说明

## 1. 源文件总览

| 文件路径 | 代码行数 | 功能描述 |
|----------|----------|----------|
| `src/voxelMapping.cpp` | ~1228行 | 主算法实现，包含ROS节点、主循环、观测模型 |
| `src/preprocess.cpp` | ~945行 | LiDAR点云预处理，特征提取 |
| `src/preprocess.h` | ~100行 | 预处理类和数据结构定义 |
| `src/IMU_Processing.hpp` | ~400行 | IMU处理、预积分、点云去畸变 |
| `include/common_lib.h` | ~300行 | 通用数据结构、常量、数学工具 |
| `include/use-ikfom.hpp` | ~200行 | 系统状态流形定义、动力学方程 |
| `include/so3_math.h` | ~150行 | SO(3)旋转群数学运算 |
| `include/voxel_map_util.hpp` | ~1500行 | VoxelMap算法实现 |
| `include/voxelmapplus_util.hpp` | ~1200行 | VoxelMap++算法实现 |
| `include/IKFoM_toolkit/esekfom/esekfom.hpp` | ~800行 | 迭代扩展卡尔曼滤波器核心实现 |

---

## 2. 主控模块：`src/voxelMapping.cpp`

### 2.1 功能概述

主控模块是整个系统的入口和核心调度中心，负责：
- ROS节点初始化和参数加载
- 传感器数据的接收和缓冲
- LiDAR与IMU数据的时间同步
- 主处理循环的调度
- 观测模型的计算
- 结果的发布

### 2.2 全局变量与配置参数

```cpp
// ===== 系统常量 =====
#define DET_RANGE       300.0f    // 检测范围（米）
#define MOV_THRESHOLD   1.5f      // 运动阈值
#define INIT_TIME       0.1       // 初始化等待时间（秒）
#define LASER_POINT_COV 0.001     // 激光点默认协方差
#define MAXN            720000    // 最大点云数量
#define PUBFRAME_PERIOD 20        // 发布周期

// ===== 数据缓冲区 =====
deque<PointCloudXYZI::Ptr> lidar_buffer;    // LiDAR点云缓冲队列
deque<double> time_buffer;                   // 时间戳缓冲队列
deque<sensor_msgs::Imu::ConstPtr> imu_buffer; // IMU消息缓冲队列

// ===== 核心数据 =====
PointCloudXYZI::Ptr feats_undistort;         // 去畸变后的点云
PointCloudXYZI::Ptr feats_undistort_down;    // 降采样后的点云
vector<M3D> var_down_lidar;                  // 降采样点的LiDAR系协方差
vector<M3D> var_down_world;                  // 降采样点的世界系协方差

// ===== 体素地图 =====
unordered_map<VOXEL_LOC, OctoTree*> voxel_map;         // VoxelMap
unordered_map<VOXEL_LOC, UnionFindNode*> voxel_map_plus; // VoxelMap++

// ===== 配置参数 =====
double voxel_size;             // 体素大小
int max_layer;                 // 最大层数
vector<int> layer_point_size;  // 各层点数阈值
double plannar_threshold;      // 平面判定阈值
int max_points_size;           // 最大点数
double sigma_num;              // 异常值剔除阈值（σ倍数）
bool b_use_voxelmap_plus;      // 是否使用VoxelMap++
double ranging_cov;            // 测距噪声
double angle_cov;              // 角度噪声
```

### 2.3 函数详解

#### 2.3.1 坐标变换函数

```cpp
// 将LiDAR系点变换到世界系（使用iKFoM状态）
void pointBodyToWorld_ikfom(PointType const *const pi,
                            PointType *const po,
                            state_ikfom &s)
// 计算: po = s.rot * (s.offset_R_L_I * pi + s.offset_T_L_I) + s.pos

// 将LiDAR系点变换到世界系（使用全局状态）
void pointLidarToWorld(PointType const *const pi, PointType *const po)
// 使用全局变量 state_point 进行变换

// 将带颜色的点从机体系变换到世界系
void RGBpointBodyToWorld(PointType const *const pi, PointType *const po)
// 同上，但保留RGB信息

// 将LiDAR系点变换到IMU（机体）系
void pointLidarToIMU(PointType const *const pi, PointType *const po)
// 计算: po = Lidar_R_wrt_IMU * pi + Lidar_T_wrt_IMU

// 批量变换整个点云
void transformLidar2World(const state_ikfom &state_point,
                          const PointCloudXYZI::Ptr &input_cloud,
                          PointCloudXYZI::Ptr &trans_cloud)

// 将LiDAR系点协方差传播到世界系
M3D transformLidarCovToWorld(Eigen::Vector3d &p_lidar,
                             const esekfom::esekf<state_ikfom, 12, input_ikfom> &kf,
                             const Eigen::Matrix3d &COV_lidar)
// 使用Jacobian矩阵进行协方差传播:
// Cov_world = J * P * J^T + R * R_LI * Cov_lidar * (R * R_LI)^T
```

#### 2.3.2 消息回调函数

```cpp
// 标准PointCloud2消息回调（Velodyne/Ouster）
void standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
// 流程:
// 1. 调用预处理模块解析点云
// 2. 将处理后的点云加入 lidar_buffer
// 3. 将时间戳加入 time_buffer

// Livox自定义消息回调
void livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg)
// 流程:
// 1. 调用预处理模块解析Livox点云
// 2. 将处理后的点云加入 lidar_buffer
// 3. 将时间戳加入 time_buffer

// IMU消息回调
void imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in)
// 流程:
// 1. 将IMU消息加入 imu_buffer
// 2. 触发信号，唤醒主循环
```

#### 2.3.3 数据同步函数

```cpp
// 同步LiDAR和IMU数据包
bool sync_packages(MeasureGroup &meas)
// 流程:
// 1. 检查 lidar_buffer 和 imu_buffer 是否非空
// 2. 取出最早的LiDAR帧
// 3. 收集该帧时间范围内的所有IMU消息
// 4. 如果IMU数据不足（时间未覆盖LiDAR帧），等待更多数据
// 5. 填充 MeasureGroup 结构并返回
//
// 输出:
//   meas.lidar: 当前帧点云
//   meas.lidar_beg_time: 帧开始时间
//   meas.lidar_end_time: 帧结束时间
//   meas.imu: 对应时间段的IMU消息队列
```

#### 2.3.4 观测模型函数

```cpp
// VoxelMap观测模型
void observation_model_share(state_ikfom &s,
                             esekfom::dyn_share_datastruct<double> &ekfom_data)
// 流程:
// 1. 将降采样点云变换到世界系
// 2. 将协方差传播到世界系
// 3. 构建带协方差的点列表 pv_list
// 4. 调用 BuildResidualListOMP() 搜索最近平面
// 5. 对每个匹配的点-平面对:
//    a. 计算残差 h = n^T * p_w + d
//    b. 计算Jacobian矩阵 H (1×23)
//    c. 计算观测噪声 R（考虑点和平面的协方差）
// 6. 组装观测方程传递给EKF

// VoxelMap++观测模型
void observation_model_share_plus(state_ikfom &s,
                                  esekfom::dyn_share_datastruct<double> &ekfom_data)
// 流程类似，但使用VoxelMap++的平面表示和残差计算方式
// 主要区别:
// - 使用3参数平面编码（ω向量）
// - 根据主方向选择不同的残差公式
// - 使用3×3平面协方差矩阵
```

#### 2.3.5 发布函数

```cpp
// 发布世界坐标系下的点云
void publish_frame_world(const ros::Publisher &pubLaserCloudFull)
// 发布到 /cloud_registered 话题

// 发布机体坐标系下的点云
void publish_frame_body(const ros::Publisher &_pub)
// 发布到 /cloud_registered_body 话题

// 发布LiDAR坐标系下的点云
void publish_frame_lidar(const ros::Publisher &_pub)
// 发布到 /cloud_registered_lidar 话题

// 发布全局地图
void publish_map(const ros::Publisher &pubLaserCloudMap)
// 发布到 /Laser_map 话题

// 发布里程计
void publish_odometry(const ros::Publisher &pubOdomAftMapped)
// 发布到 /Odometry 话题
// 同时广播 TF 变换: odom → body

// 发布运动轨迹
void publish_path(const ros::Publisher pubPath)
// 发布到 /path 话题
```

#### 2.3.6 主函数（main）

```cpp
int main(int argc, char **argv)
// 完整的初始化和主循环:
//
// === 初始化阶段 ===
// 1. 初始化ROS节点
// 2. 从参数服务器加载配置参数
// 3. 创建订阅者和发布者
// 4. 初始化预处理模块 (Preprocess)
// 5. 初始化IMU处理模块 (ImuProcess)
// 6. 初始化状态估计器 (esekf)
// 7. 初始化VoxelGrid降采样滤波器
//
// === 主循环 ===
// while (ros::ok()):
//   1. sync_packages(): 同步LiDAR和IMU数据
//   2. p_imu->Process(): IMU传播和点云去畸变
//   3. 首帧处理: 初始化VoxelMap/VoxelMap++
//   4. downSizeFilterSurf.filter(): 降采样
//   5. calcLidarCov(): 计算各点的LiDAR系协方差
//   6. kf.update_iterated_dyn_share_diagonal(): 迭代EKF更新
//   7. updateVoxelMap[OMP/Plus](): 更新体素地图
//   8. 发布里程计、路径、点云
//   9. (可选) 保存PCD文件
```

---

## 3. 预处理模块：`src/preprocess.h` + `src/preprocess.cpp`

### 3.1 类定义

```cpp
class Preprocess {
public:
    // ===== 公有接口 =====

    // 处理Livox自定义消息
    void process(const livox_ros_driver::CustomMsg::ConstPtr &msg,
                 PointCloudXYZI::Ptr &pcl_out);

    // 处理标准PointCloud2消息
    void process(const sensor_msgs::PointCloud2::ConstPtr &msg,
                 PointCloudXYZI::Ptr &pcl_out);

    // 设置预处理参数
    void set(bool feat_en, int lid_type, double bld, int pfilt_num);

    // ===== 公有成员 =====
    int lidar_type;              // LiDAR类型: 1=AVIA, 2=VELO16, 3=OUST64
    int point_filter_num;        // 点过滤因子（每N个点取1个）
    int N_SCANS;                 // 扫描线数
    int SCAN_RATE;               // 扫描频率(Hz)
    double blind;                // 盲区距离(m)
    bool feature_enabled;        // 是否启用特征提取
    bool given_offset_time;      // 是否提供时间偏移

    PointCloudXYZI pl_surf;      // 表面特征点云
    PointCloudXYZI pl_corn;      // 角点特征点云
    PointCloudXYZI pl_full;      // 完整点云

private:
    // ===== LiDAR处理函数 =====

    // Livox Avia处理器
    void avia_handler(const livox_ros_driver::CustomMsg::ConstPtr &msg);

    // Velodyne VLP-16处理器
    void velodyne_handler(const sensor_msgs::PointCloud2::ConstPtr &msg);

    // Ouster OS1-64处理器
    void oust64_handler(const sensor_msgs::PointCloud2::ConstPtr &msg);

    // ===== 特征提取函数 =====

    // 对单条扫描线进行特征提取
    void give_feature(PointCloudXYZI &pl, vector<orgtype> &types);

    // 判断一组相邻点是否构成平面
    int plane_judge(const PointCloudXYZI &pl, vector<orgtype> &types,
                    uint i_cur, uint &i_nex, Eigen::Vector3d &curr_direct);

    // 识别小尺度平面
    bool small_plane(const PointCloudXYZI &pl, vector<orgtype> &types,
                     uint i_cur, uint &i_nex, Eigen::Vector3d &curr_direct);

    // 判断相邻点之间是否存在边缘跳变
    bool edge_jump_judge(const PointCloudXYZI &pl, vector<orgtype> &types,
                         uint i, Surround nor_dir);

    // ===== 私有成员 =====
    PointCloudXYZI pl_buff[128];     // 各扫描线的点云缓冲
    vector<orgtype> typess[128];      // 各扫描线的点类型信息

    int group_size;                   // 特征分析组大小（默认8）
    double disA, disB;               // 距离判断参数
    double p2l_ratio;                // 点到线比值阈值
    double limit_maxmid;             // 最大-中间特征值比限制
    double limit_midmin;             // 中间-最小特征值比限制
    double limit_maxmin;             // 最大-最小特征值比限制
    double jump_up_limit;            // 向上跳变角度限制
    double jump_down_limit;          // 向下跳变角度限制
    double cos160;                   // cos(160°)
    double edgea, edgeb;             // 边缘判断参数
    double smallp_intersect;         // 小平面交叉阈值
    double smallp_ratio;             // 小平面比值阈值
};
```

### 3.2 数据结构

```cpp
// 点的特征类型枚举
enum Feature {
    Nor,          // 普通点
    Poss_Plane,   // 可能的平面点
    Real_Plane,   // 确认的平面点
    Edge_Jump,    // 边缘跳变点
    Edge_Plane,   // 平面边缘点
    Wire,         // 线状物体点
    ZeroPoint     // 无效点
};

// 激光雷达类型枚举
enum LID_TYPE { AVIA = 1, VELO16 = 2, OUST64 = 3 };

// 边缘跳变类型枚举
enum E_jump {
    Nr_nor,    // 正常（无跳变）
    Nr_zero,   // 距离为零
    Nr_180,    // 接近180度转向
    Nr_inf,    // 距离突变（遮挡/边缘）
    Nr_blind   // 进入盲区
};

// 相邻方向枚举
enum Surround { Prev, Next };  // 前一个点 / 后一个点

// 点的原始信息结构
struct orgtype {
    double range;          // 点到原点的距离
    double dista;          // 到下一个点的距离
    double angle[2];       // 与前后点的夹角 [前, 后]
    double intersect;      // 相邻边的交叉值
    E_jump edj[2];         // 前后跳变类型 [前, 后]
    Feature ftype;         // 特征类型
};

// Velodyne点结构（带ring和time信息）
namespace velodyne_ros {
    struct EIGEN_ALIGN16 Point {
        PCL_ADD_POINT4D;          // x, y, z, padding
        float intensity;           // 反射强度
        float time;                // 相对时间偏移
        uint16_t ring;             // 扫描线编号
    };
}

// Ouster点结构（带丰富元信息）
namespace ouster_ros {
    struct EIGEN_ALIGN16 Point {
        PCL_ADD_POINT4D;          // x, y, z, padding
        float intensity;           // 反射强度
        uint32_t t;                // 时间戳
        uint16_t reflectivity;     // 反射率
        uint8_t ring;              // 扫描线编号
        uint16_t ambient;          // 环境光
        uint32_t range;            // 距离
    };
}
```

---

## 4. IMU处理模块：`src/IMU_Processing.hpp`

### 4.1 类定义

```cpp
class ImuProcess {
public:
    // ===== 构造和重置 =====
    ImuProcess();                      // 构造函数
    ~ImuProcess();                     // 析构函数
    void Reset();                      // 重置所有状态
    void Reset(double start_timestamp, // 重置到指定时间
               const sensor_msgs::ImuConstPtr &lastimu);

    // ===== 参数设置 =====
    void set_extrinsic(const V3D &transl, const M3D &rot);   // 设置LiDAR-IMU外参
    void set_extrinsic(const V3D &transl);                    // 设置仅平移外参
    void set_extrinsic(const MD(4,4) &T);                     // 设置4x4变换矩阵外参
    void set_gyr_cov(const V3D &scaler);                      // 设置陀螺仪噪声
    void set_acc_cov(const V3D &scaler);                      // 设置加速度计噪声
    void set_gyr_bias_cov(const V3D &b_g);                    // 设置陀螺仪偏置噪声
    void set_acc_bias_cov(const V3D &b_a);                    // 设置加速度计偏置噪声

    // ===== 主处理函数 =====
    void Process(const MeasureGroup &meas,                    // 主处理入口
                 esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state,
                 PointCloudXYZI::Ptr pcl_un_);

    // ===== 公有成员 =====
    Eigen::Matrix<double, 12, 12> Q;   // 过程噪声协方差矩阵 (12×12)
    V3D cov_acc, cov_gyr;              // 实际使用的加速度和陀螺仪噪声
    V3D cov_acc_scale, cov_gyr_scale;  // 噪声的缩放因子
    V3D cov_bias_gyr, cov_bias_acc;    // 偏置漂移噪声
    double first_lidar_time;           // 第一帧LiDAR的时间戳

private:
    // ===== 内部函数 =====

    // IMU静止初始化
    void IMU_init(const MeasureGroup &meas,
                  esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state,
                  int &N);
    // 步骤:
    // 1. 收集约200帧IMU数据
    // 2. 计算加速度/陀螺仪均值和协方差
    // 3. 估计初始重力方向
    // 4. 初始化EKF状态和协方差

    // 点云运动去畸变
    void UndistortPcl(const MeasureGroup &meas,
                      esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state,
                      PointCloudXYZI &pcl_in_out);
    // 步骤:
    // 1. 前向传播: 利用IMU数据逐步计算各时刻的位姿
    // 2. 后向补偿: 将每个LiDAR点变换到扫描结束时刻

    // ===== 内部成员 =====
    PointCloudXYZI::Ptr cur_pcl_un_;        // 当前去畸变的点云
    sensor_msgs::ImuConstPtr last_imu_;     // 上一个IMU消息
    deque<sensor_msgs::ImuConstPtr> v_imu_; // IMU消息暂存队列
    vector<Pose6D> IMUpose;                 // 各IMU时刻的位姿历史
    vector<M3D> v_rot_pcl_;                 // 点云旋转矩阵列表
    M3D Lidar_R_wrt_IMU;                    // LiDAR到IMU的旋转矩阵
    V3D Lidar_T_wrt_IMU;                    // LiDAR到IMU的平移向量
    V3D mean_acc, mean_gyr;                 // 加速度和陀螺仪均值
    V3D angvel_last, acc_s_last;            // 上一时刻角速度和加速度
    double last_lidar_end_time_;            // 上一帧LiDAR的结束时间
    int init_iter_num;                      // 初始化迭代计数
    bool b_first_frame_;                    // 是否为第一帧
    bool imu_need_init_;                    // IMU是否需要初始化
};
```

---

## 5. 通用工具模块：`include/common_lib.h`

### 5.1 宏定义和常量

```cpp
#define PI_M             (3.14159265358)       // 圆周率
#define G_m_s2           (9.81)                // 重力加速度 (m/s²)
#define DIM_STATE        (18)                  // 状态维数
#define DIM_PROC_N       (12)                  // 过程噪声维数
#define CUBE_LEN         (6.0)                 // 立方体长度
#define LIDAR_SP_LEN     (2)                   // LiDAR特征长度
#define INIT_COV         (1)                   // 初始协方差值
#define NUM_MATCH_POINTS (5)                   // 平面匹配所需最小点数
#define MAX_MEAS_DIM     (10000)               // 最大测量维数
```

### 5.2 类型定义

```cpp
typedef pcl::PointXYZINormal PointType;        // 点类型（含法向量）
typedef pcl::PointCloud<PointType> PointCloudXYZI;  // 点云类型
typedef Vector3d V3D;                          // 3维double向量
typedef Matrix3d M3D;                          // 3×3 double矩阵
typedef Vector3f V3F;                          // 3维float向量
typedef Matrix3f M3F;                          // 3×3 float矩阵
```

### 5.3 核心数据结构

```cpp
// 测量数据组：包含一帧LiDAR数据和对应的IMU数据
struct MeasureGroup {
    double lidar_beg_time;                     // LiDAR帧起始时间
    double lidar_end_time;                     // LiDAR帧结束时间
    PointCloudXYZI::Ptr lidar;                 // LiDAR点云指针
    deque<sensor_msgs::Imu::ConstPtr> imu;     // IMU消息队列
};

// 6自由度位姿：记录某一时刻的完整运动状态
struct Pose6D {
    double offset_time;     // 相对于帧起始时间的偏移
    double acc[3];          // 加速度 (x, y, z)
    double gyr[3];          // 角速度 (x, y, z)
    double vel[3];          // 线速度 (x, y, z)
    double pos[3];          // 位置 (x, y, z)
    double rot[9];          // 旋转矩阵 (行优先, 3×3)
};

// 状态组：包含完整的系统状态（与iKFoM并行的状态表示）
struct StatesGroup {
    M3D rot_end;            // 旋转矩阵 (3×3)
    V3D pos_end;            // 位置向量 (3×1)
    V3D vel_end;            // 速度向量 (3×1)
    V3D bias_g;             // 陀螺仪偏置 (3×1)
    V3D bias_a;             // 加速度计偏置 (3×1)
    V3D gravity;            // 重力向量 (3×1)
    Matrix<double, DIM_STATE, DIM_STATE> cov;  // 状态协方差 (18×18)

    // 运算符重载: 支持流形上的加法和减法
    StatesGroup operator+(const Matrix<double, DIM_STATE, 1> &state_add);
    StatesGroup &operator+=(const Matrix<double, DIM_STATE, 1> &state_add);
    Matrix<double, DIM_STATE, 1> operator-(const StatesGroup &b);
};
```

### 5.4 工具函数

```cpp
// 旋转矩阵转欧拉角 (Yaw-Pitch-Roll)
Eigen::Vector3d R2ypr(const Eigen::Matrix3d &R);

// 欧拉角转旋转矩阵
Eigen::Matrix3d ypr2R(const Eigen::Vector3d &ypr);

// 重力向量转旋转矩阵
Eigen::Matrix3d g2R(const Eigen::Vector3d &g);

// 法向量拟合：最小二乘法拟合平面法向量
bool esti_normvector(Vector3d &normvec, const PointVector &point,
                     const double &threshold, const int &point_num);

// 平面拟合：PCA方法拟合平面
bool esti_plane(Vector4d &pca_result, const PointVector &point,
                const double &threshold);

// 两点距离
float calc_dist(PointType p1, PointType p2);

// 3×3矩阵的伴随矩阵
void adjugateM3D(const Matrix3d &s, Matrix3d &t);
```

---

## 6. 状态流形模块：`include/use-ikfom.hpp`

### 6.1 状态定义

```cpp
// 使用MTK（Manifold ToolKit）宏定义系统状态流形
MTK_BUILD_MANIFOLD(state_ikfom,
    ((vect3, pos))           // 位置 ∈ ℝ³ (3维)
    ((SO3, rot))             // 旋转 ∈ SO(3) (3维切空间)
    ((SO3, offset_R_L_I))    // LiDAR-IMU旋转外参 ∈ SO(3) (3维切空间)
    ((vect3, offset_T_L_I))  // LiDAR-IMU平移外参 ∈ ℝ³ (3维)
    ((vect3, vel))           // 速度 ∈ ℝ³ (3维)
    ((vect3, bg))            // 陀螺仪偏置 ∈ ℝ³ (3维)
    ((vect3, ba))            // 加速度计偏置 ∈ ℝ³ (3维)
    ((S2, grav))             // 重力 ∈ S² (2维切空间)
);
// 总误差状态维度: 3+3+3+3+3+3+3+2 = 23

// 系统输入（IMU测量）
MTK_BUILD_MANIFOLD(input_ikfom,
    ((vect3, acc))           // 加速度测量 ∈ ℝ³
    ((vect3, gyro))          // 陀螺仪测量 ∈ ℝ³
);

// 过程噪声
MTK_BUILD_MANIFOLD(process_noise_ikfom,
    ((vect3, ng))            // 陀螺仪噪声 ∈ ℝ³
    ((vect3, na))            // 加速度计噪声 ∈ ℝ³
    ((vect3, nbg))           // 陀螺仪偏置漂移 ∈ ℝ³
    ((vect3, nba))           // 加速度计偏置漂移 ∈ ℝ³
);
```

### 6.2 动力学函数

```cpp
// 状态转移函数：计算状态的时间导数
Eigen::Matrix<double, 24, 1> get_f(state_ikfom &s, const input_ikfom &in)
// 返回:
//   f[0:2]   = s.vel                                    // dp/dt = v
//   f[3:5]   = (in.gyro - s.bg)                         // dR/dt 对应的角速度
//   f[6:8]   = 0                                        // dR_LI/dt = 0
//   f[9:11]  = 0                                        // dt_LI/dt = 0
//   f[12:14] = s.rot * (in.acc - s.ba) + s.grav         // dv/dt
//   f[15:17] = 0                                        // dbg/dt = 0 (由噪声驱动)
//   f[18:20] = 0                                        // dba/dt = 0 (由噪声驱动)
//   f[21:23] = 0                                        // dg/dt = 0

// 状态转移Jacobian：∂f/∂x
Eigen::Matrix<double, 24, 23> df_dx(state_ikfom &s, const input_ikfom &in)
// 非零元素:
//   df_pos/d_vel = I₃                                  // 位置对速度的偏导
//   df_rot/d_bg = -I₃                                  // 旋转对陀螺仪偏置的偏导
//   df_vel/d_rot = -R·[a_m - ba]×                      // 速度对旋转的偏导
//   df_vel/d_ba = -R                                   // 速度对加速度偏置的偏导
//   df_vel/d_grav = I₃ (S2 Jacobian)                   // 速度对重力的偏导

// 噪声Jacobian：∂f/∂w
Eigen::Matrix<double, 24, 12> df_dw(state_ikfom &s, const input_ikfom &in)
// 非零元素:
//   df_rot/d_ng = -I₃                                  // 旋转对陀螺仪噪声
//   df_vel/d_na = -R                                   // 速度对加速度噪声
//   df_bg/d_nbg = I₃                                   // 偏置对偏置噪声
//   df_ba/d_nba = I₃

// SO(3)四元数到欧拉角转换
vect3 SO3ToEuler(const SO3 &orient)
// 将四元数转换为 [roll, pitch, yaw] 欧拉角

// 过程噪声协方差矩阵
MTK::get_cov<process_noise_ikfom>::type process_noise_cov()
// 返回 12×12 对角矩阵:
// diag(σ²_ng, σ²_ng, σ²_ng, σ²_na, ..., σ²_nbg, ..., σ²_nba, ...)
```

---

## 7. VoxelMap模块：`include/voxel_map_util.hpp`

### 7.1 数据结构

```cpp
// 体素位置（用作哈希表的键）
class VOXEL_LOC {
public:
    int64_t x, y, z;              // 体素的整数坐标
    // 哈希函数: hash = ((z * HASH_P) % MAX_N + y) * HASH_P % MAX_N + x
};

// 带协方差的3D点
struct pointWithCov {
    Eigen::Vector3d point_lidar;   // LiDAR坐标系下的点坐标
    Eigen::Vector3d point_world;   // 世界坐标系下的点坐标
    Eigen::Matrix3d cov_lidar;     // LiDAR坐标系下的协方差
    Eigen::Matrix3d cov_world;     // 世界坐标系下的协方差
};

// 平面结构
struct Plane {
    Eigen::Vector3d center;        // 平面中心
    Eigen::Vector3d normal;        // 平面法向量
    Eigen::Vector3d y_normal;      // 平面内第二方向（中等特征值）
    Eigen::Vector3d x_normal;      // 平面内第一方向（最大特征值）
    Eigen::Matrix3d covariance;    // 点集协方差矩阵
    Eigen::Matrix<double,6,6> plane_cov;  // 平面参数协方差 (6×6)
    float radius;                  // 平面半径（最大特征值的平方根）
    float min_eigen_value;         // 最小特征值（法向方差）
    float mid_eigen_value;         // 中等特征值
    float max_eigen_value;         // 最大特征值（平面范围）
    float d;                       // 平面方程常数项
    int points_size;               // 平面包含的点数
    bool is_plane;                 // 是否判定为平面
    bool is_init;                  // 是否已初始化
    bool is_update;                // 是否已更新
    int id;                        // 平面ID
};

// 点到平面匹配结构
struct ptpl {
    Eigen::Vector3d point;         // LiDAR系点坐标
    Eigen::Vector3d point_world;   // 世界系点坐标
    Eigen::Vector3d normal;        // 匹配平面的法向量
    Eigen::Vector3d center;        // 匹配平面的中心
    Eigen::Matrix<double,6,6> plane_cov;  // 平面协方差
    double d;                      // 平面方程常数项
    int layer;                     // 八叉树层数
    Eigen::Matrix3d cov_lidar;     // 点的LiDAR系协方差
    Eigen::Matrix3d cov_world;     // 点的世界系协方差
};
```

### 7.2 八叉树类

```cpp
class OctoTree {
public:
    // ===== 数据成员 =====
    vector<pointWithCov> temp_points_;   // 所有点（用于平面初始化）
    vector<pointWithCov> new_points_;    // 新增点（用于平面更新）
    Plane *plane_ptr_;                   // 指向平面对象
    int max_layer_;                      // 允许的最大层数
    int layer_;                          // 当前层编号
    int octo_state_;                     // 状态: 0=叶子(平面), 1=非叶子
    OctoTree *leaves_[8];                // 8个子节点指针
    double voxel_center_[3];             // 体素中心坐标
    float quater_length_;                // 四分之一边长
    int all_points_num_;                 // 总点数
    int new_points_num_;                 // 新增点数
    int max_points_size_;                // 最大点数阈值
    int max_cov_points_size_;            // 最大协方差计算点数
    float planer_threshold_;             // 平面判定阈值
    bool update_enable_;                 // 是否允许更新
    bool update_cov_enable_;             // 是否允许更新协方差

    // ===== 构造函数 =====
    OctoTree(int max_layer, int layer, vector<int> layer_point_size,
             int max_points_size, int max_cov_points_size,
             float planer_threshold);

    // ===== 成员函数 =====

    // 初始化平面：对点集进行PCA拟合
    void init_plane(const vector<pointWithCov> &points, Plane *plane);
    // 步骤:
    // 1. 计算点集均值
    // 2. 计算协方差矩阵
    // 3. 特征值分解
    // 4. 判断是否为平面
    // 5. 计算平面参数和协方差

    // 更新平面：加入新点后重新拟合
    void update_plane(const vector<pointWithCov> &points, Plane *plane);

    // 初始化八叉树节点
    void init_octo_tree();
    // 步骤:
    // 1. 尝试拟合平面
    // 2. 如果失败且未达最大层数，划分为8个子节点

    // 递归划分八叉树
    void cut_octo_tree();
    // 将当前节点的点分配到8个子节点中

    // 动态更新八叉树
    void UpdateOctoTree(const pointWithCov &pv);
    // 步骤:
    // 1. 如果是叶子节点: 加入新点，必要时重新拟合
    // 2. 如果是非叶子节点: 递归到对应子节点
};
```

### 7.3 全局函数

```cpp
// 构建VoxelMap（首帧初始化）
void buildVoxelMap(const vector<pointWithCov> &input_points,
                   const float voxel_size,
                   const int max_layer,
                   const vector<int> &layer_point_size,
                   const int max_points_size,
                   const int max_cov_points_size,
                   const float planer_threshold,
                   unordered_map<VOXEL_LOC, OctoTree*> &feat_map);

// 增量更新VoxelMap（单线程）
void updateVoxelMap(const vector<pointWithCov> &input_points, ...);

// 增量更新VoxelMap（OpenMP并行）
void updateVoxelMapOMP(const vector<pointWithCov> &input_points, ...);

// 搜索最近平面并构建残差列表（OpenMP并行）
void BuildResidualListOMP(const unordered_map<VOXEL_LOC, OctoTree*> &voxel_map,
                          const double voxel_size,
                          const double sigma_num,
                          const int max_layer,
                          const vector<pointWithCov> &pv_list,
                          vector<ptpl> &ptpl_list,
                          vector<Eigen::Vector3d> &non_match_list);

// 计算LiDAR测量点的协方差
void calcLidarCov(const Eigen::Vector3d &point,
                  const double ranging_cov,
                  const double angle_cov,
                  Eigen::Matrix3d &cov);

// 发布VoxelMap可视化（MarkerArray）
void pubVoxelMap(const unordered_map<VOXEL_LOC, OctoTree*> &feat_map,
                 const int publish_max_voxel_layer,
                 const ros::Publisher &pub);
```

---

## 8. VoxelMap++模块：`include/voxelmapplus_util.hpp`

### 8.1 数据结构

```cpp
// 平面结构（VoxelMap++版本）
struct Plane {
    bool is_plane;                    // 是否为有效平面
    bool is_init;                     // 是否已初始化
    int main_direction;               // 主方向: 0/1/2
    M3D plane_cov;                    // 平面参数协方差 (3×3)
    V3D n_vec;                        // 平面法向量
    bool isRootPlane;                 // 是否为根平面（并查集）
    int rgb[3];                       // 可视化颜色

    // 增量计算参数
    double xx, yy, zz;                // 二次项和: Σx², Σy², Σz²
    double xy, xz, yz;                // 交叉项和: Σxy, Σxz, Σyz
    double x, y, z;                   // 一次项和: Σx, Σy, Σz
    V3D center;                       // 平面中心
    Eigen::Matrix3d covariance;       // 点集协方差矩阵
    int points_size;                  // 包含的点数
};

// 并查集节点
class UnionFindNode {
public:
    vector<pointWithCov> temp_points_;  // 点集合
    PlanePtr plane_ptr_;                // 平面指针
    double voxel_center_[3];            // 体素中心
    int all_points_num_;                // 总点数
    int new_points_num_;                // 新增点数

    bool init_node_;                    // 是否已初始化
    bool update_enable_;                // 是否可更新
    bool is_plane;                      // 是否检测到平面
    int id;                             // 节点ID
    UnionFindNode *rootNode;            // 并查集根节点指针

    // 初始化平面
    void InitPlane(const vector<pointWithCov> &points,
                   const PlanePtr &plane,
                   UnionFindNode *node);

    // 初始化并查集节点
    void InitUnionFindNode();

    // 更新平面（添加新点并尝试融合）
    void UpdatePlane(const pointWithCov &pv,
                     VOXEL_LOC &position,
                     unordered_map<VOXEL_LOC, UnionFindNode*> &feat_map);
};
```

### 8.2 全局函数

```cpp
// 构建VoxelMap++（首帧初始化）
void BuildVoxelMap(const vector<pointWithCov> &pv_list,
                   unordered_map<VOXEL_LOC, UnionFindNode*> &feat_map);

// 增量更新VoxelMap++
void UpdateVoxelMap(const vector<pointWithCov> &pv_list,
                    unordered_map<VOXEL_LOC, UnionFindNode*> &feat_map);

// 搜索最近平面并构建残差列表（OpenMP并行）
void BuildResidualListOMP(const unordered_map<VOXEL_LOC, UnionFindNode*> &feat_map,
                          const double voxel_size,
                          const double sigma_num,
                          const vector<pointWithCov> &pv_list,
                          vector<ptpl> &ptpl_list,
                          vector<V3D> &non_match_list);

// 发布VoxelMap++可视化
void pubVoxelMap(const unordered_map<VOXEL_LOC, UnionFindNode*> &feat_map,
                 const ros::Publisher &pub);
```

---

## 9. SO(3)数学模块：`include/so3_math.h`

### 9.1 函数列表

```cpp
// 反对称矩阵
template<typename T>
Eigen::Matrix<T,3,3> skew_sym_mat(const Eigen::Matrix<T,3,1> &v);
// 返回 v 的反对称矩阵 [v]×

// 指数映射（向量 → 旋转矩阵）
template<typename T>
Eigen::Matrix<T,3,3> Exp(const Eigen::Matrix<T,3,1> &&ang);
// 使用Rodrigues公式: R = I + sin(θ)/θ·[φ]× + (1-cos(θ))/θ²·[φ]×²

template<typename T, typename Ts>
Eigen::Matrix<T,3,3> Exp(const Eigen::Matrix<T,3,1> &ang_vel, const Ts &dt);
// 计算: Exp(ω · dt)

template<typename T>
Eigen::Matrix<T,3,3> Exp(const T &v1, const T &v2, const T &v3);
// 计算: Exp([v1, v2, v3]^T)

// 对数映射（旋转矩阵 → 向量）
template<typename T>
Eigen::Matrix<T,3,1> Log(const Eigen::Matrix<T,3,3> &R);
// 返回旋转向量 φ = θ · â

// 旋转矩阵到欧拉角
template<typename T>
Eigen::Matrix<T,3,1> RotMtoEuler(const Eigen::Matrix<T,3,3> &rot);
// 返回 [roll, pitch, yaw]^T (ZYX顺序)
```

---

## 10. iESKF核心模块：`include/IKFoM_toolkit/esekfom/esekfom.hpp`

### 10.1 核心类

```cpp
template<typename state, int process_noise_dof, typename input>
class esekf {
public:
    // 状态和协方差
    state x_;                                    // 当前状态估计
    cov P_;                                      // 状态协方差矩阵

    // 初始化
    void init_dyn_share(state &x, cov &P,
                        processModel f,
                        processMatrix1 f_x,
                        processMatrix2 f_w,
                        share_datastruct &data,
                        int maximum_iter,
                        Scalar limit);

    // 预测（IMU传播）
    void predict(double dt, processNoiseCov Q,
                 const input &i_in);
    // 步骤:
    // 1. 调用 f(x, u) 计算状态导数
    // 2. 使用一阶欧拉法更新状态: x += f(x, u) * dt
    // 3. 计算 F_x 和 F_w Jacobian
    // 4. 传播协方差: P = F_x * P * F_x^T + F_w * Q * F_w^T * dt²

    // 迭代更新（对角观测噪声版本）
    void update_iterated_dyn_share_diagonal();
    // 步骤:
    // 1. 备份当前状态作为预测值
    // 2. for iter = 1 to max_iteration:
    //    a. 调用观测模型计算 h, H, R
    //    b. 计算卡尔曼增益 K
    //    c. 更新状态: x = x ⊞ K * (h - H * (x ⊟ x_pred))
    //    d. 检查收敛
    // 3. 更新协方差: P = (I - K*H) * P
};
```

---

## 11. 流形工具库（MTK）

### 11.1 模块结构

```
IKFoM_toolkit/mtk/
├── build_manifold.hpp     # 流形构建宏 MTK_BUILD_MANIFOLD
├── startIdx.hpp           # 子流形索引管理
├── src/
│   ├── mtkmath.hpp        # 数学运算工具
│   ├── SubManifold.hpp    # 子流形基类
│   └── vectview.hpp       # 向量视图（零拷贝切片）
└── types/
    ├── S2.hpp             # S² 球面流形（重力表示）
    ├── SOn.hpp            # SO(n) 旋转群
    ├── vect.hpp           # ℝⁿ 向量空间
    └── wrapped_cv_mat.hpp # OpenCV矩阵封装
```

### 11.2 关键宏

```cpp
// 构建复合流形
MTK_BUILD_MANIFOLD(name, components...)
// 自动生成:
// - 状态结构体
// - ⊞ (boxplus) 运算: 流形上的加法
// - ⊟ (boxminus) 运算: 流形上的减法
// - DOF 维度计算
// - 子流形索引映射

// S² 流形特性:
// - 2维切空间（嵌入在3维空间中）
// - boxplus: 沿切平面旋转后投影到球面
// - boxminus: 球面上两点的切向量差

// SO(3) 流形特性:
// - 3维切空间
// - boxplus: R' = R * Exp(δφ)
// - boxminus: δφ = Log(R^T * R')
```
