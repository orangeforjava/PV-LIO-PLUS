#ifndef R_VOXEL_MAP_UTIL_HPP
#define R_VOXEL_MAP_UTIL_HPP

#include "voxel_map_util.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <stack>
#include <unordered_set>

namespace r_voxel_map_ns
{
    using voxel_map_ns::Plane;
    using voxel_map_ns::VOXEL_LOC;
    using voxel_map_ns::calcLidarCov;
    using voxel_map_ns::pointWithCov;
    using voxel_map_ns::ptpl;

    static int plane_id = 0;
    static double ransac_distance_threshold = 0.10;
    static double inlier_ratio_threshold    = 0.55;
    static int ransac_iterations            = 40;
    static int min_points_threshold         = 5;
    static int validity_grid_divider        = 4;
    static int rebuild_point_threshold      = 24;

    struct GridCoord
    {
        int x = 0;
        int y = 0;

        bool operator==(const GridCoord &other) const
        {
            return x == other.x && y == other.y;
        }
    };
}

namespace std
{
    template <>
    struct hash<r_voxel_map_ns::GridCoord>
    {
        size_t operator()(const r_voxel_map_ns::GridCoord &coord) const
        {
            return (static_cast<size_t>(coord.x) << 32) ^ static_cast<size_t>(coord.y);
        }
    };
}  // namespace std

namespace r_voxel_map_ns
{
    class OctoTree
    {
    public:
        std::vector<pointWithCov> all_points_;
        std::vector<pointWithCov> plane_points_;
        std::vector<pointWithCov> pending_outliers_;
        Plane *plane_ptr_;
        int max_layer_;
        int layer_;
        OctoTree *leaves_[8];
        double voxel_center_[3];
        std::vector<int> layer_point_size_;
        float quater_length_;
        int max_points_size_;
        int max_cov_points_size_;
        float planer_threshold_;
        double voxel_size_;
        bool is_root_;
        int new_points_since_rebuild_;
        bool init_octo_;

        OctoTree(int max_layer, int layer, std::vector<int> layer_point_size,
                 int max_point_size, int max_cov_points_size, float planer_threshold,
                 float voxel_size, bool is_root = false)
            : plane_ptr_(new Plane), max_layer_(max_layer), layer_(layer), layer_point_size_(std::move(layer_point_size)), quater_length_(0.0F),
              max_points_size_(max_point_size), max_cov_points_size_(max_cov_points_size), planer_threshold_(planer_threshold),
              voxel_size_(voxel_size), is_root_(is_root), new_points_since_rebuild_(0), init_octo_(false)
        {
            voxel_center_[0] = voxel_center_[1] = voxel_center_[2] = 0.0;
            for (int i = 0; i < 8; ++i)
            {
                leaves_[i] = nullptr;
            }
        }

        ~OctoTree()
        {
            clear_children();
            delete plane_ptr_;
        }

        void BuildRecursive(const std::vector<pointWithCov> &points)
        {
            clear_children();
            *plane_ptr_               = Plane();
            plane_ptr_->update_enable = true;
            plane_points_.clear();
            pending_outliers_.clear();
            all_points_       = points;
            init_octo_        = true;
            new_points_since_rebuild_ = 0;

            if (points.size() < static_cast<size_t>(std::max(min_points_threshold, layer_point_threshold())))
            {
                return;
            }

            std::vector<pointWithCov> inliers;
            std::vector<pointWithCov> outliers;
            const bool has_candidate = ransac_plane(points, inliers, outliers);
            if (!has_candidate || inliers.size() <= inlier_ratio_threshold * points.size())
            {
                pending_outliers_ = points;
                subdivide_outliers(points);
                return;
            }

            Plane candidate_plane;
            init_plane_from_points(inliers, &candidate_plane);
            if (!candidate_plane.is_plane)
            {
                pending_outliers_ = points;
                subdivide_outliers(points);
                return;
            }

            std::vector<pointWithCov> valid_inliers;
            std::vector<pointWithCov> invalidated_points;
            if (!plane_validity_check(candidate_plane, points.size(), inliers, valid_inliers, invalidated_points))
            {
                pending_outliers_ = points;
                subdivide_outliers(points);
                return;
            }

            plane_points_ = valid_inliers;
            init_plane_from_points(plane_points_, plane_ptr_);
            plane_ptr_->is_update = true;

            pending_outliers_ = outliers;
            pending_outliers_.insert(pending_outliers_.end(), invalidated_points.begin(), invalidated_points.end());
            subdivide_outliers(pending_outliers_);
        }

        bool UpdatePoint(const pointWithCov &pv, const double sigma_num)
        {
            all_points_.push_back(pv);
            ++new_points_since_rebuild_;

            OctoTree *best_node = nullptr;
            ptpl best_ptpl;
            double best_prob = 0.0;
            find_best_match(pv, sigma_num, best_node, best_ptpl, best_prob);

            if (best_node != nullptr)
            {
                best_node->plane_points_.push_back(pv);
                best_node->recompute_plane_from_inliers();
                if (!best_node->plane_ptr_->is_plane)
                {
                    rebuild_from_all_points();
                    return false;
                }
            }
            else
            {
                pending_outliers_.push_back(pv);
            }

            const bool enough_new_points = new_points_since_rebuild_ >= std::max(rebuild_point_threshold, layer_point_threshold() * 2);
            const bool enough_outliers   = pending_outliers_.size() >= static_cast<size_t>(std::max(min_points_threshold, layer_point_threshold()));
            if (enough_new_points || enough_outliers)
            {
                rebuild_from_all_points();
            }

            return best_node != nullptr;
        }

        void FindBestMatch(const pointWithCov &pv, const double sigma_num,
                           bool &is_success, double &prob, ptpl &single_ptpl) const
        {
            is_success = false;
            prob       = 0.0;
            find_best_match_const(pv, sigma_num, single_ptpl, prob, is_success);
        }

        void CollectUpdatePlanes(const int pub_max_voxel_layer,
                                 std::vector<Plane> &plane_list) const
        {
            if (layer_ > pub_max_voxel_layer)
            {
                return;
            }

            if (plane_ptr_->is_update)
            {
                plane_list.push_back(*plane_ptr_);
            }

            for (const OctoTree *leaf : leaves_)
            {
                if (leaf != nullptr)
                {
                    leaf->CollectUpdatePlanes(pub_max_voxel_layer, plane_list);
                }
            }
        }

        void ResetUpdateFlag()
        {
            plane_ptr_->is_update = false;
            for (OctoTree *leaf : leaves_)
            {
                if (leaf != nullptr)
                {
                    leaf->ResetUpdateFlag();
                }
            }
        }

    private:
        int layer_point_threshold() const
        {
            if (layer_ >= 0 && layer_ < static_cast<int>(layer_point_size_.size()))
            {
                return std::max(min_points_threshold, layer_point_size_[layer_]);
            }
            return min_points_threshold;
        }

        void clear_children()
        {
            for (int i = 0; i < 8; ++i)
            {
                if (leaves_[i] != nullptr)
                {
                    delete leaves_[i];
                    leaves_[i] = nullptr;
                }
            }
        }

        static bool fit_plane_from_three_points(const pointWithCov &p1, const pointWithCov &p2,
                                                const pointWithCov &p3, Eigen::Vector3d &normal, double &d)
        {
            Eigen::Vector3d v1 = p2.point_world - p1.point_world;
            Eigen::Vector3d v2 = p3.point_world - p1.point_world;
            normal             = v1.cross(v2);
            const double norm  = normal.norm();
            if (norm < 1e-6)
            {
                return false;
            }
            normal.normalize();
            d = -normal.dot(p1.point_world);
            return std::isfinite(d);
        }

        bool ransac_plane(const std::vector<pointWithCov> &points,
                          std::vector<pointWithCov> &inliers,
                          std::vector<pointWithCov> &outliers) const
        {
            inliers.clear();
            outliers.clear();
            if (points.size() < 3)
            {
                outliers = points;
                return false;
            }

            const uint32_t seed = static_cast<uint32_t>(layer_ + 1) * 2654435761U
                                  ^ static_cast<uint32_t>(points.size() * 97)
                                  ^ static_cast<uint32_t>((std::llround(voxel_center_[0] * 10) & 0xffff) << 16)
                                  ^ static_cast<uint32_t>(std::llround(voxel_center_[1] * 10) & 0xffff)
                                  ^ static_cast<uint32_t>(std::llround(voxel_center_[2] * 10) & 0xffff);
            std::mt19937 rng(seed);
            std::uniform_int_distribution<size_t> dist(0, points.size() - 1);

            std::vector<size_t> best_indices;
            best_indices.reserve(points.size());

            for (int iter = 0; iter < ransac_iterations; ++iter)
            {
                size_t i1 = dist(rng);
                size_t i2 = dist(rng);
                size_t i3 = dist(rng);
                if (i1 == i2 || i1 == i3 || i2 == i3)
                {
                    continue;
                }

                Eigen::Vector3d normal;
                double d = 0.0;
                if (!fit_plane_from_three_points(points[i1], points[i2], points[i3], normal, d))
                {
                    continue;
                }

                std::vector<size_t> current_indices;
                current_indices.reserve(points.size());
                for (size_t i = 0; i < points.size(); ++i)
                {
                    const double distance = std::fabs(normal.dot(points[i].point_world) + d);
                    if (distance <= ransac_distance_threshold)
                    {
                        current_indices.push_back(i);
                    }
                }

                if (current_indices.size() > best_indices.size())
                {
                    best_indices = current_indices;
                }
            }

            if (best_indices.empty())
            {
                outliers = points;
                return false;
            }

            std::vector<bool> is_inlier(points.size(), false);
            for (size_t index : best_indices)
            {
                is_inlier[index] = true;
                inliers.push_back(points[index]);
            }
            for (size_t i = 0; i < points.size(); ++i)
            {
                if (!is_inlier[i])
                {
                    outliers.push_back(points[i]);
                }
            }
            return true;
        }

        void init_plane_from_points(const std::vector<pointWithCov> &points, Plane *plane) const
        {
            plane->plane_cov   = Eigen::Matrix<double, 6, 6>::Zero();
            plane->covariance  = Eigen::Matrix3d::Zero();
            plane->center      = Eigen::Vector3d::Zero();
            plane->normal      = Eigen::Vector3d::Zero();
            plane->points_size = static_cast<int>(points.size());
            plane->radius      = 0.0F;

            if (points.size() < 3)
            {
                plane->is_plane = false;
                plane->is_init  = true;
                return;
            }

            for (const auto &pv : points)
            {
                plane->covariance += pv.point_world * pv.point_world.transpose();
                plane->center += pv.point_world;
            }
            plane->center /= plane->points_size;
            plane->covariance /= plane->points_size;
            plane->covariance -= plane->center * plane->center.transpose();

            Eigen::EigenSolver<Eigen::Matrix3d> es(plane->covariance);
            Eigen::Matrix3cd evecs    = es.eigenvectors();
            Eigen::Vector3cd evals    = es.eigenvalues();
            Eigen::Vector3d evalsReal = evals.real();
            Eigen::Matrix3f::Index evalsMin, evalsMax;
            evalsReal.minCoeff(&evalsMin);
            evalsReal.maxCoeff(&evalsMax);
            const int evalsMid = 3 - static_cast<int>(evalsMin) - static_cast<int>(evalsMax);

            if (!plane->is_init)
            {
                plane->id      = plane_id++;
                plane->is_init = true;
            }

            if (plane->last_update_points_size == 0 || plane->points_size - plane->last_update_points_size > 20)
            {
                plane->last_update_points_size = plane->points_size;
                plane->is_update               = true;
            }

            plane->normal          = evecs.real().col(evalsMin);
            plane->y_normal        = evecs.real().col(evalsMid);
            plane->x_normal        = evecs.real().col(evalsMax);
            plane->min_eigen_value = evalsReal(evalsMin);
            plane->mid_eigen_value = evalsReal(evalsMid);
            plane->max_eigen_value = evalsReal(evalsMax);
            plane->radius          = std::sqrt(std::max(0.0, evalsReal(evalsMax)));
            plane->d               = -(plane->normal.dot(plane->center));

            if (evalsReal(evalsMin) < planer_threshold_)
            {
                Eigen::Matrix3d J_Q;
                J_Q << 1.0 / plane->points_size, 0, 0, 0, 1.0 / plane->points_size, 0, 0, 0, 1.0 / plane->points_size;

                for (size_t i = 0; i < points.size(); ++i)
                {
                    Eigen::Matrix<double, 6, 3> J;
                    Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
                    for (int m = 0; m < 3; ++m)
                    {
                        if (m != static_cast<int>(evalsMin))
                        {
                            const double denom = (plane->points_size) * (evalsReal[evalsMin] - evalsReal[m]);
                            if (std::fabs(denom) < 1e-12)
                            {
                                continue;
                            }
                            F.row(m) = (points[i].point_world - plane->center).transpose()
                                       / denom
                                       * (evecs.real().col(m) * evecs.real().col(evalsMin).transpose()
                                          + evecs.real().col(evalsMin) * evecs.real().col(m).transpose());
                        }
                    }
                    J.block<3, 3>(0, 0) = evecs.real() * F;
                    J.block<3, 3>(3, 0) = J_Q;
                    plane->plane_cov += J * points[i].cov_world * J.transpose();
                }

                plane->is_plane = true;
            }
            else
            {
                plane->is_plane = false;
            }
        }

        bool plane_validity_check(const Plane &plane, const size_t original_point_count,
                                  const std::vector<pointWithCov> &inliers,
                                  std::vector<pointWithCov> &valid_inliers,
                                  std::vector<pointWithCov> &invalid_points) const
        {
            valid_inliers.clear();
            invalid_points.clear();
            if (inliers.empty())
            {
                return false;
            }

            const double grid_resolution = voxel_size_ / std::max(1, validity_grid_divider);
            if (grid_resolution <= 0.0)
            {
                valid_inliers = inliers;
                return true;
            }

            struct CellData
            {
                std::vector<pointWithCov> points;
                int count = 0;
            };

            std::unordered_map<GridCoord, CellData> grid_map;
            grid_map.reserve(inliers.size());
            for (const auto &point : inliers)
            {
                const Eigen::Vector3d delta = point.point_world - plane.center;
                const int grid_x            = static_cast<int>(std::floor(delta.dot(plane.x_normal) / grid_resolution));
                const int grid_y            = static_cast<int>(std::floor(delta.dot(plane.y_normal) / grid_resolution));
                GridCoord coord{grid_x, grid_y};
                CellData &cell = grid_map[coord];
                cell.points.push_back(point);
                ++cell.count;
            }

            std::unordered_set<GridCoord> visited;
            visited.reserve(grid_map.size());
            int best_cluster_points = 0;
            std::vector<GridCoord> best_cluster;

            for (const auto &entry : grid_map)
            {
                const GridCoord &start = entry.first;
                if (visited.find(start) != visited.end())
                {
                    continue;
                }

                std::stack<GridCoord> stack;
                stack.push(start);
                visited.insert(start);

                std::vector<GridCoord> cluster;
                int cluster_points = 0;

                while (!stack.empty())
                {
                    GridCoord current = stack.top();
                    stack.pop();
                    cluster.push_back(current);
                    cluster_points += grid_map[current].count;

                    const std::array<GridCoord, 4> neighbors = {
                        GridCoord{current.x + 1, current.y},
                        GridCoord{current.x - 1, current.y},
                        GridCoord{current.x, current.y + 1},
                        GridCoord{current.x, current.y - 1}};
                    for (const auto &neighbor : neighbors)
                    {
                        if (grid_map.find(neighbor) == grid_map.end())
                        {
                            continue;
                        }
                        if (visited.insert(neighbor).second)
                        {
                            stack.push(neighbor);
                        }
                    }
                }

                if (cluster_points > best_cluster_points)
                {
                    best_cluster_points = cluster_points;
                    best_cluster        = cluster;
                }
            }

            if (best_cluster.empty() || best_cluster_points <= inlier_ratio_threshold * original_point_count)
            {
                invalid_points = inliers;
                return false;
            }

            std::unordered_set<GridCoord> best_cells(best_cluster.begin(), best_cluster.end());
            for (const auto &entry : grid_map)
            {
                const bool is_best = best_cells.find(entry.first) != best_cells.end();
                if (is_best)
                {
                    valid_inliers.insert(valid_inliers.end(), entry.second.points.begin(), entry.second.points.end());
                }
                else
                {
                    invalid_points.insert(invalid_points.end(), entry.second.points.begin(), entry.second.points.end());
                }
            }

            return !valid_inliers.empty();
        }

        int child_index(const Eigen::Vector3d &point) const
        {
            int xyz[3] = {0, 0, 0};
            if (point[0] > voxel_center_[0])
            {
                xyz[0] = 1;
            }
            if (point[1] > voxel_center_[1])
            {
                xyz[1] = 1;
            }
            if (point[2] > voxel_center_[2])
            {
                xyz[2] = 1;
            }
            return 4 * xyz[0] + 2 * xyz[1] + xyz[2];
        }

        void subdivide_outliers(const std::vector<pointWithCov> &outliers)
        {
            if (outliers.empty() || layer_ >= max_layer_)
            {
                return;
            }

            std::array<std::vector<pointWithCov>, 8> child_points;
            for (const auto &point : outliers)
            {
                child_points[child_index(point.point_world)].push_back(point);
            }

            for (int i = 0; i < 8; ++i)
            {
                if (child_points[i].size() < static_cast<size_t>(min_points_threshold))
                {
                    continue;
                }

                int xyz[3] = {i / 4, (i % 4) / 2, i % 2};
                leaves_[i] = new OctoTree(max_layer_, layer_ + 1, layer_point_size_,
                                          max_points_size_, max_cov_points_size_,
                                          planer_threshold_, voxel_size_, false);
                leaves_[i]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
                leaves_[i]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
                leaves_[i]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
                leaves_[i]->quater_length_   = quater_length_ / 2.0F;
                leaves_[i]->BuildRecursive(child_points[i]);
            }
        }

        void recompute_plane_from_inliers()
        {
            if (plane_points_.empty())
            {
                plane_ptr_->is_plane = false;
                return;
            }

            if (plane_points_.size() > static_cast<size_t>(max_points_size_))
            {
                plane_points_.erase(plane_points_.begin(),
                                    plane_points_.begin() + (plane_points_.size() - max_points_size_));
            }
            init_plane_from_points(plane_points_, plane_ptr_);
            plane_ptr_->is_update = true;
        }

        void rebuild_from_all_points()
        {
            const std::vector<pointWithCov> all_points = all_points_;
            BuildRecursive(all_points);
        }

        static bool evaluate_plane_candidate(const pointWithCov &pv, const Plane &plane,
                                             const int current_layer, const double sigma_num,
                                             ptpl &candidate_ptpl, double &candidate_prob)
        {
            const double radius_k                 = 3.0;
            const Eigen::Vector3d &point_world    = pv.point_world;
            const double distance_to_plane        = std::fabs(point_world.dot(plane.normal) + plane.d);
            const double distance_to_center_sq    = (point_world - plane.center).squaredNorm();
            const double projected_distance_sq    = distance_to_center_sq - distance_to_plane * distance_to_plane;
            if (projected_distance_sq < 0.0)
            {
                return false;
            }

            const double range_distance = std::sqrt(projected_distance_sq);
            if (range_distance > radius_k * plane.radius)
            {
                return false;
            }

            Eigen::Matrix<double, 1, 6> J_nq;
            J_nq.block<1, 3>(0, 0) = point_world - plane.center;
            J_nq.block<1, 3>(0, 3) = -plane.normal;
            double sigma_l         = J_nq * plane.plane_cov * J_nq.transpose();
            sigma_l += plane.normal.transpose() * pv.cov_world * plane.normal;
            if (sigma_l <= 1e-12)
            {
                return false;
            }

            if (distance_to_plane >= sigma_num * std::sqrt(sigma_l))
            {
                return false;
            }

            constexpr double kTwoPi = 6.28318530717958647692;
            candidate_prob = 1.0 / std::sqrt(kTwoPi * sigma_l)
                             * std::exp(-0.5 * distance_to_plane * distance_to_plane / sigma_l);
            candidate_ptpl.point       = pv.point_lidar;
            candidate_ptpl.point_world = pv.point_world;
            candidate_ptpl.plane_cov   = plane.plane_cov;
            candidate_ptpl.normal      = plane.normal;
            candidate_ptpl.center      = plane.center;
            candidate_ptpl.d           = plane.d;
            candidate_ptpl.layer       = current_layer;
            candidate_ptpl.cov_lidar   = pv.cov_lidar;
            candidate_ptpl.cov_world   = pv.cov_world;
            return true;
        }

        void find_best_match(const pointWithCov &pv, const double sigma_num,
                             OctoTree *&best_node, ptpl &best_ptpl, double &best_prob)
        {
            if (plane_ptr_->is_plane)
            {
                ptpl candidate_ptpl;
                double candidate_prob = 0.0;
                if (evaluate_plane_candidate(pv, *plane_ptr_, layer_, sigma_num, candidate_ptpl, candidate_prob)
                    && candidate_prob > best_prob)
                {
                    best_prob = candidate_prob;
                    best_node = this;
                    best_ptpl = candidate_ptpl;
                }
            }

            for (OctoTree *leaf : leaves_)
            {
                if (leaf != nullptr)
                {
                    leaf->find_best_match(pv, sigma_num, best_node, best_ptpl, best_prob);
                }
            }
        }

        void find_best_match_const(const pointWithCov &pv, const double sigma_num,
                                   ptpl &best_ptpl, double &best_prob,
                                   bool &is_success) const
        {
            if (plane_ptr_->is_plane)
            {
                ptpl candidate_ptpl;
                double candidate_prob = 0.0;
                if (evaluate_plane_candidate(pv, *plane_ptr_, layer_, sigma_num, candidate_ptpl, candidate_prob)
                    && candidate_prob > best_prob)
                {
                    best_prob   = candidate_prob;
                    best_ptpl   = candidate_ptpl;
                    is_success  = true;
                }
            }

            for (const OctoTree *leaf : leaves_)
            {
                if (leaf != nullptr)
                {
                    leaf->find_best_match_const(pv, sigma_num, best_ptpl, best_prob, is_success);
                }
            }
        }
    };

    inline VOXEL_LOC ComputeVoxelPosition(const pointWithCov &point, const float voxel_size)
    {
        float loc_xyz[3];
        for (int j = 0; j < 3; ++j)
        {
            loc_xyz[j] = point.point_world[j] / voxel_size;
            if (loc_xyz[j] < 0)
            {
                loc_xyz[j] -= 1.0;
            }
        }
        return VOXEL_LOC(static_cast<int64_t>(loc_xyz[0]),
                         static_cast<int64_t>(loc_xyz[1]),
                         static_cast<int64_t>(loc_xyz[2]));
    }

    inline OctoTree *CreateRootTree(const VOXEL_LOC &position, const float voxel_size,
                                    const int max_layer, const std::vector<int> &layer_point_size,
                                    const int max_points_size, const int max_cov_points_size,
                                    const float planer_threshold)
    {
        OctoTree *octo_tree                  = new OctoTree(max_layer, 0, layer_point_size, max_points_size,
                                                            max_cov_points_size, planer_threshold, voxel_size, true);
        octo_tree->quater_length_            = voxel_size / 4.0F;
        octo_tree->voxel_center_[0]          = (0.5 + position.x) * voxel_size;
        octo_tree->voxel_center_[1]          = (0.5 + position.y) * voxel_size;
        octo_tree->voxel_center_[2]          = (0.5 + position.z) * voxel_size;
        return octo_tree;
    }

    inline void buildVoxelMap(const std::vector<pointWithCov> &input_points,
                              const float voxel_size, const int max_layer,
                              const std::vector<int> &layer_point_size,
                              const int max_points_size, const int max_cov_points_size,
                              const float planer_threshold,
                              std::unordered_map<VOXEL_LOC, OctoTree *> &feat_map)
    {
        std::unordered_map<VOXEL_LOC, std::vector<pointWithCov>> grouped_points;
        for (const auto &point : input_points)
        {
            grouped_points[ComputeVoxelPosition(point, voxel_size)].push_back(point);
        }

        for (auto &entry : grouped_points)
        {
            OctoTree *octo_tree = CreateRootTree(entry.first, voxel_size, max_layer, layer_point_size,
                                                 max_points_size, max_cov_points_size, planer_threshold);
            feat_map[entry.first] = octo_tree;
            octo_tree->BuildRecursive(entry.second);
        }
    }

    inline void updateVoxelMapOMP(const std::vector<pointWithCov> &input_points,
                                  const float voxel_size, const int max_layer,
                                  const std::vector<int> &layer_point_size,
                                  const int max_points_size, const int max_cov_points_size,
                                  const float planer_threshold,
                                  std::unordered_map<VOXEL_LOC, OctoTree *> &feat_map)
    {
        std::unordered_map<VOXEL_LOC, std::vector<pointWithCov>> grouped_points;
        for (const auto &point : input_points)
        {
            grouped_points[ComputeVoxelPosition(point, voxel_size)].push_back(point);
        }

        for (auto &entry : grouped_points)
        {
            auto iter = feat_map.find(entry.first);
            if (iter == feat_map.end())
            {
                OctoTree *octo_tree = CreateRootTree(entry.first, voxel_size, max_layer, layer_point_size,
                                                     max_points_size, max_cov_points_size, planer_threshold);
                feat_map[entry.first] = octo_tree;
                octo_tree->BuildRecursive(entry.second);
                continue;
            }

            for (const auto &point : entry.second)
            {
                iter->second->UpdatePoint(point, 3.0);
            }
        }
    }

    inline void BuildResidualListOMP(const std::unordered_map<VOXEL_LOC, OctoTree *> &voxel_map,
                                     const double voxel_size, const double sigma_num,
                                     const std::vector<pointWithCov> &pv_list,
                                     std::vector<ptpl> &ptpl_list,
                                     std::vector<Eigen::Vector3d> &non_match)
    {
        ptpl_list.clear();
        non_match.clear();
        ptpl_list.reserve(pv_list.size());
        std::vector<ptpl> all_ptpl_list(pv_list.size());
        std::vector<bool> useful_ptpl(pv_list.size(), false);

#ifdef MP_EN
        omp_set_num_threads(MP_PROC_NUM);
#pragma omp parallel for
#endif
        for (size_t i = 0; i < pv_list.size(); ++i)
        {
            const pointWithCov &pv = pv_list[i];
            const VOXEL_LOC position = ComputeVoxelPosition(pv, voxel_size);
            auto iter = voxel_map.find(position);
            if (iter == voxel_map.end())
            {
                continue;
            }

            ptpl single_ptpl;
            bool is_success = false;
            double prob     = 0.0;
            iter->second->FindBestMatch(pv, sigma_num, is_success, prob, single_ptpl);

            if (!is_success)
            {
                VOXEL_LOC near_position = position;
                if (pv.point_world[0] > iter->second->voxel_center_[0] + iter->second->quater_length_)
                {
                    near_position.x += 1;
                }
                else if (pv.point_world[0] < iter->second->voxel_center_[0] - iter->second->quater_length_)
                {
                    near_position.x -= 1;
                }
                if (pv.point_world[1] > iter->second->voxel_center_[1] + iter->second->quater_length_)
                {
                    near_position.y += 1;
                }
                else if (pv.point_world[1] < iter->second->voxel_center_[1] - iter->second->quater_length_)
                {
                    near_position.y -= 1;
                }
                if (pv.point_world[2] > iter->second->voxel_center_[2] + iter->second->quater_length_)
                {
                    near_position.z += 1;
                }
                else if (pv.point_world[2] < iter->second->voxel_center_[2] - iter->second->quater_length_)
                {
                    near_position.z -= 1;
                }

                auto near_iter = voxel_map.find(near_position);
                if (near_iter != voxel_map.end())
                {
                    near_iter->second->FindBestMatch(pv, sigma_num, is_success, prob, single_ptpl);
                }
            }

            if (is_success)
            {
                useful_ptpl[i]   = true;
                all_ptpl_list[i] = single_ptpl;
            }
        }

        for (size_t i = 0; i < useful_ptpl.size(); ++i)
        {
            if (useful_ptpl[i])
            {
                ptpl_list.push_back(all_ptpl_list[i]);
            }
            else
            {
                non_match.push_back(pv_list[i].point_world);
            }
        }
    }

    inline void pubVoxelMap(const std::unordered_map<VOXEL_LOC, OctoTree *> &voxel_map,
                            const int pub_max_voxel_layer, const ros::Publisher &plane_map_pub)
    {
        double max_trace = 0.25;
        double pow_num   = 0.2;
        ros::Rate loop(500);
        float use_alpha = 0.8F;
        visualization_msgs::MarkerArray voxel_plane;
        voxel_plane.markers.reserve(1000000);
        std::vector<Plane> pub_plane_list;
        for (const auto &entry : voxel_map)
        {
            entry.second->CollectUpdatePlanes(pub_max_voxel_layer, pub_plane_list);
        }

        for (const auto &plane : pub_plane_list)
        {
            V3D plane_cov = plane.plane_cov.block<3, 3>(0, 0).diagonal();
            double trace  = std::min(max_trace, plane_cov.sum());
            trace         = std::pow(trace * (1.0 / max_trace), pow_num);
            uint8_t r, g, b;
            voxel_map_ns::mapJet(trace, 0, 1, r, g, b);
            Eigen::Vector3d plane_rgb(r / 256.0, g / 256.0, b / 256.0);
            voxel_map_ns::pubSinglePlane(voxel_plane, "plane", plane, plane.is_plane ? use_alpha : 0.0F, plane_rgb);
        }

        plane_map_pub.publish(voxel_plane);
        loop.sleep();

        for (const auto &entry : voxel_map)
        {
            entry.second->ResetUpdateFlag();
        }
    }
}

#endif
