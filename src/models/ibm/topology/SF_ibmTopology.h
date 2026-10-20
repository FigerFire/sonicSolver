#pragma once

/// @file SF_ibmTopology.h
/// @brief IBM 拓扑几何的权威存储：ghost 层、壁面截距、镜像点、法向与壁面速度。
///
/// 该对象只回答 "IBM 几何数据放在哪里"，不参与全局时间推进、通量装配或
/// MPI 通信。索引采用与 Field 一致的 (i,j,k) 含虚胞布局。

#include "SF_valueTypes.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace SF {
namespace IBM {

/// @brief 每点 IBM 壁面/镜像几何 + ghost 层号 + 有符号距离。
///
/// Ownership：由 IBM 方法（WeightBuilder / CompositeIB / mesh loader）构建，
/// 供 ghost-cell / ILW closure / 湍流壁面距离等 reader 消费。Field 不再持有。
class IBMGeometry {
public:
    /// @brief 按含虚胞网格尺寸分配所有数组并初始化为无几何状态。
    void setup(int mx, int my, int mz) {
        mx_ = mx;
        my_ = my;
        const size_t total = static_cast<size_t>(mx) * my * mz;
        ghostLayer_.assign(total, 0);
        signedDistance_.assign(total, 0.0);
        wallX_.assign(total, 0.0);
        wallY_.assign(total, 0.0);
        wallZ_.assign(total, 0.0);
        imageX_.assign(total, 0.0);
        imageY_.assign(total, 0.0);
        imageZ_.assign(total, 0.0);
        normalX_.assign(total, 0.0);
        normalY_.assign(total, 0.0);
        normalZ_.assign(total, 0.0);
        wallVelocityX_.assign(total, 0.0);
        wallVelocityY_.assign(total, 0.0);
        wallVelocityZ_.assign(total, 0.0);
    }

    /// @brief 重置全部几何为 "无几何"，保持已分配容量。
    void clear() {
        std::fill(ghostLayer_.begin(), ghostLayer_.end(), 0);
        std::fill(signedDistance_.begin(), signedDistance_.end(), 0.0);
        std::fill(wallX_.begin(), wallX_.end(), 0.0);
        std::fill(wallY_.begin(), wallY_.end(), 0.0);
        std::fill(wallZ_.begin(), wallZ_.end(), 0.0);
        std::fill(imageX_.begin(), imageX_.end(), 0.0);
        std::fill(imageY_.begin(), imageY_.end(), 0.0);
        std::fill(imageZ_.begin(), imageZ_.end(), 0.0);
        std::fill(normalX_.begin(), normalX_.end(), 0.0);
        std::fill(normalY_.begin(), normalY_.end(), 0.0);
        std::fill(normalZ_.begin(), normalZ_.end(), 0.0);
        std::fill(wallVelocityX_.begin(), wallVelocityX_.end(), 0.0);
        std::fill(wallVelocityY_.begin(), wallVelocityY_.end(), 0.0);
        std::fill(wallVelocityZ_.begin(), wallVelocityZ_.end(), 0.0);
    }

    /// @brief 是否已分配。
    bool allocated() const { return !ghostLayer_.empty(); }

    /// @brief 平铺索引，与 Field::getIdx 一致。
    inline size_t idx(int i, int j, int k) const {
        return (static_cast<size_t>(k) * my_ + j) * mx_ + i;
    }

    /// @brief IBM ghost 层号；非 ghost 点为 0。
    inline int& ghostLayer(int i, int j, int k) {
        return ghostLayer_[idx(i, j, k)];
    }
    inline int ghostLayer(int i, int j, int k) const {
        return ghostLayer_[idx(i, j, k)];
    }

    /// @brief 到 IBM 壁面的有符号距离。
    inline double& signedDistance(int i, int j, int k) {
        return signedDistance_[idx(i, j, k)];
    }
    inline double signedDistance(int i, int j, int k) const {
        return signedDistance_[idx(i, j, k)];
    }

    /// @brief 写入壁面截距、镜像点、法向与壁面速度。
    inline void setGeometry(int i, int j, int k,
                            const Vector3& wallPoint,
                            const Vector3& imagePoint,
                            const Vector3& wallNormal,
                            const Vector3& wallVelocity = Vector3()) {
        const size_t id = idx(i, j, k);
        wallX_[id] = wallPoint.x;
        wallY_[id] = wallPoint.y;
        wallZ_[id] = wallPoint.z;
        imageX_[id] = imagePoint.x;
        imageY_[id] = imagePoint.y;
        imageZ_[id] = imagePoint.z;
        normalX_[id] = wallNormal.x;
        normalY_[id] = wallNormal.y;
        normalZ_[id] = wallNormal.z;
        wallVelocityX_[id] = wallVelocity.x;
        wallVelocityY_[id] = wallVelocity.y;
        wallVelocityZ_[id] = wallVelocity.z;
    }

    /// @brief 该点是否已写入有效壁面几何（法向模长非零）。
    inline bool hasGeometry(int i, int j, int k) const {
        const size_t id = idx(i, j, k);
        const double n2 = normalX_[id] * normalX_[id]
                        + normalY_[id] * normalY_[id]
                        + normalZ_[id] * normalZ_[id];
        return n2 > 1e-24;
    }

    inline Vector3 wallPoint(int i, int j, int k) const {
        const size_t id = idx(i, j, k);
        return Vector3(wallX_[id], wallY_[id], wallZ_[id]);
    }
    inline Vector3 imagePoint(int i, int j, int k) const {
        const size_t id = idx(i, j, k);
        return Vector3(imageX_[id], imageY_[id], imageZ_[id]);
    }
    inline Vector3 wallNormal(int i, int j, int k) const {
        const size_t id = idx(i, j, k);
        return Vector3(normalX_[id], normalY_[id], normalZ_[id]);
    }
    inline Vector3 wallVelocity(int i, int j, int k) const {
        const size_t id = idx(i, j, k);
        return Vector3(wallVelocityX_[id], wallVelocityY_[id],
                       wallVelocityZ_[id]);
    }

private:
    int mx_ = 0;
    int my_ = 0;
    std::vector<int> ghostLayer_;
    std::vector<double> signedDistance_;
    std::vector<double> wallX_, wallY_, wallZ_;
    std::vector<double> imageX_, imageY_, imageZ_;
    std::vector<double> normalX_, normalY_, normalZ_;
    std::vector<double> wallVelocityX_, wallVelocityY_, wallVelocityZ_;
};

} // namespace IBM
} // namespace SF
