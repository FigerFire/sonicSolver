#pragma once

/// @file SF_couplingGraph.h
/// @brief 建立唯一 marker--Eulerian DOF 边及体积加权插值行。

#include "SF_field.h"
#include "SF_immersedConstraint.h"
#include "topology/SF_markerRegistry.h"

#include <vector>

namespace SF::IBM::Topology {

/// @brief 不含通信的 canonical marker--Eulerian 耦合图。
/// @details build 只生成 Eulerian owner 的局部 raw edge；normalize 接收 Runtime 已归并的
/// marker 权重和，因而模型不需要知道 rank、MPI 或通信计划。
class CouplingGraph {
public:
    /// @brief 为每个 marker 建立唯一的本地 Eulerian owner raw edge。
    void build(const Field& field,
               const std::vector<SurfaceMarker>& markers,
               double supportRadius);

    /// @brief 用所有 owner edge 的 canonical 权重和生成分区一致 J 行。
    void normalize(const std::vector<double>& globalNormalizations);

    /// @brief 返回本 rank owner edge 对每个 marker 的 raw 权重和。
    const std::vector<double>& rawNormalizations() const {
        return rawNormalizations_;
    }

    const std::vector<std::vector<FDM::ImmersedInterpolationWeight>>& rows()
        const { return rows_; }

private:
    std::vector<std::vector<FDM::ImmersedInterpolationWeight>> rawRows_;
    std::vector<std::vector<FDM::ImmersedInterpolationWeight>> rows_;
    std::vector<double> rawNormalizations_;
};

} // namespace SF::IBM::Topology
