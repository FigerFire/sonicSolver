#pragma once

/// @file SF_surfaceOperator.h
/// @brief 表面 Lambda_s 的 Lagrangian 积分点、J 插值和伴随传播算子。

#include "SF_field.h"
#include "constraint/SF_constraintOperator.h"
#include "geoProcessing/SF_STLGeometry.h"
#include "topology/SF_couplingGraph.h"
#include "topology/SF_markerRegistry.h"

#include <vector>

namespace SF::IBM::Forcing {

/// @brief 基于 STL 三角形重心积分和紧支撑 Wendland C2 核的表面算子。
class SurfaceOperator final : public ISurfaceConstraintOperator {
public:
    /// @brief 为当前刚体构型构造表面点与 local owner raw edge；须后续归一化。
    /// @param field Eulerian 结构网格。
    /// @param triangles 已变换到当前构型的三角形。
    /// @param center 当前质心，单位 m。
    /// @param supportRadius 核支撑半径，单位 m。
    const FDM::ImmersedSurfaceSystem& build(
        const Field& field,
        const std::vector<GeoProcessing::Triangle>& triangles,
        const Vector3& center,
        double supportRadius) override;

    /// @brief 应用 Runtime 已归并的 marker 权重和，刷新 system 内的 J 行。
    const FDM::ImmersedSurfaceSystem& normalizeDistributed(
        const std::vector<double>& normalizations) override;
    const std::vector<double>& localNormalizations() const override {
        return couplingGraph_.rawNormalizations();
    }

    const FDM::ImmersedSurfaceSystem& system() const { return system_; }

private:
    Topology::MarkerRegistry markerRegistry_;
    Topology::CouplingGraph couplingGraph_;
    FDM::ImmersedSurfaceSystem system_;
};

} // namespace SF::IBM::Forcing
