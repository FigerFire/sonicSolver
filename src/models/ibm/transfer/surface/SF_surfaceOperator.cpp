/// @file SF_surfaceOperator.cpp
/// @brief 构造表面积分点以及归一化 J/S 插值传播权重。

#include "transfer/surface/SF_surfaceOperator.h"

#include <stdexcept>

namespace SF::IBM::Forcing {
namespace {

} // namespace

const FDM::ImmersedSurfaceSystem& SurfaceOperator::build(
        const Field& field,
        const std::vector<GeoProcessing::Triangle>& triangles,
        const Vector3& center,
        double supportRadius) {
    if (!std::isfinite(supportRadius) || supportRadius <= 0.0) {
        throw std::runtime_error(
            "SurfaceOperator requires finite positive supportRadius.");
    }
    markerRegistry_.build(triangles,center);
    couplingGraph_.build(field,markerRegistry_.markers(),supportRadius);
    // build 阶段只形成 owner raw edge；调用方必须在 Runtime 完成全局归并后调用
    // normalizeDistributed，避免多 rank 在局部 support 上错误地各自归一化。
    system_={};
    return system_;
}

const FDM::ImmersedSurfaceSystem& SurfaceOperator::normalizeDistributed(
        const std::vector<double>& normalizations) {
    couplingGraph_.normalize(normalizations);
    system_={};
    system_.points.reserve(markerRegistry_.markers().size());
    for (std::size_t markerIndex=0;
         markerIndex<markerRegistry_.markers().size(); ++markerIndex) {
        const auto& marker=markerRegistry_.markers()[markerIndex];
        FDM::ImmersedSurfacePoint point;
        point.globalConstraintId=GlobalConstraintDofId::fromMarker(marker.id);
        point.position=marker.position;
        point.relativePosition=marker.relativePosition;
        point.measure=marker.measure;
        point.interpolation=couplingGraph_.rows()[markerIndex];
        system_.points.push_back(std::move(point));
    }
    if (system_.points.empty()) {
        throw std::runtime_error(
            "SurfaceOperator generated no Lagrangian surface points.");
    }
    return system_;
}

} // namespace SF::IBM::Forcing
