#pragma once

/// @file SF_markerRegistry.h
/// @brief 从 STL 表面建立稳定、唯一的串行 Lagrangian marker 注册表。

#include "SF_valueTypes.h"
#include "core/model/SF_globalEntityId.h"
#include "geoProcessing/SF_STLGeometry.h"

#include <cstdint>
#include <vector>

namespace SF::IBM::Topology {

/// @brief 一个不携带流体状态的 canonical 表面 marker。
struct SurfaceMarker {
    GlobalMarkerId id;
    Vector3 position;
    Vector3 relativePosition;
    double measure = 0.0;
};

/// @brief 负责 marker identity、几何位置和积分测度，不负责插值权重。
class MarkerRegistry {
public:
    /// @brief 以 STL 三角形索引作为当前参考拓扑下的稳定 marker id。
    void build(const std::vector<GeoProcessing::Triangle>& triangles,
               const Vector3& center);

    const std::vector<SurfaceMarker>& markers() const { return markers_; }

private:
    std::vector<SurfaceMarker> markers_;
};

} // namespace SF::IBM::Topology
