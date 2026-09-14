/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_boundaryGeometry.h
/// @brief 边界 EMPTY 维度配置；中立网格几何位于 core/mesh。

#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "core/mesh/SF_dimension.h"

#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace SF::Boundary::Geometry {

using StructuredMesh::BoundaryGeometry::BoundarySetInfo;
using StructuredMesh::BoundaryGeometry::add;
using StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet;
using StructuredMesh::BoundaryGeometry::analyzeBoundarySet;
using StructuredMesh::BoundaryGeometry::axisSize;
using StructuredMesh::BoundaryGeometry::boundaryAxisForSet;
using StructuredMesh::BoundaryGeometry::coordinateIndexForAxis;
using StructuredMesh::BoundaryGeometry::cross;
using StructuredMesh::BoundaryGeometry::dot;
using StructuredMesh::BoundaryGeometry::estimateAxisNormal;
using StructuredMesh::BoundaryGeometry::isPhysicalPoint;
using StructuredMesh::BoundaryGeometry::normalized;
using StructuredMesh::BoundaryGeometry::onAxisBoundary;
using StructuredMesh::BoundaryGeometry::point;
using StructuredMesh::BoundaryGeometry::subtract;

/// @brief 校验 EMPTY 方向仅有一个计算单元层。
inline void requireSingleCellLayerForEmpty(const Field& field,
                                           const std::string& name,
                                           int axis) {
    const int n = axisSize(field, axis);
    if (n != 2) {
        std::cerr << "[SF FATAL] EMPTY boundary set '" << name
                  << "' selects direction " << axis
                  << ", but the mesh has " << n
                  << " solution points in that direction. "
                  << "For point-based SFM, a 2D empty case must be one cell "
                  << "layer thick, i.e. exactly two solution-point planes."
                  << std::endl;
        std::exit(1);
    }
}

template <typename T>
inline void appendEmptyZones(const std::vector<BCSetting<T>>& settings,
                             std::set<std::string>& zones) {
    for (const auto& bc : settings) {
        if (bc.type == EMPTY) zones.insert(bc.name);
    }
}

inline void appendEmptyZones(const std::vector<ThermalBCSetting>& settings,
                             std::set<std::string>& zones) {
    for (const auto& bc : settings) {
        if (bc.type == ThermalBCType::Empty) zones.insert(bc.name);
    }
}

/// @brief 根据真实 EMPTY set 配置离散空间的活动方向。
///
/// EMPTY 方向由其结构化网格面推断；该函数属于边界配置，不属于中立网格查询。
template <typename... Settings>
inline void configureEmptyDimensions(const Field& field,
                                     const Settings&... settings) {
    Math::resetActiveDirections();

    std::set<std::string> emptyZones;
    (appendEmptyZones(settings, emptyZones), ...);
    if (emptyZones.empty()) return;

    int inactiveAxis = -1;
    std::array<double, 3> inactiveNormal{0.0, 0.0, 1.0};
    for (const auto& zone : emptyZones) {
        const BoundarySetInfo info = analyzeBoundarySet(field, zone);
        if (!info.valid) {
            std::cerr << "[SF FATAL] EMPTY boundary zone '" << zone
                      << "' must refer to one real mesh set located on one "
                      << "structured boundary-face direction." << std::endl;
            std::exit(1);
        }
        requireSingleCellLayerForEmpty(field, zone, info.axis);
        if (inactiveAxis >= 0 && inactiveAxis != info.axis) {
            std::cerr << "[SF FATAL] Multiple EMPTY directions are not supported "
                      << "for a 2D solve. Zones share neither the same face "
                      << "normal nor the same computational direction."
                      << std::endl;
            std::exit(1);
        }
        inactiveAxis = info.axis;
        inactiveNormal = info.normal;
    }

    Math::deactivateDirection(inactiveAxis);
    Math::setInactiveDirectionNormal(inactiveNormal);
}

} // namespace SF::Boundary::Geometry
