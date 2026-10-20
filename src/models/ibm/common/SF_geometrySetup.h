#pragma once

/// @file SF_geometrySetup.h
/// @brief IBM STL/BVH 几何加载阶段。

#include "geoProcessing/SF_STLGeometry.h"

#include <string>
#include <vector>

namespace SF::IBM::Common {

struct GeometrySetup {
    bool loaded = false;
    double loadSeconds = 0.0;
};

GeometrySetup loadGeometry(
    GeoProcessing::STLGeometry& geometry,
    const std::vector<std::string>& stlFiles,
    const std::string& caseDir);

} // namespace SF::IBM::Common
