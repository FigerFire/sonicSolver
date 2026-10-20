/// @file SF_geometrySetup.cpp
/// @brief IBM STL/BVH 几何加载阶段实现。

#include "common/SF_geometrySetup.h"

#include "common/SF_timingFormat.h"

namespace SF::IBM::Common {

GeometrySetup loadGeometry(
        GeoProcessing::STLGeometry& geometry,
        const std::vector<std::string>& stlFiles,
        const std::string& caseDir) {
    const auto start = std::chrono::steady_clock::now();
    GeometrySetup result;
    result.loaded = geometry.loadSTLFiles(stlFiles, caseDir);
    result.loadSeconds = elapsedSeconds(start);
    return result;
}

} // namespace SF::IBM::Common
