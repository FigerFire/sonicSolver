/// @file SF_markerRegistry.cpp
/// @brief canonical 表面 marker 的构造与唯一性校验。

#include "topology/SF_markerRegistry.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Topology {
namespace {

Vector3 vector(const GeoProcessing::Point& point) {
    return {point.x,point.y,point.z};
}

} // namespace

void MarkerRegistry::build(
        const std::vector<GeoProcessing::Triangle>& triangles,
        const Vector3& center) {
    markers_.clear();
    markers_.reserve(triangles.size());
    for (std::size_t index=0; index<triangles.size(); ++index) {
        const Vector3 first=vector(triangles[index].v[0]);
        const Vector3 second=vector(triangles[index].v[1]);
        const Vector3 third=vector(triangles[index].v[2]);
        const double measure=0.5*norm(cross(second-first,third-first));
        if (!std::isfinite(measure) || measure<=0.0) {
            throw std::runtime_error(
                "IBM marker registry found a degenerate STL triangle at "
                "index "+std::to_string(index)+".");
        }
        SurfaceMarker marker;
        marker.id=GlobalMarkerId::fromSurfacePrimitive(
            static_cast<std::int64_t>(index));
        marker.position=(first+second+third)*(1.0/3.0);
        marker.relativePosition=marker.position-center;
        marker.measure=measure;
        markers_.push_back(marker);
    }
    if (markers_.empty()) {
        throw std::runtime_error("IBM marker registry received no triangles.");
    }
}

} // namespace SF::IBM::Topology
