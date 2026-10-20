/// @file SF_wallMapping.cpp
/// @brief RPI 壁面沸腾闭式模型与热流分配实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_wallMapping.h"

#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

namespace {

Math::Vector3 point(const Field& field, int i, int j, int k) {
    return Math::Vector3(field.X(i, j, k),
                         field.Y(i, j, k),
                         field.Z(i, j, k));
}

double coordinate(const Field& field, int component,
                  int i, int j, int k) {
    if (component == 0) return field.X(i,j,k);
    if (component == 1) return field.Y(i,j,k);
    return field.Z(i,j,k);
}

} // namespace

std::vector<WallCellSample> mapWallPatch(const ModelContext& ctx,
                                         const WallHeatSetting& setting) {
    if (setting.patch.empty()) {
        throw std::runtime_error("RPI phaseChange wallBoiling patch is empty.");
    }

    const int axis =
        StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(ctx.field,
                                                     setting.patch);
    const int ng = ctx.field.NG();
    const int lo = ng;
    const int hi = ng + StructuredMesh::BoundaryGeometry::axisSize(ctx.field, axis) - 1;
    const auto& set = ctx.field.getSet(setting.patch);
    std::vector<WallCellSample> samples;
    samples.reserve(set.size());

    for (int id : set) {
        int i = 0, j = 0, k = 0;
        ctx.field.getIJK(id, i, j, k);
        if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(ctx.field, i, j, k)) {
            continue;
        }
        if (setting.rangeCoordinate >= 0) {
            const double location = coordinate(
                ctx.field, setting.rangeCoordinate, i, j, k);
            const double scale = std::max({
                std::abs(setting.rangeMinimum),
                std::abs(setting.rangeMaximum), 1.0});
            const double tolerance = 1.0e-12 * scale;
            if (location < setting.rangeMinimum - tolerance
                || location > setting.rangeMaximum + tolerance) {
                continue;
            }
        }

        const int coord =
            StructuredMesh::BoundaryGeometry::coordinateIndexForAxis(i, j, k, axis);
        int targetCoord = coord;
        if (coord == lo) targetCoord = coord + 1;
        else if (coord == hi) targetCoord = coord - 1;
        else {
            throw std::runtime_error(
                "RPI phaseChange patch point is not on its active boundary plane.");
        }

        std::array<int, 3> cell{i, j, k};
        cell[(size_t)axis] = targetCoord;
        const int si = cell[0];
        const int sj = cell[1];
        const int sk = cell[2];
        if (ctx.field.isSolverBoundaryPoint(si, sj, sk)) continue;
        if (ctx.field.CellFlag(si, sj, sk) != FLUID_CELL) {
            throw std::runtime_error(
                "RPI phaseChange target next to wall is not a fluid cell.");
        }

        const double distance =
            Math::distance(point(ctx.field, i, j, k),
                           point(ctx.field, si, sj, sk));
        if (!std::isfinite(distance) || distance <= 0.0) {
            throw std::runtime_error(
                "RPI phaseChange invalid wall-normal spacing.");
        }

        std::array<int, 3> lower{i, j, k};
        if (coord == hi) lower[(size_t)axis] = targetCoord;
        double face[4]{0.0,0.0,0.0,0.0};
        Math::faceMetrics(
            ctx.field, lower[0], lower[1], lower[2],
            static_cast<Math::Dir>(axis), face);
        const double area = std::sqrt(
            face[0]*face[0]+face[1]*face[1]+face[2]*face[2]);
        double areaOverVolume =
            area * ctx.field.Jac(si,sj,sk);
        if (ctx.config.eulerianEulerian.axisymmetric) {
            const int radial =
                ctx.config.eulerianEulerian.radialCoordinate;
            const auto radius = [&](int pi,int pj,int pk) {
                return radial == 0 ? ctx.field.X(pi,pj,pk)
                    : (radial == 1 ? ctx.field.Y(pi,pj,pk)
                                   : ctx.field.Z(pi,pj,pk));
            };
            const double faceRadius = radius(i,j,k);
            const double cellRadius = radius(si,sj,sk);
            if (!std::isfinite(faceRadius) || faceRadius < 0.0
                || !std::isfinite(cellRadius) || cellRadius <= 0.0) {
                throw std::runtime_error(
                    "RPI axisymmetric wall mapping found invalid radius.");
            }
            areaOverVolume *= faceRadius / cellRadius;
        }
        if (!std::isfinite(areaOverVolume)
            || areaOverVolume <= 0.0) {
            throw std::runtime_error(
                "RPI phaseChange found invalid wall area/cell-volume ratio.");
        }

        samples.push_back(
            {i,j,k,si,sj,sk,distance,areaOverVolume});
    }

    if (samples.empty()) {
        throw std::runtime_error(
            "RPI phaseChange patch '" + setting.patch
            + "' did not map to any solved fluid point.");
    }
    return samples;
}

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
