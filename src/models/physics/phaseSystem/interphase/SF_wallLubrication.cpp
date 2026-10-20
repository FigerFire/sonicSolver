/// @file SF_wallLubrication.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include "core/mesh/SF_meshBoundaryGeometry.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace SF::Physics::PhaseSystems::Interphase {
namespace {

std::vector<std::array<double, 3>> wallPoints(
        const Field& field,
        const std::string& patch) {
    const auto& indices = field.getSet(patch);
    std::vector<std::array<double, 3>> points;
    points.reserve(indices.size());
    for (int index : indices) {
        int i = 0, j = 0, k = 0;
        field.getIJK(index, i, j, k);
        if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(field, i, j, k)) continue;
        points.push_back({field.X(i,j,k), field.Y(i,j,k), field.Z(i,j,k)});
    }
    if (points.empty()) {
        throw std::runtime_error(
            "Antal wallLubrication patch '" + patch
            + "' contains no physical points.");
    }
    return points;
}

double nearestDistance(
        const Field& field,
        int i,
        int j,
        int k,
        const std::vector<std::array<double, 3>>& points) {
    double distanceSquared = std::numeric_limits<double>::max();
    for (const auto& point : points) {
        const double dx = field.X(i,j,k) - point[0];
        const double dy = field.Y(i,j,k) - point[1];
        const double dz = field.Z(i,j,k) - point[2];
        distanceSquared = std::min(
            distanceSquared, dx * dx + dy * dy + dz * dz);
    }
    const double distance = std::sqrt(distanceSquared);
    if (!std::isfinite(distance) || distance < 0.0) {
        throw std::runtime_error(
            "Antal wallLubrication found invalid wall distance.");
    }
    return distance;
}

} // namespace

void addWallLubrication(
        const PairContext& pair,
        PhaseEquationSources& sources) {
    const std::string model =
        Multiphase::normalizeModelType(pair.options.wallLubricationModel);
    if (model == "none") return;
    if (model != "antal") {
        throw std::runtime_error(
            "Eulerian wallLubrication supports only none or Antal.");
    }

    const Context& context = pair.system;
    const Field& field = context.geometry;
    const auto points = wallPoints(field, pair.options.wallPatch);
    const auto outward = StructuredMesh::BoundaryGeometry::analyzeBoundarySet(
        field, pair.options.wallPatch).normal;
    const auto& continuous = context.phases[pair.continuous];
    const auto& dispersed = context.phases[pair.dispersed];
    const double configuredDiameter =
        Multiphase::phasePairDiameter(context.config, pair.options);
    const int ghost = field.NG();
    for (int k = ghost; k < ghost + field.NZ(); ++k) {
        for (int j = ghost; j < ghost + field.NY(); ++j) {
            for (int i = ghost; i < ghost + field.NX(); ++i) {
                if (field.CellFlag(i,j,k) != FLUID_CELL
                    || field.isSolverBoundaryPoint(i,j,k)) {
                    continue;
                }
                const int cell = field.getIdx(i,j,k);
                const double distance =
                    nearestDistance(field, i, j, k, points);
                const double wallDiameter =
                    sources.wallBoilingDepartureDiameter
                        .values()[(size_t)cell];
                const double diameter = wallDiameter > 0.0
                    ? wallDiameter : configuredDiameter;
                if (!std::isfinite(diameter) || diameter <= 0.0) {
                    throw std::runtime_error(
                        "Antal wallLubrication requires positive local "
                        "bubble diameter.");
                }
                // Patch自身的边界点没有控制体体积；wall-lubrication只装配
                // 到域内单元，零距离点由相边界条件闭合。
                if (distance <= 1.0e-14) continue;
                // Antal: Cwl=max(0,C1+C2*db/y), Fd=-Cwl alpha_d
                // rho_c |Ur_parallel|^2/db n_out。
                const double coefficient = std::max(
                    0.0, pair.options.wallLubricationC1
                       + pair.options.wallLubricationC2
                         * diameter / distance);
                if (coefficient == 0.0) continue;
                std::array<double, 3> relative{};
                double normalVelocity = 0.0;
                for (int component = 0; component < 3; ++component) {
                    relative[(size_t)component] =
                        dispersed.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell]
                        - continuous.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell];
                    normalVelocity +=
                        relative[(size_t)component] * outward[(size_t)component];
                }
                double tangentialSquared = 0.0;
                for (int component = 0; component < 3; ++component) {
                    const double tangential =
                        relative[(size_t)component]
                        - normalVelocity * outward[(size_t)component];
                    tangentialSquared += tangential * tangential;
                }
                const double magnitude = coefficient
                    * dispersed.primitive.alpha.values()[(size_t)cell]
                    * continuous.primitive.density.values()[(size_t)cell]
                    * tangentialSquared / diameter;
                // 分散相受力指向域内(-n_out)，连续相取反。
                std::array<double, 3> forceOnContinuous{};
                for (int component = 0; component < 3; ++component) {
                    forceOnContinuous[(size_t)component] =
                        magnitude * outward[(size_t)component];
                }
                addConservativePairForce(
                    pair, sources, cell, forceOnContinuous);
            }
        }
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
