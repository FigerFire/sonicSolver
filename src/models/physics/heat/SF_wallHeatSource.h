/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_wallHeatSource.h
/// @brief 壁面热源到守恒能量方程的源项装配。

#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "SF_field.h"
#include "core/residual/SF_residual.h"
#include "core/interfaces/SF_log.h"
#include "SF_valueTypes.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace SF {
namespace Source {
namespace WallHeat {

inline std::string normalizeMode(std::string mode) {
    mode.erase(std::remove(mode.begin(), mode.end(), '"'), mode.end());
    mode.erase(std::remove_if(mode.begin(), mode.end(),
                              [](unsigned char c) {
                                  return c == '_' || c == '-' || std::isspace(c);
                              }),
               mode.end());
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return mode;
}

inline int axisCoordinate(int i, int j, int k, int axis) {
    if (axis == 0) return i;
    if (axis == 1) return j;
    return k;
}

inline void setAxisCoordinate(int& i, int& j, int& k, int axis, int value) {
    if (axis == 0) i = value;
    else if (axis == 1) j = value;
    else k = value;
}

inline double pointDistance(const Field& field,
                            int i0, int j0, int k0,
                            int i1, int j1, int k1) {
    const double dx = field.X(i1, j1, k1) - field.X(i0, j0, k0);
    const double dy = field.Y(i1, j1, k1) - field.Y(i0, j0, k0);
    const double dz = field.Z(i1, j1, k1) - field.Z(i0, j0, k0);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline void fatalWallHeat(const std::string& message) {
    std::cerr << "[SF FATAL] WallHeat source: " << message << std::endl;
    std::exit(1);
}

/// @brief 把一个壁面热通量源项装配到贴壁第一层流体点。
/// @param field 守恒变量场和源项容器。
/// @param setting 壁面热源设置；正 heatFlux 表示向流体加热。
inline void addSource(Field& field, Residual& residual,
                      const WallHeatSetting& setting) {
    const std::string patch = setting.patch;
    if (patch.empty() || normalizeMode(patch) == "all") {
        fatalWallHeat("setting must name one boundary patch/set, not 'all'.");
    }
    if (!setting.coupleEnergy) {
        fatalWallHeat("coupleEnergy=false is not implemented; WallHeat currently "
                      "acts through the conservative energy equation.");
    }

    const std::string mode = normalizeMode(setting.mode);
    if (!(mode.empty() || mode == "fixedheatflux" || mode == "heatflux"
          || mode == "wallheat" || mode == "wallheatflux" || mode == "q")) {
        fatalWallHeat("unsupported mode '" + setting.mode
                      + "'. Supported modes: wallHeatFlux, wallHeat, fixedHeatFlux.");
    }
    if (!std::isfinite(setting.heatFlux)) {
        fatalWallHeat("non-finite heatFlux on patch '" + patch + "'.");
    }

    const int axis = StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(field, patch);
    if (StructuredMesh::BoundaryGeometry::axisSize(field, axis) < 2) {
        fatalWallHeat("patch '" + patch
                      + "' is on an active direction with fewer than two points.");
    }

    const int ng = field.NG();
    const int lo = ng;
    const int hi = ng + StructuredMesh::BoundaryGeometry::axisSize(field, axis) - 1;
    const auto& set = field.getSet(patch);

    int applied = 0;
    int skippedIntersections = 0;
    for (int id : set) {
        int i = 0, j = 0, k = 0;
        field.getIJK(id, i, j, k);
        if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(field, i, j, k)) continue;

        const int coord = axisCoordinate(i, j, k, axis);
        int targetCoord = coord;
        if (coord == lo) targetCoord = coord + 1;
        else if (coord == hi) targetCoord = coord - 1;
        else {
            fatalWallHeat("point in patch '" + patch
                          + "' is not located on its active boundary plane.");
        }

        int si = i, sj = j, sk = k;
        setAxisCoordinate(si, sj, sk, axis, targetCoord);

        if (field.isSolverBoundaryPoint(si, sj, sk)) {
            ++skippedIntersections;
            continue;
        }
        if (field.CellFlag(si, sj, sk) != FLUID_CELL) {
            fatalWallHeat("target point next to patch '" + patch
                          + "' is not a fluid cell.");
        }

        const double distance = pointDistance(field, i, j, k, si, sj, sk);
        if (!std::isfinite(distance) || distance <= 0.0) {
            fatalWallHeat("invalid wall-normal spacing on patch '" + patch + "'.");
        }

        const double volumetricHeat = setting.heatFlux / distance;
        if (!std::isfinite(volumetricHeat)) {
            fatalWallHeat("non-finite volumetric heat source on patch '" + patch + "'.");
        }
        const int energy = field.hasStateModel()
            ? field.stateModel()->energyIndex() : E;
        residual.source(si, sj, sk, energy) += volumetricHeat;
        ++applied;
    }

    if (applied == 0) {
        fatalWallHeat("patch '" + patch
                      + "' did not map to any solved fluid point.");
    }

    static bool warnedIntersections = false;
    if (skippedIntersections > 0 && !warnedIntersections) {
        warnedIntersections = true;
        broadcast("WallHeat warning: ",
                  "skipped patch-intersection points while mapping wall heat "
                  "to solved fluid cells.");
    }
}

/// @brief 装配多个壁面热源。
/// @param field 守恒变量场和源项容器。
/// @param settings 壁面热源列表。
inline void addSource(Field& field, Residual& residual,
                      const std::vector<WallHeatSetting>& settings) {
    if (settings.empty()) {
        fatalWallHeat("sourceScheme requests WallHeat but no settings were loaded.");
    }
    for (const WallHeatSetting& setting : settings) {
        addSource(field, residual, setting);
    }
}

} // namespace WallHeat
} // namespace Source
} // namespace SF
