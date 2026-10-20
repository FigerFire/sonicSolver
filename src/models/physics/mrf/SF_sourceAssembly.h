/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_sourceAssembly.h
/// @brief MRF 平移/旋转源项共享的状态读取和守恒源装配。

#include "SF_field.h"
#include "core/residual/SF_residual.h"

#include <algorithm>

namespace SF::Source::MRF {

inline Vector3 cellPosition(const Field& field, int i, int j, int k) {
    return Vector3(field.X(i, j, k), field.Y(i, j, k), field.Z(i, j, k));
}

inline Vector3 cellVelocity(const Field& field, int i, int j, int k) {
    if (field.hasStateModel()) {
        const auto state = field.thermodynamicState(i,j,k);
        return Vector3(state.velocity[0], state.velocity[1], state.velocity[2]);
    }
    const double rho = field(i, j, k, RHO);
    return Vector3(field(i, j, k, RU) / rho,
                   field(i, j, k, RV) / rho,
                   field(i, j, k, RW) / rho);
}

inline bool appliesToZone(const Field& field, const std::string& zone,
                          int i, int j, int k) {
    const std::string setName = zone.empty() ? "all" : zone;
    const auto& set = field.getSet(setName);
    const int idx = field.getIdx(i, j, k);
    return std::find(set.begin(), set.end(), idx) != set.end();
}

inline void addBodyForceSource(Field& field, Residual& residual,
                               int i, int j, int k,
                               const Vector3& acceleration,
                               const Vector3& powerVelocity) {
    const double rho = field.hasStateModel()
        ? field.thermodynamicState(i,j,k).density
        : field(i, j, k, RHO);
    const int momentum = field.hasStateModel()
        ? field.stateModel()->momentumIndex(0) : RU;
    const int energy = field.hasStateModel()
        ? field.stateModel()->energyIndex() : E;

    residual.source(i, j, k, momentum) += rho * acceleration.x;
    residual.source(i, j, k, momentum+1) += rho * acceleration.y;
    residual.source(i, j, k, momentum+2) += rho * acceleration.z;
    residual.source(i, j, k, energy) += rho * dot(powerVelocity, acceleration);
}

inline Vector3 angularVelocity(const RotatingSetting& setting) {
    // 用户输入描述固体/STL转动；MRF源项作用在反向旋转参考系中的流体。
    return setting.axis * (-setting.omega);
}

} // namespace SF::Source::MRF
